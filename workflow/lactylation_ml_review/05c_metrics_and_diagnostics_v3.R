# 05c_metrics_and_diagnostics_v3.R
# Verified metric functions, tie-safe AP, trapezoidal PR-AUC, robust ridge logistic regression solver,
# and convergence diagnostics.

suppressPackageStartupMessages({
  library(data.table)
})

## ---- 1. Metrics ----

#' Threshold-aggregated Average Precision (tie-order invariant)
#' AP = sum_{t} (R_t - R_{t-1}) * Precision_t over distinct prediction thresholds
ap_score <- function(y, p) {
  y <- as.numeric(y); p <- as.numeric(p)
  if (length(y) == 0L) return(NA_real_)
  if (length(unique(y)) < 2L) return(NA_real_)
  P <- sum(y == 1)
  if (P == 0L || P == length(y)) {
    # If all 0s or all 1s, AP is undefined (single-class)
    return(NA_real_)
  }
  ord <- order(p, decreasing = TRUE)
  ys <- y[ord]; ps <- p[ord]
  cs <- cumsum(ys)
  u <- sort(unique(ps), decreasing = TRUE)
  # End of each tie block
  last_idx <- vapply(u, function(t) max(which(ps == t)), integer(1))
  tp <- cs[last_idx]
  n_at <- last_idx
  prec <- tp / n_at
  rec <- tp / P
  rec0 <- c(0, rec[-length(rec)])
  sum((rec - rec0) * prec)
}

#' Trapezoidal Area under the Precision-Recall Curve (PR-AUC)
#' Kept strictly separate and distinctly named from threshold-aggregated AP
pr_auc_trap <- function(y, p) {
  y <- as.numeric(y); p <- as.numeric(p)
  if (length(y) == 0L) return(NA_real_)
  if (length(unique(y)) < 2L) return(NA_real_)
  P <- sum(y == 1)
  if (P == 0L || P == length(y)) return(NA_real_)
  ord <- order(p, decreasing = TRUE)
  ys <- y[ord]; ps <- p[ord]
  u <- sort(unique(ps), decreasing = TRUE)
  last_idx <- vapply(u, function(t) max(which(ps == t)), integer(1))
  tp <- cumsum(ys)[last_idx]
  fp <- last_idx - tp
  prec <- tp / (tp + fp)
  rec <- tp / P
  prec_curve <- c(prec[1], prec)
  rec_curve <- c(0, rec)
  sum(diff(rec_curve) * (head(prec_curve, -1) + tail(prec_curve, -1)) / 2)
}

#' Area under the ROC Curve (Mann-Whitney U statistic with average tie ranks)
roc_auc <- function(y, p) {
  y <- as.numeric(y); p <- as.numeric(p)
  if (length(y) == 0L) return(NA_real_)
  if (length(unique(y)) < 2L) return(NA_real_)
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  if (n1 == 0L || n0 == 0L) return(NA_real_)
  r <- rank(p, ties.method = "average")
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

#' Brier Score (Mean Squared Error on probabilities)
brier <- function(y, p) {
  mean((as.numeric(p) - as.numeric(y))^2)
}

#' Log-Loss / Binary Cross-Entropy (clamped to [eps, 1 - eps])
logloss <- function(y, p, eps = 1e-15) {
  p_clamp <- pmin(pmax(as.numeric(p), eps), 1 - eps)
  y_num <- as.numeric(y)
  -mean(y_num * log(p_clamp) + (1 - y_num) * log(1 - p_clamp))
}

#' Standardized single-fold metric row generator
metrics_row_v3 <- function(y, p, model, cg = NA, note = "") {
  y <- as.numeric(y); p <- as.numeric(p)
  n <- length(y); pos <- sum(y == 1)
  single <- length(unique(y)) < 2L
  data.table(
    Model = model,
    ConnGroup = cg,
    N = n,
    Pos = pos,
    PosRate = if (n > 0) pos / n else NA_real_,
    AP = if (single) NA_real_ else ap_score(y, p),
    PRAUC_trap = if (single) NA_real_ else pr_auc_trap(y, p),
    ROC = if (single) NA_real_ else roc_auc(y, p),
    Brier = brier(y, p),
    LogLoss = logloss(y, p),
    Note = if (single) paste0(note, "; single-class fold: AP/PRAUC/ROC undefined (NA), Brier/LogLoss valid") else note
  )
}

## ---- 2. Robust Ridge Logistic Regression Solver ----

#' Robust Newton-Raphson / IRLS for L2-penalized Logistic Regression
#' Objective: -sum_{i=1}^n [ y_i * eta_i - log(1 + exp(eta_i)) ] + 0.5 * lambda * sum(beta[-1]^2)
#'
#' @param X Design matrix with intercept in column 1
#' @param y Binary response vector (0 or 1)
#' @param lambda Non-negative L2 penalty on non-intercept coefficients
#' @param maxit Maximum iterations (default 200)
#' @param tol Convergence tolerance on maximum absolute change in coefficients (default 1e-8)
#' @return List containing beta, converged, failed, fail_reason, iterations, irls_weight_underflow_detected, lambda
ridge_glm_v3 <- function(X, y, lambda, maxit = 200, tol = 1e-8) {
  X <- as.matrix(X)
  y <- as.numeric(y)
  p_cols <- ncol(X)
  pen <- c(0, rep(lambda, p_cols - 1))
  
  beta <- rep(0, p_cols)
  converged <- FALSE
  failed <- FALSE
  fail_reason <- NA_character_
  underflow_detected <- FALSE
  it <- 0
  
  for (it in seq_len(maxit)) {
    eta <- as.vector(X %*% beta)
    pr <- pmin(pmax(1 / (1 + exp(-eta)), 1e-12), 1 - 1e-12)
    W <- pr * (1 - pr)
    
    if (all(W < 1e-8)) {
      underflow_detected <- TRUE
      failed <- TRUE
      fail_reason <- "all_irls_weights_underflow"
      break
    }
    
    # Working response z = eta + (y - pr) / W
    z <- eta + (y - pr) / pmax(W, 1e-12)
    A <- crossprod(X, X * W) + diag(pen, nrow = p_cols, ncol = p_cols)
    b <- crossprod(X, W * z)
    
    # Invert system safely; DO NOT silently catch and return old coefficients!
    sol <- tryCatch(
      as.vector(solve(A, b)),
      error = function(e) {
        list(err = e$message)
      }
    )
    
    if (is.list(sol) && !is.null(sol$err)) {
      failed <- TRUE
      fail_reason <- paste0("linear_system_solve_error: ", sol$err)
      converged <- FALSE
      break
    }
    
    b_new <- sol
    if (!all(is.finite(b_new))) {
      failed <- TRUE
      fail_reason <- "non_finite_coefficients_encountered"
      converged <- FALSE
      break
    }
    
    if (max(abs(b_new - beta)) < tol) {
      beta <- b_new
      converged <- TRUE
      break
    }
    beta <- b_new
  }
  
  if (!converged && !failed && it >= maxit) {
    failed <- TRUE
    fail_reason <- "maximum_iterations_exceeded"
  }
  
  list(
    beta = beta,
    converged = converged,
    failed = failed,
    fail_reason = fail_reason,
    iterations = it,
    irls_weight_underflow_detected = underflow_detected,
    lambda = lambda
  )
}

#' Predict probabilities from fitted ridge logistic model
pred_ridge_v3 <- function(fit, X, eps = 1e-15) {
  if (fit$failed || !fit$converged) {
    stop("Cannot predict with non-converged or failed model: ", fit$fail_reason)
  }
  eta <- as.vector(as.matrix(X) %*% fit$beta)
  pmin(pmax(1 / (1 + exp(-eta)), eps), 1 - eps)
}
