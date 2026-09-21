# R functions: corrected metrics (threshold-aggregated AP), ridge logistic via nnls on augmented X, leakage-safe MP.
suppressPackageStartupMessages({library(data.table)})

## ---- Average Precision, threshold-aggregated (ties handled) ----
## Standard AP = sum over distinct thresholds t (desc) of (R_t - R_{t-1}) * Precision_t
ap_score <- function(y, p) {
  y <- as.numeric(y); p <- as.numeric(p)
  if (length(unique(y)) < 2L) return(NA_real_)
  P <- sum(y == 1)
  if (P == 0) return(NA_real_)
  ord <- order(p, decreasing = TRUE)
  ys <- y[ord]; ps <- p[ord]
  cs <- cumsum(ys)
  u <- sort(unique(ps), decreasing = TRUE)
  last_idx <- vapply(u, function(t) max(which(ps == t)), integer(1))  # LAST row of each tie block
  tp <- cs[last_idx]; n_at <- last_idx
  prec <- tp / n_at; rec <- tp / P
  rec0 <- c(0, rec[-length(rec)])
  sum((rec - rec0) * prec)
}
## trapezoidal PR-AUC (explicitly distinct from AP)
pr_auc_trap <- function(y, p) {
  y <- as.numeric(y); p <- as.numeric(p)
  if (length(unique(y)) < 2L) return(NA_real_)
  P <- sum(y == 1); if (P == 0) return(NA_real_)
  o <- order(p, decreasing = TRUE); ys <- y[o]; ps <- p[o]
  u <- sort(unique(ps), decreasing = TRUE)
  last_idx <- vapply(u, function(t) max(which(ps == t)), integer(1))
  tp <- cumsum(ys)[last_idx]; fp <- last_idx - tp
  prec <- tp / (tp + fp); rec <- tp / P
  prec <- c(prec[1], prec); rec <- c(0, rec)
  sum(diff(rec) * (head(prec, -1) + tail(prec, -1)) / 2)
}
roc_auc <- function(y, p) {
  if (length(unique(y)) < 2L) return(NA_real_)
  r <- rank(p); n1 <- sum(y == 1); n0 <- sum(y == 0)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
brier <- function(y, p) mean((p - y)^2)
logloss <- function(y, p, eps = 1e-15) {
  p <- pmin(pmax(p, eps), 1 - eps); -mean(y * log(p) + (1 - y) * log(1 - p))
}
## metric block, records NA with reason for single-class folds
metrics_row <- function(y, p, model, cg = NA, note = NA_character_) {
  single <- length(unique(y)) < 2L
  data.table(Model = model, ConnGroup = cg, N = length(y), Pos = sum(y),
    AP = if (single) NA_real_ else round(ap_score(y, p), 4),
    PRAUC_trap = if (single) NA_real_ else round(pr_auc_trap(y, p), 4),
    ROC = if (single) NA_real_ else round(roc_auc(y, p), 4),
    Brier = round(brier(y, p), 4),
    LogLoss = round(logloss(y, p), 4),
    Note = if (single) paste0(note, "; single-class fold -> AP/PR/ROC undefined (NA), Brier/LogLoss still valid") else note)
}

## ---- metric sanity: AP of y=c(1,1,0,0,1), p=c(.9,.8,.7,.7,.6) ----
## thresholds (last row of tie block): .9 -> n=1,tp=1 (P=1,R=1/3); .8 -> n=2,tp=2 (P=1,R=2/3);
## .7 -> n=4,tp=2 (P=0.5,R=2/3, delta=0); .6 -> n=5,tp=3 (P=0.6,R=1)
## AP = (1/3)(1) + (1/3)(1) + 0 + (1/3)(0.6) = 0.8667
stopifnot(abs(ap_score(c(1,1,0,0,1), c(.9,.8,.7,.7,.6)) - (1/3*1 + 1/3*1 + 0 + 1/3*0.6)) < 1e-12)

## ---- MP: protein detection frequency with fixed Beta(1,1) smoothing, Ref-restricted ----

## ---- module means: pre-declared modules, material-level gene medians, then mean over members ----
## Gene value per material = median over that material's samples (UNSMOOTHED log2tpm).
## Module value per material = mean over member genes; direction NOT assigned (bidirectional reactions kept unsigned).
module_mean_per_ref <- function(mat_med, member_ensembl) {
  g <- intersect(member_ensembl, rownames(mat_med))
  if (length(g) == 0) return(setNames(rep(NA_real_, ncol(mat_med)), colnames(mat_med)))
  colMeans(mat_med[g, , drop = FALSE], na.rm = TRUE)
}

## ---- z-standardization on TRAIN reference materials (unique ReferenceKey, each material weight 1) ----
## Imputation values are also train-only estimates. NOT recomputed per protein row (materials are not duplicated as weights).
fit_scaler <- function(train_refs, mat_feat_by_ref) {
  sub <- mat_feat_by_ref[, intersect(train_refs, colnames(mat_feat_by_ref)), drop = FALSE]
  mu <- rowMeans(sub, na.rm = TRUE)
  sd <- apply(sub, 1, sd, na.rm = TRUE); sd[!is.finite(sd) | sd == 0] <- 1
  impute <- apply(sub, 1, median, na.rm = TRUE)
  list(mu = mu, sd = sd, impute = impute)
}
apply_scaler <- function(scl, mat_feat_by_ref) {
  M <- mat_feat_by_ref
  for (i in seq_len(nrow(M))) { M[i, ] <- (M[i, ] - scl$mu[i]) / scl$sd[i] }
  nas <- is.na(M)
  if (any(nas)) { for (i in which(rownames(M) %in% names(scl$impute))) { M[i, nas[i, ]] <- scl$impute[i] } }
  M
}

## ---- ridge logistic fit with diagnostics (no silent fallback) ----
## L2 penalty on non-intercept coefs; Newton/IRLS inner loop in base R.
ridge_glm_v2 <- function(X, y, lambda, maxit = 200, tol = 1e-8) {
  X <- as.matrix(X)
  pen <- c(0, rep(lambda, ncol(X) - 1))
  beta <- rep(0, ncol(X)); converged <- FALSE; sep <- FALSE; it <- 0
  for (it in seq_len(maxit)) {
    eta <- as.vector(X %*% beta)
    pr <- pmin(pmax(1/(1+exp(-eta)), 1e-10), 1-1e-10)
    W <- pr * (1 - pr)
    if (max(W) < 1e-8) { sep <- TRUE; break }
    z <- eta + (y - pr) / pmax(W, 1e-12)
    A <- crossprod(X, X * W) + diag(pen)
    b <- crossprod(X, W * z)
    b_new <- tryCatch(as.vector(solve(A, b)), error = function(e) beta)
    if (max(abs(b_new - beta)) < tol) { beta <- b_new; converged <- TRUE; break }
    beta <- b_new
  }
  list(beta = beta, converged = converged, iterations = it, separation_suspected = sep, lambda = lambda)
}
pred_ridge_v2 <- function(fit, X) pmin(pmax(1/(1+exp(-(as.vector(X %*% fit$beta)))), 1e-15), 1-1e-15)

## ---- MP: protein detection frequency with fixed Beta(1,1) smoothing, Ref-restricted ----
## Returns: per-protein frequency estimated from train rows only; test rows use full outer-train estimate;
##          train rows use inner cross-fitted (leave-inner-fold) estimates to prevent self-leak.
mp_frequencies <- function(d_train, d_eval, inner_folds = NULL, a = 1, b = 1) {
  overall <- (sum(d_train$y) + a) / (nrow(d_train) + a + b)  # fallback (also the M0 rate, smoothed)
  freq <- d_train[, .(n_ref = .N, k = sum(y)), by = BaseAccession]
  freq[, freq_full := (k + a) / (n_ref + a + b)]
  freq[, freq_use := freq_full]
  freq[, method := "full_train"]
  if (!is.null(inner_folds)) {
    ## cross-fitted: for train rows in inner fold f, estimate from train rows NOT in f
    d_train[, inner_f := inner_folds]
    cf <- rbindlist(lapply(sort(unique(d_train$inner_f)), function(f) {
      tr_sub <- d_train[inner_f != f]
      ff <- tr_sub[, .(n_ref = .N, k = sum(y)), by = BaseAccession]
      ff[, freq_cf := (k + a) / (n_ref + a + b)]
      ff[, inner_f := f]
      ff
    }))
    d_train[cf, on = c("BaseAccession", "inner_f"), freq_use := i.freq_cf]
    d_train[cf, on = c("BaseAccession", "inner_f"), method := "cross_fitted"]
    d_train[, inner_f := NULL]
  }
  d_eval[, freq_use := freq[match(BaseAccession, freq$BaseAccession), freq_full]]
  d_eval[is.na(freq_use), freq_use := overall]
  d_eval[, method := ifelse(is.na(freq[match(BaseAccession, freq$BaseAccession), freq_full]), "fallback_overall", "full_train")]
  list(train = d_train, eval = d_eval, fallback_rate_eval = mean(d_eval$method == "fallback_overall"))
}
