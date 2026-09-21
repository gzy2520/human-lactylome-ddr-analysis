# 05d_training_engine_v4.R
# Reusable, strictly isolated training and evaluation engine for B-R v4.
# Shared across main outer cross-validation, inner tuning, and S2 group deletion sensitivity.

suppressPackageStartupMessages({
  library(data.table)
})

source("workflow/lactylation_ml_review/05c_metrics_and_diagnostics_v3.R")

## ---- Feature Builder Helper ----

#' Build design matrix X with column names and intercept
#' @param dt Data table containing standardized features
#' @param spec Character vector of feature names (e.g., "z_target_expr", "z_mod", "z_mp", "z_ix")
build_design_matrix_v4 <- function(dt, spec) {
  n <- nrow(dt)
  X <- matrix(1, nrow = n, ncol = 1, dimnames = list(NULL, "(Intercept)"))
  if (length(spec) == 0L) return(X)

  cols <- list()
  for (f in spec) {
    if (f == "z_ix") {
      if (!all(c("z_target_expr", "z_mod") %in% names(dt))) {
        stop("z_ix requires both z_target_expr and z_mod in table")
      }
      cols[[f]] <- dt$z_target_expr * dt$z_mod
    } else {
      if (!f %in% names(dt)) {
        stop("Feature ", f, " not found in data table")
      }
      cols[[f]] <- dt[[f]]
    }
  }
  X_mat <- cbind(X, as.matrix(as.data.table(cols)))
  colnames(X_mat) <- c("(Intercept)", spec)
  X_mat
}

## ---- Preprocessing Scaler Helper ----

#' Compute material-equal weighted mean and sd for continuous features
#' Row weight w_i = 1 / n_{m(i)} such that each material contributes total weight 1/M.
calc_weighted_moments <- function(vals, mat_keys) {
  dt <- data.table(val = as.numeric(vals), mat = as.character(mat_keys))
  mat_stats <- dt[, .(
    n = .N,
    m_mean = mean(val, na.rm = TRUE),
    m_var = if (.N > 1) mean((val - mean(val, na.rm = TRUE))^2) else 0
  ), by = mat]
  M <- nrow(mat_stats)
  if (M == 0) return(list(mu = 0, sd = 1))
  mu <- mean(mat_stats$m_mean, na.rm = TRUE)
  var_tot <- mean(mat_stats$m_var + (mat_stats$m_mean - mu)^2, na.rm = TRUE)
  sd_val <- sqrt(pmax(var_tot, 1e-12))
  if (!is.finite(sd_val) || sd_val < 1e-8) sd_val <- 1
  list(mu = mu, sd = sd_val)
}

#' Fit scalers on training data table
#' - target_expr: material-equal weighted mean and sd
#' - mact_generation_lactate: mean and sd across unique reference materials (each material weight 1)
fit_scalers_v4 <- function(dt_train, mod_col = "mact_generation_lactate") {
  # target_expr
  te_moments <- calc_weighted_moments(dt_train$target_expr, dt_train$ReferenceKey)
  imp_target <- median(dt_train$target_expr, na.rm = TRUE)
  if (!is.finite(imp_target)) imp_target <- te_moments$mu

  # module feature (material-level)
  mod_dt <- unique(dt_train[, .(ReferenceKey, mod_val = get(mod_col))])
  mu_mod <- mean(mod_dt$mod_val, na.rm = TRUE)
  sd_mod <- sd(mod_dt$mod_val, na.rm = TRUE)
  if (!is.finite(sd_mod) || sd_mod < 1e-8) sd_mod <- 1
  imp_mod <- median(mod_dt$mod_val, na.rm = TRUE)
  if (!is.finite(imp_mod)) imp_mod <- mu_mod

  list(
    mu_target = te_moments$mu,
    sd_target = te_moments$sd,
    imp_target = imp_target,
    mu_mod = mu_mod,
    sd_mod = sd_mod,
    imp_mod = imp_mod,
    n_imputed_target_train = sum(is.na(dt_train$target_expr)),
    n_imputed_mod_train = sum(is.na(dt_train[[mod_col]]))
  )
}

#' Apply scalers to transform features to z-scores
apply_scalers_v4 <- function(dt, scalers, mod_col = "mact_generation_lactate") {
  dt_out <- copy(dt)

  # target_expr
  raw_t <- dt_out$target_expr
  raw_t[is.na(raw_t)] <- scalers$imp_target
  dt_out[, z_target_expr := (raw_t - scalers$mu_target) / scalers$sd_target]

  # mod
  raw_m <- dt_out[[mod_col]]
  raw_m[is.na(raw_m)] <- scalers$imp_mod
  dt_out[, z_mod := (raw_m - scalers$mu_mod) / scalers$sd_mod]

  dt_out
}

## ---- Strictly Isolated MP Target Encoding Helper ----

#' Compute isolated MP frequencies and standardizations
#' @param dt_train Training data table
#' @param dt_eval Evaluation data table (labels must NOT be used!)
#' @param inner_conns_assign List of ConnGroups mapping to inner folds for dt_train cross-fitting
calc_isolated_mp_v4 <- function(dt_train, dt_eval, inner_assign) {
  # 1. Outer test set eval predictions: estimated from ALL training rows
  train_counts <- dt_train[, .(n = .N, k = sum(y)), by = BaseAccession]
  train_counts[, freq_train := (k + 1) / (n + 2)]

  overall_train_rate <- (sum(dt_train$y) + 1) / (nrow(dt_train) + 2)

  eval_matched <- train_counts[match(dt_eval$BaseAccession, train_counts$BaseAccession), freq_train]
  eval_is_fallback <- is.na(eval_matched)
  eval_matched[eval_is_fallback] <- overall_train_rate

  dt_eval_out <- copy(dt_eval)
  dt_eval_out[, freq_use := eval_matched]
  dt_eval_out[, mp_fallback := eval_is_fallback]

  # 2. Training rows cross-fitted strictly across inner folds of dt_train
  dt_train_out <- copy(dt_train)
  dt_train_out[, freq_use := NA_real_]
  dt_train_out[, mp_fallback := FALSE]

  # Assign inner fold to dt_train rows
  inner_f_vec <- integer(nrow(dt_train_out))
  for (f in seq_along(inner_assign)) {
    inner_f_vec[dt_train_out$ConnGroup %in% inner_assign[[f]]] <- f
  }
  dt_train_out[, inner_f := inner_f_vec]

  n_inner <- length(inner_assign)
  for (f in seq_len(n_inner)) {
    tr_sub <- dt_train_out[inner_f != f]
    ev_idx <- which(dt_train_out$inner_f == f)

    sub_counts <- tr_sub[, .(n = .N, k = sum(y)), by = BaseAccession]
    sub_counts[, freq_sub := (k + 1) / (n + 2)]
    sub_overall <- (sum(tr_sub$y) + 1) / (nrow(tr_sub) + 2)

    f_matched <- sub_counts[match(dt_train_out$BaseAccession[ev_idx], sub_counts$BaseAccession), freq_sub]
    f_fallback <- is.na(f_matched)
    f_matched[f_fallback] <- sub_overall

    dt_train_out[ev_idx, freq_use := f_matched]
    dt_train_out[ev_idx, mp_fallback := f_fallback]
  }
  dt_train_out[, inner_f := NULL]

  # 3. Standardize z_mp using material-equal weighted moments of training cross-fitted freq_use
  mp_moments <- calc_weighted_moments(dt_train_out$freq_use, dt_train_out$ReferenceKey)
  dt_train_out[, z_mp := (freq_use - mp_moments$mu) / mp_moments$sd]
  dt_eval_out[, z_mp := (freq_use - mp_moments$mu) / mp_moments$sd]

  list(
    train = dt_train_out,
    eval = dt_eval_out,
    mu_mp = mp_moments$mu,
    sd_mp = mp_moments$sd,
    fallback_rate_eval = mean(eval_is_fallback),
    fallback_rate_train = mean(dt_train_out$mp_fallback)
  )
}

## ---- Reusable Core Training and Inner Tuning Engine ----

#' Reusable Training and Prediction Function
#'
#' @param d_train Training data table (contains GroupID, ReferenceKey, ConnGroup, BaseAccession, EnsemblGeneID, y, target_expr, mact_generation_lactate)
#' @param d_eval Evaluation data table to predict on (labels in d_eval are NEVER used during fitting/predicting!)
#' @param models Vector of models to run (subset of: M0, MP, MP_cal, ME, ME2, ME_P, ME2P, XEN)
#' @param lambda_grid Candidate lambdas for ridge regression (default c(0.05, 0.1, 0.3, 1.0))
#' @param seed Random seed for inner partition (default 25)
#' @return List containing predictions on d_eval, diagnostics, fitted coefficients, scalers, selected lambdas
train_predict_cv_engine_v4 <- function(d_train,
                                       d_eval,
                                       models = c("M0", "MP", "MP_cal", "ME", "ME2", "ME_P", "ME2P", "XEN"),
                                       lambda_grid = c(0.05, 0.1, 0.3, 1.0),
                                       seed = 25) {

  # Ensure no side effects on input data tables
  tr <- copy(d_train)
  ev <- copy(d_eval)

  # Feature specifications
  model_specs <- list(
    M0 = character(0),
    MP = character(0),
    MP_cal = "z_mp",
    ME = "z_target_expr",
    ME2 = c("z_target_expr", "z_mod"),
    ME_P = c("z_target_expr", "z_mp"),
    ME2P = c("z_target_expr", "z_mod", "z_mp"),
    XEN = c("z_target_expr", "z_mod", "z_ix")
  )

  # Partition training ConnGroups into 5 inner folds
  inner_conns <- sort(unique(tr$ConnGroup))
  set.seed(seed)
  shuffled_conns <- sample(inner_conns)
  n_folds <- min(5L, length(inner_conns))
  inner_assign <- split(shuffled_conns, cut(seq_along(shuffled_conns), n_folds, labels = FALSE))

  # 1. Inner Hyperparameter Tuning (with strict isolation per inner fold)
  # Ridge models requiring tuning:
  ridge_models <- intersect(models, c("MP_cal", "ME", "ME2", "ME_P", "ME2P", "XEN"))
  selected_lambdas <- list()
  inner_perf_list <- list()
  inner_group_log <- list()

  if (length(ridge_models) > 0L) {
    for (m in ridge_models) {
      spec <- model_specs[[m]]
      # If model has only 1 lambda in grid, no tuning needed
      if (length(lambda_grid) == 1L) {
        selected_lambdas[[m]] <- lambda_grid[1]
        next
      }

      # Evaluate each candidate lambda across inner folds
      loss_by_lam <- numeric(length(lambda_grid))
      names(loss_by_lam) <- as.character(lambda_grid)

      for (lam_idx in seq_along(lambda_grid)) {
        lam_cand <- lambda_grid[lam_idx]
        fold_losses <- vector("list", n_folds)

        for (f in seq_len(n_folds)) {
          # Split inner train A and inner val V
          val_conns <- inner_assign[[f]]
          A_train <- tr[!ConnGroup %in% val_conns]
          V_val <- tr[ConnGroup %in% val_conns]

          # Fit scalers strictly on A
          scalers_A <- fit_scalers_v4(A_train)
          A_scaled <- apply_scalers_v4(A_train, scalers_A)
          V_scaled <- apply_scalers_v4(V_val, scalers_A)

          # Inner MP encoding: V uses A's labels only; A cross-fitted within A
          # Sub-partition A's ConnGroups for A's internal cross-fitting
          A_conns <- sort(unique(A_train$ConnGroup))
          n_sub <- min(4L, length(A_conns))
          set.seed(seed + f)
          A_sub_assign <- split(sample(A_conns), cut(seq_along(A_conns), n_sub, labels = FALSE))

          mp_inner <- calc_isolated_mp_v4(A_scaled, V_scaled, A_sub_assign)
          A_ready <- mp_inner$train
          V_ready <- mp_inner$eval

          # Build design matrices
          XA <- build_design_matrix_v4(A_ready, spec)
          XV <- build_design_matrix_v4(V_ready, spec)

          # Fit model on A
          fit_inner <- ridge_glm_v3(XA, A_ready$y, lam_cand)
          if (fit_inner$failed || !fit_inner$converged) {
            stop("Inner fit failed; cannot select lambda from incomplete fits")
          } else {
            preds_V <- pred_ridge_v3(fit_inner, XV)
            fold_losses[[f]] <- data.table(ConnGroup = V_ready$ConnGroup,
              row_loss = -V_ready$y * log(preds_V) - (1 - V_ready$y) * log1p(-preds_V))[,
              .(LogLoss = mean(row_loss)), by = ConnGroup]
          }
        }
        # Every validation connected group contributes exactly once, irrespective of inner-fold sizes.
        inner_group_losses <- rbindlist(fold_losses)
        stopifnot(!anyDuplicated(inner_group_losses$ConnGroup),
                  setequal(inner_group_losses$ConnGroup, inner_conns))
        loss_by_lam[lam_idx] <- mean(inner_group_losses$LogLoss)
        inner_group_losses[, `:=`(Model = m, candidate_lambda = lam_cand)]
        inner_group_log[[paste(m, lam_cand)]] <- inner_group_losses
      }

      best_lam <- as.numeric(names(which.min(loss_by_lam)))
      selected_lambdas[[m]] <- best_lam
      inner_perf_list[[m]] <- data.table(
        Model = m,
        candidate_lambda = as.numeric(names(loss_by_lam)),
        inner_mean_logloss = loss_by_lam,
        selected = (as.numeric(names(loss_by_lam)) == best_lam)
      )
    }
  }

  # 2. Final Fit on ALL Outer Training Data (d_train)
  # Fit preprocessing scalers on ALL d_train
  outer_scalers <- fit_scalers_v4(tr)
  tr_scaled <- apply_scalers_v4(tr, outer_scalers)
  ev_scaled <- apply_scalers_v4(ev, outer_scalers)

  # MP cross-fitting on ALL d_train across outer inner_assign; ev estimated from ALL d_train
  outer_mp <- calc_isolated_mp_v4(tr_scaled, ev_scaled, inner_assign)
  tr_final <- outer_mp$train
  ev_final <- outer_mp$eval

  # 3. Predict on ev for all requested models
  preds_dt <- copy(ev_final[, .(GroupID, ReferenceKey, ConnGroup, PXD, BaseAccession, EnsemblGeneID, y)])
  diag_list <- list()
  coef_list <- list()

  # M0: Prevalence baseline
  if ("M0" %in% models) {
    p_m0 <- mean(tr$y)
    preds_dt[, p_M0 := p_m0]
  }

  # MP: Direct protein frequency baseline (heuristic smoothed probability)
  if ("MP" %in% models) {
    preds_dt[, p_MP := ev_final$freq_use]
  }

  # Ridge models: MP_cal, ME, ME2, ME_P, ME2P, XEN
  for (m in ridge_models) {
    spec <- model_specs[[m]]
    lam <- selected_lambdas[[m]]

    X_tr <- build_design_matrix_v4(tr_final, spec)
    X_ev <- build_design_matrix_v4(ev_final, spec)

    fit <- ridge_glm_v3(X_tr, tr_final$y, lam)

    if (fit$failed || !fit$converged) {
      stop(sprintf("Model %s failed to converge on training data: %s", m, fit$fail_reason))
    }

    preds_dt[, (paste0("p_", m)) := pred_ridge_v3(fit, X_ev)]

    diag_list[[m]] <- data.table(
      Model = m,
      selected_lambda = lam,
      converged = fit$converged,
      failed = fit$failed,
      fail_reason = fit$fail_reason,
      iterations = fit$iterations,
      irls_weight_underflow_detected = fit$irls_weight_underflow_detected,
      n_train_rows = nrow(X_tr),
      n_features = ncol(X_tr) - 1L,
      fallback_rate_eval = outer_mp$fallback_rate_eval
    )

    coef_list[[m]] <- data.table(
      Model = m,
      selected_lambda = lam,
      Feature = colnames(X_tr),
      Coefficient = fit$beta
    )
  }

  list(
    predictions = preds_dt,
    diagnostics = rbindlist(diag_list),
    coefficients = rbindlist(coef_list),
    scalers = outer_scalers,
    mp_info = list(
      mu_mp = outer_mp$mu_mp,
      sd_mp = outer_mp$sd_mp,
      fallback_rate_eval = outer_mp$fallback_rate_eval,
      fallback_rate_train = outer_mp$fallback_rate_train
    ),
    inner_tuning = if (length(inner_perf_list) > 0) rbindlist(inner_perf_list) else data.table(),
    inner_group_losses = rbindlist(inner_group_log),
    selected_lambdas = selected_lambdas
  )
}
