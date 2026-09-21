# 08_train_BR_v3.R
# Main outer 19-fold cross-validation for B-R v3 using the unified, reusable training engine.
# Strictly isolated nested CV, corrected per-fold metrics, verified calibrations, and attribution controls.

suppressPackageStartupMessages({
  library(data.table)
})

set.seed(25)
W <- "/Users/gzy2520/Desktop/Research/kla/.worktrees/kla31-rnaseq-20260912"
setwd(W)

out_dir <- "outputs/20260921_br_v3"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Source core metrics and training engine
source("workflow/lactylation_ml_review/05c_metrics_and_diagnostics_v3.R")
source("workflow/lactylation_ml_review/05d_training_engine_v3.R")

# Load task dataset
task_file <- "outputs/20260922_br_v2/BRef_features_v2.csv.gz"
if (!file.exists(task_file)) {
  stop("Input features file not found: ", task_file)
}
d <- fread(cmd = paste("gzip -dc", shQuote(task_file)))
d[, y := as.integer(y)]

# Verify contract
cat(sprintf("Contract Check: %d rows, %d positives, %d ConnGroups, %d GroupIDs, %d ReferenceKeys, %d BaseAccessions\n",
            nrow(d), sum(d$y), uniqueN(d$ConnGroup), uniqueN(d$GroupID), uniqueN(d$ReferenceKey), uniqueN(d$BaseAccession)))
stopifnot(nrow(d) == 6744, sum(d$y) == 2708, uniqueN(d$ConnGroup) == 19)

# Check missingness
cat(sprintf("Missing values check: target_expr NAs = %d, module NAs = %d\n",
            sum(is.na(d$target_expr)), sum(is.na(d$mact_generation_lactate))))

# Save copy of contract-verified task and features to v3 output directory
fwrite(d, file.path(out_dir, "BRef_features_v3.csv.gz"))

MODELS <- c("M0", "MP", "MP_cal", "ME", "ME2", "ME_P", "ME2P", "XEN")
conns <- sort(unique(d$ConnGroup))

cat(sprintf("Starting 19-fold outer validation across models: %s\n", paste(MODELS, collapse = ", ")))

preds_list <- list()
diag_list <- list()
coef_list <- list()
scal_list <- list()
mp_cov_list <- list()
tuning_list <- list()

t_start <- Sys.time()

for (cg in conns) {
  cat(sprintf("--- Running Outer Fold %02d / 19 (ConnGroup %d) ---\n", cg, cg))
  
  tr <- d[ConnGroup != cg]
  ev <- d[ConnGroup == cg]
  
  engine_res <- train_predict_cv_engine_v3(
    d_train = tr,
    d_eval = ev,
    models = MODELS,
    lambda_grid = c(0.05, 0.1, 0.3, 1.0),
    seed = 25
  )
  
  # Out of fold predictions
  pred_sub <- engine_res$predictions
  pred_sub[, ConnGroup_out := cg]
  preds_list[[as.character(cg)]] <- pred_sub
  
  # Diagnostics
  diag_sub <- engine_res$diagnostics
  diag_sub[, ConnGroup := cg]
  diag_list[[as.character(cg)]] <- diag_sub
  
  # Coefficients
  coef_sub <- engine_res$coefficients
  coef_sub[, ConnGroup := cg]
  coef_list[[as.character(cg)]] <- coef_sub
  
  # Scalers
  sc <- engine_res$scalers
  scal_list[[as.character(cg)]] <- data.table(
    ConnGroup = cg,
    mu_target = sc$mu_target,
    sd_target = sc$sd_target,
    imp_target = sc$imp_target,
    mu_mod = sc$mu_mod,
    sd_mod = sc$sd_mod,
    imp_mod = sc$imp_mod,
    mu_mp = engine_res$mp_info$mu_mp,
    sd_mp = engine_res$mp_info$sd_mp
  )
  
  # MP coverage / fallback info
  mp_cov_list[[as.character(cg)]] <- data.table(
    ConnGroup = cg,
    n_train_rows = nrow(tr),
    n_eval_rows = nrow(ev),
    n_train_prot = uniqueN(tr$BaseAccession),
    n_eval_prot = uniqueN(ev$BaseAccession),
    fallback_rate_eval = engine_res$mp_info$fallback_rate_eval,
    fallback_rate_train = engine_res$mp_info$fallback_rate_train
  )
  
  # Inner tuning records
  if (nrow(engine_res$inner_tuning) > 0) {
    tun_sub <- copy(engine_res$inner_tuning)
    tun_sub[, ConnGroup := cg]
    tuning_list[[as.character(cg)]] <- tun_sub
  }
}

t_end <- Sys.time()
cat(sprintf("All 19 outer folds completed in %.2f seconds.\n", as.numeric(difftime(t_end, t_start, units = "secs"))))

# Combine outer predictions
P <- rbindlist(preds_list)
stopifnot(nrow(P) == 6744)
stopifnot(sum(P$y) == 2708)

# Check that every observation has exactly 1 outer prediction and no NAs
for (m in MODELS) {
  pcol <- paste0("p_", m)
  stopifnot(pcol %in% names(P))
  stopifnot(!anyNA(P[[pcol]]))
  stopifnot(all(P[[pcol]] >= 0 & P[[pcol]] <= 1))
}

# Save outputs
fwrite(P, file.path(out_dir, "BRv3_oos_preds.csv.gz"))
fwrite(rbindlist(diag_list), file.path(out_dir, "BRv3_fit_diagnostics.csv"))
fwrite(rbindlist(coef_list), file.path(out_dir, "BRv3_coefficients.csv"))
fwrite(rbindlist(scal_list), file.path(out_dir, "BRv3_scaler_params.csv"))
fwrite(rbindlist(mp_cov_list), file.path(out_dir, "BRv3_MP_coverage.csv"))
if (length(tuning_list) > 0) {
  fwrite(rbindlist(tuning_list), file.path(out_dir, "BRv3_inner_tuning_log.csv"))
}

cat("Outer predictions and fit artifacts saved to:", out_dir, "\n")

## ---- Per-Fold Metrics and Summaries ----

perf_summary_list <- list()

for (m in MODELS) {
  pcol <- paste0("p_", m)
  
  # True by=ConnGroup per-fold calculation
  per_fold <- P[, {
    sub_y <- y
    sub_p <- get(pcol)
    single <- length(unique(sub_y)) < 2L
    .(
      Model = m,
      N = length(sub_y),
      Pos = sum(sub_y),
      PosRate = sum(sub_y) / length(sub_y),
      AP = if (single) NA_real_ else ap_score(sub_y, sub_p),
      PRAUC_trap = if (single) NA_real_ else pr_auc_trap(sub_y, sub_p),
      ROC = if (single) NA_real_ else roc_auc(sub_y, sub_p),
      Brier = brier(sub_y, sub_p),
      LogLoss = logloss(sub_y, sub_p),
      Note = if (single) "single-class fold" else ""
    )
  }, by = ConnGroup]
  
  # Acceptance check: exactly 19 folds
  stopifnot(nrow(per_fold) == 19)
  stopifnot(sum(per_fold$N) == 6744)
  stopifnot(sum(per_fold$Pos) == 2708)
  
  # Save individual per-fold file
  fwrite(per_fold, file.path(out_dir, paste0("BRv3_perfold_", m, ".csv")))
  
  # Summaries
  eqw <- per_fold[, .(
    Model = m,
    Summary = "equal_weight_19folds",
    N_total = sum(N),
    AP = mean(AP, na.rm = TRUE),
    PRAUC_trap = mean(PRAUC_trap, na.rm = TRUE),
    ROC = mean(ROC, na.rm = TRUE),
    Brier = mean(Brier, na.rm = TRUE),
    LogLoss = mean(LogLoss, na.rm = TRUE)
  )]
  
  wgt <- per_fold[, .(
    Model = m,
    Summary = "n_weighted",
    N_total = sum(N),
    AP = weighted.mean(AP, N, na.rm = TRUE),
    PRAUC_trap = weighted.mean(PRAUC_trap, N, na.rm = TRUE),
    ROC = weighted.mean(ROC, N, na.rm = TRUE),
    Brier = weighted.mean(Brier, N, na.rm = TRUE),
    LogLoss = weighted.mean(LogLoss, N, na.rm = TRUE)
  )]
  
  y_all <- P$y
  p_all <- P[[pcol]]
  pl <- data.table(
    Model = m,
    Summary = "pooled_reference_only",
    N_total = length(y_all),
    AP = ap_score(y_all, p_all),
    PRAUC_trap = pr_auc_trap(y_all, p_all),
    ROC = roc_auc(y_all, p_all),
    Brier = brier(y_all, p_all),
    LogLoss = logloss(y_all, p_all)
  )
  
  perf_summary_list[[m]] <- rbind(eqw, wgt, pl)
}

perf_summary_dt <- rbindlist(perf_summary_list)
fwrite(perf_summary_dt, file.path(out_dir, "BRv3_performance_summary.csv"))

cat("\n=== B-R v3 Performance Summary ===\n")
cat("\n--- Primary Metric: Equal-Weight (19 folds) ---\n")
print(perf_summary_dt[Summary == "equal_weight_19folds", .(
  Model,
  LogLoss = round(LogLoss, 5),
  AP = round(AP, 5),
  Brier = round(Brier, 5),
  ROC = round(ROC, 5)
)])

cat("\n--- Secondary Metric: Group-Size Weighted ---\n")
print(perf_summary_dt[Summary == "n_weighted", .(
  Model,
  LogLoss = round(LogLoss, 5),
  AP = round(AP, 5),
  Brier = round(Brier, 5),
  ROC = round(ROC, 5)
)])

cat("\n--- Reference Only: Pooled Across All Rows ---\n")
print(perf_summary_dt[Summary == "pooled_reference_only", .(
  Model,
  LogLoss = round(LogLoss, 5),
  AP = round(AP, 5),
  Brier = round(Brier, 5),
  ROC = round(ROC, 5)
)])

## ---- Material-Level Calibration ----

calib_gid_list <- list()
calib_ref_list <- list()

for (m in MODELS) {
  pcol <- paste0("p_", m)
  
  # GroupID level (31 materials)
  gid_cal <- P[, .(
    Model = m,
    ReferenceKey = ReferenceKey[1],
    ConnGroup = ConnGroup[1],
    N_rows = .N,
    Pos_count = sum(y),
    obs_rate = mean(y),
    mean_pred = mean(get(pcol))
  ), by = GroupID]
  gid_cal[, calib_error := mean_pred - obs_rate]
  gid_cal[, abs_calib_error := abs(calib_error)]
  calib_gid_list[[m]] <- gid_cal
  
  # ReferenceKey level (28 RNA references)
  ref_cal <- P[, .(
    Model = m,
    n_group_ids = uniqueN(GroupID),
    N_rows = .N,
    Pos_count = sum(y),
    obs_rate = mean(y),
    mean_pred = mean(get(pcol))
  ), by = ReferenceKey]
  ref_cal[, calib_error := mean_pred - obs_rate]
  ref_cal[, abs_calib_error := abs(calib_error)]
  calib_ref_list[[m]] <- ref_cal
}

calib_gid_dt <- rbindlist(calib_gid_list)
calib_ref_dt <- rbindlist(calib_ref_list)

fwrite(calib_gid_dt, file.path(out_dir, "BRv3_material_calibration_groupid.csv"))
fwrite(calib_ref_dt, file.path(out_dir, "BRv3_material_calibration_referencekey.csv"))

calib_summary <- calib_ref_dt[, .(
  mean_abs_error = mean(abs_calib_error),
  max_abs_error = max(abs_calib_error),
  weighted_mae = weighted.mean(abs_calib_error, N_rows)
), by = Model]

fwrite(calib_summary, file.path(out_dir, "BRv3_material_calibration_summary.csv"))

cat("\n=== Material-Level Calibration MAE (across 28 ReferenceKeys) ===\n")
print(calib_summary[, .(Model, MAE = round(mean_abs_error, 5), MaxErr = round(max_abs_error, 5), WeightedMAE = round(weighted_mae, 5))])

cat("\nV3 TRAINING AND EVALUATION COMPLETED SUCCESSFULLY.\n")
