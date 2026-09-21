# 09_sensitivity_v3.R
# Two-tier sensitivity analyses for B-R v3:
# S1: Evaluation aggregation sensitivity (drop one test group, re-aggregate remaining 18 folds; no retraining).
# S2: True training sensitivity (truly drop one ConnGroup, rerun full outer 18-fold CV with inner tuning via shared engine, paired diff vs original).

suppressPackageStartupMessages({
  library(data.table)
})

set.seed(25)
W <- "/Users/gzy2520/Desktop/Research/kla/.worktrees/kla31-rnaseq-20260912"
setwd(W)

out_dir <- "outputs/20260921_br_v3"
preds_file <- file.path(out_dir, "BRv3_oos_preds.csv.gz")
task_file <- file.path(out_dir, "BRef_features_v3.csv.gz")

source("workflow/lactylation_ml_review/05c_metrics_and_diagnostics_v3.R")
source("workflow/lactylation_ml_review/05d_training_engine_v3.R")

P_orig <- fread(cmd = paste("gzip -dc", shQuote(preds_file)))
d_full <- fread(cmd = paste("gzip -dc", shQuote(task_file)))
d_full[, y := as.integer(y)]

MODELS <- c("M0", "MP", "MP_cal", "ME", "ME2", "ME_P", "ME2P", "XEN")
conns <- sort(unique(d_full$ConnGroup))

cat("=== S1: Evaluation Aggregation Sensitivity (Leave-One-Test-Group-Out) ===\n")

# Load existing per-fold files from main run to do exact fold re-averaging
perfold_orig_list <- list()
for (m in MODELS) {
  pf <- fread(file.path(out_dir, paste0("BRv3_perfold_", m, ".csv")))
  perfold_orig_list[[m]] <- pf
}

s1_list <- list()
for (drop_cg in conns) {
  for (m in MODELS) {
    pf_sub <- perfold_orig_list[[m]][ConnGroup != drop_cg]
    eqw_ap <- mean(pf_sub$AP, na.rm = TRUE)
    eqw_ll <- mean(pf_sub$LogLoss, na.rm = TRUE)
    eqw_brier <- mean(pf_sub$Brier, na.rm = TRUE)
    eqw_roc <- mean(pf_sub$ROC, na.rm = TRUE)
    
    wgt_ap <- weighted.mean(pf_sub$AP, pf_sub$N, na.rm = TRUE)
    wgt_ll <- weighted.mean(pf_sub$LogLoss, pf_sub$N, na.rm = TRUE)
    
    s1_list[[paste(drop_cg, m, sep = "_")]] <- data.table(
      dropped_test_group = drop_cg,
      Model = m,
      n_remaining_folds = nrow(pf_sub),
      n_remaining_rows = sum(pf_sub$N),
      eqw_AP = eqw_ap,
      eqw_LogLoss = eqw_ll,
      eqw_Brier = eqw_brier,
      eqw_ROC = eqw_roc,
      wgt_AP = wgt_ap,
      wgt_LogLoss = wgt_ll
    )
  }
}
s1_dt <- rbindlist(s1_list)
fwrite(s1_dt, file.path(out_dir, "S1_eval_sensitivity_drop_one_test_group.csv"))

cat("S1 completed. Range of equal-weight LogLoss across dropped test folds:\n")
s1_summary <- s1_dt[, .(
  min_eqw_LL = min(eqw_LogLoss),
  max_eqw_LL = max(eqw_LogLoss),
  range_eqw_LL = max(eqw_LogLoss) - min(eqw_LogLoss),
  min_eqw_AP = min(eqw_AP),
  max_eqw_AP = max(eqw_AP),
  range_eqw_AP = max(eqw_AP) - min(eqw_AP)
), by = Model]
print(s1_summary)

cat("\n=== S2: True Group Deletion Training Sensitivity ===\n")
cat("Rerunning full 18-fold outer validation with inner tuning for each of the 19 dropped ConnGroups...\n")

# Models for S2: focus on key models needed for MP baseline and module attribution
# M0, MP, MP_cal, ME_P, ME2P, ME, ME2
S2_MODELS <- c("M0", "MP", "MP_cal", "ME", "ME2", "ME_P", "ME2P")

s2_fold_records <- list()
s2_summary_records <- list()
checklist_records <- list()

t_s2_start <- Sys.time()

for (drop_group in conns) {
  cat(sprintf("-> S2 Deleting ConnGroup %02d / 19 from training dataset...\n", drop_group))
  
  # Truly delete ConnGroup from dataset
  d_s2 <- d_full[ConnGroup != drop_group]
  remaining_conns <- sort(unique(d_s2$ConnGroup))
  stopifnot(length(remaining_conns) == 18L)
  
  s2_preds_list <- list()
  
  # Run full outer validation on remaining 18 ConnGroups
  for (eval_cg in remaining_conns) {
    tr_s2 <- d_s2[ConnGroup != eval_cg]
    ev_s2 <- d_s2[ConnGroup == eval_cg]
    
    fit_res <- train_predict_cv_engine_v3(
      d_train = tr_s2,
      d_eval = ev_s2,
      models = S2_MODELS,
      lambda_grid = c(0.05, 0.1, 0.3, 1.0),
      seed = 25
    )
    
    sub_p <- fit_res$predictions
    sub_p[, ConnGroup_out := eval_cg]
    s2_preds_list[[as.character(eval_cg)]] <- sub_p
  }
  
  P_s2 <- rbindlist(s2_preds_list)
  stopifnot(nrow(P_s2) == nrow(d_s2))
  
  # Common evaluation row consistency check
  P_orig_sub <- P_orig[ConnGroup != drop_group]
  stopifnot(nrow(P_s2) == nrow(P_orig_sub))
  stopifnot(identical(P_s2$BaseAccession, P_orig_sub$BaseAccession))
  stopifnot(identical(P_s2$GroupID, P_orig_sub$GroupID))
  stopifnot(identical(P_s2$y, P_orig_sub$y))
  
  # Compute per-fold metrics for S2 and compare with original on the SAME remaining 18 folds
  for (m in S2_MODELS) {
    pcol <- paste0("p_", m)
    
    # S2 per-fold metrics
    s2_pf <- P_s2[, {
      sy <- y; sp <- get(pcol); single <- length(unique(sy)) < 2L
      .(
        dropped_train_group = drop_group,
        Model = m,
        N = length(sy),
        Pos = sum(sy),
        AP_new = if (single) NA_real_ else ap_score(sy, sp),
        LogLoss_new = logloss(sy, sp),
        Brier_new = brier(sy, sp),
        ROC_new = if (single) NA_real_ else roc_auc(sy, sp)
      )
    }, by = ConnGroup]
    
    # Orig per-fold metrics on same 18 folds
    orig_pf <- perfold_orig_list[[m]][ConnGroup != drop_group]
    
    merged_pf <- merge(s2_pf, orig_pf[, .(ConnGroup, AP_orig = AP, LogLoss_orig = LogLoss, Brier_orig = Brier, ROC_orig = ROC)], by = "ConnGroup")
    merged_pf[, diff_AP := AP_new - AP_orig]
    merged_pf[, diff_LogLoss := LogLoss_new - LogLoss_orig]
    merged_pf[, diff_Brier := Brier_new - Brier_orig]
    merged_pf[, diff_ROC := ROC_new - ROC_orig]
    
    s2_fold_records[[paste(drop_group, m, sep = "_")]] <- merged_pf
    
    # 18-fold summary comparison
    s2_summary_records[[paste(drop_group, m, sep = "_")]] <- data.table(
      dropped_train_group = drop_group,
      Model = m,
      n_remaining_folds = 18L,
      n_remaining_rows = nrow(P_s2),
      eqw_AP_orig = mean(merged_pf$AP_orig, na.rm = TRUE),
      eqw_AP_new = mean(merged_pf$AP_new, na.rm = TRUE),
      diff_eqw_AP = mean(merged_pf$AP_new, na.rm = TRUE) - mean(merged_pf$AP_orig, na.rm = TRUE),
      eqw_LL_orig = mean(merged_pf$LogLoss_orig, na.rm = TRUE),
      eqw_LL_new = mean(merged_pf$LogLoss_new, na.rm = TRUE),
      diff_eqw_LL = mean(merged_pf$LogLoss_new, na.rm = TRUE) - mean(merged_pf$LogLoss_orig, na.rm = TRUE),
      eqw_Brier_orig = mean(merged_pf$Brier_orig, na.rm = TRUE),
      eqw_Brier_new = mean(merged_pf$Brier_new, na.rm = TRUE),
      diff_eqw_Brier = mean(merged_pf$Brier_new, na.rm = TRUE) - mean(merged_pf$Brier_orig, na.rm = TRUE)
    )
  }
  
  checklist_records[[as.character(drop_group)]] <- data.table(
    dropped_train_group = drop_group,
    status = "COMPLETED",
    n_remaining_folds = 18L,
    n_remaining_rows = nrow(d_s2),
    models_evaluated = paste(S2_MODELS, collapse = ";"),
    timestamp = as.character(Sys.time())
  )
}

t_s2_end <- Sys.time()
cat(sprintf("S2 completed all 19 deletions in %.2f seconds.\n", as.numeric(difftime(t_s2_end, t_s2_start, units = "secs"))))

s2_folds_dt <- rbindlist(s2_fold_records)
s2_summary_dt <- rbindlist(s2_summary_records)
checklist_dt <- rbindlist(checklist_records)

# Save S2 artifacts
fwrite(s2_folds_dt, file.path(out_dir, "S2_train_sensitivity_perfold_diffs.csv"))
fwrite(s2_summary_dt, file.path(out_dir, "S2_train_sensitivity_summary.csv"))
fwrite(checklist_dt, file.path(out_dir, "S2_completion_checklist.csv"))

cat("\n=== S2 Training Sensitivity Summary (Mean Paired Diffs across 19 deletions) ===\n")
s2_overall <- s2_summary_dt[, .(
  mean_diff_LL = mean(diff_eqw_LL),
  max_abs_diff_LL = max(abs(diff_eqw_LL)),
  mean_diff_AP = mean(diff_eqw_AP),
  max_abs_diff_AP = max(abs(diff_eqw_AP))
), by = Model]
print(s2_overall)

cat("\nS2 completion checklist verified (19 / 19 groups complete):\n")
print(checklist_dt[, .(dropped_train_group, status, n_remaining_folds, n_remaining_rows)])
