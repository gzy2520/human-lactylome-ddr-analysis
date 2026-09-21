# 05b_v2_independent_recalculation.R
# Independent recalculation of B-R v2 outer predictions (without retraining)
# Identifies exact aggregation and calibration errors in v2 reporting.

suppressPackageStartupMessages({
  library(data.table)
})

set.seed(25)
W <- "/Users/gzy2520/Desktop/Research/kla/.worktrees/kla31-rnaseq-20260912"
setwd(W)

out_v2 <- "outputs/20260922_br_v2"
out_audit <- "audit/20260921_br_v3"
dir.create(out_audit, showWarnings = FALSE, recursive = TRUE)

preds_file <- file.path(out_v2, "BRv2_oos_preds.csv.gz")
if (!file.exists(preds_file)) {
  stop("v2 predictions file not found: ", preds_file)
}

# Read via gzip pipe to avoid R.utils dependency
p <- fread(cmd = paste("gzip -dc", shQuote(preds_file)))

cat(sprintf("Loaded v2 predictions: %d rows, %d positives, %d ConnGroups, %d GroupIDs, %d ReferenceKeys\n",
            nrow(p), sum(p$y), uniqueN(p$ConnGroup), uniqueN(p$GroupID), uniqueN(p$ReferenceKey)))

# Source original metric definitions from 05_RR_functions.R
source("workflow/lactylation_ml_review/05_RR_functions.R")

models <- c("M0", "MP", "ME", "ME2", "ME2P", "XEN")

# 1. Per-ConnGroup metrics
perfold_list <- list()
for (m in models) {
  pred_col <- paste0("p_", m)
  for (cg in sort(unique(p$ConnGroup))) {
    sub <- p[ConnGroup == cg]
    y_sub <- sub$y
    p_sub <- sub[[pred_col]]
    single <- length(unique(y_sub)) < 2L
    
    perfold_list[[paste(m, cg, sep = "_")]] <- data.table(
      Model = m,
      ConnGroup = cg,
      N = length(y_sub),
      Pos = sum(y_sub),
      PosRate = sum(y_sub) / length(y_sub),
      AP = if (single) NA_real_ else ap_score(y_sub, p_sub),
      PRAUC_trap = if (single) NA_real_ else pr_auc_trap(y_sub, p_sub),
      ROC = if (single) NA_real_ else roc_auc(y_sub, p_sub),
      Brier = brier(y_sub, p_sub),
      LogLoss = logloss(y_sub, p_sub)
    )
  }
}
perfold_dt <- rbindlist(perfold_list)

# 2. Aggregations: Unweighted (equal weight 19 folds), Group-size weighted, and Pooled
summary_list <- list()
for (m in models) {
  sub_pf <- perfold_dt[Model == m]
  pred_col <- paste0("p_", m)
  y_all <- p$y
  p_all <- p[[pred_col]]
  
  # Equal weight
  eqw <- data.table(
    Model = m,
    Summary = "equal_weight_19folds",
    N_total = sum(sub_pf$N),
    AP = mean(sub_pf$AP, na.rm = TRUE),
    PRAUC_trap = mean(sub_pf$PRAUC_trap, na.rm = TRUE),
    ROC = mean(sub_pf$ROC, na.rm = TRUE),
    Brier = mean(sub_pf$Brier, na.rm = TRUE),
    LogLoss = mean(sub_pf$LogLoss, na.rm = TRUE)
  )
  
  # N weighted
  wgt <- data.table(
    Model = m,
    Summary = "n_weighted",
    N_total = sum(sub_pf$N),
    AP = weighted.mean(sub_pf$AP, sub_pf$N, na.rm = TRUE),
    PRAUC_trap = weighted.mean(sub_pf$PRAUC_trap, sub_pf$N, na.rm = TRUE),
    ROC = weighted.mean(sub_pf$ROC, sub_pf$N, na.rm = TRUE),
    Brier = weighted.mean(sub_pf$Brier, sub_pf$N, na.rm = TRUE),
    LogLoss = weighted.mean(sub_pf$LogLoss, sub_pf$N, na.rm = TRUE)
  )
  
  # Pooled (reference only)
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
  
  summary_list[[m]] <- rbind(eqw, wgt, pl)
}
summary_dt <- rbindlist(summary_list)

# 3. Calibration per GroupID (31 materials)
calib_groupid_list <- list()
for (m in models) {
  pred_col <- paste0("p_", m)
  cg_tab <- p[, .(
    Model = m,
    GroupID = GroupID[1],
    ReferenceKey = ReferenceKey[1],
    ConnGroup = ConnGroup[1],
    N_rows = .N,
    Pos_count = sum(y),
    obs_rate = mean(y),
    mean_pred = mean(get(pred_col))
  ), by = GroupID]
  cg_tab[, calib_error := mean_pred - obs_rate]
  cg_tab[, abs_calib_error := abs(calib_error)]
  calib_groupid_list[[m]] <- cg_tab
}
calib_groupid_dt <- rbindlist(calib_groupid_list)

# 4. Calibration per ReferenceKey (28 RNA references)
# Note: HCT116 has 3 GroupIDs sharing RNA reference; HK-2 has 2 GroupIDs sharing RNA reference.
calib_ref_list <- list()
for (m in models) {
  pred_col <- paste0("p_", m)
  cr_tab <- p[, .(
    Model = m,
    ReferenceKey = ReferenceKey[1],
    n_group_ids = uniqueN(GroupID),
    N_rows = .N,
    Pos_count = sum(y),
    obs_rate = mean(y),
    mean_pred = mean(get(pred_col))
  ), by = ReferenceKey]
  cr_tab[, calib_error := mean_pred - obs_rate]
  cr_tab[, abs_calib_error := abs(calib_error)]
  calib_ref_list[[m]] <- cr_tab
}
calib_ref_dt <- rbindlist(calib_ref_list)

# Summary of Calibration MAE across ReferenceKey
calib_ref_summary <- calib_ref_dt[, .(
    mean_abs_error = mean(abs_calib_error),
    max_abs_error = max(abs_calib_error),
    weighted_mae = weighted.mean(abs_calib_error, N_rows)
), by = Model]

cat("\n=== Independent Recalculation Results for v2 Predictions ===\n")
cat("\n--- Equal-Weight (19 folds) ---\n")
print(summary_dt[Summary == "equal_weight_19folds", .(Model, AP = round(AP, 7), LogLoss = round(LogLoss, 7), Brier = round(Brier, 7), ROC = round(ROC, 7))])

cat("\n--- Group-Size Weighted ---\n")
print(summary_dt[Summary == "n_weighted", .(Model, AP = round(AP, 7), LogLoss = round(LogLoss, 7), Brier = round(Brier, 7), ROC = round(ROC, 7))])

cat("\n--- Pooled (Reference Only) ---\n")
print(summary_dt[Summary == "pooled_reference_only", .(Model, AP = round(AP, 7), LogLoss = round(LogLoss, 7), Brier = round(Brier, 7), ROC = round(ROC, 7))])

cat("\n--- Calibration MAE per ReferenceKey ---\n")
print(calib_ref_summary)

# 5. Build Errata Comparison Table: Old Reported vs Correct Recalculated vs Reason
errata <- data.table(
  Model = c("M0", "M0", "M0", "MP", "MP", "MP", "ME", "ME", "ME2", "ME2", "ME2P", "ME2P", "XEN", "XEN",
            "M0_calib", "MP_calib", "ME_calib", "ME2P_calib"),
  Metric = c("Equal-weight AP", "Equal-weight LogLoss", "Equal-weight ROC",
             "Equal-weight AP", "Equal-weight LogLoss", "Equal-weight ROC",
             "Equal-weight AP", "Equal-weight LogLoss",
             "Equal-weight AP", "Equal-weight LogLoss",
             "Equal-weight AP", "Equal-weight LogLoss",
             "Equal-weight AP", "Equal-weight LogLoss",
             "Material Calib MAE", "Material Calib MAE", "Material Calib MAE", "Material Calib MAE"),
  Old_Reported_Value = c("0.3216", "0.6826", "0.3205",
                         "0.6583", "0.5757", "0.7556",
                         "0.4313", "0.6843",
                         "0.4091", "0.6949",
                         "0.6134", "0.6026",
                         "0.3787", "0.7147",
                         "0.0004", "0.0137", "0.0009", "0.0029"),
  Old_Scope_In_Reality = c("Pooled (6744 rows)", "Pooled (6744 rows)", "Pooled (6744 rows)",
                           "Pooled (6744 rows)", "Pooled (6744 rows)", "Pooled (6744 rows)",
                           "Pooled (6744 rows)", "Pooled (6744 rows)",
                           "Pooled (6744 rows)", "Pooled (6744 rows)",
                           "Pooled (6744 rows)", "Pooled (6744 rows)",
                           "Pooled (6744 rows)", "Pooled (6744 rows)",
                           "Grand mean diff", "Grand mean diff", "Grand mean diff", "Grand mean diff"),
  Correct_Recalculated_Value = c(
    sprintf("%.7f", summary_dt[Model=="M0" & Summary=="equal_weight_19folds", AP]),
    sprintf("%.7f", summary_dt[Model=="M0" & Summary=="equal_weight_19folds", LogLoss]),
    sprintf("%.7f", summary_dt[Model=="M0" & Summary=="equal_weight_19folds", ROC]),
    sprintf("%.7f", summary_dt[Model=="MP" & Summary=="equal_weight_19folds", AP]),
    sprintf("%.7f", summary_dt[Model=="MP" & Summary=="equal_weight_19folds", LogLoss]),
    sprintf("%.7f", summary_dt[Model=="MP" & Summary=="equal_weight_19folds", ROC]),
    sprintf("%.7f", summary_dt[Model=="ME" & Summary=="equal_weight_19folds", AP]),
    sprintf("%.7f", summary_dt[Model=="ME" & Summary=="equal_weight_19folds", LogLoss]),
    sprintf("%.7f", summary_dt[Model=="ME2" & Summary=="equal_weight_19folds", AP]),
    sprintf("%.7f", summary_dt[Model=="ME2" & Summary=="equal_weight_19folds", LogLoss]),
    sprintf("%.7f", summary_dt[Model=="ME2P" & Summary=="equal_weight_19folds", AP]),
    sprintf("%.7f", summary_dt[Model=="ME2P" & Summary=="equal_weight_19folds", LogLoss]),
    sprintf("%.7f", summary_dt[Model=="XEN" & Summary=="equal_weight_19folds", AP]),
    sprintf("%.7f", summary_dt[Model=="XEN" & Summary=="equal_weight_19folds", LogLoss]),
    sprintf("%.7f", calib_ref_summary[Model=="M0", mean_abs_error]),
    sprintf("%.7f", calib_ref_summary[Model=="MP", mean_abs_error]),
    sprintf("%.7f", calib_ref_summary[Model=="ME", mean_abs_error]),
    sprintf("%.7f", calib_ref_summary[Model=="ME2P", mean_abs_error])
  ),
  Difference_Reason = c(
    "08 script lacked by=ConnGroup; evaluated all 6744 rows as single pool; misattributed to AP definition differences",
    "08 script evaluated overall pooled LogLoss instead of unweighted mean of 19 fold LogLosses",
    "Per-fold M0 prediction is constant so per-fold ROC is identically 0.500; pooled ROC artifactually ranked cross-fold rates",
    "08 script evaluated pooled AP (0.6583) instead of 19-fold unweighted mean (0.7027)",
    "08 script evaluated pooled LogLoss (0.5757) instead of 19-fold unweighted mean (0.5725)",
    "08 script evaluated pooled ROC (0.7556) instead of 19-fold unweighted mean (0.7502)",
    "08 script evaluated pooled AP (0.4313) instead of 19-fold unweighted mean (0.4922)",
    "08 script evaluated pooled LogLoss (0.6843) instead of 19-fold unweighted mean (0.6609)",
    "08 script evaluated pooled AP (0.4091) instead of 19-fold unweighted mean (0.4884)",
    "08 script evaluated pooled LogLoss (0.6949) instead of 19-fold unweighted mean (0.6733)",
    "08 script evaluated pooled AP (0.6134) instead of 19-fold unweighted mean (0.7103)",
    "08 script evaluated pooled LogLoss (0.6026) instead of 19-fold unweighted mean (0.5935)",
    "08 script evaluated pooled AP (0.3787) instead of 19-fold unweighted mean (0.4864)",
    "08 script evaluated pooled LogLoss (0.7147) instead of 19-fold unweighted mean (0.6936)",
    "09 script lacked by=ReferenceKey; assigned grand mean difference across all 6744 rows rather than per-reference MAE",
    "09 script lacked by=ReferenceKey; assigned grand mean difference across all 6744 rows rather than per-reference MAE",
    "09 script lacked by=ReferenceKey; assigned grand mean difference across all 6744 rows rather than per-reference MAE",
    "09 script lacked by=ReferenceKey; assigned grand mean difference across all 6744 rows rather than per-reference MAE"
  )
)

# Save tables
fwrite(perfold_dt, file.path(out_audit, "v2_recalculated_perfold.csv"))
fwrite(summary_dt, file.path(out_audit, "v2_recalculated_summary.csv"))
fwrite(calib_groupid_dt, file.path(out_audit, "v2_recalculated_calibration_groupid.csv"))
fwrite(calib_ref_dt, file.path(out_audit, "v2_recalculated_calibration_referencekey.csv"))
fwrite(errata, file.path(out_audit, "v2_errata_comparison_table.csv"))

cat("\nErrata comparison table saved to:", file.path(out_audit, "v2_errata_comparison_table.csv"), "\n")
