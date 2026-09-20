#!/usr/bin/env Rscript
# 02_gate_and_simulation.R — training admission gate + isolated simulation code test.
# REAL DATA: only audit/feasibility (no training). SIMULATION: synthetic X/y to test
# grouped-CV plumbing; simulation outputs are quarantined and NEVER reported as results.
suppressPackageStartupMessages({library(data.table)})
set.seed(25)
root <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE))[1L]), "..", ".."))
setwd(root)
out_dir <- "outputs/20260920_lactylation_ml"
sim_dir <- file.path(out_dir, "simulation_quarantine_DO_NOT_REPORT")
dir.create(sim_dir, showWarnings=FALSE, recursive=TRUE)

obs <- fread("outputs/20260920_lactylation_ml/observation_contract.csv")
# ---- Gate G1..G7 (all must pass; any FAIL blocks training) ----
gate <- data.table(
  Gate=c("G1_label_comparable", "G2_independent_units", "G3_fold_label_variation",
    "G4_no_unmanaged_sharing", "G5_feature_availability", "G6_confounding_not_deterministic", "G7_qsmooth_suitability"),
  Result=c(
    "FAIL: cross-study Kla counts / detection proportions use incompatible protocols+instruments+search spaces; no calibrated intensity/occupancy label",
    "FAIL: 28 RNA refs / 31 proteome rows collapse to 19 leakage-connected components (shared RNA refs + shared PXDs); only 21 rows have any same-sample matched modality spanning just 3 PXDs; no paired RNA-protein design",
    "FAIL: per-fold label distribution for leave-study-out cannot be guaranteed with n=1 studies and single-observation groups",
    "FAIL: HCT116_control shared x3 + HK2_untreated shared x2; PXD066054/PXD075377/PXD028488/PXD060185/PXD070007 each span multiple materials; study==material confounding",
    "FAIL: lactate-axis expression matrices absent locally (server-only); measured lactate absent; same-sample total-protein absent for most groups",
    "FAIL: detection depth (ReferenceProteinCount) spans 2193-14022; material category confounds study",
    "FAIL: qsmooth fitted jointly across all 28 refs incl. would-be test refs; cannot be refit train-only from local files"),
  BlocksTraining=c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE))
fwrite(gate, file.path(out_dir, "training_gate.csv"))
cat("GATE: all 7 gates FAIL -> FORMAL TRAINING BLOCKED (see training_gate.csv)\n")

# ---- Isolated simulation: test grouped-CV code path only ----
# Synthetic: 19 connected groups with matching size profile to real connectivity
# (1,1,2,1,1,1,2,2,1,4,5,1,1,1,2,1,1,2,1), 5 fake features.
set.seed(25)
sizes <- c(1L, 1L, 2L, 1L, 1L, 1L, 2L, 2L, 1L, 4L, 5L, 1L, 1L, 1L, 2L, 1L, 1L, 2L, 1L)
grp <- rep(sprintf("SIMG_%02d", seq_along(sizes)), times=sizes)
n <- length(grp)
X <- matrix(rnorm(n * 5), n, 5, dimnames=list(NULL, paste0("F", 1:5)))
y <- 0.8 * X[, 1] - 0.5 * X[, 2] + rnorm(n)
# Leave-one-group-out ridge (lambda fixed a priori; closed form) + mean baseline, strictly group-aware.
ridge_log <- function() {}
lambda <- 1.0
pred_ridge <- numeric(n); pred_mean <- numeric(n)
for (g in unique(grp)) {
  tr <- grp != g; te <- grp == g
  Xtr <- cbind(1, X[tr, , drop=FALSE]); Xte <- cbind(1, X[te, , drop=FALSE]); ytr <- y[tr]
  P <- diag(ncol(Xtr)); P[1, 1] <- 0
  beta <- solve(crossprod(Xtr) + lambda * P, crossprod(Xtr, ytr))
  pred_ridge[te] <- as.numeric(Xte %*% beta)
  pred_mean[te] <- mean(ytr)
}
mae <- function(a, b) mean(abs(a - b)); rmse <- function(a, b) sqrt(mean((a - b)^2))
r2_oos <- function(y, p) 1 - sum((y - p)^2) / sum((y - mean(y))^2)
sim_perf <- data.table(Model=c("SIM_baseline_mean", "SIM_ridge_lambda1"),
  MAE=c(mae(y, pred_mean), mae(y, pred_ridge)),
  RMSE=c(rmse(y, pred_mean), rmse(y, pred_ridge)),
  R2_oos=c(r2_oos(y, pred_mean), r2_oos(y, pred_ridge)))
fwrite(sim_perf, file.path(sim_dir, "SIMULATION_perf_DO_NOT_REPORT.csv"))
fwrite(data.table(Observation=seq_len(n), Group=grp, y=y, pred_mean=pred_mean, pred_ridge=pred_ridge),
  file.path(sim_dir, "SIMULATION_oos_preds_DO_NOT_REPORT.csv"))
writeLines(c("SIMULATION QUARANTINE: synthetic data only. DO NOT REPORT as study results.",
  "Purpose: verify grouped leave-one-group-out plumbing runs without leakage across groups.",
  sprintf("seed=25 lambda=%.1f groups=%d n=%d (sizes match real 19-component profile)", lambda, length(sizes), n)),
  file.path(sim_dir, "README_DO_NOT_REPORT.txt"))
cat("SIMULATION plumbing OK (quarantined).\n")
print(sim_perf)
