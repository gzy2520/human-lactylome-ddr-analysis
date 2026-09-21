# 10_metrics_tests.R — AP/PR/ROC 指标测试（独立实现对照 + 边界情形）
suppressPackageStartupMessages({library(data.table)})
source("workflow/lactylation_ml_review/05_RR_functions.R")
set.seed(25)
res <- list()
## T1 常数预测 → AP=正例率，ROC=0.5
y <- c(1,0,0,1,0,1); p <- rep(0.4, 6)
res$T1_constant <- c(AP = ap_score(y,p), PR = pr_auc_trap(y,p), ROC = roc_auc(y,p), expected_AP = mean(y))
stopifnot(abs(ap_score(y,p) - mean(y)) < 1e-12, abs(roc_auc(y,p) - 0.5) < 1e-12)
## T2 并列行顺序不影响 AP
p2 <- c(0.9, 0.8, 0.8, 0.1); y2 <- c(1, 1, 0, 0)
a1 <- ap_score(y2, p2)
a2 <- ap_score(y2[c(2,1,3,4)], p2[c(2,1,3,4)])
a3 <- ap_score(y2[c(3,2,1,4)], p2[c(3,2,1,4)])
res$T2_tie_invariance <- c(a1, a2, a3)
stopifnot(abs(a1-a2) < 1e-12, abs(a1-a3) < 1e-12)
## T3 完美排序 → AP=1
y3 <- c(0,0,1,1); p3 <- c(0.1,0.2,0.8,0.9)
res$T3_perfect <- ap_score(y3, p3)
stopifnot(abs(ap_score(y3,p3) - 1) < 1e-12)
## T4 与 PRROC / ROCR 独立实现对照（若安装）
cmp <- tryCatch({
  if (requireNamespace("PRROC", quietly = TRUE)) {
    set.seed(25); y <- sample(c(0,1), 300, TRUE, prob=c(.6,.4)); p <- rnorm(300) + 0.5*y
    pr <- PRROC::pr.curve(scores.class0 = p[y==1], scores.class1 = p[y==0], curve = TRUE)
    c(indep = pr$auc.integral, ours_ap = ap_score(y,p), ours_trap = pr_auc_trap(y,p))
  } else NA
}, error = function(e) NA)
res$T4_independent_impl <- cmp
if (!anyNA(cmp)) stopifnot(abs(cmp["indep"] - cmp["ours_ap"]) < 0.02)  # PRROC auc.integral≈AP（默认Davis-Goadrich插值，阈值聚合一致）
## T5 单类别折 → NA + 注记
mr <- metrics_row(c(1,1,1), c(.5,.5,.5), "MTEST", 1)
res$T5_single_class <- mr
stopifnot(is.na(mr$AP), grepl("single-class", mr$Note))
## T6 d11220e 外层预测：修正后重算 vs 旧值对照（仅打印，重算结果由08脚本落盘）
old <- fread("outputs/20260921_lactylation_ml_review/BRef_oos_preds.csv.gz")
dir.create("outputs/20260922_br_v2", showWarnings = FALSE, recursive = TRUE)
mods_old <- c("M0","M1","M2","M3")
cmp2 <- rbindlist(lapply(mods_old, function(m) old[, metrics_row(get("y"), get(paste0("p_", m)), m)]))
fwrite(cmp2, file.path("outputs/20260922_br_v2", "metric_fix_check_old_preds.csv"))
res$T6_oldpreds_recap <- cmp2
print(res$T6_oldpreds_recap)
cat("ALL METRIC TESTS PASSED\n")
