# 11_figures_v2.R — 校准图 + 逐折AP对比图 + 新旧指标对照图（渲染检查）
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
source("workflow/lactylation_ml_review/05_RR_functions.R")
set.seed(25)
W <- "/Users/gzy2520/Desktop/Research/kla/.worktrees/kla31-rnaseq-20260912"
setwd(W); out <- "outputs/20260922_br_v2"
P <- fread(file.path(out, "BRv2_oos_preds.csv.gz"))
## F1: 逐折 AP（MP vs ME vs M0）
per <- rbindlist(lapply(c("M0","MP","ME"), function(m) {
  P[, .(Model=m, ConnGroup=factor(ConnGroup), AP=ap_score(y, get(paste0("p_",m)))), by=ConnGroup][, .(Model, ConnGroup, AP)]
}))
per <- per[, .(Model, ConnGroup, AP)]
p1 <- ggplot(per, aes(ConnGroup, AP, fill=Model)) + geom_col(position="dodge") +
  geom_hline(yintercept=0.402, linetype="dashed") +
  labs(x="left-out connected component", y="out-of-fold AP (threshold-aggregated)",
    title="B-R v2: protein-frequency baseline (MP) vs target expression (ME) vs prevalence (M0)",
    subtitle="dashed = pooled prevalence 0.402; 19 leave-connected-component-out folds; seed 25") +
  theme_minimal(base_size=10)
ggsave(file.path(out, "fig_v2_perfold_AP.png"), p1, width=9, height=4.5, dpi=200)
ggsave(file.path(out, "fig_v2_perfold_AP.pdf"), p1, width=9, height=4.5)
fwrite(per, file.path(out, "fig_v2_perfold_AP_sourcedata.csv"))
## F2: 材料级校准（MP、ME2P）
matcal <- fread(file.path(out, "BRv2_material_level_calibration.csv"))
mc <- matcal[Model %in% c("MP","ME2P")]
p2 <- ggplot(mc, aes(mean_pred, obs_rate, color=Model)) + geom_point(size=3) +
  geom_abline(slope=1, intercept=0, linetype="dashed") +
  labs(x="mean predicted P(Kla detected | Ref)", y="observed detection rate",
    title="B-R v2 material-level calibration (out-of-fold)", subtitle="each point = one RNA reference material; diagonal = ideal") +
  theme_minimal(base_size=11)
ggsave(file.path(out, "fig_v2_material_calibration.png"), p2, width=6, height=5.5, dpi=200)
ggsave(file.path(out, "fig_v2_material_calibration.pdf"), p2, width=6, height=5.5)
## F3: 指标修正前后对照（d11220e pooled）
old <- fread("outputs/20260921_lactylation_ml_review/BRef_pooled_perf.csv")
new <- fread(file.path(out, "metric_fix_check_old_preds.csv"))
cmp <- rbind(old[, .(Model, Metric="old_rowwise_PR_AUC", Value=PRAUC)],
  new[, .(Model, Metric="new_threshold_AP", Value=AP)],
  new[, .(Model, Metric="new_trap_PR_AUC", Value=PRAUC_trap)])
fwrite(cmp, file.path(out, "metric_fix_comparison_d11220e.csv"))
cat("figures written\n")
## 渲染检查
for (f in c("fig_v2_perfold_AP.png","fig_v2_material_calibration.png")) {
  stopifnot(file.exists(file.path(out, f)), file.size(file.path(out, f)) > 10000)
}
cat("render check OK\n")
