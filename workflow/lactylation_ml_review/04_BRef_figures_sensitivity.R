#!/usr/bin/env Rscript
# 04_BRef_figures_sensitivity.R — 图 + 敏感性（删一连通组）+ 系数稳定性。
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
set.seed(25)
W <- "/Users/gzy2520/Desktop/Research/kla/.worktrees/kla31-rnaseq-20260912"
setwd(W); out <- "outputs/20260921_lactylation_ml_review"
P <- fread(file.path(out,"BRef_oos_preds.csv.gz"))
byfold <- fread(file.path(out,"BRef_perfold_perf.csv"))
mapc <- fread(file.path(out,"connectivity_recomputed.csv"))
# 折标签：成员数与PXD
connlab <- mapc[, .(members=paste(sort(GroupID),collapse="+"), pxds=paste(sort(unique(PXD)),collapse="+"), n_groups=.N), by=ConnGroup]
byfold <- merge(byfold, connlab, by="ConnGroup")

# 图1：每折 ROCAUC（M0 vs M1），附成员标注源数据
p1 <- ggplot(byfold[Model %in% c("M0","M1")], aes(x=factor(ConnGroup), y=ROCAUC, fill=Model)) +
  geom_bar(stat="identity", position="dodge") + geom_hline(yintercept=0.5, linetype="dashed") +
  labs(x="left-out connected component", y="out-of-fold ROC-AUC (within-fold)",
    title="B-R: target-gene expression vs prevalence baseline",
    subtitle="rows=(group,DDR protein)|Ref detected; y=Kla detection under that protocol; 19 leave-ConnGroup-out folds; seed=25") +
  theme_minimal(base_size=10)
ggsave(file.path(out,"fig_BRef_perfold_ROC.png"), p1, width=9, height=4.5, dpi=200)
ggsave(file.path(out,"fig_BRef_perfold_ROC.pdf"), p1, width=9, height=4.5)
fwrite(byfold, file.path(out,"fig_BRef_perfold_ROC_sourcedata.csv"))

# 图2：校准分箱（M1 pooled OOS，按预测概率十分位）
P[, bin := cut(p_M1, breaks=quantile(p_M1, probs=seq(0,1,0.1)), include.lowest=TRUE)]
cal <- P[, .(N=.N, mean_pred=mean(p_M1), obs_rate=mean(y)), by=bin]
fwrite(cal, file.path(out,"fig_BRef_calibration_M1_sourcedata.csv"))
p2 <- ggplot(cal, aes(mean_pred, obs_rate)) + geom_point(size=3) + geom_abline(slope=1, intercept=0, linetype="dashed") +
  labs(x="mean predicted P(Kla|Ref) (M1, out-of-fold)", y="observed Kla detection rate",
    title="B-R calibration (M1, pooled out-of-fold)", subtitle="10 prediction deciles; diagonal = ideal") +
  theme_minimal(base_size=11)
ggsave(file.path(out,"fig_BRef_calibration_M1.png"), p2, width=5.5, height=5, dpi=200)
ggsave(file.path(out,"fig_BRef_calibration_M1.pdf"), p2, width=5.5, height=5)

# 敏感性：删一连通组后 M1 vs M0 的 fold-mean ROCAUC 差（每次删一折重算均值）
base <- byfold[Model=="M1", mean(ROCAUC)]; base0 <- byfold[Model=="M0", mean(ROCAUC)]
sens <- rbindlist(lapply(sort(unique(byfold$ConnGroup)), function(cg) {
  s <- byfold[ConnGroup!=cg]
  data.table(dropped=cg, meanROC_M1=mean(s[Model=="M1", ROCAUC]), meanROC_M0=mean(s[Model=="M0", ROCAUC]),
    delta=mean(s[Model=="M1", ROCAUC])-mean(s[Model=="M0", ROCAUC]))
}))
fwrite(sens, file.path(out,"sensitivity_dropOneConn_M1vsM0.csv"))
cat(sprintf("full fold-mean ROC: M1=%.4f M0=%.4f delta=%.4f\n", base, base0, base-base0))
cat(sprintf("drop-one range of delta: %.4f .. %.4f\n", min(sens$delta), max(sens$delta)))
print(sens[order(delta)][c(1:.N %in% c(1:3, (.N-2):.N))])
