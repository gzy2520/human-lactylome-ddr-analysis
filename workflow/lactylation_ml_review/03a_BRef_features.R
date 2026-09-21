#!/usr/bin/env Rscript
# 03_train_BRef.R — 候选 B-R 外层留连通分组预测（19折），折内预处理，无泄漏。
# 行=(GroupID,蛋白)|Ref检出；y=Kla检出（流程下检出结果）。
# X（逐行）: ①目标基因表达 ②乳酸128均值 ③调控因子均值，均取自 UNSMOOTHED log2TPM。
# 折内步骤（仅训练折ReferenceKey）: 材料中位数聚合 -> 逐基因中心/尺度（训练估计，应用验证）。
# 模型（同一外层折比较）: M0 训练正例率基线；M1 目标基因单变量逻辑回归；
# M2 M1+乳酸均值；M3 M2+调控均值。glm(base)；若某折分离则回退M0并记录。
# 指标（外层留出 pooled + 按折）：PR-AUC（自算平均精度）、ROC-AUC（ROCR若无则Mann-Whitney）、Brier。
suppressPackageStartupMessages({library(data.table)})
set.seed(25)
W <- "/Users/gzy2520/Desktop/Research/kla/.worktrees/kla31-rnaseq-20260912"
setwd(W); out <- "outputs/20260921_lactylation_ml_review"
task <- fread(file.path(out,"BRef_task_table.csv.gz"))
task[, y := as.integer(y)]
# 基因集
lactab <- fread("outputs/20260920_lactate_metabolism/tables/membership_coverage.csv")
lac128 <- unique(lactab[InCommonMatrix==TRUE, Ensembl])
reg <- fread("outputs/20260919_repair_release/core/regulator_rna_qsmooth_28materials.csv")
regG <- unique(reg$EnsemblGeneID)
# 材料中位数聚合函数（逐ReferenceKey；读源矩阵）
expr_dir <- "outputs/20260919_expression_corrected/matrices"
med_of_ref <- function(rk, genes) {
  f <- file.path(expr_dir, paste0(rk, "_log2tpm.tsv.gz"))
  d <- fread(f, sep="\t", header=TRUE, check.names=FALSE)
  ids <- d[[1L]]; m <- as.matrix(d[, -1L, with=FALSE]); rownames(m) <- ids
  g <- intersect(genes, rownames(m))
  v <- apply(m[g, , drop=FALSE], 1, median, na.rm=TRUE)
  outv <- rep(NA_real_, length(genes)); names(outv) <- genes
  outv[g] <- v
  outv
}
need_genes <- unique(c(unique(task$EnsemblGeneID), lac128, regG))
cat(sprintf("need genes=%d (target %d + lac128 %d + reg %d)\n", length(need_genes), uniqueN(task$EnsemblGeneID), length(lac128), length(regG)))
refs <- unique(task$ReferenceKey)
mat_med <- sapply(refs, function(rk) med_of_ref(rk, need_genes))
colnames(mat_med) <- refs
cat(sprintf("material medians: %d genes x %d refs\n", nrow(mat_med), ncol(mat_med)))
# 每行特征（全量，供折内标准化用；注意：中心/尺度一律折内估计）
task[, target_expr := mat_med[cbind(match(EnsemblGeneID, rownames(mat_med)), match(ReferenceKey, colnames(mat_med)))]]
lac_med <- colMeans(mat_med[lac128, , drop=FALSE], na.rm=TRUE)
reg_med <- colMeans(mat_med[regG, , drop=FALSE], na.rm=TRUE)
task[, lac_mean := lac_med[ReferenceKey]]
task[, reg_mean := reg_med[ReferenceKey]]
fwrite(task[, .(GroupID, ReferenceKey, ConnGroup, PXD, BaseAccession, EnsemblGeneID, y, target_expr, lac_mean, reg_mean)],
  file.path(out,"BRef_features.csv.gz"))
cat("FEATURES SAVED; NA target_expr rows:", sum(is.na(task$target_expr)), "\n")
