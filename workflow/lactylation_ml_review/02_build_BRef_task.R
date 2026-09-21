#!/usr/bin/env Rscript
# 02_build_BRef_task.R — 构建候选 B-R（蛋白×材料，Ref限定DDR panel）任务表。
# 行=(GroupID, DDR蛋白) 且 RefDetected==TRUE；y=KlaDetected（该流程下检出结果，非“发生乳酸化”真值）。
# X=逐基因 UNSMOOTHED log2TPM（训练折内聚合；此处先存全量，折内聚合在03执行）。
# 特征集（冻结稳定ID）：①目标基因表达 ②乳酸模块均值（128基因，固定成员）③调控因子均值（固定panel）。
suppressPackageStartupMessages({library(data.table)})
set.seed(25)
W <- "/Users/gzy2520/Desktop/Research/kla/.worktrees/kla31-rnaseq-20260912"
setwd(W); out <- "outputs/20260921_lactylation_ml_review"

# --- 标签侧（已验收产物，不重算蛋白组） ---
p <- fread("outputs/20260919_repair_assisted/group_by_protein_pairs.tsv.gz")
pr <- p[RefDetected==TRUE, .(GroupID, BaseAccession, EnsemblGeneID, y=KlaDetected)]
cat(sprintf("B-R rows=%d positives=%d rate=%.3f groups=%d proteins=%d\n", nrow(pr), sum(pr$y), mean(pr$y), uniqueN(pr$GroupID), uniqueN(pr$BaseAccession)))

# --- 特征字典（冻结） ---
mod <- fread("config/lactate_metabolism/gene_module_membership.csv")
# 128共同矩阵基因（去重Ensembl；LDHA双ID等保留）
common_genes <- fread("outputs/20260920_lactate_metabolism/tables/membership_coverage.csv")
common128 <- unique(common_genes[InCommonMatrix==TRUE, Ensembl])
cat(sprintf("lactate common genes=%d\n", length(common128)))
# 调控因子panel（regulator RNA表，48 regulators/49 entries → 取去重Ensembl）
reg <- fread("outputs/20260919_repair_release/core/regulator_rna_qsmooth_28materials.csv")
reg_genes <- unique(reg$EnsemblGeneID)
cat(sprintf("regulator genes=%d\n", length(reg_genes)))
# 目标基因：DDR panel 371的Ensembl（来自任务表本身）
target_genes <- unique(pr$EnsemblGeneID)
# 三集合写入冻结字典
dict <- data.table(FeatureSet=c("target_gene_expr","lactate_module_mean_128","regulator_mean"),
  N=c(length(target_genes), length(common128), length(reg_genes)),
  Source=c("outputs/20260919_repair_assisted/group_by_protein_pairs.tsv.gz:EnsemblGeneID (371 DDR, 1:1 bridge)",
    "config/lactate_metabolism/gene_module_membership.csv + membership_coverage.csv InCommonMatrix (128)",
    "outputs/20260919_repair_release/core/regulator_rna_qsmooth_28materials.csv (de-duplicated Ensembl)"),
  IDKey="Ensembl (RNA); protein key UniProt BaseAccession; Symbol display-only")
fwrite(dict, file.path(out,"feature_sets_frozen.csv"))

# --- 表达侧：读源log2tpm矩阵（UNSMOOTHED），按ReferenceKey聚合材料代表值 ---
# 为折内预处理准备：此处保存“每ReferenceKey×基因”的样本级长表太重，故保存材料中位数表（全量），
# 03训练时按“仅训练折ReferenceKey”重聚合（函数化，不用全量统计量）。
idx <- fread("outputs/20260919_expression_corrected/group_index.csv")  # GroupID/ReferenceKey
# GroupID->ReferenceKey（含KLA31_15/16共享HCT116；KLA31_28共享HK2 — 以expansion为准）
exp31 <- fread("outputs/20260919_repair_qsmooth/group_expansion_31.csv")
cat("expansion rows:", nrow(exp31), "\n"); print(exp31[, .(GroupID, ReferenceKey)])
# 保存任务表（含GroupID/ReferenceKey/ConnGroup）
mapc <- fread(file.path(out,"connectivity_recomputed.csv"))
pr <- merge(pr, mapc[, .(GroupID, ReferenceKey, ConnGroup, PXD)], by="GroupID")
fwrite(pr, file.path(out,"BRef_task_table.csv.gz"))
cat("TASK TABLE SAVED\n")
# 基因集交集检查：target/lactate/regulator 在各源矩阵的可得性（抽查3个矩阵）
for (rk in c("HCT116_control","GSE132714_BPH","TCGA_PRAD_primary")) {
  f <- file.path("outputs/20260919_expression_corrected/matrices", paste0(rk,"_log2tpm.tsv.gz"))
  ids <- fread(f, select=1)[[1]]
  cat(sprintf("%s: target %d/%d; lactate128 %d/128; regul %d/%d\n", rk,
    sum(target_genes %in% ids), length(target_genes),
    sum(common128 %in% ids), sum(reg_genes %in% ids), length(reg_genes)))
}
