# 06_build_BR_v2_dataset.R — B-R v2 数据集与契约核查
# 核查6744/2708一致 → 剔除目标基因无矩阵值蛋白行 → 保存 v2 分析集合（全模型共用）。
suppressPackageStartupMessages({library(data.table)})
set.seed(25)
W <- "/Users/gzy2520/Desktop/Research/kla/.worktrees/kla31-rnaseq-20260912"
setwd(W); out <- "outputs/20260922_br_v2"; dir.create(out, showWarnings = FALSE, recursive = TRUE)
task <- fread("outputs/20260921_lactylation_ml_review/BRef_task_table.csv.gz")
cat(sprintf("contract check: rows=%d positives=%d (expect 6744/2708) groups=%d refs=%d proteins=%d\n",
  nrow(task), sum(task$y), uniqueN(task$GroupID), uniqueN(task$ReferenceKey), uniqueN(task$BaseAccession)))
stopifnot(nrow(task) == 6744, sum(task$y) == 2708)
## 目标基因可用性（各参考矩阵）
refs <- unique(task$ReferenceKey)
avail <- sapply(refs, function(rk) {
  f <- file.path("outputs/20260919_expression_corrected/matrices", paste0(rk, "_log2tpm.tsv.gz"))
  fread(f, select = 1)[[1]]
})
target_genes <- unique(task$EnsemblGeneID)
## 全参考均缺失的目标基因 → 剔除这些蛋白的所有行
miss <- target_genes[vapply(target_genes, function(g) all(!vapply(avail, function(a) g %in% a, logical(1))), logical(1))]
cat("target genes missing from ALL matrices:", length(miss), "\n"); print(miss)
n_before <- nrow(task)
task2 <- task[!(EnsemblGeneID %in% miss)]
cat(sprintf("v2 analysis set: rows=%d (dropped %d), positives=%d\n", nrow(task2), n_before - nrow(task2), sum(task2$y)))
fwrite(task2, file.path(out, "BRef_task_v2.csv.gz"))
fwrite(data.table(ReferenceKey = refs,
  n_targets_avail = vapply(refs, function(rk) sum(target_genes %in% avail[[rk]]), integer(1))),
  file.path(out, "target_coverage_by_ref.csv"))
cat("V2 DATASET SAVED\n")
