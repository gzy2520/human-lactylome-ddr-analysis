# 07_build_features_v2.R — 材料级基因中位数（仅材料内聚合，非跨材料统计）+ 冻结模块均值。
suppressPackageStartupMessages({library(data.table)})
source("workflow/lactylation_ml_review/05_RR_functions.R")
set.seed(25)
W <- "/Users/gzy2520/Desktop/Research/kla/.worktrees/kla31-rnaseq-20260912"
setwd(W); out <- "outputs/20260922_br_v2"
task <- fread(file.path(out, "BRef_task_v2.csv.gz"))
mods <- fread("config/lactylation_ml_review/modules_frozen_v2.csv")
refs <- unique(task$ReferenceKey)
need <- unique(c(task$EnsemblGeneID, mods$Ensembl))
expr_dir <- "outputs/20260919_expression_corrected/matrices"
mat_med <- sapply(refs, function(rk) {
  d <- fread(file.path(expr_dir, paste0(rk, "_log2tpm.tsv.gz")))
  m <- as.matrix(d[, -1, with = FALSE]); rownames(m) <- d[[1]]
  g <- intersect(need, rownames(m))
  v <- apply(m[g, , drop = FALSE], 1, median, na.rm = TRUE)
  v[setdiff(need, g)] <- NA_real_
  v
})
colnames(mat_med) <- refs
cat(sprintf("mat_med: %d genes x %d refs\n", nrow(mat_med), ncol(mat_med)))
## 目标基因特征（行级）
task[, target_expr := mat_med[cbind(match(EnsemblGeneID, rownames(mat_med)), match(ReferenceKey, colnames(mat_med)))]]
cat("NA target rows:", sum(is.na(task$target_expr)), "\n")
## 模块均值（材料级）
mod_names <- unique(mods$ModuleID)
for (m in mod_names) {
  v <- module_mean_per_ref(mat_med, mods[ModuleID == m, Ensembl])
  task[, (m) := v[ReferenceKey]]
  cov <- mean(!is.na(v[refs]))
  cat(sprintf("module %s: %d members, coverage=%.2f, range=[%.2f, %.2f]\n", m,
    nrow(mods[ModuleID == m]), cov, suppressWarnings(min(v, na.rm = TRUE)), suppressWarnings(max(v, na.rm = TRUE))))
}
fwrite(task, file.path(out, "BRef_features_v2.csv.gz"))
fwrite(data.table(ModuleID = mod_names,
  Members = sapply(mod_names, function(m) nrow(mods[ModuleID == m])),
  AvailInMatrices = sapply(mod_names, function(m) sum(unique(mods[ModuleID == m, Ensembl]) %in% rownames(mat_med)))),
  file.path(out, "module_coverage_v2.csv"))
cat("V2 FEATURES SAVED\n")
