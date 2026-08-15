#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args)) normalizePath(args[[1]]) else normalizePath(".")
table_root <- file.path(project_root, "results", "tables", "four_class_venn")
output_path <- file.path(table_root, "four_venn_set_counts_4x4.csv")

analyses <- data.table(
  Analysis = c(
    "all_kla_four_class_venn",
    "kla_ddr_four_class_venn",
    "reference_proteome_four_class_venn",
    "reference_proteome_ddr_four_class_venn"
  ),
  AnalysisZh = c("全部Kla蛋白", "Kla∩DDR蛋白", "普通全蛋白", "普通全蛋白∩DDR蛋白")
)
categories <- c(
  cancer_tissue = "肿瘤组织",
  normal_tissue = "非肿瘤组织",
  cancer_cells = "癌细胞系",
  normal_cells = "正常细胞系"
)

rows <- lapply(seq_len(nrow(analyses)), function(i) {
  path <- file.path(table_root, analyses$Analysis[i], "set_counts.csv")
  if (!file.exists(path)) stop("Missing set-count table: ", path)
  counts <- fread(path)
  if (!all(names(categories) %in% counts$Category)) stop("Missing category in: ", path)
  out <- data.table(Analysis = analyses$Analysis[i], AnalysisZh = analyses$AnalysisZh[i])
  for (category in names(categories)) {
    out[[categories[[category]]]] <- counts[Category == category, ProteinCount]
  }
  out
})
summary <- rbindlist(rows)
dir.create(table_root, recursive = TRUE, showWarnings = FALSE)
fwrite(summary, output_path, bom = TRUE)
message("Wrote ", output_path)
