#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args)) normalizePath(args[[1]]) else normalizePath(".")

assert <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

read_table <- function(...) {
  path <- file.path(project_root, ...)
  assert(file.exists(path), paste("Missing publication output:", path))
  fread(path)
}

catalog <- read_table("config", "sample_group_catalog.csv")
assert(nrow(catalog) == 40L, "Expected 40 catalogued sample groups.")

excluded <- read_table("config", "main_analysis_scope_exclusions.csv")
assert(
  nrow(excluded) == 3L &&
    setequal(excluded$PXD, c("PXD014870", "PXD055230", "PXD057709")),
  "The teacher-directed three-dataset exclusion contract changed."
)

grouping <- read_table("config", "four_class_sample_grouping.csv")
paired_scope <- merge(
  grouping,
  catalog[
    配置要求进入成对分析 %in% c(TRUE, "TRUE", "True", 1, "1"),
    .(PXD = 乳酸化PXD, SampleGroup = 样本组)
  ],
  by = c("PXD", "SampleGroup")
)
paired_scope <- paired_scope[
  !paste(PXD, SampleGroup, sep = "__") %in%
    paste(excluded$PXD, excluded$SampleGroup, sep = "__")
]
class_counts <- paired_scope[, .N, by = Category][order(Category)]
expected_class_counts <- data.table(
  Category = c(
    "cancer_cells", "cancer_tissue", "normal_cells", "normal_tissue"
  ),
  N = c(12L, 2L, 7L, 9L)
)
assert(
  identical(class_counts, expected_class_counts),
  "Expected paired four-class counts 12/2/7/9."
)

group_scope <- read_table(
  "results", "tables", "kla_regulator_intensity_availability_audit.csv"
)
assert(
  sum(group_scope[["定量可用"]] %in% c(TRUE, "TRUE", "True", 1, "1")) == 37L,
  "Expected 37 quantifiable Kla groups."
)

reference_rows <- read_table(
  "results", "tables", "kla_regulator_whole_proteome_heatmap_rows.csv"
)
assert(nrow(reference_rows) == 28L, "Expected 28 unique reference heatmap rows.")

paired_stats <- read_table(
  "results", "tables",
  "cell_type_kla_vs_reference_ddr_statistics_accession_only_paired_30.csv"
)
assert(nrow(paired_stats) == 30L, "Expected 30 paired analysis groups.")

membership <- read_table(
  "results", "tables", "four_class_venn",
  "kla_ddr_four_class_venn", "membership.csv"
)
assert(
  nrow(membership) == 399L && uniqueN(membership$BaseAccession) == 399L,
  "Expected 399 unique Kla-DDR proteins."
)
set_counts <- c(
  cancer_tissue = sum(membership$In_cancer_tissue),
  cancer_cells = sum(membership$In_cancer_cells),
  normal_tissue = sum(membership$In_normal_tissue),
  normal_cells = sum(membership$In_normal_cells),
  all_kla_ddr = nrow(membership)
)
assert(
  identical(as.integer(set_counts), c(178L, 381L, 183L, 292L, 399L)),
  "Expected five protein-set sizes 178/381/183/292/399."
)

pathway_counts <- read_table(
  "results", "tables", "five_set_pathway_matrix", "protein_set_counts.csv"
)
assert(
  identical(
    pathway_counts$ProteinCount,
    c(178L, 381L, 183L, 292L, 399L)
  ),
  "The 4+1 linear pathway-matrix set counts changed."
)

pathway_summary <- read_table(
  "results", "tables", "five_set_pathway_matrix",
  "pathway_state_summary_5sets_35rows.csv"
)
assert(nrow(pathway_summary) == 35L, "Expected 35 pathway summary rows.")

venn_summary <- read_table(
  "results", "tables", "four_class_venn", "four_venn_set_counts_4x4.csv"
)
assert(
  nrow(venn_summary) == 4L &&
    identical(as.integer(venn_summary$肿瘤组织), c(2714L, 178L, 8756L, 426L)) &&
    identical(as.integer(venn_summary$癌细胞系), c(4183L, 381L, 14989L, 616L)) &&
    identical(as.integer(venn_summary$非肿瘤组织), c(3423L, 183L, 18468L, 649L)) &&
    identical(as.integer(venn_summary$正常细胞系), c(3234L, 292L, 11678L, 592L)),
  "The current four-Venn 4x4 set-count summary changed."
)

required_figures <- c(
  file.path(
    "results", "figures", "kla_regulator_cross_study_relative_intensity_heatmap_en.png"
  ),
  file.path(
    "results", "figures", "kla_regulator_whole_proteome_relative_intensity_heatmap_en.png"
  ),
  file.path(
    "results", "figures", "cell_type_kla_vs_reference_ddr_fraction_en.png"
  ),
  file.path(
    "results", "figures", "four_class_venn",
    "kla_ddr_four_class_venn_30groups_en.png"
  ),
  file.path(
    "results", "figures", "five_set_pathway_matrix",
    "kla_ddr_linear_pathway_matrix_all_kla_ddr_en.png"
  )
)
missing_figures <- required_figures[!file.exists(file.path(project_root, required_figures))]
assert(
  !length(missing_figures),
  paste("Missing required figure(s):", paste(missing_figures, collapse = "; "))
)

message(
  "PASS: 40/37/30/28 scope, four classes 2/9/12/7 in display order, 399 Kla-DDR proteins, ",
  "178/381/183/292/399 sets, and selected publication figures."
)
