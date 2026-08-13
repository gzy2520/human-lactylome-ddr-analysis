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

grouping <- read_table("config", "four_class_sample_grouping.csv")
paired_scope <- merge(
  grouping,
  catalog[
    配置要求进入成对分析 %in% c(TRUE, "TRUE", "True", 1, "1"),
    .(PXD = 乳酸化PXD, SampleGroup = 样本组)
  ],
  by = c("PXD", "SampleGroup")
)
class_counts <- paired_scope[, .N, by = Category][order(Category)]
expected_class_counts <- data.table(
  Category = c(
    "cancer_cells", "cancer_tissue", "normal_cells", "normal_tissue"
  ),
  N = c(13L, 2L, 9L, 9L)
)
assert(
  identical(class_counts, expected_class_counts),
  "Expected paired four-class counts 13/2/9/9."
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
assert(nrow(reference_rows) == 30L, "Expected 30 unique reference heatmap rows.")

paired_stats <- read_table(
  "results", "tables",
  "cell_type_kla_vs_reference_ddr_statistics_accession_only_paired_33.csv"
)
assert(nrow(paired_stats) == 33L, "Expected 33 paired analysis groups.")

membership <- read_table(
  "results", "tables", "four_class_venn",
  "kla_ddr_four_class_venn", "membership.csv"
)
assert(
  nrow(membership) == 507L && uniqueN(membership$BaseAccession) == 507L,
  "Expected 507 unique Kla-DDR proteins."
)
set_counts <- c(
  normal_tissue = sum(membership$In_normal_tissue),
  normal_cells = sum(membership$In_normal_cells),
  cancer_tissue = sum(membership$In_cancer_tissue),
  cancer_cells = sum(membership$In_cancer_cells),
  all_507 = nrow(membership)
)
assert(
  identical(
    as.integer(set_counts),
    c(183L, 471L, 178L, 383L, 507L)
  ),
  "Expected five protein-set sizes 183/471/178/383/507."
)

function_inputs <- read_table(
  "results", "tables", "protein_function_inputs", "protein_metadata_and_scores.csv"
)
assignments <- read_table(
  "results", "tables", "protein_function_inputs", "pathway_assignments.csv"
)
assert(
  nrow(function_inputs) == 507L &&
    nrow(assignments) == 1175L &&
    sum(assignments$Score == 1L) == 1108L &&
    sum(assignments$Score == -1L) == 67L,
  "The fixed 507-protein pathway input contract changed."
)

bp_matrix <- read_table(
  "results", "tables", "bp_semantic_umap", "bp_semantic_binary_matrix.csv"
)
assert(
  identical(dim(bp_matrix), c(507L, 3009L)),
  "Expected a 507 x (BaseAccession + 3,008 features) BP matrix."
)

pathway_summary <- read_table(
  "results", "tables", "five_set_pathway_matrix",
  "pathway_state_summary_5sets_35rows.csv"
)
assert(nrow(pathway_summary) == 35L, "Expected 35 five-set pathway summary rows.")

venn_summary <- read_table(
  "results", "tables", "four_class_venn", "four_venn_set_counts_4x4.csv"
)
assert(
  nrow(venn_summary) == 4L &&
    identical(as.integer(venn_summary$正常组织), c(3423L, 183L, 18468L, 649L)) &&
    identical(as.integer(venn_summary$正常细胞), c(7553L, 471L, 13285L, 631L)) &&
    identical(as.integer(venn_summary$癌组织), c(2714L, 178L, 8756L, 426L)) &&
    identical(as.integer(venn_summary$癌细胞), c(4222L, 383L, 14989L, 616L)),
  "The four-Venn 4x4 set-count summary changed."
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
    "results", "figures", "bp_semantic_umap",
    "kla_ddr_raw_bp_semantic_umap_v4.png"
  ),
  file.path(
    "results", "figures", "pathway_specific_umap",
    "kla_ddr_pathway_specific_umap_3x3_v1.png"
  ),
  file.path(
    "results", "figures", "five_set_pathway_matrix",
    "kla_ddr_linear_pathway_matrix_all_507_en.png"
  )
)
missing_figures <- required_figures[!file.exists(file.path(project_root, required_figures))]
assert(
  !length(missing_figures),
  paste("Missing required figure(s):", paste(missing_figures, collapse = "; "))
)

message(
  "PASS: 40/37/33/30 group-row contract, four classes 9/9/2/13, 507 proteins, ",
  "183/471/178/383/507 sets, 3,008 BP features, and final figures."
)
