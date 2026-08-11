#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

project_root <- normalizePath(".", mustWork = TRUE)
analysis_name <- "kla_ddr_weighted_score_rank_grids_33groups_v1"
table_dir <- file.path(
  project_root,
  "reanalysis/results/tables",
  analysis_name
)
figure_dir <- file.path(
  project_root,
  "reanalysis/results/figures",
  analysis_name
)

assert <- function(condition, message) {
  if (!isTRUE(condition)) {
    stop(message, call. = FALSE)
  }
}

required_tables <- file.path(
  table_dir,
  c(
    "protein_weighted_score_audit_507.csv",
    "weighted_score_plot_data_ranked_2x5.csv",
    "weighted_score_plot_data_alphabetical_2x5.csv",
    "weighted_score_summary_by_category.csv",
    "weighted_score_coefficients.csv",
    "weighted_score_category_counts.csv",
    "weighted_score_figure_manifest.csv",
    "input_file_audit.csv",
    "session_info.txt"
  )
)
figure_stems <- c(
  "kla_ddr_weighted_score_ranked_2x5_en",
  "kla_ddr_weighted_score_ranked_2x5_zh",
  "kla_ddr_weighted_score_alphabetical_2x5_en",
  "kla_ddr_weighted_score_alphabetical_2x5_zh"
)
required_figures <- file.path(
  figure_dir,
  as.vector(outer(figure_stems, c(".png", ".pdf", ".svg"), paste0))
)
required_files <- c(required_tables, required_figures)
assert(
  all(file.exists(required_files)),
  paste(
    "Missing weighted-score output(s):",
    paste(required_files[!file.exists(required_files)], collapse = "; ")
  )
)

score_audit <- fread(
  file.path(table_dir, "protein_weighted_score_audit_507.csv")
)
ranked <- fread(
  file.path(table_dir, "weighted_score_plot_data_ranked_2x5.csv")
)
alphabetical <- fread(
  file.path(table_dir, "weighted_score_plot_data_alphabetical_2x5.csv")
)
summary_table <- fread(
  file.path(table_dir, "weighted_score_summary_by_category.csv")
)
coefficients <- fread(
  file.path(table_dir, "weighted_score_coefficients.csv")
)
category_counts <- fread(
  file.path(table_dir, "weighted_score_category_counts.csv")
)
manifest <- fread(
  file.path(table_dir, "weighted_score_figure_manifest.csv")
)

weights <- c(
  BER = 1,
  NER = 2,
  MMR = 3,
  FA = 4,
  HR = 5,
  AEJ = 6,
  NHEJ = 7
)
recalculated <- as.numeric(
  as.matrix(score_audit[, names(weights), with = FALSE]) %*% weights
)
assert(
  nrow(score_audit) == 507L &&
    uniqueN(score_audit$BaseAccession) == 507L &&
    !anyNA(score_audit$SignedScore) &&
    !anyNA(score_audit$AbsoluteScore) &&
    max(abs(score_audit$SignedScore - recalculated)) < 1e-12 &&
    max(abs(score_audit$AbsoluteScore - abs(recalculated))) < 1e-12 &&
    all(score_audit$SignedScore >= -12) &&
    all(score_audit$SignedScore <= 25),
  "The 507-protein weighted-score audit is invalid."
)

assert(
  identical(
    coefficients[IncludedInScore == TRUE, Pathway],
    names(weights)
  ) &&
    max(abs(
      coefficients[IncludedInScore == TRUE, Coefficient] -
        as.numeric(weights)
    )) < 1e-12 &&
    nrow(coefficients[IncludedInScore == FALSE]) == 2L &&
    all(is.na(coefficients[IncludedInScore == FALSE, Coefficient])),
  "The saved seven-pathway coefficients or two exclusions are invalid."
)

expected_categories <- c(
  "normal_tissue",
  "normal_cells",
  "cancer_tissue",
  "cancer_cells",
  "all_507"
)
expected_counts <- c(183L, 471L, 178L, 383L, 507L)
assert(
  identical(category_counts$Category, expected_categories) &&
    identical(category_counts$ProteinCount, expected_counts) &&
    identical(summary_table$Category, expected_categories) &&
    identical(summary_table$ProteinCount, expected_counts),
  "The five plotted protein-set counts are invalid."
)

expected_plot_rows <- 2L * sum(expected_counts)
assert(
  nrow(ranked) == expected_plot_rows &&
    nrow(alphabetical) == expected_plot_rows &&
    setequal(unique(ranked$ScoreType), c("SignedScore", "AbsoluteScore")) &&
    setequal(unique(alphabetical$ScoreType), c("SignedScore", "AbsoluteScore")),
  "The 2x5 long-form plotting tables have invalid dimensions."
)

rank_check <- ranked[
  ,
  .(
    SequentialRank = identical(OrderIndex, seq_len(.N)),
    NondecreasingScore = all(diff(PlottedScore) >= 0),
    UniqueProteins = uniqueN(BaseAccession)
  ),
  by = .(ScoreType, Category)
]
assert(
  nrow(rank_check) == 10L &&
    all(rank_check$SequentialRank) &&
    all(rank_check$NondecreasingScore) &&
    identical(
      rank_check[ScoreType == "SignedScore", UniqueProteins],
      expected_counts
    ),
  "The ascending-score ordering is not valid in every panel."
)

alphabetical_check <- alphabetical[
  ,
  .(
    SequentialIndex = identical(OrderIndex, seq_len(.N)),
    NondecreasingKey = !is.unsorted(AlphabeticalKey),
    UniqueProteins = uniqueN(BaseAccession)
  ),
  by = .(ScoreType, Category)
]
assert(
  nrow(alphabetical_check) == 10L &&
    all(alphabetical_check$SequentialIndex) &&
    all(alphabetical_check$NondecreasingKey) &&
    identical(
      alphabetical_check[ScoreType == "SignedScore", UniqueProteins],
      expected_counts
    ),
  "The GeneSymbol-alphabetical comparison ordering is invalid."
)

assert(
  nrow(manifest) == 12L &&
    setequal(manifest$File, basename(required_figures)) &&
    all(file.info(required_figures)$size > 3000) &&
    all(manifest[Format == "png", PNGDPI] == 450L),
  "The figure manifest or exported graphics are invalid."
)

cat(
  paste0(
    "PASS: revised seven-pathway scores, fixed 507 proteins, ",
    "five category counts, two score types, both orderings, ",
    "and twelve figure exports are valid.\n"
  )
)
