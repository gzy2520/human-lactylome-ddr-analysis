#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args)) normalizePath(args[[1L]]) else normalizePath(".")
analysis_name <- "five_set_pathway_matrix_revised_excel_20260816"
table_dir <- file.path(project_root, "results", "tables", analysis_name)
figure_dir <- file.path(project_root, "results", "figures", analysis_name)

assert <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

read_output <- function(filename) {
  path <- file.path(table_dir, filename)
  assert(file.exists(path), paste("Missing revised-score output:", path))
  fread(path)
}

scope <- read_output("score_workbook_scope_audit_507_to_399.csv")
assert(
  nrow(scope) == 507L &&
    uniqueN(scope$BaseAccession) == 507L &&
    sum(scope$InCurrent399 %in% c(TRUE, "TRUE", 1, "1")) == 399L,
  "Expected 507 scored proteins with exact coverage of the current 399."
)

set_counts <- read_output("protein_set_counts.csv")
assert(
  identical(
    as.integer(set_counts$ProteinCount),
    c(183L, 178L, 381L, 292L, 399L)
  ),
  "The revised-score 4+1 protein-set counts changed."
)

summary <- read_output("pathway_state_summary_5sets_35rows.csv")
assert(
  nrow(summary) == 35L &&
    all(
      summary$SuppressingCount +
        summary$UnassignedCount +
        summary$PromotingCount ==
        summary$ProteinCount
    ),
  "The revised-score 5 x 7 pathway summary is incomplete."
)

comparison <- read_output("revised_vs_previous_comparison_overview.csv")
assert(
  nrow(comparison) == 1L &&
    comparison$CurrentProteinCount == 399L &&
    comparison$ProteinsWithAnyStateChange == 12L &&
    comparison$ChangedProteinPathwayCells == 19L &&
    comparison$ProteinsWithScoreChange == 12L,
  "The revised-versus-previous score comparison changed."
)

figure_manifest <- read_output("figure_manifest.csv")
assert(
  nrow(figure_manifest) == 60L &&
    uniqueN(figure_manifest$Set) == 5L &&
    setequal(figure_manifest$Language, c("en", "zh")) &&
    setequal(figure_manifest$Format, c("png", "pdf", "svg")),
  "Expected 20 bilingual preview figures in three formats."
)
missing_figures <- figure_manifest[
  !file.exists(file.path(figure_dir, File)),
  File
]
assert(
  !length(missing_figures),
  paste(
    "Missing revised-score preview figure(s):",
    paste(missing_figures, collapse = "; ")
  )
)

message(
  "PASS: revised workbook has 507 unique proteins, covers the current 399, ",
  "retains 183/178/381/292/399 sets, changes 12 proteins/19 states, ",
  "and produces 20 bilingual matrix/summary figures."
)
