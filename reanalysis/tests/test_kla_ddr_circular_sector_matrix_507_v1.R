#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

project_root <- normalizePath(".", mustWork = TRUE)
analysis_name <- "kla_ddr_circular_sector_matrix_507_v1"
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
    "protein_order_and_seven_pathway_matrix_507.csv",
    "circular_sector_plot_data_3549.csv",
    "pathway_state_summary_507.csv",
    "pathway_state_density_by_25_rank_bin.csv",
    "pathway_sector_order_and_colors.csv",
    "circular_sector_figure_manifest.csv",
    "input_file_audit.csv",
    "session_info.txt"
  )
)
figure_stems <- c(
  "kla_ddr_circular_sector_matrix_507_en",
  "kla_ddr_circular_sector_matrix_507_zh",
  "kla_ddr_circular_sector_matrix_with_summary_507_en",
  "kla_ddr_circular_sector_matrix_with_summary_507_zh",
  "kla_ddr_linear_pathway_matrix_507_en",
  "kla_ddr_linear_pathway_matrix_507_zh",
  "kla_ddr_linear_pathway_matrix_with_summary_507_en",
  "kla_ddr_linear_pathway_matrix_with_summary_507_zh"
)
required_figures <- file.path(
  figure_dir,
  as.vector(outer(figure_stems, c(".png", ".pdf", ".svg"), paste0))
)
required_files <- c(required_tables, required_figures)
assert(
  all(file.exists(required_files)),
  paste(
    "Missing circular-sector output(s):",
    paste(required_files[!file.exists(required_files)], collapse = "; ")
  )
)

protein_order <- fread(required_tables[[1L]])
plot_data <- fread(required_tables[[2L]])
summary_table <- fread(required_tables[[3L]])
density <- fread(required_tables[[4L]])
sector_key <- fread(required_tables[[5L]])
manifest <- fread(required_tables[[6L]])

pathways <- c("BER", "NER", "MMR", "FA", "HR", "AEJ", "NHEJ")
weights <- 1:7
assert(
  nrow(protein_order) == 507L &&
    uniqueN(protein_order$BaseAccession) == 507L &&
    identical(protein_order$ProteinRank, seq_len(507L)) &&
    all(diff(protein_order$SignedScore) >= 0),
  "The 507-protein ascending score order is invalid."
)
recalculated <- as.numeric(
  as.matrix(protein_order[, pathways, with = FALSE]) %*% weights
)
assert(
  max(abs(protein_order$SignedScore - recalculated)) < 1e-12 &&
    identical(range(protein_order$SignedScore), c(-12L, 25L)),
  "The saved seven-pathway weighted scores are invalid."
)

assert(
  nrow(plot_data) == 507L * 7L &&
    uniqueN(plot_data[, .(BaseAccession, Pathway)]) == 507L * 7L &&
    setequal(unique(plot_data$Pathway), pathways) &&
    all(plot_data$State %in% c(-1L, 0L, 1L)),
  "The 3,549-cell circular matrix is invalid."
)

expected_counts <- plot_data[
  ,
  .(
    SuppressingCount = sum(State == -1L),
    UnassignedCount = sum(State == 0L),
    PromotingCount = sum(State == 1L)
  ),
  by = .(Pathway, SectorOrder)
][order(SectorOrder)]
assert(
  identical(summary_table$Pathway, pathways) &&
    identical(summary_table$SuppressingCount, expected_counts$SuppressingCount) &&
    identical(summary_table$UnassignedCount, expected_counts$UnassignedCount) &&
    identical(summary_table$PromotingCount, expected_counts$PromotingCount) &&
    all(
      summary_table$SuppressingCount +
        summary_table$UnassignedCount +
        summary_table$PromotingCount == 507L
    ),
  "The pathway state summary is invalid."
)

assert(
  sum(plot_data$State == -1L) == 67L,
  "The expected 67 dark-charcoal suppressing cells are not present."
)

assert(
  nrow(density) == 7L * 21L &&
    all(density$BinProteinCount %in% c(7L, 25L)) &&
    all(
      abs(
        density$SuppressingFraction +
          density$UnassignedFraction +
          density$PromotingFraction -
          1
      ) < 1e-12
    ),
  "The 25-rank-bin pathway density table is invalid."
)
assert(
  identical(sector_key$Pathway, pathways) &&
    max(abs(sector_key$Coefficient - as.numeric(weights))) < 1e-12 &&
    nrow(manifest) == 24L &&
    setequal(manifest$File, basename(required_figures)) &&
    all(file.info(required_figures)$size > 3000) &&
    all(manifest[Format == "png", PNGDPI] == 600L),
  "The sector key or figure exports are invalid."
)

cat(
  paste0(
    "PASS: 507 score-ordered proteins, 7 pathway sectors, 3,549 state cells, ",
    "67 dark-charcoal suppressing cells in circular and linear layouts, ",
    "25-rank-bin density tables, and 24 figure exports are valid.\n"
  )
)
