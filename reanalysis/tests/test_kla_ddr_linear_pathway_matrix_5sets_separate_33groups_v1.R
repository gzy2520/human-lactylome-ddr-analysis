#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

project_root <- normalizePath(".", mustWork = TRUE)
analysis_name <- "kla_ddr_linear_pathway_matrix_5sets_separate_33groups_v1"
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
    "protein_order_and_seven_pathway_matrix_5sets.csv",
    "linear_matrix_plot_data_12054.csv",
    "pathway_state_summary_5sets_35rows.csv",
    "protein_set_counts.csv",
    "pathway_order_and_colors.csv",
    "figure_manifest.csv",
    "input_file_audit.csv",
    "session_info.txt"
  )
)
sets <- c(
  "normal_tissue",
  "normal_cells",
  "cancer_tissue",
  "cancer_cells",
  "all_507"
)
figure_stems <- unlist(
  lapply(sets, function(set_key) {
    c(
      paste0("kla_ddr_linear_pathway_matrix_", set_key, "_en"),
      paste0("kla_ddr_linear_pathway_matrix_", set_key, "_zh"),
      paste0("kla_ddr_pathway_state_summary_", set_key, "_en"),
      paste0("kla_ddr_pathway_state_summary_", set_key, "_zh")
    )
  }),
  use.names = FALSE
)
required_figures <- file.path(
  figure_dir,
  as.vector(outer(figure_stems, c(".png", ".pdf", ".svg"), paste0))
)
required_files <- c(required_tables, required_figures)
assert(
  all(file.exists(required_files)),
  paste(
    "Missing 4+1 linear output(s):",
    paste(required_files[!file.exists(required_files)], collapse = "; ")
  )
)

protein_order <- fread(required_tables[[1L]])
plot_data <- fread(required_tables[[2L]])
summary_table <- fread(required_tables[[3L]])
set_counts <- fread(required_tables[[4L]])
pathway_key <- fread(required_tables[[5L]])
manifest <- fread(required_tables[[6L]])

expected_counts <- c(183L, 471L, 178L, 383L, 507L)
pathways <- c("BER", "NER", "MMR", "FA", "HR", "AEJ", "NHEJ")
weights <- 1:7

assert(
  identical(set_counts$Set, sets) &&
    identical(set_counts$ProteinCount, expected_counts) &&
    nrow(protein_order) == sum(expected_counts),
  "The five protein sets or counts are invalid."
)

rank_check <- protein_order[
  ,
  .(
    SequentialRank = identical(ProteinRank, seq_len(.N)),
    NondecreasingScore = all(diff(SignedScore) >= 0),
    UniqueProteins = uniqueN(BaseAccession)
  ),
  by = .(Set, SetOrder)
][order(SetOrder)]
assert(
  identical(rank_check$Set, sets) &&
    all(rank_check$SequentialRank) &&
    all(rank_check$NondecreasingScore) &&
    identical(rank_check$UniqueProteins, expected_counts),
  "Independent within-set score ordering is invalid."
)

recalculated <- as.numeric(
  as.matrix(protein_order[, pathways, with = FALSE]) %*% weights
)
assert(
  max(abs(protein_order$SignedScore - recalculated)) < 1e-12,
  "The seven-pathway weighted scores are invalid."
)

assert(
  nrow(plot_data) == 7L * sum(expected_counts) &&
    uniqueN(plot_data[, .(Set, BaseAccession, Pathway)]) ==
      7L * sum(expected_counts) &&
    setequal(unique(plot_data$Pathway), pathways) &&
    all(plot_data$State %in% c(-1L, 0L, 1L)),
  "The 12,054-cell 4+1 linear plotting matrix is invalid."
)

expected_summary <- plot_data[
  ,
  .(
    SuppressingCount = sum(State == -1L),
    UnassignedCount = sum(State == 0L),
    PromotingCount = sum(State == 1L),
    ProteinCount = uniqueN(BaseAccession)
  ),
  by = .(Set, SetOrder, Pathway, PathwayOrder)
][order(SetOrder, PathwayOrder)]
assert(
  nrow(summary_table) == 35L &&
    identical(summary_table$Set, expected_summary$Set) &&
    identical(summary_table$Pathway, expected_summary$Pathway) &&
    identical(summary_table$SuppressingCount, expected_summary$SuppressingCount) &&
    identical(summary_table$UnassignedCount, expected_summary$UnassignedCount) &&
    identical(summary_table$PromotingCount, expected_summary$PromotingCount) &&
    max(abs(
      summary_table$SuppressingFraction +
        summary_table$UnassignedFraction +
        summary_table$PromotingFraction -
        1
    )) < 1e-12,
  "The 35-row pathway-state summary is invalid."
)

assert(
  identical(pathway_key$Pathway, pathways) &&
    max(abs(pathway_key$Coefficient - as.numeric(weights))) < 1e-12 &&
    nrow(manifest) == 60L &&
    all(manifest[, .N, by = .(FigureType, Set, Language)]$N == 3L) &&
    setequal(manifest$File, basename(required_figures)) &&
    all(file.info(required_figures)$size > 3000) &&
    all(manifest[Format == "png", PNGDPI] == 600L),
  "The pathway key or 60 fully separate figure exports are invalid."
)

cat(
  paste0(
    "PASS: five independently score-ranked protein sets, 12,054 pathway cells, ",
    "35 pathway summaries, final +1/−1/0 colors, and 60 fully separate matrix/summary ",
    "exports are valid.\n"
  )
)
