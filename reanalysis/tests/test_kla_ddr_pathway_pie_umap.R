#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
script_path <- if (length(file_arg) > 0L) {
  normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE)
} else {
  normalizePath("reanalysis/tests/test_kla_ddr_pathway_pie_umap.R", mustWork = TRUE)
}
project_root <- normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = TRUE)
table_dir <- file.path(
  project_root,
  "reanalysis/results/tables/kla_ddr_pathway_pie_umap_33groups"
)
figure_dir <- file.path(
  project_root,
  "reanalysis/results/figures/kla_ddr_pathway_pie_umap_33groups"
)

required_outputs <- c(
  file.path(table_dir, "pathway_scores_507.csv"),
  file.path(table_dir, "pathway_assignments_with_evidence_1175.csv"),
  file.path(table_dir, "pathway_assignment_summary.csv"),
  file.path(table_dir, "pathway_umap_plot_data_fixed_coordinates.csv"),
  file.path(table_dir, "pathway_umap_coverage_audit.csv"),
  file.path(table_dir, "input_file_audit.csv"),
  file.path(table_dir, "pathway_color_key.csv"),
  file.path(table_dir, "session_info.txt"),
  file.path(figure_dir, "kla_ddr_pathway_signed_pie_umap_33groups.png"),
  file.path(figure_dir, "kla_ddr_pathway_signed_pie_umap_33groups.pdf"),
  file.path(figure_dir, "kla_ddr_pathway_signed_pie_umap_33groups.svg"),
  file.path(figure_dir, "kla_ddr_pathway_umap_and_signed_summary_33groups.png"),
  file.path(figure_dir, "kla_ddr_pathway_umap_and_signed_summary_33groups.pdf"),
  file.path(figure_dir, "kla_ddr_pathway_umap_and_signed_summary_33groups.svg"),
  file.path(project_root, "reanalysis/reports/UMAP_PATHWAY_PIE_33GROUP_DATA_SCOPE.md")
)
stopifnot(all(file.exists(required_outputs)))
stopifnot(all(file.info(required_outputs)$size > 0L))

coordinates <- fread(
  file.path(
    project_root,
    "reanalysis/results/tables/kla_ddr_raw_go_umap_33groups/umap_coordinates_fixed.csv"
  )
)
scores <- fread(file.path(table_dir, "pathway_scores_507.csv"))
assignments <- fread(
  file.path(table_dir, "pathway_assignments_with_evidence_1175.csv")
)
summary_counts <- fread(file.path(table_dir, "pathway_assignment_summary.csv"))
plot_data <- fread(file.path(table_dir, "pathway_umap_plot_data_fixed_coordinates.csv"))
coverage <- fread(file.path(table_dir, "pathway_umap_coverage_audit.csv"))

stopifnot(nrow(scores) == 507L)
stopifnot(uniqueN(scores$BaseAccession) == 507L)
stopifnot(setequal(scores$BaseAccession, coordinates$BaseAccession))

stopifnot(nrow(plot_data) == 507L)
stopifnot(uniqueN(plot_data$BaseAccession) == 507L)
coordinate_check <- merge(
  coordinates,
  plot_data[, .(BaseAccession, PlotUMAP1 = UMAP_1, PlotUMAP2 = UMAP_2)],
  by = "BaseAccession"
)
stopifnot(all(coordinate_check$UMAP_1 == coordinate_check$PlotUMAP1))
stopifnot(all(coordinate_check$UMAP_2 == coordinate_check$PlotUMAP2))

stopifnot(nrow(assignments) == 1175L)
stopifnot(uniqueN(assignments, by = c("ID", "Pathway", "Score")) == 1175L)
stopifnot(all(assignments$Score %in% c(-1L, 1L)))
stopifnot(!anyNA(assignments$DOI_or_SourceURL))
stopifnot(sum(assignments$Score == 1L) == 1108L)
stopifnot(sum(assignments$Score == -1L) == 67L)

stopifnot(nrow(summary_counts) == 9L)
stopifnot(sum(summary_counts$Positive) == 1108L)
stopifnot(sum(summary_counts$Negative) == 67L)
stopifnot(sum(summary_counts$TotalNonzero) == 1175L)

stopifnot(plot_data[TotalAssignmentCount > 0L, .N] == 485L)
stopifnot(plot_data[TotalAssignmentCount == 0L, .N] == 22L)
stopifnot(plot_data[NegativeAssignmentCount > 0L, .N] == 53L)
stopifnot(sum(plot_data$NegativeAssignmentCount) == 67L)

coverage_values <- setNames(coverage$Value, coverage$Item)
stopifnot(coverage_values[["Fixed UMAP proteins"]] == 507L)
stopifnot(coverage_values[["Exact BaseAccession matches"]] == 507L)
stopifnot(coverage_values[["Evidence rows matched one-to-one"]] == 1175L)

stopifnot(!any(grepl("GeneSymbol|Symbol", c("BaseAccession", "UMAP_1", "UMAP_2"))))

message("PASS: signed pathway pies use all 507 fixed BaseAccession coordinates.")
