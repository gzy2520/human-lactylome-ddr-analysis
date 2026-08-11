#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
script_path <- if (length(file_arg) > 0L) {
  normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE)
} else {
  normalizePath(
    "reanalysis/tests/test_kla_ddr_pathway_pie_umap_v2_spread.R",
    mustWork = TRUE
  )
}
project_root <- normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = TRUE)
table_dir <- file.path(
  project_root,
  "reanalysis/results/tables/kla_ddr_pathway_pie_umap_33groups_v2_spread"
)
figure_dir <- file.path(
  project_root,
  "reanalysis/results/figures/kla_ddr_pathway_pie_umap_33groups_v2_spread"
)

required_outputs <- c(
  file.path(table_dir, "umap_coordinates_v2_spread.csv"),
  file.path(table_dir, "pathway_umap_plot_data_v2_spread.csv"),
  file.path(table_dir, "pathway_assignment_long_v2.csv"),
  file.path(table_dir, "hollow_negative_wedge_polygons.csv"),
  file.path(table_dir, "pathway_assignment_summary.csv"),
  file.path(table_dir, "pie_overlap_comparison_v1_v2.csv"),
  file.path(table_dir, "umap_v2_parameters.csv"),
  file.path(table_dir, "npg_pathway_color_key.csv"),
  file.path(table_dir, "input_file_audit.csv"),
  file.path(table_dir, "session_info.txt"),
  file.path(figure_dir, "kla_ddr_pathway_hollow_negative_pie_umap_v2_spread.png"),
  file.path(figure_dir, "kla_ddr_pathway_hollow_negative_pie_umap_v2_spread.pdf"),
  file.path(figure_dir, "kla_ddr_pathway_hollow_negative_pie_umap_v2_spread.svg"),
  file.path(figure_dir, "kla_ddr_pathway_hollow_negative_umap_and_summary_v2_spread.png"),
  file.path(figure_dir, "kla_ddr_pathway_hollow_negative_umap_and_summary_v2_spread.pdf"),
  file.path(figure_dir, "kla_ddr_pathway_hollow_negative_umap_and_summary_v2_spread.svg"),
  file.path(project_root, "reanalysis/reports/UMAP_PATHWAY_PIE_33GROUP_V2_SPREAD.md")
)
stopifnot(all(file.exists(required_outputs)))
stopifnot(all(file.info(required_outputs)$size > 0L))

coordinates <- fread(file.path(table_dir, "umap_coordinates_v2_spread.csv"))
plot_data <- fread(file.path(table_dir, "pathway_umap_plot_data_v2_spread.csv"))
assignments <- fread(file.path(table_dir, "pathway_assignment_long_v2.csv"))
wedges <- fread(file.path(table_dir, "hollow_negative_wedge_polygons.csv"))
overlap <- fread(file.path(table_dir, "pie_overlap_comparison_v1_v2.csv"))
parameters <- fread(file.path(table_dir, "umap_v2_parameters.csv"))
colors <- fread(file.path(table_dir, "npg_pathway_color_key.csv"))

stopifnot(nrow(coordinates) == 507L)
stopifnot(uniqueN(coordinates$BaseAccession) == 507L)
stopifnot(all(is.finite(coordinates$UMAP_1)))
stopifnot(all(is.finite(coordinates$UMAP_2)))

stopifnot(nrow(plot_data) == 507L)
stopifnot(setequal(plot_data$BaseAccession, coordinates$BaseAccession))
stopifnot(plot_data[TotalAssignmentCount == 0L, .N] == 22L)

stopifnot(nrow(assignments) == 1175L)
stopifnot(sum(assignments$Score == 1L) == 1108L)
stopifnot(sum(assignments$Score == -1L) == 67L)
stopifnot(uniqueN(wedges$WedgeID) == 1175L)
stopifnot(uniqueN(wedges[Score == -1L]$WedgeID) == 67L)
stopifnot(all(wedges[Score == -1L]$FillKey == "NegativeHollow"))
stopifnot(all(wedges[Score == -1L]$BorderKey == wedges[Score == -1L]$Pathway))
stopifnot(all(wedges[Score == 1L]$BorderKey == "PositiveBoundary"))

stopifnot(nrow(colors) == 9L)
stopifnot(uniqueN(colors$Color) == 9L)
stopifnot(all(grepl("^#[0-9A-F]{6}$", colors$Color)))

parameter_values <- setNames(parameters$Value, parameters$Parameter)
stopifnot(parameter_values[["UMAPInput"]] == "507 x 66 protein x raw GO term binary matrix only")
stopifnot(parameter_values[["PathwayScoresUsedInUMAP"]] == "FALSE")
stopifnot(parameter_values[["NNeighbors"]] == "15")
stopifnot(parameter_values[["MinDist"]] == "3")
stopifnot(parameter_values[["Spread"]] == "10")
stopifnot(parameter_values[["NegativeSectorStyle"]] == "white hollow sector with pathway-colored border")

v1_pairs <- overlap[
  Version == "V1_fixed_compact",
  PairDistancesCloserThanPieDiameter
]
v2_pairs <- overlap[
  Version == "V2_refitted_spread",
  PairDistancesCloserThanPieDiameter
]
stopifnot(v1_pairs > 1000L)
stopifnot(v2_pairs <= 15L)
stopifnot(v2_pairs < v1_pairs)

message("PASS: V2 uses refitted spread coordinates and hollow negative sectors.")
