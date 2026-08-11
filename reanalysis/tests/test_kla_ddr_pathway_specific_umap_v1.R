#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

project_root <- normalizePath(".", mustWork = TRUE)
table_dir <- file.path(
  project_root,
  "reanalysis/results/tables/kla_ddr_pathway_specific_umap_33groups_v1"
)
figure_dir <- file.path(
  project_root,
  "reanalysis/results/figures/kla_ddr_pathway_specific_umap_33groups_v1"
)
v4_coordinate_path <- file.path(
  project_root,
  paste0(
    "reanalysis/results/tables/",
    "kla_ddr_pathway_pie_umap_33groups_v4_bp_semantic/",
    "umap_raw_coordinates_v4_bp_semantic.csv"
  )
)
assignment_path <- file.path(
  project_root,
  paste0(
    "reanalysis/results/tables/",
    "kla_ddr_pathway_pie_umap_33groups_v3_all_go/",
    "pathway_assignment_long_v3.csv"
  )
)

required_files <- c(
  file.path(table_dir, "pathway_specific_umap_plot_data.csv"),
  file.path(table_dir, "pathway_specific_umap_summary.csv"),
  file.path(table_dir, "input_file_audit.csv"),
  file.path(table_dir, "session_info.txt"),
  file.path(figure_dir, "kla_ddr_pathway_specific_umap_3x3_v1.png"),
  file.path(figure_dir, "kla_ddr_pathway_specific_umap_3x3_v1.pdf"),
  file.path(figure_dir, "kla_ddr_pathway_specific_umap_3x3_v1.svg"),
  v4_coordinate_path,
  assignment_path
)
if (any(!file.exists(required_files))) {
  stop(
    "Missing required output/input(s): ",
    paste(required_files[!file.exists(required_files)], collapse = "; ")
  )
}

assert <- function(condition, message) {
  if (!isTRUE(condition)) {
    stop(message, call. = FALSE)
  }
}

plot_data <- fread(file.path(table_dir, "pathway_specific_umap_plot_data.csv"))
summary_table <- fread(file.path(table_dir, "pathway_specific_umap_summary.csv"))
raw_coordinates <- fread(v4_coordinate_path)
assignment_long <- fread(assignment_path)

assert(
  nrow(plot_data) == 507L * 9L &&
    uniqueN(plot_data$BaseAccession) == 507L &&
    uniqueN(plot_data$DisplayLabel) == 9L,
  "Expected a complete 507-protein by nine-pathway plotting grid."
)
assert(
  plot_data[, .N, by = DisplayLabel][, all(N == 507L)] &&
    summary_table[, all(TotalProteins == 507L)] &&
    summary_table[, all(PromotingCount + SuppressingCount + UnassignedCount == 507L)],
  "Each panel must contain all 507 proteins."
)
assert(
  all(plot_data$Score %in% c(-1L, 0L, 1L)) &&
    plot_data[Score == 1L, all(Status == "Promoting (+1)")] &&
    plot_data[Score == -1L, all(Status == "Suppressing (-1)")] &&
    plot_data[Score == 0L, all(Status == "Not assigned (0)")],
  "The pathway-specific sign encoding is invalid."
)

coordinate_check <- merge(
  unique(plot_data[, .(BaseAccession, UMAP_1, UMAP_2)]),
  raw_coordinates,
  by = "BaseAccession",
  suffixes = c("_panel", "_raw")
)
assert(
  nrow(coordinate_check) == 507L &&
    isTRUE(
      all.equal(
        coordinate_check$UMAP_1_panel,
        coordinate_check$UMAP_1_raw,
        tolerance = 0
      )
    ) &&
    isTRUE(
      all.equal(
        coordinate_check$UMAP_2_panel,
        coordinate_check$UMAP_2_raw,
        tolerance = 0
      )
    ),
  "Panel coordinates are not the exact saved raw BP semantic UMAP coordinates."
)

nonzero_check <- merge(
  plot_data[Score != 0L, .(BaseAccession, DisplayLabel, Score)],
  assignment_long[, .(BaseAccession, DisplayLabel, Score)],
  by = c("BaseAccession", "DisplayLabel", "Score")
)
assert(
  nrow(nonzero_check) == nrow(assignment_long) &&
    plot_data[Score != 0L, .N] == nrow(assignment_long),
  "Nonzero panel assignments do not reproduce the signed source table."
)

individual_pngs <- list.files(
  file.path(figure_dir, "individual_pathways"),
  pattern = "_raw_umap_v1\\.png$",
  full.names = TRUE
)
individual_pdfs <- list.files(
  file.path(figure_dir, "individual_pathways"),
  pattern = "_raw_umap_v1\\.pdf$",
  full.names = TRUE
)
individual_svgs <- list.files(
  file.path(figure_dir, "individual_pathways"),
  pattern = "_raw_umap_v1\\.svg$",
  full.names = TRUE
)
assert(
  length(individual_pngs) == 9L &&
    length(individual_pdfs) == 9L &&
    length(individual_svgs) == 9L &&
    all(file.info(c(individual_pngs, individual_pdfs, individual_svgs))$size > 0),
  "Expected nine nonempty individual pathway figures in each format."
)

cat(
  paste0(
    "PASS: nine pathway-specific panels use all 507 proteins, exact raw BP ",
    "semantic UMAP coordinates, and signed solid/hollow encoding.\n"
  )
)
