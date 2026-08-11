#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(grid)
  library(png)
})

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
script_path <- if (length(file_arg) > 0L) {
  normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE)
} else {
  normalizePath(
    "reanalysis/scripts/compare_kla_ddr_umap_collision_variants.R",
    mustWork = TRUE
  )
}
project_root <- normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = TRUE)

standard_dir <- file.path(
  project_root,
  "reanalysis/results/figures/kla_ddr_pathway_pie_umap_33groups_v4_bp_semantic"
)
collision_dir <- file.path(
  project_root,
  paste0(
    "reanalysis/results/figures/",
    "kla_ddr_pathway_pie_umap_33groups_v4_bp_semantic_collision_heavy"
  )
)
standard_table_dir <- file.path(
  project_root,
  "reanalysis/results/tables/kla_ddr_pathway_pie_umap_33groups_v4_bp_semantic"
)
collision_table_dir <- file.path(
  project_root,
  paste0(
    "reanalysis/results/tables/",
    "kla_ddr_pathway_pie_umap_33groups_v4_bp_semantic_collision_heavy"
  )
)
output_dir <- file.path(
  project_root,
  paste0(
    "reanalysis/results/figures/",
    "kla_ddr_pathway_pie_umap_33groups_v4_bp_semantic_comparison"
  )
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

plot_name <- "kla_ddr_pathway_pie_umap_v4_bp_semantic.png"
standard_path <- file.path(standard_dir, plot_name)
collision_path <- file.path(collision_dir, plot_name)
standard_parameters_path <- file.path(standard_table_dir, "umap_v4_parameters.csv")
collision_parameters_path <- file.path(collision_table_dir, "umap_v4_parameters.csv")
required_paths <- c(
  standard_path,
  collision_path,
  standard_parameters_path,
  collision_parameters_path
)
if (any(!file.exists(required_paths))) {
  stop(
    "Missing comparison input(s): ",
    paste(required_paths[!file.exists(required_paths)], collapse = "; ")
  )
}

standard_image <- readPNG(standard_path)
collision_image <- readPNG(collision_path)
standard_parameters <- fread(standard_parameters_path)
collision_parameters <- fread(collision_parameters_path)

parameter_value <- function(parameters, name) {
  value <- parameters[Parameter == name, Value]
  if (length(value) != 1L) {
    stop("Expected exactly one parameter value for ", name)
  }
  value
}

standard_rho <- as.numeric(
  parameter_value(standard_parameters, "RawDisplayPairDistanceSpearman")
)
collision_rho <- as.numeric(
  parameter_value(collision_parameters, "RawDisplayPairDistanceSpearman")
)
standard_separation <- as.numeric(
  parameter_value(standard_parameters, "MinimumCenterSeparation")
)
collision_separation <- as.numeric(
  parameter_value(collision_parameters, "MinimumCenterSeparation")
)
pie_radius <- as.numeric(parameter_value(standard_parameters, "PieRadius"))

panel_labels <- c(
  paste0(
    "A  UMAP-preserving display  |  center distance ≥ ",
    format(standard_separation, digits = 3),
    "  |  edge gap = ",
    format(standard_separation - 2 * pie_radius, digits = 3),
    "  |  distance rho = ",
    format(standard_rho, digits = 4)
  ),
  paste0(
    "B  Collision-heavy display  |  center distance ≥ ",
    format(collision_separation, digits = 3),
    "  |  edge gap = ",
    format(collision_separation - 2 * pie_radius, digits = 3),
    "  |  distance rho = ",
    format(collision_rho, digits = 4)
  )
)

draw_comparison <- function() {
  grid.newpage()
  pushViewport(
    viewport(
      layout = grid.layout(
        nrow = 2L,
        ncol = 2L,
        heights = unit(c(0.07, 0.93), "npc"),
        widths = unit(c(0.5, 0.5), "npc")
      )
    )
  )
  for (column_index in 1:2) {
    pushViewport(
      viewport(layout.pos.row = 1L, layout.pos.col = column_index)
    )
    grid.text(
      panel_labels[[column_index]],
      x = unit(0.5, "npc"),
      y = unit(0.46, "npc"),
      gp = gpar(
        fontfamily = "sans",
        fontface = "bold",
        fontsize = 12,
        col = "#222222"
      )
    )
    popViewport()
  }
  images <- list(standard_image, collision_image)
  for (column_index in 1:2) {
    pushViewport(
      viewport(layout.pos.row = 2L, layout.pos.col = column_index)
    )
    grid.raster(
      images[[column_index]],
      x = unit(0.5, "npc"),
      y = unit(0.5, "npc"),
      width = unit(1, "npc"),
      height = unit(1, "npc"),
      interpolate = TRUE
    )
    grid.rect(gp = gpar(fill = NA, col = "#B8B8B8", lwd = 0.8))
    popViewport()
  }
  popViewport()
}

png_path <- file.path(
  output_dir,
  "kla_ddr_umap_collision_strategy_side_by_side.png"
)
pdf_path <- file.path(
  output_dir,
  "kla_ddr_umap_collision_strategy_side_by_side.pdf"
)

png(
  png_path,
  width = 7200,
  height = 2850,
  res = 300,
  bg = "white"
)
draw_comparison()
dev.off()

pdf(pdf_path, width = 24, height = 9.5, bg = "white", useDingbats = FALSE)
draw_comparison()
dev.off()

message("Comparison PNG: ", png_path)
message("Comparison PDF: ", pdf_path)
