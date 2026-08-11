#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
script_path <- if (length(file_arg) > 0L) {
  normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE)
} else {
  normalizePath(
    "reanalysis/tests/test_kla_ddr_pathway_pie_umap_v4_bp_semantic.R",
    mustWork = TRUE
  )
}
project_root <- normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = TRUE)

v3_table_dir <- file.path(
  project_root,
  "reanalysis/results/tables/kla_ddr_pathway_pie_umap_33groups_v3_all_go"
)
table_dir <- file.path(
  project_root,
  "reanalysis/results/tables/kla_ddr_pathway_pie_umap_33groups_v4_bp_semantic"
)
figure_dir <- file.path(
  project_root,
  "reanalysis/results/figures/kla_ddr_pathway_pie_umap_33groups_v4_bp_semantic"
)

assert <- function(condition, message) {
  if (!isTRUE(condition)) {
    stop(message, call. = FALSE)
  }
}

required_files <- c(
  file.path(
    project_root,
    "reanalysis/config/uniprot_kla_ddr_507_all_go_2026-08-10.tsv"
  ),
  file.path(v3_table_dir, "uniprot_direct_go_annotation_long.csv"),
  file.path(v3_table_dir, "protein_all_go_direct_binary_matrix.csv"),
  file.path(table_dir, "bp_direct_to_ancestor_expansion_long.csv"),
  file.path(table_dir, "bp_direct_terms_missing_from_go_db.csv"),
  file.path(table_dir, "bp_semantic_feature_dictionary.csv"),
  file.path(table_dir, "protein_bp_semantic_feature_long.csv"),
  file.path(table_dir, "protein_bp_semantic_shared_binary_matrix.csv"),
  file.path(table_dir, "umap_raw_coordinates_v4_bp_semantic.csv"),
  file.path(table_dir, "umap_display_coordinates_v4_bp_semantic.csv"),
  file.path(table_dir, "umap_raw_scaled_display_coordinate_audit.csv"),
  file.path(table_dir, "pie_spacing_before_after_collision.csv"),
  file.path(table_dir, "pathway_umap_plot_data_v4_bp_semantic.csv"),
  file.path(table_dir, "hollow_negative_wedge_polygons_v4.csv"),
  file.path(table_dir, "umap_v4_parameters.csv"),
  file.path(table_dir, "pathway_color_key_v4.csv"),
  file.path(
    figure_dir,
    "kla_ddr_pathway_pie_umap_v4_bp_semantic.png"
  ),
  file.path(
    figure_dir,
    "kla_ddr_pathway_pie_umap_v4_bp_semantic.pdf"
  ),
  file.path(
    figure_dir,
    "kla_ddr_pathway_pie_umap_v4_bp_semantic.svg"
  ),
  file.path(
    figure_dir,
    "kla_ddr_raw_bp_semantic_umap_v4.png"
  ),
  file.path(
    figure_dir,
    "kla_ddr_raw_and_pie_bp_semantic_umap_v4.png"
  )
)
assert(all(file.exists(required_files)), "At least one required V4 file is missing.")
assert(all(file.info(required_files)$size > 0), "At least one required V4 file is empty.")

direct_go <- fread(file.path(v3_table_dir, "uniprot_direct_go_annotation_long.csv"))
full_matrix <- fread(
  file.path(v3_table_dir, "protein_all_go_direct_binary_matrix.csv")
)
expansion <- fread(
  file.path(table_dir, "bp_direct_to_ancestor_expansion_long.csv")
)
missing_go_db <- fread(
  file.path(table_dir, "bp_direct_terms_missing_from_go_db.csv")
)
feature_dictionary <- fread(
  file.path(table_dir, "bp_semantic_feature_dictionary.csv")
)
semantic_long <- fread(
  file.path(table_dir, "protein_bp_semantic_feature_long.csv")
)
semantic_matrix <- fread(
  file.path(table_dir, "protein_bp_semantic_shared_binary_matrix.csv")
)
raw_coordinates <- fread(
  file.path(table_dir, "umap_raw_coordinates_v4_bp_semantic.csv")
)
display_coordinates <- fread(
  file.path(table_dir, "umap_display_coordinates_v4_bp_semantic.csv")
)
coordinate_audit <- fread(
  file.path(table_dir, "umap_raw_scaled_display_coordinate_audit.csv")
)
spacing <- fread(
  file.path(table_dir, "pie_spacing_before_after_collision.csv")
)
plot_data <- fread(
  file.path(table_dir, "pathway_umap_plot_data_v4_bp_semantic.csv")
)
wedges <- fread(
  file.path(table_dir, "hollow_negative_wedge_polygons_v4.csv")
)
parameters <- fread(file.path(table_dir, "umap_v4_parameters.csv"))
colors <- fread(file.path(table_dir, "pathway_color_key_v4.csv"))

assert(
  nrow(direct_go) == 13738L &&
    uniqueN(direct_go$GO_ID) == 3461L &&
    uniqueN(direct_go$BaseAccession) == 507L,
  "Saved full direct-GO long table is wrong."
)
assert(
  nrow(full_matrix) == 507L && ncol(full_matrix) == 3462L,
  "Saved full direct-GO matrix is not 507 x (1 + 3,461)."
)
assert(
  uniqueN(expansion$BaseAccession) == 507L &&
    all(expansion$FeatureSource %in% c("Direct", "Ancestor")),
  "BP direct-to-ancestor expansion is invalid."
)
assert(
  nrow(missing_go_db) == 4L &&
    setequal(
      missing_go_db$GO_ID,
      c("GO:0160217", "GO:0160234", "GO:0160240", "GO:0160276")
    ),
  "GO IDs unresolved against GO.db were not audited correctly."
)
assert(
  feature_dictionary[EligibleForUMAP == TRUE, .N] == 3008L,
  "Expected 3,008 shared BP semantic features."
)
assert(
  all(
    feature_dictionary[
      EligibleForUMAP == TRUE,
      ProteinCount
    ] >= 2L
  ) &&
    all(
      feature_dictionary[
        EligibleForUMAP == TRUE,
        ProteinCount
      ] <= 405L
    ),
  "BP semantic feature frequency filter is wrong."
)
assert(
  uniqueN(semantic_long$BaseAccession) == 507L &&
    uniqueN(semantic_long$FeatureGO_ID) == 3008L,
  "BP semantic long table dimensions are wrong."
)
assert(
  nrow(semantic_matrix) == 507L && ncol(semantic_matrix) == 3009L,
  "BP semantic binary matrix is not 507 x (1 + 3,008)."
)
binary <- as.matrix(semantic_matrix[, -1L])
storage.mode(binary) <- "numeric"
assert(all(binary %in% c(0, 1)), "BP semantic matrix is not binary.")
assert(
  setequal(semantic_matrix$BaseAccession, direct_go$BaseAccession),
  "BP semantic matrix protein set differs from the full-GO table."
)
assert(
  nrow(raw_coordinates) == 507L &&
    nrow(display_coordinates) == 507L &&
    uniqueN(raw_coordinates$BaseAccession) == 507L &&
    uniqueN(display_coordinates$BaseAccession) == 507L,
  "Raw or display coordinate table is incomplete."
)
assert(
  all(
    is.finite(
      as.matrix(raw_coordinates[, .(UMAP_1, UMAP_2)])
    )
  ) &&
    all(
      is.finite(
        as.matrix(
          display_coordinates[, .(UMAP_Display_1, UMAP_Display_2)]
        )
      )
    ),
  "Raw or display coordinates contain non-finite values."
)
assert(
  spacing[
    CoordinateSet == "Collision-resolved display centers",
    PairDistancesBelowThreshold
  ] == 0L &&
    spacing[
      CoordinateSet == "Collision-resolved display centers",
      MinimumNearestNeighborDistance
    ] >= 1.50 - 1e-6,
  "Enlarged V4 pies do not retain the required spacing."
)
assert(
  parameters[Parameter == "PieRadius", Value] == "0.75" &&
    parameters[Parameter == "MinimumCenterSeparation", Value] == "1.5",
  "Enlarged pie radius or center separation is wrong."
)
pair_distance_spearman <- as.numeric(
  parameters[Parameter == "RawDisplayPairDistanceSpearman", Value]
)
assert(
  is.finite(pair_distance_spearman) && pair_distance_spearman > 0.99,
  "Minimal overlap removal did not preserve the required UMAP distance structure."
)
assert(
  nrow(coordinate_audit) == 507L &&
    all(coordinate_audit$DisplayDisplacement >= 0),
  "Coordinate displacement audit is invalid."
)
assert(
  plot_data[TotalAssignmentCount == 0L, .N] == 22L,
  "Expected 22 circled-cross proteins without a scored pathway/function."
)
assert(
  uniqueN(wedges$WedgeID) == 1175L &&
    uniqueN(wedges[Score == -1L]$WedgeID) == 67L &&
    all(wedges[Score == -1L]$FillKey == "NegativeHollow") &&
    all(wedges[Score == 1L]$BorderKey == wedges[Score == 1L]$Pathway),
  "Signed hollow/filled wedge encoding is wrong."
)
assert(
  parameters[Parameter == "EmbeddingGOAspect", Value] == "BP" &&
    parameters[Parameter == "PathwayScoresUsedInUMAP", Value] == "FALSE" &&
    parameters[Parameter == "SampleDetectionUsedInUMAP", Value] == "FALSE",
  "V4 embedding inputs are not declared correctly."
)
assert(
  nrow(colors) == 9L &&
    uniqueN(colors$Color) == 9L &&
    colors[DisplayLabel == "Chromatin interaction", Color] == "#7E6148" &&
    colors[DisplayLabel == "Other support", Color] == "#6F6F6F",
  "Corrected V4 pathway colors are wrong."
)

message(
  "PASS: V4 preserves full GO archives, uses shared BP hierarchy features, ",
  "separates enlarged pies, and keeps raw/display topology roles explicit."
)
