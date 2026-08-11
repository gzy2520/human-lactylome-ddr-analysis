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
    "reanalysis/tests/test_kla_ddr_pathway_pie_umap_v3_all_go.R",
    mustWork = TRUE
  )
}
project_root <- normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = TRUE)

table_dir <- file.path(
  project_root,
  "reanalysis/results/tables/kla_ddr_pathway_pie_umap_33groups_v3_all_go"
)
figure_dir <- file.path(
  project_root,
  "reanalysis/results/figures/kla_ddr_pathway_pie_umap_33groups_v3_all_go"
)
annotation_path <- file.path(
  project_root,
  "reanalysis/config/uniprot_kla_ddr_507_all_go_2026-08-10.tsv"
)
metadata_path <- file.path(
  project_root,
  "reanalysis/config/uniprot_kla_ddr_507_all_go_2026-08-10_metadata.tsv"
)

assert <- function(condition, message) {
  if (!isTRUE(condition)) {
    stop(message, call. = FALSE)
  }
}

required_files <- c(
  annotation_path,
  metadata_path,
  file.path(table_dir, "uniprot_direct_go_annotation_long.csv"),
  file.path(table_dir, "all_go_term_dictionary_and_idf.csv"),
  file.path(table_dir, "protein_all_go_direct_binary_matrix.csv"),
  file.path(table_dir, "protein_all_go_profile_summary.csv"),
  file.path(table_dir, "duplicate_all_go_profiles.csv"),
  file.path(table_dir, "umap_coordinates_v3_all_go.csv"),
  file.path(table_dir, "pathway_umap_plot_data_v3_all_go.csv"),
  file.path(table_dir, "pathway_assignment_long_v3.csv"),
  file.path(table_dir, "hollow_negative_wedge_polygons_v3.csv"),
  file.path(table_dir, "pie_overlap_comparison_v1_v2_v3.csv"),
  file.path(table_dir, "umap_v3_parameters.csv"),
  file.path(table_dir, "pathway_color_key_v3.csv"),
  file.path(
    figure_dir,
    "kla_ddr_pathway_hollow_negative_pie_umap_v3_all_go.png"
  ),
  file.path(
    figure_dir,
    "kla_ddr_pathway_hollow_negative_pie_umap_v3_all_go.pdf"
  ),
  file.path(
    figure_dir,
    "kla_ddr_pathway_hollow_negative_pie_umap_v3_all_go.svg"
  )
)
assert(all(file.exists(required_files)), "At least one required V3 output is missing.")
assert(all(file.info(required_files)$size > 0), "At least one required V3 output is empty.")

annotation <- fread(
  annotation_path,
  sep = "\t",
  quote = "\"",
  na.strings = character(),
  check.names = FALSE
)
metadata <- fread(metadata_path, sep = "\t")
go_long <- fread(file.path(table_dir, "uniprot_direct_go_annotation_long.csv"))
term_dictionary <- fread(file.path(table_dir, "all_go_term_dictionary_and_idf.csv"))
go_binary <- fread(file.path(table_dir, "protein_all_go_direct_binary_matrix.csv"))
profiles <- fread(file.path(table_dir, "protein_all_go_profile_summary.csv"))
duplicates <- fread(file.path(table_dir, "duplicate_all_go_profiles.csv"))
coordinates <- fread(file.path(table_dir, "umap_coordinates_v3_all_go.csv"))
plot_data <- fread(file.path(table_dir, "pathway_umap_plot_data_v3_all_go.csv"))
assignments <- fread(file.path(table_dir, "pathway_assignment_long_v3.csv"))
wedges <- fread(file.path(table_dir, "hollow_negative_wedge_polygons_v3.csv"))
overlap <- fread(file.path(table_dir, "pie_overlap_comparison_v1_v2_v3.csv"))
parameters <- fread(file.path(table_dir, "umap_v3_parameters.csv"))
colors <- fread(file.path(table_dir, "pathway_color_key_v3.csv"))

assert(
  nrow(annotation) == 507L && uniqueN(annotation$BaseAccession) == 507L,
  "Raw UniProt cache is not 507 unique proteins."
)
assert(
  metadata[Key == "UniProtRelease", Value] == "2026_02" &&
    metadata[Key == "OrganismTaxonID", Value] == "9606",
  "UniProt release/taxon metadata is wrong."
)
assert(
  nrow(go_long) == 13738L &&
    uniqueN(go_long$GO_ID) == 3461L &&
    uniqueN(go_long$BaseAccession) == 507L,
  "Long GO table dimensions are wrong."
)
assert(
  setequal(unique(go_long$Ontology), c("BP", "CC", "MF")),
  "Long GO table does not contain all three GO aspects."
)
assert(
  nrow(term_dictionary) == 3461L &&
    all(term_dictionary$ProteinCount >= 1L) &&
    all(term_dictionary$IDFWeight >= 1),
  "GO dictionary or IDF weights are invalid."
)
assert(
  nrow(go_binary) == 507L && ncol(go_binary) == 3462L,
  "Saved direct-GO binary matrix is not 507 x (1 + 3,461)."
)
binary_matrix <- as.matrix(go_binary[, -1L])
storage.mode(binary_matrix) <- "numeric"
assert(all(binary_matrix %in% c(0, 1)), "Saved GO matrix is not binary.")
assert(sum(binary_matrix) == 13738L, "Saved GO matrix hit count is wrong.")
assert(
  setequal(go_binary$BaseAccession, annotation$BaseAccession),
  "Saved GO matrix protein set differs from raw UniProt cache."
)
assert(
  uniqueN(profiles$DirectGOSignature) == 506L,
  "Expected 506 unique all-GO profiles."
)
assert(
  nrow(duplicates) == 1L &&
    duplicates$ProteinsWithSignature[[1L]] == 2L &&
    duplicates$Accessions[[1L]] == "P35250;P40937",
  "Expected only the RFC2/RFC5 direct-GO duplicate profile."
)
assert(
  nrow(coordinates) == 507L &&
    uniqueN(coordinates$BaseAccession) == 507L &&
    uniqueN(coordinates[, .(UMAP_1, UMAP_2)]) == 507L &&
    all(is.finite(as.matrix(coordinates[, .(UMAP_1, UMAP_2)]))),
  "V3 UMAP coordinates are incomplete, duplicated, or non-finite."
)
assert(
  overlap[
    Version == "V3_all_GO_TFIDF_spread",
    PairDistancesCloserThanPieDiameter
  ] == 0L,
  "V3 still contains overlapping pies at the selected radius."
)
assert(
  nrow(assignments) == 1175L &&
    sum(assignments$Score == 1L) == 1108L &&
    sum(assignments$Score == -1L) == 67L,
  "Signed pathway assignment counts are wrong."
)
assert(
  uniqueN(wedges$WedgeID) == 1175L &&
    uniqueN(wedges[Score == -1L]$WedgeID) == 67L &&
    all(wedges[Score == -1L]$FillKey == "NegativeHollow"),
  "Hollow negative-sector construction is wrong."
)
assert(
  plot_data[TotalAssignmentCount == 0L, .N] == 22L,
  "Expected 22 proteins without a scored pathway/function."
)
assert(
  parameters[Parameter == "PathwayScoresUsedInUMAP", Value] == "FALSE" &&
    parameters[Parameter == "SampleDetectionUsedInUMAP", Value] == "FALSE",
  "Scores or sample detection were incorrectly declared as UMAP inputs."
)
assert(
  nrow(colors) == 9L && uniqueN(colors$Color) == 9L,
  "Pathway palette does not contain nine unique colors."
)
assert(
  colors[DisplayLabel == "Chromatin interaction", Color] == "#7E6148" &&
    colors[DisplayLabel == "Other support", Color] == "#6F6F6F",
  "Chromatin/Other colors were not separated as specified."
)

message(
  "PASS: V3 uses saved full UniProt GO tables, TF-IDF UMAP, ",
  "separated colors, hollow negative sectors, and circled-cross zero proteins."
)
