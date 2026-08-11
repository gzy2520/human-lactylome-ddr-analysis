#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

project_root <- normalizePath(".", mustWork = TRUE)
table_dir <- file.path(
  project_root,
  "reanalysis/results/tables/kla_ddr_cytoscape_pathway_clusters_33groups_v1"
)

required_files <- file.path(
  table_dir,
  c(
    "pathway_score_matrix_7.csv",
    "pathway_score_matrix_9.csv",
    "cytoscape_node_import_table.csv",
    "pathway_cluster_assignments_seed25.csv",
    "pathway_cluster_centers_k7.csv",
    "pathway_cluster_centers_k9.csv",
    "pathway_cluster_summary_k7.csv",
    "pathway_cluster_summary_k9.csv",
    "pathway_score_pattern_frequencies_k7.csv",
    "pathway_score_pattern_frequencies_k9.csv",
    "pathway_cluster_diagnostics.csv",
    "input_file_audit.csv",
    "session_info.txt"
  )
)
if (any(!file.exists(required_files))) {
  stop(
    "Missing required output(s): ",
    paste(required_files[!file.exists(required_files)], collapse = "; ")
  )
}

assert <- function(condition, message) {
  if (!isTRUE(condition)) {
    stop(message, call. = FALSE)
  }
}

nodes <- fread(file.path(table_dir, "cytoscape_node_import_table.csv"))
assignments <- fread(
  file.path(table_dir, "pathway_cluster_assignments_seed25.csv")
)
diagnostics <- fread(file.path(table_dir, "pathway_cluster_diagnostics.csv"))
k7_centers <- fread(file.path(table_dir, "pathway_cluster_centers_k7.csv"))
k9_centers <- fread(file.path(table_dir, "pathway_cluster_centers_k9.csv"))
k7_summary <- fread(file.path(table_dir, "pathway_cluster_summary_k7.csv"))
k9_summary <- fread(file.path(table_dir, "pathway_cluster_summary_k9.csv"))

score_columns <- c(
  "Score_HR",
  "Score_NHEJ",
  "Score_AEJ",
  "Score_BER",
  "Score_NER",
  "Score_MMR",
  "Score_FA",
  "Score_Chromatin_interaction",
  "Score_Other_support"
)

assert(
  nrow(nodes) == 507L &&
    uniqueN(nodes$BaseAccession) == 507L &&
    nrow(assignments) == 507L &&
    setequal(nodes$BaseAccession, assignments$BaseAccession),
  "Expected the fixed 507 unique BaseAccession proteins."
)
assert(
  all(as.matrix(nodes[, ..score_columns]) %in% c(-1L, 0L, 1L)),
  "Pathway scores must use only -1/0/+1."
)
assert(
  sum(nodes$K7_AllZero) == 189L &&
    sum(nodes$K9_AllZero) == 22L,
  "Unexpected all-zero protein counts for the 7D or 9D score matrix."
)
assert(
  identical(sort(unique(nodes$K7_Cluster)), 1:7) &&
    identical(sort(unique(nodes$K9_Cluster)), 1:9) &&
    nrow(k7_centers) == 7L &&
    nrow(k9_centers) == 9L &&
    nrow(k7_summary) == 7L &&
    nrow(k9_summary) == 9L,
  "The requested k=7 and k=9 solutions were not both generated."
)
assert(
  setequal(diagnostics$DimensionCount, c(7L, 9L)) &&
    all(diagnostics$RandomSeed == 25L) &&
    all(diagnostics$ProteinCount == 507L) &&
    diagnostics[DimensionCount == 7L, AllZeroProteinCount] == 189L &&
    diagnostics[DimensionCount == 9L, AllZeroProteinCount] == 22L,
  "Clustering diagnostics do not reproduce the fixed seed-25 scope."
)
assert(
  abs(diagnostics[DimensionCount == 7L, MeanSilhouette] - 0.6045316) < 1e-6 &&
    abs(diagnostics[DimensionCount == 9L, MeanSilhouette] - 0.3999890) < 1e-6,
  "Seed-25 k-means solutions changed unexpectedly."
)

cat(
  paste0(
    "PASS: fixed 507-protein signed pathway matrices and reproducible ",
    "seed-25 k=7/k=9 clusters are valid.\n"
  )
)
