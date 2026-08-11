#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

project_root <- normalizePath(".", mustWork = TRUE)
analysis_name <- "kla_ddr_cytoscape_pathway_clusters_33groups_v1"
table_dir <- file.path(project_root, "reanalysis/results/tables", analysis_name)
figure_dir <- file.path(project_root, "reanalysis/results/figures", analysis_name)

required_tables <- file.path(
  table_dir,
  c(
    "cytoscape_string_nodes_507.csv",
    "cytoscape_string_edges_4458.csv",
    "string_mcl_i3_cluster_sizes.csv",
    "cytoscape_cluster_color_key.csv",
    "cytoscape_network_manifest.csv",
    "cytoscape_layout_string_mcl_i3_shared_string_topology.csv",
    "cytoscape_layout_pathway_score_k7_seed25_shared_string_topology.csv",
    "cytoscape_layout_pathway_score_k9_seed25_shared_string_topology.csv",
    "kla_ddr_507_cluster_comparison_shared_string_topology_33groups_v2.cys"
  )
)
required_figures <- file.path(
  figure_dir,
  as.vector(outer(
    c(
      "kla_ddr_507_string_mcl_i3_shared_string_topology",
      "kla_ddr_507_pathway_score_k7_seed25_shared_string_topology",
      "kla_ddr_507_pathway_score_k9_seed25_shared_string_topology"
    ),
    c(".png", ".svg", ".pdf"),
    paste0
  ))
)
required_files <- c(required_tables, required_figures)
if (any(!file.exists(required_files))) {
  stop(
    "Missing required Cytoscape output(s): ",
    paste(required_files[!file.exists(required_files)], collapse = "; ")
  )
}

assert <- function(condition, message) {
  if (!isTRUE(condition)) {
    stop(message, call. = FALSE)
  }
}

nodes <- fread(file.path(table_dir, "cytoscape_string_nodes_507.csv"))
edges <- fread(file.path(table_dir, "cytoscape_string_edges_4458.csv"))
manifest <- fread(file.path(table_dir, "cytoscape_network_manifest.csv"))
color_key <- fread(file.path(table_dir, "cytoscape_cluster_color_key.csv"))

assert(
  nrow(nodes) == 507L &&
    uniqueN(nodes$BaseAccession) == 507L &&
    sum(nodes$STRING_Mapped) == 506L &&
    identical(nodes[STRING_Mapped == FALSE, BaseAccession], "A8MVJ9"),
  "The exported Cytoscape node table does not reproduce the 506+1 STRING scope."
)
assert(
  nrow(edges) == 4458L &&
    uniqueN(edges$EdgeID) == 4458L &&
    all(edges$Source %in% nodes$BaseAccession) &&
    all(edges$Target %in% nodes$BaseAccession) &&
    all(edges$STRING_score >= 0.70 & edges$STRING_score <= 1),
  "The exported high-confidence STRING edge table is invalid."
)
assert(
  sum(!is.na(nodes$MCL_cluster_i3)) == 446L &&
    uniqueN(na.omit(nodes$MCL_cluster_i3)) == 65L &&
    identical(sort(unique(nodes$K7_Cluster)), 1:7) &&
    identical(sort(unique(nodes$K9_Cluster)), 1:9),
  "Expected MCL i=3.0, k=7, and k=9 cluster labels are not all present."
)
assert(
  nrow(manifest) == 3L &&
    setequal(
      manifest$Version,
      c(
        "string_mcl_i3_shared_string_topology",
        "pathway_score_k7_seed25_shared_string_topology",
        "pathway_score_k9_seed25_shared_string_topology"
      )
    ) &&
    all(manifest$NodeCount == 507L) &&
    all(manifest$EdgeCount == 4458L),
  "The Cytoscape network manifest does not contain the three comparison versions."
)
expected_styles <- data.table(
  Version = c(
    "string_mcl_i3_shared_string_topology",
    "pathway_score_k7_seed25_shared_string_topology",
    "pathway_score_k9_seed25_shared_string_topology"
  ),
  StyleName = c(
    "Kla_DDR_STRING_MCL_i3_overview_v4",
    "Kla_DDR_pathway_score_k7_seed25_overview_v4",
    "Kla_DDR_pathway_score_k9_seed25_overview_v4"
  )
)
assert(
  nrow(merge(
    manifest[, .(Version, StyleName)],
    expected_styles,
    by = c("Version", "StyleName")
  )) == 3L,
  "The manifest does not assign a distinct expected style to each network."
)
assert(
  color_key[Version == "PathwayScore_k7_seed25", uniqueN(Cluster)] == 7L &&
    color_key[Version == "PathwayScore_k9_seed25", uniqueN(Cluster)] == 9L,
  "The saved cluster color key is incomplete."
)

k7_svg <- paste(
  readLines(
    file.path(
      figure_dir,
      "kla_ddr_507_pathway_score_k7_seed25_shared_string_topology.svg"
    ),
    warn = FALSE
  ),
  collapse = "\n"
)
k9_svg <- paste(
  readLines(
    file.path(
      figure_dir,
      "kla_ddr_507_pathway_score_k9_seed25_shared_string_topology.svg"
    ),
    warn = FALSE
  ),
  collapse = "\n"
)
assert(
  !grepl("#7f3c8d|#999999", k7_svg, ignore.case = TRUE) &&
    grepl("#7f3c8d", k9_svg, ignore.case = TRUE) &&
    grepl("#999999", k9_svg, ignore.case = TRUE),
  "The seven- and nine-pathway SVG exports do not use their distinct palettes."
)

layout_files <- required_tables[grepl("cytoscape_layout_", required_tables)]
layout_tables <- list()
for (layout_path in layout_files) {
  coordinates <- fread(layout_path)
  setorder(coordinates, BaseAccession)
  layout_tables[[basename(layout_path)]] <- coordinates
  assert(
    nrow(coordinates) == 507L &&
      uniqueN(coordinates$BaseAccession) == 507L &&
      uniqueN(coordinates[, .(X, Y)]) == 507L &&
      all(is.finite(coordinates$X)) &&
      all(is.finite(coordinates$Y)),
    paste("Invalid Cytoscape layout coordinates:", basename(layout_path))
  )
}
reference_layout <- layout_tables[[1L]]
for (coordinates in layout_tables[-1L]) {
  assert(
    identical(coordinates$BaseAccession, reference_layout$BaseAccession) &&
      max(abs(coordinates$X - reference_layout$X)) < 1e-8 &&
      max(abs(coordinates$Y - reference_layout$Y)) < 1e-8,
    "The three Cytoscape views do not use identical per-protein coordinates."
  )
}
assert(
  all(file.info(required_figures)$size > 2000) &&
    file.info(file.path(
      table_dir,
      "kla_ddr_507_cluster_comparison_shared_string_topology_33groups_v2.cys"
    ))$size > 1e6,
  "A Cytoscape figure or saved session is unexpectedly small."
)

cat(
  paste0(
    "PASS: three 507-node/4,458-edge Cytoscape networks, styles, layouts, ",
    "identical shared coordinates, figures, and the saved session are valid.\n"
  )
)
