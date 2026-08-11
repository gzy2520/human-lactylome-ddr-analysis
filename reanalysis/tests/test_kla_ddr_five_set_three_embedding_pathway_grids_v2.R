#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

project_root <- normalizePath(".", mustWork = TRUE)
analysis_name <- "kla_ddr_five_set_three_embedding_pathway_grids_33groups_v2"
table_dir <- file.path(project_root, "reanalysis/results/tables", analysis_name)

required_tables <- c(
  "embedding_coordinates_5sets_long.csv",
  "embedding_parameters_5sets.csv",
  "embedding_set_membership.csv",
  "pathway_plot_data_5sets.csv",
  "pathway_summary_5sets.csv",
  "figure_manifest_15_grids.csv",
  "parameter_tuning_metrics.csv",
  "parameter_tuning_candidate_summary.csv",
  "parameter_tuning_recommendation.csv",
  "input_file_audit.csv",
  "session_info.txt"
)
required_paths <- file.path(table_dir, required_tables)
if (any(!file.exists(required_paths))) {
  stop(
    "Missing required output(s): ",
    paste(required_paths[!file.exists(required_paths)], collapse = "; ")
  )
}

assert <- function(condition, message) {
  if (!isTRUE(condition)) {
    stop(message, call. = FALSE)
  }
}

coordinates <- fread(file.path(table_dir, "embedding_coordinates_5sets_long.csv"))
parameters <- fread(file.path(table_dir, "embedding_parameters_5sets.csv"))
membership <- fread(file.path(table_dir, "embedding_set_membership.csv"))
plot_data <- fread(file.path(table_dir, "pathway_plot_data_5sets.csv"))
summary_table <- fread(file.path(table_dir, "pathway_summary_5sets.csv"))
manifest <- fread(file.path(table_dir, "figure_manifest_15_grids.csv"))
recommendation <- fread(
  file.path(table_dir, "parameter_tuning_recommendation.csv")
)

expected_sets <- c(
  "all_507",
  "normal_tissue",
  "normal_cells",
  "cancer_tissue",
  "cancer_cells"
)
expected_counts <- c(507L, 183L, 471L, 178L, 383L)
observed_counts <- membership[
  ,
  .N,
  by = .(EmbeddingSet, SetOrder)
][order(SetOrder)]
assert(
  identical(observed_counts$EmbeddingSet, expected_sets) &&
    identical(observed_counts$N, expected_counts),
  "The five embedding sets are not 507/183/471/178/383 in requested order."
)

assert(
  nrow(coordinates) == sum(expected_counts) &&
    uniqueN(coordinates$EmbeddingSet) == 5L &&
    coordinates[
      ,
      .N,
      by = .(EmbeddingSet, SetOrder)
    ][order(SetOrder), identical(N, expected_counts)] &&
    all(is.finite(as.matrix(
      coordinates[
        ,
        .(UMAP_1, UMAP_2, TSNE_1, TSNE_2, PCA_1, PCA_2)
      ]
    ))),
  "Expected 1,722 finite rows across five independently fitted coordinate sets."
)

for (set_name in expected_sets[-1L]) {
  overlap <- merge(
    coordinates[
      EmbeddingSet == "all_507",
      .(
        BaseAccession,
        UMAP_1_all = UMAP_1,
        UMAP_2_all = UMAP_2,
        TSNE_1_all = TSNE_1,
        TSNE_2_all = TSNE_2,
        PCA_1_all = PCA_1,
        PCA_2_all = PCA_2
      )
    ],
    coordinates[
      EmbeddingSet == set_name,
      .(
        BaseAccession,
        UMAP_1_set = UMAP_1,
        UMAP_2_set = UMAP_2,
        TSNE_1_set = TSNE_1,
        TSNE_2_set = TSNE_2,
        PCA_1_set = PCA_1,
        PCA_2_set = PCA_2
      )
    ],
    by = "BaseAccession"
  )
  assert(
    nrow(overlap) == expected_counts[match(set_name, expected_sets)] &&
      any(abs(overlap$UMAP_1_all - overlap$UMAP_1_set) > 1e-10) &&
      any(abs(overlap$TSNE_1_all - overlap$TSNE_1_set) > 1e-10) &&
      any(abs(overlap$PCA_1_all - overlap$PCA_1_set) > 1e-10),
    paste0(set_name, " appears to reuse all-507 coordinates.")
  )
}

assert(
  nrow(plot_data) == sum(expected_counts) * 9L &&
    uniqueN(plot_data$EmbeddingSet) == 5L &&
    uniqueN(plot_data$DisplayLabel) == 9L &&
    plot_data[
      ,
      .N,
      by = .(EmbeddingSet, DisplayLabel)
    ][
      membership[, .(ProteinCount = .N), by = EmbeddingSet],
      on = "EmbeddingSet"
    ][, all(N == ProteinCount)],
  "Every set-pathway panel must contain exactly that set's proteins."
)
assert(
  summary_table[
    ,
    all(
      PromotingCount +
        SuppressingCount +
        UnassignedCount ==
        ProteinCount
    )
  ] &&
    plot_data[Score == 0L, all(Status == "Not assigned (0)")] &&
    plot_data[Score == 1L, all(Status == "Promoting (+1)")] &&
    plot_data[Score == -1L, all(Status == "Suppressing (-1)")],
  "Signed pathway status encoding or panel totals are inconsistent."
)

assert(
  recommendation[
    Method == "UMAP",
    NNeighbors
  ] == 12L &&
    recommendation[
      Method == "UMAP",
      MinDist
    ] == 0.35 &&
    recommendation[
      Method == "t-SNE",
      Perplexity
    ] == 25,
  "The selected UMAP/t-SNE tuning recommendation is incorrect."
)
assert(
  parameters[
    Method == "Shared" & Parameter == "OtherEmbeddingSetsUsedForCoordinates",
    all(Value == "FALSE")
  ] &&
    parameters[
      Method == "UMAP" & Parameter == "RandomSeed",
      all(Value == "25")
    ] &&
    parameters[
      Method == "t-SNE" & Parameter == "RandomSeed",
      all(Value == "25")
    ] &&
    parameters[
      Method == "PCA" & Parameter == "RandomSeedRecorded",
      all(Value == "25")
    ] &&
    parameters[
      Method == "PCA" & Parameter == "Solver",
      all(grepl("Gram matrix", Value, fixed = TRUE))
    ],
  "Independent-fit scope, seed 25, or PCA solver is not consistently recorded."
)

# The Gram-matrix PCA solver must reproduce ordinary prcomp sample scores
# exactly up to the arbitrary component signs on a subset where dgesdd works.
feature_table <- fread(
  file.path(
    project_root,
    "reanalysis/results/tables/kla_ddr_pathway_pie_umap_33groups_v4_bp_semantic",
    "protein_bp_semantic_shared_binary_matrix.csv"
  )
)
normal_tissue_ids <- membership[
  EmbeddingSet == "normal_tissue",
  BaseAccession
]
normal_tissue_matrix <- as.matrix(
  feature_table[
    match(normal_tissue_ids, BaseAccession),
    -"BaseAccession"
  ]
)
storage.mode(normal_tissue_matrix) <- "double"
normal_tissue_matrix <- normal_tissue_matrix /
  sqrt(rowSums(normal_tissue_matrix^2))
reference_pca <- prcomp(
  normal_tissue_matrix,
  center = TRUE,
  scale. = FALSE,
  rank. = 2
)
saved_pca <- coordinates[
  EmbeddingSet == "normal_tissue"
][match(normal_tissue_ids, BaseAccession)]
pca_correlations <- abs(diag(cor(
  reference_pca$x[, 1:2, drop = FALSE],
  as.matrix(saved_pca[, .(PCA_1, PCA_2)])
)))
assert(
  all(abs(pca_correlations - 1) < 1e-10),
  "Gram-matrix PCA scores are not equivalent to ordinary prcomp scores."
)

assert(
  nrow(manifest) == 15L &&
    uniqueN(manifest$EmbeddingSet) == 5L &&
    uniqueN(manifest$Method) == 3L &&
    all(manifest$IndependentFit == TRUE),
  "Expected exactly five independent sets x three methods = 15 figure entries."
)
figure_paths <- unlist(
  lapply(
    c("PNG", "PDF", "SVG"),
    function(column) file.path(project_root, manifest[[column]])
  ),
  use.names = FALSE
)
assert(
  length(figure_paths) == 45L &&
    all(file.exists(figure_paths)) &&
    all(file.info(figure_paths)$size > 0),
  "Expected 45 nonempty files for 15 figures in PNG/PDF/SVG formats."
)

cat(
  paste0(
    "PASS: 15 pathway grids use five independently fitted UMAP/t-SNE/PCA ",
    "coordinate sets, preserve 507/183/471/178/383 proteins, and record ",
    "the tuned seed-25 configuration.\n"
  )
)
