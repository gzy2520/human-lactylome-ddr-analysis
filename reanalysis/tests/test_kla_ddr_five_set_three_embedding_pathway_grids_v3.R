#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

project_root <- normalizePath(".", mustWork = TRUE)
analysis_name <- "kla_ddr_five_set_three_embedding_pathway_grids_33groups_v3"
table_dir <- file.path(project_root, "reanalysis/results/tables", analysis_name)

required_tables <- c(
  "embedding_coordinates_5sets_tuned_long.csv",
  "embedding_parameters_5sets_tuned.csv",
  "embedding_set_membership.csv",
  "pathway_plot_data_5sets_tuned.csv",
  "pathway_summary_5sets.csv",
  "figure_manifest_15_grids.csv",
  "parameter_tuning_metrics_by_set.csv",
  "parameter_tuning_ranked_candidates_by_set.csv",
  "parameter_tuning_recommendation_by_set.csv",
  "input_file_audit.csv",
  "session_info.txt",
  "parameter_tuning_session_info.txt"
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

coordinates <- fread(
  file.path(table_dir, "embedding_coordinates_5sets_tuned_long.csv")
)
parameters <- fread(file.path(table_dir, "embedding_parameters_5sets_tuned.csv"))
membership <- fread(file.path(table_dir, "embedding_set_membership.csv"))
plot_data <- fread(file.path(table_dir, "pathway_plot_data_5sets_tuned.csv"))
summary_table <- fread(file.path(table_dir, "pathway_summary_5sets.csv"))
manifest <- fread(file.path(table_dir, "figure_manifest_15_grids.csv"))
metrics <- fread(file.path(table_dir, "parameter_tuning_metrics_by_set.csv"))
ranked <- fread(
  file.path(table_dir, "parameter_tuning_ranked_candidates_by_set.csv")
)
recommendation <- fread(
  file.path(table_dir, "parameter_tuning_recommendation_by_set.csv")
)

expected_sets <- c(
  "all_507",
  "normal_tissue",
  "normal_cells",
  "cancer_tissue",
  "cancer_cells"
)
expected_counts <- c(507L, 183L, 471L, 178L, 383L)
expected_methods <- c("UMAP", "t-SNE", "PCA")

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
  uniqueN(membership, by = c("EmbeddingSet", "BaseAccession")) ==
    nrow(membership),
  "An embedding set contains duplicated BaseAccession analysis keys."
)

coordinate_counts <- coordinates[
  ,
  .N,
  by = .(EmbeddingSet, SetOrder)
][order(SetOrder)]
assert(
  nrow(coordinates) == sum(expected_counts) &&
    identical(coordinate_counts$EmbeddingSet, expected_sets) &&
    identical(coordinate_counts$N, expected_counts) &&
    all(is.finite(as.matrix(
      coordinates[
        ,
        .(UMAP_1, UMAP_2, TSNE_1, TSNE_2, PCA_1, PCA_2)
      ]
    ))),
  "Expected 1,722 finite coordinate rows across the five independent fits."
)

# Each category must be recomputed from its own proteins, not copied from all_507.
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
  expected_n <- expected_counts[match(set_name, expected_sets)]
  assert(
    nrow(overlap) == expected_n &&
      any(abs(overlap$UMAP_1_all - overlap$UMAP_1_set) > 1e-10) &&
      any(abs(overlap$TSNE_1_all - overlap$TSNE_1_set) > 1e-10) &&
      any(abs(overlap$PCA_1_all - overlap$PCA_1_set) > 1e-10),
    paste0(set_name, " appears to reuse the all_507 coordinates.")
  )
}

# Tuning covered the documented per-set candidate grid.
candidate_counts <- metrics[
  ,
  .(CandidateCount = uniqueN(Candidate)),
  by = .(EmbeddingSet, Method)
]
expected_candidate_counts <- data.table(
  Method = expected_methods,
  CandidateCount = c(32L, 14L, 5L)
)
candidate_counts <- merge(
  candidate_counts,
  expected_candidate_counts,
  by = "Method",
  suffixes = c("Observed", "Expected")
)
assert(
  nrow(candidate_counts) == 15L &&
    all(
      candidate_counts$CandidateCountObserved ==
        candidate_counts$CandidateCountExpected
    ) &&
    all(is.finite(metrics$CompositeScore)),
  "The V3 tuning grid is incomplete or contains non-finite scores."
)
assert(
  nrow(recommendation) == 15L &&
    all(recommendation[
      ,
      .N,
      by = .(EmbeddingSet, Method)
    ]$N == 1L) &&
    all(recommendation$ProteinCount == expected_counts[
      match(recommendation$EmbeddingSet, expected_sets)
    ]),
  "Expected exactly one recommendation for every set-method pair."
)
rank_one <- ranked[RankWithinSetMethod == 1L]
rank_one <- merge(
  rank_one[, .(EmbeddingSet, Method, RankedCandidate = Candidate)],
  recommendation[, .(EmbeddingSet, Method, RecommendedCandidate = Candidate)],
  by = c("EmbeddingSet", "Method")
)
assert(
  nrow(rank_one) == 15L &&
    all(rank_one$RankedCandidate == rank_one$RecommendedCandidate),
  "A formal recommendation does not match the top-ranked candidate."
)

expected_recommendations <- data.table(
  EmbeddingSet = rep(expected_sets, each = 3L),
  Method = rep(expected_methods, times = 5L),
  Candidate = c(
    "umap_nn18_md25_r15", "tsne_p20_t03", "pca_l2_all_features",
    "umap_nn08_md55_r20", "tsne_p40_t05", "pca_l2_all_features",
    "umap_nn18_md25_r15", "tsne_p25_t03",
    "pca_l2_subset_filter_02_90",
    "umap_nn12_md55_r15", "tsne_p30_t05",
    "pca_l2_subset_filter_02_90",
    "umap_nn12_md70_r20", "tsne_p25_t03",
    "pca_l2_subset_filter_02_90"
  ),
  FeatureCount = c(
    3008L, 3008L, 3008L,
    3008L, 3008L, 3008L,
    3008L, 3008L, 2968L,
    3008L, 3008L, 1688L,
    3008L, 3008L, 2418L
  )
)
recommendation_check <- merge(
  recommendation[
    ,
    .(EmbeddingSet, Method, Candidate, FeatureCount)
  ],
  expected_recommendations,
  by = c("EmbeddingSet", "Method"),
  suffixes = c("Observed", "Expected")
)
assert(
  nrow(recommendation_check) == 15L &&
    all(
      recommendation_check$CandidateObserved ==
        recommendation_check$CandidateExpected
    ) &&
    all(
      recommendation_check$FeatureCountObserved ==
        recommendation_check$FeatureCountExpected
    ),
  "The selected candidates or retained PCA feature counts have changed."
)

# Formal coordinates and long-form parameter records must use the recommendation.
coordinate_candidates <- unique(
  coordinates[
    ,
    .(
      EmbeddingSet,
      UMAPCandidate,
      TSNECandidate,
      PCAPreprocessing,
      PCAFeatureCount
    )
  ]
)
coordinate_candidates <- merge(
  coordinate_candidates,
  dcast(
    recommendation,
    EmbeddingSet ~ Method,
    value.var = "Candidate"
  ),
  by = "EmbeddingSet"
)
assert(
  nrow(coordinate_candidates) == 5L &&
    all(coordinate_candidates$UMAPCandidate == coordinate_candidates$UMAP) &&
    all(coordinate_candidates$TSNECandidate == coordinate_candidates$`t-SNE`) &&
    all(
      paste0("pca_", coordinate_candidates$PCAPreprocessing) ==
        coordinate_candidates$PCA
    ),
  "Saved coordinates do not use the formal per-set recommendations."
)
assert(
  parameters[
    Method == "Shared" & Parameter == "OtherEmbeddingSetsUsedForCoordinates",
    all(Value == "FALSE")
  ] &&
    parameters[
      Method == "Shared" & Parameter == "RandomSeedRecorded",
      all(Value == "25")
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
  "Independent-fit scope, seed 25, or the PCA solver is not recorded."
)

# Every pathway panel must contain the complete set and retain the signed encoding.
assert(
  nrow(plot_data) == sum(expected_counts) * 9L &&
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
    plot_data[Score == -1L, all(Status == "Suppressing (-1)")] &&
    all(plot_data$Score %in% c(-1L, 0L, 1L)),
  "Signed pathway status encoding or pathway totals are inconsistent."
)

# The Gram-matrix solver must reproduce prcomp scores after the selected
# cancer-tissue feature filter and row-wise L2 normalization.
feature_table <- fread(
  file.path(
    project_root,
    "reanalysis/results/tables/kla_ddr_pathway_pie_umap_33groups_v4_bp_semantic",
    "protein_bp_semantic_shared_binary_matrix.csv"
  )
)
cancer_tissue_ids <- membership[
  EmbeddingSet == "cancer_tissue",
  BaseAccession
]
cancer_tissue_matrix <- as.matrix(
  feature_table[
    match(cancer_tissue_ids, BaseAccession),
    -"BaseAccession"
  ]
)
storage.mode(cancer_tissue_matrix) <- "double"
prevalence <- colSums(cancer_tissue_matrix)
keep <- prevalence >= 2L &
  prevalence <= floor(0.90 * nrow(cancer_tissue_matrix))
cancer_tissue_matrix <- cancer_tissue_matrix[, keep, drop = FALSE]
cancer_tissue_matrix <- cancer_tissue_matrix /
  sqrt(rowSums(cancer_tissue_matrix^2))
assert(
  ncol(cancer_tissue_matrix) == 1688L,
  "The cancer-tissue PCA feature filter no longer retains 1,688 features."
)
reference_pca <- prcomp(
  cancer_tissue_matrix,
  center = TRUE,
  scale. = FALSE,
  rank. = 2
)
saved_pca <- coordinates[
  EmbeddingSet == "cancer_tissue"
][match(cancer_tissue_ids, BaseAccession)]
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
    all(manifest$IndependentFit == TRUE) &&
    all(manifest$TunedPerSet == TRUE),
  "Expected exactly five independently tuned sets x three methods."
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
    "PASS: V3 preserves 507/183/471/178/383 proteins, uses 15 per-set ",
    "tuned UMAP/t-SNE/PCA configurations at seed 25, and exports 45 ",
    "nonempty figure files with valid signed pathway states.\n"
  )
)
