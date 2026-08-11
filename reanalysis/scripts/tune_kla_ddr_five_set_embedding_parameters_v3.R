#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(uwot)
  library(Rtsne)
})

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
script_path <- if (length(file_arg) > 0L) {
  normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE)
} else {
  normalizePath(
    "reanalysis/scripts/tune_kla_ddr_five_set_embedding_parameters_v3.R",
    mustWork = TRUE
  )
}
project_root <- normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = TRUE)

v4_table_dir <- file.path(
  project_root,
  "reanalysis/results/tables/kla_ddr_pathway_pie_umap_33groups_v4_bp_semantic"
)
venn_table_dir <- file.path(
  project_root,
  "reanalysis/results/tables/four_class_venn/kla_ddr_four_class_venn"
)
feature_matrix_path <- file.path(
  v4_table_dir,
  "protein_bp_semantic_shared_binary_matrix.csv"
)
membership_path <- file.path(venn_table_dir, "membership.csv")

analysis_name <- "kla_ddr_five_set_three_embedding_pathway_grids_33groups_v3"
table_dir <- file.path(project_root, "reanalysis/results/tables", analysis_name)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

assert <- function(condition, message) {
  if (!isTRUE(condition)) {
    stop(message, call. = FALSE)
  }
}

feature_table <- fread(feature_matrix_path)
feature_accessions <- feature_table$BaseAccession
feature_matrix <- as.matrix(feature_table[, -"BaseAccession"])
storage.mode(feature_matrix) <- "double"
rownames(feature_matrix) <- feature_accessions
membership <- fread(membership_path)

assert(
  identical(dim(feature_matrix), c(507L, 3008L)) &&
    uniqueN(feature_accessions) == 507L &&
    all(feature_matrix %in% c(0, 1)),
  "Expected the fixed 507 x 3,008 binary BP semantic matrix."
)

embedding_sets <- list(
  all_507 = feature_accessions,
  normal_tissue = membership[In_normal_tissue == TRUE, BaseAccession],
  normal_cells = membership[In_normal_cells == TRUE, BaseAccession],
  cancer_tissue = membership[In_cancer_tissue == TRUE, BaseAccession],
  cancer_cells = membership[In_cancer_cells == TRUE, BaseAccession]
)
expected_counts <- c(507L, 183L, 471L, 178L, 383L)
assert(
  identical(as.integer(lengths(embedding_sets)), expected_counts),
  "The five embedding-set sizes are not 507/183/471/178/383."
)

random_seed <- 25L

cosine_distance_matrix <- function(x) {
  row_norms <- sqrt(rowSums(x^2))
  assert(all(row_norms > 0), "A protein has no retained BP semantic feature.")
  normalized <- x / row_norms
  similarity <- tcrossprod(normalized)
  similarity[similarity > 1] <- 1
  similarity[similarity < -1] <- -1
  distance <- 1 - similarity
  diag(distance) <- 0
  distance
}

gram_pca_scores <- function(x) {
  centered <- sweep(x, 2, colMeans(x), "-")
  gram <- tcrossprod(centered)
  eigendecomposition <- eigen(gram, symmetric = TRUE)
  eigenvalues <- pmax(eigendecomposition$values, 0)
  scores <- sweep(
    eigendecomposition$vectors[, 1:2, drop = FALSE],
    2,
    sqrt(eigenvalues[1:2]),
    "*"
  )
  for (component_index in 1:2) {
    anchor_index <- which.max(abs(scores[, component_index]))
    if (scores[anchor_index, component_index] < 0) {
      scores[, component_index] <- -scores[, component_index]
    }
  }
  list(
    Scores = scores,
    VariancePercent = 100 * eigenvalues[1:2] / sum(diag(gram))
  )
}

pca_preprocess <- function(x, preprocessing) {
  n <- nrow(x)
  keep <- rep(TRUE, ncol(x))
  transformed <- x
  if (grepl("subset_filter", preprocessing, fixed = TRUE)) {
    prevalence <- colSums(x)
    keep <- prevalence >= 2L & prevalence <= floor(0.90 * n)
    transformed <- transformed[, keep, drop = FALSE]
  }
  if (grepl("idf", preprocessing, fixed = TRUE)) {
    idf <- log((1 + n) / (1 + colSums(transformed))) + 1
    transformed <- sweep(transformed, 2, idf, "*")
  }
  if (grepl("l2", preprocessing, fixed = TRUE)) {
    transformed <- transformed / sqrt(rowSums(transformed^2))
  }
  list(Matrix = transformed, FeatureCount = sum(keep))
}

embedding_quality <- function(high_distance, coordinates, k = 10L) {
  n <- nrow(high_distance)
  k <- min(k, n - 1L)
  low_distance <- as.matrix(dist(coordinates))
  high_for_knn <- high_distance
  low_for_knn <- low_distance
  diag(high_for_knn) <- Inf
  diag(low_for_knn) <- Inf
  high_knn <- t(vapply(
    seq_len(n),
    function(i) order(high_for_knn[i, ])[seq_len(k)],
    integer(k)
  ))
  low_knn <- t(vapply(
    seq_len(n),
    function(i) order(low_for_knn[i, ])[seq_len(k)],
    integer(k)
  ))
  overlap <- mean(vapply(
    seq_len(n),
    function(i) length(intersect(high_knn[i, ], low_knn[i, ])) / k,
    numeric(1)
  ))

  coordinate_ranges <- apply(coordinates, 2, function(x) diff(range(x)))
  coordinate_ranges[coordinate_ranges == 0] <- 1
  normalized <- sweep(
    sweep(coordinates, 2, apply(coordinates, 2, min), "-"),
    2,
    coordinate_ranges,
    "/"
  )
  normalized_distance <- as.matrix(dist(normalized))
  diag(normalized_distance) <- Inf
  nearest_neighbor <- apply(normalized_distance, 1, min)
  upper <- upper.tri(high_distance)

  robust_ranges <- apply(
    coordinates,
    2,
    function(x) diff(quantile(x, c(0.025, 0.975), names = FALSE))
  )
  full_ranges <- apply(coordinates, 2, function(x) diff(range(x)))
  full_ranges[full_ranges == 0] <- 1
  robust_ranges[robust_ranges == 0] <- 1e-12

  data.table(
    NeighborOverlapK10 = overlap,
    GlobalDistanceSpearman = suppressWarnings(cor(
      high_distance[upper],
      low_distance[upper],
      method = "spearman"
    )),
    NormalizedNearestNeighborP05 = unname(quantile(
      nearest_neighbor,
      probs = 0.05,
      names = FALSE
    )),
    NormalizedNearestNeighborMedian = median(nearest_neighbor),
    CoreAreaFraction = prod(pmin(robust_ranges / full_ranges, 1)),
    RobustIsotropy = min(robust_ranges) / max(robust_ranges)
  )
}

safe_z <- function(x) {
  standard_deviation <- sd(x)
  if (!is.finite(standard_deviation) || standard_deviation == 0) {
    return(rep(0, length(x)))
  }
  (x - mean(x)) / standard_deviation
}

set_inputs <- lapply(names(embedding_sets), function(set_name) {
  accessions <- feature_accessions[
    feature_accessions %in% embedding_sets[[set_name]]
  ]
  x <- feature_matrix[match(accessions, feature_accessions), , drop = FALSE]
  list(
    Set = set_name,
    X = x,
    Distance = cosine_distance_matrix(x)
  )
})
names(set_inputs) <- names(embedding_sets)

metric_rows <- list()

umap_candidates <- CJ(
  NNeighbors = c(8L, 12L, 18L, 24L),
  MinDist = c(0.25, 0.40, 0.55, 0.70),
  RepulsionStrength = c(1.5, 2.0),
  unique = TRUE
)
umap_candidates[, Candidate := sprintf(
  "umap_nn%02d_md%02d_r%02d",
  NNeighbors,
  as.integer(round(MinDist * 100)),
  as.integer(round(RepulsionStrength * 10))
)]

for (candidate_index in seq_len(nrow(umap_candidates))) {
  candidate <- umap_candidates[candidate_index]
  for (set_name in names(set_inputs)) {
    input <- set_inputs[[set_name]]
    set.seed(random_seed)
    coordinates <- uwot::umap(
      X = input$X,
      n_neighbors = candidate$NNeighbors,
      n_components = 2L,
      metric = "cosine",
      n_epochs = 700L,
      scale = FALSE,
      init = "random",
      spread = 2.0,
      min_dist = candidate$MinDist,
      repulsion_strength = candidate$RepulsionStrength,
      negative_sample_rate = 10L,
      fast_sgd = FALSE,
      n_threads = 1L,
      n_sgd_threads = 1L,
      seed = random_seed,
      verbose = FALSE
    )
    metric_rows[[length(metric_rows) + 1L]] <- cbind(
      data.table(
        Method = "UMAP",
        Candidate = candidate$Candidate,
        EmbeddingSet = set_name,
        ProteinCount = nrow(input$X),
        NNeighbors = candidate$NNeighbors,
        MinDist = candidate$MinDist,
        Spread = 2.0,
        RepulsionStrength = candidate$RepulsionStrength,
        NegativeSampleRate = 10L,
        Iterations = 700L,
        Perplexity = NA_real_,
        Theta = NA_real_,
        Preprocessing = NA_character_,
        FeatureCount = ncol(input$X)
      ),
      embedding_quality(input$Distance, coordinates)
    )
  }
}

tsne_candidates <- CJ(
  Perplexity = c(10, 15, 20, 25, 30, 40, 50),
  Theta = c(0.3, 0.5),
  unique = TRUE
)
tsne_candidates[, Candidate := sprintf(
  "tsne_p%02d_t%02d",
  Perplexity,
  as.integer(round(Theta * 10))
)]

for (candidate_index in seq_len(nrow(tsne_candidates))) {
  candidate <- tsne_candidates[candidate_index]
  for (set_name in names(set_inputs)) {
    input <- set_inputs[[set_name]]
    set.seed(random_seed)
    fit <- Rtsne(
      as.dist(input$Distance),
      dims = 2L,
      perplexity = candidate$Perplexity,
      theta = candidate$Theta,
      max_iter = 1000L,
      check_duplicates = FALSE,
      pca = FALSE,
      normalize = FALSE,
      num_threads = 1L,
      verbose = FALSE
    )
    metric_rows[[length(metric_rows) + 1L]] <- cbind(
      data.table(
        Method = "t-SNE",
        Candidate = candidate$Candidate,
        EmbeddingSet = set_name,
        ProteinCount = nrow(input$X),
        NNeighbors = NA_integer_,
        MinDist = NA_real_,
        Spread = NA_real_,
        RepulsionStrength = NA_real_,
        NegativeSampleRate = NA_integer_,
        Iterations = 1000L,
        Perplexity = candidate$Perplexity,
        Theta = candidate$Theta,
        Preprocessing = NA_character_,
        FeatureCount = ncol(input$X)
      ),
      embedding_quality(input$Distance, fit$Y)
    )
  }
}

pca_candidates <- c(
  "binary_unscaled",
  "l2_all_features",
  "l2_subset_filter_02_90",
  "idf_l2_all_features",
  "idf_l2_subset_filter_02_90"
)
for (preprocessing in pca_candidates) {
  for (set_name in names(set_inputs)) {
    input <- set_inputs[[set_name]]
    prepared <- pca_preprocess(input$X, preprocessing)
    pca_result <- gram_pca_scores(prepared$Matrix)
    metric_rows[[length(metric_rows) + 1L]] <- cbind(
      data.table(
        Method = "PCA",
        Candidate = paste0("pca_", preprocessing),
        EmbeddingSet = set_name,
        ProteinCount = nrow(input$X),
        NNeighbors = NA_integer_,
        MinDist = NA_real_,
        Spread = NA_real_,
        RepulsionStrength = NA_real_,
        NegativeSampleRate = NA_integer_,
        Iterations = NA_integer_,
        Perplexity = NA_real_,
        Theta = NA_real_,
        Preprocessing = preprocessing,
        FeatureCount = prepared$FeatureCount,
        PC1VariancePercent = pca_result$VariancePercent[[1L]],
        PC2VariancePercent = pca_result$VariancePercent[[2L]]
      ),
      embedding_quality(input$Distance, pca_result$Scores)
    )
  }
}

metrics <- rbindlist(metric_rows, use.names = TRUE, fill = TRUE)
metrics[
  ,
  `:=`(
    ZNeighborOverlap = safe_z(NeighborOverlapK10),
    ZGlobalDistance = safe_z(GlobalDistanceSpearman),
    ZCrowdingP05 = safe_z(log(NormalizedNearestNeighborP05 + 1e-6)),
    ZCoreArea = safe_z(CoreAreaFraction),
    ZIsotropy = safe_z(RobustIsotropy)
  ),
  by = .(Method, EmbeddingSet)
]
metrics[
  Method %in% c("UMAP", "t-SNE"),
  CompositeScore := 0.50 * ZNeighborOverlap +
    0.15 * ZGlobalDistance +
    0.20 * ZCrowdingP05 +
    0.10 * ZCoreArea +
    0.05 * ZIsotropy
]
metrics[
  Method == "PCA",
  CompositeScore := 0.30 * ZNeighborOverlap +
    0.40 * ZGlobalDistance +
    0.15 * ZCrowdingP05 +
    0.10 * ZCoreArea +
    0.05 * ZIsotropy
]

recommendation <- metrics[
  ,
  .SD[which.max(CompositeScore)],
  by = .(EmbeddingSet, Method)
][order(match(EmbeddingSet, names(embedding_sets)), match(Method, c(
  "UMAP",
  "t-SNE",
  "PCA"
)))]
recommendation[Method == "UMAP", Iterations := 1200L]
recommendation[
  ,
  SelectionRule := fcase(
    Method %in% c("UMAP", "t-SNE"),
    paste0(
      "Per-set maximum of 0.50*z(kNN overlap) + 0.15*z(global Spearman) + ",
      "0.20*z(nearest-neighbor P05) + 0.10*z(core area) + ",
      "0.05*z(isotropy)."
    ),
    default = paste0(
      "Per-set maximum of 0.30*z(kNN overlap) + 0.40*z(global Spearman) + ",
      "0.15*z(nearest-neighbor P05) + 0.10*z(core area) + ",
      "0.05*z(isotropy)."
    )
  )
]

candidate_summary <- metrics[
  ,
  .(
    RankWithinSetMethod = frank(-CompositeScore, ties.method = "first"),
    Candidate,
    NNeighbors,
    MinDist,
    Spread,
    RepulsionStrength,
    NegativeSampleRate,
    Iterations,
    Perplexity,
    Theta,
    Preprocessing,
    FeatureCount,
    PC1VariancePercent,
    PC2VariancePercent,
    NeighborOverlapK10,
    GlobalDistanceSpearman,
    NormalizedNearestNeighborP05,
    NormalizedNearestNeighborMedian,
    CoreAreaFraction,
    RobustIsotropy,
    CompositeScore
  ),
  by = .(EmbeddingSet, ProteinCount, Method)
][order(match(EmbeddingSet, names(embedding_sets)), Method, RankWithinSetMethod)]

fwrite(metrics, file.path(table_dir, "parameter_tuning_metrics_by_set.csv"))
fwrite(
  candidate_summary,
  file.path(table_dir, "parameter_tuning_ranked_candidates_by_set.csv")
)
fwrite(
  recommendation,
  file.path(table_dir, "parameter_tuning_recommendation_by_set.csv")
)
writeLines(capture.output(sessionInfo()), file.path(
  table_dir,
  "parameter_tuning_session_info.txt"
), useBytes = TRUE)

message(
  "Per-set recommendation: ",
  file.path(table_dir, "parameter_tuning_recommendation_by_set.csv")
)
print(recommendation[
  ,
  .(
    EmbeddingSet,
    Method,
    Candidate,
    NNeighbors,
    MinDist,
    RepulsionStrength,
    Perplexity,
    Theta,
    Preprocessing,
    FeatureCount,
    CompositeScore
  )
])
