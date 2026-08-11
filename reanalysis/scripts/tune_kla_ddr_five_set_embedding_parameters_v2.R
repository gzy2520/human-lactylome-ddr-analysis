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
    "reanalysis/scripts/tune_kla_ddr_five_set_embedding_parameters_v2.R",
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

analysis_name <- "kla_ddr_five_set_three_embedding_pathway_grids_33groups_v2"
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
  similarity <- tcrossprod(x) / outer(row_norms, row_norms)
  similarity[similarity > 1] <- 1
  similarity[similarity < -1] <- -1
  distance <- 1 - similarity
  diag(distance) <- 0
  distance
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
  global_spearman <- suppressWarnings(cor(
    high_distance[upper],
    low_distance[upper],
    method = "spearman"
  ))

  data.table(
    NeighborOverlapK10 = overlap,
    GlobalDistanceSpearman = global_spearman,
    NormalizedNearestNeighborP05 = unname(quantile(
      nearest_neighbor,
      probs = 0.05,
      names = FALSE
    )),
    NormalizedNearestNeighborMedian = median(nearest_neighbor)
  )
}

set_inputs <- lapply(names(embedding_sets), function(set_name) {
  accessions <- embedding_sets[[set_name]]
  indices <- match(accessions, feature_accessions)
  x <- feature_matrix[indices, , drop = FALSE]
  list(
    Set = set_name,
    X = x,
    Distance = cosine_distance_matrix(x)
  )
})
names(set_inputs) <- names(embedding_sets)

umap_candidates <- CJ(
  NNeighbors = c(6L, 8L, 12L, 16L),
  MinDist = c(0.20, 0.35, 0.50),
  unique = TRUE
)
umap_candidates[, Candidate := sprintf(
  "umap_nn%02d_md%02d",
  NNeighbors,
  as.integer(round(MinDist * 100))
)]

umap_results <- list()
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
      n_epochs = 600L,
      scale = FALSE,
      init = "random",
      spread = 1.8,
      min_dist = candidate$MinDist,
      repulsion_strength = 1.5,
      negative_sample_rate = 10L,
      fast_sgd = FALSE,
      n_threads = 1L,
      n_sgd_threads = 1L,
      seed = random_seed,
      verbose = FALSE
    )
    quality <- embedding_quality(input$Distance, coordinates)
    umap_results[[length(umap_results) + 1L]] <- cbind(
      data.table(
        Method = "UMAP",
        Candidate = candidate$Candidate,
        EmbeddingSet = set_name,
        ProteinCount = nrow(input$X),
        NNeighbors = candidate$NNeighbors,
        MinDist = candidate$MinDist,
        Spread = 1.8,
        RepulsionStrength = 1.5,
        NegativeSampleRate = 10L,
        Iterations = 600L,
        Perplexity = NA_real_,
        Theta = NA_real_
      ),
      quality
    )
  }
}

tsne_candidates <- data.table(
  Perplexity = c(15, 25, 35, 45),
  Theta = 0.5
)
tsne_candidates[, Candidate := sprintf("tsne_p%02d", Perplexity)]

tsne_results <- list()
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
    quality <- embedding_quality(input$Distance, fit$Y)
    tsne_results[[length(tsne_results) + 1L]] <- cbind(
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
        Theta = candidate$Theta
      ),
      quality
    )
  }
}

metrics <- rbindlist(c(umap_results, tsne_results), use.names = TRUE)
metrics[
  ,
  `:=`(
    ZNeighborOverlap = as.numeric(scale(NeighborOverlapK10)),
    ZGlobalDistance = as.numeric(scale(GlobalDistanceSpearman)),
    ZCrowdingP05 = as.numeric(scale(log(
      NormalizedNearestNeighborP05 + 1e-6
    )))
  ),
  by = .(Method, EmbeddingSet)
]
metrics[
  ,
  CompositeScore := 0.70 * ZNeighborOverlap +
    0.20 * ZGlobalDistance +
    0.10 * ZCrowdingP05
]

candidate_summary <- metrics[
  ,
  .(
    MeanCompositeScore = mean(CompositeScore),
    MinCompositeScore = min(CompositeScore),
    MeanNeighborOverlapK10 = mean(NeighborOverlapK10),
    MeanGlobalDistanceSpearman = mean(GlobalDistanceSpearman),
    MeanNearestNeighborP05 = mean(NormalizedNearestNeighborP05)
  ),
  by = .(
    Method,
    Candidate,
    NNeighbors,
    MinDist,
    Spread,
    RepulsionStrength,
    NegativeSampleRate,
    Iterations,
    Perplexity,
    Theta
  )
][order(Method, -MeanCompositeScore)]

recommendation <- candidate_summary[
  ,
  .SD[which.max(MeanCompositeScore)],
  by = Method
]
recommendation[
  ,
  SelectionRule := paste0(
    "Highest mean score across five independently embedded sets; ",
    "score = 0.70*z(kNN overlap) + 0.20*z(global Spearman) + ",
    "0.10*z(5th-percentile normalized nearest-neighbor distance)."
  )
]
recommendation[
  Method == "UMAP",
  Iterations := 1000L
]

fwrite(metrics, file.path(table_dir, "parameter_tuning_metrics.csv"))
fwrite(candidate_summary, file.path(table_dir, "parameter_tuning_candidate_summary.csv"))
fwrite(recommendation, file.path(table_dir, "parameter_tuning_recommendation.csv"))

message(
  "Tuning metrics: ",
  file.path(table_dir, "parameter_tuning_metrics.csv")
)
message(
  "Recommendation: ",
  file.path(table_dir, "parameter_tuning_recommendation.csv")
)
print(recommendation)
