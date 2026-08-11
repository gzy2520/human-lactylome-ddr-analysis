#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(cluster)
  library(data.table)
})

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
script_path <- if (length(file_arg) > 0L) {
  normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE)
} else {
  normalizePath(
    "reanalysis/scripts/cluster_kla_ddr_pathway_scores_33groups.R",
    mustWork = TRUE
  )
}
project_root <- normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = TRUE)

analysis_name <- "kla_ddr_cytoscape_pathway_clusters_33groups_v1"
table_dir <- file.path(project_root, "reanalysis/results/tables", analysis_name)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

v3_dir <- file.path(
  project_root,
  "reanalysis/results/tables/kla_ddr_pathway_pie_umap_33groups_v3_all_go"
)
v4_dir <- file.path(
  project_root,
  "reanalysis/results/tables/kla_ddr_pathway_pie_umap_33groups_v4_bp_semantic"
)
feature_matrix_path <- file.path(
  v4_dir,
  "protein_bp_semantic_shared_binary_matrix.csv"
)
assignment_path <- file.path(v3_dir, "pathway_assignment_long_v3.csv")
metadata_path <- file.path(v3_dir, "pathway_umap_plot_data_v3_all_go.csv")
color_path <- file.path(v4_dir, "pathway_color_key_v4.csv")

required_inputs <- c(
  feature_matrix_path,
  assignment_path,
  metadata_path,
  color_path
)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs) > 0L) {
  stop("Missing required input(s): ", paste(missing_inputs, collapse = "; "))
}

assert <- function(condition, message) {
  if (!isTRUE(condition)) {
    stop(message, call. = FALSE)
  }
}

relative_path <- function(path) {
  sub(
    paste0("^", gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", project_root), "/?"),
    "",
    path
  )
}

adjusted_rand_index <- function(x, y) {
  tab <- table(x, y)
  choose2 <- function(z) z * (z - 1) / 2
  n <- sum(tab)
  if (n < 2L) {
    return(NA_real_)
  }
  sum_cells <- sum(choose2(tab))
  sum_rows <- sum(choose2(rowSums(tab)))
  sum_cols <- sum(choose2(colSums(tab)))
  total_pairs <- choose2(n)
  expected <- sum_rows * sum_cols / total_pairs
  maximum <- (sum_rows + sum_cols) / 2
  if (maximum == expected) {
    return(1)
  }
  (sum_cells - expected) / (maximum - expected)
}

pathway_labels <- c(
  "HR",
  "NHEJ",
  "AEJ",
  "BER",
  "NER",
  "MMR",
  "FA",
  "Chromatin interaction",
  "Other support"
)
score_column_names <- c(
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
random_seed <- 25L
nstart_final <- 250L
iter_max <- 1000L
stability_seeds <- 1:50

protein_ids <- fread(feature_matrix_path, select = "BaseAccession")$BaseAccession
assignment_long <- fread(assignment_path)
protein_metadata <- unique(
  fread(metadata_path)[, .(BaseAccession, GeneSymbol, ProteinName)]
)
pathway_colors <- fread(color_path)[order(PathwayOrder)]

assert(
  length(protein_ids) == 507L &&
    uniqueN(protein_ids) == 507L,
  "Expected exactly 507 unique BaseAccession identifiers."
)
assert(
  all(assignment_long$BaseAccession %in% protein_ids) &&
    all(assignment_long$DisplayLabel %in% pathway_labels) &&
    all(assignment_long$Score %in% c(-1L, 1L)),
  "Signed pathway assignments do not match the fixed 507-protein scope."
)
assert(
  !anyDuplicated(
    assignment_long[, .(BaseAccession, DisplayLabel)]
  ),
  "A protein-pathway pair has more than one signed assignment."
)
assert(
  nrow(protein_metadata) == 507L &&
    uniqueN(protein_metadata$BaseAccession) == 507L &&
    setequal(protein_metadata$BaseAccession, protein_ids),
  "Protein display metadata does not match the fixed 507-protein scope."
)
assert(
  identical(pathway_colors$DisplayLabel, pathway_labels),
  "The saved pathway order does not match the required 7/9 pathway definitions."
)

score_matrix <- matrix(
  0L,
  nrow = length(protein_ids),
  ncol = length(pathway_labels),
  dimnames = list(protein_ids, pathway_labels)
)
score_matrix[
  cbind(
    match(assignment_long$BaseAccession, protein_ids),
    match(assignment_long$DisplayLabel, pathway_labels)
  )
] <- assignment_long$Score
storage.mode(score_matrix) <- "double"

fit_one <- function(k) {
  x <- score_matrix[, seq_len(k), drop = FALSE]
  set.seed(random_seed)
  fit <- kmeans(
    x,
    centers = k,
    nstart = nstart_final,
    iter.max = iter_max,
    algorithm = "Hartigan-Wong"
  )

  cluster_summary <- rbindlist(lapply(seq_len(k), function(cluster_id) {
    member_index <- fit$cluster == cluster_id
    center <- fit$centers[cluster_id, ]
    ranked <- order(abs(center), decreasing = TRUE)
    selected <- ranked[abs(center[ranked]) >= 0.25]
    selected <- head(selected, 3L)
    signature <- if (length(selected) == 0L) {
      "Low-assignment pattern"
    } else {
      paste0(
        names(center)[selected],
        ifelse(center[selected] >= 0, "+", "\u2212"),
        collapse = "; "
      )
    }
    zero_fraction <- mean(rowSums(abs(x[member_index, , drop = FALSE])) == 0)
    if (zero_fraction >= 0.5) {
      signature <- paste0("Unassigned-dominant; ", signature)
    }
    data.table(
      RawCluster = cluster_id,
      ClusterSize = sum(member_index),
      AllZeroCount = sum(rowSums(abs(x[member_index, , drop = FALSE])) == 0),
      AllZeroFraction = zero_fraction,
      Signature = signature,
      PrimaryPathwayOrder = ranked[[1L]],
      PrimaryAbsoluteCenter = abs(center[ranked[[1L]]])
    )
  }))

  cluster_summary <- cluster_summary[
    order(
      Signature == "Low-assignment pattern",
      grepl("^Unassigned-dominant", Signature),
      PrimaryPathwayOrder,
      -PrimaryAbsoluteCenter,
      Signature
    )
  ]
  cluster_summary[, Cluster := seq_len(.N)]
  cluster_summary[
    ,
    ClusterCode := sprintf("K%d_C%d", k, Cluster)
  ]
  cluster_summary[
    ,
    ClusterLabel := paste0(ClusterCode, " | ", Signature)
  ]
  raw_to_canonical <- setNames(
    cluster_summary$Cluster,
    cluster_summary$RawCluster
  )
  canonical_cluster <- as.integer(
    raw_to_canonical[as.character(fit$cluster)]
  )

  centers <- as.data.table(fit$centers, keep.rownames = "RawCluster")
  centers[, RawCluster := as.integer(RawCluster)]
  centers <- merge(
    cluster_summary[
      ,
      .(
        RawCluster,
        Cluster,
        ClusterCode,
        ClusterLabel,
        ClusterSize,
        AllZeroCount,
        AllZeroFraction
      )
    ],
    centers,
    by = "RawCluster",
    sort = FALSE
  )
  setorder(centers, Cluster)

  sil <- silhouette(fit$cluster, dist(x, method = "euclidean"))
  stability_ari <- vapply(stability_seeds, function(seed_value) {
    set.seed(seed_value)
    alternative <- kmeans(
      x,
      centers = k,
      nstart = 25L,
      iter.max = iter_max,
      algorithm = "Hartigan-Wong"
    )
    adjusted_rand_index(fit$cluster, alternative$cluster)
  }, numeric(1))

  diagnostics <- data.table(
    DimensionCount = k,
    RequestedClusterCount = k,
    RandomSeed = random_seed,
    FinalNstart = nstart_final,
    IterMax = iter_max,
    ProteinCount = nrow(x),
    AllZeroProteinCount = sum(rowSums(abs(x)) == 0),
    UniqueScorePatternCount = uniqueN(as.data.table(x)),
    MeanSilhouette = mean(sil[, "sil_width"]),
    MedianSilhouette = median(sil[, "sil_width"]),
    WithinToTotalSS = fit$tot.withinss / fit$totss,
    StabilityMedianARI = median(stability_ari),
    StabilityMinARI = min(stability_ari),
    StabilityMaxARI = max(stability_ari)
  )

  pattern_frequency <- as.data.table(x)
  setnames(pattern_frequency, names(pattern_frequency), score_column_names[seq_len(k)])
  pattern_frequency <- pattern_frequency[, .N, by = names(pattern_frequency)]
  setorder(pattern_frequency, -N)
  setnames(pattern_frequency, "N", "ProteinCount")
  pattern_frequency[, PatternRank := seq_len(.N)]
  setcolorder(
    pattern_frequency,
    c("PatternRank", "ProteinCount", score_column_names[seq_len(k)])
  )

  list(
    x = x,
    cluster = canonical_cluster,
    summary = cluster_summary[order(Cluster)],
    centers = centers,
    diagnostics = diagnostics,
    pattern_frequency = pattern_frequency
  )
}

k7 <- fit_one(7L)
k9 <- fit_one(9L)

score_table <- as.data.table(score_matrix, keep.rownames = "BaseAccession")
setnames(score_table, pathway_labels, score_column_names)
score_table <- merge(
  data.table(BaseAccession = protein_ids),
  score_table,
  by = "BaseAccession",
  all.x = TRUE,
  sort = FALSE
)
score_table[, K7_AllZero := rowSums(abs(as.matrix(.SD))) == 0, .SDcols = score_column_names[1:7]]
score_table[, K9_AllZero := rowSums(abs(as.matrix(.SD))) == 0, .SDcols = score_column_names]
score_table[, K7_Cluster := k7$cluster]
score_table[, K9_Cluster := k9$cluster]

k7_labels <- setNames(k7$summary$ClusterLabel, k7$summary$Cluster)
k9_labels <- setNames(k9$summary$ClusterLabel, k9$summary$Cluster)
score_table[
  ,
  K7_ClusterLabel := unname(k7_labels[as.character(K7_Cluster)])
]
score_table[
  ,
  K9_ClusterLabel := unname(k9_labels[as.character(K9_Cluster)])
]
score_table <- merge(
  score_table,
  protein_metadata,
  by = "BaseAccession",
  all.x = TRUE,
  sort = FALSE
)
setcolorder(
  score_table,
  c(
    "BaseAccession",
    "GeneSymbol",
    "ProteinName",
    score_column_names,
    "K7_AllZero",
    "K7_Cluster",
    "K7_ClusterLabel",
    "K9_AllZero",
    "K9_Cluster",
    "K9_ClusterLabel"
  )
)

k7_matrix <- score_table[
  ,
  c("BaseAccession", "GeneSymbol", "ProteinName", score_column_names[1:7]),
  with = FALSE
]
k9_matrix <- score_table[
  ,
  c("BaseAccession", "GeneSymbol", "ProteinName", score_column_names),
  with = FALSE
]
cluster_assignments <- score_table[
  ,
  .(
    BaseAccession,
    GeneSymbol,
    ProteinName,
    K7_AllZero,
    K7_Cluster,
    K7_ClusterLabel,
    K9_AllZero,
    K9_Cluster,
    K9_ClusterLabel
  )
]

diagnostics <- rbind(k7$diagnostics, k9$diagnostics)
input_audit <- data.table(
  InputFile = relative_path(required_inputs),
  MD5 = unname(tools::md5sum(required_inputs)),
  Role = c(
    "Fixed 507 BaseAccession scope",
    "Signed manual pathway assignments",
    "Protein display metadata",
    "Pathway order and color key"
  )
)

fwrite(k7_matrix, file.path(table_dir, "pathway_score_matrix_7.csv"))
fwrite(k9_matrix, file.path(table_dir, "pathway_score_matrix_9.csv"))
fwrite(score_table, file.path(table_dir, "cytoscape_node_import_table.csv"))
fwrite(
  cluster_assignments,
  file.path(table_dir, "pathway_cluster_assignments_seed25.csv")
)
fwrite(k7$centers, file.path(table_dir, "pathway_cluster_centers_k7.csv"))
fwrite(k9$centers, file.path(table_dir, "pathway_cluster_centers_k9.csv"))
fwrite(k7$summary, file.path(table_dir, "pathway_cluster_summary_k7.csv"))
fwrite(k9$summary, file.path(table_dir, "pathway_cluster_summary_k9.csv"))
fwrite(
  k7$pattern_frequency,
  file.path(table_dir, "pathway_score_pattern_frequencies_k7.csv")
)
fwrite(
  k9$pattern_frequency,
  file.path(table_dir, "pathway_score_pattern_frequencies_k9.csv")
)
fwrite(diagnostics, file.path(table_dir, "pathway_cluster_diagnostics.csv"))
fwrite(input_audit, file.path(table_dir, "input_file_audit.csv"))

writeLines(
  c(
    capture.output(sessionInfo()),
    "",
    paste0("Random seed: ", random_seed),
    paste0("Final k-means nstart: ", nstart_final),
    paste0("k values: 7, 9")
  ),
  file.path(table_dir, "session_info.txt")
)

cat(
  "Wrote pathway-score clustering outputs to:\n",
  table_dir,
  "\n",
  sprintf(
    "k=7: all-zero=%d, mean silhouette=%.3f, cluster sizes=%s\n",
    k7$diagnostics$AllZeroProteinCount,
    k7$diagnostics$MeanSilhouette,
    paste(sort(k7$summary$ClusterSize, decreasing = TRUE), collapse = "/")
  ),
  sprintf(
    "k=9: all-zero=%d, mean silhouette=%.3f, cluster sizes=%s\n",
    k9$diagnostics$AllZeroProteinCount,
    k9$diagnostics$MeanSilhouette,
    paste(sort(k9$summary$ClusterSize, decreasing = TRUE), collapse = "/")
  ),
  sep = ""
)
