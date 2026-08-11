#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(uwot)
  library(Rtsne)
})

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
script_path <- if (length(file_arg) > 0L) {
  normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE)
} else {
  normalizePath(
    "reanalysis/scripts/plot_kla_ddr_five_set_three_embedding_pathway_grids_v3.R",
    mustWork = TRUE
  )
}
project_root <- normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = TRUE)

v3_source_table_dir <- file.path(
  project_root,
  "reanalysis/results/tables/kla_ddr_pathway_pie_umap_33groups_v3_all_go"
)
v4_source_table_dir <- file.path(
  project_root,
  "reanalysis/results/tables/kla_ddr_pathway_pie_umap_33groups_v4_bp_semantic"
)
venn_table_dir <- file.path(
  project_root,
  "reanalysis/results/tables/four_class_venn/kla_ddr_four_class_venn"
)

feature_matrix_path <- file.path(
  v4_source_table_dir,
  "protein_bp_semantic_shared_binary_matrix.csv"
)
membership_path <- file.path(venn_table_dir, "membership.csv")
set_count_path <- file.path(venn_table_dir, "set_counts.csv")
assignment_path <- file.path(
  v3_source_table_dir,
  "pathway_assignment_long_v3.csv"
)
metadata_path <- file.path(
  v3_source_table_dir,
  "pathway_umap_plot_data_v3_all_go.csv"
)
color_path <- file.path(v4_source_table_dir, "pathway_color_key_v4.csv")

analysis_name <- "kla_ddr_five_set_three_embedding_pathway_grids_33groups_v3"
table_dir <- file.path(project_root, "reanalysis/results/tables", analysis_name)
figure_dir <- file.path(project_root, "reanalysis/results/figures", analysis_name)
report_path <- file.path(
  project_root,
  "reanalysis/reports/FIVE_SET_UMAP_TSNE_PCA_PATHWAY_GRIDS_33GROUP_V3.md"
)
recommendation_path <- file.path(
  table_dir,
  "parameter_tuning_recommendation_by_set.csv"
)
tuning_metrics_path <- file.path(
  table_dir,
  "parameter_tuning_metrics_by_set.csv"
)
tuning_ranked_path <- file.path(
  table_dir,
  "parameter_tuning_ranked_candidates_by_set.csv"
)

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

required_inputs <- c(
  feature_matrix_path,
  membership_path,
  set_count_path,
  assignment_path,
  metadata_path,
  color_path,
  recommendation_path,
  tuning_metrics_path,
  tuning_ranked_path
)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs) > 0L) {
  stop(
    "Missing required input(s): ",
    paste(missing_inputs, collapse = "; "),
    "\nRun reanalysis/scripts/tune_kla_ddr_five_set_embedding_parameters_v3.R first."
  )
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

random_seed <- 25L
set_info <- data.table(
  EmbeddingSet = c(
    "all_507",
    "normal_tissue",
    "normal_cells",
    "cancer_tissue",
    "cancer_cells"
  ),
  MembershipColumn = c(
    NA_character_,
    "In_normal_tissue",
    "In_normal_cells",
    "In_cancer_tissue",
    "In_cancer_cells"
  ),
  SetLabel = c(
    "All Kla-DDR proteins",
    "Normal/non-tumor tissues",
    "Normal/non-tumor cells",
    "Cancer tissues",
    "Cancer cells"
  ),
  SetLabelZh = c(
    "全部Kla∩DDR蛋白",
    "正常/非肿瘤组织",
    "正常/非肿瘤细胞",
    "癌症组织",
    "癌症细胞"
  ),
  SetOrder = 1:5
)

feature_table <- fread(feature_matrix_path)
feature_accessions <- feature_table$BaseAccession
feature_matrix <- as.matrix(feature_table[, -"BaseAccession"])
storage.mode(feature_matrix) <- "double"
rownames(feature_matrix) <- feature_accessions
rm(feature_table)

membership <- fread(membership_path)
set_counts <- fread(set_count_path)
assignment_long <- fread(assignment_path)
protein_metadata <- unique(
  fread(metadata_path)[, .(BaseAccession, GeneSymbol, ProteinName)]
)
pathway_info <- fread(color_path)[order(PathwayOrder)]
recommendation <- fread(recommendation_path)

assert(
  identical(dim(feature_matrix), c(507L, 3008L)) &&
    uniqueN(feature_accessions) == 507L &&
    all(feature_matrix %in% c(0, 1)),
  "Expected the fixed 507 x 3,008 binary BP semantic matrix."
)
assert(
  nrow(membership) == 507L &&
    uniqueN(membership$BaseAccession) == 507L &&
    setequal(membership$BaseAccession, feature_accessions),
  "Four-category membership does not match the fixed 507 proteins."
)
assert(
  nrow(recommendation) == 15L &&
    uniqueN(recommendation$EmbeddingSet) == 5L &&
    uniqueN(recommendation$Method) == 3L &&
    all(recommendation[
      ,
      .N,
      by = .(EmbeddingSet, Method)
    ]$N == 1L),
  "Expected exactly one tuned recommendation for each set-method pair."
)
assert(
  nrow(pathway_info) == 9L &&
    all(assignment_long$Score %in% c(-1L, 1L)) &&
    all(assignment_long$BaseAccession %in% feature_accessions),
  "The signed nine-pathway inputs are invalid."
)

membership_rows <- list()
for (set_index in seq_len(nrow(set_info))) {
  set_row <- set_info[set_index]
  selected <- if (set_row$EmbeddingSet == "all_507") {
    feature_accessions
  } else {
    membership[get(set_row$MembershipColumn) == TRUE, BaseAccession]
  }
  selected <- feature_accessions[feature_accessions %in% selected]
  membership_rows[[set_index]] <- data.table(
    BaseAccession = selected,
    EmbeddingSet = set_row$EmbeddingSet,
    SetLabel = set_row$SetLabel,
    SetLabelZh = set_row$SetLabelZh,
    SetOrder = set_row$SetOrder
  )
}
embedding_set_membership <- rbindlist(membership_rows)
observed_counts <- embedding_set_membership[
  ,
  .(ProteinCount = .N),
  by = .(EmbeddingSet, SetLabel, SetLabelZh, SetOrder)
][order(SetOrder)]
assert(
  identical(observed_counts$ProteinCount, c(507L, 183L, 471L, 178L, 383L)),
  "The five embedding-set sizes are not 507/183/471/178/383."
)
expected_four_counts <- setNames(set_counts$ProteinCount, set_counts$Category)
assert(
  identical(
    observed_counts[EmbeddingSet != "all_507", ProteinCount],
    as.integer(expected_four_counts[
      observed_counts[EmbeddingSet != "all_507", EmbeddingSet]
    ])
  ),
  "The four category counts disagree with the fixed Venn membership."
)

cosine_distance_matrix <- function(x) {
  normalized <- x / sqrt(rowSums(x^2))
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

coordinate_rows <- list()
parameter_rows <- list()
for (set_index in seq_len(nrow(set_info))) {
  set_row <- set_info[set_index]
  accessions <- embedding_set_membership[
    EmbeddingSet == set_row$EmbeddingSet,
    BaseAccession
  ]
  x_binary <- feature_matrix[match(accessions, feature_accessions), , drop = FALSE]
  cosine_distance <- cosine_distance_matrix(x_binary)
  n <- nrow(x_binary)
  umap_rec <- recommendation[
    EmbeddingSet == set_row$EmbeddingSet & Method == "UMAP"
  ]
  tsne_rec <- recommendation[
    EmbeddingSet == set_row$EmbeddingSet & Method == "t-SNE"
  ]
  pca_rec <- recommendation[
    EmbeddingSet == set_row$EmbeddingSet & Method == "PCA"
  ]

  set.seed(random_seed)
  umap_coordinates <- uwot::umap(
    X = x_binary,
    n_neighbors = as.integer(umap_rec$NNeighbors),
    n_components = 2L,
    metric = "cosine",
    n_epochs = as.integer(umap_rec$Iterations),
    scale = FALSE,
    init = "random",
    spread = umap_rec$Spread,
    min_dist = umap_rec$MinDist,
    repulsion_strength = umap_rec$RepulsionStrength,
    negative_sample_rate = as.integer(umap_rec$NegativeSampleRate),
    fast_sgd = FALSE,
    n_threads = 1L,
    n_sgd_threads = 1L,
    seed = random_seed,
    verbose = FALSE
  )

  set.seed(random_seed)
  tsne_fit <- Rtsne(
    as.dist(cosine_distance),
    dims = 2L,
    perplexity = tsne_rec$Perplexity,
    theta = tsne_rec$Theta,
    max_iter = as.integer(tsne_rec$Iterations),
    check_duplicates = FALSE,
    pca = FALSE,
    normalize = FALSE,
    num_threads = 1L,
    verbose = FALSE
  )

  pca_input <- pca_preprocess(x_binary, pca_rec$Preprocessing)
  assert(
    pca_input$FeatureCount == pca_rec$FeatureCount,
    paste0("PCA feature count does not match tuning for ", set_row$EmbeddingSet, ".")
  )
  pca_result <- gram_pca_scores(pca_input$Matrix)

  assert(
    identical(dim(umap_coordinates), c(n, 2L)) &&
      identical(dim(tsne_fit$Y), c(n, 2L)) &&
      identical(dim(pca_result$Scores), c(n, 2L)) &&
      all(is.finite(umap_coordinates)) &&
      all(is.finite(tsne_fit$Y)) &&
      all(is.finite(pca_result$Scores)),
    paste0("Invalid embedding coordinates for ", set_row$EmbeddingSet, ".")
  )

  coordinate_rows[[set_index]] <- data.table(
    BaseAccession = accessions,
    EmbeddingSet = set_row$EmbeddingSet,
    SetLabel = set_row$SetLabel,
    SetLabelZh = set_row$SetLabelZh,
    SetOrder = set_row$SetOrder,
    ProteinCount = n,
    UMAP_1 = umap_coordinates[, 1],
    UMAP_2 = umap_coordinates[, 2],
    TSNE_1 = tsne_fit$Y[, 1],
    TSNE_2 = tsne_fit$Y[, 2],
    PCA_1 = pca_result$Scores[, 1],
    PCA_2 = pca_result$Scores[, 2],
    PC1VariancePercent = pca_result$VariancePercent[[1L]],
    PC2VariancePercent = pca_result$VariancePercent[[2L]],
    UMAPCandidate = umap_rec$Candidate,
    TSNECandidate = tsne_rec$Candidate,
    PCAPreprocessing = pca_rec$Preprocessing,
    PCAFeatureCount = pca_input$FeatureCount
  )

  selected_rows <- rbindlist(list(
    data.table(
      EmbeddingSet = set_row$EmbeddingSet,
      ProteinCount = n,
      Method = "Shared",
      Candidate = NA_character_,
      Parameter = c(
        "CoordinateFitScope",
        "AnalysisKey",
        "SharedBPFeatureDefinitionCount",
        "PathwayScoresUsedForCoordinates",
        "OtherEmbeddingSetsUsedForCoordinates",
        "RandomSeedRecorded"
      ),
      Value = c(
        paste0(n, " proteins from this set only"),
        "isoform-stripped UniProt BaseAccession",
        3008L,
        "FALSE",
        "FALSE",
        random_seed
      )
    ),
    data.table(
      EmbeddingSet = set_row$EmbeddingSet,
      ProteinCount = n,
      Method = "UMAP",
      Candidate = umap_rec$Candidate,
      Parameter = c(
        "Distance",
        "NNeighbors",
        "MinDist",
        "Spread",
        "RepulsionStrength",
        "NegativeSampleRate",
        "NEpochs",
        "Initialization",
        "RandomSeed",
        "TuningCompositeScore"
      ),
      Value = c(
        "cosine",
        umap_rec$NNeighbors,
        umap_rec$MinDist,
        umap_rec$Spread,
        umap_rec$RepulsionStrength,
        umap_rec$NegativeSampleRate,
        umap_rec$Iterations,
        "random",
        random_seed,
        umap_rec$CompositeScore
      )
    ),
    data.table(
      EmbeddingSet = set_row$EmbeddingSet,
      ProteinCount = n,
      Method = "t-SNE",
      Candidate = tsne_rec$Candidate,
      Parameter = c(
        "Distance",
        "Perplexity",
        "Theta",
        "MaxIterations",
        "Threads",
        "RandomSeed",
        "TuningCompositeScore"
      ),
      Value = c(
        "precomputed cosine",
        tsne_rec$Perplexity,
        tsne_rec$Theta,
        tsne_rec$Iterations,
        1L,
        random_seed,
        tsne_rec$CompositeScore
      )
    ),
    data.table(
      EmbeddingSet = set_row$EmbeddingSet,
      ProteinCount = n,
      Method = "PCA",
      Candidate = pca_rec$Candidate,
      Parameter = c(
        "Preprocessing",
        "RetainedFeatureCount",
        "Centered",
        "VarianceScaled",
        "Solver",
        "PC1VariancePercent",
        "PC2VariancePercent",
        "RandomSeedRecorded",
        "TuningCompositeScore"
      ),
      Value = c(
        pca_rec$Preprocessing,
        pca_input$FeatureCount,
        "TRUE",
        "FALSE",
        "eigendecomposition of centered protein-by-protein Gram matrix",
        pca_result$VariancePercent[[1L]],
        pca_result$VariancePercent[[2L]],
        random_seed,
        pca_rec$CompositeScore
      )
    )
  ), use.names = TRUE)
  parameter_rows[[set_index]] <- selected_rows
}

embedding_coordinates <- rbindlist(coordinate_rows)
embedding_parameters <- rbindlist(parameter_rows)
assert(
  nrow(embedding_coordinates) == sum(observed_counts$ProteinCount) &&
    all(is.finite(as.matrix(
      embedding_coordinates[
        ,
        .(UMAP_1, UMAP_2, TSNE_1, TSNE_2, PCA_1, PCA_2)
      ]
    ))),
  "The five tuned coordinate sets are incomplete."
)

pathway_levels <- pathway_info$DisplayLabel
plot_grid <- embedding_set_membership[
  ,
  .(DisplayLabel = pathway_levels),
  by = .(
    BaseAccession,
    EmbeddingSet,
    SetLabel,
    SetLabelZh,
    SetOrder
  )
]
plot_grid <- merge(
  plot_grid,
  assignment_long[, .(BaseAccession, DisplayLabel, Score)],
  by = c("BaseAccession", "DisplayLabel"),
  all.x = TRUE,
  sort = FALSE
)
plot_grid[is.na(Score), Score := 0L]
plot_grid <- merge(
  plot_grid,
  pathway_info[, .(DisplayLabel, Color, PathwayOrder)],
  by = "DisplayLabel",
  all.x = TRUE,
  sort = FALSE
)
plot_grid <- merge(
  plot_grid,
  embedding_coordinates[
    ,
    .(
      BaseAccession,
      EmbeddingSet,
      ProteinCount,
      UMAP_1,
      UMAP_2,
      TSNE_1,
      TSNE_2,
      PCA_1,
      PCA_2,
      PC1VariancePercent,
      PC2VariancePercent,
      UMAPCandidate,
      TSNECandidate,
      PCAPreprocessing,
      PCAFeatureCount
    )
  ],
  by = c("BaseAccession", "EmbeddingSet"),
  all.x = TRUE,
  sort = FALSE
)
plot_grid <- merge(
  plot_grid,
  protein_metadata,
  by = "BaseAccession",
  all.x = TRUE,
  sort = FALSE
)
plot_grid[
  ,
  Status := factor(
    fcase(
      Score == 1L, "Promoting (+1)",
      Score == -1L, "Suppressing (-1)",
      default = "Not assigned (0)"
    ),
    levels = c("Not assigned (0)", "Promoting (+1)", "Suppressing (-1)")
  )
]
plot_grid[
  ,
  PointColor := fcase(Score == 0L, "#AFAFAF", default = Color)
]
plot_grid[
  ,
  DrawOrder := fcase(Score == 0L, 1L, Score == 1L, 2L, default = 3L)
]

pathway_summary <- plot_grid[
  ,
  .(
    ProteinCount = .N,
    PromotingCount = sum(Score == 1L),
    SuppressingCount = sum(Score == -1L),
    UnassignedCount = sum(Score == 0L),
    PathwayColor = unique(Color),
    SetLabel = unique(SetLabel),
    SetLabelZh = unique(SetLabelZh),
    SetOrder = unique(SetOrder),
    PathwayOrder = unique(PathwayOrder)
  ),
  by = .(EmbeddingSet, DisplayLabel)
][order(SetOrder, PathwayOrder)]
pathway_summary[
  ,
  FacetLabel := paste0(
    DisplayLabel,
    "\n+1: ",
    PromotingCount,
    "  |  -1: ",
    SuppressingCount
  )
]
plot_grid <- merge(
  plot_grid,
  pathway_summary[, .(EmbeddingSet, DisplayLabel, FacetLabel)],
  by = c("EmbeddingSet", "DisplayLabel"),
  all.x = TRUE,
  sort = FALSE
)
setorder(plot_grid, SetOrder, PathwayOrder, DrawOrder, BaseAccession)
assert(
  nrow(plot_grid) == sum(observed_counts$ProteinCount) * 9L &&
    all(pathway_summary[
      ,
      PromotingCount + SuppressingCount + UnassignedCount == ProteinCount
    ]),
  "The five-set plotting grid is incomplete."
)

shape_values <- c(
  "Not assigned (0)" = 16,
  "Promoting (+1)" = 16,
  "Suppressing (-1)" = 1
)
size_values <- c(
  "Not assigned (0)" = 0.62,
  "Promoting (+1)" = 1.16,
  "Suppressing (-1)" = 1.62
)
alpha_values <- c(
  "Not assigned (0)" = 0.52,
  "Promoting (+1)" = 0.97,
  "Suppressing (-1)" = 1
)

save_plot_formats <- function(plot, output_directory, stem) {
  dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
  paths <- file.path(
    output_directory,
    paste0(stem, c(".png", ".pdf", ".svg"))
  )
  ggsave(
    paths[[1L]],
    plot,
    width = 12.6,
    height = 11.2,
    units = "in",
    dpi = 600,
    bg = "white"
  )
  ggsave(
    paths[[2L]],
    plot,
    width = 12.6,
    height = 11.2,
    units = "in",
    bg = "white"
  )
  grDevices::svg(
    paths[[3L]],
    width = 12.6,
    height = 11.2,
    onefile = FALSE,
    bg = "white"
  )
  print(plot)
  grDevices::dev.off()
  invisible(paths)
}

method_info <- data.table(
  Method = c("UMAP", "t-SNE", "PCA"),
  MethodKey = c("umap", "tsne", "pca"),
  XColumn = c("UMAP_1", "TSNE_1", "PCA_1"),
  YColumn = c("UMAP_2", "TSNE_2", "PCA_2"),
  MethodOrder = 1:3
)

figure_manifest_rows <- list()
for (set_index in seq_len(nrow(set_info))) {
  set_row <- set_info[set_index]
  set_data <- plot_grid[EmbeddingSet == set_row$EmbeddingSet]
  set_summary <- pathway_summary[EmbeddingSet == set_row$EmbeddingSet]
  set_data[
    ,
    FacetLabel := factor(
      FacetLabel,
      levels = set_summary[order(PathwayOrder), FacetLabel]
    )
  ]
  protein_count <- unique(set_data$ProteinCount)
  for (method_index in seq_len(nrow(method_info))) {
    method_row <- method_info[method_index]
    selected_rec <- recommendation[
      EmbeddingSet == set_row$EmbeddingSet &
        Method == method_row$Method
    ]
    x_label <- fcase(
      method_row$Method == "UMAP", "UMAP 1",
      method_row$Method == "t-SNE", "t-SNE 1",
      default = sprintf(
        "PC1 (%.2f%%)",
        unique(set_data$PC1VariancePercent)
      )
    )
    y_label <- fcase(
      method_row$Method == "UMAP", "UMAP 2",
      method_row$Method == "t-SNE", "t-SNE 2",
      default = sprintf(
        "PC2 (%.2f%%)",
        unique(set_data$PC2VariancePercent)
      )
    )
    parameter_note <- fcase(
      method_row$Method == "UMAP",
      sprintf(
        "tuned nn=%d, min_dist=%.2f, repulsion=%.1f",
        selected_rec$NNeighbors,
        selected_rec$MinDist,
        selected_rec$RepulsionStrength
      ),
      method_row$Method == "t-SNE",
      sprintf(
        "tuned perplexity=%d, theta=%.1f",
        selected_rec$Perplexity,
        selected_rec$Theta
      ),
      default = sprintf(
        "tuned PCA preprocessing=%s, %s features",
        selected_rec$Preprocessing,
        format(selected_rec$FeatureCount, big.mark = ",")
      )
    )

    pathway_plot <- ggplot(
      set_data,
      aes(
        x = .data[[method_row$XColumn]],
        y = .data[[method_row$YColumn]],
        shape = Status,
        size = Status,
        alpha = Status,
        color = I(PointColor)
      )
    ) +
      geom_point(stroke = 0.76) +
      facet_wrap(~FacetLabel, ncol = 3, scales = "fixed") +
      scale_shape_manual(
        values = shape_values,
        name = "Pathway status",
        drop = FALSE
      ) +
      scale_size_manual(values = size_values, guide = "none") +
      scale_alpha_manual(values = alpha_values, guide = "none") +
      coord_equal() +
      labs(
        title = paste0(
          set_row$SetLabel,
          " | ",
          method_row$Method,
          " pathway maps"
        ),
        subtitle = paste0(
          "Independent ",
          protein_count,
          "-protein fit | ",
          parameter_note,
          " | seed 25"
        ),
        x = x_label,
        y = y_label,
        caption = paste0(
          "Pathway-colored solid: promoting (+1); pathway-colored hollow: ",
          "suppressing (-1); medium gray: not assigned to the displayed pathway.\n",
          "Parameters and coordinates were selected independently for this set. ",
          "Absolute axes and positions should not be compared across sets."
        )
      ) +
      guides(
        shape = guide_legend(
          order = 1,
          override.aes = list(
            color = c("#8F8F8F", "#333333", "#333333"),
            size = c(2.1, 2.6, 3.0),
            alpha = 1
          )
        )
      ) +
      theme_classic(base_size = 10, base_family = "sans") +
      theme(
        panel.border = element_rect(
          color = "#B8B8B8",
          fill = NA,
          linewidth = 0.35
        ),
        strip.background = element_rect(
          fill = "#F3F3F3",
          color = "#B8B8B8",
          linewidth = 0.35
        ),
        strip.text = element_text(
          face = "bold",
          size = 9,
          lineheight = 1.08
        ),
        axis.text = element_text(size = 7.4, color = "#333333"),
        axis.title = element_text(size = 9),
        axis.line = element_blank(),
        axis.ticks = element_line(linewidth = 0.3, color = "#555555"),
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 9.2, color = "#4D4D4D"),
        plot.caption = element_text(
          size = 7.8,
          color = "#5A5A5A",
          hjust = 0,
          lineheight = 1.08
        ),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title = element_text(face = "bold", size = 8.7),
        legend.text = element_text(size = 8.1),
        plot.margin = margin(9, 10, 8, 9)
      )

    output_directory <- file.path(figure_dir, set_row$EmbeddingSet)
    stem <- sprintf(
      "%02d_%s_%02d_%s_pathway_3x3_v3",
      set_index,
      set_row$EmbeddingSet,
      method_index,
      method_row$MethodKey
    )
    output_paths <- save_plot_formats(
      pathway_plot,
      output_directory,
      stem
    )
    figure_manifest_rows[[length(figure_manifest_rows) + 1L]] <- data.table(
      EmbeddingSet = set_row$EmbeddingSet,
      SetLabel = set_row$SetLabel,
      SetOrder = set_index,
      Method = method_row$Method,
      MethodOrder = method_index,
      ProteinCount = protein_count,
      IndependentFit = TRUE,
      TunedPerSet = TRUE,
      Candidate = selected_rec$Candidate,
      PNG = relative_path(output_paths[[1L]]),
      PDF = relative_path(output_paths[[2L]]),
      SVG = relative_path(output_paths[[3L]])
    )
  }
}
figure_manifest <- rbindlist(figure_manifest_rows)

coordinate_output_path <- file.path(
  table_dir,
  "embedding_coordinates_5sets_tuned_long.csv"
)
parameter_output_path <- file.path(
  table_dir,
  "embedding_parameters_5sets_tuned.csv"
)
membership_output_path <- file.path(table_dir, "embedding_set_membership.csv")
plot_data_output_path <- file.path(table_dir, "pathway_plot_data_5sets_tuned.csv")
summary_output_path <- file.path(table_dir, "pathway_summary_5sets.csv")
manifest_output_path <- file.path(table_dir, "figure_manifest_15_grids.csv")
input_audit_path <- file.path(table_dir, "input_file_audit.csv")
session_info_path <- file.path(table_dir, "session_info.txt")

fwrite(embedding_coordinates, coordinate_output_path)
fwrite(embedding_parameters, parameter_output_path)
fwrite(embedding_set_membership, membership_output_path)
fwrite(
  plot_grid[
    ,
    .(
      BaseAccession,
      GeneSymbol,
      ProteinName,
      EmbeddingSet,
      SetLabel,
      SetLabelZh,
      SetOrder,
      ProteinCount,
      DisplayLabel,
      PathwayOrder,
      Score,
      Status,
      Color,
      UMAP_1,
      UMAP_2,
      TSNE_1,
      TSNE_2,
      PCA_1,
      PCA_2,
      PC1VariancePercent,
      PC2VariancePercent,
      UMAPCandidate,
      TSNECandidate,
      PCAPreprocessing,
      PCAFeatureCount
    )
  ],
  plot_data_output_path
)
fwrite(pathway_summary, summary_output_path)
fwrite(figure_manifest, manifest_output_path)
fwrite(
  data.table(
    InputRole = c(
      "507 x 3,008 BP semantic feature matrix",
      "fixed 33-group four-category membership",
      "fixed four-category counts",
      "signed pathway assignment",
      "protein display metadata",
      "pathway color key",
      "per-set parameter recommendations",
      "per-set tuning metrics",
      "ranked parameter candidates"
    ),
    Path = relative_path(required_inputs),
    MD5 = unname(tools::md5sum(required_inputs))
  ),
  input_audit_path
)
writeLines(capture.output(sessionInfo()), session_info_path, useBytes = TRUE)

selected_parameter_lines <- recommendation[
  order(
    match(EmbeddingSet, set_info$EmbeddingSet),
    match(Method, c("UMAP", "t-SNE", "PCA"))
  ),
  fcase(
    Method == "UMAP",
    sprintf(
      "- `%s` UMAP：n_neighbors=%d，min_dist=%.2f，repulsion=%.1f。",
      EmbeddingSet,
      NNeighbors,
      MinDist,
      RepulsionStrength
    ),
    Method == "t-SNE",
    sprintf(
      "- `%s` t-SNE：perplexity=%d，theta=%.1f。",
      EmbeddingSet,
      Perplexity,
      Theta
    ),
    default = sprintf(
      "- `%s` PCA：%s，保留%d个特征。",
      EmbeddingSet,
      Preprocessing,
      FeatureCount
    )
  )
]
report_lines <- c(
  "# 五集合逐算法独立调参的九通路图（固定33组，V3）",
  "",
  "## 与V2的区别",
  "",
  "- V2为五个集合独立计算坐标，但同一算法仍使用统一参数；V3进一步为每个集合分别选择UMAP、t-SNE和PCA配置。",
  "- 共5个集合 × 3种算法 = 15张九宫格；每张图只包含本集合蛋白，九个通路面板复用该图坐标。",
  "- V1和V2结果均保留，V3不覆盖旧产物。",
  "",
  "## 参数选择原则",
  "",
  "- UMAP和t-SNE的逐集合综合评分：10近邻保持率0.50、全局距离Spearman相关0.15、最近邻距离第5百分位0.20、核心点云面积占比0.10、稳健轴向均衡度0.05。",
  "- PCA的逐集合综合评分：10近邻保持率0.30、全局距离Spearman相关0.40、最近邻距离第5百分位0.15、核心点云面积占比0.10、稳健轴向均衡度0.05。",
  "- 核心点云面积指标用于避免少数远端蛋白把主体压缩在画布很小区域；它只用于选择算法参数，不移动任何蛋白坐标。",
  "- 所有随机算法固定seed 25并使用单线程；通路评分不参与坐标计算。",
  "",
  "## 最终逐集合配置",
  "",
  selected_parameter_lines,
  "",
  "## 图形编码与解释边界",
  "",
  "- 通路色实心：评分`+1`；通路色空心：评分`-1`；中灰：未分配至当前面板通路。",
  "- 分析键为去isoform的UniProt `BaseAccession`；GeneSymbol仅用于显示。",
  "- 不同集合使用不同参数且独立拟合，不能直接比较绝对坐标轴、方向或单个蛋白的绝对位置。",
  "",
  "## 产物",
  "",
  paste0("- 15张九宫格清单：`", relative_path(manifest_output_path), "`"),
  paste0("- 五套调参坐标：`", relative_path(coordinate_output_path), "`"),
  paste0("- 最终参数长表：`", relative_path(parameter_output_path), "`"),
  paste0("- 逐集合推荐：`", relative_path(recommendation_path), "`"),
  paste0("- 候选参数排名：`", relative_path(tuning_ranked_path), "`"),
  paste0("- 参数筛选原始指标：`", relative_path(tuning_metrics_path), "`"),
  paste0("- 输入审计：`", relative_path(input_audit_path), "`")
)
writeLines(report_lines, report_path, useBytes = TRUE)

message("Coordinates: ", coordinate_output_path)
message("Parameters: ", parameter_output_path)
message("15-grid manifest: ", manifest_output_path)
message("Figure root: ", figure_dir)
message("Report: ", report_path)
