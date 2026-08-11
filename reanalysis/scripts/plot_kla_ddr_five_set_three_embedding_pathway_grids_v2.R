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
    "reanalysis/scripts/plot_kla_ddr_five_set_three_embedding_pathway_grids_v2.R",
    mustWork = TRUE
  )
}
project_root <- normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = TRUE)

v3_table_dir <- file.path(
  project_root,
  "reanalysis/results/tables/kla_ddr_pathway_pie_umap_33groups_v3_all_go"
)
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
set_count_path <- file.path(venn_table_dir, "set_counts.csv")
assignment_path <- file.path(v3_table_dir, "pathway_assignment_long_v3.csv")
metadata_path <- file.path(v3_table_dir, "pathway_umap_plot_data_v3_all_go.csv")
color_path <- file.path(v4_table_dir, "pathway_color_key_v4.csv")

analysis_name <- "kla_ddr_five_set_three_embedding_pathway_grids_33groups_v2"
table_dir <- file.path(project_root, "reanalysis/results/tables", analysis_name)
figure_dir <- file.path(project_root, "reanalysis/results/figures", analysis_name)
report_path <- file.path(
  project_root,
  "reanalysis/reports/FIVE_SET_UMAP_TSNE_PCA_PATHWAY_GRIDS_33GROUP_V2.md"
)
tuning_recommendation_path <- file.path(
  table_dir,
  "parameter_tuning_recommendation.csv"
)
tuning_metrics_path <- file.path(table_dir, "parameter_tuning_metrics.csv")
tuning_summary_path <- file.path(
  table_dir,
  "parameter_tuning_candidate_summary.csv"
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
  tuning_recommendation_path,
  tuning_metrics_path,
  tuning_summary_path
)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs) > 0L) {
  stop(
    "Missing required input(s): ",
    paste(missing_inputs, collapse = "; "),
    "\nRun reanalysis/scripts/tune_kla_ddr_five_set_embedding_parameters_v2.R first."
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
tuning_recommendation <- fread(tuning_recommendation_path)

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
  nrow(pathway_info) == 9L &&
    uniqueN(pathway_info$DisplayLabel) == 9L &&
    all(assignment_long$Score %in% c(-1L, 1L)) &&
    all(assignment_long$BaseAccession %in% feature_accessions),
  "The signed nine-pathway assignment inputs are invalid."
)

umap_recommendation <- tuning_recommendation[Method == "UMAP"]
tsne_recommendation <- tuning_recommendation[Method == "t-SNE"]
assert(
  nrow(umap_recommendation) == 1L &&
    nrow(tsne_recommendation) == 1L &&
    umap_recommendation$NNeighbors == 12L &&
    umap_recommendation$MinDist == 0.35 &&
    tsne_recommendation$Perplexity == 25,
  "The expected five-set parameter-tuning recommendation is unavailable."
)

embedding_set_membership_rows <- list()
for (set_index in seq_len(nrow(set_info))) {
  set_row <- set_info[set_index]
  selected <- if (set_row$EmbeddingSet == "all_507") {
    feature_accessions
  } else {
    membership[
      get(set_row$MembershipColumn) == TRUE,
      BaseAccession
    ]
  }
  selected <- feature_accessions[feature_accessions %in% selected]
  embedding_set_membership_rows[[set_index]] <- data.table(
    BaseAccession = selected,
    EmbeddingSet = set_row$EmbeddingSet,
    SetLabel = set_row$SetLabel,
    SetLabelZh = set_row$SetLabelZh,
    SetOrder = set_row$SetOrder
  )
}
embedding_set_membership <- rbindlist(embedding_set_membership_rows)
observed_set_counts <- embedding_set_membership[
  ,
  .(ProteinCount = .N),
  by = .(EmbeddingSet, SetLabel, SetLabelZh, SetOrder)
][order(SetOrder)]
assert(
  identical(observed_set_counts$ProteinCount, c(507L, 183L, 471L, 178L, 383L)),
  "The five embedding-set sizes are not 507/183/471/178/383."
)
expected_four_counts <- setNames(set_counts$ProteinCount, set_counts$Category)
assert(
  identical(
    observed_set_counts[EmbeddingSet != "all_507", ProteinCount],
    as.integer(expected_four_counts[
      observed_set_counts[EmbeddingSet != "all_507", EmbeddingSet]
    ])
  ),
  "The four category counts disagree with the fixed Venn membership audit."
)

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

coordinate_rows <- list()
parameter_rows <- list()
for (set_index in seq_len(nrow(set_info))) {
  set_row <- set_info[set_index]
  accessions <- embedding_set_membership[
    EmbeddingSet == set_row$EmbeddingSet,
    BaseAccession
  ]
  x_binary <- feature_matrix[match(accessions, feature_accessions), , drop = FALSE]
  row_norms <- sqrt(rowSums(x_binary^2))
  x_l2 <- x_binary / row_norms
  cosine_distance <- cosine_distance_matrix(x_binary)
  n <- nrow(x_binary)

  set.seed(random_seed)
  umap_coordinates <- uwot::umap(
    X = x_binary,
    n_neighbors = as.integer(umap_recommendation$NNeighbors),
    n_components = 2L,
    metric = "cosine",
    n_epochs = 1000L,
    scale = FALSE,
    init = "random",
    spread = umap_recommendation$Spread,
    min_dist = umap_recommendation$MinDist,
    repulsion_strength = umap_recommendation$RepulsionStrength,
    negative_sample_rate = as.integer(
      umap_recommendation$NegativeSampleRate
    ),
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
    perplexity = tsne_recommendation$Perplexity,
    theta = tsne_recommendation$Theta,
    max_iter = 1000L,
    check_duplicates = FALSE,
    pca = FALSE,
    normalize = FALSE,
    num_threads = 1L,
    verbose = FALSE
  )

  pca_result <- gram_pca_scores(x_l2)
  pca_scores <- pca_result$Scores
  pca_variance_percent <- pca_result$VariancePercent

  assert(
    identical(dim(umap_coordinates), c(n, 2L)) &&
      identical(dim(tsne_fit$Y), c(n, 2L)) &&
      identical(dim(pca_scores), c(n, 2L)) &&
      all(is.finite(umap_coordinates)) &&
      all(is.finite(tsne_fit$Y)) &&
      all(is.finite(pca_scores)) &&
      all(is.finite(pca_variance_percent)),
    paste0("A non-finite or incorrectly sized embedding was produced for ", set_row$EmbeddingSet, ".")
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
    PCA_1 = pca_scores[, 1],
    PCA_2 = pca_scores[, 2],
    PC1VariancePercent = pca_variance_percent[[1L]],
    PC2VariancePercent = pca_variance_percent[[2L]]
  )

  parameter_rows[[set_index]] <- rbindlist(list(
    data.table(
      EmbeddingSet = set_row$EmbeddingSet,
      ProteinCount = n,
      Method = "Shared",
      Parameter = c(
        "CoordinateFitScope",
        "AnalysisKey",
        "InputFeatureCount",
        "FeatureConstruction",
        "PathwayScoresUsedForCoordinates",
        "OtherEmbeddingSetsUsedForCoordinates",
        "RandomSeedRecorded"
      ),
      Value = c(
        paste0(n, " proteins from this set only"),
        "isoform-stripped UniProt BaseAccession",
        ncol(x_binary),
        "direct UniProt BP terms plus available GO.db BP ancestors",
        "FALSE",
        "FALSE",
        random_seed
      )
    ),
    data.table(
      EmbeddingSet = set_row$EmbeddingSet,
      ProteinCount = n,
      Method = "UMAP",
      Parameter = c(
        "Input",
        "Distance",
        "NNeighbors",
        "MinDist",
        "Spread",
        "RepulsionStrength",
        "NegativeSampleRate",
        "NEpochs",
        "Initialization",
        "RandomSeed"
      ),
      Value = c(
        sprintf("%d x 3,008 binary BP semantic matrix", n),
        "cosine",
        umap_recommendation$NNeighbors,
        umap_recommendation$MinDist,
        umap_recommendation$Spread,
        umap_recommendation$RepulsionStrength,
        umap_recommendation$NegativeSampleRate,
        1000L,
        "random",
        random_seed
      )
    ),
    data.table(
      EmbeddingSet = set_row$EmbeddingSet,
      ProteinCount = n,
      Method = "t-SNE",
      Parameter = c(
        "Input",
        "Distance",
        "Perplexity",
        "Theta",
        "MaxIterations",
        "Threads",
        "RandomSeed"
      ),
      Value = c(
        sprintf("precomputed %d x %d cosine distance matrix", n, n),
        "cosine",
        tsne_recommendation$Perplexity,
        tsne_recommendation$Theta,
        1000L,
        1L,
        random_seed
      )
    ),
    data.table(
      EmbeddingSet = set_row$EmbeddingSet,
      ProteinCount = n,
      Method = "PCA",
      Parameter = c(
        "Input",
        "RowL2Normalized",
        "Centered",
        "VarianceScaled",
        "Solver",
        "DeterministicSignOrientation",
        "PC1VariancePercent",
        "PC2VariancePercent",
        "RandomSeedRecorded"
      ),
      Value = c(
        sprintf("%d x 3,008 BP semantic matrix", n),
        "TRUE",
        "TRUE",
        "FALSE",
        "eigendecomposition of the centered protein-by-protein Gram matrix",
        "largest-absolute protein score forced positive per PC",
        format(pca_variance_percent[[1L]], digits = 8),
        format(pca_variance_percent[[2L]], digits = 8),
        random_seed
      )
    )
  ), use.names = TRUE)
}

embedding_coordinates <- rbindlist(coordinate_rows)
embedding_parameters <- rbindlist(parameter_rows)
assert(
  nrow(embedding_coordinates) == sum(observed_set_counts$ProteinCount) &&
    uniqueN(embedding_coordinates$EmbeddingSet) == 5L &&
    all(is.finite(as.matrix(
      embedding_coordinates[
        ,
        .(UMAP_1, UMAP_2, TSNE_1, TSNE_2, PCA_1, PCA_2)
      ]
    ))),
  "The five independently learned coordinate sets are incomplete."
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
      PC2VariancePercent
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
  PointColor := fcase(
    Score == 0L, "#AFAFAF",
    default = Color
  )
]
plot_grid[
  ,
  DrawOrder := fcase(
    Score == 0L, 1L,
    Score == 1L, 2L,
    default = 3L
  )
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
  nrow(plot_grid) == sum(observed_set_counts$ProteinCount) * 9L &&
    all(pathway_summary[
      ,
      PromotingCount + SuppressingCount + UnassignedCount == ProteinCount
    ]),
  "The five-set nine-pathway plotting grid is incomplete."
)

shape_values <- c(
  "Not assigned (0)" = 16,
  "Promoting (+1)" = 16,
  "Suppressing (-1)" = 1
)
size_values <- c(
  "Not assigned (0)" = 0.58,
  "Promoting (+1)" = 1.10,
  "Suppressing (-1)" = 1.55
)
alpha_values <- c(
  "Not assigned (0)" = 0.52,
  "Promoting (+1)" = 0.97,
  "Suppressing (-1)" = 1
)

save_plot_formats <- function(plot, output_directory, stem, width, height) {
  dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
  paths <- file.path(
    output_directory,
    paste0(stem, c(".png", ".pdf", ".svg"))
  )
  ggsave(
    paths[[1L]],
    plot,
    width = width,
    height = height,
    units = "in",
    dpi = 600,
    bg = "white"
  )
  ggsave(
    paths[[2L]],
    plot,
    width = width,
    height = height,
    units = "in",
    bg = "white"
  )
  grDevices::svg(
    paths[[3L]],
    width = width,
    height = height,
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
  facet_levels <- set_summary[order(PathwayOrder), FacetLabel]
  set_data[, FacetLabel := factor(FacetLabel, levels = facet_levels)]
  protein_count <- unique(set_data$ProteinCount)
  pc1_variance <- unique(set_data$PC1VariancePercent)
  pc2_variance <- unique(set_data$PC2VariancePercent)

  for (method_index in seq_len(nrow(method_info))) {
    method_row <- method_info[method_index]
    x_label <- fcase(
      method_row$Method == "UMAP", "UMAP 1",
      method_row$Method == "t-SNE", "t-SNE 1",
      default = sprintf("PC1 (%.2f%%)", pc1_variance)
    )
    y_label <- fcase(
      method_row$Method == "UMAP", "UMAP 2",
      method_row$Method == "t-SNE", "t-SNE 2",
      default = sprintf("PC2 (%.2f%%)", pc2_variance)
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
      geom_point(stroke = 0.74) +
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
          "Independent fit using only these ",
          protein_count,
          " Kla-DDR proteins | 3,008 shared BP semantic features | seed 25"
        ),
        x = x_label,
        y = y_label,
        caption = paste0(
          "Pathway-colored solid: promoting (+1); pathway-colored hollow: ",
          "suppressing (-1); medium gray: not assigned to the displayed pathway.\n",
          "Coordinates were recalculated independently for this protein set. ",
          "Absolute axes and positions should not be compared across sets."
        )
      ) +
      guides(
        shape = guide_legend(
          order = 1,
          override.aes = list(
            color = c("#8F8F8F", "#333333", "#333333"),
            size = c(2.1, 2.5, 2.9),
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
        plot.subtitle = element_text(size = 9.5, color = "#4D4D4D"),
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
      "%02d_%s_%02d_%s_pathway_3x3_v2",
      set_index,
      set_row$EmbeddingSet,
      method_index,
      method_row$MethodKey
    )
    output_paths <- save_plot_formats(
      pathway_plot,
      output_directory,
      stem,
      width = 12.2,
      height = 11.2
    )
    figure_manifest_rows[[length(figure_manifest_rows) + 1L]] <- data.table(
      EmbeddingSet = set_row$EmbeddingSet,
      SetLabel = set_row$SetLabel,
      SetOrder = set_index,
      Method = method_row$Method,
      MethodOrder = method_index,
      ProteinCount = protein_count,
      IndependentFit = TRUE,
      PNG = relative_path(output_paths[[1L]]),
      PDF = relative_path(output_paths[[2L]]),
      SVG = relative_path(output_paths[[3L]])
    )
  }
}
figure_manifest <- rbindlist(figure_manifest_rows)

coordinate_output_path <- file.path(
  table_dir,
  "embedding_coordinates_5sets_long.csv"
)
parameter_output_path <- file.path(
  table_dir,
  "embedding_parameters_5sets.csv"
)
membership_output_path <- file.path(
  table_dir,
  "embedding_set_membership.csv"
)
plot_data_output_path <- file.path(
  table_dir,
  "pathway_plot_data_5sets.csv"
)
summary_output_path <- file.path(
  table_dir,
  "pathway_summary_5sets.csv"
)
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
      PC2VariancePercent
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
      "fixed 33-group four-category Kla-DDR membership",
      "fixed four-category counts",
      "signed pathway assignment long table",
      "protein display metadata",
      "pathway color key",
      "five-set parameter recommendation",
      "five-set tuning metrics",
      "five-set tuning candidate summary"
    ),
    Path = relative_path(required_inputs),
    MD5 = unname(tools::md5sum(required_inputs))
  ),
  input_audit_path
)
writeLines(capture.output(sessionInfo()), session_info_path, useBytes = TRUE)

report_lines <- c(
  "# 五个蛋白集合 × 三种独立降维 × 九通路图（固定33组，V2）",
  "",
  "## 分析范围",
  "",
  "- 第一个集合为固定33组产生的全部507个Kla∩DDR蛋白；其余四个集合分别为正常/非肿瘤组织183、正常/非肿瘤细胞471、癌症组织178、癌症细胞383。",
  "- 每个集合仅使用自身包含的蛋白重新学习UMAP、t-SNE和PCA坐标，共5 × 3 = 15张互不共用坐标的九宫格。",
  "- 分析键为去isoform的UniProt `BaseAccession`；GeneSymbol仅用于显示。",
  "- 所有集合使用相同的3,008个BP语义特征定义；通路评分不参与任何降维。",
  "",
  "## 参数筛选与正式配置",
  "",
  "- 使用五个蛋白集合上的10近邻保持率、高低维全局距离Spearman相关及归一化最近邻距离第5百分位联合评价候选参数；综合评分权重依次为0.70、0.20、0.10。",
  "- UMAP正式配置：cosine、`n_neighbors = 12`、`min_dist = 0.35`、`spread = 1.8`、`repulsion_strength = 1.5`、`negative_sample_rate = 10`、1,000轮、seed 25、单线程。",
  "- t-SNE正式配置：预计算cosine距离、perplexity 25、theta 0.5、1,000轮、seed 25、单线程。",
  "- PCA：对每个蛋白的BP二元向量先作L2归一化，再按特征中心化但不按方差缩放；使用中心化蛋白×蛋白Gram矩阵的特征分解求取与标准PCA等价的样本得分，PC符号确定性定向。该处理使PCA与UMAP/t-SNE采用的cosine语义更一致，并减弱GO注释条目数量的影响。",
  "- 三种方法均记录seed 25；PCA本身为确定性计算。",
  "",
  "## 图形解释",
  "",
  "- 通路色实心：评分`+1`；通路色空心：评分`-1`；中灰：未分配至当前面板通路。",
  "- 每张图只绘制该集合自身蛋白，不再保留集合外蛋白作为浅灰背景。",
  "- 同一张九宫格内九个面板严格复用同一套坐标；不同蛋白集合之间为独立拟合，不能比较绝对坐标轴、方向或单个蛋白的绝对位置。",
  "",
  "## 产物",
  "",
  paste0("- 15张九宫格图清单：`", relative_path(manifest_output_path), "`"),
  paste0("- 五套三种降维坐标：`", relative_path(coordinate_output_path), "`"),
  paste0("- 各集合完整参数：`", relative_path(parameter_output_path), "`"),
  paste0("- 五集合蛋白成员表：`", relative_path(membership_output_path), "`"),
  paste0("- 集合×通路绘图表：`", relative_path(plot_data_output_path), "`"),
  paste0("- 集合×通路计数：`", relative_path(summary_output_path), "`"),
  paste0("- 参数候选评价：`", relative_path(tuning_summary_path), "`"),
  paste0("- 参数筛选原始指标：`", relative_path(tuning_metrics_path), "`"),
  paste0("- 输入文件审计：`", relative_path(input_audit_path), "`")
)
writeLines(report_lines, report_path, useBytes = TRUE)

message("Coordinates: ", coordinate_output_path)
message("Parameters: ", parameter_output_path)
message("15-grid manifest: ", manifest_output_path)
message("Figure root: ", figure_dir)
message("Report: ", report_path)
