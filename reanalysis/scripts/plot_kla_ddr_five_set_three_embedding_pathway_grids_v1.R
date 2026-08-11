#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(Rtsne)
})

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
script_path <- if (length(file_arg) > 0L) {
  normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE)
} else {
  normalizePath(
    paste0(
      "reanalysis/scripts/",
      "plot_kla_ddr_five_set_three_embedding_pathway_grids_v1.R"
    ),
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
umap_coordinate_path <- file.path(
  v4_table_dir,
  "umap_raw_coordinates_v4_bp_semantic.csv"
)
umap_parameter_path <- file.path(v4_table_dir, "umap_v4_parameters.csv")
membership_path <- file.path(venn_table_dir, "membership.csv")
set_count_path <- file.path(venn_table_dir, "set_counts.csv")
assignment_path <- file.path(v3_table_dir, "pathway_assignment_long_v3.csv")
metadata_path <- file.path(v3_table_dir, "pathway_umap_plot_data_v3_all_go.csv")
color_path <- file.path(v4_table_dir, "pathway_color_key_v4.csv")

analysis_name <- "kla_ddr_five_set_three_embedding_pathway_grids_33groups_v1"
table_dir <- file.path(project_root, "reanalysis/results/tables", analysis_name)
figure_dir <- file.path(project_root, "reanalysis/results/figures", analysis_name)
report_path <- file.path(
  project_root,
  paste0(
    "reanalysis/reports/",
    "FIVE_SET_UMAP_TSNE_PCA_PATHWAY_GRIDS_33GROUP_V1.md"
  )
)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

required_inputs <- c(
  feature_matrix_path,
  umap_coordinate_path,
  umap_parameter_path,
  membership_path,
  set_count_path,
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

random_seed <- 25L
category_info <- data.table(
  Category = c(
    "all_507",
    "normal_tissue",
    "normal_cells",
    "cancer_tissue",
    "cancer_cells"
  ),
  MembershipColumn = c(
    "In_all_507",
    "In_normal_tissue",
    "In_normal_cells",
    "In_cancer_tissue",
    "In_cancer_cells"
  ),
  CategoryLabel = c(
    "All Kla-DDR proteins",
    "Normal/non-tumor tissues",
    "Normal/non-tumor cells",
    "Cancer tissues",
    "Cancer cells"
  ),
  CategoryLabelZh = c(
    "全部Kla∩DDR蛋白",
    "正常/非肿瘤组织",
    "正常/非肿瘤细胞",
    "癌症组织",
    "癌症细胞"
  ),
  CategoryOrder = 1:5
)

feature_table <- fread(feature_matrix_path)
feature_accessions <- feature_table$BaseAccession
feature_matrix <- as.matrix(feature_table[, -"BaseAccession"])
storage.mode(feature_matrix) <- "double"
rownames(feature_matrix) <- feature_accessions
rm(feature_table)

raw_umap <- fread(umap_coordinate_path)
umap_parameters <- fread(umap_parameter_path)
membership <- fread(membership_path)
membership[, In_all_507 := TRUE]
set_counts <- fread(set_count_path)
assignment_long <- fread(assignment_path)
protein_metadata <- fread(metadata_path)[
  ,
  .(BaseAccession, GeneSymbol, ProteinName)
]
pathway_info <- fread(color_path)[order(PathwayOrder)]

assert(
  identical(dim(feature_matrix), c(507L, 3008L)) &&
    length(feature_accessions) == 507L &&
    uniqueN(feature_accessions) == 507L &&
    all(feature_matrix %in% c(0, 1)),
  "Expected the fixed 507 x 3,008 binary BP semantic feature matrix."
)
assert(
  nrow(raw_umap) == 507L &&
    uniqueN(raw_umap$BaseAccession) == 507L &&
    setequal(raw_umap$BaseAccession, feature_accessions),
  "Saved raw UMAP coordinates do not match the fixed feature matrix."
)
assert(
  umap_parameters[Parameter == "RandomSeed", Value] == "25" &&
    umap_parameters[Parameter == "PathwayScoresUsedInUMAP", Value] == "FALSE" &&
    umap_parameters[Parameter == "SampleDetectionUsedInUMAP", Value] == "FALSE",
  "The saved UMAP is not the required seed-25 feature-only embedding."
)
assert(
  nrow(membership) == 507L &&
    uniqueN(membership$BaseAccession) == 507L &&
    setequal(membership$BaseAccession, feature_accessions),
  "Four-category Kla-DDR membership does not match the 507-protein scope."
)
assert(
  all(category_info$MembershipColumn %in% names(membership)) &&
    nrow(pathway_info) == 9L &&
    uniqueN(pathway_info$DisplayLabel) == 9L,
  "Required category membership columns or nine-pathway color key are missing."
)
assert(
  all(assignment_long$Score %in% c(-1L, 1L)) &&
    all(assignment_long$BaseAccession %in% feature_accessions) &&
    all(assignment_long$DisplayLabel %in% pathway_info$DisplayLabel),
  "Signed pathway assignments do not match the fixed proteins or pathways."
)

membership_long <- melt(
  membership[
    ,
    c("BaseAccession", category_info$MembershipColumn),
    with = FALSE
  ],
  id.vars = "BaseAccession",
  variable.name = "MembershipColumn",
  value.name = "DetectedInCategory"
)
membership_long <- merge(
  membership_long,
  category_info,
  by = "MembershipColumn",
  all.x = TRUE,
  sort = FALSE
)
setorder(membership_long, CategoryOrder, BaseAccession)
observed_category_counts <- membership_long[
  DetectedInCategory == TRUE,
  .(ProteinCount = .N),
  by = .(
    Category,
    CategoryLabel,
    CategoryLabelZh,
    CategoryOrder
  )
][order(CategoryOrder)]
expected_counts <- setNames(set_counts$ProteinCount, set_counts$Category)
expected_counts <- c(all_507 = 507L, expected_counts)
assert(
  identical(
    observed_category_counts$ProteinCount,
    as.integer(expected_counts[observed_category_counts$Category])
  ) &&
    identical(
      observed_category_counts$ProteinCount,
      c(507L, 183L, 471L, 178L, 383L)
    ),
  paste0(
    "The five fixed protein-set counts are not ",
    "507/183/471/178/383 in requested order."
  )
)

# All three embeddings use the same 507 x 3,008 binary BP semantic matrix.
# UMAP is reused from the audited V4 raw coordinates. t-SNE uses the same
# cosine distance as UMAP; PCA is centered but not variance-scaled so rare
# binary features are not artificially given equal variance.
raw_umap <- raw_umap[match(feature_accessions, BaseAccession)]
assert(
  identical(raw_umap$BaseAccession, feature_accessions),
  "Failed to align the saved UMAP coordinates to the feature matrix."
)

row_norms <- sqrt(rowSums(feature_matrix^2))
assert(all(row_norms > 0), "At least one protein has no retained BP semantic feature.")
cosine_similarity <- tcrossprod(feature_matrix) / outer(row_norms, row_norms)
cosine_similarity[cosine_similarity > 1] <- 1
cosine_similarity[cosine_similarity < -1] <- -1
cosine_distance_matrix <- 1 - cosine_similarity
diag(cosine_distance_matrix) <- 0
cosine_distance <- as.dist(cosine_distance_matrix)

set.seed(random_seed)
tsne_fit <- Rtsne(
  cosine_distance,
  dims = 2,
  perplexity = 30,
  theta = 0.5,
  max_iter = 1000,
  check_duplicates = FALSE,
  pca = FALSE,
  normalize = FALSE,
  num_threads = 1,
  verbose = FALSE
)
assert(
  identical(dim(tsne_fit$Y), c(507L, 2L)) &&
    all(is.finite(tsne_fit$Y)),
  "t-SNE did not return 507 finite two-dimensional coordinates."
)

set.seed(random_seed)
pca_fit <- prcomp(
  feature_matrix,
  center = TRUE,
  scale. = FALSE,
  rank. = 2
)
pca_scores <- pca_fit$x[, 1:2, drop = FALSE]
# Fix the otherwise arbitrary PC sign: the largest-absolute protein score on
# each axis is forced positive. This does not change distances or clustering.
for (component_index in 1:2) {
  anchor_index <- which.max(abs(pca_scores[, component_index]))
  if (pca_scores[anchor_index, component_index] < 0) {
    pca_scores[, component_index] <- -pca_scores[, component_index]
    pca_fit$rotation[, component_index] <-
      -pca_fit$rotation[, component_index]
  }
}
total_feature_variance <- sum(apply(feature_matrix, 2, var))
pca_variance_percent <- 100 * (pca_fit$sdev[1:2]^2) / total_feature_variance
assert(
  identical(dim(pca_scores), c(507L, 2L)) &&
    all(is.finite(pca_scores)) &&
    all(is.finite(pca_variance_percent)),
  "PCA did not return 507 finite two-dimensional coordinates and variance fractions."
)

embedding_coordinates <- data.table(
  BaseAccession = feature_accessions,
  UMAP_1 = raw_umap$UMAP_1,
  UMAP_2 = raw_umap$UMAP_2,
  TSNE_1 = tsne_fit$Y[, 1],
  TSNE_2 = tsne_fit$Y[, 2],
  PCA_1 = pca_scores[, 1],
  PCA_2 = pca_scores[, 2]
)
assert(
  all(is.finite(as.matrix(embedding_coordinates[, -"BaseAccession"]))),
  "At least one saved embedding coordinate is non-finite."
)

method_info <- data.table(
  Method = c("UMAP", "t-SNE", "PCA"),
  MethodKey = c("umap", "tsne", "pca"),
  XColumn = c("UMAP_1", "TSNE_1", "PCA_1"),
  YColumn = c("UMAP_2", "TSNE_2", "PCA_2"),
  XLabel = c(
    "UMAP 1",
    "t-SNE 1",
    sprintf("PC1 (%.2f%%)", pca_variance_percent[[1L]])
  ),
  YLabel = c(
    "UMAP 2",
    "t-SNE 2",
    sprintf("PC2 (%.2f%%)", pca_variance_percent[[2L]])
  ),
  MethodOrder = 1:3
)

embedding_parameters <- rbindlist(
  list(
    data.table(
      Method = "Shared",
      Parameter = c(
        "ProteinScope",
        "AnalysisKey",
        "InputMatrix",
        "FeatureConstruction",
        "PathwayScoresUsedForCoordinates",
        "CategoryDetectionUsedForCoordinates",
        "RandomSeedRecorded"
      ),
      Value = c(
        "507 Kla-intersection-DDR proteins from the fixed 33-group scope",
        "isoform-stripped UniProt BaseAccession",
        "507 x 3,008 binary shared BP semantic features",
        "direct UniProt BP terms plus available GO.db BP ancestors",
        "FALSE",
        "FALSE",
        random_seed
      )
    ),
    data.table(
      Method = "UMAP",
      Parameter = c(
        "CoordinateSource",
        "Distance",
        "NNeighbors",
        "MinDist",
        "Spread",
        "NEpochs",
        "Initialization",
        "RandomSeed"
      ),
      Value = c(
        relative_path(umap_coordinate_path),
        "cosine",
        "8",
        "0.2",
        "1.5",
        "1000",
        "random",
        random_seed
      )
    ),
    data.table(
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
        "precomputed 507 x 507 cosine distance matrix",
        "cosine",
        "30",
        "0.5",
        "1000",
        "1",
        random_seed
      )
    ),
    data.table(
      Method = "PCA",
      Parameter = c(
        "Input",
        "Centered",
        "VarianceScaled",
        "DeterministicSignOrientation",
        "PC1VariancePercent",
        "PC2VariancePercent",
        "RandomSeedRecorded"
      ),
      Value = c(
        "507 x 3,008 binary BP semantic matrix",
        "TRUE",
        "FALSE",
        "largest-absolute protein score forced positive per PC",
        format(pca_variance_percent[[1L]], digits = 8),
        format(pca_variance_percent[[2L]], digits = 8),
        random_seed
      )
    )
  ),
  use.names = TRUE
)

pathway_levels <- pathway_info$DisplayLabel
category_pathway_grid <- CJ(
  BaseAccession = feature_accessions,
  Category = category_info$Category,
  DisplayLabel = pathway_levels,
  unique = TRUE
)
category_pathway_grid <- merge(
  category_pathway_grid,
  membership_long[
    ,
    .(
      BaseAccession,
      Category,
      CategoryLabel,
      CategoryLabelZh,
      CategoryOrder,
      DetectedInCategory
    )
  ],
  by = c("BaseAccession", "Category"),
  all.x = TRUE,
  sort = FALSE
)
category_pathway_grid <- merge(
  category_pathway_grid,
  assignment_long[, .(BaseAccession, DisplayLabel, Score)],
  by = c("BaseAccession", "DisplayLabel"),
  all.x = TRUE,
  sort = FALSE
)
category_pathway_grid[is.na(Score), Score := 0L]
category_pathway_grid <- merge(
  category_pathway_grid,
  pathway_info[, .(DisplayLabel, Color, PathwayOrder)],
  by = "DisplayLabel",
  all.x = TRUE,
  sort = FALSE
)
category_pathway_grid <- merge(
  category_pathway_grid,
  embedding_coordinates,
  by = "BaseAccession",
  all.x = TRUE,
  sort = FALSE
)
category_pathway_grid <- merge(
  category_pathway_grid,
  protein_metadata,
  by = "BaseAccession",
  all.x = TRUE,
  sort = FALSE
)

category_pathway_grid[
  ,
  Status := factor(
    fcase(
      !DetectedInCategory, "Outside category",
      Score == 1L, "Promoting (+1)",
      Score == -1L, "Suppressing (-1)",
      default = "Detected, not assigned (0)"
    ),
    levels = c(
      "Outside category",
      "Detected, not assigned (0)",
      "Promoting (+1)",
      "Suppressing (-1)"
    )
  )
]
category_pathway_grid[
  ,
  PointColor := fcase(
    !DetectedInCategory, "#D8D8D8",
    Score == 0L, "#AFAFAF",
    default = Color
  )
]
category_pathway_grid[
  ,
  DrawOrder := fcase(
    !DetectedInCategory, 1L,
    Score == 0L, 2L,
    Score == 1L, 3L,
    default = 4L
  )
]
setorder(
  category_pathway_grid,
  CategoryOrder,
  PathwayOrder,
  DrawOrder,
  BaseAccession
)

assert(
  nrow(category_pathway_grid) == 507L * 5L * 9L &&
    all(!is.na(category_pathway_grid$DetectedInCategory)) &&
    all(is.finite(category_pathway_grid$UMAP_1)) &&
    all(is.finite(category_pathway_grid$TSNE_1)) &&
    all(is.finite(category_pathway_grid$PCA_1)),
  "The complete five-set pathway plotting grid is invalid."
)

category_pathway_summary <- category_pathway_grid[
  DetectedInCategory == TRUE,
  .(
    CategoryProteinCount = .N,
    PromotingCount = sum(Score == 1L),
    SuppressingCount = sum(Score == -1L),
    DetectedUnassignedCount = sum(Score == 0L),
    PathwayColor = unique(Color),
    CategoryLabel = unique(CategoryLabel),
    CategoryLabelZh = unique(CategoryLabelZh),
    CategoryOrder = unique(CategoryOrder),
    PathwayOrder = unique(PathwayOrder)
  ),
  by = .(Category, DisplayLabel)
][order(CategoryOrder, PathwayOrder)]
category_pathway_summary[
  ,
  FacetLabel := paste0(
    DisplayLabel,
    "\n+1: ",
    PromotingCount,
    "  |  -1: ",
    SuppressingCount
  )
]
category_pathway_grid <- merge(
  category_pathway_grid,
  category_pathway_summary[
    ,
    .(Category, DisplayLabel, FacetLabel)
  ],
  by = c("Category", "DisplayLabel"),
  all.x = TRUE,
  sort = FALSE
)
facet_levels <- category_pathway_summary[
  order(CategoryOrder, PathwayOrder),
  unique(FacetLabel)
]
category_pathway_grid[
  ,
  FacetLabel := factor(FacetLabel, levels = facet_levels)
]
setorder(
  category_pathway_grid,
  CategoryOrder,
  PathwayOrder,
  DrawOrder,
  BaseAccession
)

shape_values <- c(
  "Outside category" = 16,
  "Detected, not assigned (0)" = 16,
  "Promoting (+1)" = 16,
  "Suppressing (-1)" = 1
)
size_values <- c(
  "Outside category" = 0.34,
  "Detected, not assigned (0)" = 0.58,
  "Promoting (+1)" = 1.08,
  "Suppressing (-1)" = 1.52
)
alpha_values <- c(
  "Outside category" = 0.18,
  "Detected, not assigned (0)" = 0.52,
  "Promoting (+1)" = 0.96,
  "Suppressing (-1)" = 1
)

save_plot_formats <- function(plot, output_directory, stem, width, height) {
  dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
  png_path <- file.path(output_directory, paste0(stem, ".png"))
  pdf_path <- file.path(output_directory, paste0(stem, ".pdf"))
  svg_path <- file.path(output_directory, paste0(stem, ".svg"))
  ggsave(
    png_path,
    plot,
    width = width,
    height = height,
    units = "in",
    dpi = 600,
    bg = "white"
  )
  ggsave(
    pdf_path,
    plot,
    width = width,
    height = height,
    units = "in",
    bg = "white"
  )
  grDevices::svg(
    svg_path,
    width = width,
    height = height,
    onefile = FALSE,
    bg = "white"
  )
  print(plot)
  grDevices::dev.off()
  invisible(c(png_path, pdf_path, svg_path))
}

figure_manifest_rows <- list()
for (category_index in seq_len(nrow(category_info))) {
  category_key <- category_info$Category[[category_index]]
  category_label <- category_info$CategoryLabel[[category_index]]
  category_count <- observed_category_counts[
    Category == category_key,
    ProteinCount
  ]
  category_data <- category_pathway_grid[Category == category_key]
  legend_status <- levels(droplevels(category_data$Status))
  legend_color_values <- c(
    "Outside category" = "#D8D8D8",
    "Detected, not assigned (0)" = "#8F8F8F",
    "Promoting (+1)" = "#333333",
    "Suppressing (-1)" = "#333333"
  )
  legend_size_values <- c(
    "Outside category" = 1.8,
    "Detected, not assigned (0)" = 2.1,
    "Promoting (+1)" = 2.5,
    "Suppressing (-1)" = 2.9
  )
  for (method_index in seq_len(nrow(method_info))) {
    method_name <- method_info$Method[[method_index]]
    method_key <- method_info$MethodKey[[method_index]]
    x_column <- method_info$XColumn[[method_index]]
    y_column <- method_info$YColumn[[method_index]]
    x_label <- method_info$XLabel[[method_index]]
    y_label <- method_info$YLabel[[method_index]]
    category_method_plot <- ggplot(
      category_data,
      aes(
        x = .data[[x_column]],
        y = .data[[y_column]],
        shape = Status,
        size = Status,
        alpha = Status,
        color = I(PointColor)
      )
    ) +
      geom_point(stroke = 0.72) +
      facet_wrap(~FacetLabel, ncol = 3, scales = "fixed") +
      scale_shape_manual(
        values = shape_values,
        name = "Category/pathway status",
        drop = TRUE
      ) +
      scale_size_manual(values = size_values, guide = "none") +
      scale_alpha_manual(values = alpha_values, guide = "none") +
      coord_equal() +
      labs(
        title = paste0(category_label, " | ", method_name, " pathway maps"),
        subtitle = paste0(
          "Detected Kla-DDR proteins: ",
          category_count,
          "/507 | coordinates learned once from all 507 proteins | seed 25"
        ),
        x = x_label,
        y = y_label,
        caption = paste0(
          "Pathway-colored solid: detected and promoting (+1); pathway-colored ",
          "hollow: detected and suppressing (-1); medium gray: detected but not ",
          "assigned to the displayed pathway.\nVery light gray: not detected in ",
          "this biological category. The same ",
          method_name,
          " coordinates are reused in all five protein-set figures."
        )
      ) +
      guides(
        shape = guide_legend(
          order = 1,
          override.aes = list(
            color = unname(legend_color_values[legend_status]),
            size = unname(legend_size_values[legend_status]),
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

    output_directory <- file.path(figure_dir, category_key)
    stem <- sprintf(
      "%02d_%s_%02d_%s_pathway_3x3_v1",
      category_index,
      category_key,
      method_index,
      method_key
    )
    output_paths <- save_plot_formats(
      category_method_plot,
      output_directory,
      stem,
      width = 12.2,
      height = 11.2
    )
    figure_manifest_rows[[length(figure_manifest_rows) + 1L]] <- data.table(
      Category = category_key,
      CategoryLabel = category_label,
      CategoryOrder = category_index,
      Method = method_name,
      MethodOrder = method_index,
      DetectedProteinCount = category_count,
      PNG = relative_path(output_paths[[1L]]),
      PDF = relative_path(output_paths[[2L]]),
      SVG = relative_path(output_paths[[3L]])
    )
  }
}
figure_manifest <- rbindlist(figure_manifest_rows)

coordinate_output_path <- file.path(table_dir, "embedding_coordinates_507.csv")
parameter_output_path <- file.path(table_dir, "embedding_parameters.csv")
membership_output_path <- file.path(table_dir, "embedding_set_membership_507_long.csv")
plot_data_output_path <- file.path(table_dir, "embedding_set_pathway_plot_data.csv")
summary_output_path <- file.path(table_dir, "embedding_set_pathway_summary.csv")
manifest_output_path <- file.path(table_dir, "figure_manifest_15_grids.csv")
input_audit_path <- file.path(table_dir, "input_file_audit.csv")
session_info_path <- file.path(table_dir, "session_info.txt")

fwrite(embedding_coordinates, coordinate_output_path)
fwrite(embedding_parameters, parameter_output_path)
fwrite(
  membership_long[
    ,
    .(
      BaseAccession,
      Category,
      CategoryLabel,
      CategoryLabelZh,
      CategoryOrder,
      DetectedInCategory
    )
  ],
  membership_output_path
)
fwrite(
  category_pathway_grid[
    ,
    .(
      BaseAccession,
      GeneSymbol,
      ProteinName,
      Category,
      CategoryLabel,
      CategoryLabelZh,
      CategoryOrder,
      DetectedInCategory,
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
      PCA_2
    )
  ],
  plot_data_output_path
)
fwrite(category_pathway_summary, summary_output_path)
fwrite(figure_manifest, manifest_output_path)
fwrite(
  data.table(
    InputRole = c(
      "507 x 3,008 BP semantic feature matrix",
      "audited raw UMAP coordinates",
      "audited UMAP parameters",
      "fixed 33-group four-category Kla-DDR membership plus all-507 set",
      "fixed four-category counts",
      "signed pathway assignment long table",
      "protein display metadata",
      "pathway color key"
    ),
    Path = relative_path(required_inputs),
    MD5 = unname(tools::md5sum(required_inputs))
  ),
  input_audit_path
)
writeLines(capture.output(sessionInfo()), session_info_path, useBytes = TRUE)

report_lines <- c(
  "# 五集合 × 三种降维 × 九通路图（固定33组，V1）",
  "",
  "## 固定分析范围",
  "",
  "- 使用固定33组产生的507个Kla∩DDR蛋白。",
  "- 分析键为去isoform的UniProt `BaseAccession`；GeneSymbol仅用于显示。",
  "- 五个展示集合为全部507、正常/非肿瘤组织183、正常/非肿瘤细胞471、癌症组织178、癌症细胞383；四类成员来自现有`kla_ddr_four_class_venn/membership.csv`，同一蛋白可在多个类别中被检出。",
  "- 三种降维均使用完全相同的507 × 3,008二元BP语义特征矩阵。通路评分和类别检出均不参与坐标计算。",
  "",
  "## 降维配置",
  "",
  "- UMAP：复用已审计的V4原始坐标；cosine、`n_neighbors = 8`、`min_dist = 0.2`、`spread = 1.5`、1,000轮、seed 25。",
  "- t-SNE：同一矩阵的预计算cosine距离；perplexity 30、theta 0.5、1,000轮、单线程、seed 25。",
  "- PCA：同一二元矩阵中心化但不按方差缩放；PC符号按最大绝对蛋白得分确定性定向。PCA本身无随机性，流程仍记录seed 25。",
  "- 每种算法只计算一套507蛋白坐标，全部507和四个类别的九通路面板均复用；因此同一算法内可直接比较五个集合。",
  "",
  "## 图形编码",
  "",
  "- 通路色实心：该类别检出且评分`+1`。",
  "- 通路色空心：该类别检出且评分`-1`。",
  "- 中灰：该类别检出，但未分配至当前面板通路。",
  "- 极浅灰：未在当前生物类别中检出，仅作为全局坐标背景。",
  "",
  "## 产物",
  "",
  paste0("- 15张九宫格图清单：`", relative_path(manifest_output_path), "`"),
  paste0("- 三种降维坐标：`", relative_path(coordinate_output_path), "`"),
  paste0("- 降维配置：`", relative_path(parameter_output_path), "`"),
  paste0("- 五集合成员长表：`", relative_path(membership_output_path), "`"),
  paste0("- 22,815行集合×通路绘图表：`", relative_path(plot_data_output_path), "`"),
  paste0("- 集合×通路计数：`", relative_path(summary_output_path), "`"),
  paste0("- 输入文件审计：`", relative_path(input_audit_path), "`")
)
writeLines(report_lines, report_path, useBytes = TRUE)

message("Coordinates: ", coordinate_output_path)
message("Parameters: ", parameter_output_path)
message("15-grid manifest: ", manifest_output_path)
message("Figure root: ", figure_dir)
message("Report: ", report_path)
