#!/usr/bin/env Rscript
# generate_extended_ddr_figures.R
# Generates the full suite of extended diagnostic & mechanistic figures for DDR and DNA_repair models:
#   1. figure1_material_distributions: 28 materials box/jitter + true gold diamond
#   2. figure3_manifold_projection: 3-panel UMAP (Category, True Kla %, Pred Kla %)
#   3. figure4_residual_diagnostics: Residual vs Observed LOESS scatter + Error density
#   4. figure5_category_stratified: 4-category stratified scatter (Normal/Cancer tissue, Normal/Cancer cells)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ragg)
  library(patchwork)
  library(uwot)
})

set.seed(25)

models_root <- "outputs/20260927_ddr_fraction_models"
inputs_dir  <- "outputs/20260925_sample_fraction_inputs_final"
full_rna_path <- "outputs/20260926_full_rna_inputs/full_rna.rds"

cat("1. Loading metadata and computing deterministic UMAP coordinates...\n")
meta <- fread(file.path(inputs_dir, "sample_metadata.csv"))
setorder(meta, ReferenceKey, SampleID)

bundle <- readRDS(full_rna_path)
x <- bundle$x
stopifnot(identical(rownames(x), meta$SampleID))
rm(bundle); gc(FALSE)

# Top 2000 HVGs PCA -> UMAP
vars <- apply(x, 2, var)
top2k <- order(vars, decreasing = TRUE)[1:2000]
pca_res <- prcomp(x[, top2k], scale. = TRUE, rank. = 10)
set.seed(25)
umap_coords <- uwot::umap(pca_res$x[, 1:10], n_neighbors = 20, min_dist = 0.3, seed = 25)

df_umap <- data.table(
  SampleID = meta$SampleID,
  ReferenceKey = meta$ReferenceKey,
  Category = meta$Category,
  UMAP1 = umap_coords[, 1],
  UMAP2 = umap_coords[, 2]
)

cat_colors <- c(
  "normal_tissue" = "#377EB8",
  "cancer_tissue" = "#E41A1C",
  "normal_cells"  = "#4DAF4A",
  "cancer_cells"  = "#984EA3"
)
cat_labels <- c(
  "normal_tissue" = "Normal tissue",
  "cancer_tissue" = "Cancer tissue",
  "normal_cells"  = "Normal cells",
  "cancer_cells"  = "Cancer cells"
)

# Process both DDR and DNA_repair
for (target in c("DDR", "DNA_repair")) {
  target_dir <- file.path(models_root, target)
  fig_dir    <- file.path("outputs/20260928_report_revision", target, "figures")
  tbl_dir    <- file.path(target_dir, "tables")
  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

  cat("Processing target:", target, "...\n")
  p <- fread(file.path(tbl_dir, "predictions.csv.gz"))
  p <- merge(p, meta[, .(SampleID, Category)], by = "SampleID")

  # Split predictions: training split (n=1524) vs within-material test (n=374)
  p_tr <- p[Evaluation == "Training split"]
  p_te <- p[Evaluation == "Within-material test"]
  p_all <- p[Evaluation == "Training fit"] # full model fit

  df_split <- rbind(
    data.table(SampleID = p_tr$SampleID, ReferenceKey = p_tr$ReferenceKey, Category = p_tr$Category,
               Split = "Training (n=1,524)", Label = p_tr$Label, Predicted = p_tr$Predicted),
    data.table(SampleID = p_te$SampleID, ReferenceKey = p_te$ReferenceKey, Category = p_te$Category,
               Split = "Test (n=374)", Label = p_te$Label, Predicted = p_te$Predicted)
  )

  # -------------------------------------------------------------
  # Figure 1: Material Distributions
  # -------------------------------------------------------------
  cat("  Generating figure1_material_distributions...\n")
  mat_order <- df_split[, .(ObsLabel = Label[1]), by = ReferenceKey][order(ObsLabel), ReferenceKey]
  df_split[, ReferenceKey := factor(ReferenceKey, levels = mat_order)]
  mat_labels <- df_split[, .(Label = Label[1], Category = Category[1]), by = ReferenceKey]

  target_display <- if (target == "DDR") "Kla-DDR / DDR" else "Kla-DNA repair / DNA repair"

  g_mat <- ggplot(df_split, aes(x = ReferenceKey, y = Predicted)) +
    geom_boxplot(aes(fill = Category), alpha = 0.25, colour = "grey40", outlier.shape = NA, width = 0.55) +
    geom_jitter(aes(colour = Split), width = 0.18, size = 0.9, alpha = 0.55) +
    geom_point(data = mat_labels, aes(x = ReferenceKey, y = Label),
               shape = 23, size = 3.2, fill = "gold", colour = "black", stroke = 0.9) +
    scale_fill_manual(values = cat_colors, labels = cat_labels, name = "Material Category") +
    scale_colour_manual(values = c("Training (n=1,524)" = "#2B5C8F", "Test (n=374)" = "#E6550D"), name = "Sample Split") +
    scale_y_continuous(breaks = seq(0, 100, by = 10), limits = c(0, max(c(df_split$Predicted, mat_labels$Label), na.rm = TRUE) + 5)) +
    labs(
      title = paste("Figure 1: Individual RNA Sample Predictions vs Shared Material Labels Across 28 Materials (", target_display, ")"),
      subtitle = "Diamonds: Shared material reference percentage | Dots: Individual RNA sample predictions from split-trained RF",
      x = "Material Reference (Ordered by Measured Kla %)",
      y = paste("Predicted", target_display, "Percentage (%)"),
      caption = "Yellow diamonds represent shared material-level reference percentage.\nDots are predicted values from the split model; blue = training samples (n=1,524), orange = test samples (n=374)."
    ) +
    theme_minimal(base_size = 10, base_family = "Arial") +
    theme(
      axis.text.x = element_text(angle = 50, hjust = 1, vjust = 1, size = 8.5),
      panel.grid.minor = element_blank(),
      legend.position = "top",
      plot.title = element_text(face = "bold", size = 11.5),
      plot.subtitle = element_text(size = 9.5, colour = "grey30")
    )

  ggsave(file.path(fig_dir, "figure1_material_distributions.png"), g_mat, width = 14.5, height = 6.2, device = agg_png, dpi = 200)
  ggsave(file.path(fig_dir, "figure1_material_distributions.pdf"), g_mat, width = 14.5, height = 6.2, device = cairo_pdf)

  # -------------------------------------------------------------
  # Figure 3: UMAP Manifold Projection
  # -------------------------------------------------------------
  cat("  Generating figure3_manifold_projection...\n")
  df_dim <- merge(df_umap, p_all[, .(SampleID, Observed = Label, Predicted)], by = "SampleID")

  p_umap_cat <- ggplot(df_dim, aes(UMAP1, UMAP2, colour = Category)) +
    geom_point(size = 1.1, alpha = 0.65) +
    scale_colour_manual(values = cat_colors, labels = cat_labels, name = "Category") +
    labs(title = "A. Biological Category", x = "UMAP 1", y = "UMAP 2") +
    theme_minimal(base_size = 9.5, base_family = "Arial") +
    theme(legend.position = "bottom", plot.title = element_text(face = "bold", size = 10.5))

  max_val <- max(c(df_dim$Observed, df_dim$Predicted), na.rm = TRUE)
  mid_val <- median(df_dim$Observed, na.rm = TRUE)

  p_umap_obs <- ggplot(df_dim, aes(UMAP1, UMAP2, colour = Observed)) +
    geom_point(size = 1.1, alpha = 0.7) +
    scale_colour_gradient2(low = "#2166AC", mid = "#FFFFBF", high = "#B2182B", midpoint = mid_val, name = "Reference %") +
    labs(title = paste("B. Material Reference", target_display, "%"), x = "UMAP 1", y = "UMAP 2") +
    theme_minimal(base_size = 9.5, base_family = "Arial") +
    theme(legend.position = "bottom", plot.title = element_text(face = "bold", size = 10.5))

  p_umap_pred <- ggplot(df_dim, aes(UMAP1, UMAP2, colour = Predicted)) +
    geom_point(size = 1.1, alpha = 0.7) +
    scale_colour_gradient2(low = "#2166AC", mid = "#FFFFBF", high = "#B2182B", midpoint = mid_val, name = "Pred %") +
    labs(title = paste("C. Full-data Fitted", target_display, "%"), x = "UMAP 1", y = "UMAP 2") +
    theme_minimal(base_size = 9.5, base_family = "Arial") +
    theme(legend.position = "bottom", plot.title = element_text(face = "bold", size = 10.5))

  g_umap <- (p_umap_cat | p_umap_obs | p_umap_pred) +
    plot_annotation(
      title = paste("Figure 3: Transcriptomic Manifold Embedding (UMAP) Across 1,898 RNA Samples (", target_display, ")"),
      subtitle = "Embedding on top 2,000 highly variable genes (10 PCs). Comparison between shared material references and in-sample fitted values (not validation).",
      theme = theme(
        plot.title = element_text(face = "bold", size = 12.5, family = "Arial"),
        plot.subtitle = element_text(size = 10.5, family = "Arial", colour = "grey30")
      )
    )

  ggsave(file.path(fig_dir, "figure3_manifold_projection.png"), g_umap, width = 14.5, height = 5.2, device = agg_png, dpi = 200)
  ggsave(file.path(fig_dir, "figure3_manifold_projection.pdf"), g_umap, width = 14.5, height = 5.2, device = cairo_pdf)

  # -------------------------------------------------------------
  # Figure 4: Residual Diagnostics
  # -------------------------------------------------------------
  cat("  Generating figure4_residual_diagnostics...\n")
  df_res <- rbind(
    data.table(SampleID = p_tr$SampleID, Split = "Training split (n=1,524)",
               Observed = p_tr$Label, Predicted = p_tr$Predicted, Residual = p_tr$Predicted - p_tr$Label),
    data.table(SampleID = p_te$SampleID, Split = "Within-material test (n=374)",
               Observed = p_te$Label, Predicted = p_te$Predicted, Residual = p_te$Predicted - p_te$Label)
  )

  p_res_scatter <- ggplot(df_res, aes(x = Observed, y = Residual, colour = Split)) +
    geom_hline(yintercept = 0, linetype = "solid", colour = "grey40", linewidth = 0.6) +
    geom_hline(yintercept = c(-5, 5), linetype = "dashed", colour = "firebrick", linewidth = 0.5, alpha = 0.7) +
    geom_point(alpha = 0.45, size = 1.1) +
    geom_smooth(method = "loess", se = TRUE, linewidth = 0.8, span = 0.75) +
    scale_colour_manual(values = c("Training split (n=1,524)" = "#2B5C8F", "Within-material test (n=374)" = "#E6550D")) +
    scale_y_continuous(limits = c(-15, 15), breaks = seq(-15, 15, by = 5)) +
    labs(
      title = paste("A. Residuals vs Material Reference", target_display, "%"),
      subtitle = "Dotted red lines denote ±5 percentage point tolerance band; curves show LOESS fit",
      x = paste("Measured", target_display, "Percentage (%)"),
      y = "Residual (Predicted - Measured, pp)"
    ) +
    theme_minimal(base_size = 9.5, base_family = "Arial") +
    theme(legend.position = "top", plot.title = element_text(face = "bold", size = 11))

  p_res_density <- ggplot(df_res, aes(x = Residual, fill = Split, colour = Split)) +
    geom_vline(xintercept = 0, linetype = "solid", colour = "grey40", linewidth = 0.6) +
    geom_vline(xintercept = c(-5, 5), linetype = "dashed", colour = "firebrick", linewidth = 0.5, alpha = 0.7) +
    geom_density(alpha = 0.35, linewidth = 0.7) +
    scale_fill_manual(values = c("Training split (n=1,524)" = "#2B5C8F", "Within-material test (n=374)" = "#E6550D")) +
    scale_colour_manual(values = c("Training split (n=1,524)" = "#2B5C8F", "Within-material test (n=374)" = "#E6550D")) +
    scale_x_continuous(limits = c(-12, 12), breaks = seq(-10, 10, by = 2)) +
    labs(
      title = "B. Error Density Distribution",
      subtitle = "Error concentration near zero is descriptive; assess trends and tails separately.",
      x = "Prediction Error (Residual, pp)",
      y = "Density"
    ) +
    theme_minimal(base_size = 9.5, base_family = "Arial") +
    theme(legend.position = "top", plot.title = element_text(face = "bold", size = 11))

  g_res <- (p_res_scatter | p_res_density) +
    plot_annotation(
      title = paste("Figure 4: Model Residual Diagnostics on Training and Test Splits (", target_display, ")"),
      subtitle = "Evaluation of prediction errors across the range of target percentages and error normality",
      theme = theme(
        plot.title = element_text(face = "bold", size = 12.5, family = "Arial"),
        plot.subtitle = element_text(size = 10.5, family = "Arial", colour = "grey30")
      )
    )

  ggsave(file.path(fig_dir, "figure4_residual_diagnostics.png"), g_res, width = 12, height = 5.2, device = agg_png, dpi = 200)
  ggsave(file.path(fig_dir, "figure4_residual_diagnostics.pdf"), g_res, width = 12, height = 5.2, device = cairo_pdf)

  # -------------------------------------------------------------
  # Figure 5: Category-stratified Performance Scatter
  # -------------------------------------------------------------
  cat("  Generating figure5_category_stratified...\n")
  df_cat <- rbind(
    data.table(SampleID = p_tr$SampleID, Split = "Training split",
               Category = p_tr$Category, Label = p_tr$Label, Predicted = p_tr$Predicted),
    data.table(SampleID = p_te$SampleID, Split = "Within-material test",
               Category = p_te$Category, Label = p_te$Label, Predicted = p_te$Predicted)
  )
  df_cat[, CategoryLabel := cat_labels[Category]]

  cat_metrics <- df_cat[Split == "Within-material test", {
    e <- Predicted - Label
    list(
      N_test = .N,
      MAE = mean(abs(e)),
      R2 = 1 - sum(e^2) / sum((Label - mean(Label))^2),
      Within5 = 100 * mean(abs(e) <= 5)
    )
  }, by = CategoryLabel]
  fwrite(cat_metrics, file.path(tbl_dir, "category_test_metrics.csv"))

  cat_facet_map <- setNames(
    sprintf("%s (Test n=%d)\nMAE %.2f pp | R² %.3f | ±5 pp: %.1f%%",
            cat_metrics$CategoryLabel, cat_metrics$N_test, cat_metrics$MAE, cat_metrics$R2, cat_metrics$Within5),
    cat_metrics$CategoryLabel
  )
  df_cat[, FacetTitle := factor(cat_facet_map[CategoryLabel], levels = unname(cat_facet_map))]

  max_axis <- max(c(df_cat$Label, df_cat$Predicted), na.rm = TRUE) + 5

  g_cat <- ggplot(df_cat, aes(Label, Predicted, colour = Split)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.6) +
    geom_point(alpha = 0.5, size = 1.1) +
    facet_wrap(~FacetTitle, nrow = 1) +
    scale_colour_manual(values = c("Training split" = "#2B5C8F", "Within-material test" = "#E6550D"), name = "Subset") +
    coord_equal(xlim = c(0, max_axis), ylim = c(0, max_axis)) +
    labs(
      title = paste("Figure 5: Performance Stratified by Biological Category (", target_display, ")"),
      subtitle = "Subsets evaluated on RF_mtry135 model; metrics calculated on unseen test samples in each category",
      x = paste("Assigned Material", target_display, "Percentage (%)"),
      y = "Predicted Percentage (%)",
      caption = "Blue points: training split samples (n=1,524); Orange points: unseen within-material test samples (n=374).\nMetrics in facet banners reflect test set performance."
    ) +
    theme_minimal(base_size = 10, base_family = "Arial") +
    theme(
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = "grey80", fill = NA, linewidth = 0.5),
      strip.background = element_rect(fill = "#F0F4F8", colour = "grey80", linewidth = 0.5),
      strip.text = element_text(face = "bold", size = 8.8),
      legend.position = "top",
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 10, colour = "grey30")
    )

  ggsave(file.path(fig_dir, "figure5_category_stratified.png"), g_cat, width = 14.5, height = 5.2, device = agg_png, dpi = 200)
  ggsave(file.path(fig_dir, "figure5_category_stratified.pdf"), g_cat, width = 14.5, height = 5.2, device = cairo_pdf)

  cat("Completed all extended figures for", target, "!\n")
}

cat("ALL EXTENDED FIGURES GENERATED FOR BOTH DDR AND DNA_REPAIR SUCCESSFULLY!\n")
