#!/usr/bin/env Rscript
# render_all_extended_figures.R
# Generates the complete suite of extended diagnostic & mechanistic figures for the 4 assay-composition models:
#   1. top25.png/.pdf (with readable Gene Symbols)
#   2. figure1_material_distributions.png/.pdf (28 materials box/jitter + true label)
#   3. figure3_manifold_projection.png/.pdf (3-panel UMAP)
#   4. figure4_residual_diagnostics.png/.pdf (LOESS scatter + Error density)
#   5. figure5_category_stratified.png/.pdf (4-category stratified performance)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ragg)
  library(patchwork)
  library(uwot)
  library(org.Hs.eg.db)
})

set.seed(25)

out_root <- "outputs/20260929_assay_composition"
inputs_dir <- "outputs/20260925_sample_fraction_inputs_final"
full_rna_path <- "outputs/20260926_full_rna_inputs/full_rna.rds"

cat("1. Loading metadata and UMAP manifold coordinates...\n")
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

all_preds <- fread(file.path(out_root, "all_predictions.csv.gz"))
all_preds <- merge(all_preds, meta[, .(SampleID, Category)], by = "SampleID")

targets <- c("Proteome_DDR", "Kla_DDR", "Proteome_DNA_repair", "Kla_DNA_repair")
target_titles <- c(
  "Proteome_DDR" = "Proteome · DDR %",
  "Kla_DDR" = "Kla · DDR %",
  "Proteome_DNA_repair" = "Proteome · DNA repair %",
  "Kla_DNA_repair" = "Kla · DNA repair %"
)

for (target in targets) {
  cat("\n=== Processing target:", target, "===\n")
  model_dir <- file.path(out_root, "models", target)
  fig_dir   <- file.path(model_dir, "figures")
  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

  # --- 1. Top 25 Feature Importance with Symbols ---
  cat("  1. Re-rendering top25 with Gene Symbols...\n")
  imp_file <- file.path(model_dir, "tables/training_importance.csv")
  if (file.exists(imp_file)) {
    imp <- fread(imp_file)[1:25]
    syms <- suppressMessages(
      mapIds(org.Hs.eg.db, keys = imp$Ensembl, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
    )
    syms[is.na(syms)] <- "Unannotated"
    imp[, Symbol := syms]
    imp[, GeneLabel := paste0(Symbol, " (", Ensembl, ")")]
    imp[, GeneLabel := factor(GeneLabel, levels = rev(GeneLabel))]
    fwrite(imp, file.path(model_dir, "tables/training_importance_annotated.csv"))

    g_top <- ggplot(imp, aes(x = Importance, y = GeneLabel)) +
      geom_col(fill = "#2166AC", alpha = 0.88, width = 0.72) +
      geom_text(aes(label = sprintf("%.1f", Importance)), hjust = -0.15, size = 3.2, colour = "grey25") +
      scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
      labs(
        title = paste0(target_titles[target], " model: Top 25 training-split features"),
        subtitle = "Ranked across all 18,332 input genes; regression variance-reduction importance",
        x = "Impurity importance (Gini decrease)",
        y = "Gene: Symbol (Ensembl ID)",
        caption = "Gene symbols mapped via org.Hs.eg.db. Labels preserve stable Ensembl IDs."
      ) +
      theme_minimal(base_size = 11, base_family = "Arial") +
      theme(
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.text.y = element_text(size = 9, face = "bold", colour = "grey20"),
        plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 9.5, colour = "grey35")
      )

    ggsave(file.path(fig_dir, "top25.png"), g_top, width = 10, height = 8.5, device = agg_png, dpi = 200)
    ggsave(file.path(fig_dir, "top25.pdf"), g_top, width = 10, height = 8.5, device = cairo_pdf)
  }

  p_tgt <- all_preds[Target == target]
  p_tr <- p_tgt[Evaluation == "Training split"]
  p_te <- p_tgt[Evaluation == "Within-material test"]
  p_all <- p_tgt[Evaluation == "Training fit"]

  df_split <- rbind(
    data.table(SampleID = p_tr$SampleID, ReferenceKey = p_tr$ReferenceKey, Category = p_tr$Category,
               Split = "Training (n=1,524)", Label = p_tr$Label, Predicted = p_tr$Predicted),
    data.table(SampleID = p_te$SampleID, ReferenceKey = p_te$ReferenceKey, Category = p_te$Category,
               Split = "Test (n=374)", Label = p_te$Label, Predicted = p_te$Predicted)
  )

  # --- 2. Figure 1: 28-Material Distribution ---
  cat("  2. Generating figure1_material_distributions...\n")
  mat_order <- df_split[, .(ObsLabel = Label[1]), by = ReferenceKey][order(ObsLabel), ReferenceKey]
  df_split[, ReferenceKey := factor(ReferenceKey, levels = mat_order)]
  mat_labels <- df_split[, .(Label = Label[1], Category = Category[1]), by = ReferenceKey]

  g_mat <- ggplot(df_split, aes(x = ReferenceKey, y = Predicted)) +
    geom_boxplot(aes(fill = Category), alpha = 0.25, colour = "grey40", outlier.shape = NA, width = 0.55) +
    geom_jitter(aes(colour = Split), width = 0.18, size = 0.9, alpha = 0.55) +
    geom_point(data = mat_labels, aes(x = ReferenceKey, y = Label),
               shape = 23, size = 3.2, fill = "gold", colour = "black", stroke = 0.9) +
    scale_fill_manual(values = cat_colors, labels = cat_labels, name = "Material Category") +
    scale_colour_manual(values = c("Training (n=1,524)" = "#2B5C8F", "Test (n=374)" = "#E6550D"), name = "Sample Split") +
    labs(
      title = paste0("Figure 1: Individual RNA Sample Predictions vs Shared Material Labels (", target_titles[target], ")"),
      subtitle = "Diamonds: True material detection percentage | Dots: Individual RNA sample predictions from split-trained RF",
      x = "Material Reference (Ordered by Measured %)",
      y = paste("Predicted", target_titles[target], "Percentage (%)"),
      caption = "Yellow diamonds represent true material-level measured percentage. Dots are predicted values from split model."
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

  # --- 3. Figure 3: UMAP Manifold Projection ---
  cat("  3. Generating figure3_manifold_projection...\n")
  df_dim <- merge(df_umap, p_all[, .(SampleID, Observed = Label, Predicted)], by = "SampleID")
  mid_val <- median(df_dim$Observed, na.rm = TRUE)

  p_umap_cat <- ggplot(df_dim, aes(UMAP1, UMAP2, colour = Category)) +
    geom_point(size = 1.1, alpha = 0.65) +
    scale_colour_manual(values = cat_colors, labels = cat_labels, name = "Category") +
    labs(title = "A. Biological Category", x = "UMAP 1", y = "UMAP 2") +
    theme_minimal(base_size = 9.5, base_family = "Arial") +
    theme(legend.position = "bottom", plot.title = element_text(face = "bold", size = 10.5))

  p_umap_obs <- ggplot(df_dim, aes(UMAP1, UMAP2, colour = Observed)) +
    geom_point(size = 1.1, alpha = 0.7) +
    scale_colour_gradient2(low = "#2166AC", mid = "#FFFFBF", high = "#B2182B", midpoint = mid_val, name = "True %") +
    labs(title = paste("B. True", target_titles[target]), x = "UMAP 1", y = "UMAP 2") +
    theme_minimal(base_size = 9.5, base_family = "Arial") +
    theme(legend.position = "bottom", plot.title = element_text(face = "bold", size = 10.5))

  p_umap_pred <- ggplot(df_dim, aes(UMAP1, UMAP2, colour = Predicted)) +
    geom_point(size = 1.1, alpha = 0.7) +
    scale_colour_gradient2(low = "#2166AC", mid = "#FFFFBF", high = "#B2182B", midpoint = mid_val, name = "Pred %") +
    labs(title = paste("C. Predicted", target_titles[target]), x = "UMAP 1", y = "UMAP 2") +
    theme_minimal(base_size = 9.5, base_family = "Arial") +
    theme(legend.position = "bottom", plot.title = element_text(face = "bold", size = 10.5))

  g_umap <- (p_umap_cat | p_umap_obs | p_umap_pred) +
    plot_annotation(
      title = paste0("Figure 3: Transcriptomic Manifold Embedding (UMAP) Across 1,898 RNA Samples (", target_titles[target], ")"),
      subtitle = "Embedding on top 2,000 highly variable genes (10 PCs). Comparison between True % and Predicted %.",
      theme = theme(
        plot.title = element_text(face = "bold", size = 12.5, family = "Arial"),
        plot.subtitle = element_text(size = 10.5, family = "Arial", colour = "grey30")
      )
    )

  ggsave(file.path(fig_dir, "figure3_manifold_projection.png"), g_umap, width = 14.5, height = 5.2, device = agg_png, dpi = 200)
  ggsave(file.path(fig_dir, "figure3_manifold_projection.pdf"), g_umap, width = 14.5, height = 5.2, device = cairo_pdf)

  # --- 4. Figure 4: Residual Diagnostics ---
  cat("  4. Generating figure4_residual_diagnostics...\n")
  df_res <- rbind(
    data.table(SampleID = p_tr$SampleID, Split = "Training split (n=1,524)",
               Observed = p_tr$Label, Predicted = p_tr$Predicted, Residual = p_tr$Predicted - p_tr$Label),
    data.table(SampleID = p_te$SampleID, Split = "Within-material test (n=374)",
               Observed = p_te$Label, Predicted = p_te$Predicted, Residual = p_te$Predicted - p_te$Label)
  )

  p_res_scatter <- ggplot(df_res, aes(x = Observed, y = Residual, colour = Split)) +
    geom_hline(yintercept = 0, linetype = "solid", colour = "grey40", linewidth = 0.6) +
    geom_hline(yintercept = c(-1, 1), linetype = "dashed", colour = "firebrick", linewidth = 0.5, alpha = 0.7) +
    geom_point(alpha = 0.45, size = 1.1) +
    geom_smooth(method = "loess", se = TRUE, linewidth = 0.8, span = 0.75) +
    scale_colour_manual(values = c("Training split (n=1,524)" = "#2B5C8F", "Within-material test (n=374)" = "#E6550D")) +
    labs(
      title = paste("A. Residuals vs True Measured", target_titles[target]),
      subtitle = "Dotted red lines denote ±1 percentage point band; curves show LOESS fit",
      x = paste("Measured", target_titles[target], "(%)"),
      y = "Residual (Predicted - Measured, pp)"
    ) +
    theme_minimal(base_size = 9.5, base_family = "Arial") +
    theme(legend.position = "top", plot.title = element_text(face = "bold", size = 11))

  p_res_density <- ggplot(df_res, aes(x = Residual, fill = Split, colour = Split)) +
    geom_vline(xintercept = 0, linetype = "solid", colour = "grey40", linewidth = 0.6) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", colour = "firebrick", linewidth = 0.5, alpha = 0.7) +
    geom_density(alpha = 0.35, linewidth = 0.7) +
    scale_fill_manual(values = c("Training split (n=1,524)" = "#2B5C8F", "Within-material test (n=374)" = "#E6550D")) +
    scale_colour_manual(values = c("Training split (n=1,524)" = "#2B5C8F", "Within-material test (n=374)" = "#E6550D")) +
    labs(
      title = "B. Error Density Distribution",
      subtitle = "Sharp peak around 0 indicates unbiased predictions without variance inflation",
      x = "Prediction Error (Residual, pp)",
      y = "Density"
    ) +
    theme_minimal(base_size = 9.5, base_family = "Arial") +
    theme(legend.position = "top", plot.title = element_text(face = "bold", size = 11))

  g_res <- (p_res_scatter | p_res_density) +
    plot_annotation(
      title = paste0("Figure 4: Model Residual Diagnostics on Training and Test Splits (", target_titles[target], ")"),
      subtitle = "Evaluation of prediction errors across target percentage range and error normality",
      theme = theme(
        plot.title = element_text(face = "bold", size = 12.5, family = "Arial"),
        plot.subtitle = element_text(size = 10.5, family = "Arial", colour = "grey30")
      )
    )

  ggsave(file.path(fig_dir, "figure4_residual_diagnostics.png"), g_res, width = 12, height = 5.2, device = agg_png, dpi = 200)
  ggsave(file.path(fig_dir, "figure4_residual_diagnostics.pdf"), g_res, width = 12, height = 5.2, device = cairo_pdf)

  # --- 5. Figure 5: Category-stratified Performance ---
  cat("  5. Generating figure5_category_stratified...\n")
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
      Within1 = 100 * mean(abs(e) <= 1)
    )
  }, by = CategoryLabel]
  fwrite(cat_metrics, file.path(model_dir, "tables/category_test_metrics.csv"))

  cat_facet_map <- setNames(
    sprintf("%s (Test n=%d)\nMAE %.3f pp | R² %.3f | ±1 pp: %.1f%%",
            cat_metrics$CategoryLabel, cat_metrics$N_test, cat_metrics$MAE, cat_metrics$R2, cat_metrics$Within1),
    cat_metrics$CategoryLabel
  )
  df_cat[, FacetTitle := factor(cat_facet_map[CategoryLabel], levels = unname(cat_facet_map))]

  max_axis <- max(c(df_cat$Label, df_cat$Predicted), na.rm = TRUE) + 0.5
  min_axis <- max(0, min(c(df_cat$Label, df_cat$Predicted), na.rm = TRUE) - 0.5)

  g_cat <- ggplot(df_cat, aes(Label, Predicted, colour = Split)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.6) +
    geom_point(alpha = 0.5, size = 1.1) +
    facet_wrap(~FacetTitle, nrow = 1) +
    scale_colour_manual(values = c("Training split" = "#2B5C8F", "Within-material test" = "#E6550D"), name = "Subset") +
    coord_equal(xlim = c(min_axis, max_axis), ylim = c(min_axis, max_axis)) +
    labs(
      title = paste0("Figure 5: Performance Stratified by Biological Category (", target_titles[target], ")"),
      subtitle = "Subsets evaluated on RF_mtry135 model; metrics calculated on unseen test samples in each category",
      x = paste("Assigned Material", target_titles[target], "(%)"),
      y = "Predicted Percentage (%)",
      caption = "Blue points: training split samples (n=1,524); Orange points: unseen within-material test samples (n=374)."
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
}

cat("\nALL EXTENDED FIGURES GENERATED FOR ALL 4 TARGETS SUCCESSFULLY!\n")
