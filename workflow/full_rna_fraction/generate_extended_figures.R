#!/usr/bin/env Rscript
# generate_extended_figures.R
# Generate supplementary and extended figures for full RNA random forest model (RF_mtry135).
# Incorporates:
#   1. Supplemented 4-panel scatter plot (Training split n=1524, Test split n=374, Full fit n=1898, Source held out n=1898)
#   2. Material-level prediction distributions across all 28 materials
#   3. Top 25 feature importance barplot (Ensembl ID + Symbol)
#   4. 2D Manifold projection (UMAP & PCA) colored by Category, True Kla %, and Predicted Kla %
#   5. Residual diagnostic plots (Residual vs Observed + Density)
#   6. Category-stratified performance scatter plot (4 biological classes)

suppressPackageStartupMessages({
  library(data.table)
  library(ranger)
  library(ggplot2)
  library(ragg)
  library(org.Hs.eg.db)
  library(uwot)
  library(patchwork)
})

set.seed(25)

input_dir <- "outputs/20260926_full_rna_inputs"
model_dir <- "outputs/20260926_full_rna_models"
fig_dir   <- file.path(model_dir, "figures")
tbl_dir   <- file.path(model_dir, "tables")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tbl_dir, showWarnings = FALSE, recursive = TRUE)

cat("1. Loading RNA data and split definitions...\n")
bundle <- readRDS(file.path(input_dir, "full_rna.rds"))
x <- bundle$x
s <- bundle$metadata
y <- s$ObservedPercent
genes <- colnames(x)
rm(bundle); gc(FALSE)

frozen <- fread(file.path(tbl_dir, "split.csv"))
stopifnot(identical(frozen$SampleID, s$SampleID))
te <- which(frozen$Split == "Within-material test")
tr <- setdiff(seq_len(nrow(s)), te)
stopifnot(length(te) == 374, length(tr) == 1524)

# Calculate sample-level case weights to balance materials/donors
z_tr <- s[tr, .(ReferenceKey, DonorKey)]
z_tr[, W := 1 / .N, by = .(ReferenceKey, DonorKey)]
z_tr[, W := W / uniqueN(DonorKey), by = ReferenceKey]

z_all <- s[, .(ReferenceKey, DonorKey)]
z_all[, W := 1 / .N, by = .(ReferenceKey, DonorKey)]
z_all[, W := W / uniqueN(DonorKey), by = ReferenceKey]

cat("2. Training random forest on training split (tr, n=1524) with importance=impurity...\n")
f_tr <- ranger(x = x[tr, , drop = FALSE], y = y[tr], num.trees = 500, mtry = 135,
               min.node.size = 1, case.weights = z_tr$W, seed = 25, num.threads = 4,
               write.forest = TRUE, importance = "impurity")

pred_tr <- predict(f_tr, x[tr, , drop = FALSE])$predictions
pred_te <- predict(f_tr, x[te, , drop = FALSE])$predictions

cat("3. Training random forest on full dataset (all, n=1898) with importance=impurity...\n")
f_all <- ranger(x = x, y = y, num.trees = 500, mtry = 135,
                min.node.size = 1, case.weights = z_all$W, seed = 25, num.threads = 4,
                write.forest = TRUE, importance = "impurity")

pred_all <- predict(f_all, x)$predictions

cat("4. Reading existing Source held out predictions...\n")
existing_preds <- fread(file.path(tbl_dir, "predictions.csv.gz"))
pred_sho <- existing_preds[Model == "RF_mtry135" & Evaluation == "Source held out"]
stopifnot(nrow(pred_sho) == 1898)

# Combine evaluations into a comprehensive evaluation table
eval_rows <- list(
  data.table(Model = "RF_mtry135", Evaluation = "Training split", SampleID = s$SampleID[tr],
             ReferenceKey = s$ReferenceKey[tr], Category = s$Category[tr], ConnGroup = s$ConnGroup[tr],
             Label = y[tr], Predicted = pred_tr),
  data.table(Model = "RF_mtry135", Evaluation = "Within-material test", SampleID = s$SampleID[te],
             ReferenceKey = s$ReferenceKey[te], Category = s$Category[te], ConnGroup = s$ConnGroup[te],
             Label = y[te], Predicted = pred_te),
  data.table(Model = "RF_mtry135", Evaluation = "Full dataset fit", SampleID = s$SampleID,
             ReferenceKey = s$ReferenceKey, Category = s$Category, ConnGroup = s$ConnGroup,
             Label = y, Predicted = pred_all),
  data.table(Model = "RF_mtry135", Evaluation = "Source held out", SampleID = pred_sho$SampleID,
             ReferenceKey = pred_sho$ReferenceKey, Category = s$Category[match(pred_sho$SampleID, s$SampleID)],
             ConnGroup = pred_sho$ConnGroup, Label = pred_sho$Label, Predicted = pred_sho$Predicted)
)
p_all <- rbindlist(eval_rows)
p_all <- merge(p_all, s[, .(SampleID, DonorKey)], by = "SampleID")

# Calculate metrics (Record-level and Material-level)
p_all[, W := 1 / .N, by = .(Model, Evaluation, ReferenceKey, DonorKey)]
p_all[, W := W / uniqueN(DonorKey), by = .(Model, Evaluation, ReferenceKey)]

metrics_ext <- rbindlist(lapply(c("Record", "Material"), function(level) {
  p_all[, {
    w <- if (level == "Record") rep(1, .N) else W
    w <- w / sum(w)
    e <- Predicted - Label
    list(
      Weighting = level,
      N = .N,
      Materials = uniqueN(ReferenceKey),
      MAE = sum(w * abs(e)),
      RMSE = sqrt(sum(w * e^2)),
      R2 = 1 - sum(w * e^2) / sum(w * (Label - sum(w * Label))^2),
      Within2 = 100 * sum(w * (abs(e) <= 2)),
      Within5 = 100 * sum(w * (abs(e) <= 5))
    )
  }, by = .(Model, Evaluation)]
}))

fwrite(p_all, file.path(tbl_dir, "extended_predictions.csv.gz"))
fwrite(metrics_ext, file.path(tbl_dir, "extended_metrics.csv"))
cat("Extended metrics saved:\n")
print(metrics_ext[Weighting == "Material"])

# -------------------------------------------------------------
# Plot A: Supplemented 4-panel Scatter Plot
# -------------------------------------------------------------
cat("5. Generating Supplemented 4-panel Scatter Plot...\n")
panel_info <- metrics_ext[Model == "RF_mtry135" & Weighting == "Material"]
panel_order <- c("Training split", "Within-material test", "Full dataset fit", "Source held out")
panel_info <- panel_info[match(panel_order, Evaluation)]

panel_label <- function(e, n, ma, acc, r2) {
  sprintf("%s (n=%d)\nMAE %.2f pp | ±5 pp: %.1f%% | R² %.3f", e, n, ma, acc, r2)
}
panel_map <- setNames(
  panel_label(panel_info$Evaluation, panel_info$N, panel_info$MAE, panel_info$Within5, panel_info$R2),
  panel_info$Evaluation
)
p_all[, Panel4 := factor(panel_map[Evaluation], levels = unname(panel_map))]

g_4panel <- ggplot(p_all, aes(Label, Predicted)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.6) +
  geom_point(size = 0.85, alpha = 0.35, colour = "#2166AC") +
  facet_wrap(~Panel4, nrow = 1) +
  coord_equal(xlim = c(0, 40), ylim = c(0, 40)) +
  labs(
    title = "Nonlinear Sample-Level Prediction Across Evaluations (RF_mtry135)",
    subtitle = "18,332 common RNA genes | Random Forest (500 trees, mtry=135) | Material-balanced weights",
    x = "Assigned Material Kla Detection Percentage (%)",
    y = "Predicted Percentage (%)",
    caption = "Each dot represents an individual RNA sample. Labels are shared material-level measurements.\nPanel 1: Model trained & tested on training split; Panel 2: Model trained on split and blind-tested on held-out samples;\nPanel 3: Model fit on all 1,898 samples; Panel 4: 14-fold leave-one-source-out cross-validation."
  ) +
  theme_minimal(base_size = 10.5, base_family = "Arial") +
  theme(
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "grey80", fill = NA, linewidth = 0.5),
    strip.background = element_rect(fill = "#F0F4F8", colour = "grey80", linewidth = 0.5),
    strip.text = element_text(face = "bold", size = 9),
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 10, colour = "grey30")
  )

ggsave(file.path(fig_dir, "RF_mtry135_4panel.pdf"), g_4panel, width = 16, height = 5.2, device = cairo_pdf)
ggsave(file.path(fig_dir, "RF_mtry135_4panel.png"), g_4panel, width = 16, height = 5.2, device = agg_png, dpi = 200)

# -------------------------------------------------------------
# Plot 1: Material-level Prediction Distributions (Figure 1)
# -------------------------------------------------------------
cat("6. Generating Figure 1: Material-level Prediction Distributions...\n")
# Combine predictions on the 1,524 tr and 374 te (from the f_tr model, so predictions are from the split model)
df_split_preds <- rbind(
  data.table(SampleID = s$SampleID[tr], ReferenceKey = s$ReferenceKey[tr], Category = s$Category[tr],
             Split = "Training (n=1,524)", Label = y[tr], Predicted = pred_tr),
  data.table(SampleID = s$SampleID[te], ReferenceKey = s$ReferenceKey[te], Category = s$Category[te],
             Split = "Test (n=374)", Label = y[te], Predicted = pred_te)
)
# Order materials by Observed Kla Label ascending
mat_order <- df_split_preds[, .(ObsLabel = Label[1]), by = ReferenceKey][order(ObsLabel), ReferenceKey]
df_split_preds[, ReferenceKey := factor(ReferenceKey, levels = mat_order)]

# Material true labels table
mat_labels <- df_split_preds[, .(Label = Label[1], Category = Category[1]), by = ReferenceKey]

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

g_mat <- ggplot(df_split_preds, aes(x = ReferenceKey, y = Predicted)) +
  geom_boxplot(aes(fill = Category), alpha = 0.25, colour = "grey40", outlier.shape = NA, width = 0.55) +
  geom_jitter(aes(colour = Split), width = 0.18, size = 0.9, alpha = 0.55) +
  geom_point(data = mat_labels, aes(x = ReferenceKey, y = Label),
             shape = 23, size = 3, fill = "gold", colour = "black", stroke = 0.9) +
  scale_fill_manual(values = cat_colors, labels = cat_labels, name = "Material Category") +
  scale_colour_manual(values = c("Training (n=1,524)" = "#2B5C8F", "Test (n=374)" = "#E6550D"), name = "Sample Split") +
  scale_y_continuous(breaks = seq(0, 40, by = 5), limits = c(0, 40)) +
  labs(
    title = "Figure 1: Individual RNA Sample Predictions vs Shared Material Labels Across 28 Materials",
    subtitle = "Diamonds: True material Kla detection percentage | Dots: Individual RNA sample predictions from split-trained RF",
    x = "Material Reference (Ordered by Measured Kla %)",
    y = "Predicted Kla Percentage (%)",
    caption = "Yellow diamonds represent true material-level measured Kla detection percentage.\nDots are predicted values from the split model (f_tr); blue = training samples (n=1,524), orange = test samples (n=374)."
  ) +
  theme_minimal(base_size = 10, base_family = "Arial") +
  theme(
    axis.text.x = element_text(angle = 50, hjust = 1, vjust = 1, size = 8.5),
    panel.grid.minor = element_blank(),
    legend.position = "top",
    plot.title = element_text(face = "bold", size = 11.5),
    plot.subtitle = element_text(size = 9.5, colour = "grey30")
  )

ggsave(file.path(fig_dir, "figure1_material_distributions.pdf"), g_mat, width = 14, height = 6.2, device = cairo_pdf)
ggsave(file.path(fig_dir, "figure1_material_distributions.png"), g_mat, width = 14, height = 6.2, device = agg_png, dpi = 200)

# -------------------------------------------------------------
# Plot 2: Top 25 Feature Importance Barplot (Figure 2)
# -------------------------------------------------------------
cat("7. Generating Figure 2: Top 25 Feature Importance Barplot...\n")
imp_vec <- f_all$variable.importance
imp_sorted <- sort(imp_vec, decreasing = TRUE)
top25_ids <- names(imp_sorted)[1:25]

# Map Ensembl IDs to Symbols using org.Hs.eg.db
syms <- suppressMessages(
  mapIds(org.Hs.eg.db, keys = top25_ids, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
)
syms[is.na(syms)] <- "Unannotated"

df_imp <- data.table(
  Rank = 1:25,
  Ensembl = top25_ids,
  Symbol = syms,
  Importance = imp_sorted[1:25],
  RelativeImportance = 100 * imp_sorted[1:25] / imp_sorted[1]
)
df_imp[, GeneLabel := paste0(Symbol, " (", Ensembl, ")")]
df_imp[, GeneLabel := factor(GeneLabel, levels = rev(GeneLabel))]
fwrite(df_imp, file.path(tbl_dir, "feature_importance_top25.csv"))

g_imp <- ggplot(df_imp, aes(x = Importance, y = GeneLabel)) +
  geom_col(fill = "#2B5C8F", alpha = 0.85, width = 0.72) +
  geom_text(aes(label = sprintf("%.1f", Importance)), hjust = -0.15, size = 3, colour = "grey25") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(
    title = "Figure 2: Top 25 Feature Importance (Random Forest Impurity Importance)",
    subtitle = "Extracted from RF_mtry135 full model across 18,332 common Ensembl genes",
    x = "Variable Importance (Gini Impurity Decrease)",
    y = "Gene: Symbol (Ensembl ID)",
    caption = "Computed from 500 trees in the RF_mtry135 model. Gene symbols mapped via Ensembl annotations (org.Hs.eg.db)."
  ) +
  theme_minimal(base_size = 10, base_family = "Arial") +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(size = 8.5, face = "bold"),
    plot.title = element_text(face = "bold", size = 11.5),
    plot.subtitle = element_text(size = 9.5, colour = "grey30")
  )

ggsave(file.path(fig_dir, "figure2_feature_importance.pdf"), g_imp, width = 9.5, height = 7, device = cairo_pdf)
ggsave(file.path(fig_dir, "figure2_feature_importance.png"), g_imp, width = 9.5, height = 7, device = agg_png, dpi = 200)

# -------------------------------------------------------------
# Plot 3: 2D Manifold Projection (Figure 3: UMAP)
# -------------------------------------------------------------
cat("8. Generating Figure 3: Manifold Projection (PCA & UMAP)...\n")
vars <- apply(x, 2, var)
top2k <- order(vars, decreasing = TRUE)[1:2000]
pca_res <- prcomp(x[, top2k], scale. = TRUE, rank. = 10)
set.seed(25)
umap_res <- uwot::umap(pca_res$x[, 1:10], n_neighbors = 20, min_dist = 0.3, seed = 25)

df_dim <- data.table(
  SampleID = s$SampleID,
  ReferenceKey = s$ReferenceKey,
  Category = s$Category,
  Observed = y,
  Predicted = pred_all,
  PC1 = pca_res$x[, 1],
  PC2 = pca_res$x[, 2],
  UMAP1 = umap_res[, 1],
  UMAP2 = umap_res[, 2]
)

p_umap_cat <- ggplot(df_dim, aes(UMAP1, UMAP2, colour = Category)) +
  geom_point(size = 1.1, alpha = 0.65) +
  scale_colour_manual(values = cat_colors, labels = cat_labels, name = "Category") +
  labs(title = "A. Biological Category", x = "UMAP 1", y = "UMAP 2") +
  theme_minimal(base_size = 9.5, base_family = "Arial") +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold", size = 10.5))

p_umap_obs <- ggplot(df_dim, aes(UMAP1, UMAP2, colour = Observed)) +
  geom_point(size = 1.1, alpha = 0.7) +
  scale_colour_gradient2(low = "#2166AC", mid = "#FFFFBF", high = "#B2182B", midpoint = 15, name = "True Kla %") +
  labs(title = "B. True Measured Kla %", x = "UMAP 1", y = "UMAP 2") +
  theme_minimal(base_size = 9.5, base_family = "Arial") +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold", size = 10.5))

p_umap_pred <- ggplot(df_dim, aes(UMAP1, UMAP2, colour = Predicted)) +
  geom_point(size = 1.1, alpha = 0.7) +
  scale_colour_gradient2(low = "#2166AC", mid = "#FFFFBF", high = "#B2182B", midpoint = 15, name = "Pred Kla %") +
  labs(title = "C. Model Predicted Kla %", x = "UMAP 1", y = "UMAP 2") +
  theme_minimal(base_size = 9.5, base_family = "Arial") +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold", size = 10.5))

g_umap <- (p_umap_cat | p_umap_obs | p_umap_pred) +
  plot_annotation(
    title = "Figure 3: Transcriptomic Manifold Embedding (UMAP) Across 1,898 RNA Samples",
    subtitle = "Embedding on top 2,000 highly variable genes (10 PCs). Comparison between True Kla % and Model Predicted Kla %.",
    theme = theme(
      plot.title = element_text(face = "bold", size = 12.5, family = "Arial"),
      plot.subtitle = element_text(size = 10.5, family = "Arial", colour = "grey30")
    )
  )

ggsave(file.path(fig_dir, "figure3_manifold_projection.pdf"), g_umap, width = 14.5, height = 5.2, device = cairo_pdf)
ggsave(file.path(fig_dir, "figure3_manifold_projection.png"), g_umap, width = 14.5, height = 5.2, device = agg_png, dpi = 200)

# -------------------------------------------------------------
# Plot 4: Residual Diagnostics (Figure 4)
# -------------------------------------------------------------
cat("9. Generating Figure 4: Residual Diagnostic Plots...\n")
df_res <- rbind(
  data.table(SampleID = s$SampleID[tr], Split = "Training split (n=1,524)",
             Observed = y[tr], Predicted = pred_tr, Residual = pred_tr - y[tr]),
  data.table(SampleID = s$SampleID[te], Split = "Within-material test (n=374)",
             Observed = y[te], Predicted = pred_te, Residual = pred_te - y[te])
)

p_res_scatter <- ggplot(df_res, aes(x = Observed, y = Residual, colour = Split)) +
  geom_hline(yintercept = 0, linetype = "solid", colour = "grey40", linewidth = 0.6) +
  geom_hline(yintercept = c(-5, 5), linetype = "dashed", colour = "firebrick", linewidth = 0.5, alpha = 0.7) +
  geom_point(alpha = 0.45, size = 1.1) +
  geom_smooth(method = "loess", se = TRUE, linewidth = 0.8, span = 0.75) +
  scale_colour_manual(values = c("Training split (n=1,524)" = "#2B5C8F", "Within-material test (n=374)" = "#E6550D")) +
  scale_y_continuous(limits = c(-15, 15), breaks = seq(-15, 15, by = 5)) +
  labs(
    title = "A. Residuals vs True Measured Kla %",
    subtitle = "Dotted red lines denote ±5 percentage point tolerance band; curves show LOESS fit",
    x = "Measured Kla Detection Percentage (%)",
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
    subtitle = "Sharp peak around 0 indicates unbiased predictions without variance inflation",
    x = "Prediction Error (Residual, pp)",
    y = "Density"
  ) +
  theme_minimal(base_size = 9.5, base_family = "Arial") +
  theme(legend.position = "top", plot.title = element_text(face = "bold", size = 11))

g_res <- (p_res_scatter | p_res_density) +
  plot_annotation(
    title = "Figure 4: Model Residual Diagnostics on Training and Test Splits",
    subtitle = "Evaluation of prediction errors across the range of Kla percentages and error normality",
    theme = theme(
      plot.title = element_text(face = "bold", size = 12.5, family = "Arial"),
      plot.subtitle = element_text(size = 10.5, family = "Arial", colour = "grey30")
    )
  )

ggsave(file.path(fig_dir, "figure4_residual_diagnostics.pdf"), g_res, width = 12, height = 5.2, device = cairo_pdf)
ggsave(file.path(fig_dir, "figure4_residual_diagnostics.png"), g_res, width = 12, height = 5.2, device = agg_png, dpi = 200)

# -------------------------------------------------------------
# Plot 5: Category-stratified Performance Scatter (Figure 5)
# -------------------------------------------------------------
cat("10. Generating Figure 5: Biological Category-Stratified Scatter...\n")
df_cat <- rbind(
  data.table(SampleID = s$SampleID[tr], Split = "Training split",
             Category = s$Category[tr], Label = y[tr], Predicted = pred_tr),
  data.table(SampleID = s$SampleID[te], Split = "Within-material test",
             Category = s$Category[te], Label = y[te], Predicted = pred_te)
)
df_cat[, CategoryLabel := cat_labels[Category]]

# Calculate per-category metrics on the test split
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

g_cat <- ggplot(df_cat, aes(Label, Predicted, colour = Split)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.6) +
  geom_point(alpha = 0.5, size = 1.1) +
  facet_wrap(~FacetTitle, nrow = 1) +
  scale_colour_manual(values = c("Training split" = "#2B5C8F", "Within-material test" = "#E6550D"), name = "Subset") +
  coord_equal(xlim = c(0, 40), ylim = c(0, 40)) +
  labs(
    title = "Figure 5: Performance Stratified by Biological Category (Normal vs Cancer, Tissue vs Cell)",
    subtitle = "Subsets evaluated on RF_mtry135 model; metrics calculated on unseen test samples in each category",
    x = "Assigned Material Kla Percentage (%)",
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

ggsave(file.path(fig_dir, "figure5_category_stratified.pdf"), g_cat, width = 14, height = 5.2, device = cairo_pdf)
ggsave(file.path(fig_dir, "figure5_category_stratified.png"), g_cat, width = 14, height = 5.2, device = agg_png, dpi = 200)

cat("ALL FIGURES AND TABLES GENERATED SUCCESSFULLY!\n")
