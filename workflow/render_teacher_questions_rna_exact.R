#!/usr/bin/env Rscript
# ==============================================================================
# Master Script: Render RNA-seq Core Figures for 31 Biological Materials
# (Teacher Revision — Official MSigDB Hallmark Proliferation & Fixed DDR/Kla Panels)
#
# Deliverables Addressed:
#   1. Cell Proliferation Rate (Official MSigDB Hallmark G2M Checkpoint, 196 genes):
#      - Per-tissue boxplot with all 1,898 individual sample dots.
#      - Per-tissue ranked barplot with error bars and sample counts n.
#      - 4-category Scheme B primary boxplot (matching Figure 1b style, ANOVA).
#   2. Fixed DDR Gene Set (Project Panel, 371 genes in 18,332 space):
#      - Per-tissue boxplot with all 1,898 individual sample dots.
#      - Per-tissue ranked barplot with error bars and sample counts n.
#      - 4-category Scheme B boxplot.
#      - 31-group DDR expressed gene fraction boxplot.
#   3. Lactylation Regulators Heatmap (48 regulators across 31 materials on qsmooth):
#      - qsmooth log2(TPM + 0.5) expression heatmap (continuous gradient, 10 highlight boxes).
#      - Gene-standardized Z-score heatmap (blue-white-red, 10 highlight boxes).
#      - Unboxed clean companion versions for both.
#   4. Lactylated Protein Genes Expression (Kla Union, 5,224 genes):
#      - Per-tissue boxplot with all 1,898 individual sample dots.
#      - Per-tissue ranked barplot with error bars and sample counts n.
#      - Parallel DDR vs Kla-target gene expression comparison across 4 categories (Figure 1a style).
#
# All quantitative analyses are strictly performed on the updated 18,332-gene qsmooth matrices:
#   - outputs/20260918_qsmooth_31group_hgnc/matrices/qsmooth_B_full_log2tpm.tsv.gz (sample level)
#   - outputs/20260918_qsmooth_31group_hgnc/matrices/qsmooth_A_collapsed_log2tpm.tsv.gz (group level)
#
# Usage: Rscript workflow/render_teacher_questions_rna_exact.R [project_root] [out_dir]
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(grid)
  library(scales)
  library(ragg)
  library(matrixStats)
  library(msigdbr)
})

# Personal symbol random seed
set.seed(25)

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1L) normalizePath(args[[1L]], mustWork = TRUE) else normalizePath(".", mustWork = TRUE)
out_dir <- if (length(args) >= 2L) args[[2L]] else file.path(root, "results", "rna_teacher_questions")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

publication_font <- "Arial Unicode MS"
charcoal <- "#2F3437"
muted_text <- "#65717D"
grid_colour <- "#D9DDE3"
mean_colour <- "#C0392B"

category_order <- c("normal_tissue", "cancer_tissue", "normal_cells", "cancer_cells")
category_labels <- c(
  normal_tissue = "non-tumor\ntissues",
  cancer_tissue = "tumor\ntissues",
  normal_cells  = "normal\ncell lines",
  cancer_cells  = "cancer\ncell lines"
)
heatmap_category_labels <- c(
  normal_tissue = "non-tumor tissues",
  cancer_tissue = "tumor tissues",
  cancer_cells  = "cancer cell lines",
  normal_cells  = "normal cell lines"
)
category_fills <- c(
  normal_tissue = "#0072B2",
  cancer_tissue = "#D55E00",
  normal_cells  = "#009E73",
  cancer_cells  = "#CC79A7"
)

format_q_value <- function(q_value) {
  if (!is.finite(q_value)) return("NA")
  if (q_value < 2.2e-16) return("< 2.2e-16")
  formatC(q_value, format = "e", digits = 2)
}

read_table_gz <- function(p) {
  if (endsWith(p, ".gz")) data.table::fread(cmd = paste("gzcat", shQuote(p)))
  else data.table::fread(p)
}

# ==============================================================================
# SECTION 1: LOAD METADATA AND UPDATED QSMOOTH MATRICES (18,332 GENES)
# ==============================================================================
message(">>> Loading sample metadata, 31-group status, and 18,332-gene qsmooth matrices...")

qc <- fread(file.path(root, "outputs", "20260916_expression_extraction", "group_sample_qc.csv"))
status <- fread(file.path(root, "audit", "20260916_full_31_rna_status", "rna_31_group_status.csv"))
groups_31 <- fread(file.path(root, "data", "publication_input", "group_summary_31.csv"))

qsmooth_dir <- file.path(root, "outputs", "20260918_qsmooth_31group_hgnc")
group_exp <- fread(file.path(qsmooth_dir, "group_expansion_31.csv"))

# 1. Sample-level qsmooth matrix (18,332 genes x 1,898 samples)
qsmooth_b_path <- file.path(qsmooth_dir, "matrices", "qsmooth_B_full_log2tpm.tsv.gz")
stopifnot(file.exists(qsmooth_b_path))
qb_dt <- read_table_gz(qsmooth_b_path)
gid_col <- names(qb_dt)[1L]
full_samples <- setdiff(names(qb_dt), gid_col)
message(sprintf("Loaded qsmooth_B: %d genes x %d samples", nrow(qb_dt), length(full_samples)))

# 2. Group-level collapsed qsmooth matrix (18,332 genes x 28 references)
qsmooth_a_path <- file.path(qsmooth_dir, "matrices", "qsmooth_A_collapsed_log2tpm.tsv.gz")
stopifnot(file.exists(qsmooth_a_path))
qa_dt <- read_table_gz(qsmooth_a_path)
qa_gid_col <- names(qa_dt)[1L]
message(sprintf("Loaded qsmooth_A: %d genes x %d reference columns", nrow(qa_dt), ncol(qa_dt) - 1L))

# ==============================================================================
# SECTION 2: PROLIFERATION RATE — OFFICIAL MSIGDB HALLMARK G2M CHECKPOINT
# ==============================================================================
message(">>> Task 1: Calculating proliferation score using official MSigDB Hallmark G2M Checkpoint...")

# Query official MSigDB Hallmark collection for human
h_df <- as.data.table(msigdbr(species = "Homo sapiens", collection = "H"))
g2m_df <- h_df[gs_name == "HALLMARK_G2M_CHECKPOINT"]
g2m_genes <- intersect(g2m_df$ensembl_gene, qb_dt[[gid_col]])
message(sprintf("MSigDB Hallmark G2M Checkpoint: %d genes present in 18,332 qsmooth space", length(g2m_genes)))

# Archive official gene set definition
fwrite(g2m_df[ensembl_gene %in% g2m_genes, .(gene_symbol, ensembl_gene, entrez_id = ncbi_gene)],
       file.path(out_dir, "official_hallmark_g2m_checkpoint_genes_mapped.csv"))

# Sample-level proliferation score: mean log2(qsmooth TPM + 0.5) over 196 hallmark genes
sub_prolif <- as.matrix(qb_dt[match(g2m_genes, get(gid_col)), full_samples, with = FALSE])
sample_prolif_score <- colMeans(sub_prolif)

sample_prolif_dt <- data.table(
  SampleID = full_samples,
  ProliferationScore = sample_prolif_score
)
sample_prolif_dt <- merge(sample_prolif_dt, qc[, .(SampleID, ReferenceKey, GroupIDs)], by = "SampleID")

# Map samples to 31 group rows
sample_31_prolif_list <- list()
for (i in seq_len(nrow(group_exp))) {
  gid <- group_exp$GroupID[[i]]
  ref_k <- group_exp$ReferenceKey[[i]]
  sub_s <- copy(sample_prolif_dt[ReferenceKey == ref_k])
  sub_s[, GroupID := gid]
  sample_31_prolif_list[[i]] <- sub_s
}
sample_31_prolif_dt <- rbindlist(sample_31_prolif_list)
sample_31_prolif_dt <- merge(sample_31_prolif_dt, status[, .(GroupID, PXD, SampleGroup, Category)], by = "GroupID")
sample_31_prolif_dt <- merge(sample_31_prolif_dt, groups_31[, .(PXD, SampleGroup, ReferenceLabelEn, RowOrder)], by = c("PXD", "SampleGroup"))
sample_31_prolif_dt[, Category := factor(Category, levels = category_order)]
sample_31_prolif_dt[, CategoryLabel := factor(Category, levels = category_order, labels = unname(category_labels[category_order]))]

# Group summary statistics (Scheme B)
grp_prolif_dt <- sample_31_prolif_dt[, .(
  ProliferationMedian = median(ProliferationScore),
  ProliferationMean = mean(ProliferationScore),
  ProliferationSD = sd(ProliferationScore),
  ProliferationSEM = sd(ProliferationScore) / sqrt(.N),
  ProliferationIQR = IQR(ProliferationScore),
  N = .N
), by = .(GroupID, ReferenceKey, PXD, SampleGroup, Category, CategoryLabel, ReferenceLabelEn, RowOrder)]
grp_prolif_dt[, X := match(Category, category_order)]

fwrite(sample_31_prolif_dt, file.path(out_dir, "figure1b_hallmark_proliferation_sample_values.csv"))
fwrite(grp_prolif_dt, file.path(out_dir, "figure1b_hallmark_proliferation_31groups.csv"))

# --- 1A. Per-Tissue Proliferation Boxplot with Individual Sample Dots ---
setorder(grp_prolif_dt, Category, ProliferationMedian)
grp_prolif_dt[, DisplayLabel := sprintf("%s  %s", GroupID, ReferenceLabelEn)]
ordered_prolif_labels <- grp_prolif_dt$DisplayLabel
grp_prolif_dt[, FactorLabel := factor(DisplayLabel, levels = ordered_prolif_labels)]

sample_31_prolif_dt <- merge(sample_31_prolif_dt, grp_prolif_dt[, .(GroupID, FactorLabel)], by = "GroupID")
sample_31_prolif_dt[, FactorLabel := factor(FactorLabel, levels = ordered_prolif_labels)]

grp_prolif_summary <- sample_31_prolif_dt[, .(
  N = .N,
  Median = median(ProliferationScore),
  Mean = mean(ProliferationScore),
  Q25 = quantile(ProliferationScore, 0.25),
  Q75 = quantile(ProliferationScore, 0.75),
  MaxScore = max(ProliferationScore)
), by = .(GroupID, FactorLabel, Category, CategoryLabel)]

x_max_prolif <- 7.5

plot_prolif_by_tissue_box <- ggplot(sample_31_prolif_dt, aes(x = ProliferationScore, y = FactorLabel, fill = Category)) +
  geom_boxplot(
    aes(group = FactorLabel),
    width = 0.58, outlier.shape = NA, colour = charcoal,
    linewidth = 0.50, median.linewidth = 0.75, alpha = 0.80
  ) +
  geom_point(
    position = position_jitter(width = 0, height = 0.18, seed = 25),
    shape = 21, size = 1.9, stroke = 0.30, colour = "white", alpha = 0.65
  ) +
  geom_point(
    data = grp_prolif_summary,
    aes(x = Mean, y = FactorLabel),
    inherit.aes = FALSE, shape = 23, size = 2.4, fill = mean_colour, colour = charcoal, stroke = 0.45
  ) +
  geom_text(
    data = grp_prolif_summary,
    aes(x = 7.1, y = FactorLabel, label = paste0("n=", N)),
    inherit.aes = FALSE, family = publication_font, size = 3.6, fontface = "bold", colour = muted_text, hjust = 0
  ) +
  scale_fill_manual(
    values = category_fills,
    labels = c(
      normal_tissue = "non-tumor tissues",
      cancer_tissue = "tumor tissues",
      normal_cells  = "normal cell lines",
      cancer_cells  = "cancer cell lines"
    ),
    name = "Biological category"
  ) +
  scale_x_continuous(limits = c(1.0, 7.6), breaks = seq(1.0, 7.0, by = 1.0), expand = c(0.01, 0)) +
  labs(
    title = "Cell proliferation score distribution across 31 biological materials (sample points on qsmooth)",
    subtitle = paste0("Mean log2(qsmooth TPM + 0.5) over ", length(g2m_genes), " official MSigDB Hallmark G2M Checkpoint genes | Red diamonds = mean | n = sample count"),
    x = "Proliferation score (Hallmark G2M Checkpoint mean log2(qsmooth TPM + 0.5))",
    y = NULL
  ) +
  theme_minimal(base_family = publication_font, base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(colour = "#ECEFF1", linewidth = 0.45),
    panel.grid.major.x = element_line(colour = grid_colour, linewidth = 0.50),
    axis.text.y = element_text(size = 9.5, colour = charcoal),
    axis.text.x = element_text(size = 11.0, colour = charcoal),
    axis.title.x = element_text(size = 12.0, face = "bold", colour = charcoal, margin = margin(t = 8)),
    legend.position = "top", legend.direction = "horizontal",
    legend.title = element_text(size = 11.0, face = "bold", colour = charcoal),
    legend.text = element_text(size = 10.5, colour = charcoal),
    plot.title = element_text(size = 14.0, face = "bold", colour = charcoal, hjust = 0),
    plot.subtitle = element_text(size = 10.5, colour = muted_text, hjust = 0, margin = margin(b = 10)),
    plot.margin = margin(12, 20, 12, 12),
    plot.background = element_rect(fill = "white", colour = NA)
  )

stem_prolif_box <- file.path(out_dir, "Figure_1b_RNA_proliferation_hallmark_by_tissue_boxplot")
ggsave(paste0(stem_prolif_box, ".png"), plot_prolif_by_tissue_box, width = 14.0, height = 11.5, dpi = 300, bg = "white", device = ragg::agg_png)
ggsave(paste0(stem_prolif_box, ".pdf"), plot_prolif_by_tissue_box, width = 14.0, height = 11.5, bg = "white", device = cairo_pdf)

# --- 1B. Per-Tissue Proliferation Barplot (Ranked by Score) ---
grp_prolif_bar <- copy(grp_prolif_dt)
setorder(grp_prolif_bar, ProliferationMedian)
grp_prolif_bar[, BarLabel := factor(DisplayLabel, levels = DisplayLabel)]

plot_prolif_by_tissue_bar <- ggplot(grp_prolif_bar, aes(x = ProliferationMedian, y = BarLabel, fill = Category)) +
  geom_col(width = 0.72, alpha = 0.88, colour = charcoal, linewidth = 0.25) +
  geom_errorbar(
    aes(xmin = pmax(0, ProliferationMedian - ProliferationIQR / 2),
        xmax = ProliferationMedian + ProliferationIQR / 2),
    width = 0.35, colour = charcoal, linewidth = 0.45
  ) +
  geom_text(
    aes(label = sprintf("%.2f (n=%d)", ProliferationMedian, N)),
    hjust = -0.15, size = 3.3, family = publication_font, colour = charcoal
  ) +
  scale_fill_manual(
    values = category_fills,
    labels = c("non-tumor tissues", "tumor tissues", "normal cell lines", "cancer cell lines"),
    name = "Biological category"
  ) +
  scale_x_continuous(limits = c(0, 8.5), breaks = seq(0, 8, by = 1), expand = expansion(mult = c(0, 0.08))) +
  labs(
    title = "Cell proliferation score ranking across 31 biological materials (official MSigDB Hallmark)",
    subtitle = "Group median log2(qsmooth TPM + 0.5) over 196 Hallmark G2M Checkpoint genes | Error bars = IQR",
    x = "Proliferation score (median log2(qsmooth TPM + 0.5))",
    y = NULL
  ) +
  theme_minimal(base_family = publication_font, base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(colour = grid_colour, linewidth = 0.50),
    axis.text.y = element_text(size = 9.5, colour = charcoal),
    axis.text.x = element_text(size = 10.5, colour = charcoal),
    axis.title.x = element_text(size = 12.0, face = "bold", colour = charcoal, margin = margin(t = 6)),
    legend.position = "top", legend.direction = "horizontal",
    legend.title = element_text(size = 11.0, face = "bold", colour = charcoal),
    legend.text = element_text(size = 10.5, colour = charcoal),
    plot.title = element_text(size = 14.5, face = "bold", colour = charcoal, hjust = 0),
    plot.subtitle = element_text(size = 10.5, colour = muted_text, hjust = 0, margin = margin(b = 8)),
    plot.margin = margin(12, 20, 12, 12)
  )

stem_prolif_bar <- file.path(out_dir, "Figure_1b_RNA_proliferation_hallmark_by_tissue_barplot")
ggsave(paste0(stem_prolif_bar, ".png"), plot_prolif_by_tissue_bar, width = 14.0, height = 10.5, dpi = 300, bg = "white", device = ragg::agg_png)
ggsave(paste0(stem_prolif_bar, ".pdf"), plot_prolif_by_tissue_bar, width = 14.0, height = 10.5, bg = "white", device = cairo_pdf)

# --- 1C. Scheme B Primary Figure 1b: 31-Group Proliferation Boxplot (ANOVA) ---
grp_prolif_counts <- grp_prolif_dt[, .(N = .N, MaxScore = max(ProliferationMedian)), by = .(Category, X)]
grp_prolif_cat_stats <- grp_prolif_dt[, .(
  Mean = mean(ProliferationMedian),
  Median = median(ProliferationMedian)
), by = .(Category, X)]
grp_prolif_cat_stats[, c("x_left", "x_right") := list(X - 0.20, X + 0.20)]

aov_prolif <- summary(aov(ProliferationMedian ~ Category, data = grp_prolif_dt))[[1L]]
f_prolif <- aov_prolif["Category", "F value"]
p_prolif <- aov_prolif["Category", "Pr(>F)"]
sub_prolif_cat <- paste0("Four-category one-way ANOVA p = ", format_q_value(p_prolif), " (F = ", sprintf("%.2f", f_prolif), ")")

y_max_cat_prolif <- max(grp_prolif_dt$ProliferationMedian) + 0.85

plot_prolif_cat_box <- ggplot(grp_prolif_dt, aes(x = X, y = ProliferationMedian, fill = Category)) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5), colour = "#E5E7EB", linetype = "dashed", linewidth = 0.5) +
  geom_boxplot(
    aes(group = Category),
    width = 0.58, outlier.shape = NA, colour = charcoal,
    linewidth = 0.55, median.linewidth = 0.75, alpha = 0.82, orientation = "x", na.rm = TRUE
  ) +
  geom_segment(data = grp_prolif_cat_stats, aes(x = x_left, xend = x_right, y = Median, yend = Median), inherit.aes = FALSE, colour = charcoal, linewidth = 0.75, lineend = "round") +
  geom_segment(data = grp_prolif_cat_stats, aes(x = x_left, xend = x_right, y = Mean, yend = Mean), inherit.aes = FALSE, colour = mean_colour, linewidth = 0.80, lineend = "round") +
  geom_point(position = position_jitter(width = 0.095, height = 0, seed = 25), shape = 21, size = 3.0, stroke = 0.55, colour = "white", alpha = 0.92, na.rm = TRUE) +
  geom_text(data = grp_prolif_counts, aes(x = X, y = y_max_cat_prolif * 0.96, label = paste0("n=", N)), inherit.aes = FALSE, family = publication_font, size = 4.2, fontface = "bold", colour = muted_text, vjust = -0.15) +
  scale_fill_manual(values = category_fills, guide = "none", drop = FALSE) +
  scale_x_continuous(breaks = seq_along(category_order), labels = unname(category_labels[category_order]), limits = c(0.5, length(category_order) + 0.5), expand = c(0, 0)) +
  scale_y_continuous(limits = c(1.0, y_max_cat_prolif), breaks = scales::pretty_breaks(n = 6), expand = c(0, 0)) +
  labs(
    title = "Cell proliferation score across four biological categories (31 groups, official Hallmark)",
    subtitle = sub_prolif_cat,
    x = NULL, y = "Proliferation score (Hallmark G2M Checkpoint mean log2(TPM + 0.5))"
  ) +
  theme_minimal(base_family = publication_font, base_size = 14) +
  theme(
    panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(colour = grid_colour, linewidth = 0.50), panel.border = element_blank(),
    axis.line.x = element_line(colour = "#8C939E", linewidth = 0.60), axis.line.y = element_line(colour = "#8C939E", linewidth = 0.60),
    axis.text.x = element_text(size = 13.5, face = "bold", colour = charcoal, lineheight = 0.95, margin = margin(t = 6)),
    axis.text.y = element_text(size = 13.0, colour = charcoal),
    axis.title.y = element_text(size = 14.5, face = "bold", colour = charcoal, margin = margin(r = 10)),
    plot.title = element_text(size = 15.0, face = "bold", colour = charcoal, hjust = 0.5, margin = margin(b = 4)),
    plot.subtitle = element_text(size = 12.0, colour = muted_text, hjust = 0.5, margin = margin(b = 10)),
    plot.margin = margin(12, 16, 12, 12), plot.background = element_rect(fill = "white", colour = NA)
  ) + coord_cartesian(clip = "off")

stem_prolif_cat <- file.path(out_dir, "Figure_1b_RNA_proliferation_hallmark_31groups_boxplot")
ggsave(paste0(stem_prolif_cat, ".png"), plot_prolif_cat_box, width = 8.5, height = 7.0, dpi = 300, bg = "white", device = ragg::agg_png)
ggsave(paste0(stem_prolif_cat, ".pdf"), plot_prolif_cat_box, width = 8.5, height = 7.0, bg = "white", device = cairo_pdf)

message(">>> Task 1 complete: Figure 1b Hallmark proliferation figures saved.")

# ==============================================================================
# SECTION 3: FIXED DDR GENE SET — PER-TISSUE BOXPLOT & BARPLOT
# ==============================================================================
message(">>> Task 2: Calculating fixed DDR gene set expression across tissues and categories...")

ddr_annot <- fread(file.path(root, "outputs", "20260916_ddr_panel_31group", "kla_ddr_annotation.csv"))
ddr_genes <- unique(ddr_annot$EnsemblGeneID[!is.na(ddr_annot$EnsemblGeneID) & nzchar(ddr_annot$EnsemblGeneID)])
ddr_genes <- intersect(ddr_genes, qb_dt[[gid_col]])
message(sprintf("Fixed DDR panel: %d genes present in 18,332 space", length(ddr_genes)))

sub_ddr <- as.matrix(qb_dt[match(ddr_genes, get(gid_col)), full_samples, with = FALSE])
sample_ddr_score <- colMedians(sub_ddr)

# DDR expressed fraction (expressed threshold: log2(TPM + 0.5) >= log2(1.5))
all_sub <- as.matrix(qb_dt[, full_samples, with = FALSE])
tpm_cut_log2 <- log2(1 + 0.5)
expressed_tot <- colSums(all_sub >= tpm_cut_log2)
expressed_ddr <- colSums(sub_ddr >= tpm_cut_log2)
sample_ddr_frac <- 100 * expressed_ddr / expressed_tot

sample_ddr_dt <- data.table(
  SampleID = full_samples,
  DDRExpressionMedian = sample_ddr_score,
  DdrFractionPercentage = sample_ddr_frac
)
sample_ddr_dt <- merge(sample_ddr_dt, qc[, .(SampleID, ReferenceKey, GroupIDs)], by = "SampleID")

sample_31_ddr_list <- list()
for (i in seq_len(nrow(group_exp))) {
  gid <- group_exp$GroupID[[i]]
  ref_k <- group_exp$ReferenceKey[[i]]
  sub_s <- copy(sample_ddr_dt[ReferenceKey == ref_k])
  sub_s[, GroupID := gid]
  sample_31_ddr_list[[i]] <- sub_s
}
sample_31_ddr_dt <- rbindlist(sample_31_ddr_list)
sample_31_ddr_dt <- merge(sample_31_ddr_dt, status[, .(GroupID, PXD, SampleGroup, Category)], by = "GroupID")
sample_31_ddr_dt <- merge(sample_31_ddr_dt, groups_31[, .(PXD, SampleGroup, ReferenceLabelEn, RowOrder)], by = c("PXD", "SampleGroup"))
sample_31_ddr_dt[, Category := factor(Category, levels = category_order)]
sample_31_ddr_dt[, CategoryLabel := factor(Category, levels = category_order, labels = unname(category_labels[category_order]))]

grp_ddr_dt <- sample_31_ddr_dt[, .(
  DDRExpressionMedian = median(DDRExpressionMedian),
  DDRExpressionMean = mean(DDRExpressionMedian),
  DDRExpressionSD = sd(DDRExpressionMedian),
  DDRExpressionSEM = sd(DDRExpressionMedian) / sqrt(.N),
  DDRExpressionIQR = IQR(DDRExpressionMedian),
  DdrFractionMedian = median(DdrFractionPercentage),
  DdrFractionMean = mean(DdrFractionPercentage),
  N = .N
), by = .(GroupID, ReferenceKey, PXD, SampleGroup, Category, CategoryLabel, ReferenceLabelEn, RowOrder)]
grp_ddr_dt[, X := match(Category, category_order)]

fwrite(sample_31_ddr_dt, file.path(out_dir, "figure1a_fixed_ddr_sample_values.csv"))
fwrite(grp_ddr_dt, file.path(out_dir, "figure1a_fixed_ddr_31groups.csv"))

# --- 2A. Per-Tissue DDR Expression Boxplot with Sample Points ---
setorder(grp_ddr_dt, Category, DDRExpressionMedian)
grp_ddr_dt[, DisplayLabel := sprintf("%s  %s", GroupID, ReferenceLabelEn)]
ordered_ddr_labels <- grp_ddr_dt$DisplayLabel
grp_ddr_dt[, FactorLabel := factor(DisplayLabel, levels = ordered_ddr_labels)]

sample_31_ddr_dt <- merge(sample_31_ddr_dt, grp_ddr_dt[, .(GroupID, FactorLabel)], by = "GroupID")
sample_31_ddr_dt[, FactorLabel := factor(FactorLabel, levels = ordered_ddr_labels)]

grp_ddr_summary <- sample_31_ddr_dt[, .(
  N = .N,
  Median = median(DDRExpressionMedian),
  Mean = mean(DDRExpressionMedian),
  Q25 = quantile(DDRExpressionMedian, 0.25),
  Q75 = quantile(DDRExpressionMedian, 0.75),
  MaxScore = max(DDRExpressionMedian)
), by = .(GroupID, FactorLabel, Category, CategoryLabel)]

x_max_ddr <- 6.8

plot_ddr_by_tissue_box <- ggplot(sample_31_ddr_dt, aes(x = DDRExpressionMedian, y = FactorLabel, fill = Category)) +
  geom_boxplot(aes(group = FactorLabel), width = 0.58, outlier.shape = NA, colour = charcoal, linewidth = 0.50, median.linewidth = 0.75, alpha = 0.80) +
  geom_point(position = position_jitter(width = 0, height = 0.18, seed = 25), shape = 21, size = 1.9, stroke = 0.30, colour = "white", alpha = 0.65) +
  geom_point(data = grp_ddr_summary, aes(x = Mean, y = FactorLabel), inherit.aes = FALSE, shape = 23, size = 2.4, fill = mean_colour, colour = charcoal, stroke = 0.45) +
  geom_text(data = grp_ddr_summary, aes(x = x_max_ddr * 0.94, y = FactorLabel, label = paste0("n=", N)), inherit.aes = FALSE, family = publication_font, size = 3.6, fontface = "bold", colour = muted_text, hjust = 0) +
  scale_fill_manual(values = category_fills, labels = c(normal_tissue = "non-tumor tissues", cancer_tissue = "tumor tissues", normal_cells = "normal cell lines", cancer_cells = "cancer cell lines"), name = "Biological category") +
  scale_x_continuous(limits = c(1.2, x_max_ddr), breaks = seq(1.5, 6.5, by = 1.0), expand = c(0.01, 0)) +
  labs(
    title = "Transcriptomic DDR gene expression distribution across 31 biological materials (sample points on qsmooth)",
    subtitle = paste0("Sample median log2(qsmooth TPM + 0.5) across ", length(ddr_genes), " fixed DDR genes | Red diamonds = group mean | n = sample count"),
    x = "DDR gene expression (sample median log2(qsmooth TPM + 0.5))", y = NULL
  ) +
  theme_minimal(base_family = publication_font, base_size = 12) +
  theme(
    panel.grid.minor = element_blank(), panel.grid.major.y = element_line(colour = "#ECEFF1", linewidth = 0.45), panel.grid.major.x = element_line(colour = grid_colour, linewidth = 0.50),
    axis.text.y = element_text(size = 9.5, colour = charcoal), axis.text.x = element_text(size = 11.0, colour = charcoal),
    axis.title.x = element_text(size = 12.0, face = "bold", colour = charcoal, margin = margin(t = 8)),
    legend.position = "top", legend.direction = "horizontal", legend.title = element_text(size = 11.0, face = "bold", colour = charcoal), legend.text = element_text(size = 10.5, colour = charcoal),
    plot.title = element_text(size = 14.0, face = "bold", colour = charcoal, hjust = 0), plot.subtitle = element_text(size = 10.5, colour = muted_text, hjust = 0, margin = margin(b = 10)),
    plot.margin = margin(12, 20, 12, 12), plot.background = element_rect(fill = "white", colour = NA)
  )

stem_ddr_box <- file.path(out_dir, "Figure_1a_RNA_DDR_expression_by_tissue_boxplot")
ggsave(paste0(stem_ddr_box, ".png"), plot_ddr_by_tissue_box, width = 14.0, height = 11.5, dpi = 300, bg = "white", device = ragg::agg_png)
ggsave(paste0(stem_ddr_box, ".pdf"), plot_ddr_by_tissue_box, width = 14.0, height = 11.5, bg = "white", device = cairo_pdf)

# --- 2B. Per-Tissue DDR Expression Barplot (Ranked) ---
grp_ddr_bar <- copy(grp_ddr_dt)
setorder(grp_ddr_bar, DDRExpressionMedian)
grp_ddr_bar[, BarLabel := factor(DisplayLabel, levels = DisplayLabel)]

plot_ddr_by_tissue_bar <- ggplot(grp_ddr_bar, aes(x = DDRExpressionMedian, y = BarLabel, fill = Category)) +
  geom_col(width = 0.72, alpha = 0.88, colour = charcoal, linewidth = 0.25) +
  geom_errorbar(aes(xmin = pmax(0, DDRExpressionMedian - DDRExpressionIQR / 2), xmax = DDRExpressionMedian + DDRExpressionIQR / 2), width = 0.35, colour = charcoal, linewidth = 0.45) +
  geom_text(aes(label = sprintf("%.2f (n=%d)", DDRExpressionMedian, N)), hjust = -0.15, size = 3.3, family = publication_font, colour = charcoal) +
  scale_fill_manual(values = category_fills, labels = c("non-tumor tissues", "tumor tissues", "normal cell lines", "cancer cell lines"), name = "Biological category") +
  scale_x_continuous(limits = c(0, 6.8), breaks = seq(0, 6, by = 1), expand = expansion(mult = c(0, 0.08))) +
  labs(
    title = "Transcriptomic DDR gene expression ranking across 31 biological materials",
    subtitle = paste0("Group median log2(qsmooth TPM + 0.5) across ", length(ddr_genes), " fixed DDR genes | Error bars = IQR"),
    x = "DDR gene expression (median log2(qsmooth TPM + 0.5))", y = NULL
  ) +
  theme_minimal(base_family = publication_font, base_size = 12) +
  theme(
    panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), panel.grid.major.x = element_line(colour = grid_colour, linewidth = 0.50),
    axis.text.y = element_text(size = 9.5, colour = charcoal), axis.text.x = element_text(size = 10.5, colour = charcoal),
    axis.title.x = element_text(size = 12.0, face = "bold", colour = charcoal, margin = margin(t = 6)),
    legend.position = "top", legend.direction = "horizontal", legend.title = element_text(size = 11.0, face = "bold", colour = charcoal), legend.text = element_text(size = 10.5, colour = charcoal),
    plot.title = element_text(size = 14.5, face = "bold", colour = charcoal, hjust = 0), plot.subtitle = element_text(size = 10.5, colour = muted_text, hjust = 0, margin = margin(b = 8)),
    plot.margin = margin(12, 20, 12, 12)
  )

stem_ddr_bar <- file.path(out_dir, "Figure_1a_RNA_DDR_expression_by_tissue_barplot")
ggsave(paste0(stem_ddr_bar, ".png"), plot_ddr_by_tissue_bar, width = 14.0, height = 10.5, dpi = 300, bg = "white", device = ragg::agg_png)
ggsave(paste0(stem_ddr_bar, ".pdf"), plot_ddr_by_tissue_bar, width = 14.0, height = 10.5, bg = "white", device = cairo_pdf)

# --- 2C. Scheme B Primary Figure 1a: DDR Expression Boxplot (31 Groups, ANOVA) ---
aov_ddr_cat <- summary(aov(DDRExpressionMedian ~ Category, data = grp_ddr_dt))[[1L]]
sub_ddr_cat <- paste0("Four-category one-way ANOVA p = ", format_q_value(aov_ddr_cat["Category", "Pr(>F)"]), " (F = ", sprintf("%.2f", aov_ddr_cat["Category", "F value"]), ")")

plot_ddr_cat_box <- ggplot(grp_ddr_dt, aes(x = X, y = DDRExpressionMedian, fill = Category)) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5), colour = "#E5E7EB", linetype = "dashed", linewidth = 0.5) +
  geom_boxplot(aes(group = Category), width = 0.58, outlier.shape = NA, colour = charcoal, linewidth = 0.55, median.linewidth = 0.75, alpha = 0.82, orientation = "x", na.rm = TRUE) +
  geom_point(position = position_jitter(width = 0.095, height = 0, seed = 25), shape = 21, size = 3.0, stroke = 0.55, colour = "white", alpha = 0.92, na.rm = TRUE) +
  scale_fill_manual(values = category_fills, guide = "none", drop = FALSE) +
  scale_x_continuous(breaks = seq_along(category_order), labels = unname(category_labels[category_order]), limits = c(0.5, length(category_order) + 0.5), expand = c(0, 0)) +
  scale_y_continuous(limits = c(1.5, 6.0), breaks = seq(2, 6, by = 1), expand = c(0, 0)) +
  labs(title = "DDR gene expression across four biological categories (31 groups, fixed panel)", subtitle = sub_ddr_cat, x = NULL, y = "DDR gene expression (median log2(qsmooth TPM + 0.5))") +
  theme_minimal(base_family = publication_font, base_size = 14) +
  theme(
    panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(), panel.grid.major.y = element_line(colour = grid_colour, linewidth = 0.50),
    axis.line.x = element_line(colour = "#8C939E", linewidth = 0.60), axis.line.y = element_line(colour = "#8C939E", linewidth = 0.60),
    axis.text.x = element_text(size = 13.5, face = "bold", colour = charcoal, lineheight = 0.95, margin = margin(t = 6)),
    axis.text.y = element_text(size = 13.0, colour = charcoal), axis.title.y = element_text(size = 14.5, face = "bold", colour = charcoal, margin = margin(r = 10)),
    plot.title = element_text(size = 15.0, face = "bold", colour = charcoal, hjust = 0.5, margin = margin(b = 4)),
    plot.subtitle = element_text(size = 12.0, colour = muted_text, hjust = 0.5, margin = margin(b = 10)),
    plot.margin = margin(12, 16, 12, 12), plot.background = element_rect(fill = "white", colour = NA)
  )

stem_ddr_cat <- file.path(out_dir, "Figure_1a_RNA_DDR_expression_31groups_boxplot")
ggsave(paste0(stem_ddr_cat, ".png"), plot_ddr_cat_box, width = 8.5, height = 7.0, dpi = 300, bg = "white", device = ragg::agg_png)
ggsave(paste0(stem_ddr_cat, ".pdf"), plot_ddr_cat_box, width = 8.5, height = 7.0, bg = "white", device = cairo_pdf)

# --- 2D. 31-Group DDR Gene Fraction Boxplot ---
aov_frac_31 <- summary(aov(DdrFractionMedian ~ Category, data = grp_ddr_dt))[[1L]]
sub_frac_31 <- paste0("Four-category one-way ANOVA p = ", format_q_value(aov_frac_31["Category", "Pr(>F)"]), " (F = ", sprintf("%.2f", aov_frac_31["Category", "F value"]), ")")

plot_frac_31 <- ggplot(grp_ddr_dt, aes(x = X, y = DdrFractionMedian, fill = Category)) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5), colour = "#E5E7EB", linetype = "dashed", linewidth = 0.5) +
  geom_boxplot(aes(group = Category), width = 0.58, outlier.shape = NA, colour = charcoal, linewidth = 0.55, median.linewidth = 0.75, alpha = 0.82, orientation = "x", na.rm = TRUE) +
  geom_point(position = position_jitter(width = 0.095, height = 0, seed = 25), shape = 21, size = 3.0, stroke = 0.50, colour = "white", alpha = 0.90, na.rm = TRUE) +
  scale_fill_manual(values = category_fills, guide = "none", drop = FALSE) +
  scale_x_continuous(breaks = seq_along(category_order), labels = unname(category_labels[category_order]), limits = c(0.5, length(category_order) + 0.5), expand = c(0, 0)) +
  scale_y_continuous(limits = c(2.0, 4.0), breaks = seq(2.0, 4.0, by = 0.5), labels = function(y) paste0(y, "%"), expand = c(0, 0)) +
  labs(title = "DDR annotated gene fraction in transcriptome (31 groups, qsmooth)", subtitle = sub_frac_31, x = NULL, y = "GO-DDR expressed gene fraction (%)") +
  theme_minimal(base_family = publication_font, base_size = 14) +
  theme(
    panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(), panel.grid.major.y = element_line(colour = grid_colour, linewidth = 0.50), panel.border = element_blank(),
    axis.line.x = element_line(colour = "#8C939E", linewidth = 0.60), axis.line.y = element_line(colour = "#8C939E", linewidth = 0.60),
    axis.text.x = element_text(size = 13.5, face = "bold", colour = charcoal, lineheight = 0.95, margin = margin(t = 6)),
    axis.text.y = element_text(size = 13.0, colour = charcoal), axis.title.y = element_text(size = 14.5, face = "bold", colour = charcoal, margin = margin(r = 10)),
    plot.title = element_text(size = 15.0, face = "bold", colour = charcoal, hjust = 0.5, margin = margin(b = 4)),
    plot.subtitle = element_text(size = 12.0, colour = muted_text, hjust = 0.5, margin = margin(b = 10)),
    plot.margin = margin(12, 16, 12, 12), plot.background = element_rect(fill = "white", colour = NA)
  )

stem_frac_31 <- file.path(out_dir, "Figure_1a_RNA_DDR_fraction_31groups_boxplot")
ggsave(paste0(stem_frac_31, ".png"), plot_frac_31, width = 8.5, height = 7.0, dpi = 300, bg = "white", device = ragg::agg_png)
ggsave(paste0(stem_frac_31, ".pdf"), plot_frac_31, width = 8.5, height = 7.0, bg = "white", device = cairo_pdf)

message(">>> Task 2 complete: Figure 1a fixed DDR figures saved.")

# ==============================================================================
# SECTION 4: LACTYLATED PROTEIN GENES — PER-TISSUE BOXPLOT & BARPLOT
# ==============================================================================
message(">>> Task 3: Calculating expression of genes corresponding to lactylated proteins across tissues...")

mapping_ddr <- fread(file.path(root, "outputs", "20260916_ddr_panel_31group", "ddr_uniprot_to_ensembl.tsv"))
kla_union_genes <- unique(mapping_ddr[grepl("KlaUnion", PanelMembership) & NGeneIDs == 1L & !is.na(EnsemblGeneIDs), EnsemblGeneIDs])
kla_union_genes <- intersect(kla_union_genes, qb_dt[[gid_col]])
message(sprintf("Lactylated protein genes (KlaUnion): %d genes present in 18,332 space", length(kla_union_genes)))

sub_kla <- as.matrix(qb_dt[match(kla_union_genes, get(gid_col)), full_samples, with = FALSE])
sample_kla_score <- colMedians(sub_kla)

sample_kla_dt <- data.table(
  SampleID = full_samples,
  KlaExpressionMedian = sample_kla_score
)
sample_kla_dt <- merge(sample_kla_dt, qc[, .(SampleID, ReferenceKey, GroupIDs)], by = "SampleID")

sample_31_kla_list <- list()
for (i in seq_len(nrow(group_exp))) {
  gid <- group_exp$GroupID[[i]]
  ref_k <- group_exp$ReferenceKey[[i]]
  sub_s <- copy(sample_kla_dt[ReferenceKey == ref_k])
  sub_s[, GroupID := gid]
  sample_31_kla_list[[i]] <- sub_s
}
sample_31_kla_dt <- rbindlist(sample_31_kla_list)
sample_31_kla_dt <- merge(sample_31_kla_dt, status[, .(GroupID, PXD, SampleGroup, Category)], by = "GroupID")
sample_31_kla_dt <- merge(sample_31_kla_dt, groups_31[, .(PXD, SampleGroup, ReferenceLabelEn, RowOrder)], by = c("PXD", "SampleGroup"))
sample_31_kla_dt[, Category := factor(Category, levels = category_order)]
sample_31_kla_dt[, CategoryLabel := factor(Category, levels = category_order, labels = unname(category_labels[category_order]))]

grp_kla_dt <- sample_31_kla_dt[, .(
  KlaExpressionMedian = median(KlaExpressionMedian),
  KlaExpressionMean = mean(KlaExpressionMedian),
  KlaExpressionSD = sd(KlaExpressionMedian),
  KlaExpressionSEM = sd(KlaExpressionMedian) / sqrt(.N),
  KlaExpressionIQR = IQR(KlaExpressionMedian),
  N = .N
), by = .(GroupID, ReferenceKey, PXD, SampleGroup, Category, CategoryLabel, ReferenceLabelEn, RowOrder)]
grp_kla_dt[, X := match(Category, category_order)]

fwrite(sample_31_kla_dt, file.path(out_dir, "figure1a_lactylated_gene_sample_values.csv"))
fwrite(grp_kla_dt, file.path(out_dir, "figure1a_lactylated_gene_31groups.csv"))

# --- 4A. Per-Tissue Lactylated Gene Expression Boxplot with Sample Points ---
setorder(grp_kla_dt, Category, KlaExpressionMedian)
grp_kla_dt[, DisplayLabel := sprintf("%s  %s", GroupID, ReferenceLabelEn)]
ordered_kla_labels <- grp_kla_dt$DisplayLabel
grp_kla_dt[, FactorLabel := factor(DisplayLabel, levels = ordered_kla_labels)]

sample_31_kla_dt <- merge(sample_31_kla_dt, grp_kla_dt[, .(GroupID, FactorLabel)], by = "GroupID")
sample_31_kla_dt[, FactorLabel := factor(FactorLabel, levels = ordered_kla_labels)]

grp_kla_summary <- sample_31_kla_dt[, .(
  N = .N,
  Median = median(KlaExpressionMedian),
  Mean = mean(KlaExpressionMedian),
  Q25 = quantile(KlaExpressionMedian, 0.25),
  Q75 = quantile(KlaExpressionMedian, 0.75),
  MaxScore = max(KlaExpressionMedian)
), by = .(GroupID, FactorLabel, Category, CategoryLabel)]

x_max_kla <- 6.5

plot_kla_by_tissue_box <- ggplot(sample_31_kla_dt, aes(x = KlaExpressionMedian, y = FactorLabel, fill = Category)) +
  geom_boxplot(aes(group = FactorLabel), width = 0.58, outlier.shape = NA, colour = charcoal, linewidth = 0.50, median.linewidth = 0.75, alpha = 0.80) +
  geom_point(position = position_jitter(width = 0, height = 0.18, seed = 25), shape = 21, size = 1.9, stroke = 0.30, colour = "white", alpha = 0.65) +
  geom_point(data = grp_kla_summary, aes(x = Mean, y = FactorLabel), inherit.aes = FALSE, shape = 23, size = 2.4, fill = mean_colour, colour = charcoal, stroke = 0.45) +
  geom_text(data = grp_kla_summary, aes(x = x_max_kla * 0.94, y = FactorLabel, label = paste0("n=", N)), inherit.aes = FALSE, family = publication_font, size = 3.6, fontface = "bold", colour = muted_text, hjust = 0) +
  scale_fill_manual(values = category_fills, labels = c(normal_tissue = "non-tumor tissues", cancer_tissue = "tumor tissues", normal_cells = "normal cell lines", cancer_cells = "cancer cell lines"), name = "Biological category") +
  scale_x_continuous(limits = c(1.5, x_max_kla), breaks = seq(2.0, 6.0, by = 1.0), expand = c(0.01, 0)) +
  labs(
    title = "Transcriptomic expression distribution of lactylated protein genes across 31 materials",
    subtitle = paste0("Sample median log2(qsmooth TPM + 0.5) across ", length(kla_union_genes), " lactylated protein genes | Red diamonds = group mean | n = sample count"),
    x = "Lactylated protein gene expression (sample median log2(qsmooth TPM + 0.5))", y = NULL
  ) +
  theme_minimal(base_family = publication_font, base_size = 12) +
  theme(
    panel.grid.minor = element_blank(), panel.grid.major.y = element_line(colour = "#ECEFF1", linewidth = 0.45), panel.grid.major.x = element_line(colour = grid_colour, linewidth = 0.50),
    axis.text.y = element_text(size = 9.5, colour = charcoal), axis.text.x = element_text(size = 11.0, colour = charcoal),
    axis.title.x = element_text(size = 12.0, face = "bold", colour = charcoal, margin = margin(t = 8)),
    legend.position = "top", legend.direction = "horizontal", legend.title = element_text(size = 11.0, face = "bold", colour = charcoal), legend.text = element_text(size = 10.5, colour = charcoal),
    plot.title = element_text(size = 14.0, face = "bold", colour = charcoal, hjust = 0), plot.subtitle = element_text(size = 10.5, colour = muted_text, hjust = 0, margin = margin(b = 10)),
    plot.margin = margin(12, 20, 12, 12), plot.background = element_rect(fill = "white", colour = NA)
  )

stem_kla_box <- file.path(out_dir, "Figure_1a_RNA_lactylated_gene_expression_by_tissue_boxplot")
ggsave(paste0(stem_kla_box, ".png"), plot_kla_by_tissue_box, width = 14.0, height = 11.5, dpi = 300, bg = "white", device = ragg::agg_png)
ggsave(paste0(stem_kla_box, ".pdf"), plot_kla_by_tissue_box, width = 14.0, height = 11.5, bg = "white", device = cairo_pdf)

# --- 4B. Per-Tissue Lactylated Gene Expression Barplot (Ranked) ---
grp_kla_bar <- copy(grp_kla_dt)
setorder(grp_kla_bar, KlaExpressionMedian)
grp_kla_bar[, BarLabel := factor(DisplayLabel, levels = DisplayLabel)]

plot_kla_by_tissue_bar <- ggplot(grp_kla_bar, aes(x = KlaExpressionMedian, y = BarLabel, fill = Category)) +
  geom_col(width = 0.72, alpha = 0.88, colour = charcoal, linewidth = 0.25) +
  geom_errorbar(aes(xmin = pmax(0, KlaExpressionMedian - KlaExpressionIQR / 2), xmax = KlaExpressionMedian + KlaExpressionIQR / 2), width = 0.35, colour = charcoal, linewidth = 0.45) +
  geom_text(aes(label = sprintf("%.2f (n=%d)", KlaExpressionMedian, N)), hjust = -0.15, size = 3.3, family = publication_font, colour = charcoal) +
  scale_fill_manual(values = category_fills, labels = c("non-tumor tissues", "tumor tissues", "normal cell lines", "cancer cell lines"), name = "Biological category") +
  scale_x_continuous(limits = c(0, 6.5), breaks = seq(0, 6, by = 1), expand = expansion(mult = c(0, 0.08))) +
  labs(
    title = "Transcriptomic expression ranking of lactylated protein genes across 31 materials",
    subtitle = paste0("Group median log2(qsmooth TPM + 0.5) across ", length(kla_union_genes), " lactylated protein genes | Error bars = IQR"),
    x = "Lactylated protein gene expression (median log2(qsmooth TPM + 0.5))", y = NULL
  ) +
  theme_minimal(base_family = publication_font, base_size = 12) +
  theme(
    panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), panel.grid.major.x = element_line(colour = grid_colour, linewidth = 0.50),
    axis.text.y = element_text(size = 9.5, colour = charcoal), axis.text.x = element_text(size = 10.5, colour = charcoal),
    axis.title.x = element_text(size = 12.0, face = "bold", colour = charcoal, margin = margin(t = 6)),
    legend.position = "top", legend.direction = "horizontal", legend.title = element_text(size = 11.0, face = "bold", colour = charcoal), legend.text = element_text(size = 10.5, colour = charcoal),
    plot.title = element_text(size = 14.5, face = "bold", colour = charcoal, hjust = 0), plot.subtitle = element_text(size = 10.5, colour = muted_text, hjust = 0, margin = margin(b = 8)),
    plot.margin = margin(12, 20, 12, 12)
  )

stem_kla_bar <- file.path(out_dir, "Figure_1a_RNA_lactylated_gene_expression_by_tissue_barplot")
ggsave(paste0(stem_kla_bar, ".png"), plot_kla_by_tissue_bar, width = 14.0, height = 10.5, dpi = 300, bg = "white", device = ragg::agg_png)
ggsave(paste0(stem_kla_bar, ".pdf"), plot_kla_by_tissue_bar, width = 14.0, height = 10.5, bg = "white", device = cairo_pdf)

# --- 4C. Primary Figure 1a: DDR vs Lactylated Gene Expression Comparison (31 Groups, Two-Way ANOVA) ---
grp_dual <- merge(grp_ddr_dt[, .(GroupID, Category, CategoryLabel, X, DDRExpressionMedian)],
                  grp_kla_dt[, .(GroupID, KlaExpressionMedian)], by = "GroupID")

grp_expr_long <- melt(
  grp_dual,
  id.vars = c("GroupID", "Category", "CategoryLabel", "X"),
  measure.vars = c("DDRExpressionMedian", "KlaExpressionMedian"),
  variable.name = "GeneSet", value.name = "MedianExpression"
)
grp_expr_long[, GeneSet := factor(GeneSet, levels = c("DDRExpressionMedian", "KlaExpressionMedian"),
                                  labels = c("DDR genes", "Lactylated protein genes"))]

expr_aov_31 <- aov(MedianExpression ~ Category * GeneSet, data = grp_expr_long)
expr_aov_tab_31 <- summary(expr_aov_31)[[1L]]
cat_p_31 <- expr_aov_tab_31["Category", "Pr(>F)"]
cat_f_31 <- expr_aov_tab_31["Category", "F value"]
expr_sub_31 <- paste0("Two-way ANOVA Category factor p = ", format_q_value(cat_p_31), " (F = ", sprintf("%.2f", cat_f_31), ")")

dodge_w <- 0.72
grp_expr_summary <- grp_expr_long[, .(
  Mean = mean(MedianExpression),
  Median = median(MedianExpression),
  N = .N
), by = .(Category, CategoryLabel, GeneSet)]
grp_expr_summary[, CatIdx := as.numeric(CategoryLabel)]
grp_expr_summary[, XPos := fifelse(GeneSet == "DDR genes", CatIdx - dodge_w / 4, CatIdx + dodge_w / 4)]

y_max_e31 <- max(grp_expr_long$MedianExpression) + 0.85
y_min_e31 <- min(grp_expr_long$MedianExpression) - 0.2

plot_fig1a_dual <- ggplot(grp_expr_long, aes(x = CategoryLabel, y = MedianExpression, fill = GeneSet)) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5), colour = "#E5E7EB", linetype = "dashed", linewidth = 0.5) +
  geom_boxplot(position = position_dodge(width = dodge_w), width = 0.58, outlier.shape = NA, colour = charcoal, linewidth = 0.55, median.linewidth = 0.75, alpha = 0.82, na.rm = TRUE) +
  geom_segment(data = grp_expr_summary, aes(x = XPos - 0.16, xend = XPos + 0.16, y = Mean, yend = Mean), inherit.aes = FALSE, colour = mean_colour, linewidth = 0.80) +
  geom_point(aes(fill = GeneSet), position = position_jitterdodge(jitter.width = 0.14, dodge.width = dodge_w, seed = 25), shape = 21, size = 2.8, stroke = 0.50, colour = "white", alpha = 0.90, na.rm = TRUE) +
  geom_text(data = grp_expr_summary, aes(x = XPos, y = y_max_e31 * 0.96, label = paste0("n=", N)), inherit.aes = FALSE, size = 3.8, family = publication_font, colour = muted_text, fontface = "bold") +
  scale_fill_manual(values = c("DDR genes" = "#4E79A7", "Lactylated protein genes" = "#F28E2B")) +
  scale_y_continuous(limits = c(y_min_e31, y_max_e31), breaks = scales::pretty_breaks(n = 6), expand = expansion(mult = c(0, 0))) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE, keyheight = grid::unit(0.55, "cm"), keywidth = grid::unit(0.85, "cm"))) +
  labs(
    title = "Expression of DDR genes and lactylated protein genes across four categories (31 groups, qsmooth)",
    subtitle = expr_sub_31,
    x = NULL, y = "Group median log2(qsmooth TPM + 0.5)", fill = NULL
  ) +
  theme_minimal(base_size = 14, base_family = publication_font) +
  theme(
    panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(), panel.grid.major.y = element_line(colour = grid_colour, linewidth = 0.50), panel.border = element_blank(),
    axis.line.x = element_line(colour = "#8C939E", linewidth = 0.60), axis.line.y = element_line(colour = "#8C939E", linewidth = 0.60),
    axis.text.x = element_text(size = 13.0, colour = charcoal, face = "bold", margin = margin(t = 6)),
    axis.text.y = element_text(size = 13.0, colour = charcoal), axis.title.y = element_text(size = 15.0, face = "bold", colour = charcoal, margin = margin(r = 10)),
    plot.title = element_text(size = 14.5, face = "bold", colour = charcoal, hjust = 0.5, margin = margin(b = 4)),
    plot.subtitle = element_text(size = 12.0, colour = muted_text, hjust = 0.5, margin = margin(b = 10)),
    legend.position = "top", legend.direction = "horizontal", legend.text = element_text(size = 13.0, colour = charcoal),
    legend.key.spacing.x = grid::unit(0.35, "cm"), legend.background = element_rect(fill = "white", colour = NA), legend.margin = margin(1, 0, 4, 0),
    plot.margin = margin(10, 16, 12, 12), plot.background = element_rect(fill = "white", colour = NA)
  )

stem_dual_31 <- file.path(out_dir, "Figure_1a_RNA_DDR_and_lactylated_gene_expression_31groups_boxplot")
ggsave(paste0(stem_dual_31, ".png"), plot_fig1a_dual, width = 8.5, height = 7.0, dpi = 300, bg = "white", device = ragg::agg_png)
ggsave(paste0(stem_dual_31, ".pdf"), plot_fig1a_dual, width = 8.5, height = 7.0, bg = "white", device = cairo_pdf)

message(">>> Task 3 complete: Figure 1a lactylated protein gene expression figures saved.")

# ==============================================================================
# SECTION 5: LACTYLATION REGULATORS — QSMOOTH HEATMAP (FIGURE 3 STYLE)
# ==============================================================================
message(">>> Task 4: Building RNA regulator qsmooth heatmaps across 31 groups (Figure 3 style)...")

role_order <- c("Writer", "Eraser", "Writer-Eraser", "Reader")

role_annot_file <- file.path(root, "data", "publication_input", "Supplementary_Table_S5_Lactylation_Regulators.xlsx")
stopifnot(file.exists(role_annot_file))

s_tab <- as.data.table(readxl::read_excel(role_annot_file, sheet = "Regulator_Annotations"))
role_map_std <- unique(s_tab[Role %in% role_order & nzchar(BaseAccession), .(
  Role = trimws(as.character(Role)),
  GeneSymbol = trimws(as.character(GeneSymbol)),
  BaseAccession = trimws(as.character(BaseAccession))
)])
role_map_std[, DisplayName := ifelse(BaseAccession == "Q92830", "GCN5 (KAT2A)", GeneSymbol)]
role_map_std[, Role := factor(Role, levels = role_order)]

reg_map <- merge(role_map_std, mapping_ddr[, .(BaseAccession, EnsemblGeneIDs, EnsemblGeneViaEntrez)],
                 by = "BaseAccession", all.x = TRUE)

all_qa_genes <- qa_dt[[qa_gid_col]]

resolve_ensg <- function(ensg_str, entrez_ensg_str) {
  candidates <- c(strsplit(as.character(ensg_str), ";")[[1L]], strsplit(as.character(entrez_ensg_str), ";")[[1L]])
  candidates <- candidates[nzchar(candidates) & !is.na(candidates)]
  in_qa <- candidates[candidates %in% all_qa_genes]
  if (length(in_qa)) in_qa[1L] else candidates[1L]
}
reg_map[, EnsemblGeneID := mapply(resolve_ensg, EnsemblGeneIDs, EnsemblGeneViaEntrez)]
stopifnot(all(reg_map$EnsemblGeneID %in% all_qa_genes))
message("Verified: All 48 unique regulators (49 role entries) present natively in qsmooth_A!")

# Expand across 31 groups from qsmooth_A (0 imputation, 0 fallback)
reg_qsmooth_list <- list()
for (i in seq_len(nrow(group_exp))) {
  gid <- group_exp$GroupID[[i]]
  ref_k <- group_exp$ReferenceKey[[i]]
  
  for (j in seq_len(nrow(reg_map))) {
    acc <- reg_map$BaseAccession[[j]]
    role_val <- as.character(reg_map$Role[[j]])
    disp <- reg_map$DisplayName[[j]]
    sym <- reg_map$GeneSymbol[[j]]
    ensg <- reg_map$EnsemblGeneID[[j]]
    val <- qa_dt[get(qa_gid_col) == ensg, get(ref_k)]
    
    reg_qsmooth_list[[length(reg_qsmooth_list) + 1L]] <- data.table(
      GroupID = gid,
      ReferenceKey = ref_k,
      RegulatorBaseAccession = acc,
      Role = role_val,
      DisplayName = disp,
      GeneSymbol = sym,
      EnsemblGeneID = ensg,
      QsmoothLog2TPM = val
    )
  }
}
reg_qsmooth_dt <- rbindlist(reg_qsmooth_list)
reg_qsmooth_dt <- merge(reg_qsmooth_dt, status[, .(GroupID, PXD, SampleGroup, Category)], by = "GroupID")
reg_qsmooth_dt <- merge(reg_qsmooth_dt, groups_31[, .(PXD, SampleGroup, ReferenceLabelEn, RowOrder)], by = c("PXD", "SampleGroup"))

# Calculate Z-score per gene across 31 groups
reg_qsmooth_dt[, QsmoothZScore := (QsmoothLog2TPM - mean(QsmoothLog2TPM)) / sd(QsmoothLog2TPM), by = .(DisplayName, Role)]

fwrite(reg_qsmooth_dt, file.path(out_dir, "regulator_rna_qsmooth_31.csv"))

# Factor levels
reg_qsmooth_dt[, CategoryLabel := factor(as.character(Category), levels = category_order, labels = unname(heatmap_category_labels[category_order]))]
reg_qsmooth_dt[, RoleLabel := factor(Role, levels = role_order)]
reg_qsmooth_dt[, PlotLabel := factor(ReferenceLabelEn, levels = rev(unique(ReferenceLabelEn[order(RowOrder)])))]
reg_qsmooth_dt[, DisplayName := factor(DisplayName, levels = unique(role_map_std$DisplayName))]

highlight_genes <- c("AARS1", "ACAT2", "KRT18", "SIRT2", "PARK7", "HDAC1", "HDAC2", "BRD4", "SMARCA4", "TRIM33")

# Continuous bounding box coordinates
box_lines <- list()
for (cat_lbl in levels(reg_qsmooth_dt$CategoryLabel)) {
  sub_cat <- reg_qsmooth_dt[CategoryLabel == cat_lbl]
  n_rows <- uniqueN(sub_cat$PlotLabel)
  for (role_lbl in levels(reg_qsmooth_dt$RoleLabel)) {
    sub_panel <- sub_cat[RoleLabel == role_lbl]
    if (nrow(sub_panel) == 0L) next
    panel_genes <- levels(droplevels(sub_panel$DisplayName))
    for (g in highlight_genes) {
      if (g %in% panel_genes) {
        x_pos <- which(panel_genes == g)
        box_lines[[length(box_lines) + 1L]] <- data.frame(CategoryLabel = cat_lbl, RoleLabel = role_lbl, x = x_pos - 0.5, xend = x_pos - 0.5, y = 0.5, yend = n_rows + 0.5)
        box_lines[[length(box_lines) + 1L]] <- data.frame(CategoryLabel = cat_lbl, RoleLabel = role_lbl, x = x_pos + 0.5, xend = x_pos + 0.5, y = 0.5, yend = n_rows + 0.5)
        box_lines[[length(box_lines) + 1L]] <- data.frame(CategoryLabel = cat_lbl, RoleLabel = role_lbl, x = x_pos - 0.5, xend = x_pos + 0.5, y = 0.5, yend = 0.5)
        box_lines[[length(box_lines) + 1L]] <- data.frame(CategoryLabel = cat_lbl, RoleLabel = role_lbl, x = x_pos - 0.5, xend = x_pos + 0.5, y = n_rows + 0.5, yend = n_rows + 0.5)
      }
    }
  }
}
box_lines_df <- rbindlist(box_lines)
if (nrow(box_lines_df) > 0L) {
  box_lines_df[, CategoryLabel := factor(CategoryLabel, levels = levels(reg_qsmooth_dt$CategoryLabel))]
  box_lines_df[, RoleLabel := factor(RoleLabel, levels = levels(reg_qsmooth_dt$RoleLabel))]
}

# --- 5A. Direct Qsmooth Log2TPM Heatmap ---
rna_teal_palette <- c("#FFFFFF", "#E0F2F1", "#80CBC4", "#26A69A", "#00897B", "#004D40")
box_colour_log2 <- "#C0392B" # Crimson red for high contrast against teal

build_qsmooth_heatmap <- function(include_boxes) {
  p <- ggplot(reg_qsmooth_dt, aes(x = DisplayName, y = PlotLabel, fill = QsmoothLog2TPM)) +
    geom_tile(colour = "white", linewidth = 0.22)
  if (include_boxes && nrow(box_lines_df) > 0L) {
    p <- p + geom_segment(
      data = box_lines_df,
      aes(x = x, xend = xend, y = y, yend = yend),
      inherit.aes = FALSE,
      colour = box_colour_log2, linewidth = 1.15
    )
  }
  p +
    facet_grid(CategoryLabel ~ RoleLabel, scales = "free", space = "free") +
    scale_fill_gradientn(
      colours = rna_teal_palette,
      na.value = "#D9D9D9",
      name = "qsmooth log2(TPM + 0.5)",
      guide = guide_colourbar(
        title.position = "top",
        title.hjust = 0,
        barwidth = grid::unit(5.0, "mm"),
        barheight = grid::unit(72, "mm")
      )
    ) +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 10.5, base_family = publication_font) +
    theme(
      panel.grid = element_blank(),
      strip.text.x = element_text(face = "bold", size = 12),
      strip.text.y.right = element_text(face = "bold", size = 11.0, angle = 90, hjust = 0.5),
      strip.background = element_rect(fill = "#F2F2F2", colour = NA),
      axis.text.x = element_text(angle = 55, hjust = 1, vjust = 1, size = 10.2),
      axis.text.y = element_text(size = 9.6),
      legend.position = "right",
      legend.direction = "vertical",
      legend.title = element_text(size = 11.0, face = "bold", margin = margin(b = 8)),
      legend.text = element_text(size = 9.8),
      legend.margin = margin(0, 8, 0, 8),
      plot.margin = margin(10, 14, 10, 10)
    )
}

p_qsm_framed <- build_qsmooth_heatmap(include_boxes = TRUE)
p_qsm_unboxed <- build_qsmooth_heatmap(include_boxes = FALSE)

ggsave(file.path(out_dir, "Figure_3c_RNA_regulator_qsmooth_heatmap.png"), p_qsm_framed, width = 16.5, height = 11.5, dpi = 300, bg = "white", device = ragg::agg_png)
ggsave(file.path(out_dir, "Figure_3c_RNA_regulator_qsmooth_heatmap.pdf"), p_qsm_framed, width = 16.5, height = 11.5, bg = "white", device = cairo_pdf)
ggsave(file.path(out_dir, "Figure_3c_RNA_regulator_qsmooth_heatmap_no_frame.png"), p_qsm_unboxed, width = 16.5, height = 11.5, dpi = 300, bg = "white", device = ragg::agg_png)
ggsave(file.path(out_dir, "Figure_3c_RNA_regulator_qsmooth_heatmap_no_frame.pdf"), p_qsm_unboxed, width = 16.5, height = 11.5, bg = "white", device = cairo_pdf)

# --- 5B. Standardized Qsmooth Z-Score Heatmap ---
z_palette <- c("#2166AC", "#4393C3", "#92C5DE", "#F7F7F7", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B")
box_colour_z <- "#1B5E20" # Dark green border for clear contrast on blue-white-red

build_qsmooth_z_heatmap <- function(include_boxes) {
  p <- ggplot(reg_qsmooth_dt, aes(x = DisplayName, y = PlotLabel, fill = pmin(2.5, pmax(-2.5, QsmoothZScore)))) +
    geom_tile(colour = "white", linewidth = 0.22)
  if (include_boxes && nrow(box_lines_df) > 0L) {
    p <- p + geom_segment(
      data = box_lines_df,
      aes(x = x, xend = xend, y = y, yend = yend),
      inherit.aes = FALSE,
      colour = box_colour_z, linewidth = 1.15
    )
  }
  p +
    facet_grid(CategoryLabel ~ RoleLabel, scales = "free", space = "free") +
    scale_fill_gradientn(
      colours = z_palette,
      limits = c(-2.5, 2.5),
      na.value = "#D9D9D9",
      name = "Relative expression\n(gene Z-score)",
      guide = guide_colourbar(
        title.position = "top",
        title.hjust = 0,
        barwidth = grid::unit(5.0, "mm"),
        barheight = grid::unit(72, "mm")
      )
    ) +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 10.5, base_family = publication_font) +
    theme(
      panel.grid = element_blank(),
      strip.text.x = element_text(face = "bold", size = 12),
      strip.text.y.right = element_text(face = "bold", size = 11.0, angle = 90, hjust = 0.5),
      strip.background = element_rect(fill = "#F2F2F2", colour = NA),
      axis.text.x = element_text(angle = 55, hjust = 1, vjust = 1, size = 10.2),
      axis.text.y = element_text(size = 9.6),
      legend.position = "right",
      legend.direction = "vertical",
      legend.title = element_text(size = 11.0, face = "bold", margin = margin(b = 8)),
      legend.text = element_text(size = 9.8),
      legend.margin = margin(0, 8, 0, 8),
      plot.margin = margin(10, 14, 10, 10)
    )
}

p_z_framed <- build_qsmooth_z_heatmap(include_boxes = TRUE)
p_z_unboxed <- build_qsmooth_z_heatmap(include_boxes = FALSE)

ggsave(file.path(out_dir, "Figure_3c_RNA_regulator_qsmooth_zscore_heatmap.png"), p_z_framed, width = 16.5, height = 11.5, dpi = 300, bg = "white", device = ragg::agg_png)
ggsave(file.path(out_dir, "Figure_3c_RNA_regulator_qsmooth_zscore_heatmap.pdf"), p_z_framed, width = 16.5, height = 11.5, bg = "white", device = cairo_pdf)
ggsave(file.path(out_dir, "Figure_3c_RNA_regulator_qsmooth_zscore_heatmap_no_frame.png"), p_z_unboxed, width = 16.5, height = 11.5, dpi = 300, bg = "white", device = ragg::agg_png)
ggsave(file.path(out_dir, "Figure_3c_RNA_regulator_qsmooth_zscore_heatmap_no_frame.pdf"), p_z_unboxed, width = 16.5, height = 11.5, bg = "white", device = cairo_pdf)

message(">>> Task 4 complete: Figure 3c RNA regulator qsmooth heatmaps saved.")

# ==============================================================================
# SECTION 6: SYNCHRONIZATION TO DELIVERY DIRECTORIES
# ==============================================================================
message(">>> Synchronizing all newly generated deliverables...")

# 1. Local delivery dir
deliv_dir <- file.path(root, "outputs", "20260917_delivery_31group", "rna")
if (dir.exists(deliv_dir)) {
  file.copy(list.files(out_dir, full.names = TRUE), deliv_dir, overwrite = TRUE)
}

# 2. External desktop delivery dirs
ext_renew_top <- "/Users/gzy2520/Desktop/renew/kla"
ext_renew_rna_upper <- file.path(ext_renew_top, "RNA")
ext_renew_rna_lower <- file.path(ext_renew_top, "rna")

if (dir.exists(ext_renew_top)) {
  # Deliver all files (PDF, PNG, CSV) to RNA subfolder
  dir.create(ext_renew_rna_upper, recursive = TRUE, showWarnings = FALSE)
  file.copy(list.files(out_dir, full.names = TRUE), ext_renew_rna_upper, overwrite = TRUE)
  if (dir.exists(ext_renew_rna_lower)) {
    file.copy(list.files(out_dir, full.names = TRUE), ext_renew_rna_lower, overwrite = TRUE)
  }
  
  # Copy primary publication PNGs to ext_renew_top
  top_png_names <- c(
    # 1. Proliferation (Official Hallmark G2M Checkpoint)
    "Figure_1b_RNA_proliferation_hallmark_by_tissue_boxplot.png",
    "Figure_1b_RNA_proliferation_hallmark_by_tissue_barplot.png",
    "Figure_1b_RNA_proliferation_hallmark_31groups_boxplot.png",
    # 2. Fixed DDR Panel (371 genes)
    "Figure_1a_RNA_DDR_expression_by_tissue_boxplot.png",
    "Figure_1a_RNA_DDR_expression_by_tissue_barplot.png",
    "Figure_1a_RNA_DDR_expression_31groups_boxplot.png",
    "Figure_1a_RNA_DDR_fraction_31groups_boxplot.png",
    # 3. Lactylation Regulators Heatmaps (qsmooth)
    "Figure_3c_RNA_regulator_qsmooth_heatmap.png",
    "Figure_3c_RNA_regulator_qsmooth_heatmap_no_frame.png",
    "Figure_3c_RNA_regulator_qsmooth_zscore_heatmap.png",
    "Figure_3c_RNA_regulator_qsmooth_zscore_heatmap_no_frame.png",
    # 4. Lactylated Protein Genes Expression (5,224 genes)
    "Figure_1a_RNA_lactylated_gene_expression_by_tissue_boxplot.png",
    "Figure_1a_RNA_lactylated_gene_expression_by_tissue_barplot.png",
    "Figure_1a_RNA_DDR_and_lactylated_gene_expression_31groups_boxplot.png"
  )
  for (tp in top_png_names) {
    src_tp <- file.path(out_dir, tp)
    if (file.exists(src_tp)) file.copy(src_tp, ext_renew_top, overwrite = TRUE)
  }
}

cat("\nALL_4_REQUIRED_DELIVERABLES_SUCCESSFULLY_GENERATED_AND_DELIVERED!\n")
