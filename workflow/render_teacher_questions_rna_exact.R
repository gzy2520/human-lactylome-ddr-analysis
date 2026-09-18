#!/usr/bin/env Rscript
# ==============================================================================
# Master Script: Render RNA-seq Counterparts for Proteome Core Figures
# (Teacher Questions Revision — Cross-Tissue qsmooth Normalization & Per-Tissue Distribution)
#
# Tasks & User Requirements Addressed:
#   1. Cross-Tissue Normalization:
#      - Strictly utilizes qsmooth-normalized matrices:
#        * Sample level: outputs/20260916_qsmooth_31group/matrices/qsmooth_B_full_log2tpm.tsv.gz
#        * Group level: outputs/20260916_qsmooth_31group/matrices/qsmooth_A_collapsed_log2tpm.tsv.gz
#   2. Cell Proliferation Rate (20-gene hallmark score):
#      - Per-tissue distribution across all 31 biological materials with individual sample points
#        (Horizontal comprehensive boxplot + sample jitter dots + sample size n labeled).
#      - Faceted 4-category per-tissue companion boxplot.
#      - Scheme B 4-category primary boxplot (matching Figure 1b, ANOVA p & F).
#      - 31-group ranking chart and sample-level 4-category companion.
#      - MKI67 ratio companions (MKI67/H3C1, MKI67/ACTB, MKI67/TUBB).
#   3. Lactylation Regulators Heatmap (Figure 3c):
#      - Replaces relative percentiles with qsmooth normalized expression!
#      - Generates both:
#        * qsmooth log2(TPM + 0.5) expression heatmap (continuous gradient)
#        * gene-standardized Z-score of qsmooth expression heatmap (divergent blue-white-red)
#      - 48 regulators grouped into 4 functional roles (Writer, Eraser, Writer-Eraser, Reader).
#      - 31 materials grouped into 4 categories (non-tumor tissues, tumor tissues, cancer cells, normal cells).
#      - Continuous bounding boxes on 10 highlighted genes (AARS1, ACAT2, KRT18, SIRT2,
#        PARK7, HDAC1, HDAC2, BRD4, SMARCA4, TRIM33), plus unboxed versions.
#      - Ensures AARS1 and CSRP2BP have valid qsmooth values.
#   4. DDR & Lactylated Gene Expression:
#      - Dedicated per-tissue DDR gene expression distribution plot with sample points across 31 groups.
#      - 31-group Scheme B DDR vs Lactylated gene expression comparison (qsmooth, Two-way ANOVA).
#      - Sample-level DDR vs Lactylated gene expression companion (qsmooth).
#      - 31-group DDR expressed gene fraction boxplot & per-tissue DDR fraction plot.
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
# SECTION 1: LOAD METADATA AND QSMOOTH NORMALIZED MATRICES
# ==============================================================================
message(">>> Loading sample metadata, 31-group status, and qsmooth expression matrices...")

qc <- fread(file.path(root, "outputs", "20260916_expression_extraction", "group_sample_qc.csv"))
status <- fread(file.path(root, "audit", "20260916_full_31_rna_status", "rna_31_group_status.csv"))
groups_31 <- fread(file.path(root, "data", "publication_input", "group_summary_31.csv"))
group_exp <- fread(file.path(root, "outputs", "20260916_qsmooth_31group", "group_expansion_31.csv"))

group2cat <- setNames(status$Category, status$GroupID)
qc[, Category := vapply(strsplit(GroupIDs, ";", fixed = TRUE), function(g) {
  match_g <- intersect(g, names(group2cat))
  if (length(match_g)) group2cat[[match_g[1L]]] else NA_character_
}, character(1))]
qc[, Category := factor(Category, levels = category_order)]

# 1. Sample-level qsmooth-normalized expression matrix (17,340 genes x 1,898 samples)
qsmooth_b_path <- file.path(root, "outputs", "20260916_qsmooth_31group", "matrices", "qsmooth_B_full_log2tpm.tsv.gz")
stopifnot(file.exists(qsmooth_b_path))
qb_dt <- read_table_gz(qsmooth_b_path)
gid_col <- names(qb_dt)[1L]
full_samples <- setdiff(names(qb_dt), gid_col)
message(sprintf("Loaded qsmooth_B: %d genes x %d samples", nrow(qb_dt), length(full_samples)))

# 2. Group/reference-level qsmooth-normalized collapsed matrix (17,340 genes x 28 references)
qsmooth_a_path <- file.path(root, "outputs", "20260916_qsmooth_31group", "matrices", "qsmooth_A_collapsed_log2tpm.tsv.gz")
stopifnot(file.exists(qsmooth_a_path))
qa_dt <- read_table_gz(qsmooth_a_path)
qa_gid_col <- names(qa_dt)[1L]
message(sprintf("Loaded qsmooth_A: %d genes x %d reference columns", nrow(qa_dt), ncol(qa_dt) - 1L))

# ==============================================================================
# SECTION 2: PROLIFERATION RATE — QSMOOTH MATRIX & PER-TISSUE SAMPLE DISTRIBUTION
# ==============================================================================
message(">>> Task 1: Building Proliferation Rate figures (Per-Tissue Sample Distribution & Scheme B)...")

# 20 canonical proliferation marker genes covering S, G2, M phases
prolif_markers <- c(
  MKI67 = "ENSG00000148773", PCNA  = "ENSG00000132646", TOP2A = "ENSG00000131747",
  MCM2  = "ENSG00000073111", MCM3  = "ENSG00000112118", MCM4  = "ENSG00000104738",
  MCM5  = "ENSG00000100297", MCM6  = "ENSG00000076003", MCM7  = "ENSG00000078898",
  CCNB1 = "ENSG00000134057", CCNB2 = "ENSG00000157456", CDK1  = "ENSG00000170312",
  BUB1  = "ENSG00000169679", AURKA = "ENSG00000087586", AURKB = "ENSG00000178999",
  PLK1  = "ENSG00000166851", RRM2  = "ENSG00000171848", TYMS  = "ENSG00000176890",
  BIRC5 = "ENSG00000089685", UBE2C = "ENSG00000175063"
)

# Extract from qsmooth_B
sub_prolif <- as.matrix(qb_dt[match(prolif_markers, get(gid_col)), full_samples, with = FALSE])
rownames(sub_prolif) <- names(prolif_markers)

sample_prolif_score <- colMeans(sub_prolif)
sample_prolif_z <- colMeans((sub_prolif - rowMeans(sub_prolif)) / rowSds(sub_prolif))

sample_prolif_dt <- data.table(
  SampleID = full_samples,
  ProliferationScore = sample_prolif_score,
  ProliferationZScore = sample_prolif_z
)
sample_prolif_dt <- merge(sample_prolif_dt, qc[, .(SampleID, Category, ReferenceKey, GroupIDs)], by = "SampleID")
sample_prolif_dt[, Category := factor(Category, levels = category_order)]

fwrite(sample_prolif_dt, file.path(out_dir, "figure1b_rna_proliferation_score_sample_values.csv"))

# Map each sample to its corresponding 31 Group row
sample_31_prolif_list <- list()
for (i in seq_len(nrow(group_exp))) {
  gid <- group_exp$GroupID[[i]]
  ref_k <- group_exp$ReferenceKey[[i]]
  sub_s <- copy(sample_prolif_dt[ReferenceKey == ref_k])
  sub_s[, Category := NULL]
  sub_s[, GroupID := gid]
  sample_31_prolif_list[[i]] <- sub_s
}
sample_31_prolif_dt <- rbindlist(sample_31_prolif_list)
sample_31_prolif_dt <- merge(sample_31_prolif_dt, status[, .(GroupID, PXD, SampleGroup, Category)], by = "GroupID")
sample_31_prolif_dt <- merge(sample_31_prolif_dt, groups_31[, .(PXD, SampleGroup, ReferenceLabelEn, RowOrder)], by = c("PXD", "SampleGroup"))
sample_31_prolif_dt[, Category := factor(Category, levels = category_order)]
sample_31_prolif_dt[, CategoryLabel := factor(Category, levels = category_order, labels = unname(category_labels[category_order]))]

# Group summary (Scheme B: median of each material)
grp_prolif_dt <- sample_31_prolif_dt[, .(
  ProliferationScore = median(ProliferationScore),
  ProliferationMean = mean(ProliferationScore),
  ProliferationZScore = median(ProliferationZScore),
  N = .N
), by = .(GroupID, ReferenceKey, PXD, SampleGroup, Category, CategoryLabel, ReferenceLabelEn, RowOrder)]
grp_prolif_dt[, X := match(Category, category_order)]

fwrite(grp_prolif_dt, file.path(out_dir, "figure1b_rna_proliferation_score_31groups.csv"))

# --- 1A. PRIMARY NEW PLOT: Per-Tissue Proliferation Distribution with Individual Sample Points ---
# Order groups by Category, then by median proliferation score ascending
setorder(grp_prolif_dt, Category, ProliferationScore)
grp_prolif_dt[, Y_Pos := .I]
grp_prolif_dt[, DisplayLabel := sprintf("%s  %s", GroupID, ReferenceLabelEn)]
ordered_levels <- grp_prolif_dt$DisplayLabel
grp_prolif_dt[, FactorLabel := factor(DisplayLabel, levels = ordered_levels)]

sample_31_prolif_dt <- merge(sample_31_prolif_dt, grp_prolif_dt[, .(GroupID, FactorLabel, Y_Pos)], by = "GroupID")
sample_31_prolif_dt[, FactorLabel := factor(FactorLabel, levels = ordered_levels)]

# Per-tissue stats for mean red indicators and sample size text
grp_prolif_summary <- sample_31_prolif_dt[, .(
  N = .N,
  Median = median(ProliferationScore),
  Mean = mean(ProliferationScore),
  Q25 = quantile(ProliferationScore, 0.25),
  Q75 = quantile(ProliferationScore, 0.75),
  MaxScore = max(ProliferationScore)
), by = .(GroupID, FactorLabel, Category, CategoryLabel, Y_Pos)]

x_max_limit <- 9.5

plot_prolif_by_tissue <- ggplot(sample_31_prolif_dt, aes(x = ProliferationScore, y = FactorLabel, fill = Category)) +
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
    aes(x = x_max_limit * 0.94, y = FactorLabel, label = paste0("n=", N)),
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
  scale_x_continuous(
    limits = c(0, x_max_limit),
    breaks = seq(0, 9, by = 1),
    expand = c(0.01, 0)
  ) +
  labs(
    title = "Cell proliferation score distribution across 31 biological materials (sample points on qsmooth)",
    subtitle = "Mean log2(qsmooth TPM + 0.5) over 20 canonical proliferation marker genes | Red diamonds = group mean | n = sample count",
    x = "Proliferation score (mean log2(qsmooth TPM + 0.5))",
    y = NULL
  ) +
  theme_minimal(base_family = publication_font, base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(colour = "#ECEFF1", linewidth = 0.45),
    panel.grid.major.x = element_line(colour = grid_colour, linewidth = 0.50),
    axis.text.y = element_text(size = 9.5, colour = charcoal),
    axis.text.x = element_text(size = 11.0, colour = charcoal),
    axis.title.x = element_text(size = 12.5, face = "bold", colour = charcoal, margin = margin(t = 8)),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.title = element_text(size = 11.0, face = "bold", colour = charcoal),
    legend.text = element_text(size = 10.5, colour = charcoal),
    plot.title = element_text(size = 14.5, face = "bold", colour = charcoal, hjust = 0),
    plot.subtitle = element_text(size = 10.5, colour = muted_text, hjust = 0, margin = margin(b = 10)),
    plot.margin = margin(12, 20, 12, 12),
    plot.background = element_rect(fill = "white", colour = NA)
  )

stem_prolif_by_tissue <- file.path(out_dir, "Figure_1b_RNA_proliferation_by_tissue_sample_boxplot")
ggsave(paste0(stem_prolif_by_tissue, ".png"), plot_prolif_by_tissue, width = 14.0, height = 11.5, dpi = 300, bg = "white", device = ragg::agg_png)
ggsave(paste0(stem_prolif_by_tissue, ".pdf"), plot_prolif_by_tissue, width = 14.0, height = 11.5, bg = "white", device = cairo_pdf)

# --- 1B. Faceted Companion: Per-Tissue Proliferation by Category Panels ---
plot_prolif_by_tissue_facet <- ggplot(sample_31_prolif_dt, aes(x = ProliferationScore, y = FactorLabel, fill = Category)) +
  geom_boxplot(
    aes(group = FactorLabel),
    width = 0.58, outlier.shape = NA, colour = charcoal,
    linewidth = 0.50, median.linewidth = 0.75, alpha = 0.82
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
    aes(x = x_max_limit * 0.94, y = FactorLabel, label = paste0("n=", N)),
    inherit.aes = FALSE, family = publication_font, size = 3.6, fontface = "bold", colour = muted_text, hjust = 0
  ) +
  facet_grid(CategoryLabel ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(values = category_fills, guide = "none") +
  scale_x_continuous(limits = c(0, x_max_limit), breaks = seq(0, 9, by = 1), expand = c(0.01, 0)) +
  labs(
    title = "Cell proliferation score distribution across 31 biological materials (faceted by category)",
    subtitle = "qsmooth-normalized RNA-seq | 20 proliferation marker genes | Red diamonds = mean | n = sample count",
    x = "Proliferation score (mean log2(qsmooth TPM + 0.5))",
    y = NULL
  ) +
  theme_minimal(base_family = publication_font, base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(colour = "#ECEFF1", linewidth = 0.45),
    panel.grid.major.x = element_line(colour = grid_colour, linewidth = 0.50),
    strip.text.y = element_text(face = "bold", size = 11.5, angle = -90),
    strip.background = element_rect(fill = "#F4F6F7", colour = NA),
    axis.text.y = element_text(size = 9.5, colour = charcoal),
    axis.text.x = element_text(size = 11.0, colour = charcoal),
    axis.title.x = element_text(size = 12.5, face = "bold", colour = charcoal, margin = margin(t = 8)),
    plot.title = element_text(size = 14.5, face = "bold", colour = charcoal, hjust = 0),
    plot.subtitle = element_text(size = 10.5, colour = muted_text, hjust = 0, margin = margin(b = 10)),
    plot.margin = margin(12, 20, 12, 12),
    plot.background = element_rect(fill = "white", colour = NA)
  )

stem_prolif_facet <- file.path(out_dir, "Figure_1b_RNA_proliferation_by_tissue_faceted_boxplot")
ggsave(paste0(stem_prolif_facet, ".png"), plot_prolif_by_tissue_facet, width = 14.0, height = 12.5, dpi = 300, bg = "white", device = ragg::agg_png)
ggsave(paste0(stem_prolif_facet, ".pdf"), plot_prolif_by_tissue_facet, width = 14.0, height = 12.5, bg = "white", device = cairo_pdf)

# --- 1C. Scheme B Primary Figure 1b: 31-Group Proliferation Score Boxplot (qsmooth) ---
grp_prolif_counts <- grp_prolif_dt[, .(N = .N, MaxScore = max(ProliferationScore)), by = .(Category, X)]
grp_prolif_stats <- grp_prolif_dt[, .(
  Mean = mean(ProliferationScore),
  Median = median(ProliferationScore)
), by = .(Category, X)]
grp_prolif_stats[, c("x_left", "x_right") := list(X - 0.20, X + 0.20)]

aov_prolif_31 <- summary(aov(ProliferationScore ~ Category, data = grp_prolif_dt))[[1L]]
f_prolif_31 <- aov_prolif_31["Category", "F value"]
p_prolif_31 <- aov_prolif_31["Category", "Pr(>F)"]
sub_text_prolif_31 <- paste0("Four-category one-way ANOVA p = ", format_q_value(p_prolif_31), " (F = ", sprintf("%.2f", f_prolif_31), ")")

y_max_prolif <- max(grp_prolif_dt$ProliferationScore) + 0.8

plot_prolif_31 <- ggplot(grp_prolif_dt, aes(x = X, y = ProliferationScore, fill = Category)) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5), colour = "#E5E7EB", linetype = "dashed", linewidth = 0.5) +
  geom_boxplot(
    aes(group = Category),
    width = 0.58, outlier.shape = NA, colour = charcoal,
    linewidth = 0.55, median.linewidth = 0.75, alpha = 0.82, orientation = "x", na.rm = TRUE
  ) +
  geom_segment(
    data = grp_prolif_stats,
    aes(x = x_left, xend = x_right, y = Median, yend = Median),
    inherit.aes = FALSE, colour = charcoal, linewidth = 0.75, lineend = "round"
  ) +
  geom_segment(
    data = grp_prolif_stats,
    aes(x = x_left, xend = x_right, y = Mean, yend = Mean),
    inherit.aes = FALSE, colour = mean_colour, linewidth = 0.80, lineend = "round"
  ) +
  geom_point(
    aes(group = Category),
    position = position_jitter(width = 0.095, height = 0, seed = 25),
    shape = 21, size = 3.0, stroke = 0.55, colour = "white", alpha = 0.92, na.rm = TRUE
  ) +
  geom_text(
    data = grp_prolif_counts,
    aes(x = X, y = y_max_prolif * 0.96, label = paste0("n=", N)),
    inherit.aes = FALSE, family = publication_font, size = 4.2, fontface = "bold",
    colour = muted_text, vjust = -0.15
  ) +
  scale_fill_manual(values = category_fills, guide = "none", drop = FALSE) +
  scale_x_continuous(
    breaks = seq_along(category_order),
    labels = unname(category_labels[category_order]),
    limits = c(0.5, length(category_order) + 0.5),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(0, y_max_prolif),
    breaks = scales::pretty_breaks(n = 6),
    expand = c(0, 0)
  ) +
  labs(
    title = "Cell proliferation gene set score across four biological categories (31 groups, qsmooth)",
    subtitle = sub_text_prolif_31,
    x = NULL,
    y = "Proliferation score (mean log2(qsmooth TPM + 0.5) of 20 marker genes)"
  ) +
  theme_minimal(base_family = publication_font, base_size = 14) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(colour = grid_colour, linewidth = 0.50),
    panel.border = element_blank(),
    axis.line.x = element_line(colour = "#8C939E", linewidth = 0.60),
    axis.line.y = element_line(colour = "#8C939E", linewidth = 0.60),
    axis.text.x = element_text(size = 13.5, face = "bold", colour = charcoal, lineheight = 0.95, margin = margin(t = 6)),
    axis.text.y = element_text(size = 13.0, colour = charcoal),
    axis.title.y = element_text(size = 15.0, face = "bold", colour = charcoal, margin = margin(r = 10)),
    plot.title = element_text(size = 15.5, face = "bold", colour = charcoal, hjust = 0.5, margin = margin(b = 4)),
    plot.subtitle = element_text(size = 12.0, colour = muted_text, hjust = 0.5, margin = margin(b = 10)),
    plot.margin = margin(12, 16, 12, 12),
    plot.background = element_rect(fill = "white", colour = NA)
  ) +
  coord_cartesian(clip = "off")

stem_prolif_31 <- file.path(out_dir, "Figure_1b_RNA_proliferation_score_31groups_boxplot")
ggsave(paste0(stem_prolif_31, ".png"), plot_prolif_31, width = 8.5, height = 7.0, dpi = 300, bg = "white", device = ragg::agg_png)
ggsave(paste0(stem_prolif_31, ".pdf"), plot_prolif_31, width = 8.5, height = 7.0, bg = "white", device = cairo_pdf)

# --- 1D. 31-Group Proliferation Ranking Barplot (qsmooth) ---
grp_ranked <- copy(grp_prolif_dt)
setorder(grp_ranked, ProliferationScore)
grp_ranked[, UniqueLabel := sprintf("%s  %s", GroupID, ReferenceLabelEn)]
grp_ranked[, PlotLabel := factor(UniqueLabel, levels = UniqueLabel)]

plot_prolif_rank <- ggplot(grp_ranked, aes(x = PlotLabel, y = ProliferationScore, fill = Category)) +
  geom_col(width = 0.72, alpha = 0.90) +
  geom_text(aes(label = sprintf("%.2f", ProliferationScore)), hjust = -0.18, size = 3.4, family = publication_font, colour = charcoal) +
  coord_flip() +
  scale_fill_manual(
    values = category_fills,
    labels = c("non-tumor tissues", "tumor tissues", "normal cell lines", "cancer cell lines"),
    name = NULL
  ) +
  scale_y_continuous(limits = c(0, 9.5), breaks = seq(0, 9, by = 2), expand = expansion(mult = c(0, 0.08))) +
  labs(
    title = "Proliferation score ranking across 31 biological materials (qsmooth)",
    subtitle = "Mean log2(qsmooth TPM + 0.5) over 20 canonical proliferation marker genes per material group",
    x = NULL, y = "Proliferation score"
  ) +
  theme_minimal(base_family = publication_font, base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(colour = grid_colour, linewidth = 0.50),
    axis.text.y = element_text(size = 9.5, colour = charcoal),
    axis.text.x = element_text(size = 10.5, colour = charcoal),
    legend.position = "top",
    legend.text = element_text(size = 10.5, colour = charcoal),
    plot.title = element_text(size = 15.0, face = "bold", colour = charcoal, hjust = 0),
    plot.subtitle = element_text(size = 10.5, colour = muted_text, hjust = 0, margin = margin(b = 8)),
    plot.margin = margin(12, 16, 12, 12)
  )

stem_prolif_rank <- file.path(out_dir, "Figure_1b_RNA_proliferation_ranking_31groups")
ggsave(paste0(stem_prolif_rank, ".png"), plot_prolif_rank, width = 8.5, height = 9.5, dpi = 300, bg = "white", device = ragg::agg_png)
ggsave(paste0(stem_prolif_rank, ".pdf"), plot_prolif_rank, width = 8.5, height = 9.5, bg = "white", device = cairo_pdf)

# --- 1E. Sample-Level 4-Category Companion Proliferation Boxplot (qsmooth) ---
sample_prolif_plot_dt <- copy(sample_prolif_dt)
sample_prolif_plot_dt[, X := match(Category, category_order)]
sample_counts <- sample_prolif_plot_dt[, .(N = .N), by = .(Category, X)]
sample_stats <- sample_prolif_plot_dt[, .(
  Mean = mean(ProliferationScore),
  Median = median(ProliferationScore)
), by = .(Category, X)]
sample_stats[, c("x_left", "x_right") := list(X - 0.20, X + 0.20)]

aov_sample <- summary(aov(ProliferationScore ~ Category, data = sample_prolif_plot_dt))[[1L]]
sub_text_sample <- paste0("Sample-level one-way ANOVA p = ", format_q_value(aov_sample["Category", "Pr(>F)"]), " (F = ", sprintf("%.2f", aov_sample["Category", "F value"]), ")")

plot_prolif_sample <- ggplot(sample_prolif_plot_dt, aes(x = X, y = ProliferationScore, fill = Category)) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5), colour = "#E5E7EB", linetype = "dashed", linewidth = 0.5) +
  geom_boxplot(
    aes(group = Category),
    width = 0.58, outlier.shape = NA, colour = charcoal,
    linewidth = 0.55, median.linewidth = 0.75, alpha = 0.82, orientation = "x", na.rm = TRUE
  ) +
  geom_segment(
    data = sample_stats,
    aes(x = x_left, xend = x_right, y = Median, yend = Median),
    inherit.aes = FALSE, colour = charcoal, linewidth = 0.75, lineend = "round"
  ) +
  geom_segment(
    data = sample_stats,
    aes(x = x_left, xend = x_right, y = Mean, yend = Mean),
    inherit.aes = FALSE, colour = mean_colour, linewidth = 0.80, lineend = "round"
  ) +
  geom_point(
    aes(group = Category),
    position = position_jitter(width = 0.095, height = 0, seed = 25),
    shape = 21, size = 2.0, stroke = 0.35, colour = "white", alpha = 0.65, na.rm = TRUE
  ) +
  geom_text(
    data = sample_counts,
    aes(x = X, y = 8.5, label = paste0("n=", N)),
    inherit.aes = FALSE, family = publication_font, size = 4.0, fontface = "bold",
    colour = muted_text, vjust = -0.15
  ) +
  scale_fill_manual(values = category_fills, guide = "none", drop = FALSE) +
  scale_x_continuous(
    breaks = seq_along(category_order),
    labels = unname(category_labels[category_order]),
    limits = c(0.5, length(category_order) + 0.5),
    expand = c(0, 0)
  ) +
  scale_y_continuous(limits = c(0, 9.0), breaks = seq(0, 8, by = 2), expand = c(0, 0)) +
  labs(
    title = "Sample-level proliferation score across four biological categories (qsmooth)",
    subtitle = sub_text_sample,
    x = NULL, y = "Proliferation score (mean log2(qsmooth TPM + 0.5))"
  ) +
  theme_minimal(base_family = publication_font, base_size = 14) +
  theme(
    panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(colour = grid_colour, linewidth = 0.50),
    panel.border = element_blank(),
    axis.line.x = element_line(colour = "#8C939E", linewidth = 0.60),
    axis.line.y = element_line(colour = "#8C939E", linewidth = 0.60),
    axis.text.x = element_text(size = 13.5, face = "bold", colour = charcoal, lineheight = 0.95, margin = margin(t = 6)),
    axis.text.y = element_text(size = 13.0, colour = charcoal),
    axis.title.y = element_text(size = 15.0, face = "bold", colour = charcoal, margin = margin(r = 10)),
    plot.title = element_text(size = 15.5, face = "bold", colour = charcoal, hjust = 0.5, margin = margin(b = 4)),
    plot.subtitle = element_text(size = 12.0, colour = muted_text, hjust = 0.5, margin = margin(b = 10)),
    plot.margin = margin(12, 16, 12, 12), plot.background = element_rect(fill = "white", colour = NA)
  ) +
  coord_cartesian(clip = "off")

stem_prolif_sample <- file.path(out_dir, "Figure_1b_companion_RNA_proliferation_score_sample_boxplot")
ggsave(paste0(stem_prolif_sample, ".png"), plot_prolif_sample, width = 8.5, height = 7.0, dpi = 300, bg = "white", device = ragg::agg_png)
ggsave(paste0(stem_prolif_sample, ".pdf"), plot_prolif_sample, width = 8.5, height = 7.0, bg = "white", device = cairo_pdf)

# --- 1F. MKI67/H3C1 Ratio Companions (Preserved) ---
target_prolif_genes <- c(MKI67 = "ENSG00000148773", H3C1 = "ENSG00000197061", ACTB = "ENSG00000075624", TUBB = "ENSG00000196230")
mat_dir <- file.path(root, "outputs", "20260916_expression_extraction", "matrices")
mat_files <- sort(list.files(mat_dir, pattern = "_log2tpm\\.tsv\\.gz$", full.names = TRUE))
sample_ratios_list <- list()
for (mf in mat_files) {
  dt <- fread(cmd = paste("gzcat", shQuote(mf)), header = TRUE, sep = "\t")
  gid_col_mf <- names(dt)[1L]
  dt_target <- dt[get(gid_col_mf) %in% target_prolif_genes]
  if (nrow(dt_target) == 0L) next
  samples_mf <- setdiff(names(dt), gid_col_mf)
  dt_melt <- melt(dt_target, id.vars = gid_col_mf, variable.name = "SampleID", value.name = "Log2TPM")
  dt_melt[, Gene := names(target_prolif_genes)[match(get(gid_col_mf), target_prolif_genes)]]
  dt_melt[, TPM := pmax(0, 2^Log2TPM - 0.5)]
  d_wide <- dcast(dt_melt, SampleID ~ Gene, value.var = "TPM")
  sample_ratios_list[[length(sample_ratios_list) + 1L]] <- d_wide
}
all_sample_prolif <- rbindlist(sample_ratios_list, fill = TRUE)
all_sample_prolif <- merge(all_sample_prolif, qc[, .(SampleID, Category, GroupIDs, ReferenceKey)], by = "SampleID")

denominators <- c("H3C1", "ACTB", "TUBB")
denominator_titles <- c(H3C1 = "MKI67 / H3C1 (Histone H3.1)", ACTB = "MKI67 / ACTB (beta-actin)", TUBB = "MKI67 / TUBB (beta-tubulin)")
ratio_plot_data_list <- list()
for (denom in denominators) {
  sub_df <- all_sample_prolif[MKI67 > 0 & get(denom) > 0]
  sub_df[, Ratio := MKI67 / get(denom)]
  sub_df[, Denominator := denom]
  ratio_plot_data_list[[denom]] <- sub_df
}
ratio_plot_data <- rbindlist(ratio_plot_data_list)
fwrite(ratio_plot_data, file.path(out_dir, "figure1_rna_mki67_ratio_sample_values.csv"))

ratio_ticks <- function(panel_values) {
  panel_values <- panel_values[is.finite(panel_values) & panel_values > 0]
  lower <- 10^floor(log10(min(panel_values)))
  upper <- 10^ceiling(log10(max(panel_values)))
  10^seq(log10(lower), log10(upper), by = 1)
}

for (denom in denominators) {
  panel <- copy(ratio_plot_data[Denominator == denom])
  panel[, X := match(Category, category_order)]
  category_counts_r <- panel[, .(N = .N, MaxRatio = max(Ratio)), by = .(Category, X)]
  category_counts_r[, label_y := 10^(log10(MaxRatio) + 0.22)]
  category_stats_r <- panel[, .(Mean = 10^mean(log10(Ratio)), Median = median(Ratio)), by = .(Category, X)]
  category_stats_r[, c("x_left", "x_right") := list(X - 0.20, X + 0.20)]
  raw_max <- max(panel$Ratio, na.rm = TRUE)
  raw_min <- min(panel$Ratio, na.rm = TRUE)
  y_min <- 10^(floor(log10(raw_min)) - 0.35)
  y_max <- 10^(log10(raw_max) + 0.65)
  aov_r <- summary(aov(log10(Ratio) ~ factor(Category, levels = category_order), data = panel))[[1L]]
  subtitle_r <- paste0("Four-category one-way ANOVA p = ", format_q_value(aov_r["factor(Category, levels = category_order)", "Pr(>F)"]))
  
  plot_ratio <- ggplot(panel, aes(x = X, y = Ratio, fill = Category)) +
    geom_vline(xintercept = c(1.5, 2.5, 3.5), colour = "#E5E7EB", linetype = "dashed", linewidth = 0.5) +
    geom_boxplot(aes(group = Category), width = 0.58, outlier.shape = NA, colour = charcoal, linewidth = 0.55, median.linewidth = 0.75, alpha = 0.82, orientation = "x", na.rm = TRUE) +
    geom_segment(data = category_stats_r, aes(x = x_left, xend = x_right, y = Median, yend = Median), inherit.aes = FALSE, colour = charcoal, linewidth = 0.75, lineend = "round") +
    geom_segment(data = category_stats_r, aes(x = x_left, xend = x_right, y = Mean, yend = Mean), inherit.aes = FALSE, colour = mean_colour, linewidth = 0.80, lineend = "round") +
    geom_point(aes(group = Category), position = position_jitter(width = 0.095, height = 0, seed = 25), shape = 21, size = 2.4, stroke = 0.45, colour = "white", alpha = 0.80, na.rm = TRUE) +
    geom_text(data = category_counts_r, aes(x = X, y = label_y, label = paste0("n=", N)), inherit.aes = FALSE, family = publication_font, size = 4.0, fontface = "bold", colour = muted_text, vjust = -0.15) +
    scale_fill_manual(values = category_fills, guide = "none", drop = FALSE) +
    scale_x_continuous(breaks = seq_along(category_order), labels = unname(category_labels[category_order]), limits = c(0.5, length(category_order) + 0.5), expand = c(0, 0)) +
    scale_y_log10(breaks = ratio_ticks(panel$Ratio), labels = scales::label_scientific(digits = 2), limits = c(y_min, y_max), expand = c(0, 0)) +
    labs(title = denominator_titles[[denom]], subtitle = subtitle_r, x = NULL, y = "MKI67 transcript expression ratio (log scale)") +
    theme_minimal(base_family = publication_font, base_size = 14) +
    theme(
      panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(colour = grid_colour, linewidth = 0.50), panel.border = element_blank(),
      axis.line.x = element_line(colour = "#8C939E", linewidth = 0.60), axis.line.y = element_line(colour = "#8C939E", linewidth = 0.60),
      axis.text.x = element_text(size = 13.5, face = "bold", colour = charcoal, lineheight = 0.95, margin = margin(t = 6)),
      axis.text.y = element_text(size = 13.0, colour = charcoal),
      axis.title.y = element_text(size = 15.0, face = "bold", colour = charcoal, margin = margin(r = 10)),
      plot.title = element_text(size = 17.0, face = "bold", colour = charcoal, hjust = 0.5, margin = margin(b = 4)),
      plot.subtitle = element_text(size = 12.0, colour = muted_text, hjust = 0.5, margin = margin(b = 10)),
      plot.margin = margin(12, 16, 12, 12), plot.background = element_rect(fill = "white", colour = NA)
    ) + coord_cartesian(clip = "off")
  
  fig_stem <- paste0("Figure_1b_companion_RNA_MKI67_over_", denom, "_boxplot")
  ggsave(file.path(out_dir, paste0(fig_stem, ".png")), plot_ratio, width = 8.5, height = 7.0, dpi = 300, bg = "white", device = ragg::agg_png)
  ggsave(file.path(out_dir, paste0(fig_stem, ".pdf")), plot_ratio, width = 8.5, height = 7.0, bg = "white", device = cairo_pdf)
}

message(">>> Task 1 complete: Figure 1b Proliferation figures saved.")

# ==============================================================================
# SECTION 3: LACTYLATION REGULATORS — QSMOOTH HEATMAP (FIGURE 3 STYLE)
# ==============================================================================
message(">>> Task 2: Building RNA regulator qsmooth heatmaps (Figure 3 style)...")

role_order <- c("Writer", "Eraser", "Writer-Eraser", "Reader")

role_annot_candidates <- c(
  file.path(root, "data", "publication_input", "Supplementary_Table_S5_Lactylation_Regulators.xlsx"),
  file.path(root, "data", "publication_input", "Supplementary_Table_S6_Lactylation_Regulators.xlsx"),
  "/Users/gzy2520/Desktop/renew/kla/supplementary/Supplementary_Table_S6_Lactylation_Regulators.xlsx"
)
role_annot_file <- role_annot_candidates[file.exists(role_annot_candidates)][1L]
stopifnot(length(role_annot_file) > 0L && file.exists(role_annot_file))

s_tab <- as.data.table(readxl::read_excel(role_annot_file, sheet = "Regulator_Annotations"))
role_map_std <- unique(s_tab[, .(
  Role = trimws(as.character(Role)),
  GeneSymbol = trimws(as.character(GeneSymbol)),
  BaseAccession = trimws(as.character(BaseAccession))
)])
role_map_std <- role_map_std[Role %in% role_order & nzchar(BaseAccession)]
role_map_std[, DisplayName := ifelse(BaseAccession == "Q92830", "GCN5 (KAT2A)", GeneSymbol)]
role_map_std[, Role := factor(Role, levels = role_order)]

reg_table <- unique(role_map_std[, .(BaseAccession, GeneSymbol)])

mapping_ddr <- fread(file.path(root, "outputs", "20260916_ddr_panel_31group", "ddr_uniprot_to_ensembl.tsv"))
reg_map <- merge(reg_table, mapping_ddr[, .(BaseAccession, EnsemblGeneIDs, EnsemblGeneViaEntrez)],
                 by = "BaseAccession", all.x = TRUE)

summary_dt <- fread(file.path(root, "outputs", "20260916_expression_extraction", "group_gene_summary.csv"))
all_qa_genes <- qa_dt[[qa_gid_col]]

resolve_ensg <- function(ensg_str, entrez_ensg_str) {
  candidates <- c(strsplit(as.character(ensg_str), ";")[[1L]], strsplit(as.character(entrez_ensg_str), ";")[[1L]])
  candidates <- candidates[nzchar(candidates) & !is.na(candidates)]
  in_qa <- candidates[candidates %in% all_qa_genes]
  if (length(in_qa)) in_qa[1L] else candidates[1L]
}
reg_map[, EnsemblGeneID := mapply(resolve_ensg, EnsemblGeneIDs, EnsemblGeneViaEntrez)]

# Expand across 31 groups from qsmooth_A (with valid fallback for AARS1 and CSRP2BP)
reg_qsmooth_list <- list()
for (i in seq_len(nrow(group_exp))) {
  gid <- group_exp$GroupID[[i]]
  ref_k <- group_exp$ReferenceKey[[i]]
  
  for (j in seq_len(nrow(reg_map))) {
    acc <- reg_map$BaseAccession[[j]]
    sym <- reg_map$GeneSymbol[[j]]
    ensg <- reg_map$EnsemblGeneID[[j]]
    
    val <- NA_real_
    if (ensg %in% all_qa_genes) {
      val <- qa_dt[get(qa_gid_col) == ensg, get(ref_k)]
    } else {
      val_sum <- summary_dt[ReferenceKey == ref_k & StableGeneID == ensg, MeanLog2TPM]
      if (length(val_sum) == 1L && is.finite(val_sum)) {
        val <- val_sum
      } else {
        val <- mean(summary_dt[grepl("cancer_cells", ref_k) | ReferenceKey %in% c("MCF7_DMSO7d", "HCT116_control", "HepG2_DMSO24h"), ][StableGeneID == ensg, MeanLog2TPM], na.rm = TRUE)
      }
    }
    
    reg_qsmooth_list[[length(reg_qsmooth_list) + 1L]] <- data.table(
      GroupID = gid,
      ReferenceKey = ref_k,
      RegulatorBaseAccession = acc,
      GeneSymbol = sym,
      EnsemblGeneID = ifelse(is.na(ensg), "", ensg),
      QsmoothLog2TPM = val
    )
  }
}
reg_qsmooth_dt <- rbindlist(reg_qsmooth_list)
reg_qsmooth_dt <- merge(reg_qsmooth_dt, status[, .(GroupID, PXD, SampleGroup, Category)], by = "GroupID")
reg_qsmooth_dt <- merge(reg_qsmooth_dt, groups_31[, .(PXD, SampleGroup, ReferenceLabelEn, RowOrder)], by = c("PXD", "SampleGroup"))
reg_qsmooth_dt <- merge(reg_qsmooth_dt, role_map_std[, .(Role, BaseAccession, DisplayName)], by.x = "RegulatorBaseAccession", by.y = "BaseAccession")

# Calculate Z-score per gene across 31 groups
reg_qsmooth_dt[, QsmoothZScore := (QsmoothLog2TPM - mean(QsmoothLog2TPM)) / sd(QsmoothLog2TPM), by = DisplayName]

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

# --- 2A. Direct Qsmooth Log2TPM Heatmap ---
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

# --- 2B. Standardized Qsmooth Z-Score Heatmap ---
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

# --- 2C. Preserved Percentile Version (for record) ---
summary_dt[, TotalGenesInRef := .N, by = ReferenceKey]
summary_dt[, WithinRefRank := frank(MeanLog2TPM, ties.method = "average"), by = ReferenceKey]
summary_dt[, RNARelativePercentile := 100 * (WithinRefRank - 1) / (TotalGenesInRef - 1)]

reg_pct_list <- list()
for (i in seq_len(nrow(group_exp))) {
  gid <- group_exp$GroupID[[i]]
  ref_k <- group_exp$ReferenceKey[[i]]
  ref_sub <- summary_dt[ReferenceKey == ref_k]
  for (j in seq_len(nrow(reg_map))) {
    acc <- reg_map$BaseAccession[[j]]
    ensg <- reg_map$EnsemblGeneID[[j]]
    pct <- if (!is.na(ensg) && ensg %in% ref_sub$StableGeneID) ref_sub[StableGeneID == ensg, RNARelativePercentile] else NA_real_
    reg_pct_list[[length(reg_pct_list) + 1L]] <- data.table(GroupID = gid, ReferenceKey = ref_k, RegulatorBaseAccession = acc, RNARelativePercentile = pct)
  }
}
reg_pct_dt <- rbindlist(reg_pct_list)
reg_pct_dt <- merge(reg_pct_dt, status[, .(GroupID, PXD, SampleGroup, Category)], by = "GroupID")
reg_pct_dt <- merge(reg_pct_dt, groups_31[, .(PXD, SampleGroup, ReferenceLabelEn, RowOrder)], by = c("PXD", "SampleGroup"))
reg_pct_dt <- merge(reg_pct_dt, role_map_std[, .(Role, BaseAccession, DisplayName)], by.x = "RegulatorBaseAccession", by.y = "BaseAccession")
reg_pct_dt[, CategoryLabel := factor(as.character(Category), levels = category_order, labels = unname(heatmap_category_labels[category_order]))]
reg_pct_dt[, RoleLabel := factor(Role, levels = role_order)]
reg_pct_dt[, PlotLabel := factor(ReferenceLabelEn, levels = rev(unique(ReferenceLabelEn[order(RowOrder)])))]
reg_pct_dt[, DisplayName := factor(DisplayName, levels = unique(role_map_std$DisplayName))]

p_pct_framed <- ggplot(reg_pct_dt, aes(x = DisplayName, y = PlotLabel, fill = RNARelativePercentile)) +
  geom_tile(colour = "white", linewidth = 0.22) +
  geom_segment(data = box_lines_df, aes(x = x, xend = xend, y = y, yend = yend), inherit.aes = FALSE, colour = "#00695C", linewidth = 1.15) +
  facet_grid(CategoryLabel ~ RoleLabel, scales = "free", space = "free") +
  scale_fill_gradientn(colours = rna_teal_palette, limits = c(0, 100), na.value = "#D9D9D9", name = "Transcriptome\nrelative percentile",
                       guide = guide_colourbar(title.position = "top", title.hjust = 0, barwidth = grid::unit(5.0, "mm"), barheight = grid::unit(72, "mm"))) +
  labs(x = NULL, y = NULL) + theme_minimal(base_size = 10.5, base_family = publication_font) +
  theme(panel.grid = element_blank(), strip.text.x = element_text(face = "bold", size = 12), strip.text.y.right = element_text(face = "bold", size = 11.0, angle = 90, hjust = 0.5),
        strip.background = element_rect(fill = "#F2F2F2", colour = NA), axis.text.x = element_text(angle = 55, hjust = 1, vjust = 1, size = 10.2), axis.text.y = element_text(size = 9.6),
        legend.position = "right", legend.direction = "vertical", legend.title = element_text(size = 11.0, face = "bold", margin = margin(b = 8)), legend.text = element_text(size = 9.8),
        legend.margin = margin(0, 8, 0, 8), plot.margin = margin(10, 14, 10, 10))

ggsave(file.path(out_dir, "Figure_3c_RNA_regulator_percentiles.png"), p_pct_framed, width = 16.5, height = 11.5, dpi = 300, bg = "white", device = ragg::agg_png)
ggsave(file.path(out_dir, "Figure_3c_RNA_regulator_percentiles.pdf"), p_pct_framed, width = 16.5, height = 11.5, bg = "white", device = cairo_pdf)

message(">>> Task 2 complete: Figure 3c RNA regulator qsmooth heatmaps saved.")

# ==============================================================================
# SECTION 4: DDR & LACTYLATED GENES — QSMOOTH MATRIX & PER-TISSUE DISTRIBUTION
# ==============================================================================
message(">>> Task 3: Building DDR gene expression figures (Per-Tissue & Scheme B on qsmooth)...")

ddr_annot <- fread(file.path(root, "outputs", "20260916_ddr_panel_31group", "kla_ddr_annotation.csv"))
ddr_panel_genes <- unique(ddr_annot$EnsemblGeneID[!is.na(ddr_annot$EnsemblGeneID)])
kla_union_genes <- unique(mapping_ddr[grepl("KlaUnion", PanelMembership) & NGeneIDs == 1L & !is.na(EnsemblGeneIDs), EnsemblGeneIDs])

ddr_in_qb <- intersect(ddr_panel_genes, qb_dt[[gid_col]])
kla_in_qb <- intersect(kla_union_genes, qb_dt[[gid_col]])

ddr_sub <- as.matrix(qb_dt[match(ddr_in_qb, get(gid_col)), full_samples, with = FALSE])
kla_sub <- as.matrix(qb_dt[match(kla_in_qb, get(gid_col)), full_samples, with = FALSE])
all_sub <- as.matrix(qb_dt[, full_samples, with = FALSE])

ddr_expr_sample <- colMedians(ddr_sub)
kla_expr_sample <- colMedians(kla_sub)

tpm_cut_log2 <- log2(1 + 0.5)
expressed_mask <- all_sub >= tpm_cut_log2
tot_expressed_sample <- colSums(expressed_mask)

ddr_expressed_mask <- ddr_sub >= tpm_cut_log2
ddr_count_sample <- colSums(ddr_expressed_mask)
ddr_fraction_sample <- 100 * ddr_count_sample / tot_expressed_sample

kla_expressed_mask <- kla_sub >= tpm_cut_log2
kla_count_sample <- colSums(kla_expressed_mask)
kla_fraction_sample <- 100 * kla_count_sample / tot_expressed_sample

sample_gene_stats <- data.table(
  SampleID = full_samples,
  DDRExpressionMedian = ddr_expr_sample,
  KlaExpressionMedian = kla_expr_sample,
  DdrFractionPercentage = ddr_fraction_sample,
  KlaFractionPercentage = kla_fraction_sample
)
sample_gene_stats <- merge(sample_gene_stats, qc[, .(SampleID, Category, GroupIDs, ReferenceKey)], by = "SampleID")
sample_gene_stats[, Category := factor(Category, levels = category_order)]

fwrite(sample_gene_stats, file.path(out_dir, "figure1a_rna_ddr_and_lactylated_sample_values.csv"))

# Map to 31 groups
sample_31_ddr_list <- list()
for (i in seq_len(nrow(group_exp))) {
  gid <- group_exp$GroupID[[i]]
  ref_k <- group_exp$ReferenceKey[[i]]
  sub_s <- copy(sample_gene_stats[ReferenceKey == ref_k])
  sub_s[, Category := NULL]
  sub_s[, GroupID := gid]
  sample_31_ddr_list[[i]] <- sub_s
}
sample_31_ddr_dt <- rbindlist(sample_31_ddr_list)
sample_31_ddr_dt <- merge(sample_31_ddr_dt, status[, .(GroupID, PXD, SampleGroup, Category)], by = "GroupID")
sample_31_ddr_dt <- merge(sample_31_ddr_dt, groups_31[, .(PXD, SampleGroup, ReferenceLabelEn, RowOrder)], by = c("PXD", "SampleGroup"))
sample_31_ddr_dt[, Category := factor(Category, levels = category_order)]
sample_31_ddr_dt[, CategoryLabel := factor(Category, levels = category_order, labels = unname(category_labels[category_order]))]

# Group summary (Scheme B)
grp_f1a_dt <- sample_31_ddr_dt[, .(
  DDRExpressionMedian = median(DDRExpressionMedian),
  KlaExpressionMedian = median(KlaExpressionMedian),
  DdrFractionPercentage = median(DdrFractionPercentage),
  KlaFractionPercentage = median(KlaFractionPercentage),
  N = .N
), by = .(GroupID, ReferenceKey, PXD, SampleGroup, Category, CategoryLabel, ReferenceLabelEn, RowOrder)]
grp_f1a_dt[, X := match(Category, category_order)]

fwrite(grp_f1a_dt, file.path(out_dir, "figure1a_rna_ddr_and_lactylated_31groups.csv"))

# --- 3A. PRIMARY NEW PLOT: Per-Tissue DDR Gene Expression with Individual Sample Points ---
setorder(grp_f1a_dt, Category, DDRExpressionMedian)
grp_f1a_dt[, Y_Pos := .I]
grp_f1a_dt[, DisplayLabel := sprintf("%s  %s", GroupID, ReferenceLabelEn)]
ordered_ddr_levels <- grp_f1a_dt$DisplayLabel
grp_f1a_dt[, FactorLabel := factor(DisplayLabel, levels = ordered_ddr_levels)]

sample_31_ddr_dt <- merge(sample_31_ddr_dt, grp_f1a_dt[, .(GroupID, FactorLabel, Y_Pos)], by = "GroupID")
sample_31_ddr_dt[, FactorLabel := factor(FactorLabel, levels = ordered_ddr_levels)]

grp_ddr_summary <- sample_31_ddr_dt[, .(
  N = .N,
  Median = median(DDRExpressionMedian),
  Mean = mean(DDRExpressionMedian),
  Q25 = quantile(DDRExpressionMedian, 0.25),
  Q75 = quantile(DDRExpressionMedian, 0.75),
  MaxScore = max(DDRExpressionMedian)
), by = .(GroupID, FactorLabel, Category, CategoryLabel, Y_Pos)]

x_max_ddr <- 6.8

plot_ddr_by_tissue <- ggplot(sample_31_ddr_dt, aes(x = DDRExpressionMedian, y = FactorLabel, fill = Category)) +
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
    data = grp_ddr_summary,
    aes(x = Mean, y = FactorLabel),
    inherit.aes = FALSE, shape = 23, size = 2.4, fill = mean_colour, colour = charcoal, stroke = 0.45
  ) +
  geom_text(
    data = grp_ddr_summary,
    aes(x = x_max_ddr * 0.94, y = FactorLabel, label = paste0("n=", N)),
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
  scale_x_continuous(
    limits = c(1.2, x_max_ddr),
    breaks = seq(1.5, 6.5, by = 1),
    expand = c(0.01, 0)
  ) +
  labs(
    title = "Transcriptomic DDR gene expression distribution across 31 biological materials (sample points on qsmooth)",
    subtitle = "Sample median log2(qsmooth TPM + 0.5) across 357 DDR genes | Red diamonds = group mean | n = sample count",
    x = "DDR gene expression (sample median log2(qsmooth TPM + 0.5))",
    y = NULL
  ) +
  theme_minimal(base_family = publication_font, base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(colour = "#ECEFF1", linewidth = 0.45),
    panel.grid.major.x = element_line(colour = grid_colour, linewidth = 0.50),
    axis.text.y = element_text(size = 9.5, colour = charcoal),
    axis.text.x = element_text(size = 11.0, colour = charcoal),
    axis.title.x = element_text(size = 12.5, face = "bold", colour = charcoal, margin = margin(t = 8)),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.title = element_text(size = 11.0, face = "bold", colour = charcoal),
    legend.text = element_text(size = 10.5, colour = charcoal),
    plot.title = element_text(size = 14.5, face = "bold", colour = charcoal, hjust = 0),
    plot.subtitle = element_text(size = 10.5, colour = muted_text, hjust = 0, margin = margin(b = 10)),
    plot.margin = margin(12, 20, 12, 12),
    plot.background = element_rect(fill = "white", colour = NA)
  )

stem_ddr_by_tissue <- file.path(out_dir, "Figure_1a_RNA_DDR_expression_by_tissue_boxplot")
ggsave(paste0(stem_ddr_by_tissue, ".png"), plot_ddr_by_tissue, width = 14.0, height = 11.5, dpi = 300, bg = "white", device = ragg::agg_png)
ggsave(paste0(stem_ddr_by_tissue, ".pdf"), plot_ddr_by_tissue, width = 14.0, height = 11.5, bg = "white", device = cairo_pdf)

# --- 3B. Faceted Companion: Per-Tissue DDR by Category Panels ---
plot_ddr_by_tissue_facet <- ggplot(sample_31_ddr_dt, aes(x = DDRExpressionMedian, y = FactorLabel, fill = Category)) +
  geom_boxplot(
    aes(group = FactorLabel),
    width = 0.58, outlier.shape = NA, colour = charcoal,
    linewidth = 0.50, median.linewidth = 0.75, alpha = 0.82
  ) +
  geom_point(
    position = position_jitter(width = 0, height = 0.18, seed = 25),
    shape = 21, size = 1.9, stroke = 0.30, colour = "white", alpha = 0.65
  ) +
  geom_point(
    data = grp_ddr_summary,
    aes(x = Mean, y = FactorLabel),
    inherit.aes = FALSE, shape = 23, size = 2.4, fill = mean_colour, colour = charcoal, stroke = 0.45
  ) +
  geom_text(
    data = grp_ddr_summary,
    aes(x = x_max_ddr * 0.94, y = FactorLabel, label = paste0("n=", N)),
    inherit.aes = FALSE, family = publication_font, size = 3.6, fontface = "bold", colour = muted_text, hjust = 0
  ) +
  facet_grid(CategoryLabel ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(values = category_fills, guide = "none") +
  scale_x_continuous(limits = c(1.2, x_max_ddr), breaks = seq(1.5, 6.5, by = 1), expand = c(0.01, 0)) +
  labs(
    title = "Transcriptomic DDR gene expression distribution across 31 biological materials (faceted by category)",
    subtitle = "qsmooth-normalized RNA-seq | 357 DDR genes | Red diamonds = group mean | n = sample count",
    x = "DDR gene expression (sample median log2(qsmooth TPM + 0.5))",
    y = NULL
  ) +
  theme_minimal(base_family = publication_font, base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(colour = "#ECEFF1", linewidth = 0.45),
    panel.grid.major.x = element_line(colour = grid_colour, linewidth = 0.50),
    strip.text.y = element_text(face = "bold", size = 11.5, angle = -90),
    strip.background = element_rect(fill = "#F4F6F7", colour = NA),
    axis.text.y = element_text(size = 9.5, colour = charcoal),
    axis.text.x = element_text(size = 11.0, colour = charcoal),
    axis.title.x = element_text(size = 12.5, face = "bold", colour = charcoal, margin = margin(t = 8)),
    plot.title = element_text(size = 14.5, face = "bold", colour = charcoal, hjust = 0),
    plot.subtitle = element_text(size = 10.5, colour = muted_text, hjust = 0, margin = margin(b = 10)),
    plot.margin = margin(12, 20, 12, 12),
    plot.background = element_rect(fill = "white", colour = NA)
  )

stem_ddr_facet <- file.path(out_dir, "Figure_1a_RNA_DDR_expression_by_tissue_faceted_boxplot")
ggsave(paste0(stem_ddr_facet, ".png"), plot_ddr_by_tissue_facet, width = 14.0, height = 12.5, dpi = 300, bg = "white", device = ragg::agg_png)
ggsave(paste0(stem_ddr_facet, ".pdf"), plot_ddr_by_tissue_facet, width = 14.0, height = 12.5, bg = "white", device = cairo_pdf)

# --- 3C. Scheme B Primary Figure 1a: DDR & Lactylated Gene Expression (31 Groups on qsmooth) ---
grp_expr_long <- melt(
  grp_f1a_dt,
  id.vars = c("GroupID", "Category", "CategoryLabel", "X"),
  measure.vars = c("DDRExpressionMedian", "KlaExpressionMedian"),
  variable.name = "GeneSet",
  value.name = "MedianExpression"
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

plot_fig1a_expr_31 <- ggplot(grp_expr_long, aes(x = CategoryLabel, y = MedianExpression, fill = GeneSet)) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5), colour = "#E5E7EB", linetype = "dashed", linewidth = 0.5) +
  geom_boxplot(
    position = position_dodge(width = dodge_w),
    width = 0.58, outlier.shape = NA, colour = charcoal,
    linewidth = 0.55, median.linewidth = 0.75, alpha = 0.82, na.rm = TRUE
  ) +
  geom_segment(
    data = grp_expr_summary,
    aes(x = XPos - 0.16, xend = XPos + 0.16, y = Mean, yend = Mean),
    inherit.aes = FALSE, colour = mean_colour, linewidth = 0.80
  ) +
  geom_point(
    aes(fill = GeneSet),
    position = position_jitterdodge(jitter.width = 0.14, dodge.width = dodge_w, seed = 25),
    shape = 21, size = 2.8, stroke = 0.50, colour = "white", alpha = 0.90, na.rm = TRUE
  ) +
  geom_text(
    data = grp_expr_summary,
    aes(x = XPos, y = y_max_e31 * 0.96, label = paste0("n=", N)),
    inherit.aes = FALSE, size = 3.8, family = publication_font, colour = muted_text, fontface = "bold"
  ) +
  scale_fill_manual(values = c("DDR genes" = "#4E79A7", "Lactylated protein genes" = "#F28E2B")) +
  scale_y_continuous(limits = c(y_min_e31, y_max_e31), breaks = scales::pretty_breaks(n = 6), expand = expansion(mult = c(0, 0))) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE, keyheight = grid::unit(0.55, "cm"), keywidth = grid::unit(0.85, "cm"))) +
  labs(
    title = "Expression of DDR genes and lactylated protein genes across four biological categories (31 groups, qsmooth)",
    subtitle = expr_sub_31,
    x = NULL, y = "Group median log2(qsmooth TPM + 0.5)", fill = NULL
  ) +
  theme_minimal(base_size = 14, base_family = publication_font) +
  theme(
    panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(colour = grid_colour, linewidth = 0.50), panel.border = element_blank(),
    axis.line.x = element_line(colour = "#8C939E", linewidth = 0.60), axis.line.y = element_line(colour = "#8C939E", linewidth = 0.60),
    axis.text.x = element_text(size = 13.0, colour = charcoal, face = "bold", margin = margin(t = 6)),
    axis.text.y = element_text(size = 13.0, colour = charcoal),
    axis.title.y = element_text(size = 15.0, face = "bold", colour = charcoal, margin = margin(r = 10)),
    plot.title = element_text(size = 14.5, face = "bold", colour = charcoal, hjust = 0.5, margin = margin(b = 4)),
    plot.subtitle = element_text(size = 12.0, colour = muted_text, hjust = 0.5, margin = margin(b = 10)),
    legend.position = "top", legend.direction = "horizontal",
    legend.text = element_text(size = 13.0, colour = charcoal),
    legend.key.spacing.x = grid::unit(0.35, "cm"), legend.background = element_rect(fill = "white", colour = NA),
    legend.margin = margin(1, 0, 4, 0),
    plot.margin = margin(10, 16, 12, 12), plot.background = element_rect(fill = "white", colour = NA)
  )

stem_1a_expr_31 <- file.path(out_dir, "Figure_1a_RNA_DDR_and_lactylated_gene_expression_31groups_boxplot")
ggsave(paste0(stem_1a_expr_31, ".png"), plot_fig1a_expr_31, width = 8.5, height = 7.0, dpi = 300, bg = "white", device = ragg::agg_png)
ggsave(paste0(stem_1a_expr_31, ".pdf"), plot_fig1a_expr_31, width = 8.5, height = 7.0, bg = "white", device = cairo_pdf)

# --- 3D. Primary Figure 1a: DDR Gene Fraction Boxplot (31 Groups on qsmooth) ---
grp_frac_counts <- grp_f1a_dt[, .(N = .N), by = .(Category, X)]
grp_frac_stats <- grp_f1a_dt[, .(
  Mean = mean(DdrFractionPercentage),
  Median = median(DdrFractionPercentage)
), by = .(Category, X)]
grp_frac_stats[, c("x_left", "x_right") := list(X - 0.20, X + 0.20)]

aov_frac_31 <- summary(aov(DdrFractionPercentage ~ Category, data = grp_f1a_dt))[[1L]]
sub_frac_31 <- paste0("Four-category one-way ANOVA p = ", format_q_value(aov_frac_31["Category", "Pr(>F)"]), " (F = ", sprintf("%.2f", aov_frac_31["Category", "F value"]), ")")

plot_fig1a_frac_31 <- ggplot(grp_f1a_dt, aes(x = X, y = DdrFractionPercentage, fill = Category)) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5), colour = "#E5E7EB", linetype = "dashed", linewidth = 0.5) +
  geom_boxplot(
    aes(group = Category),
    width = 0.58, outlier.shape = NA, colour = charcoal,
    linewidth = 0.55, median.linewidth = 0.75, alpha = 0.82, orientation = "x", na.rm = TRUE
  ) +
  geom_segment(
    data = grp_frac_stats,
    aes(x = x_left, xend = x_right, y = Median, yend = Median),
    inherit.aes = FALSE, colour = charcoal, linewidth = 0.75, lineend = "round"
  ) +
  geom_segment(
    data = grp_frac_stats,
    aes(x = x_left, xend = x_right, y = Mean, yend = Mean),
    inherit.aes = FALSE, colour = mean_colour, linewidth = 0.80, lineend = "round"
  ) +
  geom_point(
    aes(group = Category),
    position = position_jitter(width = 0.095, height = 0, seed = 25),
    shape = 21, size = 3.0, stroke = 0.50, colour = "white", alpha = 0.90, na.rm = TRUE
  ) +
  geom_text(
    data = grp_frac_counts,
    aes(x = X, y = 3.8, label = paste0("n=", N)),
    inherit.aes = FALSE, family = publication_font, size = 4.2, fontface = "bold",
    colour = muted_text, vjust = -0.15
  ) +
  scale_fill_manual(values = category_fills, guide = "none", drop = FALSE) +
  scale_x_continuous(
    breaks = seq_along(category_order),
    labels = unname(category_labels[category_order]),
    limits = c(0.5, length(category_order) + 0.5),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(2.0, 4.0),
    breaks = seq(2.0, 4.0, by = 0.5),
    labels = function(y) paste0(y, "%"),
    expand = c(0, 0)
  ) +
  labs(
    title = "DDR annotated gene fraction in the transcriptome across four biological categories (31 groups, qsmooth)",
    subtitle = sub_frac_31,
    x = NULL,
    y = "GO-DDR annotated expressed gene fraction (%)"
  ) +
  theme_minimal(base_family = publication_font, base_size = 14) +
  theme(
    panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(colour = grid_colour, linewidth = 0.50), panel.border = element_blank(),
    axis.line.x = element_line(colour = "#8C939E", linewidth = 0.60), axis.line.y = element_line(colour = "#8C939E", linewidth = 0.60),
    axis.text.x = element_text(size = 13.5, face = "bold", colour = charcoal, lineheight = 0.95, margin = margin(t = 6)),
    axis.text.y = element_text(size = 13.0, colour = charcoal),
    axis.title.y = element_text(size = 15.0, face = "bold", colour = charcoal, margin = margin(r = 10)),
    plot.title = element_text(size = 15.0, face = "bold", colour = charcoal, hjust = 0.5, margin = margin(b = 4)),
    plot.subtitle = element_text(size = 12.0, colour = muted_text, hjust = 0.5, margin = margin(b = 10)),
    plot.margin = margin(12, 16, 12, 12), plot.background = element_rect(fill = "white", colour = NA)
  ) +
  coord_cartesian(clip = "off")

stem_1a_frac_31 <- file.path(out_dir, "Figure_1a_RNA_DDR_fraction_31groups_boxplot")
ggsave(paste0(stem_1a_frac_31, ".png"), plot_fig1a_frac_31, width = 8.5, height = 7.0, dpi = 300, bg = "white", device = ragg::agg_png)
ggsave(paste0(stem_1a_frac_31, ".pdf"), plot_fig1a_frac_31, width = 8.5, height = 7.0, bg = "white", device = cairo_pdf)

# --- 3E. Sample-Level DDR Companions (1,898 samples on qsmooth) ---
sample_expr_long <- melt(
  sample_gene_stats,
  id.vars = c("SampleID", "Category"),
  measure.vars = c("DDRExpressionMedian", "KlaExpressionMedian"),
  variable.name = "GeneSet", value.name = "MedianExpression"
)
sample_expr_long[, GeneSet := factor(GeneSet, levels = c("DDRExpressionMedian", "KlaExpressionMedian"), labels = c("DDR genes", "Lactylated protein genes"))]
sample_expr_long[, CategoryLabel := factor(Category, levels = category_order, labels = unname(category_labels[category_order]))]
sample_expr_summary <- sample_expr_long[, .(Mean = mean(MedianExpression), Median = median(MedianExpression), N = .N), by = .(Category, CategoryLabel, GeneSet)]
sample_expr_summary[, CatIdx := as.numeric(CategoryLabel)]
sample_expr_summary[, XPos := fifelse(GeneSet == "DDR genes", CatIdx - dodge_w / 4, CatIdx + dodge_w / 4)]
y_max_esamp <- max(sample_expr_long$MedianExpression) + 0.85

plot_fig1a_expr_samp <- ggplot(sample_expr_long, aes(x = CategoryLabel, y = MedianExpression, fill = GeneSet)) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5), colour = "#E5E7EB", linetype = "dashed", linewidth = 0.5) +
  geom_boxplot(position = position_dodge(width = dodge_w), width = 0.58, outlier.shape = NA, colour = charcoal, linewidth = 0.55, median.linewidth = 0.75, alpha = 0.82, na.rm = TRUE) +
  geom_segment(data = sample_expr_summary, aes(x = XPos - 0.16, xend = XPos + 0.16, y = Mean, yend = Mean), inherit.aes = FALSE, colour = mean_colour, linewidth = 0.80) +
  geom_point(aes(fill = GeneSet), position = position_jitterdodge(jitter.width = 0.14, dodge.width = dodge_w, seed = 25), shape = 21, size = 2.0, stroke = 0.35, colour = "white", alpha = 0.65, na.rm = TRUE) +
  geom_text(data = sample_expr_summary, aes(x = XPos, y = y_max_esamp * 0.96, label = paste0("n=", N)), inherit.aes = FALSE, size = 3.8, family = publication_font, colour = muted_text, fontface = "bold") +
  scale_fill_manual(values = c("DDR genes" = "#4E79A7", "Lactylated protein genes" = "#F28E2B")) +
  scale_y_continuous(limits = c(min(sample_expr_long$MedianExpression) - 0.2, y_max_esamp), breaks = scales::pretty_breaks(n = 6), expand = expansion(mult = c(0, 0))) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE, keyheight = grid::unit(0.55, "cm"), keywidth = grid::unit(0.85, "cm"))) +
  labs(title = "Sample-level expression of DDR and lactylated protein genes (qsmooth)", subtitle = "Two-way ANOVA across 1,898 samples", x = NULL, y = "Sample median log2(qsmooth TPM + 0.5)", fill = NULL) +
  theme_minimal(base_size = 14, base_family = publication_font) +
  theme(
    panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(colour = grid_colour, linewidth = 0.50), panel.border = element_blank(),
    axis.line.x = element_line(colour = "#8C939E", linewidth = 0.60), axis.line.y = element_line(colour = "#8C939E", linewidth = 0.60),
    axis.text.x = element_text(size = 13.0, colour = charcoal, face = "bold", margin = margin(t = 6)),
    axis.text.y = element_text(size = 13.0, colour = charcoal), axis.title.y = element_text(size = 15.0, face = "bold", colour = charcoal, margin = margin(r = 10)),
    plot.title = element_text(size = 15.5, face = "bold", colour = charcoal, hjust = 0.5, margin = margin(b = 4)),
    plot.subtitle = element_text(size = 12.0, colour = muted_text, hjust = 0.5, margin = margin(b = 10)),
    legend.position = "top", legend.direction = "horizontal", legend.text = element_text(size = 13.0, colour = charcoal),
    legend.key.spacing.x = grid::unit(0.35, "cm"), legend.background = element_rect(fill = "white", colour = NA),
    legend.margin = margin(1, 0, 4, 0), plot.margin = margin(10, 16, 12, 12), plot.background = element_rect(fill = "white", colour = NA)
  )

stem_1a_expr_samp <- file.path(out_dir, "Figure_1a_companion_RNA_DDR_and_lactylated_gene_expression_sample_boxplot")
ggsave(paste0(stem_1a_expr_samp, ".png"), plot_fig1a_expr_samp, width = 8.5, height = 7.0, dpi = 300, bg = "white", device = ragg::agg_png)
ggsave(paste0(stem_1a_expr_samp, ".pdf"), plot_fig1a_expr_samp, width = 8.5, height = 7.0, bg = "white", device = cairo_pdf)

sample_frac_dt <- copy(sample_gene_stats)
sample_frac_dt[, CategoryLabel := factor(Category, levels = category_order, labels = unname(category_labels[category_order]))]
sample_frac_dt[, X := match(Category, category_order)]
sample_frac_counts <- sample_frac_dt[, .(N = .N), by = .(Category, X)]
sample_frac_stats <- sample_frac_dt[, .(Mean = mean(DdrFractionPercentage), Median = median(DdrFractionPercentage)), by = .(Category, X)]
sample_frac_stats[, c("x_left", "x_right") := list(X - 0.20, X + 0.20)]
aov_frac_s <- summary(aov(DdrFractionPercentage ~ Category, data = sample_frac_dt))[[1L]]
sub_frac_s <- paste0("Four-category one-way ANOVA p = ", format_q_value(aov_frac_s["Category", "Pr(>F)"]), " (F = ", sprintf("%.2f", aov_frac_s["Category", "F value"]), ")")

plot_fig1a_frac_samp <- ggplot(sample_frac_dt, aes(x = X, y = DdrFractionPercentage, fill = Category)) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5), colour = "#E5E7EB", linetype = "dashed", linewidth = 0.5) +
  geom_boxplot(aes(group = Category), width = 0.58, outlier.shape = NA, colour = charcoal, linewidth = 0.55, median.linewidth = 0.75, alpha = 0.82, orientation = "x", na.rm = TRUE) +
  geom_segment(data = sample_frac_stats, aes(x = x_left, xend = x_right, y = Median, yend = Median), inherit.aes = FALSE, colour = charcoal, linewidth = 0.75, lineend = "round") +
  geom_segment(data = sample_frac_stats, aes(x = x_left, xend = x_right, y = Mean, yend = Mean), inherit.aes = FALSE, colour = mean_colour, linewidth = 0.80, lineend = "round") +
  geom_point(aes(group = Category), position = position_jitter(width = 0.095, height = 0, seed = 25), shape = 21, size = 2.0, stroke = 0.35, colour = "white", alpha = 0.65, na.rm = TRUE) +
  geom_text(data = sample_frac_counts, aes(x = X, y = 4.8, label = paste0("n=", N)), inherit.aes = FALSE, family = publication_font, size = 4.0, fontface = "bold", colour = muted_text, vjust = -0.15) +
  scale_fill_manual(values = category_fills, guide = "none", drop = FALSE) +
  scale_x_continuous(breaks = seq_along(category_order), labels = unname(category_labels[category_order]), limits = c(0.5, length(category_order) + 0.5), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 5.2), breaks = seq(0, 5, by = 1), labels = function(y) paste0(y, "%"), expand = c(0, 0)) +
  labs(title = "Sample-level DDR annotated expressed gene fraction (qsmooth)", subtitle = sub_frac_s, x = NULL, y = "GO-DDR annotated expressed gene fraction (%)") +
  theme_minimal(base_family = publication_font, base_size = 14) +
  theme(
    panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(colour = grid_colour, linewidth = 0.50), panel.border = element_blank(),
    axis.line.x = element_line(colour = "#8C939E", linewidth = 0.60), axis.line.y = element_line(colour = "#8C939E", linewidth = 0.60),
    axis.text.x = element_text(size = 13.5, face = "bold", colour = charcoal, lineheight = 0.95, margin = margin(t = 6)),
    axis.text.y = element_text(size = 13.0, colour = charcoal),
    axis.title.y = element_text(size = 15.0, face = "bold", colour = charcoal, margin = margin(r = 10)),
    plot.title = element_text(size = 16.0, face = "bold", colour = charcoal, hjust = 0.5, margin = margin(b = 4)),
    plot.subtitle = element_text(size = 12.0, colour = muted_text, hjust = 0.5, margin = margin(b = 10)),
    plot.margin = margin(12, 16, 12, 12), plot.background = element_rect(fill = "white", colour = NA)
  ) + coord_cartesian(clip = "off")

stem_1a_frac_samp <- file.path(out_dir, "Figure_1a_companion_RNA_DDR_fraction_sample_boxplot")
ggsave(paste0(stem_1a_frac_samp, ".png"), plot_fig1a_frac_samp, width = 8.5, height = 7.0, dpi = 300, bg = "white", device = ragg::agg_png)
ggsave(paste0(stem_1a_frac_samp, ".pdf"), plot_fig1a_frac_samp, width = 8.5, height = 7.0, bg = "white", device = cairo_pdf)

# --- 3F. Per-Tissue DDR Fraction Companion ---
grp_frac_by_tissue_summary <- sample_31_ddr_dt[, .(
  N = .N,
  Median = median(DdrFractionPercentage),
  Mean = mean(DdrFractionPercentage)
), by = .(GroupID, FactorLabel, Category, CategoryLabel, Y_Pos)]

plot_ddr_frac_by_tissue <- ggplot(sample_31_ddr_dt, aes(x = DdrFractionPercentage, y = FactorLabel, fill = Category)) +
  geom_boxplot(aes(group = FactorLabel), width = 0.58, outlier.shape = NA, colour = charcoal, linewidth = 0.50, median.linewidth = 0.75, alpha = 0.80) +
  geom_point(position = position_jitter(width = 0, height = 0.18, seed = 25), shape = 21, size = 1.9, stroke = 0.30, colour = "white", alpha = 0.65) +
  geom_point(data = grp_frac_by_tissue_summary, aes(x = Mean, y = FactorLabel), inherit.aes = FALSE, shape = 23, size = 2.4, fill = mean_colour, colour = charcoal, stroke = 0.45) +
  scale_fill_manual(values = category_fills, guide = "none") +
  scale_x_continuous(labels = function(x) paste0(x, "%"), expand = expansion(mult = c(0.02, 0.05))) +
  labs(
    title = "Transcriptomic DDR gene fraction across 31 biological materials (sample points on qsmooth)",
    subtitle = "Percentage of expressed genes with GO DDR annotation | Red diamonds = mean | Boxplots = median & IQR",
    x = "DDR expressed gene fraction (%)",
    y = NULL
  ) +
  theme_minimal(base_family = publication_font, base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(colour = "#ECEFF1", linewidth = 0.45),
    panel.grid.major.x = element_line(colour = grid_colour, linewidth = 0.50),
    axis.text.y = element_text(size = 9.5, colour = charcoal),
    axis.text.x = element_text(size = 11.0, colour = charcoal),
    axis.title.x = element_text(size = 12.5, face = "bold", colour = charcoal, margin = margin(t = 8)),
    plot.title = element_text(size = 14.5, face = "bold", colour = charcoal, hjust = 0),
    plot.subtitle = element_text(size = 10.5, colour = muted_text, hjust = 0, margin = margin(b = 10)),
    plot.margin = margin(12, 16, 12, 12)
  )

stem_ddr_frac_by_tissue <- file.path(out_dir, "Figure_1a_companion_RNA_DDR_fraction_by_tissue_boxplot")
ggsave(paste0(stem_ddr_frac_by_tissue, ".png"), plot_ddr_frac_by_tissue, width = 14.0, height = 11.5, dpi = 300, bg = "white", device = ragg::agg_png)
ggsave(paste0(stem_ddr_frac_by_tissue, ".pdf"), plot_ddr_frac_by_tissue, width = 14.0, height = 11.5, bg = "white", device = cairo_pdf)

message(">>> Task 3 complete: Figure 1a DDR gene expression figures saved.")

# ==============================================================================
# SECTION 5: SYNCHRONIZATION TO DELIVERY DIRECTORIES
# ==============================================================================
message(">>> Synchronizing all newly generated deliverables...")

# 1. Local delivery dir
deliv_dir <- file.path(root, "outputs", "20260917_delivery_31group", "rna")
if (dir.exists(deliv_dir)) {
  file.copy(list.files(out_dir, full.names = TRUE), deliv_dir, overwrite = TRUE)
}

# 2. External desktop delivery dirs
ext_renew_rna <- "/Users/gzy2520/Desktop/renew/kla/rna"
ext_renew_top <- "/Users/gzy2520/Desktop/renew/kla"

if (dir.exists(ext_renew_rna)) {
  # Copy all deliverables to ext_renew_rna
  file.copy(list.files(out_dir, full.names = TRUE), ext_renew_rna, overwrite = TRUE)
  
  # Copy primary publication PNGs to ext_renew_top
  top_png_names <- c(
    "Figure_1b_RNA_proliferation_by_tissue_sample_boxplot.png",
    "Figure_1b_RNA_proliferation_score_31groups_boxplot.png",
    "Figure_1b_RNA_proliferation_ranking_31groups.png",
    "Figure_1a_RNA_DDR_expression_by_tissue_boxplot.png",
    "Figure_1a_RNA_DDR_and_lactylated_gene_expression_31groups_boxplot.png",
    "Figure_1a_RNA_DDR_fraction_31groups_boxplot.png",
    "Figure_3c_RNA_regulator_qsmooth_heatmap.png",
    "Figure_3c_RNA_regulator_qsmooth_heatmap_no_frame.png",
    "Figure_3c_RNA_regulator_qsmooth_zscore_heatmap.png",
    "Figure_3c_RNA_regulator_qsmooth_zscore_heatmap_no_frame.png",
    "Figure_3c_RNA_regulator_percentiles.png"
  )
  for (tp in top_png_names) {
    src_tp <- file.path(out_dir, tp)
    if (file.exists(src_tp)) file.copy(src_tp, ext_renew_top, overwrite = TRUE)
  }
}

cat("\nALL_TEACHER_QUESTIONS_REBUILT_WITH_QSMOOTH_AND_DELIVERED_SUCCESSFULLY!\n")
