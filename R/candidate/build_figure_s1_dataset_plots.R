#!/usr/bin/env Rscript

# Render Supplementary Figure S1a and S1b:
# Fig S1a: GO-DDR annotated protein fraction (%) by individual PXD dataset
#          (Whole proteome vs Lactylome Kla), with 4 biological category ranges
#          indicated below the x-axis.
# Fig S1b: Whole-proteome MKI67 / H3C1 intensity ratio (log scale) by individual
#          PXD dataset, with 4 biological category ranges indicated below the x-axis.
# Both figures maintain strictly identical proportions, canvas dimensions, and style.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

set.seed(25)

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1L) normalizePath(args[[1L]], mustWork = TRUE) else normalizePath(".", mustWork = TRUE)
stop_if <- function(condition, message) if (!isTRUE(condition)) stop(message, call. = FALSE)

publication_font <- "Arial Unicode MS"
whole_proteome_colour <- "#4E79A7"
kla_colour <- "#F28E2B"
mean_colour <- "#C0392B"
charcoal <- "#2F3437"
muted_text <- "#65717D"
grid_colour <- "#D9DDE3"

category_order <- c("normal_tissue", "cancer_tissue", "normal_cells", "cancer_cells")
category_labels <- c(
  normal_tissue = "non-tumor tissues",
  cancer_tissue = "tumor tissues",
  normal_cells = "normal cell lines",
  cancer_cells = "cancer cell lines"
)
category_colours <- c(
  normal_tissue = "#0072B2",
  cancer_tissue = "#D55E00",
  normal_cells = "#009E73",
  cancer_cells = "#CC79A7"
)
category_strip_fills <- c(
  normal_tissue = "#DCE9E2",
  cancer_tissue = "#F0DEDE",
  normal_cells = "#E7E1EE",
  cancer_cells = "#EEE4D2"
)

candidate_dir <- normalizePath(Sys.getenv(
  "KLA_CANDIDATE_INPUT",
  unset = file.path(project_root, "data", "candidate", "escc_inclusion_20260903_pxd065830_tumor_reference", "candidate_input")
), mustWork = TRUE)

output_dir <- Sys.getenv(
  "KLA_S1_OUTPUT",
  unset = file.path(project_root, "results", "final_figures_and_tables", "Supplementary_Figure_S1")
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# =========================================================================
# 1. Figure S1a: DDR Fraction by PXD Dataset
# =========================================================================
message(">>> Generating Figure S1a (DDR fraction by PXD dataset)...")

values_s1a <- fread(file.path(candidate_dir, "figure1_sample_boxplot_values.csv"), check.names = FALSE)
stop_if(nrow(values_s1a) == 310L, "Figure S1a expects 310 source sample observations.")

# Establish group order: Category order, then study order
group_meta <- values_s1a[, .(
  Category = first(Category),
  SampleGroup = first(SampleGroup),
  RowOrder = min(RowOrder)
), by = .(PXD, SampleGroup)]
group_meta[, CategoryOrder := match(Category, category_order)]
setorder(group_meta, CategoryOrder, RowOrder, PXD)
group_meta[, GroupIndex := .I]

# Create unique label for each dataset on x-axis
# If PXD appears multiple times within category, add short disambiguation
group_meta[, PxdCountInCat := .N, by = .(Category, PXD)]
group_meta[, DatasetXLabel := fifelse(
  PxdCountInCat > 1L,
  paste0(PXD, "\n(", substr(SampleGroup, 1, 10), ")"),
  PXD
)]
# Clean specific labels for clarity
group_meta[SampleGroup == "hypertrophic scar", DatasetXLabel := paste0(PXD, "\n(scar)")]
group_meta[SampleGroup == "adjacent skin", DatasetXLabel := paste0(PXD, "\n(skin)")]
group_meta[SampleGroup == "MCF7", DatasetXLabel := paste0(PXD, "\n(MCF7)")]
group_meta[SampleGroup == "MDA-MB-468", DatasetXLabel := paste0(PXD, "\n(MDA468)")]
group_meta[SampleGroup == "T-47D", DatasetXLabel := paste0(PXD, "\n(T47D)")]
group_meta[SampleGroup == "HCT116" & RowOrder == 15, DatasetXLabel := paste0(PXD, "\n(HCT116)")]
group_meta[SampleGroup == "TALL-104", DatasetXLabel := paste0(PXD, "\n(TALL104)")]

values_s1a <- merge(values_s1a, group_meta[, .(PXD, SampleGroup, GroupIndex, DatasetXLabel)], by = c("PXD", "SampleGroup"))
setorder(values_s1a, GroupIndex)
values_s1a[, DatasetXFactor := factor(GroupIndex, levels = group_meta$GroupIndex, labels = group_meta$DatasetXLabel)]
values_s1a[, Dataset := factor(Dataset, levels = c("Whole proteome", "Lactylome (Kla)"))]

# Compute summary statistics
s1a_stats <- values_s1a[, .(
  Mean = mean(DdrFractionPercentage),
  Median = median(DdrFractionPercentage),
  N = .N
), by = .(GroupIndex, DatasetXFactor, Dataset, Category)]

dodge_w <- 0.72
s1a_stats[, CatX := as.numeric(DatasetXFactor)]
s1a_stats[, XPos := fifelse(Dataset == "Whole proteome", CatX - dodge_w / 4, CatX + dodge_w / 4)]

y_max_s1a <- max(22, ceiling((max(values_s1a$DdrFractionPercentage) + 2) / 5) * 5)

# Category ranges along x-axis for range bar
cat_ranges_s1a <- group_meta[, .(
  XMin = min(GroupIndex) - 0.45,
  XMax = max(GroupIndex) + 0.45,
  XMid = (min(GroupIndex) + max(GroupIndex)) / 2,
  CategoryLabel = unname(category_labels[first(Category)])
), by = Category]
cat_ranges_s1a[, CatOrder := match(Category, category_order)]
setorder(cat_ranges_s1a, CatOrder)

# Category boundary separators
cat_dividers_s1a <- cat_ranges_s1a$XMax[-nrow(cat_ranges_s1a)] + 0.05

p_s1a <- ggplot(values_s1a, aes(x = DatasetXFactor, y = DdrFractionPercentage, fill = Dataset)) +
  geom_vline(xintercept = cat_dividers_s1a, colour = "#CCD1D9", linetype = "dashed", linewidth = 0.55) +
  geom_boxplot(
    position = position_dodge(width = dodge_w),
    width = 0.58, outlier.shape = NA, colour = charcoal,
    linewidth = 0.65, median.linewidth = 1.1, alpha = 0.85, na.rm = TRUE
  ) +
  geom_segment(
    data = s1a_stats,
    aes(x = XPos - 0.13, xend = XPos + 0.13, y = Mean, yend = Mean),
    inherit.aes = FALSE, colour = mean_colour, linewidth = 1.25
  ) +
  geom_point(
    aes(fill = Dataset),
    position = position_jitterdodge(jitter.width = 0.12, dodge.width = dodge_w, seed = 25),
    shape = 21, size = 1.9, stroke = 0.45, colour = "white", alpha = 0.85, na.rm = TRUE
  ) +
  scale_fill_manual(values = c("Whole proteome" = whole_proteome_colour, "Lactylome (Kla)" = kla_colour)) +
  scale_y_continuous(
    limits = c(0, y_max_s1a), breaks = seq(0, y_max_s1a, by = 5),
    labels = function(y) paste0(y, "%"), expand = expansion(mult = c(0, 0))
  ) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE, keyheight = grid::unit(0.55, "cm"), keywidth = grid::unit(0.9, "cm"))) +
  labs(
    title = "GO-DDR annotated protein fraction across 31 individual datasets",
    x = NULL, y = "GO-DDR annotated protein fraction (%)", fill = NULL
  ) +
  theme_minimal(base_size = 13, base_family = publication_font) +
  theme(
    panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(colour = grid_colour, linewidth = 0.45),
    panel.border = element_blank(),
    axis.line.x = element_line(colour = "#8C939E", linewidth = 0.60),
    axis.line.y = element_line(colour = "#8C939E", linewidth = 0.60),
    axis.text.x = element_text(size = 9.5, colour = charcoal, face = "bold", angle = 50, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 11.5, colour = charcoal),
    axis.title.y = element_text(size = 14, face = "bold", colour = charcoal, margin = margin(r = 10)),
    plot.title = element_text(size = 16, face = "bold", colour = charcoal, hjust = 0.5, margin = margin(b = 4)),
    legend.position = "top", legend.direction = "horizontal",
    legend.text = element_text(size = 12.5, colour = charcoal),
    legend.margin = margin(0, 0, 4, 0),
    plot.margin = margin(8, 16, 2, 12),
    plot.background = element_rect(fill = "white", colour = NA)
  )

# Bottom category range track for S1a
p_s1a_track <- ggplot(cat_ranges_s1a) +
  geom_rect(
    aes(xmin = XMin, xmax = XMax, ymin = 0.15, ymax = 0.85, fill = Category),
    colour = charcoal, linewidth = 0.45, alpha = 0.85
  ) +
  geom_text(
    aes(x = XMid, y = 0.5, label = CategoryLabel),
    family = publication_font, size = 3.8, fontface = "bold", colour = "white"
  ) +
  scale_fill_manual(values = category_colours, guide = "none") +
  scale_x_continuous(limits = c(0.5, nrow(group_meta) + 0.5), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  theme_void() +
  theme(plot.margin = margin(0, 16, 8, 12))

fig_s1a_combined <- p_s1a / p_s1a_track + plot_layout(heights = c(12, 1.1))

# =========================================================================
# 2. Figure S1b: MKI67 / H3C1 Ratio by PXD Dataset
# =========================================================================
message(">>> Generating Figure S1b (MKI67 / H3C1 ratio by PXD dataset)...")

values_s1b <- fread(file.path(candidate_dir, "figure1_mki67_ratio_sample_values.csv"), check.names = FALSE)
values_s1b <- values_s1b[Denominator == "H3C1"]
stop_if(nrow(values_s1b) == 80L, "Figure S1b expects 80 source sample observations with H3C1 denominator.")

group_meta_b <- values_s1b[, .(
  Category = first(Category),
  SampleGroup = first(SampleGroup)
), by = .(PXD, SampleGroup)]
group_meta_b[, CategoryOrder := match(Category, category_order)]
setorder(group_meta_b, CategoryOrder, PXD, SampleGroup)
group_meta_b[, GroupIndex := .I]

group_meta_b[, DatasetXLabel := fifelse(
  PXD == "PXD066054" & Category == "normal_tissue", "PXD066054\n(BPH)",
  fifelse(
    PXD == "PXD066054" & Category == "cancer_tissue", "PXD066054\n(prostate)",
    fifelse(PXD == "PXD064038", "PXD064038\n(ESCC)", PXD)
  )
)]

values_s1b <- merge(values_s1b, group_meta_b[, .(PXD, SampleGroup, GroupIndex, DatasetXLabel)], by = c("PXD", "SampleGroup"))
setorder(values_s1b, GroupIndex)
values_s1b[, DatasetXFactor := factor(GroupIndex, levels = group_meta_b$GroupIndex, labels = group_meta_b$DatasetXLabel)]
values_s1b[, Category := factor(Category, levels = category_order)]

s1b_stats <- values_s1b[, .(
  Mean = mean(Ratio),
  Median = median(Ratio),
  N = .N,
  MaxRatio = max(Ratio)
), by = .(GroupIndex, DatasetXFactor, Category)]
s1b_stats[, label_y := 10^(log10(MaxRatio) + 0.32)]

raw_max_b <- max(values_s1b$Ratio, na.rm = TRUE)
raw_min_b <- min(values_s1b$Ratio, na.rm = TRUE)
y_min_b <- 10^(floor(log10(raw_min_b)) - 0.35)
y_max_b <- 10^(log10(raw_max_b) + 0.75)

cat_ranges_s1b <- group_meta_b[, .(
  XMin = min(GroupIndex) - 0.45,
  XMax = max(GroupIndex) + 0.45,
  XMid = (min(GroupIndex) + max(GroupIndex)) / 2,
  CategoryLabel = unname(category_labels[first(Category)])
), by = Category]
cat_ranges_s1b[, CatOrder := match(Category, category_order)]
setorder(cat_ranges_s1b, CatOrder)
cat_dividers_s1b <- cat_ranges_s1b$XMax[-nrow(cat_ranges_s1b)] + 0.05

p_s1b <- ggplot(values_s1b, aes(x = DatasetXFactor, y = Ratio, fill = Category)) +
  geom_vline(xintercept = cat_dividers_s1b, colour = "#CCD1D9", linetype = "dashed", linewidth = 0.55) +
  geom_boxplot(
    aes(group = DatasetXFactor),
    width = 0.52, outlier.shape = NA, colour = charcoal,
    linewidth = 0.70, median.linewidth = 1.25, alpha = 0.85, na.rm = TRUE
  ) +
  geom_segment(
    data = s1b_stats,
    aes(x = as.numeric(DatasetXFactor) - 0.20, xend = as.numeric(DatasetXFactor) + 0.20, y = Mean, yend = Mean),
    inherit.aes = FALSE, colour = mean_colour, linewidth = 1.35
  ) +
  geom_point(
    aes(group = DatasetXFactor),
    position = position_jitter(width = 0.12, height = 0, seed = 25),
    shape = 21, size = 2.4, stroke = 0.55, colour = "white", alpha = 0.90, na.rm = TRUE
  ) +
  geom_text(
    data = s1b_stats,
    aes(x = as.numeric(DatasetXFactor), y = label_y, label = paste0("n=", N)),
    inherit.aes = FALSE, family = publication_font, size = 3.6, fontface = "bold", colour = muted_text, vjust = -0.15
  ) +
  scale_fill_manual(values = category_colours, guide = "none") +
  scale_y_log10(
    breaks = 10^seq(floor(log10(y_min_b)), ceiling(log10(y_max_b)), by = 1),
    labels = scales::label_scientific(digits = 2),
    limits = c(y_min_b, y_max_b),
    expand = c(0, 0)
  ) +
  labs(
    title = "MKI67 / H3C1 protein intensity ratio across 11 individual datasets",
    x = NULL, y = "MKI67 / H3C1 ratio (log scale)"
  ) +
  theme_minimal(base_size = 13, base_family = publication_font) +
  theme(
    panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(colour = grid_colour, linewidth = 0.45),
    panel.border = element_blank(),
    axis.line.x = element_line(colour = "#8C939E", linewidth = 0.60),
    axis.line.y = element_line(colour = "#8C939E", linewidth = 0.60),
    axis.text.x = element_text(size = 10.5, colour = charcoal, face = "bold", angle = 50, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 11.5, colour = charcoal),
    axis.title.y = element_text(size = 14, face = "bold", colour = charcoal, margin = margin(r = 10)),
    plot.title = element_text(size = 16, face = "bold", colour = charcoal, hjust = 0.5, margin = margin(b = 4)),
    plot.margin = margin(8, 16, 2, 12),
    plot.background = element_rect(fill = "white", colour = NA)
  )

p_s1b_track <- ggplot(cat_ranges_s1b) +
  geom_rect(
    aes(xmin = XMin, xmax = XMax, ymin = 0.15, ymax = 0.85, fill = Category),
    colour = charcoal, linewidth = 0.45, alpha = 0.85
  ) +
  geom_text(
    aes(x = XMid, y = 0.5, label = CategoryLabel),
    family = publication_font, size = 3.8, fontface = "bold", colour = "white"
  ) +
  scale_fill_manual(values = category_colours, guide = "none") +
  scale_x_continuous(limits = c(0.5, nrow(group_meta_b) + 0.5), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  theme_void() +
  theme(plot.margin = margin(0, 16, 8, 12))

fig_s1b_combined <- p_s1b / p_s1b_track + plot_layout(heights = c(12, 1.1))

# =========================================================================
# 3. Save Both Figures with Exact Same Proportions (Consistent Ratio)
# =========================================================================
fig_w <- 15.0
fig_h <- 7.5

s1a_stem <- file.path(output_dir, "Figure_S1a_DDR_fraction_by_PXD")
ggsave(paste0(s1a_stem, ".png"), fig_s1a_combined, width = fig_w, height = fig_h, dpi = 300, bg = "white", device = ragg::agg_png)
ggsave(paste0(s1a_stem, ".pdf"), fig_s1a_combined, width = fig_w, height = fig_h, bg = "white", device = cairo_pdf)
message("Wrote Figure S1a to ", s1a_stem, ".png/.pdf")

s1b_stem <- file.path(output_dir, "Figure_S1b_MKI67_over_H3C1_by_PXD")
ggsave(paste0(s1b_stem, ".png"), fig_s1b_combined, width = fig_w, height = fig_h, dpi = 300, bg = "white", device = ragg::agg_png)
ggsave(paste0(s1b_stem, ".pdf"), fig_s1b_combined, width = fig_w, height = fig_h, bg = "white", device = cairo_pdf)
message("Wrote Figure S1b to ", s1b_stem, ".png/.pdf")

