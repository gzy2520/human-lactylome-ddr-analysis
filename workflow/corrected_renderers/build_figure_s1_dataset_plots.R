#!/usr/bin/env Rscript

# Render Supplementary Figure S1a and S1b:
# Fig S1a: GO-DDR annotated protein fraction (%) by individual tissue / cell line dataset
#          (Whole proteome vs Lactylome Kla), with 4 biological category banners at the top.
# Fig S1b: Whole-proteome MKI67 / H3C1 intensity ratio (log scale) across individual datasets,
#          aligned with Fig S1a (detected datasets show boxplots, undetected datasets show "ND"),
#          plus an 11-detected-only companion output.
# Both figures maintain strictly identical proportions, canvas dimensions (15.0 x 7.5 in), and style.

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
  normal_tissue = "Non-tumor tissues",
  cancer_tissue = "Tumor tissues",
  normal_cells = "Normal cell lines",
  cancer_cells = "Cancer cell lines"
)
category_colours <- c(
  normal_tissue = "#0072B2",
  cancer_tissue = "#D55E00",
  normal_cells = "#009E73",
  cancer_cells = "#CC79A7"
)

tissue_names <- c(
  "pathological rotator cuff tendon" = "Rotator cuff tendon",
  "normal human lung" = "Lung",
  "hypertrophic scar" = "Hypertrophic scar",
  "adjacent skin" = "Adjacent skin",
  "human hippocampus" = "Hippocampus",
  "normal pregnancy placenta" = "Placenta",
  "human sperm" = "Sperm",
  "BPH" = "BPH",
  "adjacent liver" = "Adjacent liver",
  "MEC and NEC ESCC groups" = "ESCC",
  "prostate cancer" = "Prostate cancer",
  "HCC" = "HCC",
  "HEK293T" = "HEK293T",
  "HMC3" = "HMC3",
  "pretreated HK-2" = "HK-2 (pretreated)",
  "HK-2 control and mannitol" = "HK-2 (mannitol)",
  "MCF10A" = "MCF10A",
  "neural stem cells" = "Neural stem cells",
  "HUVEC control and Pg infection" = "HUVEC",
  "MCF7" = "MCF7",
  "HCT116" = "HCT116",
  "HCT116 control and Roseburia co-culture" = "HCT116 (Roseburia)",
  "TALL-104" = "TALL-104",
  "HepG2 WT and SIRT1 or SIRT3 KO" = "HepG2",
  "A549" = "A549",
  "MDA-MB-468" = "MDA-MB-468",
  "T-47D" = "T-47D",
  "PC-3M" = "PC-3M",
  "glioblastoma stem cells" = "Glioblastoma stem cells",
  "RKO WT and GSK3B KO" = "RKO"
)

candidate_dir <- normalizePath(Sys.getenv(
  "KLA_CANDIDATE_INPUT",
  unset = file.path(project_root, "data", "candidate", "escc_inclusion_20260903_pxd065830_tumor_reference", "candidate_input")
), mustWork = TRUE)

sample_only <- tolower(Sys.getenv("KLA_SAMPLE_ONLY", unset = "TRUE")) %in% c("1", "true", "yes")
full_figure1_input <- Sys.getenv("KLA_FULL_FIGURE1_INPUT", unset = "")

output_dir <- Sys.getenv(
  "KLA_S1_OUTPUT",
  unset = file.path(project_root, "results", "final_figures_and_tables", "Supplementary_Figure_S1")
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# =========================================================================
# 1. Master 31-Dataset Metadata
# =========================================================================
values_s1a_raw <- fread(file.path(candidate_dir, "figure1_sample_boxplot_values.csv"), check.names = FALSE)
stop_if("ObservationType" %in% names(values_s1a_raw), "Figure S1a input lacks ObservationType.")
metadata_values_s1a <- values_s1a_raw
if (sample_only && nzchar(full_figure1_input) && file.exists(full_figure1_input)) {
  metadata_values_s1a <- fread(full_figure1_input, check.names = FALSE)
}
values_s1a <- if (sample_only) values_s1a_raw[ObservationType == "sample"] else values_s1a_raw
expected_s1a_n <- as.integer(Sys.getenv("KLA_S1A_EXPECTED_ROWS", unset = as.character(nrow(values_s1a))))
stop_if(nrow(values_s1a) == expected_s1a_n,
  paste0("Figure S1a expects ", expected_s1a_n, " observations after the configured filter."))

group_meta <- metadata_values_s1a[, .(
  Category = first(Category),
  RowOrder = min(RowOrder)
), by = .(PXD, SampleGroup)]
group_meta[, CategoryOrder := match(Category, category_order)]
setorder(group_meta, CategoryOrder, RowOrder, PXD)
group_meta[, GroupIndex := .I]

group_meta[, CleanName := unname(tissue_names[SampleGroup])]
group_meta[SampleGroup == "HCT116" & RowOrder == 15, CleanName := "HCT116 (PXD028488)"]
group_meta[SampleGroup == "HCT116" & RowOrder == 16, CleanName := "HCT116 (PXD053474)"]

cat_ranges <- group_meta[, .(
  XMin = min(GroupIndex) - 0.45,
  XMax = max(GroupIndex) + 0.45,
  XMid = (min(GroupIndex) + max(GroupIndex)) / 2,
  CategoryLabel = unname(category_labels[first(Category)])
), by = Category]
cat_ranges[, CatOrder := match(Category, category_order)]
setorder(cat_ranges, CatOrder)
cat_dividers <- cat_ranges$XMax[-nrow(cat_ranges)] + 0.05

# =========================================================================
# 2. Figure S1a: DDR Fraction by Dataset
# =========================================================================
message(">>> Generating Figure S1a (DDR fraction by dataset with top category banner)...")

v_s1a <- merge(values_s1a, group_meta[, .(PXD, SampleGroup, GroupIndex, CleanName)], by = c("PXD", "SampleGroup"))
setorder(v_s1a, GroupIndex)
v_s1a[, DatasetXFactor := factor(GroupIndex, levels = group_meta$GroupIndex, labels = group_meta$CleanName)]
v_s1a[, Dataset := factor(Dataset, levels = c("Whole proteome", "Lactylome (Kla)"))]

s1a_stats <- v_s1a[, .(
  Mean = mean(DdrFractionPercentage),
  Median = median(DdrFractionPercentage),
  N = .N
), by = .(GroupIndex, DatasetXFactor, Dataset, Category)]

dodge_w <- 0.72
s1a_stats[, CatX := as.numeric(DatasetXFactor)]
s1a_stats[, XPos := fifelse(Dataset == "Whole proteome", CatX - dodge_w / 4, CatX + dodge_w / 4)]
y_max_s1a <- max(22, ceiling((max(v_s1a$DdrFractionPercentage) + 2) / 5) * 5)

p_s1a_track <- ggplot(cat_ranges) +
  geom_rect(
    aes(xmin = XMin, xmax = XMax, ymin = 0.10, ymax = 0.90, fill = Category),
    colour = charcoal, linewidth = 0.45, alpha = 0.85
  ) +
  geom_text(
    aes(x = XMid, y = 0.5, label = CategoryLabel),
    family = publication_font, size = 4.2, fontface = "bold", colour = "white"
  ) +
  scale_fill_manual(values = category_colours, guide = "none") +
  scale_x_continuous(limits = c(0.5, nrow(group_meta) + 0.5), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  theme_void(base_family = publication_font) +
  theme(plot.margin = margin(2, 16, 2, 12))

p_s1a_main <- ggplot(v_s1a, aes(x = DatasetXFactor, y = DdrFractionPercentage, fill = Dataset)) +
  geom_vline(xintercept = cat_dividers, colour = "#CCD1D9", linetype = "dashed", linewidth = 0.55) +
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
  scale_x_discrete(limits = levels(v_s1a$DatasetXFactor), expand = expansion(add = c(0.5, 0.5))) +
  scale_y_continuous(
    limits = c(0, y_max_s1a), breaks = seq(0, y_max_s1a, by = 5),
    labels = function(y) paste0(y, "%"), expand = expansion(mult = c(0, 0))
  ) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE, keyheight = grid::unit(0.55, "cm"), keywidth = grid::unit(0.9, "cm"))) +
  labs(
    x = NULL, y = "GO-DDR annotated protein fraction (%)", fill = NULL
  ) +
  theme_minimal(base_size = 13, base_family = publication_font) +
  theme(
    panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(colour = grid_colour, linewidth = 0.45),
    panel.border = element_blank(),
    axis.line.x = element_line(colour = "#8C939E", linewidth = 0.60),
    axis.line.y = element_line(colour = "#8C939E", linewidth = 0.60),
    axis.text.x = element_text(size = 10.0, colour = charcoal, face = "bold", angle = 50, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 11.5, colour = charcoal),
    axis.title.y = element_text(size = 14, face = "bold", colour = charcoal, margin = margin(r = 10)),
    legend.position = "top", legend.direction = "horizontal",
    legend.text = element_text(size = 12.5, colour = charcoal),
    legend.margin = margin(0, 0, 4, 0),
    plot.margin = margin(2, 16, 8, 12),
    plot.background = element_rect(fill = "white", colour = NA)
  )

fig_s1a <- (p_s1a_track / p_s1a_main) +
  plot_layout(heights = c(0.9, 12)) +
  plot_annotation(
    title = "GO-DDR annotated protein fraction across 31 individual datasets",
    theme = theme(
      plot.title = element_text(family = publication_font, size = 16, face = "bold", colour = charcoal, hjust = 0.5, margin = margin(t = 6, b = 4))
    )
  )

# =========================================================================
# 3. Figure S1b: MKI67 / H3C1 Ratio by Dataset (31 Datasets Aligned with S1a)
# =========================================================================
message(">>> Generating Figure S1b (MKI67 / H3C1 ratio across all 31 datasets, ND for undetected)...")

values_s1b_raw <- fread(file.path(candidate_dir, "figure1_mki67_ratio_sample_values.csv"), check.names = FALSE)
stop_if("ObservationType" %in% names(values_s1b_raw), "Figure S1b input lacks ObservationType.")
values_s1b <- values_s1b_raw[Denominator == "H3C1"]
if (sample_only) values_s1b <- values_s1b[ObservationType == "sample"]
expected_s1b_n <- as.integer(Sys.getenv("KLA_S1B_EXPECTED_ROWS", unset = as.character(nrow(values_s1b))))
stop_if(nrow(values_s1b) == expected_s1b_n,
  paste0("Figure S1b expects ", expected_s1b_n, " observations after the configured filter."))

v_s1b <- merge(
  group_meta[, .(PXD, SampleGroup, GroupIndex, CleanName, Category)],
  values_s1b[, .(PXD, SampleGroup, Ratio, SampleID)],
  by = c("PXD", "SampleGroup"),
  all.x = TRUE
)
setorder(v_s1b, GroupIndex)
v_s1b[, DatasetXFactor := factor(GroupIndex, levels = group_meta$GroupIndex, labels = group_meta$CleanName)]
v_s1b[, Category := factor(Category, levels = category_order)]

s1b_stats <- v_s1b[!is.na(Ratio), .(
  Mean = mean(Ratio),
  Median = median(Ratio),
  N = .N,
  MaxRatio = max(Ratio)
), by = .(GroupIndex, DatasetXFactor, Category)]

raw_max_b <- max(values_s1b$Ratio, na.rm = TRUE)
raw_min_b <- min(values_s1b$Ratio, na.rm = TRUE)
y_min_b <- 10^(floor(log10(raw_min_b)) - 0.35)
y_max_b <- 10^(log10(raw_max_b) + 0.75)
s1b_stats[, label_y := 10^(log10(MaxRatio) + 0.32)]

# Datasets with ND (Not Detected)
nd_indices <- setdiff(group_meta$GroupIndex, s1b_stats$GroupIndex)
nd_groups <- copy(group_meta[GroupIndex %in% nd_indices])
nd_groups[, DatasetXFactor := factor(GroupIndex, levels = group_meta$GroupIndex, labels = group_meta$CleanName)]
nd_groups[, Y := 10^(log10(y_min_b) + 0.28)]

p_s1b_main <- ggplot() +
  geom_vline(xintercept = cat_dividers, colour = "#CCD1D9", linetype = "dashed", linewidth = 0.55) +
  geom_boxplot(
    data = v_s1b[!is.na(Ratio)],
    aes(x = DatasetXFactor, y = Ratio, fill = Category, group = DatasetXFactor),
    width = 0.52, outlier.shape = NA, colour = charcoal,
    linewidth = 0.70, median.linewidth = 1.25, alpha = 0.85, na.rm = TRUE
  ) +
  geom_segment(
    data = s1b_stats,
    aes(x = as.numeric(DatasetXFactor) - 0.20, xend = as.numeric(DatasetXFactor) + 0.20, y = Mean, yend = Mean),
    inherit.aes = FALSE, colour = mean_colour, linewidth = 1.35
  ) +
  geom_point(
    data = v_s1b[!is.na(Ratio)],
    aes(x = DatasetXFactor, y = Ratio, fill = Category, group = DatasetXFactor),
    position = position_jitter(width = 0.12, height = 0, seed = 25),
    shape = 21, size = 2.2, stroke = 0.55, colour = "white", alpha = 0.90, na.rm = TRUE
  ) +
  geom_text(
    data = s1b_stats,
    aes(x = as.numeric(DatasetXFactor), y = label_y, label = paste0("n=", N)),
    inherit.aes = FALSE, family = publication_font, size = 3.6, fontface = "bold", colour = muted_text, vjust = -0.15
  ) +
  geom_text(
    data = nd_groups,
    aes(x = as.numeric(DatasetXFactor), y = Y, label = "ND"),
    inherit.aes = FALSE, family = publication_font, size = 2.7, angle = 90,
    fontface = "italic", colour = "#98A1AA"
  ) +
  scale_fill_manual(values = category_colours, guide = "none") +
  scale_x_discrete(limits = levels(v_s1b$DatasetXFactor), expand = expansion(add = c(0.5, 0.5))) +
  scale_y_log10(
    breaks = 10^seq(floor(log10(y_min_b)), ceiling(log10(y_max_b)), by = 1),
    labels = scales::label_scientific(digits = 2),
    limits = c(y_min_b, y_max_b),
    expand = c(0, 0)
  ) +
  labs(
    x = NULL, y = "MKI67 / H3C1 ratio (log scale)"
  ) +
  theme_minimal(base_size = 13, base_family = publication_font) +
  theme(
    panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(colour = grid_colour, linewidth = 0.45),
    panel.border = element_blank(),
    axis.line.x = element_line(colour = "#8C939E", linewidth = 0.60),
    axis.line.y = element_line(colour = "#8C939E", linewidth = 0.60),
    axis.text.x = element_text(size = 10.0, colour = charcoal, face = "bold", angle = 50, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 11.5, colour = charcoal),
    axis.title.y = element_text(size = 14, face = "bold", colour = charcoal, margin = margin(r = 10)),
    plot.margin = margin(2, 16, 8, 12),
    plot.background = element_rect(fill = "white", colour = NA)
  )

p_s1b_track <- ggplot(cat_ranges) +
  geom_rect(
    aes(xmin = XMin, xmax = XMax, ymin = 0.10, ymax = 0.90, fill = Category),
    colour = charcoal, linewidth = 0.45, alpha = 0.85
  ) +
  geom_text(
    aes(x = XMid, y = 0.5, label = CategoryLabel),
    family = publication_font, size = 4.2, fontface = "bold", colour = "white"
  ) +
  scale_fill_manual(values = category_colours, guide = "none") +
  scale_x_continuous(limits = c(0.5, nrow(group_meta) + 0.5), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  theme_void(base_family = publication_font) +
  theme(plot.margin = margin(2, 16, 2, 12))

fig_s1b <- (p_s1b_track / p_s1b_main) +
  plot_layout(heights = c(0.9, 12)) +
  plot_annotation(
    title = "MKI67 / H3C1 protein intensity ratio across individual whole-proteome datasets (ND = Not Detected)",
    theme = theme(
      plot.title = element_text(family = publication_font, size = 16, face = "bold", colour = charcoal, hjust = 0.5, margin = margin(t = 6, b = 4))
    )
  )

# =========================================================================
# 4. Companion Output: Figure S1b with Only 11 Detected Datasets
# =========================================================================
detected_dataset_n <- uniqueN(values_s1b[, .(PXD, SampleGroup)])
message(">>> Generating companion Figure S1b (only ", detected_dataset_n, " sample-level detected datasets)...")

detected_meta <- copy(group_meta[GroupIndex %in% s1b_stats$GroupIndex])
setorder(detected_meta, CategoryOrder, RowOrder)
detected_meta[, DetectedIndex := .I]

v_s1b_det <- merge(values_s1b, detected_meta[, .(PXD, SampleGroup, DetectedIndex, CleanName)], by = c("PXD", "SampleGroup"))
setorder(v_s1b_det, DetectedIndex)
v_s1b_det[, DatasetXFactor := factor(DetectedIndex, levels = detected_meta$DetectedIndex, labels = detected_meta$CleanName)]
v_s1b_det[, Category := factor(Category, levels = category_order)]

s1b_det_stats <- v_s1b_det[, .(
  Mean = mean(Ratio),
  Median = median(Ratio),
  N = .N,
  MaxRatio = max(Ratio)
), by = .(DetectedIndex, DatasetXFactor, Category)]
s1b_det_stats[, label_y := 10^(log10(MaxRatio) + 0.32)]

cat_ranges_det <- detected_meta[, .(
  XMin = min(DetectedIndex) - 0.45,
  XMax = max(DetectedIndex) + 0.45,
  XMid = (min(DetectedIndex) + max(DetectedIndex)) / 2,
  CategoryLabel = unname(category_labels[first(Category)])
), by = Category]
cat_ranges_det[, CatOrder := match(Category, category_order)]
cat_ranges_det[, LabelSize := fifelse((XMax - XMin) < 2.2, 3.1, 4.2)]
cat_ranges_det[, CategoryLabel := fifelse(
  (XMax - XMin) < 2.2,
  sub(" ", "\n", CategoryLabel, fixed = TRUE),
  CategoryLabel
)]
setorder(cat_ranges_det, CatOrder)
cat_dividers_det <- cat_ranges_det$XMax[-nrow(cat_ranges_det)] + 0.05

p_s1b_det_track <- ggplot(cat_ranges_det) +
  geom_rect(
    aes(xmin = XMin, xmax = XMax, ymin = 0.10, ymax = 0.90, fill = Category),
    colour = charcoal, linewidth = 0.45, alpha = 0.85
  ) +
  geom_text(
    aes(x = XMid, y = 0.5, label = CategoryLabel, size = LabelSize),
    family = publication_font, fontface = "bold", colour = "white", show.legend = FALSE
  ) +
  scale_size_identity() +
  scale_fill_manual(values = category_colours, guide = "none") +
  scale_x_continuous(limits = c(0.5, nrow(detected_meta) + 0.5), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  theme_void(base_family = publication_font) +
  theme(plot.margin = margin(2, 16, 2, 12))

p_s1b_det_main <- ggplot(v_s1b_det, aes(x = DatasetXFactor, y = Ratio, fill = Category)) +
  geom_vline(xintercept = cat_dividers_det, colour = "#CCD1D9", linetype = "dashed", linewidth = 0.55) +
  geom_boxplot(
    aes(group = DatasetXFactor),
    width = 0.52, outlier.shape = NA, colour = charcoal,
    linewidth = 0.70, median.linewidth = 1.25, alpha = 0.85, na.rm = TRUE
  ) +
  geom_segment(
    data = s1b_det_stats,
    aes(x = as.numeric(DatasetXFactor) - 0.20, xend = as.numeric(DatasetXFactor) + 0.20, y = Mean, yend = Mean),
    inherit.aes = FALSE, colour = mean_colour, linewidth = 1.35
  ) +
  geom_point(
    aes(group = DatasetXFactor),
    position = position_jitter(width = 0.12, height = 0, seed = 25),
    shape = 21, size = 2.4, stroke = 0.55, colour = "white", alpha = 0.90, na.rm = TRUE
  ) +
  geom_text(
    data = s1b_det_stats,
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
    x = NULL, y = "MKI67 / H3C1 ratio (log scale)"
  ) +
  theme_minimal(base_size = 13, base_family = publication_font) +
  theme(
    panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(colour = grid_colour, linewidth = 0.45),
    panel.border = element_blank(),
    axis.line.x = element_line(colour = "#8C939E", linewidth = 0.60),
    axis.line.y = element_line(colour = "#8C939E", linewidth = 0.60),
    axis.text.x = element_text(size = 11.0, colour = charcoal, face = "bold", angle = 50, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 11.5, colour = charcoal),
    axis.title.y = element_text(size = 14, face = "bold", colour = charcoal, margin = margin(r = 10)),
    plot.margin = margin(2, 16, 8, 12),
    plot.background = element_rect(fill = "white", colour = NA)
  )

fig_s1b_det <- (p_s1b_det_track / p_s1b_det_main) +
  plot_layout(heights = c(0.9, 12)) +
  plot_annotation(
    title = paste0("MKI67 / H3C1 protein intensity ratio across ", detected_dataset_n, " quantified whole-proteome datasets"),
    theme = theme(
      plot.title = element_text(family = publication_font, size = 16, face = "bold", colour = charcoal, hjust = 0.5, margin = margin(t = 6, b = 4))
    )
  )

# =========================================================================
# 5. Save Both Figures with Exact Same Proportions (Consistent Ratio)
# =========================================================================
fig_w <- 15.0
fig_h <- 7.5

s1a_stem <- file.path(output_dir, "Figure_S1a_DDR_fraction_by_PXD")
ggsave(paste0(s1a_stem, ".png"), fig_s1a, width = fig_w, height = fig_h, dpi = 300, bg = "white", device = ragg::agg_png)
ggsave(paste0(s1a_stem, ".pdf"), fig_s1a, width = fig_w, height = fig_h, bg = "white", device = cairo_pdf)
message("Wrote Figure S1a to ", s1a_stem, ".png/.pdf")

# Primary S1b: 31 datasets aligned with S1a (ND for missing)
s1b_stem <- file.path(output_dir, "Figure_S1b_MKI67_over_H3C1_by_PXD")
ggsave(paste0(s1b_stem, ".png"), fig_s1b, width = fig_w, height = fig_h, dpi = 300, bg = "white", device = ragg::agg_png)
ggsave(paste0(s1b_stem, ".pdf"), fig_s1b, width = fig_w, height = fig_h, bg = "white", device = cairo_pdf)
message("Wrote Figure S1b (31 datasets aligned) to ", s1b_stem, ".png/.pdf")

# Companion S1b: 11 detected datasets only
s1b_det_stem <- file.path(output_dir, "Figure_S1b_MKI67_over_H3C1_11_detected_only")
ggsave(paste0(s1b_det_stem, ".png"), fig_s1b_det, width = fig_w, height = fig_h, dpi = 300, bg = "white", device = ragg::agg_png)
ggsave(paste0(s1b_det_stem, ".pdf"), fig_s1b_det, width = fig_w, height = fig_h, bg = "white", device = cairo_pdf)
message("Wrote companion Figure S1b (11 detected only) to ", s1b_det_stem, ".png/.pdf")
