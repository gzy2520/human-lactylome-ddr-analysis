#!/usr/bin/env Rscript

# Candidate Figure 1 renderer for dataset-level four-category boxplots.
#
# Each point is one PXD/sample-group union. The four categories are the x
# axis, and Whole proteome/Kla are the two adjacent boxes within each category.
# The approved publication renderer and its frozen inputs are not changed.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

set.seed(25)

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1L) {
  normalizePath(args[[1L]], mustWork = TRUE)
} else {
  normalizePath(".", mustWork = TRUE)
}

stop_if <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

publication_font <- "Arial Unicode MS"
whole_proteome_colour <- "#4E79A7"
kla_colour <- "#F28E2B"
mean_colour <- "#C0392B"
charcoal <- "#2F3437"
muted_text <- "#65717D"
grid_colour <- "#D9DDE3"
panel_border_colour <- "#C8CED6"

category_order <- c("normal_tissue", "cancer_tissue", "cancer_cells", "normal_cells")
category_labels <- c(
  normal_tissue = "non-tumor tissues",
  cancer_tissue = "tumor tissues",
  cancer_cells = "cancer cell lines",
  normal_cells = "normal cell lines"
)
dataset_order <- c("Whole proteome", "Lactylome (Kla)")

candidate_dir <- normalizePath(
  Sys.getenv("KLA_CANDIDATE_INPUT", unset = file.path(project_root, "data", "candidate")),
  mustWork = TRUE
)
output_dir <- normalizePath(
  Sys.getenv("KLA_CANDIDATE_OUTPUT", unset = file.path(project_root, "results", "candidate")),
  mustWork = FALSE
)
values_path <- file.path(candidate_dir, "figure1_dataset_boxplot_values.csv")
significance_path <- file.path(project_root, "R", "candidate", "boxplot_significance.R")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

stop_if(file.exists(values_path), paste0("Missing dataset-level Figure 1 input: ", values_path))
stop_if(file.exists(significance_path), paste0("Missing boxplot significance helper: ", significance_path))
values <- fread(values_path, check.names = FALSE)
source(significance_path, local = TRUE)

required_values <- c(
  "RowOrder", "PXD", "SampleGroup", "Category", "Dataset", "DatasetPXD",
  "DatasetPointID", "DatasetPointLabel", "SampleID", "SourceFile",
  "DdrFraction", "DdrFractionPercentage", "ProteinCount", "DdrProteinCount"
)
stop_if(all(required_values %in% names(values)),
  "Dataset-level Figure 1 input schema is incomplete.")
stop_if(nrow(values) > 0L, "Dataset-level Figure 1 input is empty.")
stop_if(!anyDuplicated(values$DatasetPointID),
  "Dataset-level Figure 1 input contains duplicate point IDs.")
stop_if(all(values$Dataset %in% dataset_order),
  "Dataset-level Figure 1 input contains an unknown dataset.")
stop_if(all(values$Category %in% category_order),
  "Dataset-level Figure 1 input contains an unknown category.")
stop_if(all(is.finite(values$DdrFraction) & is.finite(values$DdrFractionPercentage)),
  "Dataset-level Figure 1 fractions contain non-finite values.")
stop_if(all(values$DdrFraction >= 0 & values$DdrFraction <= 1),
  "Dataset-level Figure 1 fractions must be between 0 and 1.")
stop_if(all(values$DdrFractionPercentage >= 0 & values$DdrFractionPercentage <= 100),
  "Dataset-level Figure 1 percentages must be between 0 and 100.")

actual_counts <- values[, .(N = .N), by = .(Category, Dataset)]
expected_groups <- unique(values[, .(PXD, SampleGroup, Category)])
expected_category_counts <- expected_groups[, .(ExpectedN = .N), by = Category]
expected_counts <- expected_category_counts[
  rep(seq_len(nrow(expected_category_counts)), each = length(dataset_order))
]
expected_counts[, Dataset := rep(dataset_order, times = nrow(expected_category_counts))]
count_check <- merge(expected_counts, actual_counts, by = c("Category", "Dataset"), all = TRUE)
stop_if(nrow(count_check) == length(category_order) * length(dataset_order),
  "Dataset-level Figure 1 input does not contain eight category-dataset groups.")
stop_if(all(count_check$ExpectedN == count_check$N),
  "Dataset-level Figure 1 category-dataset counts are inconsistent.")

values[, CategoryIndex := match(Category, category_order)]
values[, CategoryLabel := factor(
  Category,
  levels = category_order,
  labels = unname(category_labels[category_order])
)]
values[, Dataset := factor(Dataset, levels = dataset_order)]

stop_if(identical(levels(values$CategoryLabel), unname(category_labels[category_order])),
  "Figure 1 category factor was not constructed.")
stop_if(identical(levels(values$Dataset), dataset_order),
  "Figure 1 dataset factor was not constructed.")

figure1_anova <- compute_category_one_way_anova(
  values,
  category_order = category_order,
  dataset_order = dataset_order,
  value_column = "DdrFractionPercentage"
)
stop_if(nrow(figure1_anova) == length(dataset_order),
  "Figure 1 one-way ANOVA results are incomplete.")
stop_if(all(is.finite(figure1_anova$PValue) & is.finite(figure1_anova$QValueBH)),
  "Figure 1 one-way ANOVA results contain a non-finite p or q value.")
fwrite(figure1_anova, file.path(output_dir, "figure1_category_one_way_anova.csv"), na = "")

format_p <- function(value) {
  if (!is.finite(value)) return("NA")
  if (value < 0.001) formatC(value, format = "e", digits = 2)
  else formatC(value, format = "f", digits = 3)
}
anova_label <- paste(vapply(seq_len(nrow(figure1_anova)), function(index) {
  row <- figure1_anova[index]
  paste0(row$Dataset, ": p=", format_p(row$PValue),
    ", q=", format_p(row$QValueBH), " ", row$Significance)
}, character(1)), collapse = "  |  ")

box_stats <- values[, .(
  Mean = mean(DdrFractionPercentage),
  Median = median(DdrFractionPercentage),
  N = .N
), by = .(Category, Dataset, CategoryIndex)]
box_stats[, Dataset := factor(Dataset, levels = dataset_order)]
box_stats[, x_center := CategoryIndex + fifelse(
  as.character(Dataset) == "Whole proteome", -0.19, 0.19
)]
box_stats[, c("x_left", "x_right") := list(x_center - 0.21, x_center + 0.21)]

y_limit <- max(14, ceiling((max(values$DdrFractionPercentage) + 1) / 2) * 2)
y_breaks <- seq(0, y_limit, by = 2)
if (tail(y_breaks, 1L) < y_limit) y_breaks <- c(y_breaks, y_limit)

figure_plot <- ggplot(
  values,
  aes(x = CategoryIndex, y = DdrFractionPercentage, fill = Dataset)
) +
  geom_boxplot(
    aes(group = Dataset),
    position = position_dodge(width = 0.76),
    width = 0.62,
    outlier.shape = NA,
    colour = charcoal,
    linewidth = 0.78,
    alpha = 0.84,
    orientation = "x",
    na.rm = TRUE
  ) +
  geom_point(
    aes(group = Dataset),
    position = position_jitterdodge(
      jitter.width = 0.045,
      jitter.height = 0,
      dodge.width = 0.76,
      seed = 25
    ),
    shape = 21,
    size = 3.15,
    stroke = 0.58,
    colour = "white",
    alpha = 0.94,
    na.rm = TRUE
  ) +
  geom_segment(
    data = box_stats,
    aes(x = x_left, xend = x_right, y = Median, yend = Median),
    inherit.aes = FALSE,
    colour = charcoal,
    linewidth = 1.15,
    lineend = "round"
  ) +
  geom_segment(
    data = box_stats,
    aes(x = x_left, xend = x_right, y = Mean, yend = Mean),
    inherit.aes = FALSE,
    colour = mean_colour,
    linewidth = 1.55,
    lineend = "round"
  ) +
  scale_fill_manual(
    values = c("Whole proteome" = whole_proteome_colour, "Lactylome (Kla)" = kla_colour),
    breaks = dataset_order
  ) +
  scale_x_continuous(
    breaks = seq_along(category_order),
    labels = unname(category_labels[category_order]),
    limits = c(0.45, length(category_order) + 0.55),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(
    limits = c(0, y_limit),
    breaks = y_breaks,
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0, 0.02))
  ) +
  guides(
    fill = guide_legend(
      nrow = 1,
      byrow = TRUE,
      keyheight = grid::unit(0.62, "cm"),
      keywidth = grid::unit(0.92, "cm")
    )
  ) +
  labs(
    title = "Dataset-level DDR fraction across four biological categories",
    subtitle = paste0("One-way ANOVA by category; one point = one PXD/sample-group union.  ", anova_label),
    x = NULL,
    y = "GO-DDR annotated protein fraction (%)",
    fill = NULL,
    caption = paste(
      "The dark horizontal line inside each box is the median; the red horizontal line is the mean.",
      "Whole-proteome points use the corresponding reference PXD, while Kla points use the source PXD.",
      "ANOVA q values are BH-adjusted across the two modality-specific omnibus tests (**** q<0.0001, *** q<0.001, ** q<0.01, * q<0.05).",
      sep = "\n"
    )
  ) +
  theme_minimal(base_size = 14, base_family = publication_font) +
  theme(
    panel.grid.major.y = element_line(colour = grid_colour, linewidth = 0.50),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.border = element_rect(colour = panel_border_colour, fill = NA, linewidth = 0.60),
    axis.text.x = element_text(size = 14.5, colour = charcoal, face = "bold", lineheight = 0.96),
    axis.text.y = element_text(size = 13.5, colour = charcoal),
    axis.title.y = element_text(size = 17.0, face = "bold", colour = charcoal, margin = margin(r = 12)),
    plot.title = element_text(size = 20.0, face = "bold", colour = charcoal, hjust = 0.5, margin = margin(b = 6)),
    plot.subtitle = element_text(size = 11.5, colour = muted_text, hjust = 0.5, margin = margin(b = 11)),
    plot.caption = element_text(size = 10.5, colour = muted_text, hjust = 0, lineheight = 1.05, margin = margin(t = 10)),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.text = element_text(size = 15.0, colour = charcoal),
    legend.key.spacing.x = grid::unit(0.38, "cm"),
    legend.background = element_rect(fill = "white", colour = NA),
    legend.margin = margin(1, 0, 10, 0),
    plot.margin = margin(12, 18, 14, 16),
    plot.background = element_rect(fill = "white", colour = NA)
  )

stop_if(length(ggplotGrob(figure_plot)$grobs) > 0L,
  "Figure 1 dataset-level boxplot grob is empty.")

stem <- file.path(output_dir, "Figure_1_DDR_fraction_candidate_category_boxplot_refined")
ggsave(
  paste0(stem, ".png"), figure_plot,
  width = 14.0, height = 8.8, dpi = 300, bg = "white", device = ragg::agg_png
)
ggsave(
  paste0(stem, ".pdf"), figure_plot,
  width = 14.0, height = 8.8, bg = "white", device = cairo_pdf
)

message(
  "Wrote dataset-level four-category Figure 1 boxplot: ", stem,
  ".png/.pdf and one-way ANOVA statistics (", nrow(expected_groups), " PXD/sample-group points per modality)."
)
