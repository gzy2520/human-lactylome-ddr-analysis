#!/usr/bin/env Rscript

# Re-render the approved-review Figure 1 layout: four biological-category
# panels, whole proteome and Kla as the two horizontal boxplots, and each dot
# as a deposited source-sample observation.  This script deliberately keeps
# the layout used in the reviewed 2026-09-02 candidate figure.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
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
panel_border_colour <- "#C8CED6"

category_order <- c("normal_tissue", "cancer_tissue", "normal_cells", "cancer_cells")
category_labels <- c(
  normal_tissue = "non-tumor tissues",
  cancer_tissue = "tumor tissues",
  normal_cells = "normal cell lines",
  cancer_cells = "cancer cell lines"
)
category_fills <- c(
  normal_tissue = "#DCE9E2",
  cancer_tissue = "#F0DEDE",
  normal_cells = "#E7E1EE",
  cancer_cells = "#EEE4D2"
)
dataset_order <- c("Whole proteome", "Lactylome (Kla)")
plot_dataset_order <- rev(dataset_order)

candidate_dir <- normalizePath(Sys.getenv(
  "KLA_CANDIDATE_INPUT", unset = file.path(project_root, "data", "candidate")
), mustWork = TRUE)
output_dir <- Sys.getenv("KLA_CANDIDATE_OUTPUT", unset = file.path(project_root, "results", "candidate"))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
values_path <- file.path(candidate_dir, "figure1_sample_boxplot_values.csv")
stop_if(file.exists(values_path), paste0("Missing Figure 1 sample input: ", values_path))
values <- fread(values_path, check.names = FALSE)

required_values <- c("PXD", "SampleGroup", "Category", "Dataset", "SampleID", "DdrFractionPercentage")
stop_if(all(required_values %in% names(values)), "Figure 1 sample input schema is incomplete.")
stop_if(all(values$Dataset %in% dataset_order), "Figure 1 contains an unknown modality.")
stop_if(all(values$Category %in% category_order), "Figure 1 contains an unknown category.")
stop_if(all(is.finite(values$DdrFractionPercentage) & values$DdrFractionPercentage >= 0 & values$DdrFractionPercentage <= 100),
  "Figure 1 fractions must be finite percentages from 0 to 100.")
stop_if(!anyDuplicated(values[, .(PXD, SampleGroup, Dataset, SampleID)]),
  "Figure 1 source-sample points must be unique within modality.")

source(file.path(project_root, "R", "candidate", "boxplot_significance.R"), local = TRUE)
anova <- compute_figure1_sample_one_way_anova(values, category_order, dataset_order)
stop_if(nrow(anova) == length(category_order) && all(is.finite(anova$PValue)),
  "Figure 1 one-way ANOVA did not produce four finite category tests.")
fwrite(anova, file.path(output_dir, "figure1_category_one_way_anova.csv"), na = "")

values[, CategoryLabel := factor(Category, levels = category_order, labels = unname(category_labels[category_order]))]
values[, Dataset := factor(Dataset, levels = plot_dataset_order)]
panel_counts <- values[, .(N = .N), by = .(CategoryLabel, Dataset)]
summary_stats <- values[, .(
  Mean = mean(DdrFractionPercentage),
  Median = median(DdrFractionPercentage)
), by = .(Category, CategoryLabel, Dataset)]
fwrite(summary_stats, file.path(output_dir, "figure1_category_boxplot_mean_median.csv"), na = "")

x_limit <- max(20, ceiling(max(values$DdrFractionPercentage) * 1.35 / 5) * 5)
annotation <- copy(anova)
annotation[, CategoryLabel := factor(Category, levels = category_order, labels = unname(category_labels[category_order]))]
annotation[, `:=`(
  x_left = x_limit * 0.805,
  x_right = x_limit * 0.825,
  x_text = x_limit * 0.842,
  y_low = 1,
  y_high = 2,
  y_mid = 1.5,
  Label = Significance
)]

figure_plot <- ggplot(values, aes(x = DdrFractionPercentage, y = Dataset, fill = Dataset)) +
  geom_boxplot(
    aes(group = Dataset), width = 0.60, outlier.shape = NA, colour = charcoal,
    linewidth = 0.82, median.linewidth = 1.35, alpha = 0.82, na.rm = TRUE
  ) +
  geom_point(
    aes(group = Dataset), position = position_jitter(width = 0, height = 0.12, seed = 25),
    shape = 21, size = 3.25, stroke = 0.62, colour = "white", alpha = 0.94, na.rm = TRUE
  ) +
  geom_point(
    data = summary_stats, aes(x = Mean, y = Dataset), inherit.aes = FALSE,
    shape = 124, size = 8.8, stroke = 1.15, colour = mean_colour
  ) +
  geom_text(
    data = panel_counts, aes(x = x_limit * 0.925, y = Dataset, label = paste0("n=", N)),
    inherit.aes = FALSE, hjust = 0, size = 4.30, family = publication_font, colour = muted_text
  ) +
  geom_segment(
    data = annotation, aes(x = x_left, xend = x_right, y = y_low, yend = y_low),
    inherit.aes = FALSE, colour = charcoal, linewidth = 0.78
  ) +
  geom_segment(
    data = annotation, aes(x = x_right, xend = x_right, y = y_low, yend = y_high),
    inherit.aes = FALSE, colour = charcoal, linewidth = 0.78
  ) +
  geom_segment(
    data = annotation, aes(x = x_left, xend = x_right, y = y_high, yend = y_high),
    inherit.aes = FALSE, colour = charcoal, linewidth = 0.78
  ) +
  geom_text(
    data = annotation, aes(x = x_text, y = y_mid, label = Label), inherit.aes = FALSE,
    hjust = 0, size = 5.0, family = publication_font, colour = charcoal
  ) +
  facet_grid(CategoryLabel ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_manual(values = c("Whole proteome" = whole_proteome_colour, "Lactylome (Kla)" = kla_colour), breaks = dataset_order) +
  scale_x_continuous(limits = c(0, x_limit), breaks = scales::pretty_breaks(n = 5), expand = expansion(mult = c(0, 0))) +
  scale_y_discrete(labels = c("Lactylome (Kla)" = "Kla", "Whole proteome" = "Whole proteome")) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE, keyheight = grid::unit(0.62, "cm"), keywidth = grid::unit(0.92, "cm"))) +
  labs(
    x = "GO-DDR annotated protein fraction (%)", y = NULL, fill = NULL,
    caption = paste(
      "Each point is one source-resolved sample observation. Dark box line = median; red vertical line = mean.",
      "Stars show BH-adjusted one-way ANOVA tests between Whole proteome and Kla within each category",
      "(**** q<0.0001, *** q<0.001, ** q<0.01, * q<0.05).",
      sep = "\n"
    )
  ) +
  theme_minimal(base_size = 14, base_family = publication_font) +
  theme(
    panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(colour = grid_colour, linewidth = 0.50),
    panel.border = element_rect(colour = panel_border_colour, fill = NA, linewidth = 0.60),
    axis.text.y = element_text(size = 15.5, colour = charcoal, face = "bold"),
    axis.text.x = element_text(size = 15, colour = charcoal),
    axis.title.x = element_text(size = 19, face = "bold", colour = charcoal, margin = margin(t = 14)),
    axis.title.y = element_blank(),
    strip.placement = "outside",
    strip.text.y.left = element_text(size = 18, face = "bold", colour = charcoal, angle = 90, lineheight = 0.95),
    strip.background = element_rect(fill = "#E7E9E7", colour = NA),
    panel.spacing.y = grid::unit(0.72, "lines"),
    legend.position = "top", legend.direction = "horizontal",
    legend.text = element_text(size = 16.5, colour = charcoal, lineheight = 1.10),
    legend.key.spacing.x = grid::unit(0.38, "cm"), legend.background = element_rect(fill = "white", colour = NA),
    legend.margin = margin(1, 0, 10, 0),
    plot.caption = element_text(size = 10.4, hjust = 0.5, colour = charcoal, margin = margin(t = 10)),
    plot.margin = margin(10, 18, 14, 12), plot.background = element_rect(fill = "white", colour = NA)
  )

apply_strip_fills <- function(plot) {
  plot_grob <- ggplotGrob(plot)
  strip_ids <- grep("^strip-l", plot_grob$layout$name)
  strip_ids <- strip_ids[order(plot_grob$layout$t[strip_ids])]
  stop_if(length(strip_ids) == length(category_order), "Figure 1 must contain four category strips.")
  for (index in seq_along(strip_ids)) {
    strip_grob <- plot_grob$grobs[[strip_ids[[index]]]]$grobs[[1L]]
    background_id <- grep("^strip.background", strip_grob$childrenOrder)
    stop_if(length(background_id) == 1L, "Unable to identify a Figure 1 category-strip background.")
    strip_grob$children[[background_id]]$gp$fill <- unname(category_fills[category_order[[index]]])
    strip_grob$children[[background_id]]$gp$col <- NA
    plot_grob$grobs[[strip_ids[[index]]]]$grobs[[1L]] <- strip_grob
  }
  plot_grob
}

plot_grob <- apply_strip_fills(figure_plot)
stem <- file.path(output_dir, "Figure_1_DDR_fraction_candidate_category_boxplot_refined")
ggsave(paste0(stem, ".png"), plot_grob, width = 15.5, height = 11.5, dpi = 300, bg = "white", device = ragg::agg_png)
ggsave(paste0(stem, ".pdf"), plot_grob, width = 15.5, height = 11.5, bg = "white", device = cairo_pdf)
message("Wrote restored-layout Figure 1 source-sample boxplot: ", stem, ".png/.pdf")
