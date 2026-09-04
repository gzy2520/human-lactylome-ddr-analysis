#!/usr/bin/env Rscript

# Render the upright Figure 1 candidate layout: four biological-category
# panels (columns), whole proteome and Kla as two vertical boxplots in each panel,
# and each dot as an individual deposited source-sample observation.

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
values[, Dataset := factor(Dataset, levels = dataset_order)]
panel_counts <- values[, .(N = .N), by = .(CategoryLabel, Dataset)]
summary_stats <- values[, .(
  Mean = mean(DdrFractionPercentage),
  Median = median(DdrFractionPercentage)
), by = .(Category, CategoryLabel, Dataset)]
fwrite(summary_stats, file.path(output_dir, "figure1_category_boxplot_mean_median.csv"), na = "")

max_val <- max(values$DdrFractionPercentage, na.rm = TRUE)
y_limit <- max(20, ceiling(max(max_val * 1.32, max_val + 4) / 5) * 5)

annotation <- copy(anova)
annotation[, CategoryLabel := factor(Category, levels = category_order, labels = unname(category_labels[category_order]))]
annotation[, `:=`(
  x_left = 1,
  x_right = 2,
  x_mid = 1.5,
  y_bar = y_limit * 0.92,
  y_tick = y_limit * 0.895,
  y_text = y_limit * 0.935,
  Label = Significance
)]

figure_plot <- ggplot(values, aes(x = Dataset, y = DdrFractionPercentage, fill = Dataset)) +
  geom_boxplot(
    aes(group = Dataset), width = 0.58, outlier.shape = NA, colour = charcoal,
    linewidth = 0.82, median.linewidth = 1.35, alpha = 0.82, na.rm = TRUE
  ) +
  geom_segment(
    data = summary_stats,
    aes(x = as.numeric(Dataset) - 0.22, xend = as.numeric(Dataset) + 0.22, y = Mean, yend = Mean),
    inherit.aes = FALSE, colour = mean_colour, linewidth = 1.55
  ) +
  geom_point(
    aes(group = Dataset), position = position_jitter(width = 0.16, height = 0, seed = 25),
    shape = 21, size = 3.0, stroke = 0.58, colour = "white", alpha = 0.88, na.rm = TRUE
  ) +
  geom_text(
    data = panel_counts, aes(x = Dataset, y = y_limit * 0.82, label = paste0("n=", N)),
    inherit.aes = FALSE, size = 4.4, family = publication_font, colour = muted_text, fontface = "bold"
  ) +
  geom_segment(
    data = annotation, aes(x = x_left, xend = x_right, y = y_bar, yend = y_bar),
    inherit.aes = FALSE, colour = charcoal, linewidth = 0.78
  ) +
  geom_segment(
    data = annotation, aes(x = x_left, xend = x_left, y = y_tick, yend = y_bar),
    inherit.aes = FALSE, colour = charcoal, linewidth = 0.78
  ) +
  geom_segment(
    data = annotation, aes(x = x_right, xend = x_right, y = y_tick, yend = y_bar),
    inherit.aes = FALSE, colour = charcoal, linewidth = 0.78
  ) +
  geom_text(
    data = annotation, aes(x = x_mid, y = y_text, label = Label), inherit.aes = FALSE,
    size = 5.6, family = publication_font, colour = charcoal, fontface = "bold"
  ) +
  facet_grid(. ~ CategoryLabel) +
  scale_fill_manual(values = c("Whole proteome" = whole_proteome_colour, "Lactylome (Kla)" = kla_colour), breaks = dataset_order) +
  scale_x_discrete(labels = c("Whole proteome" = "Whole\nproteome", "Lactylome (Kla)" = "Lactylome\n(Kla)")) +
  scale_y_continuous(limits = c(0, y_limit), breaks = scales::pretty_breaks(n = 5), labels = function(y) paste0(y, "%"), expand = expansion(mult = c(0, 0))) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE, keyheight = grid::unit(0.62, "cm"), keywidth = grid::unit(0.92, "cm"))) +
  labs(
    x = NULL, y = "GO-DDR annotated protein fraction (%)", fill = NULL,
    caption = paste(
      "Each point is one source-resolved sample observation. Dark box line = median; red horizontal line = mean.",
      "Stars show BH-adjusted one-way ANOVA tests between Whole proteome and Kla within each category",
      "(**** q<0.0001, *** q<0.001, ** q<0.01, * q<0.05).",
      sep = "\n"
    )
  ) +
  theme_minimal(base_size = 14, base_family = publication_font) +
  theme(
    panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(colour = grid_colour, linewidth = 0.50),
    panel.border = element_rect(colour = panel_border_colour, fill = NA, linewidth = 0.60),
    axis.text.x = element_text(size = 14.5, colour = charcoal, face = "bold", lineheight = 0.95),
    axis.text.y = element_text(size = 14.5, colour = charcoal),
    axis.title.y = element_text(size = 18, face = "bold", colour = charcoal, margin = margin(r = 12)),
    strip.placement = "outside",
    strip.text.x.top = element_text(size = 16, face = "bold", colour = charcoal, margin = margin(t = 6, b = 6)),
    strip.background = element_rect(fill = "#E7E9E7", colour = NA),
    panel.spacing.x = grid::unit(0.9, "lines"),
    legend.position = "top", legend.direction = "horizontal",
    legend.text = element_text(size = 15.5, colour = charcoal),
    legend.key.spacing.x = grid::unit(0.38, "cm"), legend.background = element_rect(fill = "white", colour = NA),
    legend.margin = margin(1, 0, 8, 0),
    plot.caption = element_text(size = 10.4, hjust = 0.5, colour = charcoal, margin = margin(t = 12)),
    plot.margin = margin(10, 18, 14, 14), plot.background = element_rect(fill = "white", colour = NA)
  )

apply_strip_fills <- function(plot) {
  plot_grob <- ggplotGrob(plot)
  strip_ids <- grep("^strip-t", plot_grob$layout$name)
  if (length(strip_ids) == 0L) {
    strip_ids <- grep("^strip-l", plot_grob$layout$name)
    strip_ids <- strip_ids[order(plot_grob$layout$t[strip_ids])]
  } else {
    strip_ids <- strip_ids[order(plot_grob$layout$l[strip_ids])]
  }
  for (index in seq_along(strip_ids)) {
    strip_grob <- plot_grob$grobs[[strip_ids[[index]]]]$grobs[[1L]]
    background_id <- grep("^strip.background", strip_grob$childrenOrder)
    strip_grob$children[[background_id]]$gp$fill <- unname(category_fills[category_order[[index]]])
    strip_grob$children[[background_id]]$gp$col <- NA
    plot_grob$grobs[[strip_ids[[index]]]]$grobs[[1L]] <- strip_grob
  }
  plot_grob
}

plot_grob <- apply_strip_fills(figure_plot)
stem <- file.path(output_dir, "Figure_1_DDR_fraction_candidate_category_boxplot_refined")
ggsave(paste0(stem, ".png"), plot_grob, width = 14, height = 9, dpi = 300, bg = "white", device = ragg::agg_png)
ggsave(paste0(stem, ".pdf"), plot_grob, width = 14, height = 9, bg = "white", device = cairo_pdf)
message("Wrote upright Figure 1 source-sample boxplot: ", stem, ".png/.pdf")
