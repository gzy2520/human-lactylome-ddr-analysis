#!/usr/bin/env Rscript

# Candidate Figure 1 renderer (refined visual version).
#
# The four publication categories are the plotting units.  Each category has
# two boxplots (whole proteome and Kla), and every point is one observation
# from the source-resolved input table.  This is an isolated review output;
# the approved publication renderer and its figures are not changed.

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

candidate_dir <- file.path(project_root, "data", "candidate")
output_dir <- file.path(project_root, "results", "candidate")
values_path <- file.path(candidate_dir, "figure1_sample_boxplot_values.csv")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

stop_if(file.exists(values_path), paste0("Missing Figure 1 sample input: ", values_path))
values <- fread(values_path, check.names = FALSE)

required_values <- c(
  "PXD", "SampleGroup", "Category", "Dataset", "SampleID",
  "DdrFractionPercentage"
)
stop_if(all(required_values %in% names(values)), "Figure 1 sample input schema is incomplete.")
stop_if(nrow(values) == 210L, "The category-level Figure 1 input must contain 210 observations.")
stop_if(all(values$Dataset %in% dataset_order),
  "Figure 1 sample input contains an unknown dataset.")
stop_if(all(values$Category %in% category_order),
  "Figure 1 sample input contains an unknown category.")
stop_if(all(is.finite(values$DdrFractionPercentage)),
  "Figure 1 sample fractions contain non-finite values.")
stop_if(all(values$DdrFractionPercentage >= 0 & values$DdrFractionPercentage <= 100),
  "Figure 1 sample fractions must be between 0 and 100 percent.")

expected_counts <- data.table(
  Category = rep(category_order, each = 2L),
  Dataset = rep(dataset_order, times = length(category_order)),
  ExpectedN = c(62L, 32L, 13L, 6L, 15L, 19L, 28L, 35L)
)
actual_counts <- values[, .(ActualN = .N), by = .(Category, Dataset)]
count_check <- merge(expected_counts, actual_counts, by = c("Category", "Dataset"), all = TRUE)
stop_if(nrow(count_check) == 8L, "The category-level Figure 1 input must contain eight dataset groups.")
stop_if(all(count_check$ExpectedN == count_check$ActualN),
  "The category-level Figure 1 observations do not match the confirmed sample counts.")

values[, CategoryLabel := factor(
  Category,
  levels = category_order,
  labels = unname(category_labels[category_order])
)]
values[, Dataset := factor(Dataset, levels = plot_dataset_order)]

stop_if(identical(levels(values$CategoryLabel), unname(category_labels[category_order])),
  "Figure 1 category factor was not constructed.")
stop_if(identical(levels(values$Dataset), plot_dataset_order),
  "Figure 1 dataset factor was not constructed.")

category_counts <- values[, .(N = .N), by = .(CategoryLabel, Dataset)]
category_counts[, Dataset := factor(Dataset, levels = plot_dataset_order)]

figure_plot <- ggplot(
  values,
  aes(x = DdrFractionPercentage, y = Dataset, fill = Dataset)
) +
  geom_boxplot(
    aes(group = Dataset),
    width = 0.60,
    outlier.shape = NA,
    colour = charcoal,
    linewidth = 0.82,
    fatten = 1.35,
    alpha = 0.82,
    na.rm = TRUE
  ) +
  geom_point(
    aes(group = Dataset),
    position = position_jitter(width = 0, height = 0.12, seed = 25),
    shape = 21,
    size = 3.25,
    stroke = 0.62,
    colour = "white",
    alpha = 0.94,
    na.rm = TRUE
  ) +
  geom_text(
    data = category_counts,
    aes(x = 18.0, y = Dataset, label = paste0("n=", N)),
    inherit.aes = FALSE,
    hjust = 0,
    size = 4.30,
    family = publication_font,
    colour = muted_text
  ) +
  facet_grid(
    CategoryLabel ~ .,
    scales = "free_y",
    space = "free_y",
    switch = "y"
  ) +
  scale_fill_manual(
    values = c("Whole proteome" = whole_proteome_colour, "Lactylome (Kla)" = kla_colour),
    breaks = dataset_order
  ) +
  scale_x_continuous(
    limits = c(0, 20),
    breaks = c(0, 5, 10, 15, 20),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_discrete(labels = c("Lactylome (Kla)" = "Kla", "Whole proteome" = "Whole proteome")) +
  guides(
    fill = guide_legend(
      nrow = 1,
      byrow = TRUE,
      keyheight = grid::unit(0.62, "cm"),
      keywidth = grid::unit(0.92, "cm")
    )
  ) +
  labs(
    x = "GO-DDR annotated protein fraction (%)",
    y = NULL,
    fill = NULL
  ) +
  theme_minimal(base_size = 14, base_family = publication_font) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(colour = grid_colour, linewidth = 0.50),
    panel.border = element_rect(colour = panel_border_colour, fill = NA, linewidth = 0.60),
    axis.text.y = element_text(size = 15.5, colour = charcoal, face = "bold"),
    axis.text.x = element_text(size = 15.0, colour = charcoal),
    axis.title.x = element_text(
      size = 19.0,
      face = "bold",
      colour = charcoal,
      margin = margin(t = 14)
    ),
    strip.placement = "outside",
    strip.text.y.left = element_text(
      size = 18.0,
      face = "bold",
      colour = charcoal,
      angle = 90,
      lineheight = 0.95
    ),
    strip.background = element_rect(fill = "#E7E9E7", colour = NA),
    panel.spacing.y = grid::unit(0.72, "lines"),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.text = element_text(size = 16.5, colour = charcoal, lineheight = 1.10),
    legend.key.spacing.x = grid::unit(0.38, "cm"),
    legend.background = element_rect(fill = "white", colour = NA),
    legend.margin = margin(1, 0, 10, 0),
    plot.margin = margin(10, 18, 14, 12),
    plot.background = element_rect(fill = "white", colour = NA)
  )

apply_strip_fills <- function(plot) {
  plot_grob <- ggplotGrob(plot)
  strip_ids <- grep("^strip-l", plot_grob$layout$name)
  strip_ids <- strip_ids[order(plot_grob$layout$t[strip_ids])]
  stop_if(length(strip_ids) == length(category_order),
    "Figure 1 must contain four category strips.")
  for (index in seq_along(strip_ids)) {
    strip_grob <- plot_grob$grobs[[strip_ids[[index]]]]$grobs[[1L]]
    background_id <- grep("^strip.background", strip_grob$childrenOrder)
    stop_if(length(background_id) == 1L,
      "Unable to identify a Figure 1 category-strip background.")
    strip_grob$children[[background_id]]$gp$fill <- unname(category_fills[category_order[index]])
    strip_grob$children[[background_id]]$gp$col <- NA
    plot_grob$grobs[[strip_ids[[index]]]]$grobs[[1L]] <- strip_grob
  }
  plot_grob
}

plot_grob <- apply_strip_fills(figure_plot)
stop_if(length(plot_grob$grobs) > 0L, "Figure 1 category boxplot grob is empty.")

stem <- file.path(output_dir, "Figure_1_DDR_fraction_candidate_category_boxplot_refined")
figure_width <- 15.5
figure_height <- 11.5
ggsave(
  paste0(stem, ".png"),
  plot_grob,
  width = figure_width,
  height = figure_height,
  dpi = 300,
  bg = "white",
  device = ragg::agg_png
)
ggsave(
  paste0(stem, ".pdf"),
  plot_grob,
  width = figure_width,
  height = figure_height,
  bg = "white",
  device = cairo_pdf
)

message(
  "Wrote refined candidate four-category Figure 1 boxplot: ",
  stem, ".png/.pdf (", figure_width, " x ", figure_height, " in; eight boxes)."
)
