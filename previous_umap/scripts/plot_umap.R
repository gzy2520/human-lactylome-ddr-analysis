#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(scatterpie)
  library(ggrepel)
})

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
project_root <- if (length(file_arg) > 0) {
  script_path <- sub("^--file=", "", file_arg[[1]])
  normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
} else {
  normalizePath(getwd(), mustWork = TRUE)
}

set.seed(25)

results_dir <- file.path(project_root, "results")
figures_dir <- file.path(project_root, "figures")
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

data <- read_csv(file.path(results_dir, "umap_pie_data.csv"), show_col_types = FALSE) %>%
  mutate(
    Plot_UMAP_1 = UMAP_1,
    Plot_UMAP_2 = UMAP_2
  )

plot_data <- data %>%
  arrange(Symbol)

missing_legacy <- data %>%
  filter(!coalesce(Has_Legacy_Coordinates, FALSE)) %>%
  arrange(Symbol)

write_csv(plot_data, file.path(results_dir, "umap_plot_data.csv"))
write_csv(missing_legacy, file.path(results_dir, "genes_excluded_from_legacy_umap_plot.csv"))

categories <- c("HR", "NHEJ", "BER", "NER", "MMR", "TLS", "DRR", "CP", "Other")
colors <- c(
  HR = "#56B4E9",
  NHEJ = "#009E73",
  BER = "#F0E442",
  NER = "#E69F00",
  MMR = "#D55E00",
  TLS = "#CC79A7",
  DRR = "#9C755F",
  CP = "#0072B2",
  Other = "#7F7F7F"
)

x_range <- diff(range(plot_data$Plot_UMAP_1, na.rm = TRUE))
y_range <- diff(range(plot_data$Plot_UMAP_2, na.rm = TRUE))
axis_pad <- max(c(x_range, y_range), na.rm = TRUE) * 0.08
if (!is.finite(axis_pad) || axis_pad == 0) axis_pad <- 1

pie_plot <- ggplot() +
  geom_scatterpie(
    data = plot_data,
    aes(x = Plot_UMAP_1, y = Plot_UMAP_2, group = Symbol),
    cols = categories,
    pie_scale = 0.48,
    color = "white",
    linewidth = 0.12,
    alpha = 0.95
  ) +
  geom_text_repel(
    data = plot_data,
    aes(x = Plot_UMAP_1, y = Plot_UMAP_2, label = Symbol),
    seed = 25,
    size = 2.35,
    box.padding = 0.35,
    point.padding = 0.2,
    max.overlaps = 100,
    max.iter = 5000,
    min.segment.length = 0,
    segment.alpha = 0.35,
    segment.size = 0.18
  ) +
  scale_fill_manual(values = colors, breaks = categories) +
  coord_equal(
    xlim = range(plot_data$Plot_UMAP_1, na.rm = TRUE) + c(-axis_pad, axis_pad),
    ylim = range(plot_data$Plot_UMAP_2, na.rm = TRUE) + c(-axis_pad, axis_pad),
    expand = FALSE
  ) +
  labs(
    title = "DNA Repair Pathway Classification",
    x = "1",
    y = "2",
    fill = "Pathway"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right",
    panel.grid = element_line(color = "grey88", linewidth = 0.25)
  )

primary_plot <- ggplot(
  plot_data,
  aes(
    x = Plot_UMAP_1,
    y = Plot_UMAP_2,
    color = Primary_Category,
    shape = Is_Multimechanism
  )
) +
  geom_point(size = 3.1, alpha = 0.92, stroke = 0.9) +
  geom_text_repel(
    aes(label = Symbol),
    seed = 25,
    size = 2.7,
    box.padding = 0.25,
    point.padding = 0.15,
    max.overlaps = 100,
    min.segment.length = 0
  ) +
  scale_color_manual(values = c(colors, Unclassified = "#999999"), na.value = "#999999") +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 21), labels = c(`FALSE` = "Single mechanism", `TRUE` = "Multi-mechanism")) +
  coord_equal(
    xlim = range(plot_data$Plot_UMAP_1, na.rm = TRUE) + c(-axis_pad, axis_pad),
    ylim = range(plot_data$Plot_UMAP_2, na.rm = TRUE) + c(-axis_pad, axis_pad),
    expand = FALSE
  ) +
  labs(
    title = "Primary DNA Repair Category",
    subtitle = "Outlined points indicate genes with evidence for multiple concrete repair mechanisms",
    x = "1",
    y = "2",
    color = "Primary",
    shape = "Evidence"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 9),
    legend.position = "right",
    panel.grid = element_line(color = "grey88", linewidth = 0.25)
  )

ggsave(file.path(figures_dir, "umap_repair_pies.png"), pie_plot, width = 12, height = 8, dpi = 300, limitsize = TRUE)
ggsave(file.path(figures_dir, "umap_repair_pies.pdf"), pie_plot, width = 12, height = 8, limitsize = TRUE)
ggsave(file.path(figures_dir, "umap_repair_primary.png"), primary_plot, width = 12, height = 8, dpi = 300, limitsize = TRUE)
ggsave(file.path(figures_dir, "umap_repair_primary.pdf"), primary_plot, width = 12, height = 8, limitsize = TRUE)

legacy_data <- data %>%
  filter(coalesce(Has_Legacy_Coordinates, FALSE)) %>%
  mutate(
    Plot_UMAP_1 = Legacy_UMAP_1,
    Plot_UMAP_2 = Legacy_UMAP_2
  )

if (nrow(legacy_data) > 0) {
  legacy_x_range <- diff(range(legacy_data$Plot_UMAP_1, na.rm = TRUE))
  legacy_y_range <- diff(range(legacy_data$Plot_UMAP_2, na.rm = TRUE))
  legacy_axis_pad <- max(c(legacy_x_range, legacy_y_range), na.rm = TRUE) * 0.08
  if (!is.finite(legacy_axis_pad) || legacy_axis_pad == 0) legacy_axis_pad <- 1

  legacy_plot <- ggplot() +
    geom_scatterpie(
      data = legacy_data,
      aes(x = Plot_UMAP_1, y = Plot_UMAP_2, group = Symbol),
      cols = categories,
      pie_scale = 0.55,
      color = "white",
      linewidth = 0.15,
      alpha = 0.95
    ) +
    geom_text_repel(
      data = legacy_data,
      aes(x = Plot_UMAP_1, y = Plot_UMAP_2, label = Symbol),
      seed = 25,
      size = 2.35,
      box.padding = 0.25,
      point.padding = 0.15,
      max.overlaps = 100,
      min.segment.length = 0,
      segment.alpha = 0.35,
      segment.size = 0.18
    ) +
    scale_fill_manual(values = colors, breaks = categories) +
    coord_equal(
      xlim = range(legacy_data$Plot_UMAP_1, na.rm = TRUE) + c(-legacy_axis_pad, legacy_axis_pad),
      ylim = range(legacy_data$Plot_UMAP_2, na.rm = TRUE) + c(-legacy_axis_pad, legacy_axis_pad),
      expand = FALSE
    ) +
    labs(
      title = "DNA Repair Pathway Classification (Legacy Layout)",
      x = "1",
      y = "2",
      fill = "Pathway"
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "right",
      panel.grid = element_line(color = "grey88", linewidth = 0.25)
    )

  ggsave(file.path(figures_dir, "umap_repair_pies_legacy_layout.png"), legacy_plot, width = 12, height = 8, dpi = 300, limitsize = TRUE)
  ggsave(file.path(figures_dir, "umap_repair_pies_legacy_layout.pdf"), legacy_plot, width = 12, height = 8, limitsize = TRUE)
}

message("Plotted genes in recalculated UMAP: ", nrow(plot_data))
message("Genes without legacy coordinates: ", nrow(missing_legacy))
message("Wrote figures to: ", figures_dir)
