#!/usr/bin/env Rscript

# Candidate Figure 1 renderer.
#
# This keeps the approved Figure 1 geometry: one vertical panel, four
# publication categories, and two consecutive rows per publication group.
# The historical bars are replaced by sample-level horizontal boxplots with
# the individual source observations overlaid.  The approved publication
# renderer is not changed.

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
  normal_tissue = "#DCEAF5",
  cancer_tissue = "#FCE7D4",
  normal_cells = "#DCEAF5",
  cancer_cells = "#FCE7D4"
)

candidate_dir <- normalizePath(
  Sys.getenv("KLA_CANDIDATE_INPUT", unset = file.path(project_root, "data", "candidate")),
  mustWork = TRUE
)
output_dir <- normalizePath(
  Sys.getenv("KLA_CANDIDATE_OUTPUT", unset = file.path(project_root, "results", "candidate")),
  mustWork = FALSE
)
values_path <- file.path(candidate_dir, "figure1_sample_boxplot_values.csv")
count_path <- file.path(candidate_dir, "biological_sample_count_record.csv")
expected_group_count <- as.integer(Sys.getenv("KLA_CANDIDATE_EXPECTED_GROUPS", unset = "30"))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

stop_if(file.exists(values_path), paste0("Missing Figure 1 sample input: ", values_path))
stop_if(file.exists(count_path), paste0("Missing sample-count record: ", count_path))

values <- fread(values_path, check.names = FALSE)
counts <- fread(count_path, check.names = FALSE)

required_values <- c(
  "RowOrder", "PXD", "SampleGroup", "Category", "DisplayGroup", "Dataset",
  "SampleID", "ConditionLabel", "SampleClass", "ObservationType", "SourceMode",
  "SourceFile", "ReferencePXD", "ProteinCount", "DdrProteinCount",
  "DdrFractionPercentage"
)
stop_if(all(required_values %in% names(values)), "Figure 1 sample input schema is incomplete.")
stop_if(nrow(counts) == expected_group_count,
  paste0("The sample-count record must cover ", expected_group_count, " publication groups."))
stop_if(!anyDuplicated(counts[, .(PXD, SampleGroup)]), "The sample-count record contains duplicate groups.")
stop_if(all(values$Dataset %in% c("Whole proteome", "Lactylome (Kla)")),
  "Figure 1 sample input contains an unknown dataset.")
stop_if(all(values$Category %in% category_order), "Figure 1 sample input contains an unknown category.")
stop_if(all(is.finite(values$DdrFractionPercentage)), "Figure 1 sample fractions contain non-finite values.")
stop_if(all(values$DdrFractionPercentage >= 0 & values$DdrFractionPercentage <= 100),
  "Figure 1 sample fractions must be between 0 and 100 percent.")
stop_if(!anyDuplicated(values[, .(PXD, SampleGroup, Dataset, SampleID)]),
  "Figure 1 sample observations are duplicated.")

expected_n <- melt(
  counts,
  id.vars = c("RowOrder", "PXD", "SampleGroup"),
  measure.vars = c("KlaSampleCount", "ReferenceSampleCount"),
  variable.name = "CountType",
  value.name = "ExpectedN"
)
expected_n[, Dataset := fifelse(CountType == "KlaSampleCount", "Lactylome (Kla)", "Whole proteome")]
actual_n <- values[, .(ActualN = .N), by = .(RowOrder, PXD, SampleGroup, Dataset)]
aggregate_n <- values[, .(
  AggregateOnly = all(ObservationType == "aggregate"),
  ObservationTypes = paste(sort(unique(ObservationType)), collapse = ";")
), by = .(RowOrder, PXD, SampleGroup, Dataset)]
n_check <- merge(
  expected_n[, .(RowOrder, PXD, SampleGroup, Dataset, ExpectedN)],
  actual_n,
  by = c("RowOrder", "PXD", "SampleGroup", "Dataset"),
  all = TRUE
)
n_check <- merge(
  n_check,
  aggregate_n,
  by = c("RowOrder", "PXD", "SampleGroup", "Dataset"),
  all.x = TRUE
)
stop_if(nrow(n_check) == expected_group_count * 2L,
  "Figure 1 sample input does not contain both datasets for all groups.")
count_ok <- !is.na(n_check$ExpectedN) & !is.na(n_check$ActualN) & (
  n_check$ExpectedN == n_check$ActualN |
    (n_check$ActualN == 1L & n_check$ExpectedN > 1L & n_check$AggregateOnly)
)
stop_if(all(count_ok), paste(
  "Figure 1 sample counts disagree with the count record; only a single",
  "explicit aggregate observation may represent multiple source replicates."
))

values[, CategoryLabel := factor(
  Category,
  levels = category_order,
  labels = unname(category_labels[category_order])
)]
values[, Dataset := factor(Dataset, levels = c("Whole proteome", "Lactylome (Kla)"))]
values[, DatasetShort := fifelse(Dataset == "Whole proteome", "whole proteome", "Kla")]
values[, RowID := paste(PXD, SampleGroup, Dataset, sep = "__")]
values[, RowLabel := fifelse(
  Dataset == "Whole proteome",
  paste0(DisplayGroup, " · whole proteome Ref:", ReferencePXD),
  paste0(DisplayGroup, " · Kla:", PXD)
)]
row_order <- unique(values[order(RowOrder, Dataset), RowID])
row_labels <- unique(values[order(RowOrder, Dataset), .(RowID, RowLabel)])
values[, PlotRow := factor(RowID, levels = rev(row_order))]
values[, N := .N, by = .(RowID)]
row_counts <- unique(values[, .(CategoryLabel, PlotRow, N)])
row_stats <- values[, .(
  Mean = mean(DdrFractionPercentage),
  Median = median(DdrFractionPercentage)
), by = PlotRow]
fwrite(
  row_stats[, .(PlotRow, Mean, Median)],
  file.path(output_dir, "figure1_original_boxplot_mean_median.csv"),
  na = ""
)

stop_if(identical(levels(values$CategoryLabel), unname(category_labels[category_order])),
  "Figure 1 category factor was not constructed.")
stop_if(length(row_order) == expected_group_count * 2L,
  paste0("Figure 1 must contain ", expected_group_count * 2L, " dataset rows."))

figure_plot <- ggplot(
  values,
  aes(x = DdrFractionPercentage, y = PlotRow, fill = Dataset)
) +
  geom_boxplot(
    aes(group = PlotRow),
    width = 0.68,
    outlier.shape = NA,
    colour = charcoal,
    linewidth = 0.55,
    fatten = 1.35,
    orientation = "y",
    na.rm = TRUE
  ) +
  geom_point(
    aes(group = PlotRow),
    position = position_jitter(width = 0, height = 0.08, seed = 25),
    shape = 21,
    size = 2.05,
    stroke = 0.42,
    colour = "white",
    na.rm = TRUE
  ) +
  geom_point(
    data = row_stats,
    aes(x = Mean, y = PlotRow),
    inherit.aes = FALSE,
    shape = 124,
    size = 8.2,
    colour = "#C0392B",
    na.rm = TRUE
  ) +
  geom_text(
    data = row_counts,
    aes(x = 17.3, y = PlotRow, label = paste0("n=", N)),
    inherit.aes = FALSE,
    hjust = 0,
    size = 3.0,
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
    breaks = c("Whole proteome", "Lactylome (Kla)")
  ) +
  scale_x_continuous(
    limits = c(0, 19),
    breaks = c(0, 5, 10, 15),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_discrete(labels = setNames(row_labels$RowLabel, row_labels$RowID)) +
  guides(
    fill = guide_legend(
      ncol = 1,
      byrow = TRUE,
      keyheight = grid::unit(0.84, "cm"),
      keywidth = grid::unit(1.02, "cm")
    )
  ) +
  labs(
    x = "GO-DDR annotated protein fraction (%)",
    y = NULL,
    fill = NULL,
    caption = paste(
      "Each point is a source-resolved observation; one aggregate point is retained when the source only provides an aggregate profile.",
      "The box center line is the median and the red vertical line is the mean.",
      sep = "\n"
    )
  ) +
  theme_minimal(base_size = 10, base_family = publication_font) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(colour = grid_colour, linewidth = 0.45),
    panel.border = element_rect(colour = panel_border_colour, fill = NA, linewidth = 0.42),
    axis.text.y = element_text(size = 7.4, colour = charcoal, lineheight = 0.92),
    axis.text.x = element_text(size = 11.5, colour = charcoal),
    axis.title.x = element_text(
      size = 16.5,
      face = "bold",
      colour = charcoal,
      margin = margin(t = 12)
    ),
    strip.placement = "outside",
    strip.text.y.left = element_text(
      size = 15,
      face = "bold",
      colour = charcoal,
      angle = 90
    ),
    strip.background = element_rect(fill = "#DCEAF5", colour = NA),
    panel.spacing.y = grid::unit(0.72, "lines"),
    legend.position = "inside",
    legend.position.inside = c(0.80, 0.992),
    legend.justification.inside = c(1, 1),
    legend.direction = "vertical",
    legend.text = element_text(size = 17.5, colour = charcoal, lineheight = 1.12),
    legend.key.spacing.y = grid::unit(0.24, "cm"),
    legend.background = element_rect(
      fill = scales::alpha("white", 0.92),
      colour = panel_border_colour,
      linewidth = 0.45
    ),
    legend.margin = margin(11, 14, 11, 14),
    plot.margin = margin(8, 16, 12, 12),
    plot.background = element_rect(fill = "white", colour = NA)
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
    strip_grob$children[[background_id]]$gp$fill <- unname(category_fills[category_order[index]])
    strip_grob$children[[background_id]]$gp$col <- NA
    plot_grob$grobs[[strip_ids[[index]]]]$grobs[[1L]] <- strip_grob
  }
  plot_grob
}

plot_grob <- apply_strip_fills(figure_plot)
stop_if(length(plot_grob$grobs) > 0L, "Figure 1 boxplot grob is empty.")

stem <- file.path(output_dir, "Figure_1_DDR_fraction_candidate_sample_boxplot")
figure_height <- max(18, nrow(row_labels) * 0.30 + 4.2)
ggsave(
  paste0(stem, ".png"),
  plot_grob,
  width = 15.5,
  height = figure_height,
  dpi = 300,
  bg = "white",
  device = ragg::agg_png
)
ggsave(
  paste0(stem, ".pdf"),
  plot_grob,
  width = 15.5,
  height = figure_height,
  bg = "white",
  device = cairo_pdf
)

message(
  "Wrote candidate original-layout Figure 1 boxplot: ",
  stem, ".png/.pdf (15.5 x ", sprintf("%.1f", figure_height), " in)."
)
