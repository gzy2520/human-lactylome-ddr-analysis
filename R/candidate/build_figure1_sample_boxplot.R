#!/usr/bin/env Rscript

# Candidate-only Figure 1 renderer.
#
# The approved publication figure is not changed here.  This figure shows the
# distribution of source-defined sample-level Kla-DDR fractions within each
# publication group.  The blue diamond is the frozen whole-proteome, group-level
# reference used in the manuscript; it is not treated as a sample replicate.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
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
kla_border_colour <- "#C36E12"
charcoal <- "#2F3437"
muted_text <- "#4B5563"
grid_colour <- "#D9DDE3"

candidate_dir <- file.path(project_root, "data", "candidate")
output_dir <- file.path(project_root, "results", "candidate")
values_path <- file.path(candidate_dir, "sample_boxplot_values.csv")
reconciliation_path <- file.path(candidate_dir, "sample_boxplot_reconciliation.csv")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

stop_if(file.exists(values_path), paste0("Missing sample-level input: ", values_path))
stop_if(file.exists(reconciliation_path), paste0("Missing reconciliation input: ", reconciliation_path))

values <- fread(values_path, check.names = FALSE)
reconciliation <- fread(reconciliation_path, check.names = FALSE)

required_columns <- c(
  "RowOrder", "PXD", "SampleGroup", "Category", "DisplayGroup", "SampleID",
  "ConditionLabel", "SampleClass", "ObservationType", "KlaProteinCount",
  "KlaDdrProteinCount", "KlaDdrFractionPercentage", "ReferenceFraction"
)
stop_if(all(required_columns %in% names(values)), "Sample-level input is missing required columns.")
stop_if(nrow(values) > 0L, "Sample-level input is empty.")
stop_if(all(is.finite(values$KlaDdrFractionPercentage)), "Sample-level fractions contain non-finite values.")
stop_if(all(values$KlaDdrFractionPercentage >= 0 & values$KlaDdrFractionPercentage <= 100),
  "Sample-level fractions must be percentages between 0 and 100.")
stop_if(!anyDuplicated(values[, .(PXD, SampleGroup, SampleID)]),
  "A sample-level input row is duplicated.")
stop_if(nrow(reconciliation) == 30L, "Sample-level reconciliation must cover 30 publication groups.")
stop_if(all(reconciliation$GroupUnionStatus == "PASS"),
  "Sample-level inputs do not reconcile to the frozen publication groups.")

category_order <- c("normal_tissue", "cancer_tissue", "cancer_cells", "normal_cells")
category_labels <- c(
  normal_tissue = "Non-tumor tissues",
  cancer_tissue = "Tumor tissues",
  cancer_cells = "Cancer cell lines",
  normal_cells = "Normal cell lines"
)
category_fills <- c(
  normal_tissue = "#DCEAF5",
  cancer_tissue = "#FCE7D4",
  cancer_cells = "#FCE7D4",
  normal_cells = "#DCEAF5"
)
stop_if(all(values$Category %in% category_order), "Unknown publication category in sample-level input.")

group_counts <- values[, .(N = .N), by = .(RowOrder, PXD, SampleGroup, Category, DisplayGroup)]
group_meta <- unique(values[, .(
  RowOrder,
  PXD,
  SampleGroup,
  Category,
  DisplayGroup,
  ReferenceFraction
)])
stop_if(nrow(group_meta) == 30L, "The sample-level plot must contain exactly 30 publication groups.")
stop_if(!anyDuplicated(group_meta[, .(PXD, SampleGroup)]), "Publication groups are not unique.")

group_meta <- merge(
  group_meta,
  group_counts,
  by = c("RowOrder", "PXD", "SampleGroup", "Category", "DisplayGroup"),
  all.x = TRUE,
  sort = FALSE
)
setorder(group_meta, RowOrder)
group_meta[, RowLabel := paste0(DisplayGroup, " · ", PXD)]
group_meta[, PlotRow := factor(RowLabel, levels = rev(RowLabel))]
values <- merge(
  values,
  group_meta[, .(PXD, SampleGroup, RowLabel, PlotRow)],
  by = c("PXD", "SampleGroup"),
  all.x = TRUE,
  sort = FALSE
)

make_panel <- function(category, show_x_title = FALSE, show_x_text = TRUE, show_legend = FALSE) {
  panel_meta <- group_meta[Category == category]
  panel_values <- values[Category == category]
  panel_values[, PlotRow := factor(RowLabel, levels = rev(panel_meta$RowLabel))]
  panel_values[, StripLabel := category_labels[[category]]]
  panel_reference <- panel_meta[, .(
    PlotRow = factor(RowLabel, levels = rev(panel_meta$RowLabel)),
    ReferenceFraction,
    StripLabel = category_labels[[category]]
  )]
  panel_meta[, StripLabel := category_labels[[category]]]

  ggplot(panel_values, aes(x = KlaDdrFractionPercentage, y = PlotRow)) +
    geom_boxplot(
      aes(group = PlotRow, fill = "Kla sample-level"),
      width = 0.56,
      outlier.shape = NA,
      colour = kla_border_colour,
      linewidth = 0.55,
      na.rm = TRUE
    ) +
    geom_point(
      aes(fill = "Kla sample-level"),
      position = position_jitter(width = 0, height = 0.085, seed = 25),
      shape = 21,
      size = 2.65,
      stroke = 0.45,
      colour = "white",
      na.rm = TRUE
    ) +
    geom_point(
      data = panel_reference,
      aes(x = ReferenceFraction, y = PlotRow, shape = "Whole proteome reference"),
      inherit.aes = FALSE,
      size = 3.75,
      stroke = 0.65,
      fill = whole_proteome_colour,
      colour = "white"
    ) +
    geom_text(
      data = panel_meta,
      aes(
        x = 17.35,
        y = PlotRow,
        label = paste0("n=", N)
      ),
      inherit.aes = FALSE,
      hjust = 0,
      size = 4.05,
      family = publication_font,
      colour = charcoal
    ) +
    scale_fill_manual(
      values = c("Kla sample-level" = kla_colour),
      breaks = "Kla sample-level",
      name = NULL
    ) +
    scale_shape_manual(
      values = c("Whole proteome reference" = 23),
      breaks = "Whole proteome reference",
      name = NULL
    ) +
    scale_x_continuous(
      limits = c(0, 18.5),
      breaks = c(0, 5, 10, 15),
      expand = expansion(mult = c(0, 0))
    ) +
    scale_y_discrete() +
    facet_grid(
      StripLabel ~ .,
      scales = "free_y",
      space = "free_y",
      switch = "y"
    ) +
    labs(
      x = if (show_x_title) "GO-DDR annotated Kla protein fraction (%)" else NULL,
      y = NULL
    ) +
    theme_minimal(base_family = publication_font, base_size = 12) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(colour = grid_colour, linewidth = 0.45),
      panel.border = element_rect(colour = "#9AA5B1", fill = NA, linewidth = 0.65),
      panel.background = element_rect(fill = "white", colour = NA),
      axis.text.y = element_text(size = 11.5, colour = charcoal, lineheight = 0.92),
      axis.text.x = if (show_x_text) element_text(size = 13.2, colour = charcoal) else element_blank(),
      axis.ticks.x = if (show_x_text) element_line(colour = "#8D99A6", linewidth = 0.35) else element_blank(),
      axis.line.x = if (show_x_text) element_line(colour = "#8D99A6", linewidth = 0.35) else element_blank(),
      axis.title.x = element_text(
        size = 16.2,
        face = "bold",
        colour = charcoal,
        margin = margin(t = 10)
      ),
      plot.title = element_text(
        size = 18,
        face = "bold",
        colour = charcoal,
        hjust = 0,
        margin = margin(b = 6)
      ),
      strip.placement = "outside",
      strip.text.y.left = element_text(
        size = 13.5,
        face = "bold",
        colour = charcoal,
        angle = 90
      ),
      strip.background = element_rect(
        fill = unname(category_fills[[category]]),
        colour = NA
      ),
      strip.switch.pad.grid = grid::unit(0.16, "lines"),
      legend.position = if (show_legend) "bottom" else "none",
      legend.text = element_text(size = 13.2, colour = charcoal),
      legend.key.height = grid::unit(0.42, "cm"),
      legend.key.width = grid::unit(0.85, "cm"),
      plot.margin = margin(7, 10, if (show_x_title) 9 else 3, 10),
      plot.background = element_rect(fill = "white", colour = NA)
    )
}

# The two tissue categories form one top block.  The two cell-line categories
# retain separate full-width rows below, preserving the original group order.
tissue_block <- patchwork::wrap_plots(
  make_panel("normal_tissue", show_x_text = TRUE),
  make_panel("cancer_tissue", show_x_text = TRUE),
  nrow = 1,
  widths = c(1.95, 1.05),
  guides = "keep"
) +
  plot_annotation(
    title = "Tissues",
    theme = theme(
      plot.title = element_text(
        family = publication_font,
        size = 18,
        face = "bold",
        colour = charcoal,
        hjust = 0,
        margin = margin(b = 4)
      )
    )
  )

final_plot <- tissue_block /
  make_panel("cancer_cells", show_x_text = TRUE) /
  make_panel("normal_cells", show_x_title = TRUE, show_x_text = TRUE, show_legend = TRUE)

final_plot <- final_plot + plot_layout(guides = "keep", heights = c(1.0, 1.0, 1.0))

final_plot <- final_plot +
  plot_annotation(
    title = "Sample-level Kla–DDR fractions across the 30 publication groups",
    subtitle = paste(
      "Orange boxes and points show the distribution of source-defined sample observations within each group.",
      "Blue diamonds show the frozen whole-proteome group reference."
    ),
    caption = paste(
      "n is the number of observations used for each boxplot. A single observation is retained when the deposited",
      "processed data do not support independent biological sample identities; technical channels and pooled runs",
      "are not counted as biological replicates."
    ),
    theme = theme(
      plot.title = element_text(
        family = publication_font,
        size = 24,
        face = "bold",
        colour = charcoal,
        hjust = 0,
        margin = margin(b = 5)
      ),
      plot.subtitle = element_text(
        family = publication_font,
        size = 14,
        colour = muted_text,
        hjust = 0,
        lineheight = 1.08,
        margin = margin(b = 8)
      ),
      plot.caption = element_text(
        family = publication_font,
        size = 11.2,
        colour = muted_text,
        hjust = 0,
        lineheight = 1.05,
        margin = margin(t = 8)
      ),
      plot.margin = margin(12, 16, 12, 16)
    )
  )

plot_grob <- patchwork::patchworkGrob(final_plot)
stop_if(length(plot_grob$grobs) > 0L, "Sample-level boxplot could not be assembled.")

stem <- file.path(output_dir, "Figure_1_candidate_sample_boxplot")
ggsave(
  paste0(stem, ".png"),
  plot_grob,
  width = 23,
  height = 25,
  dpi = 300,
  bg = "white",
  device = ragg::agg_png
)
ggsave(
  paste0(stem, ".pdf"),
  plot_grob,
  width = 23,
  height = 25,
  bg = "white",
  device = cairo_pdf
)
message("Wrote candidate sample-level boxplot: ", stem, ".png/.pdf")
