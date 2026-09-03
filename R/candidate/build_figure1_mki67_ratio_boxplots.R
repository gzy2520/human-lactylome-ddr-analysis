#!/usr/bin/env Rscript

# Render three isolated Figure 1 candidate plots using whole-proteome
# sample-level MKI67 expression normalized by ACTB, TUBB or H3C1.
# The approved publication figures and their source tables are not changed.

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
charcoal <- "#2F3437"
muted_text <- "#65717D"
grid_colour <- "#D9DDE3"
panel_border_colour <- "#C8CED6"

category_order <- c("normal_tissue", "cancer_tissue", "normal_cells", "cancer_cells")
category_labels <- c(
  normal_tissue = "non-tumor\ntissues",
  cancer_tissue = "tumor\ntissues",
  normal_cells = "normal\ncell lines",
  cancer_cells = "cancer\ncell lines"
)
category_fills <- c(
  normal_tissue = "#DCE9E2",
  cancer_tissue = "#F0DEDE",
  normal_cells = "#E7E1EE",
  cancer_cells = "#EEE4D2"
)
denominator_order <- c("ACTB", "TUBB", "H3C1")
denominator_titles <- c(
  ACTB = "MKI67 / ACTB",
  TUBB = "MKI67 / TUBB",
  H3C1 = "MKI67 / H3C1"
)
denominator_display <- c(
  ACTB = "ACTB (beta-actin)",
  TUBB = "TUBB (beta-tubulin)",
  H3C1 = "H3C1 (histone H3.1)"
)

render_options <- if (length(args) >= 2L) args[-1L] else character()
show_significance <- !("--no-significance" %in% render_options)
render_options <- setdiff(render_options, "--no-significance")
denominators_to_render <- if (length(render_options)) render_options else denominator_order
stop_if(all(denominators_to_render %in% denominator_order),
  "The requested MKI67 ratio denominator is not supported.")

candidate_dir <- normalizePath(Sys.getenv(
  "KLA_MKI67_INPUT", unset = file.path(project_root, "data", "candidate")
), mustWork = TRUE)
output_dir_name <- if (show_significance) "mki67_ratio_boxplot" else "mki67_ratio_boxplot_no_significance"
output_dir <- Sys.getenv(
  "KLA_MKI67_OUTPUT",
  unset = file.path(project_root, "results", "candidate", output_dir_name)
)
values_path <- file.path(candidate_dir, "figure1_mki67_ratio_sample_values.csv")
significance_path <- file.path(project_root, "R", "candidate", "boxplot_significance.R")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

stop_if(file.exists(values_path), paste0("Missing MKI67 ratio input: ", values_path))
stop_if(file.exists(significance_path), paste0("Missing significance helper: ", significance_path))
source(significance_path, local = TRUE)

values <- fread(values_path, check.names = FALSE)
required_values <- c(
  "PXD", "SampleGroup", "Category", "Dataset", "SampleID", "ConditionLabel",
  "ObsKey", "SourceFile", "Denominator", "RatioLabel", "Ratio"
)
stop_if(all(required_values %in% names(values)), "MKI67 ratio input schema is incomplete.")
stop_if(all(values$Dataset == "Whole proteome"), "The ratio plot input contains a non-whole-proteome observation.")
stop_if(all(values$Category %in% category_order), "The ratio plot input contains an unknown category.")
stop_if(all(values$Denominator %in% denominator_order), "The ratio plot input contains an unknown denominator.")
stop_if(all(is.finite(values$Ratio) & values$Ratio > 0), "MKI67 ratios must be finite and positive.")
stop_if(uniqueN(values[, .(ObsKey, Denominator)]) == nrow(values),
  "The ratio plot input contains duplicate observation-denominator records.")
stop_if(all(nzchar(trimws(values$SourceFile))), "A ratio source path is empty.")
stop_if(all(!startsWith(values$SourceFile, "/")), "A ratio source path is absolute.")

values[, X := match(Category, category_order)]
values[, Category := factor(Category, levels = category_order)]
values[, Denominator := factor(Denominator, levels = denominator_order)]

ratio_test_values <- copy(values[, .(
  Category = as.character(Category),
  Denominator = as.character(Denominator),
  Ratio
)])
global_significance <- compute_ratio_global_significance(
  ratio_test_values,
  denominator_order,
  category_order
)
stop_if(nrow(global_significance) == length(denominator_order),
  "MKI67 ratio omnibus significance results are incomplete.")
fwrite(global_significance, file.path(candidate_dir, "figure1_mki67_ratio_significance.csv"))

ratio_ticks <- function(panel_values) {
  panel_values <- panel_values[is.finite(panel_values) & panel_values > 0]
  lower <- 10^floor(log10(min(panel_values)))
  upper <- 10^ceiling(log10(max(panel_values)))
  10^seq(log10(lower), log10(upper), by = 1)
}

format_q_value <- function(q_value) {
  if (!is.finite(q_value)) return("NA")
  formatC(q_value, format = "e", digits = 2)
}

make_plot <- function(denominator) {
  panel <- copy(values[as.character(Denominator) == denominator])
  panel[, Denominator := as.character(Denominator)]
  stop_if(nrow(panel) > 0L, paste0("No data available for ", denominator))
  category_counts <- panel[, .(N = .N, MaxRatio = max(Ratio)), by = .(Category, X)]
  category_counts[, label_y := MaxRatio * 1.18]
  category_stats <- panel[, .(
    Mean = mean(Ratio),
    Median = median(Ratio)
  ), by = .(Category, X)]
  category_stats[, c("x_left", "x_right") := list(X - 0.20, X + 0.20)]

  raw_max <- max(panel$Ratio, na.rm = TRUE)
  raw_min <- min(panel$Ratio, na.rm = TRUE)
  y_min <- 10^(floor(log10(raw_min)) - 0.35)
  preliminary_max <- max(c(raw_max * 1.35, category_counts$label_y * 1.10))
  panel_global <- global_significance[Denominator == denominator]
  global_q <- if (nrow(panel_global)) panel_global$QValueBH[[1L]] else NA_real_
  y_max <- max(preliminary_max, raw_max * 1.50)

  background <- data.table(
    xmin = seq_along(category_order) - 0.5,
    xmax = seq_along(category_order) + 0.5,
    ymin = y_min,
    ymax = y_max,
    Category = factor(category_order, levels = category_order)
  )

  subtitle_text <- if (show_significance) {
    paste0(
      "Four-category one-way ANOVA q=", format_q_value(global_q),
      "\nWhole-proteome source observations; ", denominator_display[[denominator]],
      " used as the denominator"
    )
  } else {
    paste0(
      "Whole-proteome source observations; ", denominator_display[[denominator]],
      " used as the denominator"
    )
  }
  caption_text <- paste(
    "Each point is one source-resolved whole-proteome observation.",
    "Ratios use exact UniProt accessions MKI67=P46013 and",
    switch(denominator, ACTB = "ACTB=P60709", TUBB = "TUBB=P07437", H3C1 = "H3C1=P68431"),
    "with no imputation.",
    "The dark horizontal line inside each box is the median and the red horizontal line is the mean.",
    sep = "\n"
  )
  if (show_significance) {
    caption_text <- paste(
      caption_text,
      "The four-category one-way ANOVA is performed on log10 ratios and BH-adjusted across the three denominators.",
      "**** q<0.0001, *** q<0.001, ** q<0.01, * q<0.05.",
      sep = "\n"
    )
  }

  plot <- ggplot(panel, aes(x = X, y = Ratio, fill = Category)) +
    geom_rect(
      data = background,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = Category),
      inherit.aes = FALSE,
      colour = NA,
      alpha = 0.24
    ) +
    geom_boxplot(
      aes(group = Category),
      width = 0.58,
      outlier.shape = NA,
      colour = charcoal,
      linewidth = 0.90,
      median.linewidth = 1.45,
      alpha = 0.90,
      orientation = "x",
      na.rm = TRUE
    ) +
    geom_segment(
      data = category_stats,
      aes(x = x_left, xend = x_right, y = Median, yend = Median),
      inherit.aes = FALSE,
      colour = charcoal,
      linewidth = 1.15,
      lineend = "round"
    ) +
    geom_segment(
      data = category_stats,
      aes(x = x_left, xend = x_right, y = Mean, yend = Mean),
      inherit.aes = FALSE,
      colour = "#C0392B",
      linewidth = 1.55,
      lineend = "round"
    ) +
    geom_point(
      aes(group = Category),
      position = position_jitter(width = 0.095, height = 0, seed = 25),
      shape = 21,
      size = 3.25,
      stroke = 0.72,
      colour = "white",
      alpha = 0.95,
      na.rm = TRUE
    ) +
    geom_text(
      data = category_counts,
      aes(x = X, y = label_y, label = paste0("n=", N)),
      inherit.aes = FALSE,
      family = publication_font,
      size = 4.55,
      fontface = "bold",
      colour = muted_text,
      vjust = -0.15
    ) +
    scale_fill_manual(values = category_fills, guide = "none", drop = FALSE) +
    scale_x_continuous(
      breaks = seq_along(category_order),
      labels = unname(category_labels[category_order]),
      limits = c(0.5, length(category_order) + 0.5),
      expand = c(0, 0)
    ) +
    scale_y_log10(
      breaks = ratio_ticks(panel$Ratio),
      labels = scales::label_scientific(digits = 2),
      limits = c(y_min, y_max),
      expand = c(0, 0)
    ) +
    labs(
      title = denominator_titles[[denominator]],
      subtitle = subtitle_text,
      x = NULL,
      y = "MKI67 protein intensity ratio (log scale)",
      caption = caption_text
    ) +
    theme_minimal(base_family = publication_font, base_size = 15) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(colour = grid_colour, linewidth = 0.52),
      panel.border = element_rect(colour = panel_border_colour, fill = NA, linewidth = 0.68),
      axis.text.x = element_text(size = 14.8, face = "bold", colour = charcoal, lineheight = 0.95),
      axis.text.y = element_text(size = 14.3, colour = charcoal),
      axis.title.y = element_text(size = 18.5, face = "bold", colour = charcoal, margin = margin(r = 12)),
      plot.title = element_text(
        size = 24.0,
        face = "bold",
        colour = charcoal,
        hjust = 0.5,
        margin = margin(b = 5)
      ),
      plot.subtitle = element_text(size = 14.2, colour = muted_text, hjust = 0.5, margin = margin(b = 12)),
      plot.caption = element_text(size = 10.8, colour = muted_text, hjust = 0, lineheight = 1.05, margin = margin(t = 12)),
      plot.margin = margin(14, 20, 16, 18),
      plot.background = element_rect(fill = "white", colour = NA)
    ) +
    coord_cartesian(clip = "off")

  output_stem <- file.path(output_dir, paste0("Figure_1_MKI67_over_", denominator, "_boxplot"))
  ggsave(paste0(output_stem, ".png"), plot, width = 11.2, height = 8.6, dpi = 300,
    bg = "white", device = ragg::agg_png)
  ggsave(paste0(output_stem, ".pdf"), plot, width = 11.2, height = 8.6,
    bg = "white", device = cairo_pdf)
}

for (denominator in denominators_to_render) make_plot(denominator)

message("Wrote isolated Figure 1 MKI67 ratio boxplots to ", output_dir)
