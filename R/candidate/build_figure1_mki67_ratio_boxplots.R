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
  normal_tissue = "#0072B2",
  cancer_tissue = "#D55E00",
  normal_cells  = "#009E73",
  cancer_cells  = "#CC79A7"
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
  category_counts[, label_y := 10^(log10(MaxRatio) + 0.22)]
  category_stats <- panel[, .(
    Mean = mean(Ratio),
    Median = median(Ratio)
  ), by = .(Category, X)]
  category_stats[, c("x_left", "x_right") := list(X - 0.20, X + 0.20)]

  raw_max <- max(panel$Ratio, na.rm = TRUE)
  raw_min <- min(panel$Ratio, na.rm = TRUE)
  y_min <- 10^(floor(log10(raw_min)) - 0.35)
  panel_global <- global_significance[Denominator == denominator]
  global_q <- if (nrow(panel_global)) panel_global$QValueBH[[1L]] else NA_real_
  y_max <- 10^(log10(raw_max) + 0.65)

  background <- data.table(
    xmin = seq_along(category_order) - 0.5,
    xmax = seq_along(category_order) + 0.5,
    ymin = y_min,
    ymax = y_max,
    Category = factor(category_order, levels = category_order)
  )

  subtitle_text <- if (show_significance) {
    paste0("Four-category one-way ANOVA q = ", format_q_value(global_q))
  } else {
    NULL
  }
  caption_text <- NULL

  plot <- ggplot(panel, aes(x = X, y = Ratio, fill = Category)) +
    geom_vline(xintercept = c(1.5, 2.5, 3.5), colour = "#E5E7EB", linetype = "dashed", linewidth = 0.5) +
    geom_boxplot(
      aes(group = Category),
      width = 0.58,
      outlier.shape = NA,
      colour = charcoal,
      linewidth = 0.55,
      median.linewidth = 0.75,
      alpha = 0.82,
      orientation = "x",
      na.rm = TRUE
    ) +
    geom_segment(
      data = category_stats,
      aes(x = x_left, xend = x_right, y = Median, yend = Median),
      inherit.aes = FALSE,
      colour = charcoal,
      linewidth = 0.75,
      lineend = "round"
    ) +
    geom_segment(
      data = category_stats,
      aes(x = x_left, xend = x_right, y = Mean, yend = Mean),
      inherit.aes = FALSE,
      colour = "#C0392B",
      linewidth = 0.80,
      lineend = "round"
    ) +
    geom_point(
      aes(group = Category),
      position = position_jitter(width = 0.095, height = 0, seed = 25),
      shape = 21,
      size = 2.6,
      stroke = 0.50,
      colour = "white",
      alpha = 0.90,
      na.rm = TRUE
    ) +
    geom_text(
      data = category_counts,
      aes(x = X, y = label_y, label = paste0("n=", N)),
      inherit.aes = FALSE,
      family = publication_font,
      size = 4.2,
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
    theme_minimal(base_family = publication_font, base_size = 14) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(colour = grid_colour, linewidth = 0.50),
      panel.border = element_blank(),
      axis.line.x = element_line(colour = "#8C939E", linewidth = 0.60),
      axis.line.y = element_line(colour = "#8C939E", linewidth = 0.60),
      axis.text.x = element_text(size = 13.5, face = "bold", colour = charcoal, lineheight = 0.95, margin = margin(t = 6)),
      axis.text.y = element_text(size = 13.0, colour = charcoal),
      axis.title.y = element_text(size = 15.5, face = "bold", colour = charcoal, margin = margin(r = 10)),
      plot.title = element_text(
        size = 18.0,
        face = "bold",
        colour = charcoal,
        hjust = 0.5,
        margin = margin(b = 4)
      ),
      plot.subtitle = element_text(size = 12.0, colour = muted_text, hjust = 0.5, margin = margin(b = 10)),
      plot.caption = element_text(size = 9.8, colour = muted_text, hjust = 0, lineheight = 1.10, margin = margin(t = 10)),
      plot.margin = margin(12, 16, 12, 12),
      plot.background = element_rect(fill = "white", colour = NA)
    ) +
    coord_cartesian(clip = "off")

  output_stem <- file.path(output_dir, paste0("Figure_1_MKI67_over_", denominator, "_boxplot"))
  ggsave(paste0(output_stem, ".png"), plot, width = 8.5, height = 7.0, dpi = 300,
    bg = "white", device = ragg::agg_png)
  ggsave(paste0(output_stem, ".pdf"), plot, width = 8.5, height = 7.0,
    bg = "white", device = cairo_pdf)
}

for (denominator in denominators_to_render) make_plot(denominator)

message("Wrote isolated Figure 1 MKI67 ratio boxplots to ", output_dir)
