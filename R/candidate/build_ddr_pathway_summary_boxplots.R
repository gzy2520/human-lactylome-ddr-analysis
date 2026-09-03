#!/usr/bin/env Rscript

# Candidate dataset-level DDR pathway-summary renderers.
#
# There are seven pathways for each modality (Kla and whole proteome). Each
# plot has four category positions and, within each position, adjacent
# Up/positive and Down/negative boxes. No direction is plotted on opposite
# sides of zero. Every point is one PXD/sample-group union.

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
mean_colour <- "#C0392B"
median_colour <- "#2F3437"
muted_text <- "#65717D"
grid_colour <- "#D9DDE3"
panel_border_colour <- "#C8CED6"
negative_fill <- "#7B838B"

category_order <- c("normal_tissue", "cancer_tissue", "cancer_cells", "normal_cells")
category_labels <- c(
  normal_tissue = "non-tumor tissues",
  cancer_tissue = "tumor tissues",
  cancer_cells = "cancer cell lines",
  normal_cells = "normal cell lines"
)
pathway_order <- c("BER", "NER", "MMR", "FA", "HR", "AEJ", "NHEJ")
dataset_order <- c("Lactylome (Kla)", "Whole proteome")
direction_order <- c("Up/positive", "Down/negative")

candidate_dir <- normalizePath(
  Sys.getenv("KLA_CANDIDATE_INPUT", unset = file.path(project_root, "data", "candidate")),
  mustWork = TRUE
)
output_dir <- normalizePath(
  Sys.getenv(
    "KLA_CANDIDATE_OUTPUT",
    unset = file.path(project_root, "results", "candidate", "pathway_summary_dataset_boxplot")
  ),
  mustWork = FALSE
)
values_path <- file.path(candidate_dir, "pathway_summary_dataset_boxplot_values.csv")
pathway_display_path <- file.path(
  normalizePath(
    Sys.getenv("KLA_PUBLICATION_INPUT", unset = file.path(project_root, "data", "publication_input")),
    mustWork = TRUE
  ),
  "pathway_display.csv"
)
significance_path <- file.path(project_root, "R", "candidate", "boxplot_significance.R")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

stop_if(file.exists(values_path), paste0("Missing dataset-level pathway input: ", values_path))
stop_if(file.exists(pathway_display_path), paste0("Missing pathway display input: ", pathway_display_path))
stop_if(file.exists(significance_path), paste0("Missing boxplot significance helper: ", significance_path))

values <- fread(values_path, check.names = FALSE)
pathway_display <- fread(pathway_display_path, check.names = FALSE)
source(significance_path, local = TRUE)

required_values <- c(
  "RowOrder", "PXD", "SampleGroup", "Category", "Dataset", "DatasetPXD",
  "DatasetPointID", "DatasetPointLabel", "SampleID", "SourceFile",
  "Pathway", "PositiveProteinCount", "NegativeProteinCount",
  "AnyPathwayProteinCount", "PathwayScoreMappedProteinCount",
  "DdrProteinCount", "PositiveFraction", "NegativeFraction", "SignedFraction",
  "PathwayScoreCoverage", "PathwayScoreBasis"
)
stop_if(all(required_values %in% names(values)),
  "Dataset-level pathway input schema is incomplete.")
stop_if(nrow(values) > 0L, "Dataset-level pathway input is empty.")
stop_if(all(values$Dataset %in% dataset_order),
  "Dataset-level pathway input contains an unexpected dataset.")
stop_if(all(values$Category %in% category_order),
  "Dataset-level pathway input contains an unknown category.")
stop_if(all(values$Pathway %in% pathway_order),
  "Dataset-level pathway input contains an unknown pathway.")
stop_if(!anyDuplicated(values[, .(DatasetPointID, Pathway)]),
  "Dataset-level pathway input contains duplicate point/pathway rows.")
stop_if(uniqueN(values$DatasetPointID) * length(pathway_order) == nrow(values),
  "Dataset-level pathway input does not have seven pathway rows per point.")
stop_if(all(is.finite(values$PositiveFraction) & is.finite(values$NegativeFraction)),
  "Dataset-level pathway fractions contain non-finite values.")
stop_if(all(values$PositiveFraction >= 0 & values$PositiveFraction <= 1),
  "A positive dataset-level pathway fraction is outside 0-1.")
stop_if(all(values$NegativeFraction >= 0 & values$NegativeFraction <= 1),
  "A negative dataset-level pathway fraction is outside 0-1.")
stop_if(all(is.finite(values$PathwayScoreCoverage) & values$PathwayScoreCoverage >= 0 &
  values$PathwayScoreCoverage <= 1),
  "A pathway-score coverage value is outside 0-1.")

pathway_colours <- setNames(pathway_display$Color, pathway_display$Pathway)
missing_colours <- setdiff(pathway_order, names(pathway_colours))
if (length(missing_colours)) {
  pathway_colours <- c(
    pathway_colours,
    setNames(rep("#4E79A7", length(missing_colours)), missing_colours)
  )
}
pathway_colours <- pathway_colours[pathway_order]

pathway_anova <- compute_pathway_two_way_anova(
  values,
  category_order = category_order,
  pathway_order = pathway_order,
  dataset_order = dataset_order
)
stop_if(nrow(pathway_anova) == length(dataset_order) * length(pathway_order) * 3L,
  "Dataset-level pathway two-way ANOVA results are incomplete.")
stop_if(all(is.finite(pathway_anova$PValue) & is.finite(pathway_anova$QValueBH)),
  "Dataset-level pathway ANOVA results contain a non-finite p or q value.")
fwrite(pathway_anova, file.path(output_dir, "pathway_summary_two_way_anova.csv"), na = "")

format_p <- function(value) {
  if (!is.finite(value)) return("NA")
  if (value < 0.001) formatC(value, format = "e", digits = 2)
  else formatC(value, format = "f", digits = 3)
}

make_anova_label <- function(dataset, pathway) {
  stats <- pathway_anova[Dataset == dataset & Pathway == pathway]
  term_labels <- c(
    "CategoryFactor" = "category",
    "DirectionFactor" = "direction",
    "CategoryFactor:DirectionFactor" = "interaction"
  )
  paste(vapply(seq_len(nrow(stats)), function(index) {
    row <- stats[index]
    term <- if (row$Term %in% names(term_labels)) term_labels[[row$Term]] else row$Term
    paste0(term, " q=", format_p(row$QValueBH), " ", row$Significance)
  }, character(1)), collapse = "\n")
}

make_panel <- function(dataset, pathway) {
  panel <- copy(values[Dataset == dataset & Pathway == pathway])
  panel[, CategoryIndex := match(Category, category_order)]

  plot_data <- rbindlist(list(
    panel[, .(
      DatasetPointID, DatasetPointLabel, Category, CategoryIndex,
      Direction = direction_order[[1L]],
      ValuePercent = PositiveFraction * 100
    )],
    panel[, .(
      DatasetPointID, DatasetPointLabel, Category, CategoryIndex,
      Direction = direction_order[[2L]],
      ValuePercent = NegativeFraction * 100
    )]
  ), fill = TRUE)
  plot_data[, Direction := factor(Direction, levels = direction_order)]

  box_stats <- plot_data[, .(
    Mean = mean(ValuePercent),
    Median = median(ValuePercent),
    N = .N
  ), by = .(Category, Direction, CategoryIndex)]
  box_stats[, Direction := factor(Direction, levels = direction_order)]
  box_stats[, x_center := CategoryIndex + fifelse(
    as.character(Direction) == direction_order[[1L]], -0.19, 0.19
  )]
  box_stats[, c("x_left", "x_right") := list(x_center - 0.21, x_center + 0.21)]

  y_max <- max(plot_data$ValuePercent)
  y_limit <- max(10, ceiling((y_max + 2) / 5) * 5)
  y_breaks <- seq(0, y_limit, by = 5)
  if (tail(y_breaks, 1L) < y_limit) y_breaks <- c(y_breaks, y_limit)

  coverage_text <- sprintf(
    "%.1f-%.1f%%",
    min(panel$PathwayScoreCoverage) * 100,
    max(panel$PathwayScoreCoverage) * 100
  )
  anova_label <- make_anova_label(dataset, pathway)

  figure_plot <- ggplot(
    plot_data,
    aes(x = CategoryIndex, y = ValuePercent, fill = Direction)
  ) +
    geom_boxplot(
      aes(group = Direction),
      position = position_dodge(width = 0.76),
      width = 0.62,
      outlier.shape = NA,
      colour = median_colour,
      linewidth = 0.78,
      alpha = 0.86,
      orientation = "x",
      na.rm = TRUE
    ) +
    geom_point(
      aes(group = Direction),
      position = position_jitterdodge(
        jitter.width = 0.045,
        jitter.height = 0,
        dodge.width = 0.76,
        seed = 25
      ),
      shape = 21,
      size = 3.05,
      stroke = 0.56,
      colour = "white",
      alpha = 0.94,
      na.rm = TRUE
    ) +
    geom_segment(
      data = box_stats,
      aes(x = x_left, xend = x_right, y = Median, yend = Median),
      inherit.aes = FALSE,
      colour = median_colour,
      linewidth = 1.12,
      lineend = "round"
    ) +
    geom_segment(
      data = box_stats,
      aes(x = x_left, xend = x_right, y = Mean, yend = Mean),
      inherit.aes = FALSE,
      colour = mean_colour,
      linewidth = 1.52,
      lineend = "round"
    ) +
    scale_fill_manual(
      values = c(
        "Up/positive" = unname(pathway_colours[[pathway]]),
        "Down/negative" = negative_fill
      ),
      breaks = direction_order
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
      title = paste0(pathway, " · ", dataset),
      subtitle = paste0(
        "Up/positive and Down/negative are adjacent within each category; ",
        "one point = one PXD/sample-group union.\n",
        "Two-way ANOVA (BH-adjusted):\n", anova_label
      ),
      x = NULL,
      y = "Fraction of modality-specific DDR proteins (%)",
      fill = NULL,
      caption = paste(
        "The dark horizontal line inside each box is the median; the red horizontal line is the mean.",
        "Fractions use all modality-specific DDR proteins as the denominator.",
        paste0("Whole-proteome pathway states use S4 BaseAccession matches only; mapped coverage range: ", coverage_text, "."),
        "BH adjustment is across 42 terms from 14 plots: Category, Direction and Category × Direction (**** q<0.0001, *** q<0.001, ** q<0.01, * q<0.05).",
        sep = "\n"
      )
    ) +
    theme_minimal(base_size = 14, base_family = publication_font) +
    theme(
      panel.grid.major.y = element_line(colour = grid_colour, linewidth = 0.50),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.border = element_rect(colour = panel_border_colour, fill = NA, linewidth = 0.60),
      axis.text.x = element_text(size = 13.5, colour = median_colour, face = "bold", lineheight = 0.96),
      axis.text.y = element_text(size = 13.5, colour = median_colour),
      axis.title.y = element_text(size = 17.0, face = "bold", colour = median_colour, margin = margin(r = 12)),
      plot.title = element_text(size = 20.0, face = "bold", colour = unname(pathway_colours[[pathway]]), hjust = 0.5, margin = margin(b = 6)),
      plot.subtitle = element_text(size = 11.0, colour = muted_text, hjust = 0.5, margin = margin(b = 11)),
      plot.caption = element_text(size = 10.0, colour = muted_text, hjust = 0, lineheight = 1.05, margin = margin(t = 10)),
      legend.position = "top",
      legend.direction = "horizontal",
      legend.text = element_text(size = 14.5, colour = median_colour),
      legend.key.spacing.x = grid::unit(0.38, "cm"),
      legend.background = element_rect(fill = "white", colour = NA),
      legend.margin = margin(1, 0, 10, 0),
      plot.margin = margin(12, 18, 14, 16),
      plot.background = element_rect(fill = "white", colour = NA)
    )

  safe_dataset <- if (dataset == "Lactylome (Kla)") "Kla" else "whole_proteome"
  stem <- paste0(
    "Figure_2_DDR_pathway_summary_",
    pathway,
    "_",
    safe_dataset,
    "_boxplot"
  )
  png_path <- file.path(output_dir, paste0(stem, ".png"))
  pdf_path <- file.path(output_dir, paste0(stem, ".pdf"))
  ggsave(
    png_path, figure_plot, width = 13.5, height = 8.8, dpi = 300,
    bg = "white", device = ragg::agg_png
  )
  ggsave(
    pdf_path, figure_plot, width = 13.5, height = 8.8, bg = "white",
    device = cairo_pdf
  )
  data.table(
    Dataset = dataset,
    Pathway = pathway,
    OutputStem = stem,
    PNG = paste0(stem, ".png"),
    PDF = paste0(stem, ".pdf"),
    PointN = uniqueN(panel$DatasetPointID),
    CategoryN = uniqueN(panel$Category),
    PathwayScoreCoverageMin = min(panel$PathwayScoreCoverage),
    PathwayScoreCoverageMax = max(panel$PathwayScoreCoverage)
  )
}

manifest <- rbindlist(lapply(dataset_order, function(dataset) {
  rbindlist(lapply(pathway_order, function(pathway) make_panel(dataset, pathway)), fill = TRUE)
}), fill = TRUE)
manifest[, PathwayOrder := match(Pathway, pathway_order)]
setorder(manifest, Dataset, PathwayOrder)
manifest[, PathwayOrder := NULL]
fwrite(manifest, file.path(output_dir, "pathway_summary_dataset_boxplot_manifest.csv"), na = "")

stop_if(nrow(manifest) == length(dataset_order) * length(pathway_order),
  "Dataset-level pathway boxplot manifest is incomplete.")
stop_if(all(file.exists(file.path(output_dir, manifest$PNG))), "A pathway PNG output is missing.")
stop_if(all(file.exists(file.path(output_dir, manifest$PDF))), "A pathway PDF output is missing.")

message(
  "Wrote ", nrow(manifest), " dataset-level DDR pathway-summary boxplots (",
  length(pathway_order), " pathways x ", length(dataset_order), " modalities) to ", output_dir,
  ". Two-way ANOVA terms and BH-adjusted q values are in pathway_summary_two_way_anova.csv."
)
