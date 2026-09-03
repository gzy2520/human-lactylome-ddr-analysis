#!/usr/bin/env Rscript

# Render seven source-sample Kla-DDR pathway summaries. Each pathway receives
# one upright figure with four biological-category columns and paired Pro plus
# Inh bar plots with SEM error bars.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

set.seed(25)

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1L) normalizePath(args[[1L]], mustWork = TRUE) else normalizePath(".", mustWork = TRUE)
stop_if <- function(condition, message) if (!isTRUE(condition)) stop(message, call. = FALSE)

publication_font <- "Arial Unicode MS"
charcoal <- "#2F3437"
muted_text <- "#65717D"
grid_colour <- "#D9DDE3"
panel_border_colour <- "#C8CED6"
down_colour <- "#98A1AA"

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
pathway_order <- c("BER", "NER", "MMR", "FA", "HR", "AEJ", "NHEJ")
pathway_colours <- c(
  BER = "#54BED4", NER = "#F59E83", MMR = "#8C9ABD", FA = "#8ED5C4",
  HR = "#4C669D", AEJ = "#00A98F", NHEJ = "#E94F3D"
)
direction_order <- c("Pro", "Inh")

candidate_dir <- normalizePath(Sys.getenv(
  "KLA_CANDIDATE_INPUT", unset = file.path(project_root, "data", "candidate")
), mustWork = TRUE)
output_dir <- Sys.getenv(
  "KLA_CANDIDATE_OUTPUT",
  unset = file.path(project_root, "results", "candidate", "pathway_summary_by_pathway")
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
values_path <- file.path(candidate_dir, "figure1_pathway_summary_sample_boxplot_values.csv")
stop_if(file.exists(values_path), paste0("Missing sample-level pathway input: ", values_path))
values <- fread(values_path, check.names = FALSE)

required_values <- c(
  "PXD", "SampleGroup", "Category", "Dataset", "SampleID", "Pathway",
  "PositiveFraction", "NegativeFraction", "KlaDdrProteinCount", "SourceFile"
)
stop_if(all(required_values %in% names(values)), "Pathway summary input schema is incomplete.")
stop_if(all(as.character(values$Dataset) == "Lactylome (Kla)"),
  "This requested seven-panel pathway summary must use source-resolved Kla observations.")
stop_if(setequal(unique(values$Category), category_order), "Pathway summary categories are incomplete.")
stop_if(setequal(unique(values$Pathway), pathway_order), "Pathway summary pathways are incomplete.")
stop_if(!anyDuplicated(values[, .(PXD, SampleGroup, SampleID, Pathway)]),
  "A sample/pathway observation is duplicated.")
stop_if(all(is.finite(values$PositiveFraction) & is.finite(values$NegativeFraction) &
            values$PositiveFraction >= 0 & values$NegativeFraction >= 0),
  "Pathway fractions must be finite and non-negative.")

source(file.path(project_root, "R", "candidate", "boxplot_significance.R"), local = TRUE)
pathway_anova <- compute_pathway_sample_two_way_anova(values, category_order, pathway_order, direction_order = direction_order)
stop_if(nrow(pathway_anova) == length(pathway_order) * 3L && all(is.finite(pathway_anova$PValue)),
  "Pathway two-way ANOVA did not produce 21 finite term tests.")
fwrite(pathway_anova, file.path(output_dir, "pathway_summary_two_way_anova.csv"), na = "")

make_subtitle <- function(pathway) {
  terms <- pathway_anova[Pathway == pathway]
  q_for <- function(term) {
    row <- terms[Term == term]
    if (!nrow(row)) "NA" else paste0(row$Significance[[1L]], " (q=", formatC(row$QValueBH[[1L]], format = "e", digits = 2), ")")
  }
  paste0(
    "Two-way ANOVA, BH-adjusted: category ", q_for("CategoryFactor"),
    "; direction ", q_for("DirectionFactor"),
    "; category × direction ", q_for("CategoryFactor:DirectionFactor")
  )
}

apply_strip_fills <- function(plot) {
  plot_grob <- ggplotGrob(plot)
  strip_ids <- grep("^strip-t", plot_grob$layout$name)
  if (length(strip_ids) == 0L) {
    strip_ids <- grep("^strip-l", plot_grob$layout$name)
    strip_ids <- strip_ids[order(plot_grob$layout$t[strip_ids])]
  } else {
    strip_ids <- strip_ids[order(plot_grob$layout$l[strip_ids])]
  }
  stop_if(length(strip_ids) == length(category_order), "Each pathway plot must contain four category strips.")
  for (index in seq_along(strip_ids)) {
    strip_grob <- plot_grob$grobs[[strip_ids[[index]]]]$grobs[[1L]]
    background_id <- grep("^strip.background", strip_grob$childrenOrder)
    stop_if(length(background_id) == 1L, "Unable to identify a pathway category-strip background.")
    strip_grob$children[[background_id]]$gp$fill <- unname(category_fills[category_order[[index]]])
    strip_grob$children[[background_id]]$gp$col <- NA
    plot_grob$grobs[[strip_ids[[index]]]]$grobs[[1L]] <- strip_grob
  }
  plot_grob
}

manifest <- rbindlist(lapply(pathway_order, function(pathway) {
  pathway_values <- copy(values[Pathway == pathway])
  long <- rbindlist(list(
    pathway_values[, .(PXD, SampleGroup, Category, SampleID, Direction = direction_order[[1L]], ValuePercent = PositiveFraction * 100)],
    pathway_values[, .(PXD, SampleGroup, Category, SampleID, Direction = direction_order[[2L]], ValuePercent = NegativeFraction * 100)]
  ))
  long[, CategoryLabel := factor(Category, levels = category_order, labels = unname(category_labels[category_order]))]
  long[, Direction := factor(Direction, levels = direction_order)]

  stats <- long[, .(
    N = .N,
    Mean = mean(ValuePercent),
    SD = if (.N > 1L) sd(ValuePercent) else 0,
    SEM = if (.N > 1L) sd(ValuePercent) / sqrt(.N) else 0
  ), by = .(Category, CategoryLabel, Direction)]
  stats[, ErrorMin := pmax(0, Mean - SEM)]
  stats[, ErrorMax := Mean + SEM]

  y_limit <- max(15, ceiling(max(c(stats$ErrorMax, long$ValuePercent)) * 1.15 / 5) * 5)
  counts <- stats[Direction == direction_order[[1L]], .(CategoryLabel, N, LabelY = y_limit * 0.92)]

  figure_plot <- ggplot() +
    geom_col(
      data = stats, aes(x = Direction, y = Mean, fill = Direction),
      width = 0.58, colour = charcoal, linewidth = 0.65, alpha = 0.65
    ) +
    geom_point(
      data = long, aes(x = Direction, y = ValuePercent, fill = Direction),
      position = position_jitter(width = 0.16, height = 0, seed = 25),
      shape = 21, size = 2.6, stroke = 0.55, alpha = 0.85, colour = charcoal
    ) +
    geom_errorbar(
      data = stats, aes(x = Direction, ymin = ErrorMin, ymax = ErrorMax),
      width = 0.22, linewidth = 0.90, colour = charcoal
    ) +
    geom_text(
      data = counts, aes(x = 2.35, y = LabelY, label = paste0("n=", N)),
      inherit.aes = FALSE, hjust = 1, size = 4.2, family = publication_font, colour = muted_text
    ) +
    facet_grid(. ~ CategoryLabel) +
    scale_fill_manual(values = c("Pro" = pathway_colours[[pathway]], "Inh" = down_colour)) +
    scale_colour_manual(values = c("Pro" = pathway_colours[[pathway]], "Inh" = down_colour)) +
    scale_y_continuous(
      limits = c(0, y_limit), breaks = scales::pretty_breaks(n = 5),
      labels = function(y) paste0(y, "%"), expand = expansion(mult = c(0, 0))
    ) +
    guides(
      fill = guide_legend(title = NULL, nrow = 1, byrow = TRUE),
      colour = "none"
    ) +
    labs(
      title = paste0(pathway, " pathway"), subtitle = make_subtitle(pathway),
      x = NULL, y = "Relative portion of Kla-DDR proteins (%)",
      caption = paste(
        "Bars represent group mean; error bars indicate ± SEM; points represent individual source-resolved Kla observations.",
        "Pro = positive fraction; Inh = negative fraction. Fractions use each sample's Kla-DDR protein count as denominator.",
        sep = "\n"
      )
    ) +
    theme_minimal(base_size = 14, base_family = publication_font) +
    theme(
      panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(colour = grid_colour, linewidth = 0.50),
      panel.border = element_rect(colour = panel_border_colour, fill = NA, linewidth = 0.60),
      axis.text.x = element_text(size = 15, colour = charcoal, face = "bold"),
      axis.text.y = element_text(size = 14.5, colour = charcoal),
      axis.title.y = element_text(size = 17, face = "bold", colour = charcoal, margin = margin(r = 12)),
      plot.title = element_text(size = 22, face = "bold", colour = pathway_colours[[pathway]], hjust = 0.5, margin = margin(b = 2)),
      plot.subtitle = element_text(size = 11, colour = muted_text, hjust = 0.5, margin = margin(b = 10)),
      strip.placement = "outside",
      strip.text.x.top = element_text(size = 15.5, face = "bold", colour = charcoal, margin = margin(t = 6, b = 6)),
      strip.background = element_rect(fill = "#E7E9E7", colour = NA),
      panel.spacing.x = grid::unit(0.9, "lines"),
      legend.position = "top", legend.text = element_text(size = 13.5, colour = charcoal),
      legend.margin = margin(0, 0, 6, 0),
      plot.caption = element_text(size = 10.5, hjust = 0.5, colour = muted_text, margin = margin(t = 12)),
      plot.margin = margin(12, 18, 14, 14), plot.background = element_rect(fill = "white", colour = NA)
    )

  plot_grob <- apply_strip_fills(figure_plot)
  stem_barplot <- paste0("Figure_2_DDR_pathway_summary_", pathway, "_barplot")
  stem_boxplot <- paste0("Figure_2_DDR_pathway_summary_", pathway, "_boxplot")

  ggsave(file.path(output_dir, paste0(stem_barplot, ".png")), plot_grob, width = 14, height = 9, dpi = 300, bg = "white", device = ragg::agg_png)
  ggsave(file.path(output_dir, paste0(stem_barplot, ".pdf")), plot_grob, width = 14, height = 9, bg = "white", device = cairo_pdf)

  ggsave(file.path(output_dir, paste0(stem_boxplot, ".png")), plot_grob, width = 14, height = 9, dpi = 300, bg = "white", device = ragg::agg_png)
  ggsave(file.path(output_dir, paste0(stem_boxplot, ".pdf")), plot_grob, width = 14, height = 9, bg = "white", device = cairo_pdf)

  data.table(
    Pathway = pathway, Dataset = "Lactylome (Kla)",
    PNG = paste0(stem_boxplot, ".png"), PDF = paste0(stem_boxplot, ".pdf"),
    BarplotPNG = paste0(stem_barplot, ".png"), BarplotPDF = paste0(stem_barplot, ".pdf"),
    CategoryPanels = 4L, BoxesPerFigure = 8L, InputPoints = uniqueN(long$SampleID)
  )
}), fill = TRUE)

fwrite(manifest, file.path(output_dir, "pathway_summary_by_pathway_manifest.csv"), na = "")
message("Wrote seven upright four-category pathway barplots with SEM error bars to ", output_dir)
