#!/usr/bin/env Rscript

# Candidate-only Figure 1 layout.
#
# The approved publication figure is not changed here.  This renderer keeps
# the selected publication values and historical group order, but lays them out as
# three horizontal rows: two tissue panels in the first row, followed by the
# normal- and cancer-cell-line panels.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

set.seed(25)

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args)) {
  normalizePath(args[[1L]], mustWork = TRUE)
} else {
  normalizePath(".", mustWork = TRUE)
}

stop_if <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

input_dir <- normalizePath(
  Sys.getenv("KLA_PUBLICATION_INPUT", unset = file.path(project_root, "data", "publication_input")),
  mustWork = TRUE
)
candidate_dir <- normalizePath(
  Sys.getenv("KLA_CANDIDATE_DESIGN_INPUT", unset = file.path(project_root, "data", "candidate")),
  mustWork = TRUE
)
output_dir <- normalizePath(
  Sys.getenv("KLA_CANDIDATE_OUTPUT", unset = file.path(project_root, "results", "candidate")),
  mustWork = FALSE
)
expected_group_count <- as.integer(Sys.getenv("KLA_PUBLICATION_EXPECTED_GROUPS", unset = "30"))
expected_category_counts <- as.integer(strsplit(
  Sys.getenv("KLA_PUBLICATION_CATEGORY_COUNTS", unset = "normal_tissue=9;cancer_tissue=2;cancer_cells=12;normal_cells=7"),
  "[;]"
)[[1L]] |> vapply(function(item) sub("^.*=", "", item), character(1)))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

publication_font <- "Arial Unicode MS"
whole_proteome_colour <- "#4E79A7"
kla_colour <- "#F28E2B"
charcoal <- "#2F3437"
grid_colour <- "#D9DDE3"
line_colour <- "#9AA5B1"

category_order <- c("normal_tissue", "cancer_tissue", "cancer_cells", "normal_cells")
category_labels <- c(
  normal_tissue = "Non-tumor tissues",
  cancer_tissue = "Tumor tissues",
  normal_cells = "Normal cell lines",
  cancer_cells = "Cancer cell lines"
)
category_fills <- c(
  normal_tissue = "#DCEAF5",
  cancer_tissue = "#FCE7D4",
  normal_cells = "#DCEAF5",
  cancer_cells = "#FCE7D4"
)
match_linetypes <- c(
  same_study_sample = "solid",
  same_biospecimen = "solid",
  same_material = "solid",
  external_tissue = "dashed",
  external_adjacent = "dashed",
  external_disease = "dashed",
  external_cell_line = "dashed",
  same_study_caveat = "dotted",
  process_control = "dotted"
)

groups <- fread(file.path(input_dir, "group_summary_30.csv"), check.names = FALSE)
design <- fread(file.path(candidate_dir, "sample_design_30.csv"), na.strings = c("", "NA"))
stop_if(nrow(groups) == expected_group_count,
  paste0("Publication input must contain ", expected_group_count, " groups."))
stop_if(nrow(design) == expected_group_count,
  paste0("Candidate design must contain ", expected_group_count, " groups."))

groups[, GroupKey := paste(PXD, SampleGroup, sep = "__")]
design[, GroupKey := paste(PXD, SampleGroup, sep = "__")]
stop_if(!anyDuplicated(groups$GroupKey), "Frozen publication groups are not unique.")
stop_if(!anyDuplicated(design$GroupKey), "Candidate design rows are not unique.")
stop_if(setequal(groups$GroupKey, design$GroupKey), "Candidate design does not cover the selected publication groups.")

data <- merge(
  groups[, .(
    RowOrder,
    PXD,
    SampleGroup,
    Category,
    DisplayGroup = KlaLabelEn,
    KlaFraction = KlaDdrFraction * 100,
    KlaDdr = KlaDdrProteinCount,
    KlaTotal = KlaProteinCount,
    ReferenceFraction = ReferenceDdrFraction * 100,
    ReferenceDdr = ReferenceDdrProteinCount,
    ReferenceTotal = ReferenceProteinCount
  )],
  design[, .(
    PXD,
    SampleGroup,
    KlaN,
    ReferenceN,
    MatchClass
  )],
  by = c("PXD", "SampleGroup"),
  all.x = TRUE,
  sort = FALSE
)
setorder(data, RowOrder)
stop_if(!anyNA(data$ReferenceFraction), "Every group must have a reference fraction.")
stop_if(all(data$MatchClass %in% names(match_linetypes)), "Unknown reference match class.")
stop_if(all(data$Category %in% category_order), "Unknown publication category.")

format_n <- function(x) ifelse(is.na(x), "NR", as.character(x))

data[, DeltaPP := KlaFraction - ReferenceFraction]
data[, RowLabel := paste0(DisplayGroup, " · ", PXD)]
data[, NLabel := paste0("n=", format_n(KlaN), "/", format_n(ReferenceN))]
data[, DeltaLabel := sprintf("%+.1f pp", DeltaPP)]
data[, ReferenceLabel := sprintf("%d/%d (%.1f%%)", ReferenceDdr, ReferenceTotal, ReferenceFraction)]
data[, KlaLabel := sprintf("%d/%d (%.1f%%)", KlaDdr, KlaTotal, KlaFraction)]

stop_if(identical(as.integer(table(factor(data$Category, levels = category_order))), expected_category_counts),
  "The selected publication category counts do not match the configured category contract."
)

make_category_panel <- function(category, show_x_title = FALSE) {
  panel_data <- data[Category == category]
  setorder(panel_data, RowOrder)
  panel_data[, PlotRow := factor(as.character(RowOrder), levels = rev(as.character(RowOrder)))]
  reference_points <- panel_data[, .(PlotRow, ReferenceFraction, ReferenceLabel)]
  kla_points <- panel_data[, .(PlotRow, KlaFraction, KlaLabel)]

  ggplot(panel_data, aes(y = PlotRow)) +
    geom_segment(
      aes(x = ReferenceFraction, xend = KlaFraction, yend = PlotRow, linetype = MatchClass),
      linewidth = 0.75,
      colour = line_colour,
      lineend = "round"
    ) +
    geom_point(
      aes(x = ReferenceFraction, fill = "Whole proteome"),
      shape = 21,
      size = 3.8,
      stroke = 0.65,
      colour = "white"
    ) +
    geom_point(
      aes(x = KlaFraction, fill = "Kla"),
      shape = 21,
      size = 3.8,
      stroke = 0.65,
      colour = "white"
    ) +
    geom_text(
      data = reference_points,
      aes(x = ReferenceFraction, label = ReferenceLabel),
      hjust = 1.08,
      size = 3.0,
      position = position_nudge(y = 0.15),
      family = publication_font,
      colour = whole_proteome_colour
    ) +
    geom_text(
      data = kla_points,
      aes(x = KlaFraction, label = KlaLabel),
      hjust = -0.08,
      size = 3.0,
      position = position_nudge(y = -0.15),
      family = publication_font,
      colour = "#C36E12"
    ) +
    geom_text(
      aes(x = 16.0, label = NLabel),
      hjust = 0,
      size = 3.35,
      family = publication_font,
      colour = charcoal
    ) +
    geom_text(
      aes(x = 20.5, label = DeltaLabel),
      hjust = 0,
      size = 3.35,
      family = publication_font,
      colour = charcoal,
      fontface = "bold"
    ) +
    geom_vline(xintercept = 15, colour = "#C8CED6", linewidth = 0.45) +
    scale_fill_manual(
      values = c("Whole proteome" = whole_proteome_colour, "Kla" = kla_colour),
      breaks = c("Whole proteome", "Kla")
    ) +
    scale_linetype_manual(values = match_linetypes, guide = "none") +
    scale_x_continuous(
      limits = c(0, 27),
      breaks = c(0, 5, 10, 15),
      expand = expansion(mult = c(0, 0))
    ) +
    scale_y_discrete(labels = setNames(panel_data$RowLabel, as.character(panel_data$RowOrder))) +
    labs(
      title = category_labels[[category]],
      x = if (show_x_title) "GO-DDR annotated protein fraction (%)" else NULL,
      y = NULL
    ) +
    theme_minimal(base_family = publication_font, base_size = 12) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(colour = grid_colour, linewidth = 0.45),
      axis.text.y = element_text(size = 9.6, colour = charcoal, lineheight = 0.92),
      axis.text.x = element_text(size = 12.5, colour = charcoal),
      axis.title.x = element_text(size = 16.5, face = "bold", colour = charcoal, margin = margin(t = 10)),
      plot.title = element_text(size = 17.5, face = "bold", colour = charcoal, hjust = 0, margin = margin(b = 6)),
      plot.margin = margin(8, 10, if (show_x_title) 8 else 3, 10),
      plot.background = element_rect(fill = "white", colour = NA)
    )
}

final_plot <- patchwork::wrap_plots(
  make_category_panel("normal_tissue"),
  make_category_panel("cancer_tissue"),
  make_category_panel("cancer_cells"),
  make_category_panel("normal_cells", show_x_title = TRUE),
  design = "AB\nCC\nDD",
  widths = c(1.12, 0.88),
  heights = c(1.35, 0.96, 1.25),
  guides = "collect"
) +
  plot_annotation(
    title = paste0("DDR-annotated protein fractions across the ", expected_group_count, " publication groups"),
    subtitle = paste(
      "Each row compares the whole-proteome reference (blue) with the linked Kla measurement (orange).",
      "The two tissue categories share the upper block; the two cell-line categories remain separate below."
    ),
    caption = paste(
      "Solid, dashed and dotted connectors retain the original same-study, external and caveat/process reference classes.",
      "n gives the available Kla/reference sample counts; delta is Kla minus reference in percentage points."
    ),
    theme = theme(
      plot.title = element_text(
        family = publication_font,
        size = 22,
        face = "bold",
        colour = charcoal,
        hjust = 0,
        margin = margin(b = 5)
      ),
      plot.subtitle = element_text(
        family = publication_font,
        size = 13,
        colour = "#4B5563",
        hjust = 0,
        lineheight = 1.08,
        margin = margin(b = 8)
      ),
      plot.caption = element_text(
        family = publication_font,
        size = 10.8,
        colour = "#4B5563",
        hjust = 0,
        lineheight = 1.05,
        margin = margin(t = 8)
      ),
      plot.margin = margin(12, 14, 12, 14)
    )
  )

final_plot <- final_plot & theme(
  legend.position = "bottom",
  legend.text = element_text(size = 13.5, colour = charcoal),
  legend.title = element_blank(),
  legend.key.height = grid::unit(0.42, "cm"),
  legend.key.width = grid::unit(0.85, "cm")
)

plot_grob <- patchwork::patchworkGrob(final_plot)
stop_if(length(plot_grob$grobs) > 0L, "Three-row candidate plot could not be assembled.")

stem <- file.path(output_dir, "Figure_1_candidate_three_row_tissue_block")
ggsave(
  paste0(stem, ".png"), plot_grob,
  width = 20.5,
  height = 21.5,
  dpi = 300,
  bg = "white",
  device = ragg::agg_png
)
ggsave(
  paste0(stem, ".pdf"), plot_grob,
  width = 20.5,
  height = 21.5,
  bg = "white",
  device = cairo_pdf
)
message("Wrote candidate three-row Figure 1: ", stem, ".png/.pdf")
