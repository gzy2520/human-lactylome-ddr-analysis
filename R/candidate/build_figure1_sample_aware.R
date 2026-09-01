#!/usr/bin/env Rscript

# Candidate only: make a sample-aware redesign of Figure 1 without changing
# the frozen publication workflow or its approved figure outputs.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

set.seed(25)

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args)) {
  normalizePath(args[[1L]], mustWork = TRUE)
} else {
  normalizePath(".", mustWork = TRUE)
}

input_dir <- file.path(project_root, "data", "publication_input")
design_path <- file.path(project_root, "data", "candidate", "sample_design_30.csv")
output_dir <- file.path(project_root, "results", "candidate")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

publication_font <- "Arial Unicode MS"
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

stop_if <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}
format_n <- function(x) ifelse(is.na(x), "NR", as.character(x))

groups <- fread(file.path(input_dir, "group_summary_30.csv"))
design <- fread(design_path, na.strings = c("", "NA"))
stop_if(nrow(groups) == 30L, "Frozen publication input must contain 30 groups.")
stop_if(nrow(design) == 30L, "Candidate sample design must contain 30 groups.")
stop_if(!anyDuplicated(groups[, .(PXD, SampleGroup)]), "Publication groups are not unique.")
stop_if(!anyDuplicated(design[, .(PXD, SampleGroup)]), "Candidate design rows are not unique.")

groups[, GroupKey := paste(PXD, SampleGroup, sep = "__")]
design[, GroupKey := paste(PXD, SampleGroup, sep = "__")]
stop_if(setequal(groups$GroupKey, design$GroupKey), "Candidate design does not cover exactly the 30 frozen groups.")

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
    KlaSampleDesign,
    ReferenceSampleDesign,
    Aggregation,
    MatchClass
  )],
  by = c("PXD", "SampleGroup"),
  all.x = TRUE,
  sort = FALSE
)
setorder(data, RowOrder)
stop_if(!anyNA(data$ReferenceFraction), "Candidate Fig. 1 requires a reference fraction for every group.")
stop_if(all(data$MatchClass %in% names(match_linetypes)), "Candidate design contains an unknown match class.")

data[, Category := factor(Category, levels = category_order)]
data[, CategoryLabel := factor(
  Category,
  levels = category_order,
  labels = unname(category_labels[category_order])
)]
data[, MatchType := factor(MatchClass, levels = names(match_linetypes))]
data[, RowLabel := paste0(DisplayGroup, " · ", PXD)]
data[, NLabel := paste0("n=", format_n(KlaN), "/", format_n(ReferenceN))]
data[, DeltaPP := KlaFraction - ReferenceFraction]
data[, DeltaLabel := sprintf("%+.1f pp", DeltaPP)]
data[, ReferenceLabel := sprintf("%d/%d (%.1f%%)", ReferenceDdr, ReferenceTotal, ReferenceFraction)]
data[, KlaLabel := sprintf("%d/%d (%.1f%%)", KlaDdr, KlaTotal, KlaFraction)]
data[, PlotRow := factor(as.character(RowOrder), levels = rev(as.character(RowOrder)))]

reference_points <- data[, .(
  CategoryLabel,
  PlotRow,
  ReferenceFraction,
  ReferenceLabel
)]
kla_points <- data[, .(
  CategoryLabel,
  PlotRow,
  KlaFraction,
  KlaLabel
)]

candidate_plot <- ggplot(data, aes(y = PlotRow)) +
  geom_segment(
    aes(
      x = ReferenceFraction,
      xend = KlaFraction,
      yend = PlotRow,
      linetype = MatchType
    ),
    linewidth = 0.65,
    colour = "#9AA5B1"
  ) +
  geom_point(
    aes(x = ReferenceFraction),
    shape = 21,
    size = 3.1,
    stroke = 0.55,
    fill = "#4E79A7",
    colour = "white"
  ) +
  geom_point(
    aes(x = KlaFraction),
    shape = 21,
    size = 3.1,
    stroke = 0.55,
    fill = "#F28E2B",
    colour = "white"
  ) +
  geom_text(
    data = reference_points,
    aes(x = ReferenceFraction, label = ReferenceLabel),
    hjust = 1.08,
    size = 2.65,
    position = position_nudge(y = 0.13),
    family = publication_font,
    colour = "#4E79A7"
  ) +
  geom_text(
    data = kla_points,
    aes(x = KlaFraction, label = KlaLabel),
    hjust = -0.08,
    size = 2.65,
    position = position_nudge(y = -0.13),
    family = publication_font,
    colour = "#C36E12"
  ) +
  geom_text(
    aes(x = 16.0, label = NLabel),
    hjust = 0,
    size = 3.0,
    family = publication_font,
    colour = "#30343B"
  ) +
  geom_text(
    aes(x = 20.5, label = DeltaLabel),
    hjust = 0,
    size = 3.0,
    family = publication_font,
    colour = "#20252B",
    fontface = "bold"
  ) +
  geom_vline(xintercept = 15, colour = "#C8CED6", linewidth = 0.45) +
  facet_grid(
    CategoryLabel ~ .,
    scales = "free_y",
    space = "free_y",
    switch = "y"
  ) +
  scale_linetype_manual(values = match_linetypes, guide = "none") +
  scale_x_continuous(
    limits = c(0, 27),
    breaks = c(0, 5, 10, 15),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_discrete(
    labels = setNames(data$RowLabel, as.character(data$RowOrder))
  ) +
  labs(
    x = "GO-DDR annotated protein fraction (%)",
    y = NULL,
    title = "DDR annotation across the 30 publication sample groups",
    subtitle = "Blue = whole proteome reference; orange = Kla. Solid/dashed/dotted lines indicate same-study, external, and caveat/process references.\nRight: n (Kla/reference) and delta (Kla minus reference).",
    caption = "Each point is a publication sample-group summary. NR indicates that a common biological sample count was not recoverable from the processed evidence; it is not imputed. Per-group Kla sample structure and the reference relation (also encoded by line style) are tabulated in the supplementary sample design."
  ) +
  theme_minimal(base_size = 10, base_family = publication_font) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(colour = "#D9DDE3", linewidth = 0.45),
    axis.text.y = element_text(size = 8.0, colour = "#30343B", lineheight = 0.92),
    axis.text.x = element_text(size = 12, colour = "#30343B"),
    axis.title.x = element_text(size = 17.5, face = "bold", colour = "#20252B", margin = margin(t = 12)),
    plot.title = element_text(size = 18, face = "bold", colour = "#20252B", margin = margin(b = 5)),
    plot.subtitle = element_text(size = 10, colour = "#4B5563", lineheight = 1.05, margin = margin(b = 9)),
    plot.caption = element_text(size = 8.8, colour = "#4B5563", hjust = 0, margin = margin(t = 9)),
    strip.placement = "outside",
    strip.text.y.left = element_text(size = 16, face = "bold", colour = "#30343B", angle = 90),
    strip.background = element_rect(fill = "#DCEAF5", colour = NA),
    panel.spacing.y = grid::unit(0.9, "lines"),
    plot.margin = margin(10, 14, 14, 14)
  )

plot_grob <- ggplotGrob(candidate_plot)
strip_ids <- grep("^strip-l", plot_grob$layout$name)
strip_ids <- strip_ids[order(plot_grob$layout$t[strip_ids])]
stop_if(length(strip_ids) == length(category_order), "Candidate Fig. 1 category strips are incomplete.")
for (index in seq_along(strip_ids)) {
  strip_grob <- plot_grob$grobs[[strip_ids[[index]]]]$grobs[[1L]]
  background_id <- grep("^strip.background", strip_grob$childrenOrder)
  stop_if(length(background_id) == 1L, "Could not identify candidate category strip background.")
  fill <- unname(category_fills[category_order[index]])
  strip_grob$children[[background_id]]$gp$fill <- fill
  strip_grob$children[[background_id]]$gp$col <- NA
  plot_grob$grobs[[strip_ids[[index]]]]$grobs[[1L]] <- strip_grob
}

stem <- file.path(output_dir, "Figure_1_candidate_sample_aware")
ggsave(paste0(stem, ".png"), plot_grob, width = 18.5, height = 20.5, dpi = 300, bg = "white", device = ragg::agg_png)
ggsave(paste0(stem, ".pdf"), plot_grob, width = 18.5, height = 20.5, bg = "white", device = cairo_pdf)
message("Wrote candidate Figure 1: ", stem, ".png/.pdf")
