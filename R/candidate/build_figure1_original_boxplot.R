#!/usr/bin/env Rscript

# Candidate Figure 1 renderer.
#
# This preserves the reviewed PXD-axis layout: one row per publication
# sample-group dataset, four biological-category strips, and a right-hand
# statistics column.  The two source modalities are now horizontal boxplots
# dodged within the same row, with source-resolved sample observations overlaid.
# The approved publication renderer is not changed.

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
mean_colour <- "#C0392B"
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
# Keep the reviewed palette: tissue/cell-line and normal/tumor distinctions
# remain visible through the established strip colors and explicit labels.
category_fills <- c(
  normal_tissue = "#DCEAF5",
  cancer_tissue = "#FCE7D4",
  normal_cells = "#DCEAF5",
  cancer_cells = "#FCE7D4"
)
dataset_order <- c("Whole proteome", "Lactylome (Kla)")

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
  "DdrFractionPercentage", "FrozenKlaDdrFraction", "FrozenReferenceDdrFraction"
)
stop_if(all(required_values %in% names(values)), "Figure 1 sample input schema is incomplete.")
stop_if(nrow(counts) == expected_group_count,
  paste0("The sample-count record must cover ", expected_group_count, " publication groups."))
stop_if(!anyDuplicated(counts[, .(PXD, SampleGroup)]), "The sample-count record contains duplicate groups.")
stop_if(all(values$Dataset %in% dataset_order),
  "Figure 1 sample input contains an unknown dataset.")
stop_if(all(values$Category %in% category_order),
  "Figure 1 sample input contains an unknown category.")
stop_if(all(is.finite(values$DdrFractionPercentage)),
  "Figure 1 sample fractions contain non-finite values.")
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
values[, Dataset := factor(Dataset, levels = dataset_order)]
values[, GroupID := paste(PXD, SampleGroup, sep = "__")]

group_meta <- values[, .(
  RowOrder = first(RowOrder),
  PXD = first(PXD),
  SampleGroup = first(SampleGroup),
  DisplayGroup = first(DisplayGroup),
  Category = as.character(first(Category)),
  CategoryLabel = as.character(first(CategoryLabel)),
  ReferencePXD = first(ReferencePXD),
  WholeProteomeN = sum(as.character(Dataset) == dataset_order[[1L]]),
  KlaN = sum(as.character(Dataset) == dataset_order[[2L]]),
  FrozenKlaFraction = first(FrozenKlaDdrFraction),
  FrozenReferenceFraction = first(FrozenReferenceDdrFraction)
), by = GroupID]
setorder(group_meta, RowOrder, PXD, SampleGroup)
group_meta[, PlotRow := factor(GroupID, levels = rev(GroupID))]
group_meta[, CategoryLabel := factor(
  CategoryLabel,
  levels = unname(category_labels[category_order])
)]
group_meta[, RowLabel := paste0(DisplayGroup, " · ", PXD)]
group_meta[, NLabel := paste0("n=", KlaN, "/", WholeProteomeN)]
group_meta[, DeltaPP := FrozenKlaFraction - FrozenReferenceFraction]
group_meta[, DeltaLabel := sprintf("%+.1f pp", DeltaPP)]
stop_if(nrow(group_meta) == expected_group_count,
  paste0("Figure 1 must contain ", expected_group_count, " PXD/sample-group rows."))
stop_if(all(group_meta$WholeProteomeN > 0L & group_meta$KlaN > 0L),
  "Every Figure 1 PXD/sample-group row must contain both modalities.")

plot_row_levels <- levels(group_meta$PlotRow)
values[, PlotRow := factor(GroupID, levels = plot_row_levels)]
values[, N := .N, by = .(GroupID, Dataset)]

classification_audit <- group_meta[, .(
  GroupID,
  RowOrder,
  PXD,
  SampleGroup,
  DisplayGroup,
  Category,
  CategoryLabel,
  ReferencePXD,
  WholeProteomeN,
  KlaN,
  FrozenKlaFraction,
  FrozenReferenceFraction,
  DeltaPP,
  NLabel,
  DeltaLabel
)]
fwrite(
  classification_audit,
  file.path(output_dir, "figure1_original_boxplot_classification.csv"),
  na = ""
)
escc_classification <- classification_audit[
  PXD == "PXD064038" & SampleGroup == "MEC and NEC ESCC groups"
]
if (any(values$PXD == "PXD064038" & values$SampleGroup == "MEC and NEC ESCC groups")) {
  stop_if(nrow(escc_classification) == 1L &&
            escc_classification$Category == "cancer_tissue" &&
            escc_classification$CategoryLabel == "tumor tissues",
    "PXD064038 ESCC MEC/NEC is not classified as tumor tissue.")
}

row_stats <- values[, .(
  RowOrder = first(RowOrder),
  PXD = first(PXD),
  SampleGroup = first(SampleGroup),
  Category = as.character(first(Category)),
  CategoryLabel = as.character(first(CategoryLabel)),
  PlotRow = first(PlotRow),
  N = .N,
  Mean = mean(DdrFractionPercentage),
  Median = median(DdrFractionPercentage)
), by = .(GroupID, Dataset)]
row_stats[, Dataset := as.character(Dataset)]
row_stats[, CategoryLabel := factor(
  CategoryLabel,
  levels = unname(category_labels[category_order])
)]
fwrite(
  row_stats[, .(
    GroupID, RowOrder, PXD, SampleGroup, Category, CategoryLabel, Dataset,
    N, Mean, Median
  )],
  file.path(output_dir, "figure1_original_boxplot_mean_median.csv"),
  na = ""
)

source(file.path(project_root, "R", "candidate", "boxplot_significance.R"), local = TRUE)
original_anova <- compute_figure1_original_dataset_one_way_anova(
  values,
  dataset_order = dataset_order
)
original_anova[, CategoryLabel := factor(
  Category,
  levels = category_order,
  labels = unname(category_labels[category_order])
)]
original_anova <- merge(
  original_anova,
  group_meta[, .(RowOrder, PXD, SampleGroup, GroupID, PlotRow)],
  by = c("RowOrder", "PXD", "SampleGroup"),
  all.x = TRUE,
  sort = FALSE
)
fwrite(
  original_anova,
  file.path(output_dir, "figure1_original_dataset_one_way_anova.csv"),
  na = ""
)
anova_annotation <- original_anova[is.finite(QValueBH) & QValueBH < 0.05]

label_map <- stats::setNames(group_meta$RowLabel, group_meta$GroupID)
figure_plot <- ggplot(values, aes(x = DdrFractionPercentage, y = PlotRow, fill = Dataset)) +
  geom_boxplot(
    aes(group = interaction(GroupID, Dataset)),
    position = position_dodge(width = 0.62),
    width = 0.30,
    outlier.shape = NA,
    colour = charcoal,
    linewidth = 0.55,
    median.linewidth = 1.35,
    orientation = "y",
    na.rm = TRUE
  ) +
  geom_point(
    aes(group = interaction(GroupID, Dataset)),
    position = position_jitterdodge(
      jitter.width = 0,
      jitter.height = 0.055,
      dodge.width = 0.62,
      seed = 25
    ),
    shape = 21,
    size = 2.25,
    stroke = 0.45,
    colour = "white",
    na.rm = TRUE
  ) +
  geom_point(
    data = row_stats,
    aes(x = Mean, y = PlotRow, group = Dataset),
    inherit.aes = FALSE,
    position = position_dodge(width = 0.62),
    shape = 124,
    size = 6.9,
    colour = mean_colour,
    na.rm = TRUE
  ) +
  geom_text(
    data = group_meta,
    aes(x = 16.0, y = PlotRow, label = NLabel),
    inherit.aes = FALSE,
    hjust = 0,
    size = 3.0,
    family = publication_font,
    colour = muted_text
  ) +
  geom_text(
    data = group_meta,
    aes(x = 20.5, y = PlotRow, label = DeltaLabel),
    inherit.aes = FALSE,
    hjust = 0,
    size = 3.0,
    family = publication_font,
    colour = charcoal,
    fontface = "bold"
  ) +
  geom_text(
    data = anova_annotation,
    aes(x = 25.5, y = PlotRow, label = Significance),
    inherit.aes = FALSE,
    hjust = 0,
    size = 3.1,
    family = publication_font,
    colour = charcoal,
    fontface = "bold"
  ) +
  geom_vline(xintercept = 15, colour = panel_border_colour, linewidth = 0.45) +
  facet_grid(
    CategoryLabel ~ .,
    scales = "free_y",
    space = "free_y",
    switch = "y"
  ) +
  scale_fill_manual(
    values = c("Whole proteome" = whole_proteome_colour, "Lactylome (Kla)" = kla_colour),
    breaks = dataset_order,
    guide = "none"
  ) +
  scale_x_continuous(
    limits = c(0, 27),
    breaks = c(0, 5, 10, 15),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_discrete(
    labels = label_map,
    expand = expansion(mult = c(0.03, 0.03))
  ) +
  labs(
    x = "GO-DDR annotated protein fraction (%)",
    y = NULL,
    title = paste0("DDR annotation across the ", expected_group_count, " publication sample groups"),
    subtitle = paste(
      "Blue = whole-proteome observations; orange = Kla observations. Dark box line = median; red line = mean.",
      "Right: n (Kla/whole) and delta (Kla minus whole-proteome group reference). Stars: BH-adjusted one-way ANOVA within each PXD/sample-group row.",
      sep = "\n"
    ),
    caption = paste(
      "Each point is a source-resolved sample observation; one aggregate point is retained when the source only provides an aggregate profile.",
      "PXD064038 ESCC MEC/NEC is in the tumor tissues category and uses the 94 ESCC-T observations from PXD065830 as its whole-proteome reference.",
      "The n column is Kla/whole-proteome; rows with one observation in each modality have no estimable ANOVA p value.",
      sep = "\n"
    )
  ) +
  theme_minimal(base_size = 10, base_family = publication_font) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(colour = grid_colour, linewidth = 0.45),
    panel.border = element_rect(colour = panel_border_colour, fill = NA, linewidth = 0.42),
    axis.text.y = element_text(size = 8.0, colour = charcoal, lineheight = 0.92),
    axis.text.x = element_text(size = 12, colour = charcoal),
    axis.title.x = element_text(
      size = 17.5,
      face = "bold",
      colour = charcoal,
      margin = margin(t = 12)
    ),
    plot.title = element_text(size = 18, face = "bold", colour = charcoal, margin = margin(b = 5)),
    plot.subtitle = element_text(size = 10, colour = muted_text, lineheight = 1.05, margin = margin(b = 9)),
    plot.caption = element_text(size = 8.8, colour = muted_text, hjust = 0, margin = margin(t = 9)),
    strip.placement = "outside",
    strip.text.y.left = element_text(size = 16, face = "bold", colour = charcoal, angle = 90),
    strip.background = element_rect(fill = "#E7E9E7", colour = NA),
    panel.spacing.y = grid::unit(0.9, "lines"),
    plot.margin = margin(10, 14, 14, 14),
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
    strip_grob$children[[background_id]]$gp$fill <- unname(category_fills[category_order[[index]]])
    strip_grob$children[[background_id]]$gp$col <- NA
    plot_grob$grobs[[strip_ids[[index]]]]$grobs[[1L]] <- strip_grob
  }
  plot_grob
}

plot_grob <- apply_strip_fills(figure_plot)
stop_if(length(plot_grob$grobs) > 0L, "Figure 1 boxplot grob is empty.")

stem <- file.path(output_dir, "Figure_1_DDR_fraction_candidate_sample_boxplot")
figure_height <- max(20.5, nrow(group_meta) * 0.50 + 4.5)
ggsave(
  paste0(stem, ".png"),
  plot_grob,
  width = 18.5,
  height = figure_height,
  dpi = 300,
  bg = "white",
  device = ragg::agg_png
)
ggsave(
  paste0(stem, ".pdf"),
  plot_grob,
  width = 18.5,
  height = figure_height,
  bg = "white",
  device = cairo_pdf
)

message(
  "Wrote candidate original-layout Figure 1 boxplot: ",
  stem, ".png/.pdf (18.5 x ", sprintf("%.1f", figure_height), " in)."
)
