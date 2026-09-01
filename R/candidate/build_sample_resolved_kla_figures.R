#!/usr/bin/env Rscript

# Candidate-only rendering of the committed sample-resolved Kla–DDR inputs.
# This does not alter the frozen publication figures or workflow.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
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

candidate_dir <- file.path(project_root, "data", "candidate")
output_dir <- file.path(project_root, "results", "candidate", "sample_resolved")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

publication_font <- "Arial Unicode MS"
whole_proteome_colour <- "#4E79A7"
kla_colour <- "#F28E2B"
charcoal <- "#2F3437"
grid_colour <- "#D9DDE3"
pathway_order <- c("BER", "NER", "MMR", "FA", "HR", "NHEJ", "AEJ")
condition_order <- c(
  "Adjacent skin",
  "Hypertrophic scar",
  "Human hippocampus",
  "Benign prostatic hyperplasia",
  "Prostate cancer",
  "RKO GSK3B WT",
  "RKO GSK3B KO",
  "HK-2 control",
  "HK-2 mannitol"
)
study_labels <- c(
  PXD046800 = "PXD046800 · human scar / adjacent skin · Tier A",
  PXD050470 = "PXD050470 · human hippocampus · Tier A",
  PXD066054 = "PXD066054 · prostate tissue · Tier A",
  PXD078013 = "PXD078013 · RKO cells · Tier B",
  PXD078736 = "PXD078736 · HK-2 cells · Tier B"
)

summary_path <- file.path(candidate_dir, "sample_resolved_kla_summary.csv")
paired_path <- file.path(candidate_dir, "sample_resolved_matched_modalities.csv")
pathway_path <- file.path(candidate_dir, "sample_resolved_pathway_profile.csv")
reconciliation_path <- file.path(candidate_dir, "sample_resolved_source_reconciliation.csv")
for (item_path in c(summary_path, paired_path, pathway_path, reconciliation_path)) {
  stop_if(file.exists(item_path), paste0("Missing candidate input: ", item_path))
}

sample_summary <- fread(summary_path)
paired <- fread(paired_path)
pathway_profile <- fread(pathway_path)
reconciliation <- fread(reconciliation_path)

stop_if(nrow(sample_summary) == 31L, "Sample-resolved summary must contain 31 source samples.")
stop_if(nrow(paired) == 21L, "Sample-resolved matched table must contain 21 same-sample pairs.")
stop_if(nrow(pathway_profile) == 31L * length(pathway_order), "Sample-resolved pathway profile is incomplete.")
stop_if(
  all(reconciliation[, KlaProteinCountMatchesFrozen & KlaDdrProteinCountMatchesFrozen & WholeProteomeProteinCountMatchesFrozen & WholeProteomeDdrProteinCountMatchesFrozen]),
  "Source reconciliation must pass before rendering candidate figures."
)

sample_summary[, ConditionLabel := factor(ConditionLabel, levels = condition_order)]
sample_summary[, StudyLabel := factor(study_labels[PXD], levels = unname(study_labels))]
sample_summary[, TierLabel := factor(
  EvidenceTier,
  levels = c("A_same_source_same_sample", "B_kla_sample_only"),
  labels = c(
    "Tier A — same-source sample identifiers in both modalities",
    "Tier B — Kla sample identifiers; external references excluded"
  )
)]
sample_summary[, GroupKey := paste(PXD, ConditionLabel, sep = "__")]
group_levels <- unique(sample_summary[order(PXD, ConditionLabel), GroupKey])
sample_summary[, GroupKey := factor(GroupKey, levels = rev(group_levels))]
sample_summary[, GroupLabel := paste0(as.character(ConditionLabel), " · ", PXD)]

group_statistics <- sample_summary[, .(
  SampleCount = .N,
  Minimum = min(KlaDdrFraction * 100),
  Mean = mean(KlaDdrFraction * 100),
  Maximum = max(KlaDdrFraction * 100)
), by = .(GroupKey, GroupLabel, TierLabel)]
group_statistics[, NLabel := paste0("n=", SampleCount)]
display_maximum <- max(c(sample_summary$KlaDdrFraction * 100, group_statistics$Maximum))
sample_x_limit <- ceiling((display_maximum + 1.5) / 2) * 2

sample_fraction_plot <- ggplot(sample_summary, aes(y = GroupKey, x = KlaDdrFraction * 100)) +
  geom_segment(
    data = group_statistics,
    aes(y = GroupKey, yend = GroupKey, x = Minimum, xend = Maximum),
    inherit.aes = FALSE,
    linewidth = 1.0,
    colour = "#AAB3BD",
    lineend = "round"
  ) +
  geom_point(
    position = position_jitter(height = 0.13, width = 0, seed = 25),
    shape = 21,
    size = 4.0,
    stroke = 0.75,
    fill = kla_colour,
    colour = "white"
  ) +
  geom_point(
    data = group_statistics,
    aes(y = GroupKey, x = Mean),
    inherit.aes = FALSE,
    shape = 23,
    size = 4.7,
    stroke = 0.75,
    fill = kla_colour,
    colour = charcoal
  ) +
  geom_text(
    data = group_statistics,
    aes(y = GroupKey, x = sample_x_limit - 0.1, label = NLabel),
    inherit.aes = FALSE,
    hjust = 1,
    family = publication_font,
    size = 4.2,
    colour = charcoal
  ) +
  facet_grid(TierLabel ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_y_discrete(labels = setNames(group_statistics$GroupLabel, group_statistics$GroupKey)) +
  scale_x_continuous(
    limits = c(0, sample_x_limit),
    breaks = pretty_breaks(n = 6),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    title = "Sample-resolved Kla–DDR coverage",
    subtitle = "Each circle is one deposited source sample/run; horizontal lines show the observed range and diamonds the arithmetic mean. No cross-study statistical pooling is used.",
    x = "Kla proteins annotated to DDR (%)",
    y = NULL,
    caption = "Tier B source samples are shown only as Kla measurements. Their external whole-proteome references are deliberately not displayed as paired observations."
  ) +
  theme_minimal(base_family = publication_font, base_size = 14) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(colour = grid_colour, linewidth = 0.45),
    axis.text.y = element_text(size = 12.5, colour = charcoal),
    axis.text.x = element_text(size = 13, colour = charcoal),
    axis.title.x = element_text(size = 17, face = "bold", margin = margin(t = 10)),
    plot.title = element_text(size = 22, face = "bold", colour = charcoal, margin = margin(b = 5)),
    plot.subtitle = element_text(size = 13, colour = "#4B5563", lineheight = 1.05, margin = margin(b = 10)),
    plot.caption = element_text(size = 11.5, colour = "#4B5563", hjust = 0, lineheight = 1.05, margin = margin(t = 10)),
    strip.placement = "outside",
    strip.background = element_rect(fill = "#F2F2F2", colour = NA),
    strip.text.y.left = element_text(size = 13.5, face = "bold", colour = charcoal, angle = 90),
    panel.spacing.y = grid::unit(0.8, "lines"),
    plot.margin = margin(12, 15, 12, 15)
  )

paired[, ConditionLabel := factor(ConditionLabel, levels = condition_order)]
paired[, StudyLabel := factor(study_labels[PXD], levels = unname(study_labels[1:3]))]
paired[, SampleLabel := paste0(SampleID, "  ·  ", as.character(ConditionLabel))]
paired[, PairRow := factor(SampleKey, levels = rev(unique(paired[order(PXD, ConditionLabel, SampleID), SampleKey])))]
paired[, DeltaLabel := sprintf("%+.1f pp", DeltaPercentagePoints)]
paired_x_limit <- ceiling((max(c(paired$KlaDdrFraction, paired$WholeProteomeDdrFraction)) * 100 + 1.8) / 2) * 2

paired_plot <- ggplot(paired, aes(y = PairRow)) +
  geom_segment(
    aes(
      x = WholeProteomeDdrFraction * 100,
      xend = KlaDdrFraction * 100,
      yend = PairRow
    ),
    linewidth = 0.95,
    colour = "#9AA5B1",
    lineend = "round"
  ) +
  geom_point(
    aes(x = WholeProteomeDdrFraction * 100, fill = "Whole proteome"),
    shape = 21,
    size = 3.9,
    stroke = 0.65,
    colour = "white"
  ) +
  geom_point(
    aes(x = KlaDdrFraction * 100, fill = "Kla"),
    shape = 21,
    size = 3.9,
    stroke = 0.65,
    colour = "white"
  ) +
  geom_text(
    aes(x = paired_x_limit - 0.1, label = DeltaLabel),
    hjust = 1,
    family = publication_font,
    size = 3.9,
    colour = charcoal
  ) +
  facet_grid(StudyLabel ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_manual(values = c("Whole proteome" = whole_proteome_colour, "Kla" = kla_colour)) +
  scale_y_discrete(labels = setNames(paired$SampleLabel, paired$SampleKey)) +
  scale_x_continuous(
    limits = c(0, paired_x_limit),
    breaks = pretty_breaks(n = 6),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    title = "Within-sample comparison where both modalities were deposited",
    subtitle = "Only same-source samples with matching identifiers are connected. The label at right is Kla minus whole-proteome coverage.",
    x = "DDR-annotated protein fraction (%)",
    y = NULL,
    fill = NULL,
    caption = "Both values follow the corresponding frozen Figure 1 membership tables. They are descriptive sample-level coverage measurements, not independent protein-level statistical tests."
  ) +
  theme_minimal(base_family = publication_font, base_size = 14) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(colour = grid_colour, linewidth = 0.45),
    axis.text.y = element_text(size = 11.6, colour = charcoal),
    axis.text.x = element_text(size = 13, colour = charcoal),
    axis.title.x = element_text(size = 17, face = "bold", margin = margin(t = 10)),
    plot.title = element_text(size = 22, face = "bold", colour = charcoal, margin = margin(b = 5)),
    plot.subtitle = element_text(size = 13, colour = "#4B5563", lineheight = 1.05, margin = margin(b = 10)),
    plot.caption = element_text(size = 11.5, colour = "#4B5563", hjust = 0, lineheight = 1.05, margin = margin(t = 10)),
    legend.position = "bottom",
    legend.text = element_text(size = 13),
    strip.placement = "outside",
    strip.background = element_rect(fill = "#F2F2F2", colour = NA),
    strip.text.y.left = element_text(size = 13.5, face = "bold", colour = charcoal, angle = 90),
    panel.spacing.y = grid::unit(0.8, "lines"),
    plot.margin = margin(12, 15, 12, 15)
  )

sample_coverage_figure <- sample_fraction_plot / paired_plot +
  plot_layout(heights = c(0.9, 1.25)) +
  plot_annotation(tag_levels = "A")

pathway_profile[, ConditionLabel := factor(ConditionLabel, levels = condition_order)]
pathway_profile[, StudyLabel := factor(study_labels[PXD], levels = unname(study_labels))]
pathway_profile[, Pathway := factor(Pathway, levels = pathway_order)]
pathway_profile[, SampleLabel := paste0(SampleID, "  ·  ", as.character(ConditionLabel))]
pathway_profile[, SampleRow := factor(SampleKey, levels = rev(unique(pathway_profile[order(PXD, ConditionLabel, SampleID), SampleKey])))]
pathway_profile[, PositiveFraction := PositiveProteinCount / KlaDdrProteinCount]
pathway_profile[, NegativeFraction := NegativeProteinCount / KlaDdrProteinCount]

pathway_figure <- ggplot(pathway_profile, aes(x = Pathway, y = SampleRow)) +
  geom_tile(
    aes(fill = PositiveFraction),
    colour = "white",
    linewidth = 0.75,
    width = 0.94,
    height = 0.90
  ) +
  geom_point(
    data = pathway_profile[NegativeFraction > 0],
    aes(size = NegativeFraction),
    shape = 22,
    fill = charcoal,
    colour = "white",
    stroke = 0.25
  ) +
  facet_grid(StudyLabel ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_gradientn(
    colours = c("#FFFFFF", "#E8F3FA", "#9ECAE1", "#4292C6", "#08519C"),
    values = rescale(c(0, 0.08, 0.16, 0.24, 0.32)),
    limits = c(0, 0.32),
    oob = squish,
    name = "Positive pathway\nmembers / Kla-DDR proteins"
  ) +
  scale_size_continuous(
    range = c(2.5, 5.2),
    limits = c(0, max(pathway_profile$NegativeFraction)),
    name = "Negative\nmember fraction"
  ) +
  scale_y_discrete(labels = setNames(pathway_profile$SampleLabel, pathway_profile$SampleKey)) +
  labs(
    title = "Sample-resolved Kla–DDR pathway repertoire",
    subtitle = "Tile intensity shows positively scored pathway proteins as a fraction of each sample's Kla-DDR proteins. Charcoal squares mark pathways containing negatively scored members.",
    x = NULL,
    y = NULL,
    caption = "The seven curated pathway assignments are non-exclusive; a protein can contribute to more than one pathway. This panel uses the frozen S4 pathway scores without reclassification."
  ) +
  theme_minimal(base_family = publication_font, base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 11.2, colour = charcoal),
    axis.text.x = element_text(size = 14, face = "bold", colour = charcoal),
    plot.title = element_text(size = 22, face = "bold", colour = charcoal, margin = margin(b = 5)),
    plot.subtitle = element_text(size = 13, colour = "#4B5563", lineheight = 1.05, margin = margin(b = 10)),
    plot.caption = element_text(size = 11.5, colour = "#4B5563", hjust = 0, lineheight = 1.05, margin = margin(t = 10)),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.title = element_text(size = 11.5),
    legend.text = element_text(size = 11),
    strip.placement = "outside",
    strip.background = element_rect(fill = "#F2F2F2", colour = NA),
    strip.text.y.left = element_text(size = 13.5, face = "bold", colour = charcoal, angle = 90),
    panel.spacing.y = grid::unit(0.75, "lines"),
    plot.margin = margin(12, 15, 12, 15)
  )

coverage_stem <- file.path(output_dir, "Figure_S_candidate_sample_resolved_Kla_DDR")
pathway_stem <- file.path(output_dir, "Figure_S_candidate_sample_pathway_repertoire")
ggsave(paste0(coverage_stem, ".png"), sample_coverage_figure, width = 16.5, height = 22, dpi = 300, bg = "white", device = ragg::agg_png)
ggsave(paste0(coverage_stem, ".pdf"), sample_coverage_figure, width = 16.5, height = 22, bg = "white", device = cairo_pdf)
ggsave(paste0(pathway_stem, ".png"), pathway_figure, width = 15.5, height = 18.5, dpi = 300, bg = "white", device = ragg::agg_png)
ggsave(paste0(pathway_stem, ".pdf"), pathway_figure, width = 15.5, height = 18.5, bg = "white", device = cairo_pdf)

message("Wrote candidate sample-resolved figures to ", output_dir)
