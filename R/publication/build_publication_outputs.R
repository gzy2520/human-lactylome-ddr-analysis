#!/usr/bin/env Rscript

# Rebuild only the figures reported in DDR_Kla_manuscript_V3.docx from the
# frozen, publication-scoped inputs.  All biological joins use BaseAccession;
# gene symbols are retained solely as display annotations in supplementary data.

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(ggVennDiagram)
  library(readxl)
  library(tidyr)
})

set.seed(25)

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args)) {
  normalizePath(args[[1]], mustWork = TRUE)
} else {
  normalizePath(".", mustWork = TRUE)
}
input_dir <- file.path(project_root, "data", "publication_input")
figure_dir <- file.path(project_root, "results", "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
publication_font <- "Arial Unicode MS"

assert <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

assert(
  publication_font %in% systemfonts::system_fonts()$family,
  paste("Required publication font is unavailable:", publication_font)
)

input_path <- function(filename) file.path(input_dir, filename)
md5_file <- function(path) digest::digest(file = path, algo = "md5", serialize = FALSE)

required_inputs <- c(
  "group_summary_30.csv",
  "regulator_kla_percentiles_30.csv",
  "regulator_reference_percentiles_30.csv",
  "pathway_display.csv",
  "Supplementary_Table_S4_Pathway_Protein_Ranking.xlsx",
  "Supplementary_Table_S5_Lactylation_Regulators.xlsx",
  "venn_all_kla.csv",
  "venn_kla_ddr.csv",
  "venn_reference.csv",
  "venn_reference_ddr.csv"
)
missing_inputs <- required_inputs[!file.exists(input_path(required_inputs))]
assert(!length(missing_inputs), paste("Missing publication input:", paste(missing_inputs, collapse = ", ")))

input_manifest <- fread(input_path("INPUT_MANIFEST.csv"))
assert(identical(names(input_manifest), c("File", "Bytes", "MD5")), "The frozen input manifest has an unexpected schema.")
actual_input_files <- list.files(input_dir, recursive = TRUE, full.names = FALSE, no.. = TRUE)
assert(
  setequal(actual_input_files, c("INPUT_MANIFEST.csv", input_manifest$File)),
  "The frozen publication input contains an unmanifested or missing file."
)
for (index in seq_len(nrow(input_manifest))) {
  filename <- input_manifest$File[[index]]
  path <- input_path(filename)
  assert(file.exists(path), paste("Frozen input is missing:", filename))
  assert(file.info(path)$size == input_manifest$Bytes[[index]], paste("Frozen input byte count changed:", filename))
  assert(md5_file(path) == input_manifest$MD5[[index]], paste("Frozen input checksum changed:", filename))
}

category_order <- c("normal_tissue", "cancer_tissue", "cancer_cells", "normal_cells")
category_labels <- c(
  normal_tissue = "non-tumor tissues",
  cancer_tissue = "tumor tissues",
  cancer_cells = "cancer cell lines",
  normal_cells = "normal cell lines"
)
category_counts <- c(
  normal_tissue = 9L,
  cancer_tissue = 2L,
  cancer_cells = 12L,
  normal_cells = 7L
)
pathway_order <- c("BER", "NER", "MMR", "FA", "HR", "AEJ", "NHEJ")
pathway_weights <- stats::setNames(seq_along(pathway_order), pathway_order)
role_order <- c("Writer", "Eraser", "Writer-Eraser", "Reader")

pathway_display <- fread(input_path("pathway_display.csv"))
assert(
  identical(pathway_display$Pathway[order(pathway_display$PathwayOrder)], pathway_order),
  "The pathway display input must define BER, NER, MMR, FA, HR, AEJ and NHEJ in manuscript order."
)
pathway_colours <- stats::setNames(pathway_display$Color, pathway_display$Pathway)
pathway_legend_labels <- c(
  stats::setNames(paste(pathway_order, "positive factor"), pathway_order),
  "suppressing" = "negative regulator",
  "unassigned" = "unassigned"
)

is_true <- function(x) as.character(x) %in% c("TRUE", "True", "true", "1", "T")

groups <- fread(input_path("group_summary_30.csv")) |>
  as_tibble() |>
  arrange(RowOrder) |>
  mutate(
    GroupKey = paste(PXD, SampleGroup, sep = "__"),
    Category = factor(Category, levels = category_order),
    CategoryLabel = unname(category_labels[as.character(Category)]),
    PlotLabel = ifelse(
      !is.na(KlaLabelEn) & nzchar(KlaLabelEn),
      KlaLabelEn,
      SampleGroup
    )
  )
assert(nrow(groups) == 30L, "The frozen publication input must start with exactly 30 groups.")
assert(!anyDuplicated(groups$GroupKey), "Publication groups must be unique PXD/sample-group pairs.")
observed_category_counts <- table(factor(as.character(groups$Category), levels = category_order))
assert(
  identical(as.integer(observed_category_counts), as.integer(category_counts[category_order])),
  "Publication group counts must be 9/2/12/7 in manuscript category order."
)

save_figure <- function(plot, stem, width, height) {
  ggsave(file.path(figure_dir, paste0(stem, ".png")), plot,
    width = width, height = height, dpi = 300, bg = "white", device = ragg::agg_png
  )
  ggsave(file.path(figure_dir, paste0(stem, ".pdf")), plot,
    width = width, height = height, device = cairo_pdf, bg = "white"
  )
}

# Figure 1: independent Kla and whole-proteome DDR fractions for the 30 groups.
fraction_data <- bind_rows(
  groups |>
    transmute(Category, CategoryLabel, RowOrder, PlotLabel,
      Measurement = "Kla proteome", Fraction = KlaDdrFraction),
  groups |>
    transmute(Category, CategoryLabel, RowOrder, PlotLabel,
      Measurement = "Whole proteome", Fraction = ReferenceDdrFraction)
) |>
  mutate(
    PlotLabel = factor(PlotLabel, levels = rev(unique(groups$PlotLabel))),
    Measurement = factor(Measurement, levels = c("Kla proteome", "Whole proteome"))
  )

fraction_plot <- ggplot(fraction_data, aes(x = PlotLabel, y = Fraction, fill = Measurement)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.72) +
  coord_flip() +
  facet_grid(CategoryLabel ~ ., scales = "free_y", space = "free_y") +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1), expand = expansion(mult = c(0, 0.04))) +
  scale_fill_manual(values = c("Kla proteome" = "#C44E52", "Whole proteome" = "#4C72B0")) +
  labs(x = NULL, y = "DDR protein fraction", fill = NULL) +
  theme_classic(base_size = 10, base_family = publication_font) +
  theme(
    strip.background = element_rect(fill = "#F0F0F0", colour = NA),
    strip.text = element_text(face = "bold"),
    axis.text.y = element_text(size = 7),
    legend.position = "top",
    panel.spacing.y = grid::unit(0.25, "lines")
  )
save_figure(fraction_plot, "Figure_1_DDR_fraction", 10, 12)

build_role_map <- function() {
  read_excel(input_path("Supplementary_Table_S5_Lactylation_Regulators.xlsx"), sheet = "Regulator_Annotations") |>
    transmute(
      Role = trimws(as.character(Role)),
      GeneSymbol = trimws(as.character(GeneSymbol)),
      BaseAccession = trimws(as.character(BaseAccession))
    ) |>
    filter(Role %in% role_order, !is.na(BaseAccession), nzchar(BaseAccession)) |>
    distinct(Role, BaseAccession, .keep_all = TRUE) |>
    mutate(
      DisplayName = ifelse(BaseAccession == "Q92830", "GCN5 (KAT2A)", GeneSymbol),
      Role = factor(Role, levels = role_order)
    )
}

role_map <- build_role_map()
assert(!anyDuplicated(role_map[c("Role", "BaseAccession")]), "Each regulator role/accession pair must be unique.")

draw_percentile_heatmap <- function(data, value_column, stem) {
  values <- data |>
    mutate(
      GroupKey = paste(PXD, SampleGroup, sep = "__"),
      Value = .data[[value_column]]
    ) |>
    inner_join(groups |> select(GroupKey, RowOrder, PlotLabel, Category), by = "GroupKey", relationship = "many-to-one") |>
    inner_join(role_map |> select(Role, BaseAccession, DisplayName), by = c("RegulatorBaseAccession" = "BaseAccession"), relationship = "many-to-many") |>
    filter(!is.na(Value)) |>
    distinct(Role, DisplayName, PlotLabel, .keep_all = TRUE) |>
    mutate(
      PlotLabel = factor(PlotLabel, levels = groups$PlotLabel),
      Role = factor(Role, levels = role_order)
    )
  assert(nrow(values) > 0L, paste("No values were available for", stem))
  plot <- ggplot(values, aes(x = PlotLabel, y = DisplayName, fill = Value)) +
    geom_tile(colour = "white", linewidth = 0.1) +
    facet_grid(Role ~ ., scales = "free_y", space = "free_y") +
    scale_fill_gradient(low = "#F7FBFF", high = "#08519C", limits = c(0, 100), name = "Percentile") +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 9, base_family = publication_font) +
    theme(
      panel.grid = element_blank(),
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 60, hjust = 1, size = 6),
      axis.text.y = element_text(size = 7),
      panel.spacing.y = grid::unit(0.3, "lines")
    )
  save_figure(plot, stem, 14, 11)
}

kla_percentiles <- fread(input_path("regulator_kla_percentiles_30.csv")) |>
  as_tibble() |>
  transmute(PXD, SampleGroup, RegulatorBaseAccession, RelativeKlaPercentile)
draw_percentile_heatmap(kla_percentiles, "RelativeKlaPercentile", "Figure_3b_Kla_regulator_percentiles")

reference_percentiles <- fread(input_path("regulator_reference_percentiles_30.csv")) |>
  as_tibble() |>
  transmute(PXD, SampleGroup, RegulatorBaseAccession, WholeProteomeRelativePercentile)
draw_percentile_heatmap(reference_percentiles, "WholeProteomeRelativePercentile", "Figure_3a_reference_regulator_percentiles")

venn_colours <- c("#0072B2", "#E69F00", "#CC79A7", "#009E73")

draw_exact_venn <- function(filename, stem) {
  membership <- fread(input_path(filename)) |>
    as_tibble()
  membership_columns <- paste0("In_", category_order)
  assert(all(membership_columns %in% names(membership)), paste("Invalid Venn input:", filename))
  sets <- lapply(category_order, function(category) {
    membership$BaseAccession[is_true(membership[[paste0("In_", category)]])]
  })
  names(sets) <- unname(category_labels[category_order])
  venn_plot <- ggVennDiagram(
    sets,
    label = "count",
    label_geom = "text",
    label_size = 5,
    set_color = venn_colours,
    set_size = 1,
    edge_size = 0.7
  ) +
    scale_fill_gradient(low = "#F7F7F7", high = "#BDBDBD", guide = "none") +
    theme_void(base_size = 10, base_family = publication_font) +
    theme(plot.margin = margin(25, 85, 25, 85))
  save_figure(venn_plot, stem, 11, 8.5)
}

draw_exact_venn("venn_reference_ddr.csv", "Figure_2a_whole_proteome_DDR_Venn")
draw_exact_venn("venn_kla_ddr.csv", "Figure_2b_Kla_DDR_Venn")
draw_exact_venn("venn_reference.csv", "Supplementary_Figure_S1a_whole_proteome_Venn")
draw_exact_venn("venn_all_kla.csv", "Supplementary_Figure_S1b_Kla_proteome_Venn")

# Four signed pathway matrices reported in the manuscript. The frozen S4
# ranking workbook is the data source; it is never recalculated by this code.
kla_ddr_membership <- fread(input_path("venn_kla_ddr.csv"))
assert(nrow(kla_ddr_membership) == 399L && uniqueN(kla_ddr_membership$BaseAccession) == 399L,
  "The final Kla-DDR union must contain 399 BaseAccessions."
)
matrix_specs <- list(
  list(key = "normal_tissue", sheet = "NonTumorTissues", label = "non-tumor tissues", expected = 183L),
  list(key = "cancer_tissue", sheet = "TumorTissues", label = "tumor tissues", expected = 178L),
  list(key = "cancer_cells", sheet = "CancerCellLines", label = "cancer cell lines", expected = 381L),
  list(key = "normal_cells", sheet = "NormalCellLines", label = "normal cell lines", expected = 292L)
)

pathway_panels <- list()
for (spec in matrix_specs) {
  panel <- as.data.table(read_excel(input_path("Supplementary_Table_S4_Pathway_Protein_Ranking.xlsx"), sheet = spec$sheet))
  panel <- panel[!is.na(BaseAccession) & nzchar(trimws(BaseAccession))]
  panel[, BaseAccession := trimws(as.character(BaseAccession))]
  assert(all(c("BaseAccession", "SignedScore", pathway_order) %in% names(panel)), paste("Frozen S4 is missing required columns in", spec$sheet))
  assert(!anyDuplicated(panel$BaseAccession), paste("Frozen S4 has duplicated BaseAccessions in", spec$sheet))
  assert(nrow(panel) == spec$expected, paste("Unexpected pathway-matrix size for", spec$label))
  expected_ids <- kla_ddr_membership$BaseAccession[is_true(kla_ddr_membership[[paste0("In_", spec$key)]])]
  assert(setequal(panel$BaseAccession, expected_ids), paste("Frozen S4 membership does not match", spec$label))
  score_matrix <- as.matrix(panel[, ..pathway_order])
  storage.mode(score_matrix) <- "numeric"
  assert(all(score_matrix %in% c(-1, 0, 1)), paste("Frozen S4 pathway states must be -1, 0 or +1 in", spec$sheet))
  assert(identical(as.numeric(panel$SignedScore), as.numeric(score_matrix %*% pathway_weights[pathway_order])), paste("Frozen S4 signed scores changed in", spec$sheet))
  assert(identical(order(panel$SignedScore, panel$BaseAccession), seq_len(nrow(panel))), paste("Frozen S4 order changed in", spec$sheet))
  long <- melt(
    panel[, c("BaseAccession", pathway_order), with = FALSE],
    id.vars = "BaseAccession", variable.name = "Pathway", value.name = "State"
  )
  long[, `:=`(
    Pathway = factor(Pathway, levels = pathway_order),
    Rank = match(BaseAccession, panel$BaseAccession),
    FillClass = fifelse(State == 1, as.character(Pathway), fifelse(State == -1, "suppressing", "unassigned"))
  )]
  expected_summary_states <- CJ(Pathway = pathway_order, State = c(-1, 1))
  expected_summary_states[, Pathway := factor(Pathway, levels = pathway_order)]
  pathway_summary <- long[, .(ProteinCount = .N), by = .(Pathway, State, FillClass)] |>
    merge(expected_summary_states, all.y = TRUE, by = c("Pathway", "State")) |>
    as.data.table()
  pathway_summary[is.na(ProteinCount), ProteinCount := 0L]
  pathway_summary[, `:=`(
    Panel = spec$label,
    Direction = fifelse(State == 1, "promoting", "suppressing"),
    SignedFraction = fifelse(State == 1, ProteinCount / nrow(panel), -ProteinCount / nrow(panel)),
    FillClass = fifelse(State == 1, as.character(Pathway), "suppressing")
  )]
  pathway_panels[[spec$key]] <- list(
    label = spec$label,
    matrix = long,
    summary = pathway_summary
  )
}

draw_pathway_matrices <- function(keys, stem) {
  matrix_data <- rbindlist(lapply(keys, function(key) {
    data <- copy(pathway_panels[[key]]$matrix)
    data[, Panel := pathway_panels[[key]]$label]
    data
  }))
  matrix_data[, Panel := factor(Panel, levels = vapply(keys, function(key) pathway_panels[[key]]$label, character(1)))]
  matrix_plot <- ggplot(matrix_data, aes(x = Rank, y = Pathway, fill = FillClass)) +
    geom_tile(colour = NA) +
    scale_fill_manual(
      values = c(pathway_colours, "suppressing" = "#303030", "unassigned" = "#E6E6E6"),
      breaks = c(pathway_order, "suppressing", "unassigned"),
      labels = pathway_legend_labels,
      name = NULL
    ) +
    facet_wrap(~Panel, nrow = 1, scales = "free_x") +
    scale_x_continuous(expand = c(0, 0), breaks = NULL) +
    labs(x = "Kla-DDR proteins ranked by signed pathway score", y = NULL) +
    theme_minimal(base_size = 10, base_family = publication_font) +
    theme(
      panel.grid = element_blank(),
      strip.text = element_text(face = "bold"),
      legend.position = "top"
    )
  save_figure(matrix_plot, stem, 14, 3.6)
}

draw_pathway_summary <- function(key, stem) {
  summary_data <- pathway_panels[[key]]$summary
  summary_plot <- ggplot(summary_data, aes(x = Pathway, y = SignedFraction, fill = FillClass)) +
    geom_col(width = 0.72) +
    geom_hline(yintercept = 0, colour = "#666666", linewidth = 0.3) +
    scale_fill_manual(
      values = c(pathway_colours, "suppressing" = "#303030"),
      breaks = c(pathway_order, "suppressing"),
      labels = pathway_legend_labels[c(pathway_order, "suppressing")],
      name = NULL
    ) +
    scale_y_continuous(labels = scales::label_percent(accuracy = 1)) +
    labs(x = NULL, y = "Fraction of Kla-DDR proteins", title = pathway_panels[[key]]$label) +
    theme_classic(base_size = 10, base_family = publication_font) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "top"
    )
  save_figure(summary_plot, stem, 7, 4.5)
}

draw_pathway_matrices(c("cancer_tissue", "normal_tissue"), "Figure_2c_DDR_pathway_matrices_tissues")
draw_pathway_summary("cancer_tissue", "Figure_2d_DDR_pathway_summary_tumor_tissue")
draw_pathway_summary("normal_tissue", "Figure_2e_DDR_pathway_summary_non_tumor_tissue")
draw_pathway_matrices(c("cancer_cells", "normal_cells"), "Supplementary_Figure_S2a_DDR_pathway_matrices_cell_lines")
draw_pathway_summary("cancer_cells", "Supplementary_Figure_S2b_DDR_pathway_summary_cancer_cell_lines")
draw_pathway_summary("normal_cells", "Supplementary_Figure_S2c_DDR_pathway_summary_normal_cell_lines")

message("PASS: rebuilt only manuscript figures from the frozen 30-group publication input.")
