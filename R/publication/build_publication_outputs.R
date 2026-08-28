#!/usr/bin/env Rscript

# Rebuild only the figures reported in DDR_Kla_manuscript_V3.docx from the
# frozen, publication-scoped inputs.  All biological joins use BaseAccession;
# gene symbols are retained solely as display annotations in supplementary data.

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(eulerr)
  library(ggplot2)
  library(patchwork)
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
input_dir <- normalizePath(
  Sys.getenv("KLA_PUBLICATION_INPUT", unset = file.path(project_root, "data", "publication_input")),
  mustWork = TRUE
)
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

# Figure 1 retains the established publication layout: a reference bar is
# shown once for a shared whole-proteome source, immediately above the linked
# Kla bar(s).  Only the input table is new; the colour and layout contract is
# the one used before the repository was reduced to the final manuscript set.
reference_rows <- groups |>
  mutate(
    ReferenceKey = paste(
      ReferencePXD, ReferenceEvidenceFile, ReferenceProteinCount,
      ReferenceDdrProteinCount, sep = "||"
    )
  ) |>
  group_by(Category, ReferenceKey) |>
  arrange(RowOrder, .by_group = TRUE) |>
  summarise(
    RowOrder = first(RowOrder),
    PXD = first(ReferencePXD),
    DisplayLabel = paste0(first(ReferenceLabelEn), " · ", first(ReferencePXD)),
    Fraction = first(ReferenceDdrFraction),
    Ddr = first(ReferenceDdrProteinCount),
    Total = first(ReferenceProteinCount),
    BarType = "reference",
    .groups = "drop"
  )
kla_rows <- groups |>
  transmute(
    Category, RowOrder, PXD,
    DisplayLabel = paste0(KlaLabelEn, " · ", PXD),
    Fraction = KlaDdrFraction,
    Ddr = KlaDdrProteinCount,
    Total = KlaProteinCount,
    BarType = "kla"
  )
fraction_data <- bind_rows(reference_rows, kla_rows) |>
  mutate(
    BarPriority = if_else(BarType == "reference", 0L, 1L)
  ) |>
  arrange(match(as.character(Category), category_order), RowOrder, BarPriority, DisplayLabel) |>
  mutate(
    BarOrder = row_number(),
    PlotRow = factor(as.character(BarOrder), levels = rev(as.character(BarOrder))),
    CategoryLabel = factor(
      as.character(Category), levels = category_order,
      labels = unname(category_labels[category_order])
    ),
    Dataset = factor(
      BarType, levels = c("reference", "kla"),
      labels = c("Whole proteome", "Lactylome (Kla)")
    ),
    BarLabel = sprintf("%s/%s (%.1f%%)", Ddr, Total, 100 * Fraction)
  )
assert(max(fraction_data$Fraction, na.rm = TRUE) <= 0.15, "The established Figure 1 axis is 0-15%; a value exceeds that range.")
fraction_plot <- ggplot(fraction_data, aes(x = Fraction * 100, y = PlotRow, fill = Dataset)) +
  geom_col(width = 0.68, colour = "white", linewidth = 0.2) +
  geom_text(aes(label = BarLabel), hjust = -0.04, size = 2.45, family = publication_font, colour = "#30343B") +
  facet_grid(CategoryLabel ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_manual(values = c("Whole proteome" = "#4E79A7", "Lactylome (Kla)" = "#F28E2B")) +
  scale_x_continuous(limits = c(0, 15), breaks = c(0, 5, 10, 15), expand = expansion(mult = c(0, 0))) +
  scale_y_discrete(labels = setNames(fraction_data$DisplayLabel, as.character(fraction_data$BarOrder))) +
  guides(fill = guide_legend(ncol = 1, byrow = TRUE, keyheight = grid::unit(0.84, "cm"), keywidth = grid::unit(1.02, "cm"))) +
  labs(x = "GO-DDR annotated protein fraction (%)", y = NULL, fill = NULL) +
  theme_minimal(base_size = 10, base_family = publication_font) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(colour = "#D9DDE3", linewidth = 0.45),
    axis.text.y = element_text(size = 6.8, colour = "#30343B"),
    axis.text.x = element_text(size = 10, colour = "#30343B"),
    axis.title.x = element_text(size = 16.5, face = "bold", colour = "#20252B", margin = margin(t = 12)),
    strip.placement = "outside",
    strip.text.y.left = element_text(size = 15, face = "bold", colour = "#30343B", angle = 90),
    strip.background = element_rect(fill = "#DCEAF5", colour = NA),
    panel.spacing.y = grid::unit(0.72, "lines"),
    legend.position = "inside",
    legend.position.inside = c(0.965, 0.992),
    legend.justification.inside = c(1, 1),
    legend.direction = "vertical",
    legend.text = element_text(size = 17.5, colour = "#20252B", lineheight = 1.12),
    legend.key.spacing.y = grid::unit(0.24, "cm"),
    legend.background = element_rect(fill = scales::alpha("white", 0.92), colour = "#C8CED6", linewidth = 0.45),
    legend.margin = margin(11, 14, 11, 14),
    plot.margin = margin(8, 16, 12, 12)
  )
save_figure(fraction_plot, "Figure_1_DDR_fraction", 13.5, max(12, nrow(fraction_data) * 0.28 + 4.2))

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

draw_percentile_heatmap <- function(data, value_column, stem, measurement_label) {
  values <- data |>
    mutate(
      GroupKey = paste(PXD, SampleGroup, sep = "__"),
      Value = .data[[value_column]]
    ) |>
    inner_join(
      groups |> select(GroupKey, RowOrder, KlaLabelEn, ReferenceLabelEn, Category),
      by = "GroupKey", relationship = "many-to-one"
    ) |>
    inner_join(role_map |> select(Role, BaseAccession, DisplayName), by = c("RegulatorBaseAccession" = "BaseAccession"), relationship = "many-to-many") |>
    filter(!is.na(Value)) |>
    mutate(
      PlotLabel = if (identical(value_column, "WholeProteomeRelativePercentile")) ReferenceLabelEn else KlaLabelEn
    ) |>
    distinct(Role, DisplayName, PlotLabel, .keep_all = TRUE) |>
    mutate(
      CategoryLabel = factor(
        as.character(Category), levels = c("normal_tissue", "cancer_tissue", "normal_cells", "cancer_cells"),
        labels = c("Normal/non-tumor tissues", "Cancer tissues", "Normal/non-tumor cells", "Cancer cells")
      ),
      RoleLabel = factor(Role, levels = role_order),
      PlotLabel = factor(PlotLabel, levels = rev(unique(PlotLabel[order(RowOrder)]))),
      DisplayName = factor(DisplayName, levels = unique(role_map$DisplayName))
    )
  assert(nrow(values) > 0L, paste("No values were available for", stem))
  plot <- ggplot(values, aes(x = DisplayName, y = PlotLabel, fill = Value)) +
    geom_tile(colour = "white", linewidth = 0.22) +
    facet_grid(CategoryLabel ~ RoleLabel, scales = "free", space = "free") +
    scale_fill_gradientn(
      colours = c("#FFFFFF", "#FFF3E0", "#FDBB84", "#FC8D59", "#B2182B"),
      values = scales::rescale(c(0, 20, 50, 80, 100)),
      limits = c(0, 100), na.value = "#D9D9D9", name = measurement_label
    ) +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 8.5, base_family = publication_font) +
    theme(
      panel.grid = element_blank(),
      strip.text.x = element_text(face = "bold", size = 9),
      strip.text.y = element_text(face = "bold", size = 8),
      strip.background = element_rect(fill = "#F2F2F2", colour = NA),
      axis.text.x = element_text(angle = 55, hjust = 1, vjust = 1, size = 7),
      axis.text.y = element_text(size = 7.2),
      legend.position = "bottom",
      legend.key.width = grid::unit(35, "mm"),
      plot.margin = margin(8, 10, 8, 8)
    )
  save_figure(plot, stem, 15.5, 11.5)
}

kla_percentiles <- fread(input_path("regulator_kla_percentiles_30.csv")) |>
  as_tibble() |>
  transmute(PXD, SampleGroup, RegulatorBaseAccession, RelativeKlaPercentile)
draw_percentile_heatmap(
  kla_percentiles, "RelativeKlaPercentile", "Figure_3b_Kla_regulator_percentiles",
  "Within-sample Kla\npercentile"
)

reference_percentiles <- fread(input_path("regulator_reference_percentiles_30.csv")) |>
  as_tibble() |>
  transmute(PXD, SampleGroup, RegulatorBaseAccession, WholeProteomeRelativePercentile)
draw_percentile_heatmap(
  reference_percentiles, "WholeProteomeRelativePercentile", "Figure_3a_reference_regulator_percentiles",
  "Whole-proteome\npercentile"
)

draw_exact_venn <- function(filename, stem, title) {
  membership <- fread(input_path(filename)) |>
    as_tibble()
  membership_columns <- paste0("In_", category_order)
  assert(all(membership_columns %in% names(membership)), paste("Invalid Venn input:", filename))
  sets <- lapply(category_order, function(category) {
    membership$BaseAccession[is_true(membership[[paste0("In_", category)]])]
  })
  labels <- unname(category_labels[category_order])
  names(sets) <- labels
  fit <- eulerr::euler(sets)
  render <- function(device, path, width, height, resolution = NULL) {
    if (identical(device, "png")) {
      grDevices::png(path, width = width, height = height, res = resolution, type = "cairo")
    } else {
      grDevices::cairo_pdf(path, width = width, height = height, family = publication_font)
    }
    on.exit(grDevices::dev.off(), add = TRUE)
    graphics::par(family = publication_font, mar = c(1, 1, 7.5, 1))
    diagram <- plot(
      fit,
      labels = FALSE,
      legend = list(labels = labels, side = "bottom", nrow = 1, ncol = 4, byrow = TRUE, cex = 0.9),
      quantities = list(cex = 0.9),
      fills = list(fill = c("#E69F00", "#D55E00", "#009E73", "#CC79A7"), alpha = 0.52),
      edges = list(col = "#4B4B4B", lwd = 1.2),
      main = list(label = title, cex = 0.85),
      sub = "Deduplicated by UniProt BaseAccession; ellipse areas are fitted proportional to set sizes",
      sub.cex = 0.78,
      quantities.cex = 1.0
    )
    print(diagram)
  }
  render("png", file.path(figure_dir, paste0(stem, ".png")), 2400, 2300, 300)
  render("pdf", file.path(figure_dir, paste0(stem, ".pdf")), 8.5, 8.0)
}

draw_exact_venn("venn_reference_ddr.csv", "Figure_2a_whole_proteome_DDR_Venn", "Reference whole-proteome DDR proteins across four tissue and cell categories")
draw_exact_venn("venn_kla_ddr.csv", "Figure_2b_Kla_DDR_Venn", "Kla and DDR proteins across four tissue and cell categories")
draw_exact_venn("venn_reference.csv", "Supplementary_Figure_S1a_whole_proteome_Venn", "Reference whole-proteome proteins across four tissue and cell categories")
draw_exact_venn("venn_all_kla.csv", "Supplementary_Figure_S1b_Kla_proteome_Venn", "All Kla proteins across four tissue and cell categories")

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

zero_fill <- "#F1F3F5"
suppressing_fill <- "#2F3437"
guide_colour <- "#4B5563"

make_pathway_legend <- function() {
  legend_data <- data.table(
    X = 1:3,
    Label = c(
      "+1 promoting: solid pathway colour",
      "-1 suppressing: dark charcoal",
      "0 unassigned: light gray"
    ),
    Fill = c("#3C5488", suppressing_fill, zero_fill)
  )
  ggplot() +
    geom_rect(
      data = legend_data,
      aes(xmin = X - 0.31, xmax = X + 0.31, ymin = 0.38, ymax = 0.68, fill = Fill),
      colour = NA
    ) +
    geom_text(data = legend_data, aes(x = X, y = 0.18, label = Label), family = publication_font, size = 3.8, colour = "#374151") +
    scale_fill_identity() +
    coord_cartesian(xlim = c(0.45, 3.55), ylim = c(0.05, 0.78), clip = "off") +
    theme_void(base_family = publication_font) +
    theme(plot.margin = margin(0, 10, 0, 10))
}

legacy_matrix_panel <- function(key) {
  panel <- copy(pathway_panels[[key]]$matrix)
  n <- uniqueN(panel$BaseAccession)
  panel[, `:=`(
    PathwayOrder = match(as.character(Pathway), pathway_order),
    Y = 8 - match(as.character(Pathway), pathway_order),
    XMin = Rank - 0.5,
    XMax = Rank + 0.5
  )]
  panel[, `:=`(YMin = Y - 0.42, YMax = Y + 0.42)]
  ggplot() +
    geom_rect(data = panel[State == 0], aes(xmin = XMin, xmax = XMax, ymin = YMin, ymax = YMax), fill = zero_fill, colour = NA) +
    geom_rect(data = panel[State == 1], aes(xmin = XMin, xmax = XMax, ymin = YMin, ymax = YMax, fill = as.character(Pathway)), colour = NA) +
    geom_rect(data = panel[State == -1], aes(xmin = XMin, xmax = XMax, ymin = YMin, ymax = YMax), fill = suppressing_fill, colour = NA) +
    scale_fill_manual(values = pathway_colours, guide = "none") +
    scale_x_continuous(limits = c(0.5, n + 0.5), breaks = unique(as.integer(round(c(1, (n + 1) / 2, n)))), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0.5, 7.5), breaks = 7:1, labels = pathway_order, expand = c(0, 0)) +
    labs(
      title = paste0(pathway_panels[[key]]$label, "  |  n = ", n),
      x = "Protein rank (ascending signed pathway score)", y = NULL
    ) +
    theme_minimal(base_family = publication_font, base_size = 10) +
    theme(
      panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(colour = "#D1D5DB", linewidth = 0.34),
      axis.text.y = element_text(face = "bold", colour = "#374151", size = 9.8),
      axis.text.x = element_text(colour = "#4B5563", size = 8.5),
      axis.title.x = element_text(colour = "#374151", size = 9.5),
      plot.title = element_text(face = "bold", size = 11.5, colour = "#111827"),
      plot.margin = margin(7, 12, 7, 10)
    )
}

legacy_summary_panel <- function(key) {
  panel <- copy(pathway_panels[[key]]$matrix)
  n <- uniqueN(panel$BaseAccession)
  summary_data <- panel[, .(
    SuppressingCount = sum(State == -1),
    UnassignedCount = sum(State == 0),
    PromotingCount = sum(State == 1)
  ), by = Pathway]
  summary_data[, PathwayOrder := match(as.character(Pathway), pathway_order)]
  setorder(summary_data, PathwayOrder)
  summary_data[, `:=`(
    Y = 8 - PathwayOrder,
    PromotingFraction = PromotingCount / n,
    SuppressingFraction = SuppressingCount / n,
    PositiveLabel = sprintf("%d\n(%.1f%%)", PromotingCount, 100 * PromotingCount / n),
    NegativeLabel = sprintf("%d\n(%.1f%%)", SuppressingCount, 100 * SuppressingCount / n),
    ZeroLabel = paste0("0: ", UnassignedCount)
  )]
  ggplot() +
    geom_rect(data = summary_data, aes(xmin = 0, xmax = PromotingFraction, ymin = Y - 0.29, ymax = Y + 0.29, fill = as.character(Pathway)), colour = NA) +
    geom_rect(data = summary_data, aes(xmin = -SuppressingFraction, xmax = 0, ymin = Y - 0.29, ymax = Y + 0.29), fill = suppressing_fill, colour = NA) +
    geom_vline(xintercept = 0, colour = guide_colour, linewidth = 0.4) +
    geom_text(data = summary_data, aes(x = PromotingFraction + 0.008, y = Y, label = PositiveLabel), family = publication_font, hjust = 0, lineheight = 0.92, size = 2.75, colour = "#374151") +
    geom_text(data = summary_data, aes(x = pmin(-SuppressingFraction - 0.005, -0.005), y = Y, label = NegativeLabel), family = publication_font, hjust = 1, lineheight = 0.92, size = 2.75, colour = "#374151") +
    geom_text(data = summary_data, aes(x = 0.48, y = Y, label = ZeroLabel), family = publication_font, hjust = 0, size = 2.8, colour = "#6B7280") +
    scale_fill_manual(values = pathway_colours, guide = "none") +
    scale_x_continuous(limits = c(-0.13, 0.62), breaks = c(-0.1, 0, 0.1, 0.2, 0.3, 0.4), labels = function(x) paste0(abs(round(100 * x)), "%"), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0.5, 7.5), breaks = 7:1, labels = pathway_order, expand = c(0, 0)) +
    labs(title = paste0(pathway_panels[[key]]$label, "\nn = ", n), x = "Protein fraction", y = NULL) +
    theme_minimal(base_family = publication_font, base_size = 9) +
    theme(
      panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(colour = "#E5E7EB", linewidth = 0.32),
      axis.text.y = element_text(face = "bold", colour = "#374151", size = 8.8),
      axis.text.x = element_text(colour = "#4B5563", size = 7.8),
      axis.title.x = element_text(colour = "#374151", size = 8.5),
      plot.title = element_text(face = "bold", size = 10.5, colour = "#111827"),
      plot.margin = margin(8, 13, 8, 8)
    )
}

draw_pathway_matrices <- function(keys, stem) {
  plot <- wrap_plots(lapply(keys, legacy_matrix_panel), nrow = 1) /
    make_pathway_legend() + plot_layout(heights = c(1, 0.17))
  save_figure(plot, stem, 15.8, 6.7)
}

draw_pathway_summary <- function(key, stem) {
  plot <- legacy_summary_panel(key) /
    make_pathway_legend() + plot_layout(heights = c(1, 0.17))
  save_figure(plot, stem, 8.8, 7.2)
}

draw_pathway_matrices(c("cancer_tissue", "normal_tissue"), "Figure_2c_DDR_pathway_matrices_tissues")
draw_pathway_summary("cancer_tissue", "Figure_2d_DDR_pathway_summary_tumor_tissue")
draw_pathway_summary("normal_tissue", "Figure_2e_DDR_pathway_summary_non_tumor_tissue")
draw_pathway_matrices(c("cancer_cells", "normal_cells"), "Supplementary_Figure_S2a_DDR_pathway_matrices_cell_lines")
draw_pathway_summary("cancer_cells", "Supplementary_Figure_S2b_DDR_pathway_summary_cancer_cell_lines")
draw_pathway_summary("normal_cells", "Supplementary_Figure_S2c_DDR_pathway_summary_normal_cell_lines")

message("PASS: rebuilt only manuscript figures from the frozen 30-group publication input.")
