#!/usr/bin/env Rscript

# Rebuild only the figures reported in DDR_Kla_manuscript_V3.docx from the
# frozen, publication-scoped inputs.  All biological joins use BaseAccession;
# gene symbols are retained solely as display annotations in supplementary data.

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggVennDiagram)
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
fraction_category_order <- c("normal_tissue", "cancer_tissue", "normal_cells", "cancer_cells")
pathway_display_order <- c("BER", "NER", "MMR", "FA", "HR", "NHEJ", "AEJ")
venn_category_order <- c("normal_tissue", "cancer_tissue", "cancer_cells", "normal_cells")
fraction_display_names <- c(
  "pathological rotator cuff tendon" = "Pathological rotator cuff tendon",
  "normal human lung" = "Normal human lung",
  "hypertrophic scar" = "Hypertrophic scar",
  "adjacent skin" = "Adjacent skin",
  "human hippocampus" = "Human hippocampus",
  "normal pregnancy placenta" = "Normal pregnancy placenta",
  "human sperm" = "Human sperm",
  "BPH" = "Benign prostatic hyperplasia",
  "adjacent liver" = "Adjacent liver",
  "prostate cancer" = "Prostate cancer",
  "HCC" = "HCC",
  "HEK293T" = "HEK293T",
  "HMC3" = "HMC3",
  "pretreated HK-2" = "Pretreated HK-2",
  "MCF10A" = "MCF10A",
  "neural stem cells" = "Neural stem cells",
  "HUVEC control and Pg infection" = "HUVEC control/Pg infection",
  "HK-2 control and mannitol" = "HK-2 control/mannitol",
  "MCF7" = "MCF7",
  "HCT116" = "HCT116",
  "TALL-104" = "TALL-104",
  "HepG2 WT and SIRT1 or SIRT3 KO" = "HepG2 WT/SIRT1-KO/SIRT3-KO",
  "A549" = "A549",
  "MDA-MB-468" = "MDA-MB-468",
  "T-47D" = "T-47D",
  "PC-3M" = "PC-3M",
  "HCT116 control and Roseburia co-culture" = "HCT116 control/Roseburia co-culture",
  "glioblastoma stem cells" = "Glioblastoma stem cells",
  "RKO WT and GSK3B KO" = "RKO WT/GSK3B-KO"
)

pathway_display <- fread(input_path("pathway_display.csv"))
assert(
  identical(pathway_display$Pathway[order(pathway_display$PathwayOrder)], pathway_order),
  "The pathway display input must define BER, NER, MMR, FA, HR, AEJ and NHEJ in manuscript order."
)
pathway_colours <- stats::setNames(pathway_display$Color, pathway_display$Pathway)

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
assert(
  setequal(unique(groups$SampleGroup), names(fraction_display_names)),
  "Figure 1 requires the historical display label for every publication sample group."
)
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
    SampleGroup = first(SampleGroup),
    DisplayLabel = paste0(unname(fraction_display_names[first(SampleGroup)]), " · whole proteome Ref:", first(ReferencePXD)),
    Fraction = first(ReferenceDdrFraction),
    Ddr = first(ReferenceDdrProteinCount),
    Total = first(ReferenceProteinCount),
    BarType = "reference",
    .groups = "drop"
  )
kla_rows <- groups |>
  transmute(
    Category, RowOrder, PXD, SampleGroup,
    DisplayLabel = paste0(unname(fraction_display_names[SampleGroup]), " · Kla:", PXD),
    Fraction = KlaDdrFraction,
    Ddr = KlaDdrProteinCount,
    Total = KlaProteinCount,
    BarType = "kla"
  )
fraction_data <- bind_rows(reference_rows, kla_rows) |>
  mutate(
    BarPriority = if_else(BarType == "reference", 0L, 1L),
    FigureOrder = case_when(
      as.character(Category) == "cancer_tissue" & SampleGroup == "HCC" ~ 1L,
      as.character(Category) == "cancer_tissue" & SampleGroup == "prostate cancer" ~ 2L,
      TRUE ~ as.integer(RowOrder) + 2L
    )
  ) |>
  arrange(match(as.character(Category), fraction_category_order), FigureOrder, BarPriority, DisplayLabel) |>
  mutate(
    BarOrder = row_number(),
    PlotRow = factor(as.character(BarOrder), levels = rev(as.character(BarOrder))),
    CategoryLabel = factor(
      as.character(Category), levels = fraction_category_order,
      labels = unname(category_labels[fraction_category_order])
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
apply_figure_one_strip_fills <- function(plot) {
  plot_grob <- ggplotGrob(plot)
  strip_ids <- grep("^strip-l", plot_grob$layout$name)
  strip_ids <- strip_ids[order(plot_grob$layout$t[strip_ids])]
  strip_fills <- c("#DCEAF5", "#FCE7D4", "#DCEAF5", "#FCE7D4")
  assert(length(strip_ids) == length(strip_fills), "Figure 1 must contain four category strips.")
  for (index in seq_along(strip_ids)) {
    strip_grob <- plot_grob$grobs[[strip_ids[[index]]]]$grobs[[1]]
    background_id <- grep("^strip.background", strip_grob$childrenOrder)
    assert(length(background_id) == 1L, "Unable to identify a Figure 1 category-strip background.")
    strip_grob$children[[background_id]]$gp$fill <- strip_fills[[index]]
    strip_grob$children[[background_id]]$gp$col <- NA
    plot_grob$grobs[[strip_ids[[index]]]]$grobs[[1]] <- strip_grob
  }
  plot_grob
}
save_figure(
  apply_figure_one_strip_fills(fraction_plot),
  "Figure_1_DDR_fraction", 13.5, max(12, nrow(fraction_data) * 0.28 + 4.2)
)

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

draw_percentile_heatmap <- function(data, value_column, stem, measurement_label, colours) {
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
        as.character(Category), levels = category_order,
        labels = unname(category_labels[category_order])
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
      colours = colours,
      values = scales::rescale(c(0, 20, 50, 80, 100)),
      limits = c(0, 100),
      na.value = "#D9D9D9",
      name = measurement_label,
      guide = guide_colourbar(
        title.position = "left",
        title.hjust = 1,
        barwidth = grid::unit(76, "mm"),
        barheight = grid::unit(4.2, "mm")
      )
    ) +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 10.5, base_family = publication_font) +
    theme(
      panel.grid = element_blank(),
      strip.text.x = element_text(face = "bold", size = 12),
      strip.text.y.right = element_text(face = "bold", size = 11, angle = 180),
      strip.background = element_rect(fill = "#F2F2F2", colour = NA),
      axis.text.x = element_text(angle = 55, hjust = 1, vjust = 1, size = 10.2),
      axis.text.y = element_text(size = 9.6),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.title = element_text(size = 10.8),
      legend.text = element_text(size = 9.8),
      plot.margin = margin(10, 12, 10, 10)
    )
  save_figure(plot, stem, 15.5, 11.5)
}

kla_percentiles <- fread(input_path("regulator_kla_percentiles_30.csv")) |>
  as_tibble() |>
  transmute(PXD, SampleGroup, RegulatorBaseAccession, RelativeKlaPercentile)
draw_percentile_heatmap(
  kla_percentiles, "RelativeKlaPercentile", "Figure_3b_Kla_regulator_percentiles",
  "Lactylome (Kla) percentile",
  c("#FFFFFF", "#E8F3FA", "#9ECAE1", "#4292C6", "#08519C")
)

reference_percentiles <- fread(input_path("regulator_reference_percentiles_30.csv")) |>
  as_tibble() |>
  transmute(PXD, SampleGroup, RegulatorBaseAccession, WholeProteomeRelativePercentile)
draw_percentile_heatmap(
  reference_percentiles, "WholeProteomeRelativePercentile", "Figure_3a_reference_regulator_percentiles",
  "Whole-proteome percentile",
  c("#FFFFFF", "#FFF3E0", "#FDBB84", "#FC8D59", "#B2182B")
)

draw_exact_venn <- function(filename, stem, title) {
  membership <- fread(input_path(filename)) |>
    as_tibble()
  membership_columns <- paste0("In_", venn_category_order)
  assert(all(membership_columns %in% names(membership)), paste("Invalid Venn input:", filename))
  venn_labels <- c(
    normal_tissue = "non-tumor\ntissues",
    cancer_tissue = "tumor\ntissues",
    cancer_cells = "cancer cell\nlines",
    normal_cells = "normal\ncell lines"
  )
  venn_colours <- c(
    normal_tissue = "#0072B2",
    cancer_tissue = "#E69F00",
    cancer_cells = "#CC79A7",
    normal_cells = "#009E73"
  )
  lighten_region_colour <- function(region_id, strength = 0.42) {
    set_ids <- as.integer(strsplit(region_id, "/", fixed = TRUE)[[1]])
    rgb_values <- grDevices::col2rgb(unname(venn_colours)[set_ids])
    base_rgb <- rowMeans(rgb_values)
    mixed_rgb <- (1 - strength) * 255 + strength * base_rgb
    grDevices::rgb(mixed_rgb[1], mixed_rgb[2], mixed_rgb[3], maxColorValue = 255)
  }
  exact_sets <- lapply(venn_category_order, function(category) {
    membership$BaseAccession[is_true(membership[[paste0("In_", category)]])]
  })
  names(exact_sets) <- unname(venn_labels[venn_category_order])
  plot <- ggVennDiagram(
    exact_sets,
    label = "count",
    label_geom = "text",
    label_size = 5.5,
    label_color = "#1D1D1F",
    label_alpha = 1,
    label_font = publication_font,
    set_color = unname(venn_colours[venn_category_order]),
    set_size = 5.5,
    edge_size = 0.7,
    force_upset = FALSE
  )
  region_ids <- as.character(plot$layers[[1]]$data$id)
  plot$layers[[1]]$data$region_fill <- vapply(region_ids, lighten_region_colour, character(1))
  plot$layers[[1]]$mapping <- aes(
    x = .data$X,
    y = .data$Y,
    fill = .data$region_fill,
    group = .data$id
  )
  plot <- plot +
    scale_fill_identity(guide = "none") +
    labs(title = title) +
    theme_void(base_family = publication_font) +
    theme(
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.margin = margin(40, 120, 40, 120),
      text = element_text(family = publication_font),
      plot.title = element_text(hjust = 0.5, size = 20, colour = "#111111", margin = margin(b = 7))
    )
  plot$coordinates$clip <- "off"
  save_figure(plot, stem, 11, 8.5)
}

draw_exact_venn("venn_reference_ddr.csv", "Figure_2a_whole_proteome_DDR_Venn", "DDR Proteome")
draw_exact_venn("venn_kla_ddr.csv", "Figure_2b_Kla_DDR_Venn", "DDR Lactylome")
draw_exact_venn("venn_reference.csv", "Supplementary_Figure_S1a_whole_proteome_Venn", "Whole Proteome")
draw_exact_venn("venn_all_kla.csv", "Supplementary_Figure_S1b_Kla_proteome_Venn", "Whole Lactylome")

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

legacy_matrix_panel <- function(key, show_x_title = FALSE) {
  panel <- copy(pathway_panels[[key]]$matrix)
  n <- uniqueN(panel$BaseAccession)
  panel[, `:=`(
    PathwayOrder = match(as.character(Pathway), pathway_display_order),
    Y = 8 - match(as.character(Pathway), pathway_display_order),
    XMin = Rank - 1,
    XMax = Rank
  )]
  panel[, `:=`(YMin = Y - 0.42, YMax = Y + 0.42)]
  category_label <- pathway_panels[[key]]$label
  label_x <- n + max(6, ceiling(n * 0.025))
  ggplot() +
    geom_rect(data = panel[State == 0], aes(xmin = XMin, xmax = XMax, ymin = YMin, ymax = YMax), fill = zero_fill, colour = NA) +
    geom_rect(data = panel[State == 1], aes(xmin = XMin, xmax = XMax, ymin = YMin, ymax = YMax, fill = as.character(Pathway)), colour = NA) +
    geom_rect(data = panel[State == -1], aes(xmin = XMin, xmax = XMax, ymin = YMin, ymax = YMax), fill = suppressing_fill, colour = NA) +
    scale_fill_manual(values = pathway_colours, guide = "none") +
    geom_text(data = data.frame(X = label_x, Y = 4, Label = category_label), aes(x = X, y = Y, label = Label), inherit.aes = FALSE, angle = 90, hjust = 0.5, family = publication_font, size = 4.8, colour = "#20252B") +
    scale_x_continuous(limits = c(0, label_x + 2), breaks = unique(as.integer(round(c(0, n / 2, n)))), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0.5, 7.5), breaks = 7:1, labels = pathway_display_order, expand = c(0, 0)) +
    labs(
      x = if (show_x_title) "Lactylated DDR Protein (#)" else NULL, y = NULL
    ) +
    theme_minimal(base_family = publication_font, base_size = 10) +
    theme(
      panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(colour = "#D1D5DB", linewidth = 0.34),
      axis.text.y = element_text(face = "bold", colour = "#374151", size = 9.8),
      axis.text.x = element_text(colour = "#4B5563", size = 8.5),
      axis.title.x = element_text(colour = "#20252B", size = 12, margin = margin(t = 8)),
      plot.margin = margin(7, 24, 4, 10)
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
  summary_data[, PathwayOrder := match(as.character(Pathway), pathway_display_order)]
  setorder(summary_data, PathwayOrder)
  summary_data[, `:=`(
    Y = 8 - PathwayOrder,
    PromotingFraction = PromotingCount / n,
    SuppressingFraction = SuppressingCount / n,
    PositiveLabel = sprintf("%d\n(%.1f%%)", PromotingCount, 100 * PromotingCount / n),
    NegativeLabel = sprintf("%d\n(%.1f%%)", SuppressingCount, 100 * SuppressingCount / n)
  )]
  ggplot() +
    geom_rect(data = summary_data, aes(xmin = 0, xmax = PromotingFraction, ymin = Y - 0.29, ymax = Y + 0.29, fill = as.character(Pathway)), colour = NA) +
    geom_rect(data = summary_data, aes(xmin = -SuppressingFraction, xmax = 0, ymin = Y - 0.29, ymax = Y + 0.29), fill = suppressing_fill, colour = NA) +
    geom_vline(xintercept = 0, colour = guide_colour, linewidth = 0.4) +
    geom_text(data = summary_data, aes(x = PromotingFraction + 0.008, y = Y, label = PositiveLabel), family = publication_font, hjust = 0, lineheight = 0.92, size = 2.75, colour = "#374151") +
    geom_text(data = summary_data, aes(x = pmin(-SuppressingFraction - 0.005, -0.005), y = Y, label = NegativeLabel), family = publication_font, hjust = 1, lineheight = 0.92, size = 2.75, colour = "#374151") +
    scale_fill_manual(values = pathway_colours, guide = "none") +
    scale_x_continuous(limits = c(-0.13, 0.45), breaks = c(-0.1, 0, 0.1, 0.2, 0.3, 0.4), labels = function(x) paste0(abs(round(100 * x)), "%"), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0.5, 7.5), breaks = 7:1, labels = pathway_display_order, expand = c(0, 0)) +
    labs(title = pathway_panels[[key]]$label, x = "Relative Portion (%)", y = NULL) +
    theme_minimal(base_family = publication_font, base_size = 9) +
    theme(
      panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(colour = "#E5E7EB", linewidth = 0.32),
      axis.text.y = element_text(face = "bold", colour = "#374151", size = 8.8),
      axis.text.x = element_text(colour = "#4B5563", size = 7.8),
      axis.title.x = element_text(colour = "#20252B", size = 11, margin = margin(t = 8)),
      plot.title = element_text(face = "plain", hjust = 1, size = 14, colour = "#20252B"),
      plot.margin = margin(8, 18, 8, 8)
    )
}

draw_pathway_matrices <- function(keys, stem) {
  plot <- wrap_plots(
    lapply(seq_along(keys), function(index) {
      legacy_matrix_panel(keys[[index]], show_x_title = index == length(keys))
    }),
    ncol = 1
  )
  save_figure(plot, stem, 10.8, 7.2)
}

draw_pathway_summary <- function(key, stem) {
  save_figure(legacy_summary_panel(key), stem, 8.0, 5.4)
}

draw_pathway_matrices(c("cancer_tissue", "normal_tissue"), "Figure_2c_DDR_pathway_matrices_tissues")
draw_pathway_summary("cancer_tissue", "Figure_2d_DDR_pathway_summary_tumor_tissue")
draw_pathway_summary("normal_tissue", "Figure_2e_DDR_pathway_summary_non_tumor_tissue")
draw_pathway_matrices(c("cancer_cells", "normal_cells"), "Supplementary_Figure_S2a_DDR_pathway_matrices_cell_lines")
draw_pathway_summary("cancer_cells", "Supplementary_Figure_S2b_DDR_pathway_summary_cancer_cell_lines")
draw_pathway_summary("normal_cells", "Supplementary_Figure_S2c_DDR_pathway_summary_normal_cell_lines")

message("PASS: rebuilt only manuscript figures from the frozen 30-group publication input.")
