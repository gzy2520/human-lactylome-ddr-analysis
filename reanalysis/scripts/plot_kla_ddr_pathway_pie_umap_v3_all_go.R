#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(Matrix)
  library(patchwork)
  library(readxl)
  library(uwot)
})

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
script_path <- if (length(file_arg) > 0L) {
  normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE)
} else {
  normalizePath(
    "reanalysis/scripts/plot_kla_ddr_pathway_pie_umap_v3_all_go.R",
    mustWork = TRUE
  )
}
project_root <- normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = TRUE)

workbook_path <- file.path(
  project_root,
  "data/identifier/260810乳酸化DDR基因评分表.xlsx"
)
uniprot_path <- file.path(
  project_root,
  "reanalysis/config/uniprot_kla_ddr_507_all_go_2026-08-10.tsv"
)
uniprot_metadata_path <- file.path(
  project_root,
  "reanalysis/config/uniprot_kla_ddr_507_all_go_2026-08-10_metadata.tsv"
)
raw_ddr_go_path <- file.path(
  project_root,
  paste0(
    "reanalysis/results/tables/kla_ddr_raw_go_umap_33groups/",
    "protein_raw_go_term_binary_matrix.csv"
  )
)
v1_coordinates_path <- file.path(
  project_root,
  paste0(
    "reanalysis/results/tables/kla_ddr_raw_go_umap_33groups/",
    "umap_coordinates_fixed.csv"
  )
)
v2_coordinates_path <- file.path(
  project_root,
  paste0(
    "reanalysis/results/tables/kla_ddr_pathway_pie_umap_33groups_v2_spread/",
    "umap_coordinates_v2_spread.csv"
  )
)
table_dir <- file.path(
  project_root,
  "reanalysis/results/tables/kla_ddr_pathway_pie_umap_33groups_v3_all_go"
)
figure_dir <- file.path(
  project_root,
  "reanalysis/results/figures/kla_ddr_pathway_pie_umap_33groups_v3_all_go"
)
report_path <- file.path(
  project_root,
  "reanalysis/reports/UMAP_PATHWAY_PIE_33GROUP_V3_ALL_GO.md"
)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

required_inputs <- c(
  workbook_path,
  uniprot_path,
  uniprot_metadata_path,
  raw_ddr_go_path,
  v1_coordinates_path,
  v2_coordinates_path
)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs) > 0L) {
  stop("Missing required input(s): ", paste(missing_inputs, collapse = "; "))
}

stop_if_false <- function(condition, message) {
  if (!isTRUE(condition)) {
    stop(message, call. = FALSE)
  }
}

base_accession <- function(x) {
  x <- trimws(as.character(x))
  x <- sub("^(sp|tr)\\|", "", x)
  x <- sub("\\|.*$", "", x)
  x <- sub("^.*:", "", x)
  sub("-[0-9]+$", "", x)
}

relative_path <- function(path) {
  sub(
    paste0("^", gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", project_root), "/?"),
    "",
    path
  )
}

pathway_info <- data.table(
  Pathway = c(
    "HR",
    "NHEJ",
    "AEJ",
    "BER",
    "NER",
    "MMR",
    "FA",
    "Chromatin Interaction",
    "Others (Transcription, RNA processing and proteostasis)"
  ),
  DisplayLabel = c(
    "HR",
    "NHEJ",
    "AEJ",
    "BER",
    "NER",
    "MMR",
    "FA",
    "Chromatin interaction",
    "Other support"
  ),
  ColorName = c(
    "Chambray",
    "Cinnabar",
    "PersianGreen",
    "Shakespeare",
    "Apricot",
    "WildBlueYonder",
    "MonteCarlo",
    "RomanCoffee",
    "NeutralCharcoal"
  ),
  Color = c(
    "#3C5488",
    "#E64B35",
    "#00A087",
    "#4DBBD5",
    "#F39B7F",
    "#8491B4",
    "#91D1C2",
    "#7E6148",
    "#6F6F6F"
  ),
  PathwayOrder = seq_len(9L)
)
pathway_columns <- pathway_info$Pathway

scoring_raw <- as.data.table(
  read_excel(workbook_path, sheet = "评分表", .name_repair = "minimal")
)
scoring <- scoring_raw[
  !is.na(BaseAccession) & nzchar(trimws(as.character(BaseAccession)))
]
scoring[, BaseAccession := base_accession(BaseAccession)]
scoring[, ID := as.integer(ID)]
for (column in pathway_columns) {
  set(scoring, j = column, value = as.integer(scoring[[column]]))
}

raw_ddr_go <- fread(raw_ddr_go_path)
v1_coordinates <- fread(v1_coordinates_path)
v2_coordinates <- fread(v2_coordinates_path)
uniprot <- fread(
  uniprot_path,
  sep = "\t",
  quote = "\"",
  na.strings = character(),
  check.names = FALSE
)
setnames(uniprot, "BaseAccession", "BaseAccession", skip_absent = TRUE)
uniprot[, BaseAccession := base_accession(BaseAccession)]

stop_if_false(
  nrow(scoring) == 507L && uniqueN(scoring$BaseAccession) == 507L,
  "The scoring workbook does not contain 507 unique BaseAccessions."
)
stop_if_false(
  nrow(uniprot) == 507L && uniqueN(uniprot$BaseAccession) == 507L,
  "The UniProt annotation cache does not contain 507 unique BaseAccessions."
)
stop_if_false(
  setequal(scoring$BaseAccession, uniprot$BaseAccession),
  "The scoring and UniProt protein sets do not match exactly."
)
stop_if_false(
  setequal(scoring$BaseAccession, raw_ddr_go$BaseAccession),
  "The scoring set is not the same 507-protein Kla-intersection-DDR set."
)
stop_if_false(
  setequal(v1_coordinates$BaseAccession, scoring$BaseAccession) &&
    setequal(v2_coordinates$BaseAccession, scoring$BaseAccession),
  "V1/V2 coordinates do not match the 507-protein set."
)

aspect_spec <- data.table(
  Ontology = c("BP", "CC", "MF"),
  SourceColumn = c(
    "Gene Ontology (biological process)",
    "Gene Ontology (cellular component)",
    "Gene Ontology (molecular function)"
  )
)
stop_if_false(
  all(aspect_spec$SourceColumn %in% names(uniprot)),
  "Expected UniProt GO aspect columns are absent."
)

parse_go_cell <- function(accession, ontology, value) {
  if (is.na(value) || !nzchar(value)) {
    return(NULL)
  }
  pieces <- strsplit(value, "; ", fixed = TRUE)[[1L]]
  matched <- grepl("\\[GO:[0-9]{7}\\]$", pieces)
  pieces <- pieces[matched]
  if (length(pieces) == 0L) {
    return(NULL)
  }
  data.table(
    BaseAccession = accession,
    Ontology = ontology,
    GO_ID = sub("^.*\\[(GO:[0-9]{7})\\]$", "\\1", pieces),
    GO_Term = sub("\\s*\\[GO:[0-9]{7}\\]$", "", pieces)
  )
}

go_long <- rbindlist(
  lapply(seq_len(nrow(aspect_spec)), function(aspect_index) {
    ontology <- aspect_spec$Ontology[[aspect_index]]
    source_column <- aspect_spec$SourceColumn[[aspect_index]]
    rbindlist(
      lapply(seq_len(nrow(uniprot)), function(protein_index) {
        parse_go_cell(
          uniprot$BaseAccession[[protein_index]],
          ontology,
          uniprot[[source_column]][[protein_index]]
        )
      }),
      use.names = TRUE
    )
  }),
  use.names = TRUE
)
go_long <- unique(go_long, by = c("BaseAccession", "Ontology", "GO_ID"))
setorder(go_long, BaseAccession, Ontology, GO_ID)

stop_if_false(nrow(go_long) == 13738L, "Expected 13,738 direct protein-GO annotations.")
stop_if_false(uniqueN(go_long$GO_ID) == 3461L, "Expected 3,461 unique direct GO terms.")
stop_if_false(
  uniqueN(go_long[Ontology == "BP"]$BaseAccession) == 507L,
  "At least one protein lacks a direct UniProt biological-process annotation."
)

protein_levels <- sort(scoring$BaseAccession)
term_levels <- sort(unique(go_long$GO_ID))
go_binary <- sparseMatrix(
  i = match(go_long$BaseAccession, protein_levels),
  j = match(go_long$GO_ID, term_levels),
  x = 1,
  dims = c(length(protein_levels), length(term_levels)),
  dimnames = list(protein_levels, term_levels)
)
go_binary@x[] <- 1

term_document_frequency <- Matrix::colSums(go_binary > 0)
term_idf <- log((nrow(go_binary) + 1) / (term_document_frequency + 1)) + 1
go_tfidf <- go_binary %*% Diagonal(x = as.numeric(term_idf))

stop_if_false(
  identical(dim(go_binary), c(507L, 3461L)),
  "All-GO matrix is not 507 proteins by 3,461 terms."
)
stop_if_false(sum(go_binary) == 13738L, "All-GO matrix does not contain 13,738 hits.")
stop_if_false(all(Matrix::rowSums(go_binary) >= 1L), "At least one protein lacks GO input.")

signature_table <- go_long[
  ,
  .(
    DirectGOCount = .N,
    BPTermCount = sum(Ontology == "BP"),
    CCTermCount = sum(Ontology == "CC"),
    MFTermCount = sum(Ontology == "MF"),
    DirectGOSignature = paste(sort(GO_ID), collapse = ";")
  ),
  by = BaseAccession
]
signature_counts <- signature_table[
  ,
  .(
    ProteinsWithSignature = .N,
    Accessions = paste(sort(BaseAccession), collapse = ";")
  ),
  by = DirectGOSignature
][order(-ProteinsWithSignature, Accessions)]
duplicate_signatures <- signature_counts[ProteinsWithSignature > 1L]

stop_if_false(
  uniqueN(signature_table$DirectGOSignature) == 506L,
  "Expected 506 unique full-GO direct annotation profiles."
)
stop_if_false(
  nrow(duplicate_signatures) == 1L &&
    duplicate_signatures$Accessions[[1L]] == "P35250;P40937",
  "The expected RFC2/RFC5 duplicate all-GO profile was not recovered."
)

umap_parameters <- list(
  n_neighbors = 15L,
  min_dist = 3.0,
  spread = 10.0,
  repulsion_strength = 1.8,
  negative_sample_rate = 10L,
  n_epochs = 750L,
  seed = 25L
)

set.seed(umap_parameters$seed)
v3_matrix <- uwot::umap(
  X = as.matrix(go_tfidf),
  n_neighbors = umap_parameters$n_neighbors,
  n_components = 2L,
  metric = "cosine",
  n_epochs = umap_parameters$n_epochs,
  scale = FALSE,
  init = "random",
  spread = umap_parameters$spread,
  min_dist = umap_parameters$min_dist,
  repulsion_strength = umap_parameters$repulsion_strength,
  negative_sample_rate = umap_parameters$negative_sample_rate,
  fast_sgd = FALSE,
  n_threads = 1L,
  n_sgd_threads = 1L,
  seed = umap_parameters$seed,
  verbose = FALSE
)
stop_if_false(
  identical(dim(v3_matrix), c(507L, 2L)) && all(is.finite(v3_matrix)),
  "V3 UMAP did not produce 507 finite two-dimensional coordinates."
)

v3_coordinates <- data.table(
  BaseAccession = protein_levels,
  UMAP_1 = v3_matrix[, 1L],
  UMAP_2 = v3_matrix[, 2L]
)
scoring <- scoring[match(protein_levels, BaseAccession)]
plot_data <- merge(
  v3_coordinates,
  scoring,
  by = "BaseAccession",
  all.x = TRUE,
  sort = FALSE
)
plot_data <- merge(
  plot_data,
  signature_table[
    ,
    .(BaseAccession, DirectGOCount, BPTermCount, CCTermCount, MFTermCount)
  ],
  by = "BaseAccession",
  all.x = TRUE,
  sort = FALSE
)
setorder(plot_data, BaseAccession)

score_matrix <- as.matrix(scoring[, ..pathway_columns])
stop_if_false(!anyNA(score_matrix), "The pathway score matrix contains missing values.")
stop_if_false(
  all(score_matrix %in% c(-1L, 0L, 1L)),
  "Pathway scores outside {-1, 0, 1} were detected."
)

assignment_long <- melt(
  scoring,
  id.vars = c("ID", "BaseAccession", "GeneSymbol", "ProteinName", "Note"),
  measure.vars = pathway_columns,
  variable.name = "Pathway",
  value.name = "Score",
  variable.factor = FALSE
)[Score != 0L]
assignment_long <- merge(
  assignment_long,
  pathway_info[
    ,
    .(Pathway, DisplayLabel, ColorName, Color, PathwayOrder)
  ],
  by = "Pathway",
  all.x = TRUE,
  sort = FALSE
)
setorder(assignment_long, BaseAccession, PathwayOrder)
assignment_long[
  ,
  Direction := fifelse(Score == 1L, "Promoting (+1)", "Suppressing (-1)")
]

plot_data[
  ,
  PositiveAssignmentCount := rowSums(.SD == 1L),
  .SDcols = pathway_columns
]
plot_data[
  ,
  NegativeAssignmentCount := rowSums(.SD == -1L),
  .SDcols = pathway_columns
]
plot_data[
  ,
  TotalAssignmentCount := PositiveAssignmentCount + NegativeAssignmentCount
]

stop_if_false(nrow(assignment_long) == 1175L, "Expected 1,175 nonzero assignments.")
stop_if_false(sum(assignment_long$Score == 1L) == 1108L, "Expected 1,108 promoting assignments.")
stop_if_false(sum(assignment_long$Score == -1L) == 67L, "Expected 67 suppressing assignments.")
stop_if_false(plot_data[TotalAssignmentCount == 0L, .N] == 22L, "Expected 22 all-zero proteins.")

pie_radius <- 0.95
arc_points <- 20L

make_protein_wedges <- function(accession, center_x, center_y) {
  assignments <- assignment_long[BaseAccession == accession]
  sector_count <- nrow(assignments)
  if (sector_count == 0L) {
    return(NULL)
  }
  pieces <- vector("list", sector_count)
  for (sector_index in seq_len(sector_count)) {
    start_angle <- pi / 2 - 2 * pi * (sector_index - 1L) / sector_count
    end_angle <- pi / 2 - 2 * pi * sector_index / sector_count
    theta <- seq(start_angle, end_angle, length.out = arc_points)
    assignment <- assignments[sector_index]
    pieces[[sector_index]] <- data.table(
      BaseAccession = accession,
      WedgeID = paste(accession, sector_index, sep = "__"),
      Pathway = assignment$Pathway,
      DisplayLabel = assignment$DisplayLabel,
      Score = assignment$Score,
      Direction = assignment$Direction,
      X = c(center_x, center_x + pie_radius * cos(theta), center_x),
      Y = c(center_y, center_y + pie_radius * sin(theta), center_y),
      FillKey = if (assignment$Score == 1L) assignment$Pathway else "NegativeHollow",
      BorderKey = if (assignment$Score == 1L) "PositiveBoundary" else assignment$Pathway
    )
  }
  rbindlist(pieces)
}

wedge_polygons <- rbindlist(
  lapply(seq_len(nrow(plot_data)), function(index) {
    make_protein_wedges(
      plot_data$BaseAccession[[index]],
      plot_data$UMAP_1[[index]],
      plot_data$UMAP_2[[index]]
    )
  })
)
stop_if_false(
  uniqueN(wedge_polygons$WedgeID) == 1175L,
  "Wedge construction did not yield one sector per nonzero assignment."
)

fill_values <- c(
  setNames(pathway_info$Color, pathway_info$Pathway),
  NegativeHollow = "#FFFFFF"
)
border_values <- c(
  PositiveBoundary = "#FFFFFF",
  setNames(pathway_info$Color, pathway_info$Pathway)
)

summary_counts <- assignment_long[
  ,
  .(
    Positive = sum(Score == 1L),
    Negative = sum(Score == -1L),
    TotalNonzero = .N
  ),
  by = .(Pathway, DisplayLabel, ColorName, Color, PathwayOrder)
]
setorder(summary_counts, PathwayOrder)

nearest_neighbor_metrics <- function(coordinates, version, radius = pie_radius) {
  coordinate_matrix <- as.matrix(coordinates[, .(UMAP_1, UMAP_2)])
  distances <- as.matrix(stats::dist(coordinate_matrix))
  diag(distances) <- Inf
  nearest <- apply(distances, 1L, min)
  upper_distances <- distances[upper.tri(distances)]
  data.table(
    Version = version,
    PieRadius = radius,
    PieDiameter = 2 * radius,
    MinimumNearestNeighborDistance = min(nearest),
    P10NearestNeighborDistance = unname(quantile(nearest, 0.10)),
    MedianNearestNeighborDistance = median(nearest),
    MeanNearestNeighborDistance = mean(nearest),
    ProteinsWithNearestNeighborCloserThanPieDiameter = sum(nearest < 2 * radius),
    PairDistancesCloserThanPieDiameter = sum(upper_distances < 2 * radius),
    UMAP1Range = diff(range(coordinate_matrix[, 1L])),
    UMAP2Range = diff(range(coordinate_matrix[, 2L]))
  )
}

overlap_comparison <- rbind(
  nearest_neighbor_metrics(v1_coordinates, "V1_DDR_GO_fixed_compact"),
  nearest_neighbor_metrics(v2_coordinates, "V2_DDR_GO_refitted_spread"),
  nearest_neighbor_metrics(v3_coordinates, "V3_all_GO_TFIDF_spread")
)

unassigned_data <- plot_data[TotalAssignmentCount == 0L]

pie_plot <- ggplot() +
  geom_point(
    data = unassigned_data,
    aes(x = UMAP_1, y = UMAP_2, shape = "No scored pathway/function"),
    color = "#A6A6A6",
    size = 2.35,
    stroke = 0.45
  ) +
  geom_polygon(
    data = wedge_polygons,
    aes(x = X, y = Y, group = WedgeID, fill = FillKey, color = BorderKey),
    linewidth = 0.32,
    linejoin = "round"
  ) +
  scale_fill_manual(
    values = fill_values,
    breaks = pathway_info$Pathway,
    labels = setNames(pathway_info$DisplayLabel, pathway_info$Pathway),
    name = "Pathway / function",
    drop = FALSE
  ) +
  scale_color_manual(values = border_values, guide = "none") +
  scale_shape_manual(
    values = c("No scored pathway/function" = 13),
    name = NULL
  ) +
  coord_equal() +
  labs(
    title = "Kla-DDR pathway assignments in a full-GO functional UMAP",
    subtitle = paste0(
      "507 proteins | 3,461 UniProt BP/CC/MF terms | TF-IDF-weighted cosine distance"
    ),
    x = "UMAP 1",
    y = "UMAP 2",
    caption = paste0(
      "Filled sectors: promoting score +1. White sectors with pathway-colored borders: ",
      "suppressing score -1. Gray circled crosses: no scored pathway/function.\n",
      "Sector area denotes membership, not abundance. Coordinates reflect full direct ",
      "GO annotation similarity; pathway scores were not used to fit the UMAP."
    )
  ) +
  guides(
    fill = guide_legend(
      override.aes = list(color = NA),
      order = 1,
      ncol = 1
    ),
    shape = guide_legend(order = 2)
  ) +
  theme_classic(base_size = 10, base_family = "sans") +
  theme(
    axis.line = element_line(linewidth = 0.4, color = "#333333"),
    axis.ticks = element_line(linewidth = 0.35, color = "#333333"),
    axis.text = element_text(color = "#333333", size = 8),
    axis.title = element_text(size = 9),
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 9, color = "#4D4D4D"),
    plot.caption = element_text(size = 7.4, color = "#5A5A5A", hjust = 0),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 9),
    legend.text = element_text(size = 8),
    legend.key.height = grid::unit(0.34, "cm"),
    plot.margin = margin(8, 8, 7, 8)
  )

bar_long <- melt(
  summary_counts,
  id.vars = c("Pathway", "DisplayLabel", "ColorName", "Color", "PathwayOrder"),
  measure.vars = c("Positive", "Negative"),
  variable.name = "Direction",
  value.name = "Count"
)
bar_long[, PlotCount := fifelse(Direction == "Negative", -Count, Count)]
bar_long[
  ,
  DirectionLabel := fifelse(
    Direction == "Positive",
    "Promoting (+1)",
    "Suppressing (-1)"
  )
]
bar_long[
  ,
  DisplayLabel := factor(DisplayLabel, levels = rev(pathway_info$DisplayLabel))
]

bar_plot <- ggplot(
  bar_long,
  aes(x = PlotCount, y = DisplayLabel, fill = DirectionLabel)
) +
  geom_vline(xintercept = 0, linewidth = 0.4, color = "#555555") +
  geom_col(width = 0.68) +
  geom_text(
    data = bar_long[Direction == "Positive"],
    aes(label = Count),
    hjust = -0.25,
    size = 2.8
  ) +
  geom_text(
    data = bar_long[Direction == "Negative"],
    aes(label = Count),
    hjust = 1.25,
    size = 2.8
  ) +
  scale_fill_manual(
    values = c(
      "Promoting (+1)" = "#3C5488",
      "Suppressing (-1)" = "#E64B35"
    ),
    name = "Evidence direction"
  ) +
  scale_x_continuous(
    limits = c(-35, 345),
    breaks = c(-20, 0, 100, 200, 300),
    labels = function(x) abs(x),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    title = "Signed pathway/function assignment counts",
    x = "Number of protein-pathway assignments",
    y = NULL
  ) +
  theme_classic(base_size = 9.5, base_family = "sans") +
  theme(
    plot.title = element_text(face = "bold", size = 10.5),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(color = "#333333", size = 8),
    axis.title.x = element_text(size = 9),
    legend.position = "top",
    legend.justification = "left",
    legend.title = element_text(face = "bold", size = 8.5),
    legend.text = element_text(size = 8),
    plot.margin = margin(4, 14, 6, 8)
  )

combined_plot <- pie_plot / bar_plot +
  plot_layout(heights = c(3.45, 1.65))

save_plot_formats <- function(plot, stem, width, height) {
  png_path <- file.path(figure_dir, paste0(stem, ".png"))
  pdf_path <- file.path(figure_dir, paste0(stem, ".pdf"))
  svg_path <- file.path(figure_dir, paste0(stem, ".svg"))
  ggsave(
    png_path,
    plot,
    width = width,
    height = height,
    units = "in",
    dpi = 600,
    bg = "white",
    limitsize = TRUE
  )
  ggsave(
    pdf_path,
    plot,
    width = width,
    height = height,
    units = "in",
    bg = "white",
    limitsize = TRUE
  )
  grDevices::svg(
    svg_path,
    width = width,
    height = height,
    onefile = FALSE,
    bg = "white"
  )
  print(plot)
  grDevices::dev.off()
  invisible(c(png_path, pdf_path, svg_path))
}

pie_paths <- save_plot_formats(
  pie_plot,
  "kla_ddr_pathway_hollow_negative_pie_umap_v3_all_go",
  width = 10.6,
  height = 7.8
)
combined_paths <- save_plot_formats(
  combined_plot,
  "kla_ddr_pathway_hollow_negative_umap_and_summary_v3_all_go",
  width = 10.6,
  height = 11.0
)

term_dictionary <- unique(go_long[, .(Ontology, GO_ID, GO_Term)])
term_frequency <- go_long[
  ,
  .(ProteinCount = uniqueN(BaseAccession)),
  by = .(Ontology, GO_ID)
]
term_dictionary <- merge(
  term_dictionary,
  term_frequency,
  by = c("Ontology", "GO_ID"),
  all.x = TRUE
)
term_dictionary[
  ,
  IDFWeight := log((length(protein_levels) + 1) / (ProteinCount + 1)) + 1
]
setorder(term_dictionary, Ontology, GO_ID)

aspect_summary <- go_long[
  ,
  .(
    ProteinGOHits = .N,
    UniqueGOTerms = uniqueN(GO_ID),
    AnnotatedProteins = uniqueN(BaseAccession)
  ),
  by = Ontology
]
aspect_summary[, AspectOrder := match(Ontology, c("BP", "CC", "MF"))]
setorder(aspect_summary, AspectOrder)
aspect_summary[, AspectOrder := NULL]

annotation_profile_summary <- data.table(
  Metric = c(
    "Proteins",
    "ProteinGOHits",
    "UniqueDirectGOTerms",
    "UniqueDirectGOProfiles",
    "ProteinsInDuplicateProfiles",
    "MaximumDuplicateProfileSize",
    "MinimumTermsPerProtein",
    "MedianTermsPerProtein",
    "MeanTermsPerProtein",
    "MaximumTermsPerProtein"
  ),
  Value = c(
    nrow(go_binary),
    sum(go_binary),
    ncol(go_binary),
    uniqueN(signature_table$DirectGOSignature),
    signature_counts[ProteinsWithSignature > 1L, sum(ProteinsWithSignature)],
    max(signature_counts$ProteinsWithSignature),
    min(signature_table$DirectGOCount),
    median(signature_table$DirectGOCount),
    mean(signature_table$DirectGOCount),
    max(signature_table$DirectGOCount)
  )
)

parameter_table <- data.table(
  Parameter = c(
    "AnalysisUnit",
    "ProteinAnalysisKey",
    "UMAPInput",
    "GOAspects",
    "GOAnnotationType",
    "GOFeatureWeighting",
    "PathwayScoresUsedInUMAP",
    "SampleDetectionUsedInUMAP",
    "DistanceMetric",
    "NNeighbors",
    "MinDist",
    "Spread",
    "RepulsionStrength",
    "NegativeSampleRate",
    "NEpochs",
    "Initialization",
    "RandomSeed",
    "NearestNeighborThreads",
    "SGDThreads",
    "PieRadius",
    "PositiveSectorStyle",
    "NegativeSectorStyle",
    "AllZeroProteinStyle",
    "Palette"
  ),
  Value = as.character(c(
    "one unique protein per pie",
    "isoform-stripped UniProt BaseAccession",
    "507 x 3,461 protein x direct GO term matrix",
    "BP;CC;MF",
    "direct UniProt annotations; no ancestor propagation",
    "IDF = log((507 + 1) / (protein frequency + 1)) + 1",
    "FALSE",
    "FALSE",
    "cosine on TF-IDF-weighted binary features",
    umap_parameters$n_neighbors,
    umap_parameters$min_dist,
    umap_parameters$spread,
    umap_parameters$repulsion_strength,
    umap_parameters$negative_sample_rate,
    umap_parameters$n_epochs,
    "random",
    umap_parameters$seed,
    1L,
    1L,
    pie_radius,
    "solid pathway color with white sector boundary",
    "white hollow sector with pathway-colored border",
    "light-gray circled cross (pch 13)",
    "NPG-derived colors plus neutral charcoal"
  ))
)

input_audit <- data.table(
  InputRole = c(
    "pathway scoring workbook",
    "UniProt full direct GO annotation cache",
    "UniProt retrieval metadata",
    "33-group-derived 507-protein scope verification",
    "V1 coordinate comparison only",
    "V2 coordinate comparison only"
  ),
  Path = relative_path(required_inputs),
  MD5 = unname(tools::md5sum(required_inputs))
)

binary_dense <- as.data.table(as.matrix(go_binary))
binary_dense[, BaseAccession := protein_levels]
setcolorder(binary_dense, c("BaseAccession", term_levels))

fwrite(
  go_long,
  file.path(table_dir, "uniprot_direct_go_annotation_long.csv")
)
fwrite(
  term_dictionary,
  file.path(table_dir, "all_go_term_dictionary_and_idf.csv")
)
fwrite(
  binary_dense,
  file.path(table_dir, "protein_all_go_direct_binary_matrix.csv")
)
fwrite(
  signature_table,
  file.path(table_dir, "protein_all_go_profile_summary.csv")
)
fwrite(
  duplicate_signatures,
  file.path(table_dir, "duplicate_all_go_profiles.csv")
)
fwrite(
  aspect_summary,
  file.path(table_dir, "go_aspect_annotation_summary.csv")
)
fwrite(
  annotation_profile_summary,
  file.path(table_dir, "all_go_annotation_profile_summary.csv")
)
fwrite(v3_coordinates, file.path(table_dir, "umap_coordinates_v3_all_go.csv"))
fwrite(
  plot_data,
  file.path(table_dir, "pathway_umap_plot_data_v3_all_go.csv")
)
fwrite(
  assignment_long,
  file.path(table_dir, "pathway_assignment_long_v3.csv")
)
fwrite(
  wedge_polygons,
  file.path(table_dir, "hollow_negative_wedge_polygons_v3.csv")
)
fwrite(summary_counts, file.path(table_dir, "pathway_assignment_summary.csv"))
fwrite(
  overlap_comparison,
  file.path(table_dir, "pie_overlap_comparison_v1_v2_v3.csv")
)
fwrite(parameter_table, file.path(table_dir, "umap_v3_parameters.csv"))
fwrite(pathway_info, file.path(table_dir, "pathway_color_key_v3.csv"))
fwrite(input_audit, file.path(table_dir, "input_file_audit.csv"))
writeLines(
  capture.output(sessionInfo()),
  file.path(table_dir, "session_info.txt"),
  useBytes = TRUE
)

v3_overlap_pairs <- overlap_comparison[
  Version == "V3_all_GO_TFIDF_spread",
  PairDistancesCloserThanPieDiameter
]
v3_close_proteins <- overlap_comparison[
  Version == "V3_all_GO_TFIDF_spread",
  ProteinsWithNearestNeighborCloserThanPieDiameter
]
uniprot_metadata <- fread(uniprot_metadata_path, sep = "\t")
uniprot_release <- uniprot_metadata[Key == "UniProtRelease", Value]
uniprot_release_date <- uniprot_metadata[Key == "UniProtReleaseDate", Value]

report_lines <- c(
  "# 33组Kla∩DDR通路饼图UMAP V3：全GO功能空间",
  "",
  "## 为什么改为全GO",
  "",
  "- 正式分析对象仍是33组数据得到的507个Kla∩DDR蛋白；一个饼代表一个去isoform的UniProt `BaseAccession`。",
  "- V1/V2仅使用66个DDR相关GO term，205种输入模式中最大重复组包含102个蛋白，难以支撑可读的功能布局。",
  paste0(
    "- V3从UniProt ", uniprot_release, "（", uniprot_release_date,
    "）取得这507个蛋白的全部直接BP/CC/MF注释：13,738条蛋白–GO关系、3,461个不同GO term。"
  ),
  "- 全GO直接注释形成506种模式；仅RFC2（P35250）与RFC5（P40937）完全相同。二者属于同一Replication factor C复合体近缘小亚基，保留接近具有生物学合理性。",
  "",
  "## UMAP输入和边界",
  "",
  "- 使用507 × 3,461蛋白×GO二值矩阵；不对GO term预先归入HR/NHEJ等通路。",
  "- 纳入BP、CC和MF直接注释，不扩展GO祖先节点，避免人为重复加入宽泛上位term。",
  "- 对每个GO term使用逆文档频率权重：`IDF = log((507 + 1) / (注释蛋白数 + 1)) + 1`，再以cosine距离拟合UMAP；常见的nucleus、cytoplasm等term因此影响较低。",
  "- 9列通路评分和33组检出模式均不参与UMAP，只在坐标生成后用于饼图编码或范围核验。",
  "- UMAP邻近表示全GO注释相似性，不代表物理互作、共同表达、乳酸化丰度或通路激活强度。",
  "",
  "## 扇形、全0蛋白与配色",
  "",
  "- `+1`：对应通路颜色的实心扇形；`-1`：白色空心扇形并以对应通路颜色描边。",
  "- 9列全为0的22个蛋白改为浅灰色空心圆内带叉（circled cross），明确表示未获当前评分表的通路/功能分配，不再画成类似第10类通路的实心灰点。",
  "- Chromatin interaction保留棕色`#7E6148`；Other support改为中性炭灰`#6F6F6F`，消除原两个棕色过近的问题。",
  "- 配色以NPG/Nature Reviews Cancer风格为基础，并使用中性炭灰扩展；导出600 dpi PNG和PDF/SVG矢量图。",
  "",
  "## 可读性",
  "",
  paste0(
    "- 饼半径为", pie_radius, "坐标单位；V3有", v3_overlap_pairs,
    "对蛋白中心距离小于一个饼直径，涉及", v3_close_proteins, "个蛋白。"
  ),
  "",
  "## 主要输出",
  "",
  "- 推荐图：`reanalysis/results/figures/kla_ddr_pathway_pie_umap_33groups_v3_all_go/kla_ddr_pathway_hollow_negative_pie_umap_v3_all_go.{png,pdf,svg}`",
  "- 组合图：`reanalysis/results/figures/kla_ddr_pathway_pie_umap_33groups_v3_all_go/kla_ddr_pathway_hollow_negative_umap_and_summary_v3_all_go.{png,pdf,svg}`",
  "- 坐标：`reanalysis/results/tables/kla_ddr_pathway_pie_umap_33groups_v3_all_go/umap_coordinates_v3_all_go.csv`",
  "- 全GO长表：`reanalysis/results/tables/kla_ddr_pathway_pie_umap_33groups_v3_all_go/uniprot_direct_go_annotation_long.csv`",
  "- 全GO二值矩阵：`reanalysis/results/tables/kla_ddr_pathway_pie_umap_33groups_v3_all_go/protein_all_go_direct_binary_matrix.csv`",
  "- 注释与IDF：`reanalysis/results/tables/kla_ddr_pathway_pie_umap_33groups_v3_all_go/all_go_term_dictionary_and_idf.csv`",
  "- 参数：`reanalysis/results/tables/kla_ddr_pathway_pie_umap_33groups_v3_all_go/umap_v3_parameters.csv`"
)
writeLines(report_lines, report_path, useBytes = TRUE)

message("V3 UMAP proteins: ", nrow(v3_coordinates))
message("Direct protein-GO annotations: ", nrow(go_long))
message("Unique direct GO terms: ", uniqueN(go_long$GO_ID))
message("Unique direct GO profiles: ", uniqueN(signature_table$DirectGOSignature))
message("V3 overlap pairs at pie diameter: ", v3_overlap_pairs)
message("Hollow negative wedges: ", sum(assignment_long$Score == -1L))
message("V3 pie UMAP: ", pie_paths[[1L]])
message("V3 combined figure: ", combined_paths[[1L]])
message("V3 report: ", report_path)
