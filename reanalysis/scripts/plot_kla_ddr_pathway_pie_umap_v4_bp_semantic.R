#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(GO.db)
  library(ggplot2)
  library(Matrix)
  library(patchwork)
  library(uwot)
})

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
script_path <- if (length(file_arg) > 0L) {
  normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE)
} else {
  normalizePath(
    "reanalysis/scripts/plot_kla_ddr_pathway_pie_umap_v4_bp_semantic.R",
    mustWork = TRUE
  )
}
project_root <- normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = TRUE)
variant_tag <- trimws(Sys.getenv("KLA_V4_VARIANT_TAG", ""))
variant_tag <- gsub("[^A-Za-z0-9_-]+", "_", variant_tag)
variant_suffix <- if (nzchar(variant_tag)) paste0("_", variant_tag) else ""

v3_table_dir <- file.path(
  project_root,
  "reanalysis/results/tables/kla_ddr_pathway_pie_umap_33groups_v3_all_go"
)
go_long_path <- file.path(v3_table_dir, "uniprot_direct_go_annotation_long.csv")
full_go_matrix_path <- file.path(
  v3_table_dir,
  "protein_all_go_direct_binary_matrix.csv"
)
plot_data_path <- file.path(
  v3_table_dir,
  "pathway_umap_plot_data_v3_all_go.csv"
)
assignment_path <- file.path(v3_table_dir, "pathway_assignment_long_v3.csv")
color_path <- file.path(v3_table_dir, "pathway_color_key_v3.csv")
v3_coordinates_path <- file.path(v3_table_dir, "umap_coordinates_v3_all_go.csv")
uniprot_path <- file.path(
  project_root,
  "reanalysis/config/uniprot_kla_ddr_507_all_go_2026-08-10.tsv"
)
uniprot_metadata_path <- file.path(
  project_root,
  "reanalysis/config/uniprot_kla_ddr_507_all_go_2026-08-10_metadata.tsv"
)
table_dir <- file.path(
  project_root,
  paste0(
    "reanalysis/results/tables/kla_ddr_pathway_pie_umap_33groups_v4_bp_semantic",
    variant_suffix
  )
)
figure_dir <- file.path(
  project_root,
  paste0(
    "reanalysis/results/figures/kla_ddr_pathway_pie_umap_33groups_v4_bp_semantic",
    variant_suffix
  )
)
report_path <- file.path(
  project_root,
  paste0(
    "reanalysis/reports/UMAP_PATHWAY_PIE_33GROUP_V4_BP_SEMANTIC",
    toupper(variant_suffix),
    ".md"
  )
)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

required_inputs <- c(
  uniprot_path,
  uniprot_metadata_path,
  go_long_path,
  full_go_matrix_path,
  plot_data_path,
  assignment_path,
  color_path,
  v3_coordinates_path
)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs) > 0L) {
  stop("Missing required input(s): ", paste(missing_inputs, collapse = "; "))
}

assert <- function(condition, message) {
  if (!isTRUE(condition)) {
    stop(message, call. = FALSE)
  }
}

relative_path <- function(path) {
  sub(
    paste0("^", gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", project_root), "/?"),
    "",
    path
  )
}

go_long <- fread(go_long_path)
plot_source <- fread(plot_data_path)
assignment_long <- fread(assignment_path)
pathway_info <- fread(color_path)
v3_coordinates <- fread(v3_coordinates_path)
uniprot_metadata <- fread(uniprot_metadata_path, sep = "\t")

protein_levels <- sort(unique(go_long$BaseAccession))
assert(length(protein_levels) == 507L, "Expected 507 proteins in the direct-GO table.")
assert(
  setequal(protein_levels, plot_source$BaseAccession) &&
    all(assignment_long$BaseAccession %in% protein_levels) &&
    uniqueN(assignment_long$BaseAccession) == 485L,
  "The GO, pathway-score, and nonzero-assignment protein sets do not match."
)
assert(
  nrow(go_long) == 13738L && uniqueN(go_long$GO_ID) == 3461L,
  "The saved full direct-GO table does not match the V3 audit."
)
assert(
  nrow(pathway_info) == 9L &&
    pathway_info[DisplayLabel == "Other support", Color] == "#6F6F6F",
  "The corrected nine-color pathway key is missing."
)

# Biological-process annotations define the functional pathway layout. All
# direct BP/CC/MF annotations remain preserved in the V3 full-GO tables.
bp_direct <- go_long[Ontology == "BP"]
assert(
  uniqueN(bp_direct$BaseAccession) == 507L,
  "At least one protein lacks a direct BP annotation."
)

bp_ancestor_map <- as.list(GOBPANCESTOR)
bp_root <- "GO:0008150"
bp_terms_missing_from_go_db <- bp_direct[
  !GO_ID %in% names(bp_ancestor_map),
  .(
    DirectAnnotationCount = .N,
    AnnotatedProteinCount = uniqueN(BaseAccession),
    UniProtGO_Term = paste(sort(unique(GO_Term)), collapse = ";")
  ),
  by = GO_ID
][order(GO_ID)]
bp_expansion <- rbindlist(
  lapply(seq_len(nrow(bp_direct)), function(index) {
    direct_id <- bp_direct$GO_ID[[index]]
    ancestors <- bp_ancestor_map[[direct_id]]
    feature_ids <- unique(c(direct_id, ancestors))
    feature_ids <- feature_ids[
      !is.na(feature_ids) &
        feature_ids != "all" &
        feature_ids != bp_root
    ]
    data.table(
      BaseAccession = bp_direct$BaseAccession[[index]],
      DirectGO_ID = direct_id,
      FeatureGO_ID = feature_ids,
      FeatureSource = ifelse(
        feature_ids == direct_id,
        "Direct",
        "Ancestor"
      )
    )
  }),
  use.names = TRUE
)
bp_expansion <- unique(
  bp_expansion,
  by = c("BaseAccession", "DirectGO_ID", "FeatureGO_ID")
)
setorder(bp_expansion, BaseAccession, DirectGO_ID, FeatureGO_ID)

bp_feature_long <- unique(
  bp_expansion[, .(BaseAccession, FeatureGO_ID)],
  by = c("BaseAccession", "FeatureGO_ID")
)
bp_feature_frequency <- bp_feature_long[
  ,
  .(ProteinCount = uniqueN(BaseAccession)),
  by = FeatureGO_ID
]
minimum_protein_count <- 2L
maximum_protein_count <- floor(0.80 * length(protein_levels))
eligible_features <- bp_feature_frequency[
  ProteinCount >= minimum_protein_count &
    ProteinCount <= maximum_protein_count
]
bp_semantic_long <- bp_feature_long[
  eligible_features,
  on = "FeatureGO_ID",
  nomatch = 0L
]
setorder(bp_semantic_long, BaseAccession, FeatureGO_ID)

feature_levels <- sort(unique(bp_semantic_long$FeatureGO_ID))
bp_semantic_binary <- sparseMatrix(
  i = match(bp_semantic_long$BaseAccession, protein_levels),
  j = match(bp_semantic_long$FeatureGO_ID, feature_levels),
  x = 1,
  dims = c(length(protein_levels), length(feature_levels)),
  dimnames = list(protein_levels, feature_levels)
)
bp_semantic_binary@x[] <- 1

assert(
  nrow(bp_semantic_binary) == 507L &&
    ncol(bp_semantic_binary) > 1000L,
  "BP semantic matrix is incomplete."
)
assert(
  all(Matrix::rowSums(bp_semantic_binary) >= 1L),
  "At least one protein lacks an eligible shared BP semantic feature."
)

umap_parameters <- list(
  n_neighbors = 8L,
  min_dist = 0.2,
  spread = 1.5,
  repulsion_strength = 1.5,
  negative_sample_rate = 10L,
  n_epochs = 1000L,
  seed = 25L,
  initialization = "random"
)

set.seed(umap_parameters$seed)
raw_umap_matrix <- uwot::umap(
  X = as.matrix(bp_semantic_binary),
  n_neighbors = umap_parameters$n_neighbors,
  n_components = 2L,
  metric = "cosine",
  n_epochs = umap_parameters$n_epochs,
  scale = FALSE,
  init = umap_parameters$initialization,
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
assert(
  identical(dim(raw_umap_matrix), c(507L, 2L)) &&
    all(is.finite(raw_umap_matrix)),
  "BP semantic UMAP did not produce 507 finite two-dimensional coordinates."
)

raw_coordinates <- data.table(
  BaseAccession = protein_levels,
  UMAP_1 = raw_umap_matrix[, 1L],
  UMAP_2 = raw_umap_matrix[, 2L]
)

# Use a uniform scale followed by deterministic local collision resolution.
# The raw UMAP coordinates are retained as the inferential embedding; the
# adjusted coordinates are used only as pie centers.
scaled_matrix <- scale(
  raw_umap_matrix,
  center = colMeans(raw_umap_matrix),
  scale = FALSE
)
display_target_range <- as.numeric(
  Sys.getenv("KLA_V4_DISPLAY_TARGET_RANGE", "100")
)
assert(
  length(display_target_range) == 1L &&
    is.finite(display_target_range) &&
    display_target_range > 0,
  "KLA_V4_DISPLAY_TARGET_RANGE must be one positive finite number."
)
uniform_scale_factor <- display_target_range / max(
  diff(range(scaled_matrix[, 1L])),
  diff(range(scaled_matrix[, 2L]))
)
scaled_matrix <- scaled_matrix * uniform_scale_factor

pie_radius <- 0.75
minimum_center_separation <- as.numeric(
  Sys.getenv("KLA_V4_MINIMUM_CENTER_SEPARATION", "1.50")
)
assert(
  length(minimum_center_separation) == 1L &&
    is.finite(minimum_center_separation) &&
    minimum_center_separation >= 2 * pie_radius,
  "KLA_V4_MINIMUM_CENTER_SEPARATION must be at least the pie diameter (1.50)."
)

resolve_circle_collisions <- function(
  coordinates,
  minimum_separation,
  max_iterations = 2500L
) {
  result <- coordinates
  point_count <- nrow(result)
  pair_indices <- which(
    upper.tri(matrix(FALSE, point_count, point_count)),
    arr.ind = TRUE
  )
  completed_iterations <- 0L
  for (iteration in seq_len(max_iterations)) {
    delta_x <- result[pair_indices[, 2L], 1L] -
      result[pair_indices[, 1L], 1L]
    delta_y <- result[pair_indices[, 2L], 2L] -
      result[pair_indices[, 1L], 2L]
    pair_distance <- sqrt(delta_x^2 + delta_y^2)
    overlapping <- which(pair_distance < minimum_separation)
    if (length(overlapping) == 0L) {
      completed_iterations <- iteration - 1L
      break
    }
    ordered_overlaps <- overlapping[order(pair_distance[overlapping])]
    for (pair_index in ordered_overlaps) {
      first <- pair_indices[pair_index, 1L]
      second <- pair_indices[pair_index, 2L]
      vector <- result[second, ] - result[first, ]
      distance <- sqrt(sum(vector^2))
      if (distance >= minimum_separation) {
        next
      }
      if (distance < 1e-10) {
        angle <- ((first * 37L + second * 101L) %% 360L) * pi / 180
        unit_vector <- c(cos(angle), sin(angle))
      } else {
        unit_vector <- vector / distance
      }
      push <- (minimum_separation - distance) / 2 + 1e-6
      result[first, ] <- result[first, ] - push * unit_vector
      result[second, ] <- result[second, ] + push * unit_vector
    }
    completed_iterations <- iteration
  }
  attr(result, "iterations") <- completed_iterations
  result
}

display_matrix <- resolve_circle_collisions(
  scaled_matrix,
  minimum_center_separation
)
collision_iterations <- attr(display_matrix, "iterations")

display_coordinates <- data.table(
  BaseAccession = protein_levels,
  UMAP_Display_1 = display_matrix[, 1L],
  UMAP_Display_2 = display_matrix[, 2L]
)
coordinate_audit <- merge(
  raw_coordinates,
  display_coordinates,
  by = "BaseAccession",
  sort = FALSE
)
coordinate_audit[
  ,
  `:=`(
    UMAP_Scaled_1 = scaled_matrix[, 1L],
    UMAP_Scaled_2 = scaled_matrix[, 2L]
  )
]
coordinate_audit[
  ,
  DisplayDisplacement := sqrt(
    (UMAP_Display_1 - UMAP_Scaled_1)^2 +
      (UMAP_Display_2 - UMAP_Scaled_2)^2
  )
]

nearest_neighbor_summary <- function(x, y, label, threshold) {
  matrix_xy <- cbind(x, y)
  distances <- as.matrix(stats::dist(matrix_xy))
  diag(distances) <- Inf
  nearest <- apply(distances, 1L, min)
  upper_distances <- distances[upper.tri(distances)]
  data.table(
    CoordinateSet = label,
    Threshold = threshold,
    MinimumNearestNeighborDistance = min(nearest),
    P10NearestNeighborDistance = unname(quantile(nearest, 0.10)),
    MedianNearestNeighborDistance = median(nearest),
    PairDistancesBelowThreshold = sum(upper_distances < threshold),
    ProteinsWithNearestNeighborBelowThreshold = sum(nearest < threshold),
    XRange = diff(range(x)),
    YRange = diff(range(y))
  )
}

spacing_summary <- rbind(
  nearest_neighbor_summary(
    coordinate_audit$UMAP_Scaled_1,
    coordinate_audit$UMAP_Scaled_2,
    "Uniformly scaled raw UMAP",
    minimum_center_separation
  ),
  nearest_neighbor_summary(
    coordinate_audit$UMAP_Display_1,
    coordinate_audit$UMAP_Display_2,
    "Collision-resolved display centers",
    minimum_center_separation
  )
)
raw_pair_distances <- as.numeric(stats::dist(raw_umap_matrix))
display_pair_distances <- as.numeric(stats::dist(display_matrix))
pair_distance_spearman <- cor(
  raw_pair_distances,
  display_pair_distances,
  method = "spearman"
)

plot_data <- merge(
  display_coordinates,
  plot_source[, !c("UMAP_1", "UMAP_2"), with = FALSE],
  by = "BaseAccession",
  all.x = TRUE,
  sort = FALSE
)
setorder(plot_data, BaseAccession)
assert(
  plot_data[TotalAssignmentCount == 0L, .N] == 22L,
  "Expected 22 proteins without a scored pathway/function."
)
assert(
  nrow(assignment_long) == 1175L &&
    sum(assignment_long$Score == -1L) == 67L,
  "Pathway assignment counts changed unexpectedly."
)

arc_points <- 20L
make_wedges <- function(accession, center_x, center_y) {
  assignments <- assignment_long[
    BaseAccession == accession
  ][order(PathwayOrder)]
  sector_count <- nrow(assignments)
  if (sector_count == 0L) {
    return(NULL)
  }
  rbindlist(lapply(seq_len(sector_count), function(sector_index) {
    start_angle <- pi / 2 - 2 * pi * (sector_index - 1L) / sector_count
    end_angle <- pi / 2 - 2 * pi * sector_index / sector_count
    theta <- seq(start_angle, end_angle, length.out = arc_points)
    assignment <- assignments[sector_index]
    data.table(
      BaseAccession = accession,
      WedgeID = paste(accession, sector_index, sep = "__"),
      Pathway = assignment$Pathway,
      Score = assignment$Score,
      X = c(center_x, center_x + pie_radius * cos(theta), center_x),
      Y = c(center_y, center_y + pie_radius * sin(theta), center_y),
      FillKey = if (assignment$Score == 1L) {
        assignment$Pathway
      } else {
        "NegativeHollow"
      },
      BorderKey = if (assignment$Score == 1L) {
        assignment$Pathway
      } else {
        assignment$Pathway
      }
    )
  }))
}

wedge_polygons <- rbindlist(
  lapply(seq_len(nrow(plot_data)), function(index) {
    make_wedges(
      plot_data$BaseAccession[[index]],
      plot_data$UMAP_Display_1[[index]],
      plot_data$UMAP_Display_2[[index]]
    )
  })
)
assert(
  uniqueN(wedge_polygons$WedgeID) == 1175L,
  "Wedge construction did not yield one sector per nonzero assignment."
)

fill_values <- c(
  setNames(pathway_info$Color, pathway_info$Pathway)
)
border_values <- setNames(pathway_info$Color, pathway_info$Pathway)
unassigned_data <- plot_data[TotalAssignmentCount == 0L]

pie_plot <- ggplot() +
  geom_point(
    data = unassigned_data,
    aes(
      x = UMAP_Display_1,
      y = UMAP_Display_2,
      shape = "No scored pathway/function"
    ),
    color = "#A6A6A6",
    size = 2.85,
    stroke = 0.50
  ) +
  geom_polygon(
    data = wedge_polygons[Score == 1L],
    aes(
      x = X,
      y = Y,
      group = WedgeID,
      fill = Pathway,
      color = Pathway
    ),
    linewidth = 0.42,
    linejoin = "round"
  ) +
  geom_polygon(
    data = wedge_polygons[Score == -1L],
    aes(x = X, y = Y, group = WedgeID, color = Pathway),
    fill = "#FFFFFF",
    linewidth = 0.46,
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
    title = "Kla-DDR pathway map derived from a BP semantic UMAP",
    subtitle = paste0(
      "507 proteins | direct UniProt BP plus available GO.db ancestors | ",
      ncol(bp_semantic_binary), " shared features"
    ),
    x = "UMAP display 1",
    y = "UMAP display 2",
    caption = paste0(
      "Filled sectors: promoting +1; white sectors with pathway-colored borders: ",
      "suppressing -1; gray circled crosses: no scored pathway/function.\n",
      "Raw UMAP used BP semantic features only. Pie centers were collision-resolved ",
      "for readability; use the separately saved raw coordinates for local topology."
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

raw_point_plot <- ggplot(
  raw_coordinates,
  aes(x = UMAP_1, y = UMAP_2)
) +
  geom_point(
    color = "#3C5488",
    fill = "#FFFFFF",
    shape = 21,
    size = 1.7,
    stroke = 0.42,
    alpha = 0.90
  ) +
  coord_equal() +
  labs(
    title = "Raw BP semantic UMAP before pie collision adjustment",
    subtitle = paste0(
      "n_neighbors = ", umap_parameters$n_neighbors,
      " | min_dist = ", umap_parameters$min_dist,
      " | cosine distance"
    ),
    x = "UMAP 1",
    y = "UMAP 2",
    caption = "This is the inferential embedding; no display displacement is applied."
  ) +
  theme_classic(base_size = 10, base_family = "sans") +
  theme(
    plot.title = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 8.5, color = "#4D4D4D"),
    plot.caption = element_text(size = 7.5, hjust = 0, color = "#5A5A5A"),
    axis.text = element_text(size = 8, color = "#333333"),
    axis.title = element_text(size = 9)
  )

coordinate_comparison_plot <- raw_point_plot + pie_plot +
  plot_layout(widths = c(1, 1.65)) +
  plot_annotation(
    title = "Raw functional embedding and collision-resolved pathway display"
  )

save_formats <- function(plot, stem, width, height) {
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
    bg = "white"
  )
  ggsave(
    pdf_path,
    plot,
    width = width,
    height = height,
    units = "in",
    bg = "white"
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

pie_paths <- save_formats(
  pie_plot,
  "kla_ddr_pathway_pie_umap_v4_bp_semantic",
  width = 10.6,
  height = 7.8
)
raw_paths <- save_formats(
  raw_point_plot,
  "kla_ddr_raw_bp_semantic_umap_v4",
  width = 7.2,
  height = 6.2
)
comparison_paths <- save_formats(
  coordinate_comparison_plot,
  "kla_ddr_raw_and_pie_bp_semantic_umap_v4",
  width = 15.5,
  height = 7.8
)

go_term_name <- function(go_id) {
  term_object <- GOTERM[[go_id]]
  if (is.null(term_object)) {
    return(NA_character_)
  }
  Term(term_object)
}
bp_feature_dictionary <- copy(bp_feature_frequency)
bp_feature_dictionary[
  ,
  GO_Term := vapply(FeatureGO_ID, go_term_name, character(1L))
]
bp_feature_dictionary[
  ,
  `:=`(
    EligibleForUMAP = FeatureGO_ID %in% eligible_features$FeatureGO_ID,
    InclusionReason = fifelse(
      ProteinCount < minimum_protein_count,
      "excluded: singleton feature cannot define between-protein similarity",
      fifelse(
        ProteinCount > maximum_protein_count,
        "excluded: near-root feature annotated to >80% of proteins",
        "included: shared BP semantic feature"
      )
    )
  )
]
setorder(bp_feature_dictionary, -ProteinCount, FeatureGO_ID)

bp_binary_dense <- as.data.table(as.matrix(bp_semantic_binary))
bp_binary_dense[, BaseAccession := protein_levels]
setcolorder(bp_binary_dense, c("BaseAccession", feature_levels))

method_parameters <- data.table(
  Parameter = c(
    "ProteinSet",
    "ProteinAnalysisKey",
    "FullGOArchive",
    "EmbeddingGOAspect",
    "HierarchyExpansion",
    "ExcludedRoot",
    "MinimumFeatureProteinCount",
    "MaximumFeatureProteinCount",
    "DistanceMetric",
    "NNeighbors",
    "MinDist",
    "Spread",
    "RepulsionStrength",
    "NegativeSampleRate",
    "NEpochs",
    "Initialization",
    "RandomSeed",
    "PathwayScoresUsedInUMAP",
    "SampleDetectionUsedInUMAP",
    "DisplayTargetRange",
    "PieRadius",
    "MinimumCenterSeparation",
    "CollisionIterations",
    "RawDisplayPairDistanceSpearman"
  ),
  Value = as.character(c(
    "507 Kla-intersection-DDR proteins from the fixed 33-group scope",
    "isoform-stripped UniProt BaseAccession",
    "all direct BP/CC/MF annotations retained in V3 tables",
    "BP",
    "each direct BP term plus all available GO.db BP ancestors",
    bp_root,
    minimum_protein_count,
    maximum_protein_count,
    "cosine on binary shared BP semantic features",
    umap_parameters$n_neighbors,
    umap_parameters$min_dist,
    umap_parameters$spread,
    umap_parameters$repulsion_strength,
    umap_parameters$negative_sample_rate,
    umap_parameters$n_epochs,
    umap_parameters$initialization,
    umap_parameters$seed,
    "FALSE",
    "FALSE",
    display_target_range,
    pie_radius,
    minimum_center_separation,
    collision_iterations,
    format(pair_distance_spearman, digits = 8)
  ))
)

input_audit <- data.table(
  InputRole = c(
    "raw UniProt full-GO response",
    "UniProt retrieval metadata",
    "saved direct protein-GO long table",
    "saved full 507 x 3,461 direct-GO matrix",
    "pathway score/plot data",
    "signed pathway assignment long table",
    "corrected pathway color key",
    "V3 all-aspect coordinate comparison"
  ),
  Path = relative_path(required_inputs),
  MD5 = unname(tools::md5sum(required_inputs))
)

fwrite(
  bp_expansion,
  file.path(table_dir, "bp_direct_to_ancestor_expansion_long.csv")
)
fwrite(
  bp_terms_missing_from_go_db,
  file.path(table_dir, "bp_direct_terms_missing_from_go_db.csv")
)
fwrite(
  bp_feature_dictionary,
  file.path(table_dir, "bp_semantic_feature_dictionary.csv")
)
fwrite(
  bp_semantic_long,
  file.path(table_dir, "protein_bp_semantic_feature_long.csv")
)
fwrite(
  bp_binary_dense,
  file.path(table_dir, "protein_bp_semantic_shared_binary_matrix.csv")
)
fwrite(
  raw_coordinates,
  file.path(table_dir, "umap_raw_coordinates_v4_bp_semantic.csv")
)
fwrite(
  display_coordinates,
  file.path(table_dir, "umap_display_coordinates_v4_bp_semantic.csv")
)
fwrite(
  coordinate_audit,
  file.path(table_dir, "umap_raw_scaled_display_coordinate_audit.csv")
)
fwrite(
  spacing_summary,
  file.path(table_dir, "pie_spacing_before_after_collision.csv")
)
fwrite(
  plot_data,
  file.path(table_dir, "pathway_umap_plot_data_v4_bp_semantic.csv")
)
fwrite(
  wedge_polygons,
  file.path(table_dir, "hollow_negative_wedge_polygons_v4.csv")
)
fwrite(
  method_parameters,
  file.path(table_dir, "umap_v4_parameters.csv")
)
fwrite(
  pathway_info,
  file.path(table_dir, "pathway_color_key_v4.csv")
)
fwrite(input_audit, file.path(table_dir, "input_file_audit.csv"))
writeLines(
  capture.output(sessionInfo()),
  file.path(table_dir, "session_info.txt"),
  useBytes = TRUE
)

uniprot_release <- uniprot_metadata[Key == "UniProtRelease", Value]
uniprot_release_date <- uniprot_metadata[Key == "UniProtReleaseDate", Value]
display_overlap_pairs <- spacing_summary[
  CoordinateSet == "Collision-resolved display centers",
  PairDistancesBelowThreshold
]

report_lines <- c(
  "# 33组Kla∩DDR通路饼图UMAP V4：BP层级语义与防碰撞展示",
  "",
  "## V3问题与V4修正",
  "",
  "- V3直接混合BP/CC/MF并使用IDF；3,461个直接GO term中1,864个只命中一个蛋白，稀有term权重使个体区分压过共享功能，因此图上几乎没有清楚的功能簇。",
  "- V4仍完整保留UniProt原始全GO表、13,738条蛋白–GO长表以及507 × 3,461直接GO二值矩阵，但最终通路布局只使用BP，因为当前饼图表达的是DNA修复通路/生物过程。",
  "- CC和MF未被删除，仍保存在全GO归档表中；它们不再参与最终坐标，以免宽泛定位或分子功能跨过程连接不同BP簇。",
  "",
  "## BP语义特征",
  "",
  paste0(
    "- 注释来源：UniProt ", uniprot_release, "（", uniprot_release_date,
    "）；GO层级来源：GO.db ", as.character(packageVersion("GO.db")), "。"
  ),
  "- 每个直接BP term扩展到当前GO.db版本可提供的BP祖先节点，删除BP根`GO:0008150`。",
  paste0(
    "- 有", nrow(bp_terms_missing_from_go_db),
    "个UniProt直接BP ID不在当前GO.db祖先映射中；这些直接term仍保留，",
    "但无法补充祖先。明细保存在`bp_direct_terms_missing_from_go_db.csv`。"
  ),
  paste0(
    "- 只有命中至少", minimum_protein_count, "个、至多",
    maximum_protein_count, "个蛋白的共享特征参与UMAP；单蛋白term仍保存在原始表中，",
    "但不能贡献蛋白间相似性；覆盖>80%蛋白的近根term被排除。"
  ),
  paste0(
    "- 最终UMAP矩阵为507 × ", ncol(bp_semantic_binary),
    "，使用binary cosine、`n_neighbors = 8`、`min_dist = 0.2`、",
    "`spread = 1.5`、1,000轮、种子25、单线程。"
  ),
  "- 通路评分和33组检出模式均不参与UMAP。",
  "",
  "## 饼图防碰撞",
  "",
  "- 原始UMAP坐标单独保存，是功能结构解释使用的坐标。",
  paste0(
    "- 主饼图先对原始坐标统一缩放，再执行确定性防碰撞；饼半径为",
    format(pie_radius, digits = 4),
    "、直径为", format(2 * pie_radius, digits = 4),
    "、最小中心距为", format(minimum_center_separation, digits = 4),
    "，饼边缘最小净间距为",
    format(minimum_center_separation - 2 * pie_radius, digits = 4),
    "。间距越大，展示坐标对原始UMAP局部结构的改变越强。"
  ),
  paste0(
    "- 防碰撞后小于设定中心距的蛋白对为", display_overlap_pairs,
    "；原始与展示坐标全部蛋白对距离的Spearman相关为",
    format(pair_distance_spearman, digits = 4), "。"
  ),
  "- 防碰撞会改变稠密区域中的局部邻居关系；主图只用于读取通路组合和宽尺度几何。局部功能拓扑必须使用另存的原始UMAP坐标/点图解释。",
  "",
  "## 编码与配色",
  "",
  "- `+1`为通路颜色实心扇形，并用完全相同的通路颜色描边以覆盖polygon闭合缝；扇形之间不加白色分隔线。`-1`为白色空心扇形并用通路色描边。",
  "- 22个全0蛋白使用浅灰圈叉；Chromatin interaction为`#7E6148`，Other support为`#6F6F6F`。",
  "",
  "## 保存的GO表和矩阵",
  "",
  "- UniProt原始全GO返回表：`reanalysis/config/uniprot_kla_ddr_507_all_go_2026-08-10.tsv`",
  "- UniProt下载元数据：`reanalysis/config/uniprot_kla_ddr_507_all_go_2026-08-10_metadata.tsv`",
  "- 全部直接蛋白–GO长表：`reanalysis/results/tables/kla_ddr_pathway_pie_umap_33groups_v3_all_go/uniprot_direct_go_annotation_long.csv`",
  "- 507 × 3,461全部直接GO矩阵：`reanalysis/results/tables/kla_ddr_pathway_pie_umap_33groups_v3_all_go/protein_all_go_direct_binary_matrix.csv`",
  "- BP直接term到祖先节点展开表：`reanalysis/results/tables/kla_ddr_pathway_pie_umap_33groups_v4_bp_semantic/bp_direct_to_ancestor_expansion_long.csv`",
  "- V4实际使用的BP语义矩阵：`reanalysis/results/tables/kla_ddr_pathway_pie_umap_33groups_v4_bp_semantic/protein_bp_semantic_shared_binary_matrix.csv`",
  "",
  "## 图形输出",
  "",
  "- 推荐饼图：`reanalysis/results/figures/kla_ddr_pathway_pie_umap_33groups_v4_bp_semantic/kla_ddr_pathway_pie_umap_v4_bp_semantic.{png,pdf,svg}`",
  "- 未位移原始点图：`reanalysis/results/figures/kla_ddr_pathway_pie_umap_33groups_v4_bp_semantic/kla_ddr_raw_bp_semantic_umap_v4.{png,pdf,svg}`",
  "- 原始与展示对照：`reanalysis/results/figures/kla_ddr_pathway_pie_umap_33groups_v4_bp_semantic/kla_ddr_raw_and_pie_bp_semantic_umap_v4.{png,pdf,svg}`"
)
writeLines(report_lines, report_path, useBytes = TRUE)

message("V4 proteins: ", nrow(raw_coordinates))
message("Eligible shared BP semantic features: ", ncol(bp_semantic_binary))
message("Collision iterations: ", collision_iterations)
message("Display overlap pairs: ", display_overlap_pairs)
message("Raw/display pair-distance Spearman: ", format(pair_distance_spearman, digits = 8))
message("V4 pie UMAP: ", pie_paths[[1L]])
message("V4 raw UMAP: ", raw_paths[[1L]])
message("V4 comparison: ", comparison_paths[[1L]])
message("V4 report: ", report_path)
