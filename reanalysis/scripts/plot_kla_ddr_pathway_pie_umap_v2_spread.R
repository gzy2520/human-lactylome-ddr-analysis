#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
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
    "reanalysis/scripts/plot_kla_ddr_pathway_pie_umap_v2_spread.R",
    mustWork = TRUE
  )
}
project_root <- normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = TRUE)

workbook_path <- file.path(
  project_root,
  "data/identifier/260810乳酸化DDR基因评分表.xlsx"
)
raw_go_matrix_path <- file.path(
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
table_dir <- file.path(
  project_root,
  "reanalysis/results/tables/kla_ddr_pathway_pie_umap_33groups_v2_spread"
)
figure_dir <- file.path(
  project_root,
  "reanalysis/results/figures/kla_ddr_pathway_pie_umap_33groups_v2_spread"
)
report_path <- file.path(
  project_root,
  "reanalysis/reports/UMAP_PATHWAY_PIE_33GROUP_V2_SPREAD.md"
)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

required_inputs <- c(workbook_path, raw_go_matrix_path, v1_coordinates_path)
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
  # NPG palette inspired by plots in Nature Reviews Cancer, with nine colors
  # selected from the ten-color set to avoid using two closely related reds.
  ColorName = c(
    "Chambray",
    "Cinnabar",
    "PersianGreen",
    "Shakespeare",
    "Apricot",
    "WildBlueYonder",
    "MonteCarlo",
    "RomanCoffee",
    "Sandrift"
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
    "#B09C85"
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

raw_go_binary <- fread(raw_go_matrix_path)
v1_coordinates <- fread(v1_coordinates_path)

stop_if_false(nrow(scoring) == 507L, "The workbook does not contain 507 populated scoring rows.")
stop_if_false(
  uniqueN(scoring$BaseAccession) == 507L,
  "Scoring BaseAccessions are not unique."
)
stop_if_false(
  nrow(raw_go_binary) == 507L && uniqueN(raw_go_binary$BaseAccession) == 507L,
  "The raw-GO binary matrix does not contain 507 unique proteins."
)
stop_if_false(
  ncol(raw_go_binary) == 67L,
  "The raw-GO binary matrix is not 507 proteins by 66 GO terms."
)
stop_if_false(
  setequal(scoring$BaseAccession, raw_go_binary$BaseAccession),
  "The scoring proteins and raw-GO matrix proteins do not match exactly."
)
stop_if_false(
  setequal(v1_coordinates$BaseAccession, raw_go_binary$BaseAccession),
  "The V1 coordinates and raw-GO matrix proteins do not match exactly."
)

raw_go_binary <- raw_go_binary[match(sort(BaseAccession), BaseAccession)]
scoring <- scoring[match(raw_go_binary$BaseAccession, BaseAccession)]
v1_coordinates <- v1_coordinates[match(raw_go_binary$BaseAccession, BaseAccession)]

raw_go_matrix <- as.matrix(raw_go_binary[, -1L])
storage.mode(raw_go_matrix) <- "numeric"
stop_if_false(all(raw_go_matrix %in% c(0, 1)), "The raw-GO UMAP matrix is not binary.")
stop_if_false(all(rowSums(raw_go_matrix) >= 1L), "At least one protein lacks a GO term.")
stop_if_false(sum(raw_go_matrix) == 1029L, "Expected 1,029 binary protein-GO hits.")

score_matrix <- as.matrix(scoring[, ..pathway_columns])
stop_if_false(!anyNA(score_matrix), "The 507 scoring rows contain missing values.")
stop_if_false(
  all(score_matrix %in% c(-1L, 0L, 1L)),
  "Pathway scores outside {-1, 0, 1} were detected."
)

# V2 prioritizes readable pie placement while retaining the same unsupervised
# raw-GO input. Pathway scores are not included in the embedding.
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
v2_matrix <- uwot::umap(
  X = raw_go_matrix,
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
  identical(dim(v2_matrix), c(507L, 2L)) && all(is.finite(v2_matrix)),
  "V2 UMAP did not produce 507 finite two-dimensional coordinates."
)

v2_coordinates <- data.table(
  BaseAccession = raw_go_binary$BaseAccession,
  UMAP_1 = v2_matrix[, 1L],
  UMAP_2 = v2_matrix[, 2L]
)
plot_data <- merge(
  v2_coordinates,
  scoring,
  by = "BaseAccession",
  all.x = TRUE,
  sort = FALSE
)
setorder(plot_data, BaseAccession)

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
  nearest_neighbor_metrics(v1_coordinates, "V1_fixed_compact"),
  nearest_neighbor_metrics(v2_coordinates, "V2_refitted_spread")
)

unassigned_data <- plot_data[TotalAssignmentCount == 0L]

pie_plot <- ggplot() +
  geom_point(
    data = unassigned_data,
    aes(x = UMAP_1, y = UMAP_2, shape = "No scored assignment"),
    fill = "#D9D9D9",
    color = "#8C8C8C",
    size = 2.1,
    stroke = 0.35
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
    values = c("No scored assignment" = 21),
    name = NULL
  ) +
  coord_equal() +
  labs(
    title = "Kla-DDR pathway assignments on a readability-optimized UMAP",
    subtitle = paste0(
      "507 proteins | raw GO-term input only | expanded local spacing for pie readability"
    ),
    x = "UMAP 1",
    y = "UMAP 2",
    caption = paste0(
      "Filled sectors: promoting score +1. White sectors with pathway-colored borders: ",
      "suppressing score -1.\n",
      "Sector area denotes membership, not abundance. Small positional differences ",
      "among identical GO profiles are display-only."
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
  id.vars = c(
    "Pathway", "DisplayLabel", "ColorName", "Color", "PathwayOrder"
  ),
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
  grDevices::svg(svg_path, width = width, height = height, onefile = FALSE, bg = "white")
  print(plot)
  grDevices::dev.off()
  invisible(c(png_path, pdf_path, svg_path))
}

pie_paths <- save_plot_formats(
  pie_plot,
  "kla_ddr_pathway_hollow_negative_pie_umap_v2_spread",
  width = 10.6,
  height = 7.8
)
combined_paths <- save_plot_formats(
  combined_plot,
  "kla_ddr_pathway_hollow_negative_umap_and_summary_v2_spread",
  width = 10.6,
  height = 11.0
)

parameter_table <- data.table(
  Parameter = c(
    "AnalysisUnit",
    "ProteinAnalysisKey",
    "UMAPInput",
    "PathwayScoresUsedInUMAP",
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
    "Palette",
    "PaletteSource"
  ),
  Value = c(
    "one unique protein per pie",
    "isoform-stripped UniProt BaseAccession",
    "507 x 66 protein x raw GO term binary matrix only",
    "FALSE",
    "cosine",
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
    "NPG nine-color selection",
    "https://github.com/nanxstats/ggsci/blob/master/R/palettes.R"
  )
)

input_audit <- data.table(
  InputRole = c(
    "pathway scoring workbook",
    "507 x 66 raw-GO binary UMAP input",
    "V1 coordinate comparison only"
  ),
  Path = sub(
    paste0("^", gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", project_root), "/?"),
    "",
    required_inputs
  ),
  MD5 = unname(tools::md5sum(required_inputs))
)

fwrite(v2_coordinates, file.path(table_dir, "umap_coordinates_v2_spread.csv"))
fwrite(
  plot_data,
  file.path(table_dir, "pathway_umap_plot_data_v2_spread.csv")
)
fwrite(
  assignment_long,
  file.path(table_dir, "pathway_assignment_long_v2.csv")
)
fwrite(
  wedge_polygons,
  file.path(table_dir, "hollow_negative_wedge_polygons.csv")
)
fwrite(summary_counts, file.path(table_dir, "pathway_assignment_summary.csv"))
fwrite(overlap_comparison, file.path(table_dir, "pie_overlap_comparison_v1_v2.csv"))
fwrite(parameter_table, file.path(table_dir, "umap_v2_parameters.csv"))
fwrite(pathway_info, file.path(table_dir, "npg_pathway_color_key.csv"))
fwrite(input_audit, file.path(table_dir, "input_file_audit.csv"))
writeLines(
  capture.output(sessionInfo()),
  file.path(table_dir, "session_info.txt"),
  useBytes = TRUE
)

v1_overlap_pairs <- overlap_comparison[
  Version == "V1_fixed_compact",
  PairDistancesCloserThanPieDiameter
]
v2_overlap_pairs <- overlap_comparison[
  Version == "V2_refitted_spread",
  PairDistancesCloserThanPieDiameter
]

report_lines <- c(
  "# 33组Kla∩DDR通路饼图UMAP V2：扩展间距与负向空心扇形",
  "",
  "## 本版修改",
  "",
  "- 不再使用V1固定坐标；从同一507 × 66原始GO-term二值矩阵重新拟合UMAP。",
  "- 通路评分仍不参与UMAP，仅在坐标生成后用于绘制饼图。",
  "- 新参数：cosine距离、`n_neighbors = 15`、`min_dist = 3`、`spread = 10`、`repulsion_strength = 1.8`、`negative_sample_rate = 10`、750轮、随机种子25、单线程。",
  paste0(
    "- 饼图半径为", pie_radius, "坐标单位。以该直径计算，V1中有",
    v1_overlap_pairs, "对蛋白距离小于一个饼图直径；V2降为",
    v2_overlap_pairs, "对。"
  ),
  "- 102个仅命中`GO:0006974`的蛋白拥有完全相同的输入特征；它们在V2中的相对摊开仅用于显示，不能解释为生物学距离。",
  "",
  "## 扇形编码",
  "",
  "- `+1`促进性评分：使用对应通路颜色的实心扇形，扇形之间以细白线分隔。",
  "- `-1`抑制性评分：使用白色空心扇形，并以对应通路颜色描边。",
  "- `0`：不绘制该通路扇形；9列全为0的22个蛋白显示为灰色圆点。",
  "- 每个非零通路扇形等权；扇形面积不表示表达量、蛋白丰度或证据强弱。",
  "",
  "## 配色与版式",
  "",
  "- 使用NPG离散色板中9个颜色，该色板源自Nature Reviews Cancer图形配色。",
  "- 调色板来源：`https://github.com/nanxstats/ggsci/blob/master/R/palettes.R`。",
  "- 图中文字使用统一无衬线字体，并同时导出600 dpi PNG及PDF/SVG矢量文件，便于按Nature双栏宽度缩放。",
  "- Nature作者指南：`https://www.nature.com/nature/for-authors/final-submission`。",
  "",
  "## 输出",
  "",
  "- 推荐饼图：`reanalysis/results/figures/kla_ddr_pathway_pie_umap_33groups_v2_spread/kla_ddr_pathway_hollow_negative_pie_umap_v2_spread.{png,pdf,svg}`",
  "- 组合图：`reanalysis/results/figures/kla_ddr_pathway_pie_umap_33groups_v2_spread/kla_ddr_pathway_hollow_negative_umap_and_summary_v2_spread.{png,pdf,svg}`",
  "- V2坐标：`reanalysis/results/tables/kla_ddr_pathway_pie_umap_33groups_v2_spread/umap_coordinates_v2_spread.csv`",
  "- 参数：`reanalysis/results/tables/kla_ddr_pathway_pie_umap_33groups_v2_spread/umap_v2_parameters.csv`",
  "- 重叠比较：`reanalysis/results/tables/kla_ddr_pathway_pie_umap_33groups_v2_spread/pie_overlap_comparison_v1_v2.csv`"
)
writeLines(report_lines, report_path, useBytes = TRUE)

message("V2 UMAP proteins: ", nrow(v2_coordinates))
message("V1 overlap pairs at pie diameter: ", v1_overlap_pairs)
message("V2 overlap pairs at pie diameter: ", v2_overlap_pairs)
message("Hollow negative wedges: ", sum(assignment_long$Score == -1L))
message("V2 pie UMAP: ", pie_paths[[1L]])
message("V2 combined figure: ", combined_paths[[1L]])
message("V2 report: ", report_path)
