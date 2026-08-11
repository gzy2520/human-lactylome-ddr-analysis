#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(readxl)
  library(scatterpie)
})

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
script_path <- if (length(file_arg) > 0L) {
  normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE)
} else {
  normalizePath("reanalysis/scripts/plot_kla_ddr_pathway_pie_umap.R", mustWork = TRUE)
}
project_root <- normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = TRUE)

workbook_path <- file.path(
  project_root,
  "data/identifier/260810乳酸化DDR基因评分表.xlsx"
)
coordinates_path <- file.path(
  project_root,
  "reanalysis/results/tables/kla_ddr_raw_go_umap_33groups/umap_coordinates_fixed.csv"
)
raw_go_parameters_path <- file.path(
  project_root,
  "reanalysis/results/tables/kla_ddr_raw_go_umap_33groups/umap_parameters.csv"
)
table_dir <- file.path(
  project_root,
  "reanalysis/results/tables/kla_ddr_pathway_pie_umap_33groups"
)
figure_dir <- file.path(
  project_root,
  "reanalysis/results/figures/kla_ddr_pathway_pie_umap_33groups"
)
report_path <- file.path(
  project_root,
  "reanalysis/reports/UMAP_PATHWAY_PIE_33GROUP_DATA_SCOPE.md"
)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

required_inputs <- c(workbook_path, coordinates_path, raw_go_parameters_path)
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

lighten_color <- function(hex_color, white_fraction = 0.58) {
  rgb_value <- grDevices::col2rgb(hex_color)
  mixed <- rgb_value * (1 - white_fraction) + 255 * white_fraction
  grDevices::rgb(mixed[1L], mixed[2L], mixed[3L], maxColorValue = 255)
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
  Key = c(
    "HR",
    "NHEJ",
    "AEJ",
    "BER",
    "NER",
    "MMR",
    "FA",
    "Chromatin",
    "Other"
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
  Color = c(
    "#0072B2",
    "#D55E00",
    "#CC79A7",
    "#009E73",
    "#E69F00",
    "#56B4E9",
    "#F0E442",
    "#6A3D9A",
    "#7F7F7F"
  )
)
pathway_columns <- pathway_info$Pathway

scoring_raw <- as.data.table(
  read_excel(workbook_path, sheet = "评分表", .name_repair = "minimal")
)
reference_raw <- as.data.table(
  read_excel(workbook_path, sheet = "参考文献表", .name_repair = "minimal")
)
coordinates <- fread(coordinates_path)
raw_go_parameters <- fread(raw_go_parameters_path)

stop_if_false(
  all(
    c("ID", "BaseAccession", "GeneSymbol", "ProteinName", "Note", pathway_columns) %in%
      names(scoring_raw)
  ),
  "The scoring sheet is missing required columns."
)
stop_if_false(
  all(
    c(
      "RefID", "Symbol", "Pathway", "RepresentativeResearchArticle",
      "DOI_or_SourceURL", "Score", "EvidenceNote"
    ) %in% names(reference_raw)
  ),
  "The reference sheet is missing required columns."
)
stop_if_false(
  identical(names(coordinates), c("BaseAccession", "UMAP_1", "UMAP_2")),
  "The fixed-coordinate table has unexpected columns."
)

scoring <- scoring_raw[
  !is.na(BaseAccession) & nzchar(trimws(as.character(BaseAccession)))
]
scoring[, BaseAccession := base_accession(BaseAccession)]
scoring[, ID := as.integer(ID)]
for (column in pathway_columns) {
  set(scoring, j = column, value = as.integer(scoring[[column]]))
}

stop_if_false(nrow(scoring_raw) == 514L, "Expected 514 rows below the scoring-sheet header.")
stop_if_false(
  sum(is.na(scoring_raw$BaseAccession) | !nzchar(trimws(scoring_raw$BaseAccession))) == 7L,
  "Expected seven trailing ID-only template rows in the scoring sheet."
)
stop_if_false(nrow(scoring) == 507L, "Expected 507 populated scoring rows.")
stop_if_false(
  uniqueN(scoring$BaseAccession) == 507L,
  "Scoring-sheet BaseAccessions are not unique after isoform stripping."
)
stop_if_false(
  nrow(coordinates) == 507L && uniqueN(coordinates$BaseAccession) == 507L,
  "The fixed coordinate table does not contain 507 unique BaseAccessions."
)
stop_if_false(
  setequal(scoring$BaseAccession, coordinates$BaseAccession),
  "The 507 scoring-sheet proteins do not exactly match the 507 fixed UMAP coordinates."
)

score_matrix <- as.matrix(scoring[, ..pathway_columns])
stop_if_false(!anyNA(score_matrix), "The 507 populated scoring rows contain missing pathway scores.")
stop_if_false(
  all(score_matrix %in% c(-1L, 0L, 1L)),
  "Pathway scores outside {-1, 0, 1} were detected."
)

score_long <- melt(
  scoring,
  id.vars = c("ID", "BaseAccession", "GeneSymbol", "ProteinName", "Note"),
  measure.vars = pathway_columns,
  variable.name = "Pathway",
  value.name = "Score",
  variable.factor = FALSE
)
assignment_long <- score_long[Score != 0L]
assignment_long[, Direction := fifelse(Score == 1L, "Promoting (+1)", "Suppressing (-1)")]

reference <- copy(reference_raw)
setnames(reference, "RefID", "ID")
reference[, ID := as.integer(ID)]
reference[, Score := as.integer(Score)]

stop_if_false(nrow(assignment_long) == 1175L, "Expected 1,175 nonzero protein-pathway assignments.")
stop_if_false(nrow(reference) == 1175L, "Expected 1,175 pathway-reference rows.")
stop_if_false(
  uniqueN(reference, by = c("ID", "Pathway", "Score")) == 1175L,
  "Reference rows are not unique by ID + Pathway + Score."
)

assignment_with_evidence <- merge(
  assignment_long,
  reference[
    ,
    .(
      ID,
      Pathway,
      Score,
      ReferenceSymbol = Symbol,
      RepresentativeResearchArticle,
      DOI_or_SourceURL,
      EvidenceNote
    )
  ],
  by = c("ID", "Pathway", "Score"),
  all.x = TRUE,
  sort = FALSE
)
stop_if_false(
  nrow(assignment_with_evidence) == 1175L &&
    !anyNA(assignment_with_evidence$DOI_or_SourceURL),
  "At least one nonzero scoring assignment lacks a matching evidence row."
)
stop_if_false(
  setequal(
    paste(reference$ID, reference$Pathway, reference$Score, sep = "\r"),
    paste(
      assignment_with_evidence$ID,
      assignment_with_evidence$Pathway,
      assignment_with_evidence$Score,
      sep = "\r"
    )
  ),
  "The scoring and reference sheets are not in one-to-one assignment agreement."
)

plot_data <- merge(
  coordinates,
  scoring,
  by = "BaseAccession",
  all.x = TRUE,
  sort = FALSE
)
setorder(plot_data, BaseAccession)
stop_if_false(nrow(plot_data) == 507L, "Coordinate-to-score merge changed the protein count.")

for (index in seq_len(nrow(pathway_info))) {
  pathway <- pathway_info$Pathway[[index]]
  key <- pathway_info$Key[[index]]
  plot_data[, (paste0(key, "_Pos")) := as.integer(get(pathway) == 1L)]
  plot_data[, (paste0(key, "_Neg")) := as.integer(get(pathway) == -1L)]
}

positive_columns <- paste0(pathway_info$Key, "_Pos")
negative_columns <- paste0(pathway_info$Key, "_Neg")
pie_columns <- as.vector(rbind(positive_columns, negative_columns))

plot_data[, PositiveAssignmentCount := rowSums(.SD), .SDcols = positive_columns]
plot_data[, NegativeAssignmentCount := rowSums(.SD), .SDcols = negative_columns]
plot_data[
  ,
  TotalAssignmentCount := PositiveAssignmentCount + NegativeAssignmentCount
]
plot_data[
  ,
  AssignmentStatus := fifelse(
    TotalAssignmentCount == 0L,
    "No scored assignment",
    "At least one scored assignment"
  )
]
plot_data[, PieRadius := 0.82]

stop_if_false(
  plot_data[TotalAssignmentCount == 0L, .N] == 22L,
  "Expected 22 proteins with no nonzero pathway assignment."
)
stop_if_false(
  plot_data[NegativeAssignmentCount > 0L, .N] == 53L,
  "Expected 53 proteins with at least one suppressing assignment."
)
stop_if_false(
  sum(plot_data$NegativeAssignmentCount) == 67L,
  "Expected 67 suppressing protein-pathway assignments."
)
stop_if_false(
  plot_data[TotalAssignmentCount > 0L, .N] == 485L,
  "Expected 485 proteins with at least one nonzero assignment."
)

summary_counts <- rbindlist(
  lapply(seq_len(nrow(pathway_info)), function(index) {
    pathway <- pathway_info$Pathway[[index]]
    data.table(
      Pathway = pathway,
      Key = pathway_info$Key[[index]],
      DisplayLabel = pathway_info$DisplayLabel[[index]],
      Positive = sum(scoring[[pathway]] == 1L),
      Negative = sum(scoring[[pathway]] == -1L)
    )
  })
)
summary_counts[, TotalNonzero := Positive + Negative]
summary_counts[, DisplayOrder := seq_len(.N)]

fill_values <- c(
  setNames(pathway_info$Color, positive_columns),
  setNames(
    vapply(pathway_info$Color, lighten_color, character(1L)),
    negative_columns
  )
)
fill_breaks <- positive_columns
fill_labels <- setNames(pathway_info$DisplayLabel, positive_columns)

assigned_data <- plot_data[TotalAssignmentCount > 0L]
unassigned_data <- plot_data[TotalAssignmentCount == 0L]

pie_plot <- ggplot() +
  geom_point(
    data = unassigned_data,
    aes(x = UMAP_1, y = UMAP_2, shape = AssignmentStatus),
    color = "#BDBDBD",
    fill = "#D9D9D9",
    size = 2.2,
    stroke = 0.35
  ) +
  geom_scatterpie(
    data = assigned_data,
    aes(x = UMAP_1, y = UMAP_2, r = PieRadius),
    cols = pie_columns,
    color = "white",
    linewidth = 0.12,
    alpha = 0.97
  ) +
  scale_fill_manual(
    values = fill_values,
    breaks = fill_breaks,
    labels = fill_labels,
    name = "Pathway / function",
    drop = FALSE
  ) +
  scale_shape_manual(
    values = c("No scored assignment" = 21),
    name = NULL
  ) +
  coord_equal() +
  labs(
    title = "Kla-DDR pathway assignments on the fixed raw-GO UMAP",
    subtitle = paste0(
      "33 sample groups | 507 proteins | 485 assigned | 22 without a scored assignment"
    ),
    x = "UMAP 1",
    y = "UMAP 2",
    caption = paste0(
      "Equal pie sectors indicate pathway/function membership, not abundance. ",
      "Solid sectors: score +1; pale sectors: score -1. Coordinates were not refitted."
    )
  ) +
  guides(
    fill = guide_legend(
      override.aes = list(alpha = 1, color = NA),
      order = 1,
      ncol = 1
    ),
    shape = guide_legend(order = 2)
  ) +
  theme_classic(base_size = 11, base_family = "sans") +
  theme(
    axis.line = element_line(linewidth = 0.45, colour = "#333333"),
    axis.ticks = element_line(linewidth = 0.4, colour = "#333333"),
    axis.text = element_text(colour = "#333333"),
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 10, colour = "#4D4D4D"),
    plot.caption = element_text(size = 8.2, colour = "#5A5A5A", hjust = 0),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 8.7),
    legend.key.height = grid::unit(0.36, "cm"),
    plot.margin = margin(10, 10, 8, 10)
  )

bar_long <- melt(
  summary_counts,
  id.vars = c("Pathway", "Key", "DisplayLabel", "DisplayOrder"),
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
bar_long[, DisplayLabel := factor(DisplayLabel, levels = rev(pathway_info$DisplayLabel))]

bar_plot <- ggplot(bar_long, aes(x = PlotCount, y = DisplayLabel, fill = DirectionLabel)) +
  geom_vline(xintercept = 0, linewidth = 0.45, color = "#555555") +
  geom_col(width = 0.68) +
  geom_text(
    data = bar_long[Direction == "Positive"],
    aes(label = Count),
    hjust = -0.25,
    size = 3
  ) +
  geom_text(
    data = bar_long[Direction == "Negative"],
    aes(label = Count),
    hjust = 1.25,
    size = 3
  ) +
  scale_fill_manual(
    values = c(
      "Promoting (+1)" = "#2C7FB8",
      "Suppressing (-1)" = "#B2182B"
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
  theme_classic(base_size = 10.5, base_family = "sans") +
  theme(
    plot.title = element_text(face = "bold", size = 11.5),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(colour = "#333333", size = 9),
    legend.position = "top",
    legend.justification = "left",
    legend.title = element_text(face = "bold", size = 9.5),
    legend.text = element_text(size = 9),
    plot.margin = margin(5, 16, 8, 10)
  )

combined_plot <- pie_plot / bar_plot +
  plot_layout(heights = c(3.35, 1.8))

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
  "kla_ddr_pathway_signed_pie_umap_33groups",
  width = 10.6,
  height = 7.0
)
combined_paths <- save_plot_formats(
  combined_plot,
  "kla_ddr_pathway_umap_and_signed_summary_33groups",
  width = 10.6,
  height = 10.4
)

scoring_output <- scoring[
  ,
  c(
    "ID", "BaseAccession", "GeneSymbol", "ProteinName",
    pathway_columns, "Note"
  ),
  with = FALSE
]
setorder(scoring_output, BaseAccession)
setorder(assignment_with_evidence, BaseAccession, Pathway)

coverage_audit <- data.table(
  Item = c(
    "Fixed UMAP proteins",
    "Populated scoring rows",
    "Exact BaseAccession matches",
    "Proteins with at least one assignment",
    "Proteins with no nonzero assignment",
    "Proteins with at least one suppressing assignment",
    "Promoting protein-pathway assignments",
    "Suppressing protein-pathway assignments",
    "Total nonzero protein-pathway assignments",
    "Evidence rows matched one-to-one",
    "Pathway/function categories"
  ),
  Value = c(
    nrow(coordinates),
    nrow(scoring),
    sum(coordinates$BaseAccession %chin% scoring$BaseAccession),
    plot_data[TotalAssignmentCount > 0L, .N],
    plot_data[TotalAssignmentCount == 0L, .N],
    plot_data[NegativeAssignmentCount > 0L, .N],
    sum(plot_data$PositiveAssignmentCount),
    sum(plot_data$NegativeAssignmentCount),
    nrow(assignment_with_evidence),
    sum(!is.na(assignment_with_evidence$DOI_or_SourceURL)),
    nrow(pathway_info)
  )
)

input_audit <- data.table(
  InputRole = c(
    "pathway scoring and evidence workbook",
    "fixed raw-GO UMAP coordinates",
    "raw-GO UMAP parameters"
  ),
  Path = sub(
    paste0("^", gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", project_root), "/?"),
    "",
    required_inputs
  ),
  MD5 = unname(tools::md5sum(required_inputs))
)

fwrite(scoring_output, file.path(table_dir, "pathway_scores_507.csv"))
fwrite(
  assignment_with_evidence,
  file.path(table_dir, "pathway_assignments_with_evidence_1175.csv")
)
fwrite(summary_counts, file.path(table_dir, "pathway_assignment_summary.csv"))
fwrite(
  plot_data,
  file.path(table_dir, "pathway_umap_plot_data_fixed_coordinates.csv")
)
fwrite(coverage_audit, file.path(table_dir, "pathway_umap_coverage_audit.csv"))
fwrite(input_audit, file.path(table_dir, "input_file_audit.csv"))
fwrite(pathway_info, file.path(table_dir, "pathway_color_key.csv"))
writeLines(
  capture.output(sessionInfo()),
  file.path(table_dir, "session_info.txt"),
  useBytes = TRUE
)

report_lines <- c(
  "# 33组Kla∩DDR通路评分与固定UMAP展示审计",
  "",
  "## 数据合并",
  "",
  "- 评分来源：`data/identifier/260810乳酸化DDR基因评分表.xlsx`的“评分表”工作表。",
  "- 工作表含507个有`BaseAccession`的正式评分行，另有7个仅保留ID的空白模板行；空白行不参与分析。",
  "- 507个评分蛋白与`umap_coordinates_fixed.csv`中的507个蛋白按去isoform的UniProt `BaseAccession`逐一完全匹配。",
  "- `GeneSymbol`和蛋白名称仅保留用于显示与审计，不作为合并、去重或分析键。",
  "- UMAP坐标直接读取既有固定坐标，未重新拟合。",
  "",
  "## 通路评分",
  "",
  "- 9个通路/功能类别：HR、NHEJ、AEJ、BER、NER、MMR、FA、Chromatin Interaction和Other support。",
  "- 分值仅为`-1/0/+1`：`+1`表示促进性证据，`-1`表示抑制性证据，`0`表示未分配该通路。",
  "- 共485个蛋白至少有1个非零分配；22个蛋白9列均为0。",
  "- 共1,175个非零蛋白–通路分配，其中1,108个为`+1`，67个为`-1`；53个蛋白至少含1个`-1`分配。",
  "- “参考文献表”恰好有1,175行，并可按`ID + Pathway + Score`与每个非零分配一对一匹配。",
  "",
  "## 展示方案",
  "",
  "- 主图保留固定UMAP坐标，并把每个有分配的蛋白画成等权饼图；扇区面积表示通路成员关系，不表示表达量、强度或评分大小。",
  "- 实色扇区表示`+1`，同色浅色扇区表示`-1`；全零蛋白显示为灰点。",
  "- 由于饼图中的浅色负向扇区在密集区域较小，推荐使用“UMAP饼图 + 正负分配计数条形图”的组合图作为主展示；同时保留纯UMAP饼图。",
  "- Chromatin Interaction和Other support是功能类别，不应在文字中误称为经典DNA修复通路。",
  "",
  "## 输出",
  "",
  "- 推荐组合图：`reanalysis/results/figures/kla_ddr_pathway_pie_umap_33groups/kla_ddr_pathway_umap_and_signed_summary_33groups.{png,pdf,svg}`",
  "- 纯UMAP饼图：`reanalysis/results/figures/kla_ddr_pathway_pie_umap_33groups/kla_ddr_pathway_signed_pie_umap_33groups.{png,pdf,svg}`",
  "- 固定坐标绘图数据：`reanalysis/results/tables/kla_ddr_pathway_pie_umap_33groups/pathway_umap_plot_data_fixed_coordinates.csv`",
  "- 带文献证据的1,175条分配：`reanalysis/results/tables/kla_ddr_pathway_pie_umap_33groups/pathway_assignments_with_evidence_1175.csv`"
)
writeLines(report_lines, report_path, useBytes = TRUE)

message("Fixed UMAP proteins: ", nrow(plot_data))
message("Assigned proteins: ", plot_data[TotalAssignmentCount > 0L, .N])
message("Unassigned proteins: ", plot_data[TotalAssignmentCount == 0L, .N])
message("Nonzero assignments: ", nrow(assignment_with_evidence))
message("Suppressing assignments: ", sum(plot_data$NegativeAssignmentCount))
message("Pie UMAP: ", pie_paths[[1L]])
message("Recommended combined figure: ", combined_paths[[1L]])
message("Audit report: ", report_path)
