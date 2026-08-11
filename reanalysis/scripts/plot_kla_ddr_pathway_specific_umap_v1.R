#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
script_path <- if (length(file_arg) > 0L) {
  normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE)
} else {
  normalizePath(
    "reanalysis/scripts/plot_kla_ddr_pathway_specific_umap_v1.R",
    mustWork = TRUE
  )
}
project_root <- normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = TRUE)

v3_table_dir <- file.path(
  project_root,
  "reanalysis/results/tables/kla_ddr_pathway_pie_umap_33groups_v3_all_go"
)
v4_table_dir <- file.path(
  project_root,
  "reanalysis/results/tables/kla_ddr_pathway_pie_umap_33groups_v4_bp_semantic"
)
assignment_path <- file.path(v3_table_dir, "pathway_assignment_long_v3.csv")
metadata_path <- file.path(v3_table_dir, "pathway_umap_plot_data_v3_all_go.csv")
coordinate_path <- file.path(v4_table_dir, "umap_raw_coordinates_v4_bp_semantic.csv")
color_path <- file.path(v4_table_dir, "pathway_color_key_v4.csv")

analysis_name <- "kla_ddr_pathway_specific_umap_33groups_v1"
table_dir <- file.path(project_root, "reanalysis/results/tables", analysis_name)
figure_dir <- file.path(project_root, "reanalysis/results/figures", analysis_name)
individual_figure_dir <- file.path(figure_dir, "individual_pathways")
report_path <- file.path(
  project_root,
  "reanalysis/reports/UMAP_PATHWAY_SPECIFIC_33GROUP_V1.md"
)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(individual_figure_dir, recursive = TRUE, showWarnings = FALSE)

required_inputs <- c(
  assignment_path,
  metadata_path,
  coordinate_path,
  color_path
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

assignment_long <- fread(assignment_path)
protein_metadata <- fread(metadata_path)[
  ,
  .(BaseAccession, GeneSymbol, ProteinName)
]
raw_coordinates <- fread(coordinate_path)
pathway_info <- fread(color_path)[order(PathwayOrder)]

assert(
  nrow(raw_coordinates) == 507L &&
    uniqueN(raw_coordinates$BaseAccession) == 507L,
  "Expected 507 unique proteins in the raw BP semantic UMAP coordinates."
)
assert(
  nrow(pathway_info) == 9L &&
    uniqueN(pathway_info$DisplayLabel) == 9L,
  "Expected nine unique pathway/function display labels."
)
assert(
  all(assignment_long$Score %in% c(-1L, 1L)) &&
    uniqueN(
      assignment_long,
      by = c("BaseAccession", "DisplayLabel")
    ) == nrow(assignment_long),
  "Signed pathway assignments must be unique nonzero protein-pathway pairs."
)
assert(
  setequal(raw_coordinates$BaseAccession, protein_metadata$BaseAccession),
  "Raw UMAP and protein metadata sets do not match."
)
assert(
  all(assignment_long$BaseAccession %in% raw_coordinates$BaseAccession) &&
    all(assignment_long$DisplayLabel %in% pathway_info$DisplayLabel),
  "Assignment proteins or pathway labels do not match the fixed inputs."
)

pathway_levels <- pathway_info$DisplayLabel
full_grid <- CJ(
  BaseAccession = raw_coordinates$BaseAccession,
  DisplayLabel = pathway_levels,
  unique = TRUE
)
plot_data <- merge(
  full_grid,
  assignment_long[, .(BaseAccession, DisplayLabel, Score)],
  by = c("BaseAccession", "DisplayLabel"),
  all.x = TRUE,
  sort = FALSE
)
plot_data[is.na(Score), Score := 0L]
plot_data <- merge(
  plot_data,
  raw_coordinates,
  by = "BaseAccession",
  all.x = TRUE,
  sort = FALSE
)
plot_data <- merge(
  plot_data,
  protein_metadata,
  by = "BaseAccession",
  all.x = TRUE,
  sort = FALSE
)
plot_data <- merge(
  plot_data,
  pathway_info[, .(DisplayLabel, Color, PathwayOrder)],
  by = "DisplayLabel",
  all.x = TRUE,
  sort = FALSE
)

plot_data[
  ,
  Status := factor(
    fifelse(
      Score == 1L,
      "Promoting (+1)",
      fifelse(Score == -1L, "Suppressing (-1)", "Not assigned (0)")
    ),
    levels = c("Not assigned (0)", "Promoting (+1)", "Suppressing (-1)")
  )
]
plot_data[
  ,
  PointColor := fifelse(Score == 0L, "#C9C9C9", Color)
]
plot_data[
  ,
  DrawOrder := fifelse(Score == 0L, 1L, fifelse(Score == 1L, 2L, 3L))
]
setorder(plot_data, PathwayOrder, DrawOrder, BaseAccession)

assert(
  nrow(plot_data) == 507L * 9L &&
    plot_data[, uniqueN(BaseAccession), by = DisplayLabel][, all(V1 == 507L)] &&
    all(is.finite(plot_data$UMAP_1)) &&
    all(is.finite(plot_data$UMAP_2)),
  "The complete pathway-specific plotting grid is invalid."
)

pathway_summary <- plot_data[
  ,
  .(
    PromotingCount = sum(Score == 1L),
    SuppressingCount = sum(Score == -1L),
    UnassignedCount = sum(Score == 0L),
    TotalProteins = .N,
    Color = unique(Color),
    PathwayOrder = unique(PathwayOrder)
  ),
  by = DisplayLabel
][order(PathwayOrder)]
pathway_summary[
  ,
  FacetLabel := paste0(
    DisplayLabel,
    "\n+1: ",
    PromotingCount,
    "  |  -1: ",
    SuppressingCount
  )
]
plot_data <- merge(
  plot_data,
  pathway_summary[, .(DisplayLabel, FacetLabel)],
  by = "DisplayLabel",
  all.x = TRUE,
  sort = FALSE
)
plot_data[
  ,
  FacetLabel := factor(
    FacetLabel,
    levels = pathway_summary$FacetLabel
  )
]
setorder(plot_data, PathwayOrder, DrawOrder, BaseAccession)

shape_values <- c(
  "Not assigned (0)" = 16,
  "Promoting (+1)" = 16,
  "Suppressing (-1)" = 1
)
size_values_overview <- c(
  "Not assigned (0)" = 0.62,
  "Promoting (+1)" = 1.15,
  "Suppressing (-1)" = 1.55
)
alpha_values <- c(
  "Not assigned (0)" = 0.42,
  "Promoting (+1)" = 0.96,
  "Suppressing (-1)" = 1
)

overview_plot <- ggplot(
  plot_data,
  aes(
    x = UMAP_1,
    y = UMAP_2,
    shape = Status,
    size = Status,
    alpha = Status,
    color = I(PointColor)
  )
) +
  geom_point(stroke = 0.72) +
  facet_wrap(~FacetLabel, ncol = 3, scales = "fixed") +
  scale_shape_manual(
    values = shape_values,
    name = "Pathway assignment",
    drop = FALSE
  ) +
  scale_size_manual(
    values = size_values_overview,
    guide = "none"
  ) +
  scale_alpha_manual(
    values = alpha_values,
    guide = "none"
  ) +
  coord_equal() +
  labs(
    title = "Pathway-specific distributions on the raw BP semantic UMAP",
    subtitle = paste0(
      "507 proteins in the Kla-DDR intersection from the fixed 33-group scope | ",
      "same raw UMAP coordinates in every panel"
    ),
    x = "UMAP 1",
    y = "UMAP 2",
    caption = paste0(
      "Solid pathway-colored points: promoting (+1); hollow points with ",
      "pathway-colored outlines: suppressing (-1); light-gray points: not assigned.\n",
      "All panels use the exact same raw BP semantic UMAP coordinates; ",
      "no collision displacement was applied."
    )
  ) +
  guides(
    shape = guide_legend(
      order = 1,
      override.aes = list(
        color = c("#BDBDBD", "#333333", "#333333"),
        size = c(2.2, 2.6, 3.0),
        alpha = 1
      )
    )
  ) +
  theme_classic(base_size = 10, base_family = "sans") +
  theme(
    panel.border = element_rect(color = "#B8B8B8", fill = NA, linewidth = 0.35),
    strip.background = element_rect(fill = "#F3F3F3", color = "#B8B8B8", linewidth = 0.35),
    strip.text = element_text(face = "bold", size = 9, lineheight = 1.08),
    axis.text = element_text(size = 7.4, color = "#333333"),
    axis.title = element_text(size = 9),
    axis.line = element_blank(),
    axis.ticks = element_line(linewidth = 0.3, color = "#555555"),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 9.5, color = "#4D4D4D"),
    plot.caption = element_text(size = 8, color = "#5A5A5A", hjust = 0),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_text(face = "bold", size = 9),
    legend.text = element_text(size = 8.5),
    plot.margin = margin(9, 10, 8, 9)
  )

save_plot_formats <- function(plot, output_directory, stem, width, height) {
  png_path <- file.path(output_directory, paste0(stem, ".png"))
  pdf_path <- file.path(output_directory, paste0(stem, ".pdf"))
  svg_path <- file.path(output_directory, paste0(stem, ".svg"))
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

overview_paths <- save_plot_formats(
  overview_plot,
  figure_dir,
  "kla_ddr_pathway_specific_umap_3x3_v1",
  width = 12.2,
  height = 11.2
)

individual_paths <- vector("list", nrow(pathway_summary))
for (pathway_index in seq_len(nrow(pathway_summary))) {
  pathway_label <- pathway_summary$DisplayLabel[[pathway_index]]
  pathway_data <- plot_data[DisplayLabel == pathway_label]
  pathway_color <- pathway_summary$Color[[pathway_index]]
  individual_plot <- ggplot(
    pathway_data,
    aes(
      x = UMAP_1,
      y = UMAP_2,
      shape = Status,
      size = Status,
      alpha = Status,
      color = I(PointColor)
    )
  ) +
    geom_point(stroke = 0.82) +
    scale_shape_manual(
      values = shape_values,
      name = "Pathway assignment",
      drop = FALSE
    ) +
    scale_size_manual(
      values = c(
        "Not assigned (0)" = 1.10,
        "Promoting (+1)" = 2.10,
        "Suppressing (-1)" = 2.60
      ),
      guide = "none"
    ) +
    scale_alpha_manual(values = alpha_values, guide = "none") +
    coord_equal() +
    labs(
      title = pathway_label,
      subtitle = paste0(
        "Promoting +1: ",
        pathway_summary$PromotingCount[[pathway_index]],
        " | Suppressing -1: ",
        pathway_summary$SuppressingCount[[pathway_index]],
        " | Not assigned: ",
        pathway_summary$UnassignedCount[[pathway_index]]
      ),
      x = "UMAP 1",
      y = "UMAP 2",
      caption = paste0(
        "Raw BP semantic UMAP coordinates; no collision displacement. ",
        "Pathway color: ",
        pathway_color,
        "."
      )
    ) +
    guides(
      shape = guide_legend(
        override.aes = list(
          color = c("#BDBDBD", pathway_color, pathway_color),
          size = c(2.2, 2.8, 3.2),
          alpha = 1
        )
      )
    ) +
    theme_classic(base_size = 10, base_family = "sans") +
    theme(
      axis.text = element_text(size = 8, color = "#333333"),
      axis.title = element_text(size = 9.5),
      axis.line = element_line(linewidth = 0.4, color = "#333333"),
      axis.ticks = element_line(linewidth = 0.35, color = "#333333"),
      plot.title = element_text(face = "bold", size = 14, color = pathway_color),
      plot.subtitle = element_text(size = 9, color = "#4D4D4D"),
      plot.caption = element_text(size = 7.7, hjust = 0, color = "#5A5A5A"),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.title = element_text(face = "bold", size = 9),
      legend.text = element_text(size = 8.5),
      plot.margin = margin(9, 9, 8, 9)
    )
  file_label <- gsub("[^A-Za-z0-9]+", "_", pathway_label)
  file_label <- gsub("^_|_$", "", file_label)
  stem <- sprintf("%02d_%s_raw_umap_v1", pathway_index, file_label)
  individual_paths[[pathway_index]] <- save_plot_formats(
    individual_plot,
    individual_figure_dir,
    stem,
    width = 7.2,
    height = 6.4
  )
}

plot_table_path <- file.path(table_dir, "pathway_specific_umap_plot_data.csv")
summary_table_path <- file.path(table_dir, "pathway_specific_umap_summary.csv")
input_audit_path <- file.path(table_dir, "input_file_audit.csv")
session_info_path <- file.path(table_dir, "session_info.txt")
fwrite(
  plot_data[
    ,
    .(
      BaseAccession,
      GeneSymbol,
      ProteinName,
      DisplayLabel,
      PathwayOrder,
      Score,
      Status,
      Color,
      UMAP_1,
      UMAP_2
    )
  ],
  plot_table_path
)
fwrite(pathway_summary, summary_table_path)
fwrite(
  data.table(
    InputRole = c(
      "signed pathway assignment long table",
      "protein display metadata",
      "raw BP semantic UMAP coordinates",
      "pathway color key"
    ),
    Path = relative_path(required_inputs),
    MD5 = unname(tools::md5sum(required_inputs))
  ),
  input_audit_path
)
writeLines(capture.output(sessionInfo()), session_info_path, useBytes = TRUE)

report_lines <- c(
  "# Kla-DDR通路分面UMAP（33组，V1）",
  "",
  "## 固定范围",
  "",
  "- 使用固定33组范围产生的507个Kla∩DDR蛋白。",
  "- 分析与合并键为去isoform的UniProt `BaseAccession`；GeneSymbol仅用于显示。",
  "- 所有面板复用V4 BP semantic UMAP的原始坐标，未执行防碰撞位移。",
  "- UMAP坐标由直接UniProt BP term及GO.db中可获得的BP祖先构建；通路评分不参与UMAP。",
  "",
  "## 图形编码",
  "",
  "- 每个面板只突出一个通路/功能，并保留全部507个蛋白作为共同背景。",
  "- `+1`：通路颜色实心点。",
  "- `-1`：白色空心点，边框使用该通路颜色。",
  "- `0`：浅灰色背景点，表示未分配至当前面板的通路/功能。",
  "- 九宫格和九张单独通路图使用相同原始UMAP坐标范围，因此可以直接比较空间分布。",
  "",
  "## 产物",
  "",
  paste0("- 九宫格主图：`", relative_path(overview_paths[[1L]]), "`"),
  paste0("- 九张单独通路图目录：`", relative_path(individual_figure_dir), "`"),
  paste0("- 完整4,563行绘图表：`", relative_path(plot_table_path), "`"),
  paste0("- 通路计数表：`", relative_path(summary_table_path), "`"),
  paste0("- 输入文件审计：`", relative_path(input_audit_path), "`")
)
writeLines(report_lines, report_path, useBytes = TRUE)

message("Overview PNG: ", overview_paths[[1L]])
message("Individual pathway directory: ", individual_figure_dir)
message("Plot data: ", plot_table_path)
message("Summary table: ", summary_table_path)
message("Report: ", report_path)
