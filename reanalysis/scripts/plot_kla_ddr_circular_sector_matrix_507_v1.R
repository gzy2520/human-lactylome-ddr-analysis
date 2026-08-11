#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(readxl)
})

full_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", full_args, value = TRUE)
script_path <- if (length(file_arg) > 0L) {
  normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE)
} else {
  normalizePath(
    "reanalysis/scripts/plot_kla_ddr_circular_sector_matrix_507_v1.R",
    mustWork = TRUE
  )
}
project_root <- normalizePath(
  file.path(dirname(script_path), "..", ".."),
  mustWork = TRUE
)

analysis_name <- "kla_ddr_circular_sector_matrix_507_v1"
table_dir <- file.path(project_root, "reanalysis/results/tables", analysis_name)
figure_dir <- file.path(project_root, "reanalysis/results/figures", analysis_name)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

score_workbook_path <- file.path(
  project_root,
  "data/identifier/260810乳酸化DDR基因评分表.xlsx"
)
color_key_path <- file.path(
  project_root,
  paste0(
    "reanalysis/results/tables/",
    "kla_ddr_pathway_pie_umap_33groups_v4_bp_semantic/",
    "pathway_color_key_v4.csv"
  )
)

assert <- function(condition, message) {
  if (!isTRUE(condition)) {
    stop(message, call. = FALSE)
  }
}

required_inputs <- c(score_workbook_path, color_key_path)
assert(
  all(file.exists(required_inputs)),
  paste(
    "Missing required input(s):",
    paste(required_inputs[!file.exists(required_inputs)], collapse = "; ")
  )
)

weights <- c(
  BER = 1,
  NER = 2,
  MMR = 3,
  FA = 4,
  HR = 5,
  AEJ = 6,
  NHEJ = 7
)
pathway_order <- names(weights)

score_raw <- as.data.table(
  read_excel(score_workbook_path, sheet = "评分表")
)
required_columns <- c(
  "BaseAccession",
  "GeneSymbol",
  "ProteinName",
  pathway_order,
  "score"
)
assert(
  all(required_columns %in% names(score_raw)),
  "The revised score workbook is missing required columns."
)

scores <- score_raw[
  !is.na(BaseAccession) & nzchar(trimws(BaseAccession))
]
scores[, BaseAccession := trimws(BaseAccession)]
assert(
  nrow(scores) == 507L &&
    uniqueN(scores$BaseAccession) == 507L,
  "Expected exactly 507 unique BaseAccession identifiers."
)
assert(
  all(as.matrix(scores[, pathway_order, with = FALSE]) %in% c(-1, 0, 1)),
  "The seven pathway columns contain values outside -1/0/+1."
)

recalculated_score <- as.numeric(
  as.matrix(scores[, pathway_order, with = FALSE]) %*% weights
)
assert(
  !anyNA(scores$score) &&
    max(abs(as.numeric(scores$score) - recalculated_score)) < 1e-12,
  "The workbook score does not match the seven-pathway weighted formula."
)
scores[, SignedScore := recalculated_score]
setorder(scores, SignedScore, BaseAccession)
scores[, ProteinRank := seq_len(.N)]

color_key <- fread(color_key_path)
pathway_info <- data.table(
  Pathway = pathway_order,
  Coefficient = as.numeric(weights),
  SectorOrder = seq_along(pathway_order)
)
pathway_info <- merge(
  pathway_info,
  color_key[, .(Pathway = DisplayLabel, Color, ColorName)],
  by = "Pathway",
  all.x = TRUE,
  sort = FALSE
)
setorder(pathway_info, SectorOrder)
assert(
  nrow(pathway_info) == 7L &&
    !anyNA(pathway_info$Color) &&
    identical(pathway_info$Pathway, pathway_order),
  "The seven pathway colors could not be recovered."
)

matrix_long <- melt(
  scores[
    ,
    c(
      "BaseAccession",
      "GeneSymbol",
      "ProteinName",
      "SignedScore",
      "ProteinRank",
      pathway_order
    ),
    with = FALSE
  ],
  id.vars = c(
    "BaseAccession",
    "GeneSymbol",
    "ProteinName",
    "SignedScore",
    "ProteinRank"
  ),
  measure.vars = pathway_order,
  variable.name = "Pathway",
  value.name = "State"
)
matrix_long[, Pathway := as.character(Pathway)]
matrix_long <- merge(
  matrix_long,
  pathway_info,
  by = "Pathway",
  all.x = TRUE,
  sort = FALSE
)
setorder(matrix_long, SectorOrder, ProteinRank)

inner_hole <- 80
sector_half_width <- 0.43
matrix_long[, SectorCenter := SectorOrder]
matrix_long[, XMin := SectorCenter - sector_half_width]
matrix_long[, XMax := SectorCenter + sector_half_width]
matrix_long[, YMin := inner_hole + ProteinRank - 1]
matrix_long[, YMax := inner_hole + ProteinRank]
matrix_long[, StateLabel := fifelse(
  State == 1,
  "Promoting (+1)",
  fifelse(State == -1, "Suppressing (-1)", "Unassigned (0)")
)]

linear_matrix <- copy(matrix_long)
linear_matrix[, LinearY := 8 - SectorOrder]
linear_matrix[, LinearXMin := ProteinRank - 0.5]
linear_matrix[, LinearXMax := ProteinRank + 0.5]
linear_matrix[, LinearYMin := LinearY - 0.42]
linear_matrix[, LinearYMax := LinearY + 0.42]

pathway_summary <- matrix_long[
  ,
  .(
    SuppressingCount = sum(State == -1),
    UnassignedCount = sum(State == 0),
    PromotingCount = sum(State == 1),
    TotalProteins = .N
  ),
  by = .(
    Pathway,
    SectorOrder,
    Coefficient,
    Color,
    ColorName
  )
][order(SectorOrder)]
pathway_summary[
  ,
  `:=`(
    SuppressingFraction = SuppressingCount / TotalProteins,
    UnassignedFraction = UnassignedCount / TotalProteins,
    PromotingFraction = PromotingCount / TotalProteins
  )
]

rank_bin_size <- 25L
matrix_long[, RankBin := ceiling(ProteinRank / rank_bin_size)]
rank_bin_density <- matrix_long[
  ,
  .(
    BinStartRank = min(ProteinRank),
    BinEndRank = max(ProteinRank),
    BinProteinCount = uniqueN(BaseAccession),
    SuppressingCount = sum(State == -1),
    UnassignedCount = sum(State == 0),
    PromotingCount = sum(State == 1),
    SuppressingFraction = mean(State == -1),
    UnassignedFraction = mean(State == 0),
    PromotingFraction = mean(State == 1),
    MinimumSignedScore = min(SignedScore),
    MaximumSignedScore = max(SignedScore)
  ),
  by = .(Pathway, SectorOrder, RankBin)
][order(SectorOrder, RankBin)]

fwrite(
  scores[
    ,
    c(
      "BaseAccession",
      "GeneSymbol",
      "ProteinName",
      "SignedScore",
      "ProteinRank",
      pathway_order
    ),
    with = FALSE
  ],
  file.path(table_dir, "protein_order_and_seven_pathway_matrix_507.csv")
)
fwrite(
  matrix_long,
  file.path(table_dir, "circular_sector_plot_data_3549.csv")
)
fwrite(
  pathway_summary,
  file.path(table_dir, "pathway_state_summary_507.csv")
)
fwrite(
  rank_bin_density,
  file.path(table_dir, "pathway_state_density_by_25_rank_bin.csv")
)
fwrite(
  pathway_info,
  file.path(table_dir, "pathway_sector_order_and_colors.csv")
)

input_audit <- data.table(
  Input = c("Revised DDR score workbook", "Existing pathway color key"),
  Path = c(
    file.path("data/identifier", basename(score_workbook_path)),
    paste0(
      "reanalysis/results/tables/",
      "kla_ddr_pathway_pie_umap_33groups_v4_bp_semantic/",
      basename(color_key_path)
    )
  ),
  MD5 = unname(tools::md5sum(required_inputs)),
  Role = c(
    "507 proteins, seven signed pathway states, weighted score",
    "stable pathway colors reused from previous figures"
  )
)
fwrite(input_audit, file.path(table_dir, "input_file_audit.csv"))

zero_fill <- "#F1F3F5"
suppressing_fill <- "#2F3437"
guide_color <- "#4B5563"
outer_rank <- inner_hole + nrow(scores)
label_radius <- outer_rank + 34
plot_radius <- label_radius + 13
sector_start <- -pi / 2 - pi / 7

pathway_labels <- copy(pathway_info)
pathway_labels[, Label := Pathway]

make_circle_plot <- function(language) {
  is_zh <- identical(language, "zh")
  base_family <- if (is_zh) "PingFang SC" else "Helvetica"

  ggplot() +
    geom_rect(
      data = matrix_long[State == 0],
      aes(xmin = XMin, xmax = XMax, ymin = YMin, ymax = YMax),
      fill = zero_fill,
      colour = NA
    ) +
    geom_rect(
      data = matrix_long[State == 1],
      aes(
        xmin = XMin,
        xmax = XMax,
        ymin = YMin,
        ymax = YMax,
        fill = Color
      ),
      colour = NA
    ) +
    geom_rect(
      data = matrix_long[State == -1],
      aes(
        xmin = XMin,
        xmax = XMax,
        ymin = YMin,
        ymax = YMax
      ),
      fill = suppressing_fill,
      colour = NA
    ) +
    geom_hline(
      yintercept = c(inner_hole, outer_rank),
      colour = "#9CA3AF",
      linewidth = 0.25
    ) +
    geom_text(
      data = pathway_labels,
      aes(x = SectorOrder, y = label_radius, label = Label, colour = Color),
      family = base_family,
      fontface = "bold",
      size = 4.1,
      lineheight = 0.92,
      show.legend = FALSE
    ) +
    scale_fill_identity() +
    scale_colour_identity() +
    scale_x_continuous(limits = c(0.5, 7.5), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, plot_radius), expand = c(0, 0)) +
    coord_polar(
      theta = "x",
      start = sector_start,
      direction = 1,
      clip = "off"
    ) +
    theme_void(base_family = base_family) +
    theme(
      plot.margin = margin(18, 28, 18, 28)
    )
}

make_manual_legend <- function(language) {
  is_zh <- identical(language, "zh")
  base_family <- if (is_zh) "PingFang SC" else "Helvetica"
  legend_labels <- if (is_zh) {
    c(
      "+1 促进：通路实色",
      "−1 抑制：深炭灰",
      "0 未分配：浅灰"
    )
  } else {
    c(
      "+1 promoting: solid pathway color",
      "−1 suppressing: dark charcoal",
      "0 unassigned: light gray"
    )
  }
  legend_data <- data.table(
    X = 1:3,
    Label = legend_labels,
    Fill = c("#3C5488", suppressing_fill, zero_fill),
    Border = c("#3C5488", suppressing_fill, zero_fill)
  )

  ggplot() +
    geom_rect(
      data = legend_data,
      aes(
        xmin = X - 0.31,
        xmax = X + 0.31,
        ymin = 0.38,
        ymax = 0.68,
        fill = Fill,
        colour = Border
      ),
      linewidth = 0.45
    ) +
    geom_text(
      data = legend_data,
      aes(x = X, y = 0.18, label = Label),
      family = base_family,
      size = 3.7,
      colour = "#374151"
    ) +
    scale_fill_identity() +
    scale_colour_identity() +
    coord_cartesian(xlim = c(0.45, 3.55), ylim = c(0.05, 0.78), clip = "off") +
    theme_void(base_family = base_family) +
    theme(plot.margin = margin(0, 10, 0, 10))
}

make_summary_plot <- function(language) {
  is_zh <- identical(language, "zh")
  base_family <- if (is_zh) "PingFang SC" else "Helvetica"
  positive <- pathway_summary[
    ,
    .(
      Pathway,
      SectorOrder,
      Y = 8 - SectorOrder,
      Color,
      Fraction = PromotingFraction,
      Count = PromotingCount,
      Direction = "promoting"
    )
  ]
  negative <- pathway_summary[
    ,
    .(
      Pathway,
      SectorOrder,
      Y = 8 - SectorOrder,
      Color,
      Fraction = -SuppressingFraction,
      Count = SuppressingCount,
      Direction = "suppressing"
    )
  ]
  labels <- pathway_summary[
    ,
    .(
      Pathway,
      SectorOrder,
      Y = 8 - SectorOrder,
      ZeroLabel = if (is_zh) {
        paste0("0：", UnassignedCount)
      } else {
        paste0("0: ", UnassignedCount)
      }
    )
  ]

  positive[, CountLabel := sprintf("%d (%.1f%%)", Count, 100 * Fraction)]
  negative[, CountLabel := sprintf("%d (%.1f%%)", Count, 100 * abs(Fraction))]
  ggplot() +
    geom_rect(
      data = positive,
      aes(
        xmin = 0,
        xmax = Fraction,
        ymin = Y - 0.29,
        ymax = Y + 0.29,
        fill = Color
      )
    ) +
    geom_rect(
      data = negative,
      aes(
        xmin = Fraction,
        xmax = 0,
        ymin = Y - 0.29,
        ymax = Y + 0.29,
        colour = Color
      ),
      fill = suppressing_fill,
      linewidth = 0.6
    ) +
    geom_vline(
      xintercept = 0,
      colour = guide_color,
      linewidth = 0.4
    ) +
    geom_text(
      data = positive,
      aes(
        x = Fraction + 0.008,
        y = Y,
        label = CountLabel
      ),
      family = base_family,
      hjust = 0,
      size = 3.3,
      colour = "#374151"
    ) +
    geom_text(
      data = negative,
      aes(
        x = pmin(Fraction - 0.004, -0.004),
        y = Y,
        label = CountLabel
      ),
      family = base_family,
      hjust = 1,
      size = 3.3,
      colour = "#374151"
    ) +
    geom_text(
      data = labels,
      aes(
        x = 0.54,
        y = Y,
        label = ZeroLabel
      ),
      family = base_family,
      hjust = 0,
      size = 3.15,
      colour = "#6B7280"
    ) +
    scale_fill_identity() +
    scale_colour_identity() +
    scale_x_continuous(
      limits = c(-0.16, 0.66),
      breaks = c(-0.1, 0, 0.1, 0.2, 0.3, 0.4),
      labels = function(x) paste0(abs(round(100 * x)), "%"),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      limits = c(0.5, 7.5),
      breaks = 7:1,
      labels = pathway_order,
      expand = c(0, 0)
    ) +
    labs(
      title = if (is_zh) {
        "各通路状态比例"
      } else {
        "State proportions by pathway"
      },
      subtitle = if (is_zh) {
        "左：−1抑制；右：+1促进；数字为n（占507的比例）"
      } else {
        "Left: −1 suppressing; right: +1 promoting; labels are n (% of 507)"
      },
      x = if (is_zh) "蛋白比例" else "Fraction of proteins",
      y = NULL
    ) +
    theme_minimal(base_family = base_family, base_size = 10) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(colour = "#E5E7EB", linewidth = 0.35),
      axis.text.y = element_text(
        face = "bold",
        colour = "#374151",
        size = 10.5
      ),
      axis.text.x = element_text(colour = "#4B5563"),
      axis.title.x = element_text(colour = "#374151"),
      plot.title = element_text(face = "bold", size = 13, colour = "#111827"),
      plot.subtitle = element_text(size = 9.2, colour = "#6B7280"),
      plot.margin = margin(20, 45, 18, 14)
    )
}

make_figure <- function(language, include_summary = FALSE) {
  is_zh <- identical(language, "zh")
  base_family <- if (is_zh) "PingFang SC" else "Helvetica"
  circle <- make_circle_plot(language)
  legend <- make_manual_legend(language)
  circle_with_legend <- circle / legend +
    plot_layout(heights = c(1, 0.115))

  title <- if (is_zh) {
    "507个Kla∩DDR蛋白的7通路环形扇区矩阵"
  } else {
    "Seven-pathway circular sector matrix of 507 Kla-DDR proteins"
  }
  subtitle <- if (is_zh) {
    paste0(
      "每条同心弧代表1个蛋白；得分仅用于决定由内向外的升序排列；",
      "扇区颜色仅表示原始−1/0/+1通路状态"
    )
  } else {
    paste0(
      "Each concentric arc is one protein; score is used only for the ascending ",
      "inner-to-outer order; sector colors encode the original −1/0/+1 pathway states"
    )
  }
  caption <- if (is_zh) {
    paste0(
      "得分不参与颜色、弧线宽度或状态比例计算。局部颜色集中可描述为“局部比例较高”，",
      "不能仅凭本图称为统计学富集。"
    )
  } else {
    paste0(
      "Score does not determine color, arc width, or state proportion. Local color ",
      "concentration is descriptive and is not a statistical enrichment test."
    )
  }

  body <- if (include_summary) {
    circle_with_legend | make_summary_plot(language)
  } else {
    circle_with_legend
  }
  widths <- if (include_summary) c(1.65, 1) else 1

  body +
    plot_layout(widths = widths) +
    plot_annotation(
      title = title,
      subtitle = subtitle,
      caption = caption,
      theme = theme(
        text = element_text(family = base_family, colour = "#111827"),
        plot.title = element_text(face = "bold", size = 17),
        plot.subtitle = element_text(size = 10.5, colour = "#4B5563"),
        plot.caption = element_text(size = 8.8, colour = "#6B7280", hjust = 0),
        plot.margin = margin(10, 16, 8, 16)
      )
    )
}

make_linear_plot <- function(language) {
  is_zh <- identical(language, "zh")
  base_family <- if (is_zh) "PingFang SC" else "Helvetica"

  ggplot() +
    geom_rect(
      data = linear_matrix[State == 0],
      aes(
        xmin = LinearXMin,
        xmax = LinearXMax,
        ymin = LinearYMin,
        ymax = LinearYMax
      ),
      fill = zero_fill,
      colour = NA
    ) +
    geom_rect(
      data = linear_matrix[State == 1],
      aes(
        xmin = LinearXMin,
        xmax = LinearXMax,
        ymin = LinearYMin,
        ymax = LinearYMax,
        fill = Color
      ),
      colour = NA
    ) +
    geom_rect(
      data = linear_matrix[State == -1],
      aes(
        xmin = LinearXMin,
        xmax = LinearXMax,
        ymin = LinearYMin,
        ymax = LinearYMax
      ),
      fill = suppressing_fill,
      colour = NA
    ) +
    scale_fill_identity() +
    scale_colour_identity() +
    scale_x_continuous(
      limits = c(0.5, 507.5),
      breaks = c(1, 127, 254, 380, 507),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      limits = c(0.5, 7.5),
      breaks = 7:1,
      labels = pathway_order,
      expand = c(0, 0)
    ) +
    labs(
      x = if (is_zh) {
        "蛋白排名（按得分升序；得分仅用于排序）"
      } else {
        "Protein rank (ascending score; score used for ordering only)"
      },
      y = NULL
    ) +
    theme_minimal(base_family = base_family, base_size = 10) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(colour = "#D1D5DB", linewidth = 0.34),
      axis.text.y = element_text(
        face = "bold",
        colour = "#374151",
        size = 10.5
      ),
      axis.text.x = element_text(colour = "#4B5563", size = 9),
      axis.title.x = element_text(colour = "#374151", size = 10),
      plot.margin = margin(10, 12, 8, 10)
    )
}

make_linear_figure <- function(language, include_summary = FALSE) {
  is_zh <- identical(language, "zh")
  base_family <- if (is_zh) "PingFang SC" else "Helvetica"
  matrix_with_legend <- make_linear_plot(language) /
    make_manual_legend(language) +
    plot_layout(heights = c(1, 0.24))

  title <- if (is_zh) {
    "507个Kla∩DDR蛋白的7通路线性状态矩阵"
  } else {
    "Linear seven-pathway state matrix of 507 Kla-DDR proteins"
  }
  subtitle <- if (is_zh) {
    paste0(
      "环形扇区图的展开版本；每列为同一个蛋白；",
      "得分仅决定从左到右的升序排列"
    )
  } else {
    paste0(
      "Unrolled version of the circular sector matrix; each column is one protein; ",
      "score determines only the left-to-right order"
    )
  }
  caption <- if (is_zh) {
    paste0(
      "颜色仅来自原始−1/0/+1通路状态。连续色块可用于观察局部比例，",
      "不能仅凭本图称为统计学富集。"
    )
  } else {
    paste0(
      "Colors encode only the original −1/0/+1 pathway states. ",
      "Contiguous blocks show local proportions, not statistical enrichment."
    )
  }

  body <- if (include_summary) {
    matrix_with_legend | make_summary_plot(language)
  } else {
    matrix_with_legend
  }
  widths <- if (include_summary) c(2.3, 1) else 1

  body +
    plot_layout(widths = widths) +
    plot_annotation(
      title = title,
      subtitle = subtitle,
      caption = caption,
      theme = theme(
        text = element_text(family = base_family, colour = "#111827"),
        plot.title = element_text(face = "bold", size = 17),
        plot.subtitle = element_text(size = 10.5, colour = "#4B5563"),
        plot.caption = element_text(size = 8.8, colour = "#6B7280", hjust = 0),
        plot.margin = margin(10, 16, 8, 16)
      )
    )
}

save_figure <- function(plot, stem, width, height) {
  output_paths <- file.path(
    figure_dir,
    paste0(stem, c(".png", ".pdf", ".svg"))
  )
  ggsave(
    output_paths[[1L]],
    plot,
    width = width,
    height = height,
    units = "in",
    dpi = 600,
    device = ragg::agg_png,
    background = "white"
  )
  ggsave(
    output_paths[[2L]],
    plot,
    width = width,
    height = height,
    units = "in",
    device = cairo_pdf,
    bg = "white"
  )
  ggsave(
    output_paths[[3L]],
    plot,
    width = width,
    height = height,
    units = "in",
    device = grDevices::svg,
    bg = "white"
  )
  output_paths
}

circle_en <- make_figure("en", include_summary = FALSE)
circle_zh <- make_figure("zh", include_summary = FALSE)
summary_en <- make_figure("en", include_summary = TRUE)
summary_zh <- make_figure("zh", include_summary = TRUE)
linear_en <- make_linear_figure("en", include_summary = FALSE)
linear_zh <- make_linear_figure("zh", include_summary = FALSE)
linear_summary_en <- make_linear_figure("en", include_summary = TRUE)
linear_summary_zh <- make_linear_figure("zh", include_summary = TRUE)

figure_paths <- c(
  save_figure(
    circle_en,
    "kla_ddr_circular_sector_matrix_507_en",
    width = 10.2,
    height = 11.3
  ),
  save_figure(
    circle_zh,
    "kla_ddr_circular_sector_matrix_507_zh",
    width = 10.2,
    height = 11.3
  ),
  save_figure(
    summary_en,
    "kla_ddr_circular_sector_matrix_with_summary_507_en",
    width = 16.8,
    height = 10.4
  ),
  save_figure(
    summary_zh,
    "kla_ddr_circular_sector_matrix_with_summary_507_zh",
    width = 16.8,
    height = 10.4
  ),
  save_figure(
    linear_en,
    "kla_ddr_linear_pathway_matrix_507_en",
    width = 15.8,
    height = 6.4
  ),
  save_figure(
    linear_zh,
    "kla_ddr_linear_pathway_matrix_507_zh",
    width = 15.8,
    height = 6.4
  ),
  save_figure(
    linear_summary_en,
    "kla_ddr_linear_pathway_matrix_with_summary_507_en",
    width = 20,
    height = 7.2
  ),
  save_figure(
    linear_summary_zh,
    "kla_ddr_linear_pathway_matrix_with_summary_507_zh",
    width = 20,
    height = 7.2
  )
)

figure_manifest <- data.table(
  FigureType = rep(
    c(
      "circle_only",
      "circle_with_pathway_summary",
      "linear_only",
      "linear_with_pathway_summary"
    ),
    each = 6L
  ),
  Language = rep(rep(c("en", "zh"), each = 3L), times = 4L),
  Format = rep(c("png", "pdf", "svg"), times = 8L),
  File = basename(figure_paths),
  PNGDPI = fifelse(grepl("\\.png$", figure_paths), 600L, NA_integer_)
)
fwrite(
  figure_manifest,
  file.path(table_dir, "circular_sector_figure_manifest.csv")
)

capture.output(
  sessionInfo(),
  file = file.path(table_dir, "session_info.txt")
)

cat(
  "Created circular and linear seven-pathway matrices for 507 proteins.\n",
  "Primary Chinese circular figure:\n",
  file.path(
    figure_dir,
    "kla_ddr_circular_sector_matrix_with_summary_507_zh.png"
  ),
  "\nPrimary Chinese linear figure:\n",
  file.path(
    figure_dir,
    "kla_ddr_linear_pathway_matrix_with_summary_507_zh.png"
  ),
  "\n",
  sep = ""
)
