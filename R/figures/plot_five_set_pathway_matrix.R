#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(readxl)
})

full_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", full_args, value = TRUE)
script_path <- if (length(file_arg) > 0L) {
  normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE)
} else {
  normalizePath(
    "R/figures/plot_five_set_pathway_matrix.R",
    mustWork = TRUE
  )
}
project_root <- normalizePath(
  file.path(dirname(script_path), "..", ".."),
  mustWork = TRUE
)

analysis_name <- Sys.getenv(
  "KLA_SCORE_ANALYSIS_NAME",
  unset = "five_set_pathway_matrix"
)
table_dir <- file.path(project_root, "results/tables", analysis_name)
figure_dir <- file.path(project_root, "results/figures", analysis_name)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

score_workbook_relative <- Sys.getenv(
  "KLA_SCORE_WORKBOOK_PATH",
  unset = "data/identifier/260810乳酸化DDR基因评分表.xlsx"
)
score_workbook_path <- if (grepl("^/", score_workbook_relative)) {
  score_workbook_relative
} else {
  file.path(project_root, score_workbook_relative)
}
membership_path <- file.path(
  project_root,
  "results/tables/four_class_venn/",
  "kla_ddr_four_class_venn/membership.csv"
)
color_key_path <- file.path(
  project_root,
  "results", "tables", "protein_function_inputs",
  "pathway_colors.csv"
)
display_override_relative <- Sys.getenv(
  "KLA_SCORE_PATHWAY_DISPLAY",
  unset = ""
)
display_override_path <- if (!nzchar(display_override_relative)) {
  ""
} else if (grepl("^/", display_override_relative)) {
  display_override_relative
} else {
  file.path(project_root, display_override_relative)
}
display_input_path <- if (nzchar(display_override_path)) {
  display_override_path
} else {
  color_key_path
}

assert <- function(condition, message) {
  if (!isTRUE(condition)) {
    stop(message, call. = FALSE)
  }
}

required_inputs <- c(
  score_workbook_path,
  membership_path,
  display_input_path
)
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
score_columns <- names(weights)

set_info <- data.table(
  Set = c(
    "normal_tissue",
    "cancer_tissue",
    "cancer_cells",
    "normal_cells",
    "all_kla_ddr"
  ),
  MembershipColumn = c(
    "In_normal_tissue",
    "In_cancer_tissue",
    "In_cancer_cells",
    "In_normal_cells",
    "In_all_kla_ddr"
  ),
  SetOrder = 1:5,
  SetLabelEn = c(
    "non-tumor tissues",
    "tumor tissues",
    "cancer cell lines",
    "normal cell lines",
    "All Kla-DDR proteins"
  ),
  SetLabelZh = c(
    "非肿瘤组织",
    "肿瘤组织",
    "癌细胞系",
    "正常细胞系",
    "全部Kla∩DDR蛋白"
  )
)

score_raw <- as.data.table(
  read_excel(score_workbook_path, sheet = "评分表")
)
required_score_columns <- c(
  "BaseAccession",
  "GeneSymbol",
  "ProteinName",
  score_columns,
  "score"
)
assert(
  all(required_score_columns %in% names(score_raw)),
  "The revised score workbook is missing required columns."
)
scores <- score_raw[
  !is.na(BaseAccession) & nzchar(trimws(BaseAccession))
]
scores[, BaseAccession := trimws(BaseAccession)]
assert(
  uniqueN(scores$BaseAccession) == nrow(scores) &&
    all(as.matrix(scores[, score_columns, with = FALSE]) %in% c(-1, 0, 1)),
  "The score workbook does not contain a valid unique seven-pathway matrix."
)
recalculated_score <- as.numeric(
  as.matrix(scores[, score_columns, with = FALSE]) %*% weights
)
assert(
  !anyNA(scores$score) &&
    max(abs(as.numeric(scores$score) - recalculated_score)) < 1e-12,
  "The workbook score does not match the seven-pathway weighted formula."
)
scores[, SignedScore := recalculated_score]
scored_input_count <- nrow(scores)

membership <- fread(membership_path)
membership_columns <- set_info$MembershipColumn[set_info$Set != "all_kla_ddr"]
assert(
  nrow(membership) > 0L &&
    uniqueN(membership$BaseAccession) == nrow(membership) &&
    all(membership_columns %in% names(membership)) &&
    all(membership$BaseAccession %in% scores$BaseAccession),
  "One or more proteins in the current membership lack pathway scores."
)
score_scope_audit <- scores[
  ,
  .(
    BaseAccession,
    GeneSymbol,
    ProteinName,
    WorkbookScore = as.numeric(score),
    RecalculatedSignedScore = SignedScore,
    InCurrent399 = BaseAccession %in% membership$BaseAccession
  )
]
assert(
  scored_input_count == 507L &&
    sum(score_scope_audit$InCurrent399) == 399L,
  "Expected 507 scored proteins with complete coverage of the current 399."
)
scores <- scores[BaseAccession %in% membership$BaseAccession]
membership[, In_all_kla_ddr := TRUE]
membership_long <- melt(
  membership[
    ,
    c("BaseAccession", set_info$MembershipColumn),
    with = FALSE
  ],
  id.vars = "BaseAccession",
  variable.name = "MembershipColumn",
  value.name = "Included"
)
membership_long <- merge(
  membership_long,
  set_info,
  by = "MembershipColumn",
  all.x = TRUE,
  sort = FALSE
)
membership_long <- membership_long[Included == TRUE]

if (nzchar(display_override_path)) {
  pathway_info <- fread(display_override_path)
  assert(
    all(c("Pathway", "PathwayOrder", "Color", "ColorName") %in%
      names(pathway_info)),
    "The pathway display override is missing required columns."
  )
  pathway_info[, Coefficient := as.numeric(weights[Pathway])]
  setorder(pathway_info, PathwayOrder)
} else {
  color_key <- fread(color_key_path)
  pathway_info <- data.table(
    Pathway = score_columns,
    Coefficient = as.numeric(weights),
    PathwayOrder = seq_along(score_columns)
  )
  pathway_info <- merge(
    pathway_info,
    color_key[, .(Pathway = DisplayLabel, Color, ColorName)],
    by = "Pathway",
    all.x = TRUE,
    sort = FALSE
  )
  setorder(pathway_info, PathwayOrder)
}
pathway_order <- pathway_info$Pathway
assert(
  nrow(pathway_info) == 7L &&
    !anyNA(pathway_info$Color) &&
    setequal(pathway_info$Pathway, score_columns) &&
    !anyNA(pathway_info$Coefficient),
  "The seven fixed pathway colors could not be recovered."
)

protein_by_set <- merge(
  membership_long,
  scores[
    ,
    c(
      "BaseAccession",
      "GeneSymbol",
      "ProteinName",
      "SignedScore",
      pathway_order
    ),
    with = FALSE
  ],
  by = "BaseAccession",
  all.x = TRUE,
  sort = FALSE
)
setorder(protein_by_set, SetOrder, SignedScore, BaseAccession)
protein_by_set[
  ,
  ProteinRank := seq_len(.N),
  by = Set
]

observed_set_counts <- protein_by_set[
  ,
  .(
    ProteinCount = .N,
    MinimumSignedScore = min(SignedScore),
    MaximumSignedScore = max(SignedScore)
  ),
  by = .(
    Set,
    SetOrder,
    SetLabelEn,
    SetLabelZh
  )
][order(SetOrder)]
expected_counts <- c(
  vapply(
    membership_columns,
    function(column) sum(membership[[column]] %in% c(TRUE, "TRUE", 1, "1")),
    integer(1),
    USE.NAMES = FALSE
  ),
  nrow(membership)
)
assert(
  identical(observed_set_counts$ProteinCount, expected_counts),
  "The five protein-set counts do not match the current membership table."
)

matrix_long <- melt(
  protein_by_set,
  id.vars = c(
    "BaseAccession",
    "GeneSymbol",
    "ProteinName",
    "SignedScore",
    "ProteinRank",
    "Set",
    "SetOrder",
    "SetLabelEn",
    "SetLabelZh"
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
matrix_long[, PathwayY := 8 - PathwayOrder]
matrix_long[, XMin := ProteinRank - 0.5]
matrix_long[, XMax := ProteinRank + 0.5]
matrix_long[, YMin := PathwayY - 0.42]
matrix_long[, YMax := PathwayY + 0.42]
setorder(matrix_long, SetOrder, PathwayOrder, ProteinRank)

pathway_summary <- matrix_long[
  ,
  .(
    SuppressingCount = sum(State == -1),
    UnassignedCount = sum(State == 0),
    PromotingCount = sum(State == 1),
    ProteinCount = .N
  ),
  by = .(
    Set,
    SetOrder,
    SetLabelEn,
    SetLabelZh,
    Pathway,
    PathwayOrder,
    Color
  )
][order(SetOrder, PathwayOrder)]
pathway_summary[
  ,
  `:=`(
    SuppressingFraction = SuppressingCount / ProteinCount,
    UnassignedFraction = UnassignedCount / ProteinCount,
    PromotingFraction = PromotingCount / ProteinCount
  )
]

fwrite(
  protein_by_set[
    ,
    c(
      "Set",
      "SetOrder",
      "SetLabelEn",
      "SetLabelZh",
      "BaseAccession",
      "GeneSymbol",
      "ProteinName",
      "SignedScore",
      "ProteinRank",
      pathway_order
    ),
    with = FALSE
  ],
  file.path(table_dir, "protein_order_and_seven_pathway_matrix_5sets.csv")
)
fwrite(
  matrix_long,
  file.path(table_dir, "linear_matrix_plot_data.csv")
)
fwrite(
  pathway_summary,
  file.path(table_dir, "pathway_state_summary_5sets_35rows.csv")
)
fwrite(
  observed_set_counts,
  file.path(table_dir, "protein_set_counts.csv")
)
fwrite(
  score_scope_audit,
  file.path(table_dir, "score_workbook_scope_audit_507_to_399.csv")
)
fwrite(
  pathway_info,
  file.path(table_dir, "pathway_order_and_colors.csv")
)

input_audit <- data.table(
  Input = c(
    "Revised DDR score workbook",
    "Current four-category membership",
    "Pathway display configuration"
  ),
  Path = c(
    file.path("data/identifier", basename(score_workbook_path)),
    paste0(
      "results/tables/four_class_venn/",
      "kla_ddr_four_class_venn/",
      basename(membership_path)
    ),
    sub(paste0("^", project_root, "/?"), "", display_input_path)
  ),
  SHA256 = vapply(
    required_inputs,
    function(path) digest::digest(
      file = path,
      algo = "sha256",
      serialize = FALSE
    ),
    character(1)
  ),
  Role = c(
    paste0(
      scored_input_count,
      " scored proteins; filtered to ",
      nrow(scores),
      " current proteins; seven signed pathway states and score ordering key"
    ),
    "non-tumor tissues, tumor tissues, cancer cell lines, normal cell lines membership",
    "stable pathway colors reused from previous figures"
  )
)
fwrite(input_audit, file.path(table_dir, "input_file_audit.csv"))

zero_fill <- "#F1F3F5"
suppressing_fill <- "#2F3437"
guide_color <- "#4B5563"

matrix_breaks <- function(n) {
  unique(as.integer(round(c(1, (n + 1) / 2, n))))
}

make_matrix_panel <- function(set_key, language, show_x_title = FALSE) {
  is_zh <- identical(language, "zh")
  base_family <- "Arial Unicode MS"
  panel <- matrix_long[Set == set_key]
  n <- unique(panel[, uniqueN(BaseAccession)])
  set_row <- set_info[Set == set_key]
  ggplot() +
    geom_rect(
      data = panel[State == 0],
      aes(xmin = XMin, xmax = XMax, ymin = YMin, ymax = YMax),
      fill = zero_fill,
      colour = NA
    ) +
    geom_rect(
      data = panel[State == 1],
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
      data = panel[State == -1],
      aes(xmin = XMin, xmax = XMax, ymin = YMin, ymax = YMax),
      fill = suppressing_fill,
      colour = NA
    ) +
    scale_fill_identity() +
    scale_x_continuous(
      limits = c(0.5, n + 0.5),
      breaks = matrix_breaks(n),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      limits = c(0.5, 7.5),
      breaks = 7:1,
      labels = pathway_order,
      expand = c(0, 0)
    ) +
    labs(
      title = NULL,
      x = if (show_x_title) {
        if (is_zh) {
          "蛋白排名（各面板内按得分升序；得分仅用于排序）"
        } else {
          "Protein rank (ascending score within each panel; ordering only)"
        }
      } else {
        NULL
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
        size = 9.8
      ),
      axis.text.x = element_text(colour = "#4B5563", size = 12.5),
      axis.title.x = element_text(colour = "#374151", size = 10.5),
      plot.margin = margin(7, 12, 7, 10)
    )
}

make_matrix_figure <- function(set_key, language) {
  make_matrix_panel(set_key, language, show_x_title = TRUE)
}

make_summary_panel <- function(set_key, language, show_y = TRUE) {
  is_zh <- identical(language, "zh")
  base_family <- "Arial Unicode MS"
  panel <- copy(pathway_summary[Set == set_key])
  panel[, Y := 8 - PathwayOrder]
  panel[, PositiveLabel := sprintf(
    "%d\n(%.1f%%)",
    PromotingCount,
    100 * PromotingFraction
  )]
  panel[, NegativeLabel := sprintf(
    "%d\n(%.1f%%)",
    SuppressingCount,
    100 * SuppressingFraction
  )]
  panel[, ZeroLabel := paste0("0: ", UnassignedCount)]
  ggplot() +
    geom_rect(
      data = panel,
      aes(
        xmin = 0,
        xmax = PromotingFraction,
        ymin = Y - 0.29,
        ymax = Y + 0.29,
        fill = Color
      ),
      colour = NA
    ) +
    geom_rect(
      data = panel,
      aes(
        xmin = -SuppressingFraction,
        xmax = 0,
        ymin = Y - 0.29,
        ymax = Y + 0.29
      ),
      fill = suppressing_fill,
      colour = NA
    ) +
    geom_vline(
      xintercept = 0,
      colour = guide_color,
      linewidth = 0.4
    ) +
    geom_text(
      data = panel,
      aes(
        x = PromotingFraction + 0.008,
        y = Y,
        label = PositiveLabel
      ),
      family = base_family,
      hjust = 0,
      lineheight = 0.92,
      size = 2.75,
      colour = "#374151"
    ) +
    geom_text(
      data = panel,
      aes(
        x = pmin(-SuppressingFraction - 0.005, -0.005),
        y = Y,
        label = NegativeLabel
      ),
      family = base_family,
      hjust = 1,
      lineheight = 0.92,
      size = 2.75,
      colour = "#374151"
    ) +
    geom_text(
      data = panel,
      aes(x = 0.48, y = Y, label = ZeroLabel),
      family = base_family,
      hjust = 0,
      size = 2.8,
      colour = "#6B7280"
    ) +
    scale_fill_identity() +
    scale_x_continuous(
      limits = c(-0.13, 0.62),
      breaks = c(-0.1, 0, 0.1, 0.2, 0.3, 0.4),
      labels = function(x) paste0(abs(round(100 * x)), "%"),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      limits = c(0.5, 7.5),
      breaks = 7:1,
      labels = if (show_y) pathway_order else rep("", 7),
      expand = c(0, 0)
    ) +
    labs(
      title = NULL,
      x = if (is_zh) "蛋白比例" else "Protein fraction",
      y = NULL
    ) +
    theme_minimal(base_family = base_family, base_size = 9) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(colour = "#E5E7EB", linewidth = 0.32),
      axis.text.y = element_text(
        face = "bold",
        colour = "#374151",
        size = 8.8
      ),
      axis.text.x = element_text(colour = "#4B5563", size = 11),
      axis.title.x = element_text(colour = "#374151", size = 9.5),
      plot.margin = margin(8, 13, 8, 8)
    )
}

make_summary_figure <- function(set_key, language) {
  make_summary_panel(set_key, language, show_y = TRUE)
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

figure_paths <- character()
manifest_rows <- list()
manifest_index <- 0L
summary_height_inches <- 5.2

for (set_key in set_info$Set) {
  for (language in c("en", "zh")) {
    matrix_stem <- paste0(
      "kla_ddr_linear_pathway_matrix_",
      set_key,
      "_",
      language
    )
    summary_stem <- paste0(
      "kla_ddr_pathway_state_summary_",
      set_key,
      "_",
      language
    )
    matrix_paths <- save_figure(
      make_matrix_figure(set_key, language),
      matrix_stem,
      width = 15.8,
      height = 6.7
    )
    summary_paths <- save_figure(
      make_summary_figure(set_key, language),
      summary_stem,
      width = 8.8,
      height = summary_height_inches
    )
    figure_paths <- c(figure_paths, matrix_paths, summary_paths)

    manifest_index <- manifest_index + 1L
    manifest_rows[[manifest_index]] <- data.table(
      FigureType = "linear_matrix",
      Set = set_key,
      Language = language,
      Format = c("png", "pdf", "svg"),
      File = basename(matrix_paths),
      WidthInches = 15.8,
      HeightInches = 6.7,
      PNGDPI = c(600L, NA_integer_, NA_integer_)
    )
    manifest_index <- manifest_index + 1L
    manifest_rows[[manifest_index]] <- data.table(
      FigureType = "pathway_summary",
      Set = set_key,
      Language = language,
      Format = c("png", "pdf", "svg"),
      File = basename(summary_paths),
      WidthInches = 8.8,
      HeightInches = summary_height_inches,
      PNGDPI = c(600L, NA_integer_, NA_integer_)
    )
  }
}

figure_manifest <- rbindlist(manifest_rows)
fwrite(
  figure_manifest,
  file.path(table_dir, "figure_manifest.csv")
)

capture.output(
  sessionInfo(),
  file = file.path(table_dir, "session_info.txt")
)

cat(
  "Created five separate linear matrices and five separate summaries per language.\n",
  "Example Chinese normal-tissue matrix:\n",
  file.path(
    figure_dir,
    "kla_ddr_linear_pathway_matrix_normal_tissue_zh.png"
  ),
  "\nExample Chinese normal-tissue summary:\n",
  file.path(
    figure_dir,
    "kla_ddr_pathway_state_summary_normal_tissue_zh.png"
  ),
  "\n",
  sep = ""
)
