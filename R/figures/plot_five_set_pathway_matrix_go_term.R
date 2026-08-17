#!/usr/bin/env Rscript

# Draw the 4+1 current 30-group linear pathway matrices and summaries from
# direct GO-term counts. Direct term counts determine protein order; the matrix
# itself encodes pathway presence/absence.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

full_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", full_args, value = TRUE)
script_path <- if (length(file_arg) > 0L) {
  normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE)
} else {
  normalizePath(
    "R/figures/plot_five_set_pathway_matrix_go_term.R",
    mustWork = TRUE
  )
}
project_root <- normalizePath(
  file.path(dirname(script_path), "..", ".."),
  mustWork = TRUE
)

analysis_name <- "five_set_pathway_matrix_go_term_30groups"
table_dir <- file.path(project_root, "results/tables", analysis_name)
figure_dir <- file.path(project_root, "results/figures", analysis_name)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

score_path <- file.path(
  project_root,
  "results/tables/go_term_pathway_scoring_30groups/",
  "protein_pathway_direct_term_count_matrix.csv"
)
membership_path <- file.path(
  project_root,
  "results/tables/four_class_venn/",
  "kla_ddr_four_class_venn/membership.csv"
)
display_path <- file.path(project_root, "config/seven_pathway_display.csv")
required_inputs <- c(score_path, membership_path, display_path)

assert <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

assert(
  all(file.exists(required_inputs)),
  paste(
    "Missing required input(s):",
    paste(required_inputs[!file.exists(required_inputs)], collapse = "; ")
  )
)

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

pathway_info <- fread(display_path)
setorder(pathway_info, PathwayOrder)
pathway_order <- pathway_info$Pathway
assert(
  identical(pathway_order, c("BER", "NER", "MMR", "FA", "HR", "NHEJ", "AEJ")) &&
    !anyNA(pathway_info$Color),
  "The seven-pathway display contract changed."
)

scores <- fread(score_path)
required_score_columns <- c(
  "BaseAccession",
  "GeneSymbolAudit",
  "ProteinNameAudit",
  pathway_order,
  "Others",
  "SevenPathwayTermScore",
  "DirectGOCount"
)
assert(
  all(required_score_columns %in% names(scores)) &&
    nrow(scores) == 399L &&
    uniqueN(scores$BaseAccession) == 399L,
  "Expected a complete 399-protein direct-GO pathway count matrix."
)
assert(
  all(as.matrix(scores[, ..pathway_order]) >= 0L) &&
    all(scores$Others >= 0L) &&
    identical(
      as.integer(scores$SevenPathwayTermScore),
      as.integer(rowSums(scores[, ..pathway_order]))
    ),
  "The direct-GO term counts or seven-pathway term score are invalid."
)

membership <- fread(membership_path)
membership[, In_all_kla_ddr := TRUE]
membership_columns <- set_info$MembershipColumn[
  set_info$Set != "all_kla_ddr"
]
assert(
  nrow(membership) == 399L &&
    uniqueN(membership$BaseAccession) == 399L &&
    all(membership_columns %in% names(membership)) &&
    setequal(membership$BaseAccession, scores$BaseAccession),
  "The pathway scores and current four-class membership are not identical."
)

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
membership_long <- membership_long[Included %in% c(TRUE, "TRUE", 1, "1")]

protein_by_set <- merge(
  membership_long,
  scores[
    ,
    c(
      "BaseAccession",
      "GeneSymbolAudit",
      "ProteinNameAudit",
      pathway_order,
      "Others",
      "SevenPathwayTermScore",
      "DirectGOCount"
    ),
    with = FALSE
  ],
  by = "BaseAccession",
  all.x = TRUE,
  sort = FALSE
)

# Term counts define score and order. BaseAccession is the only tie-breaker so
# the display does not manufacture pathway-specific blocks within equal scores.
setorderv(
  protein_by_set,
  c(
    "SetOrder",
    "SevenPathwayTermScore",
    "BaseAccession"
  )
)
protein_by_set[, ProteinRank := seq_len(.N), by = Set]

observed_set_counts <- protein_by_set[
  ,
  .(
    ProteinCount = .N,
    MinimumSevenPathwayTermScore = min(SevenPathwayTermScore),
    MaximumSevenPathwayTermScore = max(SevenPathwayTermScore),
    ProteinsWithoutSevenPathwayTerm = sum(SevenPathwayTermScore == 0L)
  ),
  by = .(
    Set,
    SetOrder,
    SetLabelEn,
    SetLabelZh
  )
][order(SetOrder)]
assert(
  identical(
    observed_set_counts$ProteinCount,
    c(183L, 178L, 381L, 292L, 399L)
  ),
  "Expected five current protein-set sizes 183/178/381/292/399."
)

matrix_long <- melt(
  protein_by_set,
  id.vars = c(
    "BaseAccession",
    "GeneSymbolAudit",
    "ProteinNameAudit",
    "Others",
    "SevenPathwayTermScore",
    "DirectGOCount",
    "ProteinRank",
    "Set",
    "SetOrder",
    "SetLabelEn",
    "SetLabelZh"
  ),
  measure.vars = pathway_order,
  variable.name = "Pathway",
  value.name = "DirectTermCount"
)
matrix_long[, Pathway := as.character(Pathway)]
matrix_long <- merge(
  matrix_long,
  pathway_info,
  by = "Pathway",
  all.x = TRUE,
  sort = FALSE
)
matrix_long[, `:=`(
  Present = DirectTermCount > 0L,
  PathwayY = 8L - PathwayOrder,
  XMin = ProteinRank - 0.5,
  XMax = ProteinRank + 0.5,
  YMin = 8L - PathwayOrder - 0.42,
  YMax = 8L - PathwayOrder + 0.42
)]
setorder(matrix_long, SetOrder, PathwayOrder, ProteinRank)

pathway_summary <- matrix_long[
  ,
  .(
    AssignedProteinCount = sum(Present),
    UnassignedProteinCount = sum(!Present),
    ProteinCount = .N,
    TotalDirectTermAssignments = sum(DirectTermCount),
    MeanDirectTermsAmongAssigned = ifelse(
      any(Present),
      as.numeric(mean(DirectTermCount[Present])),
      0.0
    ),
    MedianDirectTermsAmongAssigned = ifelse(
      any(Present),
      as.numeric(median(DirectTermCount[Present])),
      0.0
    )
  ),
  by = .(
    Set,
    SetOrder,
    SetLabelEn,
    SetLabelZh,
    Pathway,
    PathwayOrder,
    Color,
    ColorName
  )
][order(SetOrder, PathwayOrder)]
pathway_summary[, `:=`(
  AssignedFraction = AssignedProteinCount / ProteinCount,
  UnassignedFraction = UnassignedProteinCount / ProteinCount
)]
assert(nrow(pathway_summary) == 35L, "Expected a 5 x 7 pathway summary.")

fwrite(
  protein_by_set[
    ,
    c(
      "Set",
      "SetOrder",
      "SetLabelEn",
      "SetLabelZh",
      "BaseAccession",
      "GeneSymbolAudit",
      "ProteinNameAudit",
      "SevenPathwayTermScore",
      "Others",
      "DirectGOCount",
      "ProteinRank",
      pathway_order
    ),
    with = FALSE
  ],
  file.path(table_dir, "protein_order_and_direct_term_counts_5sets.csv")
)
fwrite(
  matrix_long,
  file.path(table_dir, "linear_matrix_plot_data.csv")
)
fwrite(
  pathway_summary,
  file.path(table_dir, "pathway_presence_summary_5sets_35rows.csv")
)
fwrite(
  observed_set_counts,
  file.path(table_dir, "protein_set_counts.csv")
)
fwrite(
  pathway_info,
  file.path(table_dir, "pathway_order_and_colors.csv")
)

input_audit <- data.table(
  Input = c(
    "GO-term-derived protein pathway count matrix",
    "Current four-category membership",
    "Seven-pathway display configuration"
  ),
  Path = c(
    sub(paste0("^", project_root, "/?"), "", score_path),
    sub(paste0("^", project_root, "/?"), "", membership_path),
    sub(paste0("^", project_root, "/?"), "", display_path)
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
    "direct GO-term counts, Others count, seven-pathway term-score ordering key",
    "non-tumor tissues, tumor tissues, cancer cell lines, normal cell lines membership",
    "fixed pathway order and publication colors"
  )
)
fwrite(input_audit, file.path(table_dir, "input_file_audit.csv"))

zero_fill <- "#EEF1F4"
guide_color <- "#D1D5DB"
base_family <- "Arial Unicode MS"

matrix_breaks <- function(n) {
  unique(as.integer(round(c(1, (n + 1) / 2, n))))
}

make_matrix_panel <- function(set_key, language) {
  is_zh <- identical(language, "zh")
  panel <- matrix_long[Set == set_key]
  n <- uniqueN(panel$BaseAccession)
  ggplot() +
    geom_rect(
      data = panel[Present == FALSE],
      aes(xmin = XMin, xmax = XMax, ymin = YMin, ymax = YMax),
      fill = zero_fill,
      colour = NA
    ) +
    geom_rect(
      data = panel[Present == TRUE],
      aes(
        xmin = XMin,
        xmax = XMax,
        ymin = YMin,
        ymax = YMax,
        fill = Color
      ),
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
      subtitle = NULL,
      caption = NULL,
      x = if (is_zh) {
        "蛋白排名（按七通路直接GO term总数升序）"
      } else {
        "Protein rank (ascending total direct GO terms across seven pathways)"
      },
      y = NULL
    ) +
    theme_minimal(base_family = base_family, base_size = 10) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(
        colour = guide_color,
        linewidth = 0.34
      ),
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

make_summary_panel <- function(set_key, language) {
  is_zh <- identical(language, "zh")
  panel <- copy(pathway_summary[Set == set_key])
  panel[, `:=`(
    Y = 8L - PathwayOrder,
    AssignedLabel = sprintf(
      "%d (%.1f%%)",
      AssignedProteinCount,
      100 * AssignedFraction
    )
  )]
  ggplot(panel) +
    geom_rect(
      aes(
        xmin = 0,
        xmax = 1,
        ymin = Y - 0.29,
        ymax = Y + 0.29
      ),
      fill = zero_fill,
      colour = NA
    ) +
    geom_rect(
      aes(
        xmin = 0,
        xmax = AssignedFraction,
        ymin = Y - 0.29,
        ymax = Y + 0.29,
        fill = Color
      ),
      colour = NA
    ) +
    geom_text(
      aes(
        x = pmin(AssignedFraction + 0.018, 0.86),
        y = Y,
        label = AssignedLabel
      ),
      family = base_family,
      hjust = 0,
      size = 3.25,
      colour = "#374151"
    ) +
    scale_fill_identity() +
    scale_x_continuous(
      limits = c(0, 1.02),
      breaks = seq(0, 1, by = 0.2),
      labels = function(x) paste0(round(100 * x), "%"),
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
      subtitle = NULL,
      caption = NULL,
      x = if (is_zh) {
        "至少含1个对应直接GO term的蛋白比例"
      } else {
        "Proteins with at least one corresponding direct GO term"
      },
      y = NULL
    ) +
    theme_minimal(base_family = base_family, base_size = 9) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(
        colour = "#E5E7EB",
        linewidth = 0.32
      ),
      axis.text.y = element_text(
        face = "bold",
        colour = "#374151",
        size = 9.2
      ),
      axis.text.x = element_text(colour = "#4B5563", size = 11),
      axis.title.x = element_text(colour = "#374151", size = 9.5),
      plot.margin = margin(8, 13, 8, 8)
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

manifest_rows <- list()
manifest_index <- 0L
summary_height_inches <- 5.2
for (set_key in set_info$Set) {
  for (language in c("en", "zh")) {
    matrix_stem <- paste0(
      "kla_ddr_go_term_linear_pathway_matrix_",
      set_key,
      "_",
      language
    )
    summary_stem <- paste0(
      "kla_ddr_go_term_pathway_summary_",
      set_key,
      "_",
      language
    )
    matrix_paths <- save_figure(
      make_matrix_panel(set_key, language),
      matrix_stem,
      width = 15.8,
      height = 6.7
    )
    summary_paths <- save_figure(
      make_summary_panel(set_key, language),
      summary_stem,
      width = 8.8,
      height = summary_height_inches
    )

    manifest_index <- manifest_index + 1L
    manifest_rows[[manifest_index]] <- data.table(
      FigureType = "go_term_linear_matrix",
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
      FigureType = "go_term_pathway_summary",
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

fwrite(
  rbindlist(manifest_rows),
  file.path(table_dir, "figure_manifest.csv")
)
capture.output(
  sessionInfo(),
  file = file.path(table_dir, "session_info.txt")
)

message(
  "Created five separate GO-term-derived linear matrices and summaries ",
  "per language in: ",
  figure_dir
)
