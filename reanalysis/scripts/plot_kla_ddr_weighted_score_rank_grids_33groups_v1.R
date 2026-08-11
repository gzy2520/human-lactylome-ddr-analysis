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
    "reanalysis/scripts/plot_kla_ddr_weighted_score_rank_grids_33groups_v1.R",
    mustWork = TRUE
  )
}
project_root <- normalizePath(
  file.path(dirname(script_path), "..", ".."),
  mustWork = TRUE
)

analysis_name <- "kla_ddr_weighted_score_rank_grids_33groups_v1"
table_dir <- file.path(project_root, "reanalysis/results/tables", analysis_name)
figure_dir <- file.path(project_root, "reanalysis/results/figures", analysis_name)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

score_workbook_path <- file.path(
  project_root,
  "data/identifier/260810乳酸化DDR基因评分表.xlsx"
)
membership_path <- file.path(
  project_root,
  "reanalysis/results/tables/four_class_venn/kla_ddr_four_class_venn/membership.csv"
)

assert <- function(condition, message) {
  if (!isTRUE(condition)) {
    stop(message, call. = FALSE)
  }
}

required_inputs <- c(score_workbook_path, membership_path)
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
excluded_non_pathway_columns <- c(
  "Chromatin Interaction",
  "Others (Transcription, RNA processing and proteostasis)"
)

score_raw <- as.data.table(
  read_excel(score_workbook_path, sheet = "评分表")
)
required_score_columns <- c(
  "ID",
  "BaseAccession",
  "GeneSymbol",
  "ProteinName",
  names(weights),
  excluded_non_pathway_columns,
  "score"
)
assert(
  all(required_score_columns %in% names(score_raw)),
  "The revised DDR score sheet is missing required columns."
)

blank_accession_rows <- score_raw[
  is.na(BaseAccession) | !nzchar(trimws(BaseAccession))
]
scores <- score_raw[
  !is.na(BaseAccession) & nzchar(trimws(BaseAccession))
]
scores[, BaseAccession := trimws(BaseAccession)]

assert(
  nrow(score_raw) == 514L &&
    nrow(blank_accession_rows) == 7L &&
    nrow(scores) == 507L &&
    uniqueN(scores$BaseAccession) == 507L,
  "The score workbook does not contain the expected 507 proteins plus seven blank rows."
)
assert(
  all(as.matrix(scores[, names(weights), with = FALSE]) %in% c(-1, 0, 1)),
  "The seven canonical pathway columns contain values outside -1/0/+1."
)

recalculated_score <- as.numeric(
  as.matrix(scores[, names(weights), with = FALSE]) %*% weights
)
assert(
  !anyNA(scores$score) &&
    all(is.finite(scores$score)) &&
    identical(as.numeric(scores$score), recalculated_score),
  paste0(
    "The cached Excel score values do not match ",
    "BER*1 + NER*2 + MMR*3 + FA*4 + HR*5 + AEJ*6 + NHEJ*7."
  )
)

scores[, SignedScore := as.numeric(score)]
scores[, AbsoluteScore := abs(SignedScore)]
scores[, ComponentwiseAbsoluteScore_NotPlotted := as.numeric(
  as.matrix(abs(scores[, names(weights), with = FALSE])) %*% weights
)]

membership <- fread(membership_path)
membership_columns <- c(
  "In_normal_tissue",
  "In_normal_cells",
  "In_cancer_tissue",
  "In_cancer_cells"
)
assert(
  nrow(membership) == 507L &&
    uniqueN(membership$BaseAccession) == 507L &&
    all(membership_columns %in% names(membership)) &&
    setequal(membership$BaseAccession, scores$BaseAccession),
  "The four-category membership table does not match the fixed 507-protein score set."
)

category_info <- data.table(
  Category = c(
    "normal_tissue",
    "normal_cells",
    "cancer_tissue",
    "cancer_cells",
    "all_507"
  ),
  MembershipColumn = c(membership_columns, "In_all_507"),
  CategoryOrder = 1:5,
  CategoryLabelEn = c(
    "Normal/non-tumor tissues",
    "Normal/non-tumor cells",
    "Cancer tissues",
    "Cancer cells",
    "All Kla-DDR proteins"
  ),
  CategoryLabelZh = c(
    "正常/非肿瘤组织",
    "正常/非肿瘤细胞",
    "癌症组织",
    "癌症细胞",
    "全部Kla∩DDR蛋白"
  )
)

membership[, In_all_507 := TRUE]
membership_long <- melt(
  membership[
    ,
    c("BaseAccession", category_info$MembershipColumn),
    with = FALSE
  ],
  id.vars = "BaseAccession",
  variable.name = "MembershipColumn",
  value.name = "DetectedInCategory"
)
membership_long <- merge(
  membership_long,
  category_info,
  by = "MembershipColumn",
  all.x = TRUE,
  sort = FALSE
)
membership_long <- membership_long[DetectedInCategory == TRUE]

protein_columns <- c(
  "BaseAccession",
  "GeneSymbol",
  "ProteinName",
  names(weights),
  excluded_non_pathway_columns,
  "SignedScore",
  "AbsoluteScore",
  "ComponentwiseAbsoluteScore_NotPlotted"
)
category_scores <- merge(
  membership_long,
  scores[, ..protein_columns],
  by = "BaseAccession",
  all.x = TRUE,
  sort = FALSE
)

observed_category_counts <- category_scores[
  ,
  .(ProteinCount = .N),
  by = .(
    Category,
    CategoryOrder,
    CategoryLabelEn,
    CategoryLabelZh
  )
][order(CategoryOrder)]
assert(
  identical(
    observed_category_counts$ProteinCount,
    c(183L, 471L, 178L, 383L, 507L)
  ),
  "The five protein-set counts are not 183/471/178/383/507."
)

score_long <- melt(
  category_scores,
  id.vars = setdiff(
    names(category_scores),
    c("SignedScore", "AbsoluteScore")
  ),
  measure.vars = c("SignedScore", "AbsoluteScore"),
  variable.name = "ScoreType",
  value.name = "PlottedScore"
)
score_long[, ScoreTypeOrder := fifelse(ScoreType == "SignedScore", 1L, 2L)]
score_long[, ScoreTypeLabelEn := fifelse(
  ScoreType == "SignedScore",
  "Original signed score",
  "Absolute score |signed score|"
)]
score_long[, ScoreTypeLabelZh := fifelse(
  ScoreType == "SignedScore",
  "原始有向得分",
  "得分绝对值 |原始得分|"
)]

ranked_plot_data <- copy(score_long)
setorder(
  ranked_plot_data,
  ScoreTypeOrder,
  CategoryOrder,
  PlottedScore,
  BaseAccession
)
ranked_plot_data[
  ,
  OrderIndex := seq_len(.N),
  by = .(ScoreType, Category)
]
ranked_plot_data[, OrderingMethod := "score_ascending"]

alphabetical_plot_data <- copy(score_long)
alphabetical_plot_data[, HasValidGeneSymbol := (
  !is.na(GeneSymbol) &
    nzchar(trimws(GeneSymbol)) &
    !grepl(
      "^\\[?no gene symbol",
      trimws(GeneSymbol),
      ignore.case = TRUE
    )
)]
alphabetical_plot_data[, AlphabeticalKey := fifelse(
  HasValidGeneSymbol,
  toupper(trimws(GeneSymbol)),
  BaseAccession
)]
setorder(
  alphabetical_plot_data,
  ScoreTypeOrder,
  CategoryOrder,
  AlphabeticalKey,
  BaseAccession
)
alphabetical_plot_data[
  ,
  OrderIndex := seq_len(.N),
  by = .(ScoreType, Category)
]
alphabetical_plot_data[, OrderingMethod := "gene_symbol_alphabetical"]

score_summary <- category_scores[
  ,
  .(
    ProteinCount = .N,
    Minimum = min(SignedScore),
    FirstQuartile = as.numeric(quantile(SignedScore, 0.25)),
    Median = median(SignedScore),
    ThirdQuartile = as.numeric(quantile(SignedScore, 0.75)),
    Maximum = max(SignedScore),
    NegativeCount = sum(SignedScore < 0),
    ZeroCount = sum(SignedScore == 0),
    PositiveCount = sum(SignedScore > 0)
  ),
  by = .(
    Category,
    CategoryOrder,
    CategoryLabelEn,
    CategoryLabelZh
  )
][order(CategoryOrder)]

weight_table <- data.table(
  Pathway = names(weights),
  Coefficient = as.numeric(weights),
  IncludedInScore = TRUE
)
weight_table <- rbind(
  weight_table,
  data.table(
    Pathway = excluded_non_pathway_columns,
    Coefficient = NA_real_,
    IncludedInScore = FALSE
  ),
  fill = TRUE
)

protein_score_audit <- scores[
  ,
  c(
    "ID",
    protein_columns
  ),
  with = FALSE
]
setorder(protein_score_audit, BaseAccession)
setorder(ranked_plot_data, ScoreTypeOrder, CategoryOrder, OrderIndex)
setorder(alphabetical_plot_data, ScoreTypeOrder, CategoryOrder, OrderIndex)

fwrite(
  protein_score_audit,
  file.path(table_dir, "protein_weighted_score_audit_507.csv")
)
fwrite(
  ranked_plot_data,
  file.path(table_dir, "weighted_score_plot_data_ranked_2x5.csv")
)
fwrite(
  alphabetical_plot_data,
  file.path(table_dir, "weighted_score_plot_data_alphabetical_2x5.csv")
)
fwrite(
  score_summary,
  file.path(table_dir, "weighted_score_summary_by_category.csv")
)
fwrite(
  weight_table,
  file.path(table_dir, "weighted_score_coefficients.csv")
)
fwrite(
  observed_category_counts,
  file.path(table_dir, "weighted_score_category_counts.csv")
)

input_audit <- data.table(
  Input = c("Revised DDR score workbook", "Fixed four-category membership"),
  Path = c(
    file.path("data/identifier", basename(score_workbook_path)),
    file.path(
      "reanalysis/results/tables/four_class_venn/kla_ddr_four_class_venn",
      basename(membership_path)
    )
  ),
  MD5 = unname(tools::md5sum(required_inputs)),
  ActiveProteinCount = c(nrow(scores), nrow(membership)),
  IdentifierKey = "isoform-stripped UniProt BaseAccession"
)
fwrite(input_audit, file.path(table_dir, "input_file_audit.csv"))

palette <- c(
  SignedScore = "#276FBF",
  AbsoluteScore = "#D97706"
)

theme_score_panel <- function(base_family) {
  theme_minimal(base_family = base_family, base_size = 9) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(
        colour = "#E3E7EB",
        linewidth = 0.32
      ),
      axis.line.x = element_line(colour = "#4B5563", linewidth = 0.35),
      axis.ticks.x = element_line(colour = "#4B5563", linewidth = 0.3),
      axis.text = element_text(colour = "#374151"),
      axis.title = element_text(colour = "#1F2937", size = 9),
      plot.title = element_text(
        colour = "#111827",
        face = "bold",
        size = 9.5,
        hjust = 0
      ),
      plot.subtitle = element_text(colour = "#6B7280", size = 8),
      plot.margin = margin(5, 7, 5, 7)
    )
}

score_breaks <- function(n) {
  unique(as.integer(round(c(1, (n + 1) / 2, n))))
}

make_panel <- function(
    data,
    category_index,
    score_type,
    language,
    ordering_method,
    show_x,
    show_y,
    show_category_title) {
  category_key <- category_info$Category[[category_index]]
  panel <- data[
    Category == category_key & ScoreType == score_type
  ]
  n <- nrow(panel)
  is_signed <- identical(score_type, "SignedScore")
  is_zh <- identical(language, "zh")
  base_family <- if (is_zh) "PingFang SC" else "Helvetica"
  category_title <- if (is_zh) {
    category_info$CategoryLabelZh[[category_index]]
  } else {
    category_info$CategoryLabelEn[[category_index]]
  }
  x_title <- if (ordering_method == "score_ascending") {
    if (is_zh) "蛋白排名（得分升序）" else "Protein rank (ascending score)"
  } else {
    if (is_zh) {
      "蛋白序号（GeneSymbol字母顺序）"
    } else {
      "Protein index (GeneSymbol alphabetical)"
    }
  }
  y_title <- if (is_signed) {
    if (is_zh) "原始有向得分" else "Original signed score"
  } else {
    if (is_zh) "得分绝对值" else "Absolute score"
  }

  plot <- ggplot(panel, aes(x = OrderIndex, y = PlottedScore)) +
    geom_line(
      colour = unname(palette[[score_type]]),
      linewidth = 0.48,
      lineend = "round"
    ) +
    geom_point(
      colour = unname(palette[[score_type]]),
      size = 0.62,
      alpha = 0.78
    )
  if (is_signed) {
    plot <- plot + geom_hline(
      yintercept = 0,
      colour = "#4B5563",
      linewidth = 0.38,
      linetype = "22"
    )
  }

  plot +
    scale_x_continuous(
      breaks = score_breaks(n),
      limits = c(1, n),
      expand = expansion(mult = c(0.012, 0.012))
    ) +
    scale_y_continuous(
      breaks = if (is_signed) seq(-10, 25, 5) else seq(0, 25, 5),
      limits = if (is_signed) c(-13, 26) else c(0, 26),
      expand = expansion(mult = c(0, 0))
    ) +
    labs(
      title = if (show_category_title) category_title else NULL,
      subtitle = if (show_category_title) paste0("n = ", n) else NULL,
      x = if (show_x) x_title else NULL,
      y = if (show_y) y_title else NULL
    ) +
    theme_score_panel(base_family) +
    theme(
      axis.text.x = if (show_x) {
        element_text(size = 7.2)
      } else {
        element_blank()
      },
      axis.ticks.x = if (show_x) {
        element_line(colour = "#4B5563", linewidth = 0.3)
      } else {
        element_blank()
      },
      axis.text.y = if (show_y) {
        element_text(size = 7.5)
      } else {
        element_blank()
      },
      axis.ticks.y = if (show_y) {
        element_line(colour = "#4B5563", linewidth = 0.3)
      } else {
        element_blank()
      },
      axis.line.y = if (show_y) {
        element_line(colour = "#4B5563", linewidth = 0.35)
      } else {
        element_blank()
      }
    )
}

make_grid <- function(data, language, ordering_method) {
  is_zh <- identical(language, "zh")
  base_family <- if (is_zh) "PingFang SC" else "Helvetica"
  top_row <- lapply(seq_len(nrow(category_info)), function(i) {
    make_panel(
      data = data,
      category_index = i,
      score_type = "SignedScore",
      language = language,
      ordering_method = ordering_method,
      show_x = FALSE,
      show_y = i == 1L,
      show_category_title = TRUE
    )
  })
  bottom_row <- lapply(seq_len(nrow(category_info)), function(i) {
    make_panel(
      data = data,
      category_index = i,
      score_type = "AbsoluteScore",
      language = language,
      ordering_method = ordering_method,
      show_x = TRUE,
      show_y = i == 1L,
      show_category_title = FALSE
    )
  })

  title <- if (is_zh) {
    "Kla∩DDR蛋白的通路难度加权得分"
  } else {
    "Difficulty-weighted pathway scores of Kla-DDR proteins"
  }
  subtitle <- if (ordering_method == "score_ascending") {
    if (is_zh) {
      paste0(
        "仅使用7条DDR通路；每个面板内独立按得分升序排列；",
        "相同得分按BaseAccession排序"
      )
    } else {
      paste0(
        "Seven DDR pathways only; proteins independently ranked by score ",
        "within each panel; ties ordered by BaseAccession"
      )
    }
  } else {
    if (is_zh) {
      paste0(
        "仅使用7条DDR通路；按GeneSymbol字母顺序展示；",
        "BaseAccession仍为分析键"
      )
    } else {
      paste0(
        "Seven DDR pathways only; display order follows GeneSymbol; ",
        "BaseAccession remains the analysis key"
      )
    }
  }
  caption <- if (is_zh) {
    paste0(
      "权重：BER=1，NER=2，MMR=3，FA=4，HR=5，AEJ=6，NHEJ=7；",
      "第二行为原始加权得分的绝对值。"
    )
  } else {
    paste0(
      "Weights: BER=1, NER=2, MMR=3, FA=4, HR=5, AEJ=6, NHEJ=7; ",
      "the second row is the absolute value of the signed weighted score."
    )
  }

  wrap_plots(c(top_row, bottom_row), ncol = 5, byrow = TRUE) +
    plot_annotation(
      title = title,
      subtitle = subtitle,
      caption = caption,
      theme = theme(
        text = element_text(family = base_family, colour = "#111827"),
        plot.title = element_text(face = "bold", size = 16),
        plot.subtitle = element_text(size = 10.5, colour = "#4B5563"),
        plot.caption = element_text(size = 8.5, colour = "#6B7280", hjust = 0),
        plot.margin = margin(10, 12, 8, 12)
      )
    )
}

save_grid <- function(plot, stem) {
  png_path <- file.path(figure_dir, paste0(stem, ".png"))
  pdf_path <- file.path(figure_dir, paste0(stem, ".pdf"))
  svg_path <- file.path(figure_dir, paste0(stem, ".svg"))

  ggsave(
    png_path,
    plot,
    width = 20,
    height = 8.4,
    units = "in",
    dpi = 450,
    device = ragg::agg_png,
    background = "white"
  )
  ggsave(
    pdf_path,
    plot,
    width = 20,
    height = 8.4,
    units = "in",
    device = cairo_pdf,
    bg = "white"
  )
  ggsave(
    svg_path,
    plot,
    width = 20,
    height = 8.4,
    units = "in",
    device = grDevices::svg,
    bg = "white"
  )
  invisible(c(png_path, pdf_path, svg_path))
}

ranked_grid_en <- make_grid(
  ranked_plot_data,
  language = "en",
  ordering_method = "score_ascending"
)
ranked_grid_zh <- make_grid(
  ranked_plot_data,
  language = "zh",
  ordering_method = "score_ascending"
)
alphabetical_grid_en <- make_grid(
  alphabetical_plot_data,
  language = "en",
  ordering_method = "gene_symbol_alphabetical"
)
alphabetical_grid_zh <- make_grid(
  alphabetical_plot_data,
  language = "zh",
  ordering_method = "gene_symbol_alphabetical"
)

figure_paths <- c(
  save_grid(
    ranked_grid_en,
    "kla_ddr_weighted_score_ranked_2x5_en"
  ),
  save_grid(
    ranked_grid_zh,
    "kla_ddr_weighted_score_ranked_2x5_zh"
  ),
  save_grid(
    alphabetical_grid_en,
    "kla_ddr_weighted_score_alphabetical_2x5_en"
  ),
  save_grid(
    alphabetical_grid_zh,
    "kla_ddr_weighted_score_alphabetical_2x5_zh"
  )
)

figure_manifest <- data.table(
  Ordering = rep(
    c("score_ascending", "gene_symbol_alphabetical"),
    each = 6L
  ),
  Language = rep(rep(c("en", "zh"), each = 3L), times = 2L),
  Format = rep(c("png", "pdf", "svg"), times = 4L),
  File = basename(figure_paths),
  WidthInches = 20,
  HeightInches = 8.4,
  PNGDPI = fifelse(grepl("\\.png$", figure_paths), 450L, NA_integer_)
)
fwrite(
  figure_manifest,
  file.path(table_dir, "weighted_score_figure_manifest.csv")
)

capture.output(
  sessionInfo(),
  file = file.path(table_dir, "session_info.txt")
)

cat(
  "Created ranked and alphabetical 2x5 weighted-score grids for 507 proteins.\n",
  "Primary ranked figure:\n",
  file.path(
    figure_dir,
    "kla_ddr_weighted_score_ranked_2x5_zh.png"
  ),
  "\n",
  sep = ""
)
