#!/usr/bin/env Rscript

# Render the 4 split DDR pathway matrices as standalone publication figures:
# Supplementary Figure S3a: DDR pathway matrix — tumor tissues (192 proteins)
# Supplementary Figure S3b: DDR pathway matrix — non-tumor tissues (183 proteins)
# Supplementary Figure S4a: DDR pathway matrix — cancer cell lines (381 proteins)
# Supplementary Figure S4b: DDR pathway matrix — normal cell lines (292 proteins)

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
  library(ggplot2)
  library(dplyr)
})

set.seed(25)

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1L) normalizePath(args[[1L]], mustWork = TRUE) else normalizePath(".", mustWork = TRUE)
stop_if <- function(condition, message) if (!isTRUE(condition)) stop(message, call. = FALSE)

publication_font <- "Arial Unicode MS"
pathway_order <- c("BER", "NER", "MMR", "FA", "HR", "AEJ", "NHEJ")
pathway_weights <- stats::setNames(seq_along(pathway_order), pathway_order)
pathway_display_order <- c("BER", "NER", "MMR", "FA", "HR", "NHEJ", "AEJ")
pathway_colours <- c(
  BER = "#54BED4", NER = "#F59E83", MMR = "#8C9ABD", FA = "#8ED5C4",
  HR = "#4C669D", AEJ = "#00A98F", NHEJ = "#E94F3D"
)
zero_fill <- "#F1F3F5"
suppressing_fill <- "#2F3437"
guide_colour <- "#4B5563"

input_dir <- normalizePath(Sys.getenv(
  "KLA_PUBLICATION_INPUT",
  unset = file.path(project_root, "data", "candidate", "escc_inclusion_20260903_pxd065830_tumor_reference", "publication_input")
), mustWork = TRUE)

s3_dir <- file.path(project_root, "results", "final_figures_and_tables", "Supplementary_Figure_S3")
s4_dir <- file.path(project_root, "results", "final_figures_and_tables", "Supplementary_Figure_S4")
dir.create(s3_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(s4_dir, recursive = TRUE, showWarnings = FALSE)

# Load 31-group S4/S5 pathway protein ranking table
ranking_path <- file.path(input_dir, "Supplementary_Table_S4_Pathway_Protein_Ranking.xlsx")
stop_if(file.exists(ranking_path), paste("Missing pathway ranking workbook:", ranking_path))

kla_ddr_membership <- fread(file.path(input_dir, "venn_kla_ddr.csv"))

matrix_specs <- list(
  cancer_tissue = list(
    key = "cancer_tissue", sheet = "TumorTissues",
    label = "tumor tissues",
    out_dir = s3_dir,
    stem = "Figure_S3a_DDR_pathway_matrix_tumor_tissues",
    title = "DDR pathway matrix — tumor tissues (192 proteins)"
  ),
  normal_tissue = list(
    key = "normal_tissue", sheet = "NonTumorTissues",
    label = "non-tumor tissues",
    out_dir = s3_dir,
    stem = "Figure_S3b_DDR_pathway_matrix_non_tumor_tissues",
    title = "DDR pathway matrix — non-tumor tissues (183 proteins)"
  ),
  cancer_cells = list(
    key = "cancer_cells", sheet = "CancerCellLines",
    label = "cancer cell lines",
    out_dir = s4_dir,
    stem = "Figure_S4a_DDR_pathway_matrix_cancer_cell_lines",
    title = "DDR pathway matrix — cancer cell lines (381 proteins)"
  ),
  normal_cells = list(
    key = "normal_cells", sheet = "NormalCellLines",
    label = "normal cell lines",
    out_dir = s4_dir,
    stem = "Figure_S4b_DDR_pathway_matrix_normal_cell_lines",
    title = "DDR pathway matrix — normal cell lines (292 proteins)"
  )
)

for (spec in matrix_specs) {
  panel <- as.data.table(read_excel(ranking_path, sheet = spec$sheet))
  panel <- panel[!is.na(BaseAccession) & nzchar(trimws(BaseAccession))]
  panel[, BaseAccession := trimws(as.character(BaseAccession))]
  
  n <- nrow(panel)
  message("Processing ", spec$label, ": ", n, " proteins")
  
  long <- melt(
    panel[, c("BaseAccession", pathway_order), with = FALSE],
    id.vars = "BaseAccession", variable.name = "Pathway", value.name = "State"
  )
  long[, `:=`(
    Pathway = factor(Pathway, levels = pathway_order),
    Rank = match(BaseAccession, panel$BaseAccession)
  )]
  long[, `:=`(
    PathwayOrder = match(as.character(Pathway), pathway_display_order),
    Y = 8 - match(as.character(Pathway), pathway_display_order),
    XMin = Rank - 1,
    XMax = Rank
  )]
  long[, `:=`(YMin = Y - 0.42, YMax = Y + 0.42)]
  
  label_x <- n + max(6, ceiling(n * 0.025))
  
  p <- ggplot() +
    geom_rect(data = long[State == 0], aes(xmin = XMin, xmax = XMax, ymin = YMin, ymax = YMax), fill = zero_fill, colour = NA) +
    geom_rect(data = long[State == 1], aes(xmin = XMin, xmax = XMax, ymin = YMin, ymax = YMax, fill = as.character(Pathway)), colour = NA) +
    geom_rect(data = long[State == -1], aes(xmin = XMin, xmax = XMax, ymin = YMin, ymax = YMax), fill = suppressing_fill, colour = NA) +
    scale_fill_manual(values = pathway_colours, guide = "none") +
    geom_text(
      data = data.frame(X = label_x, Y = 4, Label = spec$label),
      aes(x = X, y = Y, label = Label), inherit.aes = FALSE,
      angle = 90, hjust = 0.5, family = publication_font, size = 4.8, colour = "#20252B"
    ) +
    scale_x_continuous(
      limits = c(0, label_x + 2),
      breaks = unique(as.integer(round(c(0, n / 2, n)))),
      expand = c(0, 0)
    ) +
    scale_y_continuous(limits = c(0.5, 7.5), breaks = 7:1, labels = pathway_display_order, expand = c(0, 0)) +
    labs(
      title = spec$title,
      x = "Lactylated DDR Protein (#)", y = NULL
    ) +
    theme_minimal(base_family = publication_font, base_size = 11) +
    theme(
      panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(colour = "#D1D5DB", linewidth = 0.38),
      axis.text.y = element_text(face = "bold", colour = "#374151", size = 10.5),
      axis.text.x = element_text(colour = "#4B5563", size = 9.5),
      axis.title.x = element_text(colour = "#20252B", size = 12.5, margin = margin(t = 8)),
      plot.title = element_text(face = "bold", size = 14, colour = "#20252B", hjust = 0.5, margin = margin(b = 6)),
      plot.margin = margin(10, 24, 10, 12),
      plot.background = element_rect(fill = "white", colour = NA)
    )
  
  out_stem <- file.path(spec$out_dir, spec$stem)
  ggsave(paste0(out_stem, ".png"), p, width = 11.0, height = 4.8, dpi = 300, bg = "white", device = ragg::agg_png)
  ggsave(paste0(out_stem, ".pdf"), p, width = 11.0, height = 4.8, bg = "white", device = cairo_pdf)
  message("Saved: ", out_stem, ".png/.pdf")
}

