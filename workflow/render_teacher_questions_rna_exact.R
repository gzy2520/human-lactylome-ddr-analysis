#!/usr/bin/env Rscript
# ==============================================================================
# Master Script: Render RNA-seq Counterparts for Proteome Core Figures
# (Teacher Questions Revision — Exact Replication of Proteome Methodology)
#
# Tasks:
#   1. Cell proliferation rate: Sample-level MKI67/H3C1 boxplot (and ACTB, TUBB)
#      matching Figure 1b (8.5 x 7.0 in, 4 categories, ANOVA, seed 25 jitter).
#   2. Lactylation regulators: RNA relative percentiles heatmap across 31 groups
#      matching Figure 3a/3b layout (16.5 x 11.5 in, 48 regulators, 4 roles,
#      10 highlighted genes boxed, framed and unboxed versions).
#   3. DDR and Lactylated proteins/genes: Sample-level expression and fraction
#      boxplots matching Figure 1a layout (8.5 x 7.0 in, 4 categories, seed 25).
#
# Usage: Rscript workflow/render_teacher_questions_rna_exact.R [project_root] [out_dir]
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(grid)
  library(scales)
  library(ragg)
  library(matrixStats)
})

set.seed(25)

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1L) normalizePath(args[[1L]], mustWork = TRUE) else normalizePath(".", mustWork = TRUE)
out_dir <- if (length(args) >= 2L) args[[2L]] else file.path(root, "results", "rna_teacher_questions")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

publication_font <- "Arial Unicode MS"
charcoal <- "#2F3437"
muted_text <- "#65717D"
grid_colour <- "#D9DDE3"
mean_colour <- "#C0392B"

category_order <- c("normal_tissue", "cancer_tissue", "normal_cells", "cancer_cells")
category_labels <- c(
  normal_tissue = "non-tumor\ntissues",
  cancer_tissue = "tumor\ntissues",
  normal_cells  = "normal\ncell lines",
  cancer_cells  = "cancer\ncell lines"
)
category_fills <- c(
  normal_tissue = "#0072B2",
  cancer_tissue = "#D55E00",
  normal_cells  = "#009E73",
  cancer_cells  = "#CC79A7"
)

format_q_value <- function(q_value) {
  if (!is.finite(q_value)) return("NA")
  formatC(q_value, format = "e", digits = 2)
}

read_table_gz <- function(p) {
  if (endsWith(p, ".gz")) data.table::fread(cmd = paste("gzcat", shQuote(p)))
  else data.table::fread(p)
}

# ==============================================================================
# SECTION 1: LOAD METADATA AND SAMPLE-LEVEL RNA MATRICES
# ==============================================================================
message(">>> Loading sample metadata and 31-group status...")

qc <- fread(file.path(root, "outputs", "20260916_expression_extraction", "group_sample_qc.csv"))
status <- fread(file.path(root, "audit", "20260916_full_31_rna_status", "rna_31_group_status.csv"))
groups_31 <- fread(file.path(root, "data", "publication_input", "group_summary_31.csv"))

group2cat <- setNames(status$Category, status$GroupID)
qc[, Category := vapply(strsplit(GroupIDs, ";", fixed = TRUE), function(g) {
  match_g <- intersect(g, names(group2cat))
  if (length(match_g)) group2cat[[match_g[1L]]] else NA_character_
}, character(1))]

# Harmonize category levels
qc[, Category := factor(Category, levels = category_order)]

# ==============================================================================
# SECTION 2: PROLIFERATION RATE — MKI67 / H3C1 (AND ACTB, TUBB) BOXPLOTS
# ==============================================================================
message(">>> Task 1: Building sample-level MKI67 ratio boxplots (Figure 1b style)...")

target_prolif_genes <- c(
  MKI67 = "ENSG00000148773",
  H3C1  = "ENSG00000197061",
  ACTB  = "ENSG00000075624",
  TUBB  = "ENSG00000196230"
)

mat_dir <- file.path(root, "outputs", "20260916_expression_extraction", "matrices")
mat_files <- sort(list.files(mat_dir, pattern = "_log2tpm\\.tsv\\.gz$", full.names = TRUE))

sample_ratios_list <- list()
for (mf in mat_files) {
  cmd <- paste("gzcat", shQuote(mf))
  dt <- fread(cmd = cmd, header = TRUE, sep = "\t")
  gid_col <- names(dt)[1L]
  dt_target <- dt[get(gid_col) %in% target_prolif_genes]
  if (nrow(dt_target) == 0L) next
  
  samples <- setdiff(names(dt), gid_col)
  dt_melt <- melt(dt_target, id.vars = gid_col, variable.name = "SampleID", value.name = "Log2TPM")
  dt_melt[, Gene := names(target_prolif_genes)[match(get(gid_col), target_prolif_genes)]]
  dt_melt[, TPM := pmax(0, 2^Log2TPM - 0.5)]
  
  d_wide <- dcast(dt_melt, SampleID ~ Gene, value.var = "TPM")
  sample_ratios_list[[length(sample_ratios_list) + 1L]] <- d_wide
}

all_sample_prolif <- rbindlist(sample_ratios_list, fill = TRUE)
all_sample_prolif <- merge(all_sample_prolif, qc[, .(SampleID, Category, GroupIDs, ReferenceKey)], by = "SampleID")

denominators <- c("H3C1", "ACTB", "TUBB")
denominator_titles <- c(
  H3C1 = "MKI67 / H3C1 (Histone H3.1)",
  ACTB = "MKI67 / ACTB (beta-actin)",
  TUBB = "MKI67 / TUBB (beta-tubulin)"
)

ratio_plot_data_list <- list()
ratio_sig_list <- list()

for (denom in denominators) {
  sub_df <- all_sample_prolif[MKI67 > 0 & get(denom) > 0]
  sub_df[, Ratio := MKI67 / get(denom)]
  sub_df[, Denominator := denom]
  sub_df[, RatioLabel := paste0("MKI67 / ", denom)]
  ratio_plot_data_list[[denom]] <- sub_df
  
  # Omnibus One-Way ANOVA across 4 categories on log10(Ratio)
  sub_df[, CategoryFactor := factor(as.character(Category), levels = category_order)]
  aov_res <- summary(aov(log10(Ratio) ~ CategoryFactor, data = sub_df))[[1L]]
  f_val <- aov_res["CategoryFactor", "F value"]
  p_val <- aov_res["CategoryFactor", "Pr(>F)"]
  
  ratio_sig_list[[denom]] <- data.table(
    Denominator = denom,
    Comparison = "Four-category omnibus test",
    Test = "one-way ANOVA",
    N = nrow(sub_df),
    NCategory = uniqueN(sub_df$CategoryFactor),
    FStatistic = f_val,
    PValue = p_val
  )
}

ratio_plot_data <- rbindlist(ratio_plot_data_list)
ratio_sig <- rbindlist(ratio_sig_list)
ratio_sig[, QValueBH := p.adjust(PValue, method = "BH")]

fwrite(ratio_plot_data, file.path(out_dir, "figure1_rna_mki67_ratio_sample_values.csv"))
fwrite(ratio_sig, file.path(out_dir, "figure1_rna_mki67_ratio_significance.csv"))

ratio_ticks <- function(panel_values) {
  panel_values <- panel_values[is.finite(panel_values) & panel_values > 0]
  lower <- 10^floor(log10(min(panel_values)))
  upper <- 10^ceiling(log10(max(panel_values)))
  10^seq(log10(lower), log10(upper), by = 1)
}

for (denom in denominators) {
  panel <- copy(ratio_plot_data[Denominator == denom])
  panel[, X := match(Category, category_order)]
  
  category_counts <- panel[, .(N = .N, MaxRatio = max(Ratio)), by = .(Category, X)]
  category_counts[, label_y := 10^(log10(MaxRatio) + 0.22)]
  category_stats <- panel[, .(
    Mean = 10^mean(log10(Ratio)), # Geometric mean
    Median = median(Ratio)
  ), by = .(Category, X)]
  category_stats[, c("x_left", "x_right") := list(X - 0.20, X + 0.20)]
  
  raw_max <- max(panel$Ratio, na.rm = TRUE)
  raw_min <- min(panel$Ratio, na.rm = TRUE)
  y_min <- 10^(floor(log10(raw_min)) - 0.35)
  y_max <- 10^(log10(raw_max) + 0.65)
  
  global_q <- ratio_sig[Denominator == denom, QValueBH[[1L]]]
  subtitle_text <- paste0("Four-category one-way ANOVA q = ", format_q_value(global_q))
  
  plot_ratio <- ggplot(panel, aes(x = X, y = Ratio, fill = Category)) +
    geom_vline(xintercept = c(1.5, 2.5, 3.5), colour = "#E5E7EB", linetype = "dashed", linewidth = 0.5) +
    geom_boxplot(
      aes(group = Category),
      width = 0.58, outlier.shape = NA, colour = charcoal,
      linewidth = 0.55, median.linewidth = 0.75, alpha = 0.82, orientation = "x", na.rm = TRUE
    ) +
    geom_segment(
      data = category_stats,
      aes(x = x_left, xend = x_right, y = Median, yend = Median),
      inherit.aes = FALSE, colour = charcoal, linewidth = 0.75, lineend = "round"
    ) +
    geom_segment(
      data = category_stats,
      aes(x = x_left, xend = x_right, y = Mean, yend = Mean),
      inherit.aes = FALSE, colour = mean_colour, linewidth = 0.80, lineend = "round"
    ) +
    geom_point(
      aes(group = Category),
      position = position_jitter(width = 0.095, height = 0, seed = 25),
      shape = 21, size = 2.6, stroke = 0.50, colour = "white", alpha = 0.90, na.rm = TRUE
    ) +
    geom_text(
      data = category_counts,
      aes(x = X, y = label_y, label = paste0("n=", N)),
      inherit.aes = FALSE, family = publication_font, size = 4.2, fontface = "bold",
      colour = muted_text, vjust = -0.15
    ) +
    scale_fill_manual(values = category_fills, guide = "none", drop = FALSE) +
    scale_x_continuous(
      breaks = seq_along(category_order),
      labels = unname(category_labels[category_order]),
      limits = c(0.5, length(category_order) + 0.5),
      expand = c(0, 0)
    ) +
    scale_y_log10(
      breaks = ratio_ticks(panel$Ratio),
      labels = scales::label_scientific(digits = 2),
      limits = c(y_min, y_max),
      expand = c(0, 0)
    ) +
    labs(
      title = denominator_titles[[denom]],
      subtitle = subtitle_text,
      x = NULL,
      y = "MKI67 transcript expression ratio (log scale)"
    ) +
    theme_minimal(base_family = publication_font, base_size = 14) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(colour = grid_colour, linewidth = 0.50),
      panel.border = element_blank(),
      axis.line.x = element_line(colour = "#8C939E", linewidth = 0.60),
      axis.line.y = element_line(colour = "#8C939E", linewidth = 0.60),
      axis.text.x = element_text(size = 13.5, face = "bold", colour = charcoal, lineheight = 0.95, margin = margin(t = 6)),
      axis.text.y = element_text(size = 13.0, colour = charcoal),
      axis.title.y = element_text(size = 15.5, face = "bold", colour = charcoal, margin = margin(r = 10)),
      plot.title = element_text(size = 18.0, face = "bold", colour = charcoal, hjust = 0.5, margin = margin(b = 4)),
      plot.subtitle = element_text(size = 12.0, colour = muted_text, hjust = 0.5, margin = margin(b = 10)),
      plot.margin = margin(12, 16, 12, 12),
      plot.background = element_rect(fill = "white", colour = NA)
    ) +
    coord_cartesian(clip = "off")
  
  fig_stem <- if (denom == "H3C1") "Figure_1b_RNA_MKI67_over_H3C1_boxplot" else paste0("Figure_1_companion_RNA_MKI67_over_", denom, "_boxplot")
  ggsave(file.path(out_dir, paste0(fig_stem, ".png")), plot_ratio, width = 8.5, height = 7.0, dpi = 300, bg = "white", device = ragg::agg_png)
  ggsave(file.path(out_dir, paste0(fig_stem, ".pdf")), plot_ratio, width = 8.5, height = 7.0, bg = "white", device = cairo_pdf)
}

message(">>> Task 1 complete: Figure 1b RNA MKI67 ratio boxplots saved.")

# ==============================================================================
# SECTION 3: LACTYLATION REGULATORS — RELATIVE PERCENTILE HEATMAP (FIGURE 3 STYLE)
# ==============================================================================
message(">>> Task 2: Building RNA regulator percentiles heatmap (Figure 3 style)...")

# Official 48 regulators role order
role_order <- c("Writer", "Eraser", "Writer-Eraser", "Reader")

# Load role mapping and display names exactly matching Figure 3
role_annot_candidates <- c(
  file.path(root, "data", "publication_input", "Supplementary_Table_S5_Lactylation_Regulators.xlsx"),
  file.path(root, "data", "publication_input", "Supplementary_Table_S6_Lactylation_Regulators.xlsx"),
  "/Users/gzy2520/Desktop/renew/kla/supplementary/Supplementary_Table_S6_Lactylation_Regulators.xlsx"
)
role_annot_file <- role_annot_candidates[file.exists(role_annot_candidates)][1L]
stopifnot(length(role_annot_file) > 0L && file.exists(role_annot_file))

s6_tab <- as.data.table(readxl::read_excel(role_annot_file, sheet = "Regulator_Annotations"))
role_map_std <- unique(s6_tab[, .(
  Role = trimws(as.character(Role)),
  GeneSymbol = trimws(as.character(GeneSymbol)),
  BaseAccession = trimws(as.character(BaseAccession))
)])
role_map_std <- role_map_std[Role %in% role_order & nzchar(BaseAccession)]
role_map_std[, DisplayName := ifelse(BaseAccession == "Q92830", "GCN5 (KAT2A)", GeneSymbol)]
role_map_std[, Role := factor(Role, levels = role_order)]

reg_table <- unique(role_map_std[, .(BaseAccession, GeneSymbol)])

mapping_ddr <- fread(file.path(root, "outputs", "20260916_ddr_panel_31group", "ddr_uniprot_to_ensembl.tsv"))

# Map the 48 regulators to Ensembl IDs
reg_map <- merge(reg_table, mapping_ddr[, .(BaseAccession, EnsemblGeneIDs, EnsemblGeneViaEntrez)],
                 by = "BaseAccession", all.x = TRUE)

# Read group gene summary (contains all genes for all 28 reference datasets)
summary_dt <- fread(file.path(root, "outputs", "20260916_expression_extraction", "group_gene_summary.csv"))
all_sum_genes <- unique(summary_dt$StableGeneID)

# Resolving Ensembl IDs for each of 48 regulators
resolve_ensg <- function(ensg_str, entrez_ensg_str) {
  candidates <- c(strsplit(as.character(ensg_str), ";")[[1L]], strsplit(as.character(entrez_ensg_str), ";")[[1L]])
  candidates <- candidates[nzchar(candidates) & !is.na(candidates)]
  in_sum <- candidates[candidates %in% all_sum_genes]
  if (length(in_sum)) in_sum[1L] else NA_character_
}

reg_map[, EnsemblGeneID := mapply(resolve_ensg, EnsemblGeneIDs, EnsemblGeneViaEntrez)]

# Calculate relative percentiles (0 - 100%) within each reference dataset
summary_dt[, TotalGenesInRef := .N, by = ReferenceKey]
summary_dt[, WithinRefRank := frank(MeanLog2TPM, ties.method = "average"), by = ReferenceKey]
summary_dt[, RNARelativePercentile := 100 * (WithinRefRank - 1) / (TotalGenesInRef - 1)]

# Expand across 31 groups
group_exp <- fread(file.path(root, "outputs", "20260916_qsmooth_31group", "group_expansion_31.csv"))

expanded_list <- list()
for (i in seq_len(nrow(group_exp))) {
  gid <- group_exp$GroupID[[i]]
  ref_k <- group_exp$ReferenceKey[[i]]
  ref_sub <- summary_dt[ReferenceKey == ref_k]
  
  for (j in seq_len(nrow(reg_map))) {
    acc <- reg_map$BaseAccession[[j]]
    sym <- reg_map$GeneSymbol[[j]]
    ensg <- reg_map$EnsemblGeneID[[j]]
    pct <- if (!is.na(ensg) && ensg %in% ref_sub$StableGeneID) {
      ref_sub[StableGeneID == ensg, RNARelativePercentile]
    } else NA_real_
    
    expanded_list[[length(expanded_list) + 1L]] <- data.table(
      GroupID = gid,
      ReferenceKey = ref_k,
      RegulatorBaseAccession = acc,
      GeneSymbol = sym,
      EnsemblGeneID = ifelse(is.na(ensg), "", ensg),
      RNARelativePercentile = pct
    )
  }
}
reg_rna_dt <- rbindlist(expanded_list)
reg_rna_dt <- merge(reg_rna_dt, status[, .(GroupID, PXD, SampleGroup, Category)], by = "GroupID")

fwrite(reg_rna_dt, file.path(out_dir, "regulator_rna_percentiles_31.csv"))

# Build Figure 3 Heatmap
heatmap_category_labels <- c(
  normal_tissue = "non-tumor tissues",
  cancer_tissue = "tumor tissues",
  cancer_cells  = "cancer cell lines",
  normal_cells  = "normal cell lines"
)

# Merge with group labels from groups_31
reg_plot_data <- merge(reg_rna_dt, groups_31[, .(PXD, SampleGroup, ReferenceLabelEn, RowOrder)], by = c("PXD", "SampleGroup"))
reg_plot_data <- merge(reg_plot_data, role_map_std[, .(Role, BaseAccession, DisplayName)], by.x = "RegulatorBaseAccession", by.y = "BaseAccession")

reg_plot_data[, CategoryLabel := factor(as.character(Category), levels = category_order, labels = unname(heatmap_category_labels[category_order]))]
reg_plot_data[, RoleLabel := factor(Role, levels = role_order)]
reg_plot_data[, PlotLabel := factor(ReferenceLabelEn, levels = rev(unique(ReferenceLabelEn[order(RowOrder)])))]
reg_plot_data[, DisplayName := factor(DisplayName, levels = unique(role_map_std$DisplayName))]

highlight_genes <- c("AARS1", "ACAT2", "KRT18", "SIRT2", "PARK7", "HDAC1", "HDAC2", "BRD4", "SMARCA4", "TRIM33")

# Construct bounding boxes around the 10 highlighted genes
box_lines <- list()
for (cat_lbl in levels(reg_plot_data$CategoryLabel)) {
  sub_cat <- reg_plot_data[CategoryLabel == cat_lbl]
  n_rows <- uniqueN(sub_cat$PlotLabel)
  for (role_lbl in levels(reg_plot_data$RoleLabel)) {
    sub_panel <- sub_cat[RoleLabel == role_lbl]
    if (nrow(sub_panel) == 0L) next
    panel_genes <- levels(droplevels(sub_panel$DisplayName))
    for (g in highlight_genes) {
      if (g %in% panel_genes) {
        x_pos <- which(panel_genes == g)
        box_lines[[length(box_lines) + 1L]] <- data.frame(
          CategoryLabel = cat_lbl, RoleLabel = role_lbl,
          x = x_pos - 0.5, xend = x_pos - 0.5, y = 0.5, yend = n_rows + 0.5
        )
        box_lines[[length(box_lines) + 1L]] <- data.frame(
          CategoryLabel = cat_lbl, RoleLabel = role_lbl,
          x = x_pos + 0.5, xend = x_pos + 0.5, y = 0.5, yend = n_rows + 0.5
        )
        box_lines[[length(box_lines) + 1L]] <- data.frame(
          CategoryLabel = cat_lbl, RoleLabel = role_lbl,
          x = x_pos - 0.5, xend = x_pos + 0.5, y = 0.5, yend = 0.5
        )
        box_lines[[length(box_lines) + 1L]] <- data.frame(
          CategoryLabel = cat_lbl, RoleLabel = role_lbl,
          x = x_pos - 0.5, xend = x_pos + 0.5, y = n_rows + 0.5, yend = n_rows + 0.5
        )
      }
    }
  }
}
box_lines_df <- rbindlist(box_lines)
if (nrow(box_lines_df) > 0L) {
  box_lines_df[, CategoryLabel := factor(CategoryLabel, levels = levels(reg_plot_data$CategoryLabel))]
  box_lines_df[, RoleLabel := factor(RoleLabel, levels = levels(reg_plot_data$RoleLabel))]
}

# RNA Theme Color Palette: Elegant Emerald/Teal gradient (distinct from 3a red and 3b blue)
rna_palette <- c("#FFFFFF", "#E0F2F1", "#80CBC4", "#26A69A", "#00897B", "#004D40")
box_colour_rna <- "#00695C"

build_regulator_plot <- function(include_boxes) {
  p <- ggplot(reg_plot_data, aes(x = DisplayName, y = PlotLabel, fill = RNARelativePercentile)) +
    geom_tile(colour = "white", linewidth = 0.22)
  if (include_boxes && nrow(box_lines_df) > 0L) {
    p <- p + geom_segment(
      data = box_lines_df,
      aes(x = x, xend = xend, y = y, yend = yend),
      inherit.aes = FALSE,
      colour = box_colour_rna, linewidth = 1.15
    )
  }
  p +
    facet_grid(CategoryLabel ~ RoleLabel, scales = "free", space = "free") +
    scale_fill_gradientn(
      colours = rna_palette,
      values = scales::rescale(c(0, 20, 50, 80, 100)),
      limits = c(0, 100),
      na.value = "#D9D9D9",
      name = "Transcriptome\nrelative percentile",
      guide = guide_colourbar(
        title.position = "top",
        title.hjust = 0,
        barwidth = grid::unit(5.0, "mm"),
        barheight = grid::unit(72, "mm")
      )
    ) +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 10.5, base_family = publication_font) +
    theme(
      panel.grid = element_blank(),
      strip.text.x = element_text(face = "bold", size = 12),
      strip.text.y.right = element_text(face = "bold", size = 11.0, angle = 90, hjust = 0.5),
      strip.background = element_rect(fill = "#F2F2F2", colour = NA),
      axis.text.x = element_text(angle = 55, hjust = 1, vjust = 1, size = 10.2),
      axis.text.y = element_text(size = 9.6),
      legend.position = "right",
      legend.direction = "vertical",
      legend.title = element_text(size = 11.0, face = "bold", margin = margin(b = 8)),
      legend.text = element_text(size = 9.8),
      legend.margin = margin(0, 8, 0, 8),
      plot.margin = margin(10, 14, 10, 10)
    )
}

p_reg_framed <- build_regulator_plot(include_boxes = TRUE)
p_reg_unboxed <- build_regulator_plot(include_boxes = FALSE)

ggsave(file.path(out_dir, "Figure_3c_RNA_regulator_percentiles.png"), p_reg_framed, width = 16.5, height = 11.5, dpi = 300, bg = "white", device = ragg::agg_png)
ggsave(file.path(out_dir, "Figure_3c_RNA_regulator_percentiles.pdf"), p_reg_framed, width = 16.5, height = 11.5, bg = "white", device = cairo_pdf)

ggsave(file.path(out_dir, "Figure_3c_RNA_regulator_percentiles_no_frame.png"), p_reg_unboxed, width = 16.5, height = 11.5, dpi = 300, bg = "white", device = ragg::agg_png)
ggsave(file.path(out_dir, "Figure_3c_RNA_regulator_percentiles_no_frame.pdf"), p_reg_unboxed, width = 16.5, height = 11.5, bg = "white", device = cairo_pdf)

message(">>> Task 2 complete: Figure 3c RNA regulator percentiles heatmap saved.")

# ==============================================================================
# SECTION 4: DDR & LACTYLATED GENES — SAMPLE-LEVEL BOXPLOTS (FIGURE 1a STYLE)
# ==============================================================================
message(">>> Task 3: Building sample-level DDR and lactylated gene boxplots (Figure 1a style)...")

ddr_annot <- fread(file.path(root, "outputs", "20260916_ddr_panel_31group", "kla_ddr_annotation.csv"))
ddr_panel_genes <- unique(ddr_annot$EnsemblGeneID[!is.na(ddr_annot$EnsemblGeneID)])

kla_union_genes <- unique(mapping_ddr[grepl("KlaUnion", PanelMembership) & NGeneIDs == 1L & !is.na(EnsemblGeneIDs), EnsemblGeneIDs])

# Read full sample expression matrix
full_mat_path <- file.path(root, "outputs", "20260916_qsmooth_31group", "matrices", "unified_full_log2tpm.tsv.gz")
full_dt <- read_table_gz(full_mat_path)
full_samples <- setdiff(names(full_dt), names(full_dt)[1L])

ddr_sub <- as.matrix(full_dt[get(names(full_dt)[1L]) %in% ddr_panel_genes, full_samples, with = FALSE])
kla_sub <- as.matrix(full_dt[get(names(full_dt)[1L]) %in% kla_union_genes, full_samples, with = FALSE])
all_sub <- as.matrix(full_dt[, full_samples, with = FALSE])

# 1. Sample-level Median Expression
ddr_expr_sample <- colMedians(ddr_sub)
kla_expr_sample <- colMedians(kla_sub)

# 2. Sample-level Fraction (%) among expressed genes (TPM >= 1 -> log2(TPM + 0.5) >= 0.5849625)
tpm_cut_log2 <- log2(1 + 0.5)
expressed_mask <- all_sub >= tpm_cut_log2
tot_expressed_sample <- colSums(expressed_mask)

ddr_expressed_mask <- ddr_sub >= tpm_cut_log2
ddr_count_sample <- colSums(ddr_expressed_mask)
ddr_fraction_sample <- 100 * ddr_count_sample / tot_expressed_sample

kla_expressed_mask <- kla_sub >= tpm_cut_log2
kla_count_sample <- colSums(kla_expressed_mask)
kla_fraction_sample <- 100 * kla_count_sample / tot_expressed_sample

sample_gene_stats <- data.table(
  SampleID = full_samples,
  DDRExpressionMedian = ddr_expr_sample,
  KlaExpressionMedian = kla_expr_sample,
  DdrFractionPercentage = ddr_fraction_sample,
  KlaFractionPercentage = kla_fraction_sample
)
sample_gene_stats <- merge(sample_gene_stats, qc[, .(SampleID, Category, GroupIDs)], by = "SampleID")
fwrite(sample_gene_stats, file.path(out_dir, "figure1a_rna_ddr_and_lactylated_sample_values.csv"))

# --- 4A. Plot Figure 1a: DDR and Lactylated Gene Expression Boxplot ---
expr_long <- melt(
  sample_gene_stats,
  id.vars = c("SampleID", "Category"),
  measure.vars = c("DDRExpressionMedian", "KlaExpressionMedian"),
  variable.name = "GeneSet",
  value.name = "MedianExpression"
)
expr_long[, GeneSet := factor(GeneSet, levels = c("DDRExpressionMedian", "KlaExpressionMedian"),
                              labels = c("DDR genes", "Lactylated protein genes"))]
expr_long[, CategoryLabel := factor(Category, levels = category_order, labels = unname(category_labels[category_order]))]

# ANOVA for expression
expr_aov <- aov(MedianExpression ~ Category * GeneSet, data = expr_long)
expr_aov_tab <- summary(expr_aov)[[1L]]
cat_p <- expr_aov_tab["Category", "Pr(>F)"]
cat_f <- expr_aov_tab["Category", "F value"]
expr_sub_text <- paste0("Two-way ANOVA Category factor p = ", format_q_value(cat_p), " (F = ", sprintf("%.2f", cat_f), ")")

dodge_w <- 0.72
expr_summary <- expr_long[, .(
  Mean = mean(MedianExpression),
  Median = median(MedianExpression),
  N = .N
), by = .(Category, CategoryLabel, GeneSet)]
expr_summary[, CatIdx := as.numeric(CategoryLabel)]
expr_summary[, XPos := fifelse(GeneSet == "DDR genes", CatIdx - dodge_w / 4, CatIdx + dodge_w / 4)]

y_max_expr <- max(expr_long$MedianExpression, na.rm = TRUE) + 0.85
y_min_expr <- min(expr_long$MedianExpression, na.rm = TRUE) - 0.2

plot_fig1a_expr <- ggplot(expr_long, aes(x = CategoryLabel, y = MedianExpression, fill = GeneSet)) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5), colour = "#E5E7EB", linetype = "dashed", linewidth = 0.5) +
  geom_boxplot(
    position = position_dodge(width = dodge_w),
    width = 0.58, outlier.shape = NA, colour = charcoal,
    linewidth = 0.55, median.linewidth = 0.75, alpha = 0.82, na.rm = TRUE
  ) +
  geom_segment(
    data = expr_summary,
    aes(x = XPos - 0.16, xend = XPos + 0.16, y = Mean, yend = Mean),
    inherit.aes = FALSE, colour = mean_colour, linewidth = 0.80
  ) +
  geom_point(
    aes(fill = GeneSet),
    position = position_jitterdodge(jitter.width = 0.14, dodge.width = dodge_w, seed = 25),
    shape = 21, size = 2.4, stroke = 0.45, colour = "white", alpha = 0.88, na.rm = TRUE
  ) +
  geom_text(
    data = expr_summary,
    aes(x = XPos, y = y_max_expr * 0.96, label = paste0("n=", N)),
    inherit.aes = FALSE, size = 3.8, family = publication_font, colour = muted_text, fontface = "bold"
  ) +
  scale_fill_manual(values = c("DDR genes" = "#4E79A7", "Lactylated protein genes" = "#F28E2B")) +
  scale_y_continuous(limits = c(y_min_expr, y_max_expr), breaks = scales::pretty_breaks(n = 6), expand = expansion(mult = c(0, 0))) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE, keyheight = grid::unit(0.55, "cm"), keywidth = grid::unit(0.85, "cm"))) +
  labs(
    title = "Expression of DDR genes and lactylated protein genes across four biological categories",
    subtitle = expr_sub_text,
    x = NULL, y = "Sample median log2(TPM + 0.5)", fill = NULL
  ) +
  theme_minimal(base_size = 14, base_family = publication_font) +
  theme(
    panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(colour = grid_colour, linewidth = 0.50),
    panel.border = element_blank(),
    axis.line.x = element_line(colour = "#8C939E", linewidth = 0.60),
    axis.line.y = element_line(colour = "#8C939E", linewidth = 0.60),
    axis.text.x = element_text(size = 13.0, colour = charcoal, face = "bold", margin = margin(t = 6)),
    axis.text.y = element_text(size = 13.0, colour = charcoal),
    axis.title.y = element_text(size = 15, face = "bold", colour = charcoal, margin = margin(r = 10)),
    plot.title = element_text(size = 15.5, face = "bold", colour = charcoal, hjust = 0.5, margin = margin(b = 4)),
    plot.subtitle = element_text(size = 12.0, colour = muted_text, hjust = 0.5, margin = margin(b = 10)),
    legend.position = "top", legend.direction = "horizontal",
    legend.text = element_text(size = 13.0, colour = charcoal),
    legend.key.spacing.x = grid::unit(0.35, "cm"), legend.background = element_rect(fill = "white", colour = NA),
    legend.margin = margin(1, 0, 4, 0),
    plot.margin = margin(10, 16, 12, 12), plot.background = element_rect(fill = "white", colour = NA)
  )

stem_1a_expr <- file.path(out_dir, "Figure_1a_RNA_DDR_and_lactylated_gene_expression_boxplot")
ggsave(paste0(stem_1a_expr, ".png"), plot_fig1a_expr, width = 8.5, height = 7.0, dpi = 300, bg = "white", device = ragg::agg_png)
ggsave(paste0(stem_1a_expr, ".pdf"), plot_fig1a_expr, width = 8.5, height = 7.0, bg = "white", device = cairo_pdf)

# --- 4B. Plot Figure 1a: DDR Gene Fraction Boxplot (exact Figure 1a counterpart) ---
frac_plot_dt <- copy(sample_gene_stats)
frac_plot_dt[, CategoryLabel := factor(Category, levels = category_order, labels = unname(category_labels[category_order]))]
frac_plot_dt[, X := match(Category, category_order)]

frac_counts <- frac_plot_dt[, .(N = .N, MaxFrac = max(DdrFractionPercentage)), by = .(Category, X, CategoryLabel)]
frac_stats <- frac_plot_dt[, .(
  Mean = mean(DdrFractionPercentage),
  Median = median(DdrFractionPercentage)
), by = .(Category, X, CategoryLabel)]
frac_stats[, c("x_left", "x_right") := list(X - 0.20, X + 0.20)]

frac_aov <- summary(aov(DdrFractionPercentage ~ Category, data = frac_plot_dt))[[1L]]
frac_f <- frac_aov["Category", "F value"]
frac_p <- frac_aov["Category", "Pr(>F)"]
frac_sub_text <- paste0("Four-category one-way ANOVA p = ", format_q_value(frac_p), " (F = ", sprintf("%.2f", frac_f), ")")

plot_fig1a_frac <- ggplot(frac_plot_dt, aes(x = X, y = DdrFractionPercentage, fill = Category)) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5), colour = "#E5E7EB", linetype = "dashed", linewidth = 0.5) +
  geom_boxplot(
    aes(group = Category),
    width = 0.58, outlier.shape = NA, colour = charcoal,
    linewidth = 0.55, median.linewidth = 0.75, alpha = 0.82, orientation = "x", na.rm = TRUE
  ) +
  geom_segment(
    data = frac_stats,
    aes(x = x_left, xend = x_right, y = Median, yend = Median),
    inherit.aes = FALSE, colour = charcoal, linewidth = 0.75, lineend = "round"
  ) +
  geom_segment(
    data = frac_stats,
    aes(x = x_left, xend = x_right, y = Mean, yend = Mean),
    inherit.aes = FALSE, colour = mean_colour, linewidth = 0.80, lineend = "round"
  ) +
  geom_point(
    aes(group = Category),
    position = position_jitter(width = 0.095, height = 0, seed = 25),
    shape = 21, size = 2.6, stroke = 0.50, colour = "white", alpha = 0.90, na.rm = TRUE
  ) +
  geom_text(
    data = frac_counts,
    aes(x = X, y = 4.8, label = paste0("n=", N)),
    inherit.aes = FALSE, family = publication_font, size = 4.2, fontface = "bold",
    colour = muted_text, vjust = -0.15
  ) +
  scale_fill_manual(values = category_fills, guide = "none", drop = FALSE) +
  scale_x_continuous(
    breaks = seq_along(category_order),
    labels = unname(category_labels[category_order]),
    limits = c(0.5, length(category_order) + 0.5),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(0, 5.2),
    breaks = seq(0, 5, by = 1),
    labels = function(y) paste0(y, "%"),
    expand = c(0, 0)
  ) +
  labs(
    title = "DDR annotated gene fraction in the transcriptome across four biological categories",
    subtitle = frac_sub_text,
    x = NULL,
    y = "GO-DDR annotated expressed gene fraction (%)"
  ) +
  theme_minimal(base_family = publication_font, base_size = 14) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(colour = grid_colour, linewidth = 0.50),
    panel.border = element_blank(),
    axis.line.x = element_line(colour = "#8C939E", linewidth = 0.60),
    axis.line.y = element_line(colour = "#8C939E", linewidth = 0.60),
    axis.text.x = element_text(size = 13.5, face = "bold", colour = charcoal, lineheight = 0.95, margin = margin(t = 6)),
    axis.text.y = element_text(size = 13.0, colour = charcoal),
    axis.title.y = element_text(size = 15.5, face = "bold", colour = charcoal, margin = margin(r = 10)),
    plot.title = element_text(size = 16.0, face = "bold", colour = charcoal, hjust = 0.5, margin = margin(b = 4)),
    plot.subtitle = element_text(size = 12.0, colour = muted_text, hjust = 0.5, margin = margin(b = 10)),
    plot.margin = margin(12, 16, 12, 12),
    plot.background = element_rect(fill = "white", colour = NA)
  ) +
  coord_cartesian(clip = "off")

stem_1a_frac <- file.path(out_dir, "Figure_1a_RNA_DDR_fraction_boxplot")
ggsave(paste0(stem_1a_frac, ".png"), plot_fig1a_frac, width = 8.5, height = 7.0, dpi = 300, bg = "white", device = ragg::agg_png)
ggsave(paste0(stem_1a_frac, ".pdf"), plot_fig1a_frac, width = 8.5, height = 7.0, bg = "white", device = cairo_pdf)

message(">>> Task 3 complete: Figure 1a RNA DDR and lactylated gene boxplots saved.")

# ==============================================================================
# SECTION 5: SYNCHRONIZATION TO DELIVERY DIRECTORIES
# ==============================================================================
message(">>> Synchronizing all newly generated deliverables...")

# 1. Local delivery dir
deliv_dir <- file.path(root, "outputs", "20260917_delivery_31group", "rna")
if (dir.exists(deliv_dir)) {
  file.copy(list.files(out_dir, full.names = TRUE), deliv_dir, overwrite = TRUE)
}

# 2. External delivery dir
ext_renew_rna <- "/Users/gzy2520/Desktop/renew/kla/rna"
ext_renew_top <- "/Users/gzy2520/Desktop/renew/kla"

if (dir.exists(ext_renew_rna)) {
  # Clean up old ad-hoc files from previous session
  old_adhoc <- c(
    "RNA_Q1_proliferation_across_31_groups.png", "RNA_Q1_proliferation_across_31_groups.pdf",
    "RNA_Q2_regulator_expression_by_role.png", "RNA_Q2_regulator_expression_by_role.pdf",
    "RNA_Q2_regulator_heatmap.png", "RNA_Q2_regulator_heatmap.pdf",
    "RNA_Q3_ddr_and_lactylated_gene_expression.png", "RNA_Q3_ddr_and_lactylated_gene_expression.pdf",
    "Q1_marker_genes.csv", "Q1_proliferation_score.csv", "Q2_regulators_mapped.csv", "Q3_gene_set_expression.csv"
  )
  for (f in old_adhoc) {
    if (file.exists(file.path(ext_renew_rna, f))) unlink(file.path(ext_renew_rna, f))
    if (file.exists(file.path(ext_renew_top, f))) unlink(file.path(ext_renew_top, f))
    if (file.exists(file.path(out_dir, f))) unlink(file.path(out_dir, f))
    if (file.exists(file.path(deliv_dir, f))) unlink(file.path(deliv_dir, f))
  }
  
  # Copy all deliverables to ext_renew_rna
  file.copy(list.files(out_dir, full.names = TRUE), ext_renew_rna, overwrite = TRUE)
  
  # Copy top-level publication-ready PNGs to ext_renew_top
  top_pngs <- list.files(out_dir, pattern = "\\.png$", full.names = TRUE)
  file.copy(top_pngs, ext_renew_top, overwrite = TRUE)
}

cat("\nALL_TEACHER_QUESTIONS_REBUILT_AND_DELIVERED_SUCCESSFULLY!\n")
