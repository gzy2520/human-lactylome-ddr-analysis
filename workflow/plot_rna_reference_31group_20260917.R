#!/usr/bin/env Rscript
# RNA reference panels for the 31-group set (Publication Refined).
#
# Companion to the proteome figures: the proteome records which proteins were captured,
# these panels show what the transcriptome says about the same material classes, using the
# qsmooth-smoothed log2(TPM + 0.5) reference profile and the DDR panel lifted onto it.
#
# Panels:
#   RNA_1   DDR-panel expression across the 31 material classes, faceted by DDR pathway,
#           with biological category header tracks and material labels.
#   RNA_2a  Decomposition of lactylome non-detection (11,067 pairs):
#           proteome-only 2x2, with the RNA split confined to the cell both assays missed.
#   RNA_2b  Transcript expression distribution behind each detection state with medians.
#   RNA_2   Unified publication-ready composite panel combining 2a and 2b.
#
# Usage: plot_rna_reference_31group_20260917.R [project_root] [out_dir]

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(grid)
  library(gtable)
  library(patchwork)
  library(scales)
})

set.seed(25)

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1L) normalizePath(args[[1L]], mustWork = TRUE) else normalizePath(".", mustWork = TRUE)
out_dir <- if (length(args) >= 2L) args[[2L]] else file.path(root, "results", "rna_reference_31group")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

publication_font <- "Arial Unicode MS"
text_dark <- "#20252B"
text_body <- "#30343B"

read_table_gz <- function(p) {
  if (endsWith(p, ".gz")) data.table::fread(cmd = paste("gzcat", shQuote(p)))
  else data.table::fread(p)
}

save_panel <- function(plot, stem, width, height) {
  ggsave(file.path(out_dir, paste0(stem, ".png")), plot, width = width, height = height,
         dpi = 300, bg = "white", device = ragg::agg_png)
  ggsave(file.path(out_dir, paste0(stem, ".pdf")), plot, width = width, height = height,
         device = cairo_pdf, bg = "white")
}

## ---- Inputs --------------------------------------------------------------

mat <- read_table_gz(file.path(root, "outputs", "20260916_qsmooth_31group", "matrices",
                               "qsmooth_A_collapsed_log2tpm.tsv.gz"))
expansion <- fread(file.path(root, "outputs", "20260916_qsmooth_31group", "group_expansion_31.csv"))
ledger <- fread(file.path(root, "audit", "20260916_full_31_rna_status", "rna_31_group_status.csv"))
annotation <- fread(file.path(root, "outputs", "20260916_ddr_panel_31group", "kla_ddr_annotation.csv"))
panel_expr <- read_table_gz(file.path(root, "outputs", "20260916_ddr_panel_31group",
                                      "kla_ddr_expression_31groups.tsv.gz"))
pairs <- read_table_gz(file.path(root, "outputs", "20260917_rna_assisted_ddr",
                                 "group_by_protein_pairs.tsv.gz"))

rna <- as.matrix(mat[, -1L]); rownames(rna) <- mat[[1L]]
rna <- rna[, match(expansion$ReferenceKey, colnames(rna)), drop = FALSE]
colnames(rna) <- expansion$GroupID

category_levels <- c("normal_tissue", "cancer_tissue", "cancer_cells", "normal_cells")
group_meta <- ledger[, .(GroupID, PXD, SampleGroup, Category)]
group_meta[, Category := factor(Category, levels = category_levels)]
setorder(group_meta, Category, GroupID)

col_labels <- c(
  KLA31_01 = "KLA31_01: Rotator tendon",
  KLA31_02 = "KLA31_02: Normal lung",
  KLA31_03 = "KLA31_03: Hypertrophic scar",
  KLA31_04 = "KLA31_04: Adjacent skin",
  KLA31_05 = "KLA31_05: Hippocampus",
  KLA31_06 = "KLA31_06: Placenta",
  KLA31_07 = "KLA31_07: Sperm",
  KLA31_08 = "KLA31_08: BPH",
  KLA31_09 = "KLA31_09: Adjacent liver",
  KLA31_10 = "KLA31_10: ESCC",
  KLA31_11 = "KLA31_11: Prostate cancer",
  KLA31_12 = "KLA31_12: HCC",
  KLA31_13 = "KLA31_13: MCF7",
  KLA31_14 = "KLA31_14: HCT116 (a)",
  KLA31_15 = "KLA31_15: HCT116 (b)",
  KLA31_16 = "KLA31_16: HCT116 (c)",
  KLA31_17 = "KLA31_17: TALL-104",
  KLA31_18 = "KLA31_18: HepG2",
  KLA31_19 = "KLA31_19: A549",
  KLA31_20 = "KLA31_20: MDA-MB-468",
  KLA31_21 = "KLA31_21: T-47D",
  KLA31_22 = "KLA31_22: PC-3M",
  KLA31_23 = "KLA31_23: GSC (MES28)",
  KLA31_24 = "KLA31_24: RKO",
  KLA31_25 = "KLA31_25: HEK293T",
  KLA31_26 = "KLA31_26: HMC3",
  KLA31_27 = "KLA31_27: HK-2 (a)",
  KLA31_28 = "KLA31_28: HK-2 (b)",
  KLA31_29 = "KLA31_29: MCF10A",
  KLA31_30 = "KLA31_30: NSC",
  KLA31_31 = "KLA31_31: HUVEC"
)

pathways <- c("HR", "BER", "NER", "FA", "NHEJ", "MMR", "AEJ")
pathway_colours <- c(
  BER = "#4DBBD5", 
  NER = "#E64B35", 
  MMR = "#8491B4", 
  FA = "#00A087", 
  HR = "#3C5488", 
  AEJ = "#F39B7F", 
  NHEJ = "#91D1C2", 
  unassigned = "#6C757D"
)

panel <- annotation[BaseAccession %in% panel_expr$BaseAccession &
                      EnsemblGeneID %in% rownames(rna),
                    c("BaseAccession", "EnsemblGeneID", pathways), with = FALSE]

panel[, PathwayLabel := apply(.SD, 1L, function(v) {
  hit <- pathways[!is.na(v) & as.character(v) %in% c("1", "-1")]
  if (length(hit)) hit[[1L]] else "unassigned"
}), .SDcols = pathways]

pathway_order <- c("HR", "BER", "NER", "FA", "NHEJ", "MMR", "AEJ", "unassigned")
panel[, PathwayLabel := factor(PathwayLabel, levels = pathway_order)]

sub <- rna[panel$EnsemblGeneID, group_meta$GroupID, drop = FALSE]
scaled <- t(scale(t(sub)))
scaled[!is.finite(scaled)] <- 0

gene_order <- panel[, .(Gene = EnsemblGeneID, Pathway = PathwayLabel,
                        Mean = rowMeans(sub))][order(factor(Pathway, levels = pathway_order), -Mean)]

long <- data.table(
  Gene = rep(panel$EnsemblGeneID, times = ncol(scaled)),
  Pathway = rep(panel$PathwayLabel, times = ncol(scaled)),
  GroupID = rep(colnames(scaled), each = nrow(scaled)),
  Z = as.vector(scaled)
)
long <- merge(long, group_meta[, .(GroupID, Category)], by = "GroupID")
long[, Category := factor(Category, levels = category_levels)]
long[, Pathway := factor(Pathway, levels = pathway_order)]
long[, Gene := factor(Gene, levels = rev(gene_order$Gene))]
long[, GroupID := factor(GroupID, levels = group_meta$GroupID)]
setorder(long, Category, GroupID)

## ---- RNA_1: DDR panel expression heatmap ---------------------------------

p_heat <- ggplot(long, aes(GroupID, Gene, fill = Z)) +
  geom_raster() +
  geom_vline(xintercept = c(9.5, 12.5, 24.5), colour = "white", linewidth = 1.0) +
  facet_grid(Pathway ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_gradientn(
    colours = c("#2166AC", "#4393C3", "#92C5DE", "#E0F3F8", "#F7F7F7", "#FEE0D2", "#FC9272", "#DE2D26", "#A50F15"),
    limits = c(-2.5, 2.5), oob = scales::squish, name = "Row Z-score"
  ) +
  scale_x_discrete(expand = c(0, 0), labels = col_labels) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(title = "DDR panel expression across the 31 material classes",
       subtitle = sprintf("%d Kla \u2229 DDR genes, qsmooth log2(TPM + 0.5), row-scaled; columns: 31 biological materials, rows: DDR pathways", nrow(panel)),
       x = NULL, y = NULL) +
  theme_minimal(base_size = 9, base_family = publication_font) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7.5, colour = text_body, face = "plain"),
    axis.text.y = element_blank(), 
    axis.ticks.x = element_line(colour = "#B0B7C0", linewidth = 0.4),
    panel.grid = element_blank(),
    panel.spacing.y = unit(1.8, "pt"),
    strip.text.y.left = element_text(angle = 0, size = 8.5, face = "bold", colour = "white", hjust = 0.5),
    strip.background = element_rect(fill = "#555555", colour = "white", linewidth = 0.5),
    strip.placement = "outside",
    plot.title = element_text(size = 14, face = "bold", colour = text_dark, margin = margin(b = 4)),
    plot.subtitle = element_text(size = 9.5, colour = "#5A626A", margin = margin(b = 8)),
    legend.position = "right",
    legend.title = element_text(size = 9, face = "bold", colour = text_dark),
    legend.text = element_text(size = 8, colour = text_body),
    legend.key.height = unit(2.2, "cm"),
    legend.key.width = unit(0.4, "cm"),
    plot.margin = margin(12, 16, 10, 8)
  )

# Convert to gtable and recolor strips dynamically
g <- ggplotGrob(p_heat)
strips <- grep("strip-l", g$layout$name)
for (idx in strips) {
  strip_box <- g$grobs[[idx]]$grobs[[1]]
  txt_name <- grep("text", names(strip_box$children), value = TRUE)
  lbl <- strip_box$children[[txt_name]]$children[[1]]$label
  rect_name <- grep("background", names(strip_box$children), value = TRUE)
  if (lbl %in% names(pathway_colours)) {
    strip_box$children[[rect_name]]$gp$fill <- pathway_colours[[lbl]]
    strip_box$children[[rect_name]]$gp$col <- "white"
  }
  if (lbl == "AEJ") {
    strip_box$children[[txt_name]]$children[[1]]$gp$fontsize <- 6.5
  } else if (lbl == "MMR") {
    strip_box$children[[txt_name]]$children[[1]]$gp$fontsize <- 7.5
  }
  g$grobs[[idx]]$grobs[[1]] <- strip_box
}

# Add top category banner matching the 4 biological classes
top_anno_grob <- grobTree(
  rectGrob(x = unit(4.5/31, "npc"), width = unit(9/31, "npc"), y = unit(0.5, "npc"), height = unit(0.82, "npc"),
           gp = gpar(fill = "#3A7D65", col = "white", lwd = 1.5)),
  rectGrob(x = unit(10.5/31, "npc"), width = unit(3/31, "npc"), y = unit(0.5, "npc"), height = unit(0.82, "npc"),
           gp = gpar(fill = "#B03A2E", col = "white", lwd = 1.5)),
  rectGrob(x = unit(18/31, "npc"), width = unit(12/31, "npc"), y = unit(0.5, "npc"), height = unit(0.82, "npc"),
           gp = gpar(fill = "#D35400", col = "white", lwd = 1.5)),
  rectGrob(x = unit(27.5/31, "npc"), width = unit(7/31, "npc"), y = unit(0.5, "npc"), height = unit(0.82, "npc"),
           gp = gpar(fill = "#5B5EA6", col = "white", lwd = 1.5)),
  textGrob("Non-tumor tissues (n=9)", x = unit(4.5/31, "npc"), y = unit(0.5, "npc"),
           gp = gpar(col = "white", fontface = "bold", fontsize = 8, fontfamily = publication_font)),
  textGrob("Tumor tissues (n=3)", x = unit(10.5/31, "npc"), y = unit(0.5, "npc"),
           gp = gpar(col = "white", fontface = "bold", fontsize = 8, fontfamily = publication_font)),
  textGrob("Cancer cell lines (n=12)", x = unit(18/31, "npc"), y = unit(0.5, "npc"),
           gp = gpar(col = "white", fontface = "bold", fontsize = 8, fontfamily = publication_font)),
  textGrob("Normal cell lines (n=7)", x = unit(27.5/31, "npc"), y = unit(0.5, "npc"),
           gp = gpar(col = "white", fontface = "bold", fontsize = 8, fontfamily = publication_font))
)

panel_col <- g$layout$l[g$layout$name == "panel-1-1"]
top_panel_row <- min(g$layout$t[grepl("^panel", g$layout$name)])

g <- gtable_add_rows(g, heights = unit(0.65, "cm"), pos = top_panel_row - 1)
g <- gtable_add_grob(g, top_anno_grob, t = top_panel_row, b = top_panel_row, l = panel_col, r = panel_col, z = 10, name = "top-category-annotation")

# Save RNA_1
ggsave(file.path(out_dir, "RNA_1_DDR_panel_expression_31groups.png"), g, width = 11.0, height = 13.5,
       dpi = 300, bg = "white", device = ragg::agg_png)
ggsave(file.path(out_dir, "RNA_1_DDR_panel_expression_31groups.pdf"), g, width = 11.0, height = 13.5,
       device = cairo_pdf, bg = "white")

## ---- RNA_2: what non-detection means -------------------------------------

# The proteome alone already gives a 2x2: the lactylome (Kla) against the matched
# non-enriched whole-proteome reference (Ref). Its Kla-/Ref+ cell - protein present, not
# lactylated - needs no transcriptome at all. What the transcriptome adds is confined to the
# Kla-/Ref- cell, where both assays missed the protein: it separates "the protein is genuinely
# not there" from "the transcript is there and both assays still missed it".
CELL_NEITHER <- "Kla- / Ref-\nmissed by both assays"
pairs[, Cell := fifelse(KlaDetected & RefDetected, "Kla+ / Ref+",
                 fifelse(KlaDetected & !RefDetected, "Kla+ / Ref-",
                 fifelse(!KlaDetected & RefDetected, "Kla- / Ref+\nprotein present,\nnot lactylated",
                         CELL_NEITHER)))]
pairs[, Transcript := fifelse(ExprPct > 0.5, "transcript above group median",
                                                "transcript below group median")]

cell_levels <- c("Kla+ / Ref+", "Kla+ / Ref-",
                 "Kla- / Ref+\nprotein present,\nnot lactylated", CELL_NEITHER)
segment_levels <- c("proteome result", "transcript below group median", "transcript above group median")

bars <- rbind(
  pairs[Cell != CELL_NEITHER, .(Segment = "proteome result", N = .N), by = Cell],
  pairs[Cell == CELL_NEITHER, .(Segment = as.character(Transcript), N = .N),
        by = .(Cell, Transcript)][, .(Cell, Segment, N)]
)
bars[, Segment := as.character(Segment)]
bars[, Cell := factor(Cell, levels = cell_levels)]
bars[, Segment := factor(Segment, levels = segment_levels)]
counts <- bars[order(Cell, Segment)]
counts[, Label := fifelse(Segment == "proteome result", format(N, big.mark = ","), format(N, big.mark = ","))]

cell_labels <- c("Kla+ / Ref+", "Kla+ / Ref-",
                 "Kla- / Ref+\nprotein present,\nnot lactylated", CELL_NEITHER)
counts[, DisplayLabel := factor(Cell, levels = cell_levels,
                                labels = c("Kla+ / Ref+", "Kla+ / Ref-",
                                           "Kla- / Ref+\nprotein present,\nnot lactylated",
                                           "Kla- / Ref-\nmissed by\nboth assays"))]

p_2a <- ggplot(counts, aes(DisplayLabel, N, fill = Segment)) +
  geom_col(width = 0.66) +
  geom_text(data = counts[Segment == "proteome result"],
            aes(y = N, label = Label), vjust = -0.5, size = 3.5, fontface = "bold",
            family = publication_font, colour = text_dark, show.legend = FALSE) +
  geom_text(data = counts[Cell == CELL_NEITHER & Segment != "proteome result"],
            aes(y = N, label = Label), position = position_stack(vjust = 0.5),
            size = 3.1, fontface = "bold", family = publication_font,
            colour = "white", show.legend = FALSE) +
  scale_fill_manual(values = c("proteome result" = "#8C97A3",
                               "transcript below group median" = "#C9CFD6",
                               "transcript above group median" = "#B2182B"),
                    breaks = c("transcript below group median", "transcript above group median"),
                    drop = FALSE,
                    name = "transcript level within the Kla- / Ref- cell") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "Detection status of the DDR panel across the 31 groups",
       subtitle = sprintf(paste0("%s pairs, by lactylome (Kla) and matched whole-proteome ",
                                 "reference (Ref);\nonly the last bar is split by transcript level."),
                          format(nrow(pairs), big.mark = ",")),
       x = NULL, y = "pairs") +
  theme_minimal(base_size = 9, base_family = publication_font) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 7, colour = text_body, lineheight = 1.05),
        legend.position = "top", legend.title = element_text(size = 7.2, colour = text_body),
        legend.text = element_text(size = 7.2, colour = text_body),
        plot.title = element_text(size = 12, face = "bold", colour = text_dark),
        plot.subtitle = element_text(size = 7.6, colour = text_body))

# density panel: the same pairs, split three ways by what the proteome said
dens <- pairs[, .(Status = fifelse(KlaDetected, "Captured as Kla",
                            fifelse(RefDetected, "In reference proteome only", "Neither")),
                  RNA = RNA)]
dens[, Status := factor(Status, levels = c("Captured as Kla", "In reference proteome only", "Neither"))]
status_dens_colors <- c(
  "Captured as Kla" = "#E67E22",
  "In reference proteome only" = "#2980B9",
  "Neither" = "#95A5A6"
)
medians <- dens[, .(Median = median(RNA)), by = Status]

p_2b <- ggplot(dens, aes(RNA, fill = Status, colour = Status)) +
  geom_hline(yintercept = seq(0.05, 0.20, 0.05), colour = "#EDF0F2", linewidth = 0.45) +
  geom_density(alpha = 0.22, linewidth = 0.95, adjust = 1.2) +
  geom_vline(data = medians, aes(xintercept = Median, colour = Status),
             linetype = "dashed", linewidth = 0.65, show.legend = FALSE) +
  annotate("label", x = 7.6, y = 0.21, hjust = 0, vjust = 1,
           label = sprintf("Median log\u2082(TPM + 0.5):\n  \u25cf Kla: %.2f\n  \u25cf Ref proteome: %.2f\n  \u25cf Neither: %.2f",
                           medians[Status == "Captured as Kla"]$Median,
                           medians[Status == "In reference proteome only"]$Median,
                           medians[Status == "Neither"]$Median),
           size = 2.9, family = publication_font, colour = text_dark, fill = "#F8F9FA",
           lineheight = 1.2) +
  scale_fill_manual(values = status_dens_colors, name = NULL) +
  scale_colour_manual(values = status_dens_colors, name = NULL) +
  scale_x_continuous(breaks = seq(-2, 12, 2), limits = c(-1.5, 12)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = expression(log[2]*"(TPM + 0.5)"), y = "Density",
       title = "Expression Distribution by Detection State",
       subtitle = sprintf("Transcript abundance across %s pairs; dashed lines mark the group medians",
                         format(nrow(pairs), big.mark = ","))) +
  theme_minimal(base_size = 10, base_family = publication_font) +
  theme(
    legend.position = "top",
    legend.direction = "horizontal",
    legend.text = element_text(size = 8.5, colour = text_dark),
    legend.key.size = unit(0.4, "cm"),
    legend.margin = margin(0, 0, 4, 0),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.line.x = element_line(colour = "#4A5158", linewidth = 0.5),
    axis.line.y = element_line(colour = "#4A5158", linewidth = 0.5),
    axis.text.x = element_text(size = 8.5, colour = text_body),
    axis.text.y = element_text(size = 8.5, colour = text_body),
    axis.title.x = element_text(size = 9.5, face = "bold", colour = text_dark, margin = margin(t = 6)),
    axis.title.y = element_text(size = 9.5, face = "bold", colour = text_dark, margin = margin(r = 8)),
    plot.title = element_text(size = 12.5, face = "bold", colour = text_dark, margin = margin(b = 3)),
    plot.subtitle = element_text(size = 8.5, colour = "#5A626A", margin = margin(b = 8)),
    plot.margin = margin(12, 14, 10, 10)
  )

p_combined <- p_2a + p_2b + plot_layout(widths = c(1.05, 1.15)) +
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size = 14, face = "bold", family = publication_font, colour = text_dark))

save_panel(p_2a, "RNA_2a_non_detection_meaning", 7.8, 4.9)
save_panel(p_2b, "RNA_2b_expression_by_detection_state", 6.8, 5.0)
save_panel(p_combined, "RNA_2_lactylome_transcriptome_coupling", 12.8, 5.2)

## ---- Companion tables ----------------------------------------------------

fwrite(panel[, .(BaseAccession, EnsemblGeneID, PathwayLabel,
                 Annotation = PathwayLabel)],
       file.path(out_dir, "RNA_1_ddr_panel_genes.csv"))
fwrite(counts[, .(PairStatus = as.character(Cell), Segment = as.character(Segment), Pairs = N)],
       file.path(out_dir, "RNA_2_non_detection_counts.csv"))
writeLines(trimws(capture.output(sessionInfo()), which = "right"),
           file.path(out_dir, "sessionInfo.txt"))

writeLines(c(
  "# RNA reference panels (2026-09-17)", "",
  "Companion panels to the 31-group proteome figures. They read the qsmooth-smoothed",
  "log2(TPM + 0.5) reference profile and the Kla n DDR panel lifted onto it.", "",
  "RNA_1  DDR-panel expression across the 31 material classes: faceted by DDR pathway,",
  "       genes ordered by mean expression within each block, columns ordered by material",
  "       category, row-scaled so pattern rather than absolute level is visible.",
  "RNA_2a Detection status of every (group x panel protein) pair, classified by the lactylome",
  "       (Kla) and the matched whole-proteome reference (Ref). The proteome alone gives the",
  "       2x2; only the Kla-/Ref- cell is split by transcript level.",
  "RNA_2b Transcript-level density behind each detection state, with group medians.",
  "RNA_2  Composite of 2a and 2b.", "",
  "Scope: descriptive. The panels report where each pair falls and what the transcriptome",
  "adds to the one cell the proteome cannot resolve; no mechanism is inferred.",
  "The RNA reference is a material-class profile drawn from different studies than the",
  "proteome, so groups are compared as classes, not as paired samples.",
  "Expression is judged within each group rather than against an absolute cut-off; the",
  "sweeps in outputs/20260917_rna_assisted_ddr/ show how the splits move with the cut."
), file.path(out_dir, "README.md"))

# Copy to delivery folder
delivery_dir <- file.path(root, "outputs", "20260917_delivery_31group", "rna")
if (dir.exists(delivery_dir)) {
  file.copy(list.files(out_dir, full.names = TRUE), delivery_dir, overwrite = TRUE)
}

cat(sprintf("ALL_RNA_PLOTS_BEAUTIFIED_AND_SAVED -> %s\n", out_dir))
