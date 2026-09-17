#!/usr/bin/env Rscript
# RNA reference panels for the 31-group set.
#
# Companion to the proteome figures: the proteome records which proteins were captured,
# these panels show what the transcriptome says about the same material classes, using the
# qsmooth-smoothed log2(TPM + 0.5) reference profile and the DDR panel lifted onto it.
#
# Two figures:
#   RNA_1  DDR-panel expression across the 31 material classes, rows grouped by DDR pathway
#   RNA_2  what the lactylome's non-detection means: absent transcript versus present but
#          not captured, and the expression distribution behind each case
#
# Colours and font follow the publication renderers so the panels sit alongside them.
#
# Usage: plot_rna_reference_31group_20260917.R <project_root> [out_dir]
args <- commandArgs(TRUE)
stopifnot(length(args) >= 1L)
root <- normalizePath(args[[1L]], mustWork = TRUE)
out_dir <- if (length(args) >= 2L) args[[2L]] else file.path(root, "results", "rna_reference_31group")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
set.seed(25)

publication_font <- "Arial Unicode MS"
text_dark <- "#20252B"; text_body <- "#30343B"
strip_fills <- c(normal_tissue = "#DCEAF5", cancer_tissue = "#F6D6D6",
                 cancer_cells = "#FCE7D4", normal_cells = "#DCEFE4")
kla_colour <- "#F28E2B"; reference_colour <- "#4E79A7"
pathway_colours <- c(BER = "#4DBBD5", NER = "#F39B7F", MMR = "#8491B4", FA = "#91D1C2",
                     HR = "#3C5488", AEJ = "#00A087", NHEJ = "#E64B35")
category_levels <- c("normal_tissue", "cancer_tissue", "cancer_cells", "normal_cells")

# same devices as the publication renderers: both resolve the system font correctly,
# where the default pdf() device rejects it
save_panel <- function(plot, stem, width, height) {
  ggsave(file.path(out_dir, paste0(stem, ".png")), plot, width = width, height = height,
         dpi = 300, bg = "white", device = ragg::agg_png)
  ggsave(file.path(out_dir, paste0(stem, ".pdf")), plot, width = width, height = height,
         device = cairo_pdf, bg = "white")
}

## ---- inputs --------------------------------------------------------------

mat <- fread(file.path(root, "outputs", "20260916_qsmooth_31group", "matrices",
                       "qsmooth_A_collapsed_log2tpm.tsv.gz"))
expansion <- fread(file.path(root, "outputs", "20260916_qsmooth_31group",
                             "group_expansion_31.csv"))
ledger <- fread(file.path(root, "audit", "20260916_full_31_rna_status", "rna_31_group_status.csv"))
annotation <- fread(file.path(root, "outputs", "20260916_ddr_panel_31group",
                              "kla_ddr_annotation.csv"))
pairs <- fread(file.path(root, "outputs", "20260917_rna_assisted_ddr",
                         "group_by_protein_pairs.tsv.gz"))

rna <- as.matrix(mat[, -1L]); rownames(rna) <- mat[[1L]]
rna <- rna[, match(expansion$ReferenceKey, colnames(rna)), drop = FALSE]
colnames(rna) <- expansion$GroupID

group_meta <- ledger[, .(GroupID, PXD, SampleGroup, Category)]
group_meta[, Category := factor(Category, levels = category_levels)]
setorder(group_meta, Category, GroupID)

pathways <- c("BER", "NER", "MMR", "FA", "HR", "AEJ", "NHEJ")

## ---- RNA_1: DDR panel expression -----------------------------------------

# the Kla n DDR panel is the set the overlay resolved into the profile
panel_expr <- fread(file.path(root, "outputs", "20260916_ddr_panel_31group",
                              "kla_ddr_expression_31groups.tsv.gz"))
panel <- annotation[BaseAccession %in% panel_expr$BaseAccession &
                      EnsemblGeneID %in% rownames(rna),
                    c("BaseAccession", "EnsemblGeneID", pathways), with = FALSE]
# the S4 states are +1 / 0 / -1, where 0 means "no directional assignment", not "member";
# a gene can carry states in several pathways, so its first is used as the display block
panel[, PathwayLabel := apply(.SD, 1L, function(v) {
  hit <- pathways[!is.na(v) & as.character(v) %in% c("1", "-1")]
  if (length(hit)) hit[[1L]] else "unassigned"
}), .SDcols = pathways]
panel[, PathwayLabel := factor(PathwayLabel, levels = c(pathways, "unassigned"))]

sub <- rna[panel$EnsemblGeneID, group_meta$GroupID, drop = FALSE]
scaled <- t(scale(t(sub)))                      # per-gene z-score, so patterns, not level
scaled[!is.finite(scaled)] <- 0

long <- data.table(
  Gene = rep(panel$EnsemblGeneID, times = ncol(scaled)),
  Pathway = rep(panel$PathwayLabel, times = ncol(scaled)),
  GroupID = rep(colnames(scaled), each = nrow(scaled)),
  Z = as.vector(scaled)
)
long <- merge(long, group_meta[, .(GroupID, Category)], by = "GroupID")
long[, Category := factor(Category, levels = category_levels)]
long[, Pathway := factor(Pathway, levels = c(pathways, "unassigned"))]
setorder(long, Category, GroupID)

# order genes by pathway then by mean expression, so the block structure reads cleanly
gene_order <- panel[, .(Gene = EnsemblGeneID, Pathway = PathwayLabel,
                        Mean = rowMeans(sub))][order(factor(Pathway, levels = c(pathways, "unassigned")), -Mean)]
long[, Gene := factor(Gene, levels = rev(gene_order$Gene))]
group_order <- unique(long[order(Category, GroupID), .(GroupID)])
long[, GroupID := factor(GroupID, levels = group_order$GroupID)]

# rows are faceted by pathway so the blocks are labelled, and genes are ordered by mean
# expression inside each block; columns are ordered by material category
long[, Gene := factor(Gene, levels = rev(gene_order$Gene))]
long[, Category := factor(Category, levels = category_levels)]
setorder(long, Category, GroupID)
long[, GroupID := factor(GroupID, levels = unique(long$GroupID))]

pathway_strip <- setNames(pathway_colours[pathways], pathways)
rna1 <- ggplot(long, aes(GroupID, Gene, fill = Z)) +
  geom_raster() +
  facet_grid(Pathway ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_gradient2(low = "#2F6690", mid = "#F7F7F7", high = "#B2182B", midpoint = 0,
                       limits = c(-2.5, 2.5), oob = scales::squish, name = "row z") +
  scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
  labs(title = "DDR panel expression across the 31 material classes",
       subtitle = sprintf(paste0("%d Kla \u2229 DDR genes, qsmooth log2(TPM + 0.5), row-scaled; ",
                                 "columns grouped by material category, rows by DDR pathway"),
                          nrow(panel)),
       x = NULL, y = NULL) +
  theme_minimal(base_size = 8, base_family = publication_font) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6, colour = text_body),
        axis.text.y = element_blank(), panel.grid = element_blank(),
        panel.spacing.y = unit(1.2, "pt"),
        strip.text.y.left = element_text(angle = 0, size = 7.5, face = "bold",
                                         colour = "white", hjust = 0.5),
        strip.background = element_rect(fill = "#4A5158", colour = NA),
        strip.placement = "outside",
        plot.title = element_text(size = 14, face = "bold", colour = text_dark),
        plot.subtitle = element_text(size = 9, colour = text_body),
        legend.position = "right", legend.title = element_text(size = 8.5),
        legend.text = element_text(size = 8),
        plot.margin = margin(10, 14, 10, 6))
save_panel(rna1, "RNA_1_DDR_panel_expression_31groups", 9.5, 12)

## ---- RNA_2: what non-detection means -------------------------------------

status <- pairs[, .(Status = fifelse(KlaDetected, "Captured as Kla",
                              fifelse(RefDetected & ExprPct > 0.5, "Present, expressed, not Kla",
                              fifelse(ExprPct > 0.5, "Expressed, not detected",
                                      "Low transcript"))))]
status[, Status := factor(Status, levels = c("Captured as Kla", "Present, expressed, not Kla",
                                             "Expressed, not detected", "Low transcript"))]
counts <- status[, .N, by = Status][order(Status)]
counts[, Pct := 100 * N / sum(N)]
counts[, Label := sprintf("%s\n%d (%.1f%%)", Status, N, Pct)]

rna2a <- ggplot(counts, aes(Status, N, fill = Status)) +
  geom_col(width = 0.68) +
  geom_text(aes(label = sprintf("%d  (%.1f%%)", N, Pct)), vjust = -0.35, size = 3,
            family = publication_font, colour = text_dark) +
  scale_fill_manual(values = c("Captured as Kla" = kla_colour,
                               "Present, expressed, not Kla" = "#B2182B",
                               "Expressed, not detected" = reference_colour,
                               "Low transcript" = "#B8BEC7")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.16))) +
  labs(title = "What lactylome non-detection means",
       subtitle = sprintf("%d (group \u00d7 Kla \u2229 DDR protein) pairs; expression taken within each group",
                          nrow(pairs)),
       x = NULL, y = "pairs") +
  theme_minimal(base_size = 9, base_family = publication_font) +
  theme(legend.position = "none", panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 7.5, colour = text_body, lineheight = 1.1),
        plot.title = element_text(size = 12, face = "bold", colour = text_dark),
        plot.subtitle = element_text(size = 8, colour = text_body))

dens <- pairs[, .(Status = fifelse(KlaDetected, "Captured as Kla",
                            fifelse(RefDetected, "In reference proteome only", "Neither")),
                  RNA = RNA)]
dens[, Status := factor(Status, levels = c("Captured as Kla", "In reference proteome only", "Neither"))]
rna2b <- ggplot(dens, aes(RNA, fill = Status, colour = Status)) +
  geom_density(alpha = 0.42, linewidth = 0.55, adjust = 1.1) +
  scale_fill_manual(values = c("Captured as Kla" = kla_colour,
                               "In reference proteome only" = reference_colour,
                               "Neither" = "#B8BEC7"), name = NULL) +
  scale_colour_manual(values = c("Captured as Kla" = kla_colour,
                                 "In reference proteome only" = reference_colour,
                                 "Neither" = "#8A929C"), guide = "none") +
  labs(x = expression(log[2]*"(TPM + 0.5)"), y = "density",
       title = "Expression behind each detection state") +
  theme_minimal(base_size = 9, base_family = publication_font) +
  theme(legend.position = "top", legend.text = element_text(size = 7.5, colour = text_body),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 12, face = "bold", colour = text_dark),
        axis.title = element_text(size = 8.5, colour = text_body))

save_panel(rna2a, "RNA_2a_non_detection_meaning", 6.4, 4.4)
save_panel(rna2b, "RNA_2b_expression_by_detection_state", 6.4, 4.0)

## ---- companion tables ----------------------------------------------------

fwrite(panel[, .(BaseAccession, EnsemblGeneID, PathwayLabel,
                 Annotation = PathwayLabel)],
       file.path(out_dir, "RNA_1_ddr_panel_genes.csv"))
fwrite(counts[, .(Status, Pairs = N, Percent = round(Pct, 2))],
       file.path(out_dir, "RNA_2_non_detection_counts.csv"))
writeLines(trimws(capture.output(sessionInfo()), which = "right"),
           file.path(out_dir, "sessionInfo.txt"))
writeLines(c(
  "# RNA reference panels (2026-09-17)", "",
  "Companion panels to the 31-group proteome figures. They read the qsmooth-smoothed",
  "log2(TPM + 0.5) reference profile and the Kla n DDR panel lifted onto it.", "",
  "RNA_1  DDR-panel expression across the material classes, rows grouped by DDR pathway,",
  "       row-scaled so patterns rather than absolute level are visible.",
  "RNA_2  the state of every (group x panel protein) pair: captured as Kla, present in the",
  "       reference proteome but not lactylated, expressed but not detected, or low transcript.", "",
  "The RNA reference is a material-class profile from different studies than the proteome,",
  "so these panels compare classes, not paired samples.",
  "Expression is judged within each group, not against an absolute cut-off; the sweeps in",
  "outputs/20260917_rna_assisted_ddr/ show how the shares move with the cut."
), file.path(out_dir, "README.md"))
cat(sprintf("RNA_PLOTS_DONE -> %s\n", out_dir))
