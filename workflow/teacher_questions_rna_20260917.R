#!/usr/bin/env Rscript
if (Sys.getenv("KLA_ALLOW_LEGACY_RNA") != "1") {
  stop("Historical RNA renderer (20-gene panel/31-row expansion). Use rebuild_rna_release_20260919.sh. For historical reproduction only, set KLA_ALLOW_LEGACY_RNA=1 and choose a separate empty output directory.", call. = FALSE)
}
# The three questions the supervisor asked, answered from the transcriptome side.
#
#   1. cell proliferation rate
#   2. expression of the lactylation writer / eraser / reader enzymes
#   3. expression of DDR genes and of the genes behind the lactylated proteins
#
# Questions 1 and 2 are the priorities. Everything here is descriptive: values per material
# class, no test between groups and no mechanism claimed.
#
# Usage: teacher_questions_rna_20260917.R <project_root> [out_dir]
args <- commandArgs(TRUE)
stopifnot(length(args) >= 1L)
root <- normalizePath(args[[1L]], mustWork = TRUE)
out_dir <- if (length(args) >= 2L) args[[2L]] else file.path(root, "results", "rna_teacher_questions")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
set.seed(25)

publication_font <- "Arial Unicode MS"
text_dark <- "#20252B"; text_body <- "#30343B"
category_levels <- c("normal_tissue", "cancer_tissue", "cancer_cells", "normal_cells")
category_colours <- c(normal_tissue = "#4E79A7", cancer_tissue = "#B2182B",
                      cancer_cells = "#F28E2B", normal_cells = "#00A087")
role_colours <- c(Writer = "#2166AC", "Writer-Eraser" = "#762A83",
                  Eraser = "#B2182B", Reader = "#1B7837")

save_panel <- function(plot, stem, width, height) {
  ggsave(file.path(out_dir, paste0(stem, ".png")), plot, width = width, height = height,
         dpi = 300, bg = "white", device = ragg::agg_png)
  ggsave(file.path(out_dir, paste0(stem, ".pdf")), plot, width = width, height = height,
         device = cairo_pdf, bg = "white")
}

## ---- inputs --------------------------------------------------------------

# Cross-tissue matrices carry the HGNC rename fix (2026-09-18); see build_symbol_lookup().
qsmooth_dir <- file.path("outputs", "20260918_qsmooth_31group_hgnc")
mat <- fread(file.path(root, qsmooth_dir, "matrices",
                       "qsmooth_A_collapsed_log2tpm.tsv.gz"))
expansion <- fread(file.path(root, qsmooth_dir, "group_expansion_31.csv"))
ledger <- fread(file.path(root, "audit", "20260916_full_31_rna_status", "rna_31_group_status.csv"))
mapping <- fread(file.path(root, "outputs", "20260916_ddr_panel_31group",
                           "ddr_uniprot_to_ensembl.tsv"))
regulators <- fread(file.path(root, "data", "publication_input",
                              "regulator_kla_percentiles_31.csv"))
panel_expr <- fread(file.path(root, "outputs", "20260916_ddr_panel_31group",
                              "kla_ddr_expression_31groups.tsv.gz"))

rna <- as.matrix(mat[, -1L]); rownames(rna) <- mat[[1L]]
rna <- rna[, match(expansion$ReferenceKey, colnames(rna)), drop = FALSE]
colnames(rna) <- expansion$GroupID

group_meta <- unique(ledger[, .(GroupID, PXD, SampleGroup, Category)])
group_meta[, Category := factor(Category, levels = category_levels)]
setorder(group_meta, Category, GroupID)
group_meta[, SampleLabel := gsub("_", " ", substr(SampleGroup, 1, 26))]
# the three HCT116 rows and the two HK-2 rows share a matrix, so their sample labels repeat;
# the axis label carries the group id to keep every row distinct
group_meta[, ShortLabel := sprintf("%s  %s", GroupID, SampleLabel)]

# symbol -> Ensembl, through the official NCBI tables only
sym_tab <- fread(file.path(root, "metadata", "annotation", "human_symbol_to_entrez.tsv"),
                 header = FALSE); names(sym_tab) <- c("symbol", "entrez")
sym_tab[, entrez := as.character(entrez)]
dup <- sym_tab[, .N, by = symbol][N > 1L, symbol]
sym_tab <- sym_tab[!symbol %in% dup][!duplicated(symbol)]
symbol_to_ensg <- function(symbols) {
  e <- sym_tab$entrez[match(symbols, sym_tab$symbol)]
  m <- unique(mapping[, .(entrez = EntrezGeneIDs, ensg = EnsemblGeneIDs)][entrez != ""])
  m <- m[, .(ensg = unlist(strsplit(ensg[1L], ";"))), by = entrez][!duplicated(entrez)]
  unname(m$ensg[match(e, m$entrez)])
}

## ---- Q1: proliferation ---------------------------------------------------

# A documented proliferation marker panel. MKI67 is shown on its own as well, because it is
# the same gene the proteome side uses for its ratio, so the two can be read side by side.
PROLIF_MARKERS <- c("MKI67", "PCNA", "TOP2A", "MCM2", "MCM3", "MCM4", "MCM5", "MCM6", "MCM7",
                    "CCNB1", "CCNB2", "CDK1", "BUB1", "AURKA", "AURKB", "PLK1", "RRM2",
                    "TYMS", "BIRC5", "UBE2C")
marker_ensg <- symbol_to_ensg(PROLIF_MARKERS)
marker_tab <- data.table(Symbol = PROLIF_MARKERS, EnsemblGeneID = marker_ensg)
marker_tab[, Available := !is.na(EnsemblGeneID) & EnsemblGeneID %in% rownames(rna)]
cat(sprintf("proliferation markers mapped into the profile: %d / %d\n",
            sum(marker_tab$Available), nrow(marker_tab)))
if (any(!marker_tab$Available)) print(marker_tab[Available == FALSE])

use <- marker_tab[Available == TRUE, EnsemblGeneID]
prof <- rna[use, , drop = FALSE]
prolif <- data.table(GroupID = colnames(prof),
                     ProlifScore = colMeans(prof),
                     MKI67 = as.numeric(rna[marker_tab[Symbol == "MKI67", EnsemblGeneID], ]),
                     NMarkers = nrow(prof))
prolif <- merge(group_meta, prolif, by = "GroupID")
setorder(prolif, -ProlifScore)

cat("\nQ1 proliferation score, highest and lowest:\n")
print(prolif[, .(GroupID, ShortLabel, Category,
                 ProlifScore = round(ProlifScore, 2), MKI67 = round(MKI67, 2))][c(1:5, 27:31)])

p_q1 <- ggplot(prolif, aes(reorder(ShortLabel, ProlifScore), ProlifScore, fill = Category)) +
  geom_col(width = 0.72) +
  coord_flip() +
  scale_fill_manual(values = category_colours, name = NULL,
                    labels = c("normal tissue", "tumour tissue", "cancer cells", "normal cells")) +
  labs(title = "Q1  Proliferation across the 31 material classes",
       subtitle = sprintf(paste0("mean log2(TPM + 0.5) over %d proliferation markers; ",
                                 "MKI67 is the gene the proteome side uses as well"), nrow(prof)),
       x = NULL, y = "mean marker expression") +
  theme_minimal(base_size = 8.5, base_family = publication_font) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 7, colour = text_body),
        legend.position = "top", legend.text = element_text(size = 7.5, colour = text_body),
        plot.title = element_text(size = 12.5, face = "bold", colour = text_dark),
        plot.subtitle = element_text(size = 7.8, colour = text_body))
save_panel(p_q1, "RNA_Q1_proliferation_across_31_groups", 7.2, 8.0)

## ---- Q2: lactylation writers, erasers and readers -------------------------

reg <- unique(regulators[, .(Role, GeneSymbol, BaseAccession = RegulatorBaseAccession,
                             Display = RegulatorDisplayName)])
reg <- merge(reg, mapping[, .(BaseAccession, EnsemblGeneIDs, NGeneIDs)],
             by = "BaseAccession", all.x = TRUE)
reg[, Available := !is.na(NGeneIDs) & NGeneIDs == 1L & EnsemblGeneIDs %in% rownames(rna)]
cat(sprintf("\nlactylation regulators mapped into the profile: %d / %d\n",
            sum(reg$Available), nrow(reg)))
print(reg[, .(mapped = sum(Available), total = .N), by = Role])

reg_ok <- reg[Available == TRUE]
reg_expr <- rna[reg_ok$EnsemblGeneIDs, , drop = FALSE]
reg_long <- data.table(
  EnsemblGeneID = rep(reg_ok$EnsemblGeneIDs, times = ncol(reg_expr)),
  Symbol = rep(reg_ok$GeneSymbol, times = ncol(reg_expr)),
  Role = rep(reg_ok$Role, times = ncol(reg_expr)),
  GroupID = rep(colnames(reg_expr), each = nrow(reg_expr)),
  Expression = as.vector(reg_expr)
)
reg_long <- merge(reg_long, group_meta[, .(GroupID, ShortLabel, Category)], by = "GroupID")
reg_long[, Role := factor(Role, levels = c("Writer", "Writer-Eraser", "Eraser", "Reader"))]
reg_long[, Category := factor(Category, levels = category_levels)]
setorder(reg_long, Role, Symbol, Category, GroupID)

save_panel(
  ggplot(reg_long, aes(GroupID, Expression, fill = Category)) +
    geom_boxplot(outlier.size = 0.5, linewidth = 0.35, width = 0.62) +
    facet_wrap(~ Role, ncol = 1, scales = "free_y", strip.position = "right") +
    scale_fill_manual(values = category_colours, name = NULL,
                      labels = c("normal tissue", "tumour tissue", "cancer cells", "normal cells")) +
    labs(title = "Q2  Lactylation enzyme expression across the 31 material classes",
         subtitle = sprintf("log2(TPM + 0.5) per group, %d regulators grouped by role",
                             nrow(reg_ok)),
         x = NULL, y = "log2(TPM + 0.5)") +
    theme_minimal(base_size = 8, base_family = publication_font) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6, colour = text_body),
          panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
          strip.text.y = element_text(size = 8, face = "bold", angle = 0),
          legend.position = "top", legend.text = element_text(size = 7.5, colour = text_body),
          plot.title = element_text(size = 12.5, face = "bold", colour = text_dark),
          plot.subtitle = element_text(size = 7.8, colour = text_body)),
  "RNA_Q2_regulator_expression_by_role", 10.5, 9.5)

# heatmap: gene x group, rows grouped by role
reg_z <- t(scale(t(reg_expr)))
reg_z[!is.finite(reg_z)] <- 0
hm <- data.table(
  EnsemblGeneID = rep(rownames(reg_z), times = ncol(reg_z)),
  Symbol = rep(reg_ok$GeneSymbol, times = ncol(reg_z)),
  Role = rep(reg_ok$Role, times = ncol(reg_z)),
  GroupID = rep(colnames(reg_z), each = nrow(reg_z)),
  Z = as.vector(reg_z)
)
hm <- merge(hm, group_meta[, .(GroupID, Category)], by = "GroupID")
hm[, Role := factor(Role, levels = c("Writer", "Writer-Eraser", "Eraser", "Reader"))]
row_order <- unique(reg_ok[order(Role, GeneSymbol), GeneSymbol])
hm[, Symbol := factor(Symbol, levels = rev(unique(row_order)))]
hm[, GroupID := factor(GroupID, levels = group_meta$GroupID)]

save_panel(
  ggplot(hm, aes(GroupID, Symbol, fill = Z)) +
    geom_raster() +
    facet_grid(Role ~ ., scales = "free_y", space = "free_y", switch = "y") +
    scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", midpoint = 0,
                         limits = c(-2, 2), oob = scales::squish, name = "row z") +
    scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
    geom_vline(xintercept = c(9.5, 12.5, 24.5), colour = "white", linewidth = 0.8) +
    labs(title = "Q2  Lactylation regulators, gene by gene",
         subtitle = "row-scaled log2(TPM + 0.5); columns are the 31 groups in category order",
         x = NULL, y = NULL) +
    theme_minimal(base_size = 8, base_family = publication_font) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6, colour = text_body),
          axis.text.y = element_text(size = 6.5, colour = text_body),
          panel.grid = element_blank(), panel.spacing.y = unit(1.2, "pt"),
          strip.text.y.left = element_text(angle = 0, size = 7.5, face = "bold", colour = "white", hjust = 0.5),
          strip.background = element_rect(fill = "#4A5158", colour = NA),
          strip.placement = "outside",
          plot.title = element_text(size = 12.5, face = "bold", colour = text_dark),
          plot.subtitle = element_text(size = 7.8, colour = text_body)),
  "RNA_Q2_regulator_heatmap", 9.0, 8.5)

## ---- Q3: DDR genes and the genes behind the lactylated proteins -----------

# three gene sets, all read on the same expression scale: the panel used for the DDR figures,
# the wider DDR GO universe, and every gene behind a lactylated protein
set_genes <- list(
  "Kla n DDR panel" = panel_expr$EnsemblGeneID,
  "DDR GO universe" = mapping[grepl("DdrGoUniverse", PanelMembership) & NGeneIDs == 1L &
                                EnsemblGeneIDs %in% rownames(rna), EnsemblGeneIDs],
  "Lactylated (Kla) genes" = mapping[grepl("KlaUnion", PanelMembership) & NGeneIDs == 1L &
                                       EnsemblGeneIDs %in% rownames(rna), EnsemblGeneIDs]
)
set_genes <- lapply(set_genes, unique)
cat("\nQ3 gene sets:\n")
for (nm in names(set_genes)) cat(sprintf("  %-24s %d genes\n", nm, length(set_genes[[nm]])))

q3_long <- rbindlist(lapply(names(set_genes), function(nm) {
  g <- set_genes[[nm]]
  data.table(GeneSet = nm, GroupID = colnames(rna),
             MedianExpression = apply(rna[g, , drop = FALSE], 2L, median),
             MeanExpression = colMeans(rna[g, , drop = FALSE]))
}))
q3_long <- merge(q3_long, group_meta[, .(GroupID, ShortLabel, Category)], by = "GroupID")
q3_long[, GeneSet := factor(GeneSet, levels = names(set_genes))]
q3_long[, Category := factor(Category, levels = category_levels)]
setorder(q3_long, GeneSet, Category, GroupID)

save_panel(
  ggplot(q3_long, aes(reorder(ShortLabel, MedianExpression), MedianExpression, fill = Category)) +
    geom_col(width = 0.72) + coord_flip() +
    facet_wrap(~ GeneSet, ncol = 1, scales = "free_x", strip.position = "right") +
    scale_fill_manual(values = category_colours, name = NULL,
                      labels = c("normal tissue", "tumour tissue", "cancer cells", "normal cells")) +
    labs(title = "Q3  Expression of DDR genes and of the genes behind lactylated proteins",
         subtitle = "median log2(TPM + 0.5) of each gene set per material class",
         x = NULL, y = "median log2(TPM + 0.5)") +
    theme_minimal(base_size = 8, base_family = publication_font) +
    theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
          axis.text.y = element_text(size = 6.5, colour = text_body),
          strip.text.y = element_text(size = 8, face = "bold", angle = 0),
          legend.position = "top", legend.text = element_text(size = 7.5, colour = text_body),
          plot.title = element_text(size = 12, face = "bold", colour = text_dark),
          plot.subtitle = element_text(size = 7.8, colour = text_body)),
  "RNA_Q3_ddr_and_lactylated_gene_expression", 8.2, 10.5)

fwrite(q3_long, file.path(out_dir, "Q3_gene_set_expression.csv"))
fwrite(reg_ok[, .(Role, GeneSymbol, EnsemblGeneID = EnsemblGeneIDs)],
       file.path(out_dir, "Q2_regulators_mapped.csv"))
fwrite(prolif, file.path(out_dir, "Q1_proliferation_score.csv"))
fwrite(marker_tab, file.path(out_dir, "Q1_marker_genes.csv"))
fwrite(reg_long, file.path(out_dir, "Q2_regulator_expression.tsv.gz"), sep = "\t")
writeLines(trimws(capture.output(sessionInfo()), which = "right"),
           file.path(out_dir, "sessionInfo.txt"))
cat("\nTEACHER_QUESTIONS_DONE\n")
