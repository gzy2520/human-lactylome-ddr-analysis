#!/usr/bin/env Rscript

# Validate the ANOVA contracts for the restored source-sample figure layouts.

suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args)) normalizePath(args[[1L]], mustWork = TRUE) else normalizePath(".", mustWork = TRUE)
stop_if <- function(condition, message) if (!isTRUE(condition)) stop(message, call. = FALSE)
candidate_dir <- normalizePath(Sys.getenv(
  "KLA_CANDIDATE_INPUT", unset = file.path(project_root, "data", "candidate")
), mustWork = TRUE)

source(file.path(project_root, "R", "candidate", "boxplot_significance.R"), local = TRUE)
category_order <- c("normal_tissue", "cancer_tissue", "normal_cells", "cancer_cells")
dataset_order <- c("Whole proteome", "Lactylome (Kla)")
pathway_order <- c("BER", "NER", "MMR", "FA", "HR", "AEJ", "NHEJ")

figure1 <- fread(file.path(candidate_dir, "figure1_sample_boxplot_values.csv"))
figure1_anova <- compute_figure1_sample_one_way_anova(figure1, category_order, dataset_order)
stop_if(nrow(figure1_anova) == length(category_order), "Figure 1 must contain one one-way ANOVA per category panel.")
stop_if(all(figure1_anova$Test == "one-way ANOVA"), "Figure 1 does not use one-way ANOVA.")
stop_if(all(figure1_anova$NWholeProteome > 0L & figure1_anova$NKla > 0L),
  "Every Figure 1 category must contain both modalities.")
stop_if(all(is.finite(figure1_anova$PValue) & is.finite(figure1_anova$QValueBH)),
  "Figure 1 one-way ANOVA contains a non-finite p or q value.")
stop_if(all(figure1_anova$Significance %in% c("****", "***", "**", "*", "ns")),
  "Figure 1 ANOVA star labels are invalid.")

original_anova <- compute_figure1_original_dataset_one_way_anova(figure1, dataset_order)
expected_original_rows <- uniqueN(figure1[, .(PXD, SampleGroup)])
stop_if(nrow(original_anova) == expected_original_rows,
  "The PXD-axis Figure 1 ANOVA must contain one row per PXD/sample-group.")
stop_if(all(original_anova$Test == "one-way ANOVA"),
  "The PXD-axis Figure 1 does not use one-way ANOVA.")
stop_if(all(original_anova$NWholeProteome > 0L & original_anova$NKla > 0L),
  "Every PXD-axis Figure 1 row must contain both modalities.")
stop_if(all(original_anova$Significance %in% c("****", "***", "**", "*", "ns", "NA")),
  "The PXD-axis Figure 1 ANOVA labels are invalid.")
if (any(figure1$PXD == "PXD064038" & figure1$SampleGroup == "MEC and NEC ESCC groups")) {
  escc_anova <- original_anova[
    PXD == "PXD064038" & SampleGroup == "MEC and NEC ESCC groups"
  ]
  stop_if(nrow(escc_anova) == 1L && escc_anova$NWholeProteome == 94L && escc_anova$NKla == 6L,
    "The PXD064038 PXD-axis Figure 1 ANOVA row has the wrong sample counts.")
}

pathway <- fread(file.path(candidate_dir, "figure1_pathway_summary_sample_boxplot_values.csv"))
pathway_anova <- compute_pathway_sample_two_way_anova(pathway, category_order, pathway_order)
stop_if(nrow(pathway_anova) == length(pathway_order) * 3L,
  "Seven pathway figures must contain three two-way ANOVA terms each.")
stop_if(all(pathway_anova$Dataset == "Lactylome (Kla)"), "Pathway summary must be Kla-only.")
stop_if(all(pathway_anova$Test == "two-way ANOVA"), "Pathway summary does not use two-way ANOVA.")
stop_if(setequal(pathway_anova$Term, c("CategoryFactor", "DirectionFactor", "CategoryFactor:DirectionFactor")),
  "Pathway two-way ANOVA terms are incomplete.")
stop_if(all(pathway_anova$NPoint > 0L & pathway_anova$N == pathway_anova$NPoint * 2L),
  "Pathway ANOVA must use a positive and a down fraction from every source sample.")
stop_if(all(is.finite(pathway_anova$PValue) & is.finite(pathway_anova$QValueBH)),
  "Pathway two-way ANOVA contains a non-finite p or q value.")

message("PASS: restored Figure 1 uses category and PXD-axis one-way ANOVA; seven Kla pathways use 21 two-way ANOVA terms.")
