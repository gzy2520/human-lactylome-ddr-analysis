#!/usr/bin/env Rscript

# Validate one-way and two-way ANOVA statistics used by the isolated review
# boxplots. This check reads inputs and does not write figures.

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1L) {
  normalizePath(args[[1L]], mustWork = TRUE)
} else {
  normalizePath(".", mustWork = TRUE)
}
stop_if <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

candidate_dir <- normalizePath(
  Sys.getenv("KLA_CANDIDATE_INPUT", unset = file.path(project_root, "data", "candidate")),
  mustWork = TRUE
)
expected_group_count <- as.integer(Sys.getenv("KLA_CANDIDATE_EXPECTED_GROUPS", unset = "30"))
source(file.path(project_root, "R", "candidate", "boxplot_significance.R"), local = TRUE)

category_order <- c("normal_tissue", "cancer_tissue", "cancer_cells", "normal_cells")
dataset_order <- c("Whole proteome", "Lactylome (Kla)")
pathway_order <- c("BER", "NER", "MMR", "FA", "HR", "AEJ", "NHEJ")

figure1 <- fread(file.path(candidate_dir, "figure1_dataset_boxplot_values.csv"))
figure1_anova <- compute_category_one_way_anova(
  figure1, category_order, dataset_order, "DdrFractionPercentage"
)
stop_if(nrow(figure1_anova) == 2L,
  "Figure 1 one-way ANOVA results must contain one row per modality.")
stop_if(setequal(figure1_anova$Dataset, dataset_order),
  "Figure 1 one-way ANOVA results do not cover both modalities.")
stop_if(all(figure1_anova$Test == "one-way ANOVA"),
  "Figure 1 does not use one-way ANOVA.")
stop_if(all(figure1_anova$N == expected_group_count),
  "Figure 1 ANOVA point counts changed.")
stop_if(all(figure1_anova$NCategory == length(category_order)),
  "Figure 1 ANOVA category counts changed.")
stop_if(all(is.finite(figure1_anova$PValue) & is.finite(figure1_anova$QValueBH)),
  "Figure 1 one-way ANOVA results contain a non-finite p or q value.")
stop_if(all(figure1_anova$Significance %in% c("****", "***", "**", "*", "ns")),
  "Figure 1 one-way ANOVA significance labels are invalid.")

pathway <- fread(file.path(candidate_dir, "pathway_summary_dataset_boxplot_values.csv"))
pathway_anova <- compute_pathway_two_way_anova(
  pathway, category_order, pathway_order,
  c("Lactylome (Kla)", "Whole proteome")
)
stop_if(nrow(pathway_anova) == 2L * length(pathway_order) * 3L,
  "Pathway two-way ANOVA results must contain three terms for each plot.")
stop_if(setequal(unique(pathway_anova$Dataset), c("Lactylome (Kla)", "Whole proteome")),
  "Pathway two-way ANOVA results do not cover both modalities.")
stop_if(setequal(unique(pathway_anova$Pathway), pathway_order),
  "Pathway two-way ANOVA results do not cover seven pathways.")
stop_if(all(pathway_anova$Test == "two-way ANOVA"),
  "Pathway summary does not use two-way ANOVA.")
stop_if(all(pathway_anova$N == expected_group_count * 2L),
  "Pathway two-way ANOVA observation counts changed.")
stop_if(all(pathway_anova$NPoint == expected_group_count),
  "Pathway two-way ANOVA point counts changed.")
stop_if(setequal(unique(pathway_anova$Term),
  c("CategoryFactor", "DirectionFactor", "CategoryFactor:DirectionFactor")),
  "Pathway two-way ANOVA terms are incomplete.")
stop_if(all(is.finite(pathway_anova$PValue) & is.finite(pathway_anova$QValueBH)),
  "Pathway two-way ANOVA results contain a non-finite p or q value.")
stop_if(all(pathway_anova$Significance %in% c("****", "***", "**", "*", "ns")),
  "Pathway two-way ANOVA significance labels are invalid.")
stop_if(all(pathway_anova[Term == "CategoryFactor", Df] == length(category_order) - 1L),
  "Pathway two-way ANOVA category degrees of freedom changed.")
stop_if(all(pathway_anova[Term == "DirectionFactor", Df] == 1L),
  "Pathway two-way ANOVA direction degrees of freedom changed.")
stop_if(all(pathway_anova[Term == "CategoryFactor:DirectionFactor", Df] ==
  (length(category_order) - 1L)),
  "Pathway two-way ANOVA interaction degrees of freedom changed.")

message(
  "PASS: Figure 1 uses one-way ANOVA and the 14 pathway plots use three-term two-way ANOVA."
)
