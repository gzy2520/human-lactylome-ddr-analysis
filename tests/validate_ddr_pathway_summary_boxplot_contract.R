#!/usr/bin/env Rscript

# Validate the dataset-level DDR pathway-summary boxplot input.
# This check reads data and never writes figures.

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
input_path <- file.path(candidate_dir, "pathway_summary_dataset_boxplot_values.csv")
values <- fread(input_path, check.names = FALSE)

required_values <- c(
  "RowOrder", "PXD", "SampleGroup", "Category", "Dataset", "DatasetPXD",
  "DatasetPointID", "DatasetPointLabel", "SampleID", "SourceFile",
  "Pathway", "PositiveProteinCount", "NegativeProteinCount",
  "AnyPathwayProteinCount", "PathwayScoreMappedProteinCount",
  "DdrProteinCount", "KlaDdrProteinCount", "ReferenceDdrProteinCount",
  "PathwayScoreCoverage", "PathwayScoreBasis",
  "PositiveFraction", "NegativeFraction", "SignedFraction"
)
stop_if(all(required_values %in% names(values)),
  "Dataset-level pathway-summary input is missing required columns.")
pathway_order <- c("BER", "NER", "MMR", "FA", "HR", "AEJ", "NHEJ")
category_order <- c("normal_tissue", "cancer_tissue", "cancer_cells", "normal_cells")
dataset_order <- c("Lactylome (Kla)", "Whole proteome")

stop_if(nrow(values) == expected_group_count * 2L * length(pathway_order),
  "Dataset-level pathway-summary input has an unexpected row count.")
stop_if(uniqueN(values$DatasetPointID) == expected_group_count * 2L,
  "Dataset-level pathway-summary input has an unexpected point count.")
stop_if(setequal(unique(values$Dataset), dataset_order),
  "Dataset-level pathway-summary input does not contain exactly two modalities.")
stop_if(setequal(unique(values$Category), category_order),
  "Dataset-level pathway-summary input does not contain exactly four categories.")
stop_if(setequal(unique(values$Pathway), pathway_order),
  "Dataset-level pathway-summary input does not contain exactly seven pathways.")
stop_if(!anyDuplicated(values[, .(DatasetPointID, Pathway)]),
  "Dataset-level pathway-summary input contains duplicate point/pathway rows.")
stop_if(all(values[, .N, by = .(DatasetPointID, Pathway)]$N == 1L),
  "Each dataset-level point must contribute one row per pathway.")
stop_if(all(values$ObservationType == "dataset_union"),
  "Pathway-summary points are not dataset unions.")
stop_if(all(nzchar(trimws(values$SourceFile))),
  "Pathway-summary input contains an empty source path.")
stop_if(all(!startsWith(values$SourceFile, "/")),
  "Pathway-summary source paths must be relative.")
stop_if(all(is.finite(values$PositiveFraction) & is.finite(values$NegativeFraction) &
  is.finite(values$SignedFraction) & is.finite(values$PathwayScoreCoverage)),
  "Pathway-summary fractions contain non-finite values.")
stop_if(all(values$PositiveFraction >= 0 & values$PositiveFraction <= 1),
  "A positive pathway fraction is outside 0-1.")
stop_if(all(values$NegativeFraction >= 0 & values$NegativeFraction <= 1),
  "A negative pathway fraction is outside 0-1.")
stop_if(all(values$SignedFraction >= -1 & values$SignedFraction <= 1),
  "A signed pathway fraction is outside -1 to 1.")
stop_if(all(abs(values$SignedFraction - values$PositiveFraction + values$NegativeFraction) < 1e-12),
  "Signed pathway fractions do not match their components.")
stop_if(all(values$PathwayScoreMappedProteinCount <= values$DdrProteinCount),
  "Mapped pathway proteins exceed the modality-specific DDR denominator.")
stop_if(all(abs(values$PathwayScoreCoverage -
  values$PathwayScoreMappedProteinCount / values$DdrProteinCount) < 1e-12),
  "Pathway-score coverage does not match its numerator and denominator.")

category_counts <- unique(values[, .(PXD, SampleGroup, Category)])[, .N, by = Category]
stop_if(all(category_counts$Category %in% category_order) &&
  sum(category_counts$N) == expected_group_count,
  "Pathway-summary category point counts do not match the expected group scope.")

message(
  "PASS: dataset-level pathway-summary input covers ", expected_group_count,
  " PXD/sample-group points per modality across seven pathways."
)
