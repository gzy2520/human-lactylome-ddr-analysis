#!/usr/bin/env Rscript

# Validate the dataset-level, four-category Figure 1 boxplot input.
# This check reads candidate tables and never writes publication outputs.

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
values <- fread(file.path(candidate_dir, "figure1_dataset_boxplot_values.csv"), check.names = FALSE)

required_values <- c(
  "RowOrder", "PXD", "SampleGroup", "Category", "Dataset", "DatasetPXD",
  "DatasetPointID", "DatasetPointLabel", "SampleID", "ObservationType",
  "SourceFile", "ProteinCount", "DdrProteinCount", "DdrFraction",
  "DdrFractionPercentage"
)
stop_if(all(required_values %in% names(values)),
  "Dataset-level Figure 1 input is missing required columns.")
stop_if(nrow(values) == expected_group_count * 2L,
  "Dataset-level Figure 1 input must contain two modality rows per publication group.")
stop_if(uniqueN(values[, .(PXD, SampleGroup)]) == expected_group_count,
  "Dataset-level Figure 1 input does not cover the expected publication groups.")
stop_if(uniqueN(values[, .(PXD, SampleGroup, Dataset)]) == expected_group_count * 2L,
  "Dataset-level Figure 1 input does not contain both modalities for every group.")
stop_if(!anyDuplicated(values$DatasetPointID),
  "Dataset-level Figure 1 input contains duplicate point IDs.")
stop_if(!anyDuplicated(values[, .(PXD, SampleGroup, Dataset)]),
  "A category-level boxplot point is duplicated.")
stop_if(setequal(unique(values$Dataset), c("Whole proteome", "Lactylome (Kla)")),
  "Dataset-level Figure 1 input does not contain exactly two modalities.")
stop_if(setequal(unique(values$Category),
  c("normal_tissue", "cancer_tissue", "cancer_cells", "normal_cells")),
  "Dataset-level Figure 1 input does not contain exactly four categories.")
stop_if(all(values$ObservationType == "dataset_union"),
  "Category-level Figure 1 points are not dataset unions.")
stop_if(all(nzchar(trimws(values$PXD)) & nzchar(trimws(values$SampleGroup))),
  "A dataset-level Figure 1 point has an empty PXD/group key.")
stop_if(all(nzchar(trimws(values$SourceFile))),
  "A dataset-level Figure 1 source path is empty.")
stop_if(all(!startsWith(values$SourceFile, "/")),
  "A dataset-level Figure 1 source path is absolute.")
stop_if(all(is.finite(values$DdrFraction) & is.finite(values$DdrFractionPercentage)),
  "Dataset-level Figure 1 fractions contain non-finite values.")
stop_if(all(values$DdrFraction >= 0 & values$DdrFraction <= 1),
  "A dataset-level Figure 1 fraction is outside 0-1.")
stop_if(all(values$DdrFractionPercentage >= 0 & values$DdrFractionPercentage <= 100),
  "A dataset-level Figure 1 percentage is outside 0-100.")

groups <- unique(values[, .(PXD, SampleGroup, Category)])
expected_category_counts <- groups[, .(ExpectedN = .N), by = Category]
actual_counts <- values[, .(ActualN = .N), by = .(Category, Dataset)]
expected_counts <- expected_category_counts[
  rep(seq_len(nrow(expected_category_counts)), each = 2L)
]
expected_counts[, Dataset := rep(c("Whole proteome", "Lactylome (Kla)"), times = nrow(expected_category_counts))]
count_check <- merge(expected_counts, actual_counts, by = c("Category", "Dataset"), all = TRUE)
stop_if(nrow(count_check) == 8L && all(count_check$ExpectedN == count_check$ActualN),
  "Dataset-level Figure 1 category counts are inconsistent.")

message(
  "PASS: dataset-level Figure 1 input contains ", expected_group_count,
  " PXD/sample-group points per modality across four categories."
)
