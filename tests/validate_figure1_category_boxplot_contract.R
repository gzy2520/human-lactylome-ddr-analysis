#!/usr/bin/env Rscript

# Validate the four-category, eight-box Figure 1 candidate input.
# This check reads the committed candidate table and never writes figures.

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

values_path <- file.path(project_root, "data", "candidate", "figure1_sample_boxplot_values.csv")
values <- fread(values_path, check.names = FALSE)

required_values <- c(
  "PXD", "SampleGroup", "Category", "Dataset", "SampleID",
  "SourceFile", "DdrFraction", "DdrFractionPercentage"
)
stop_if(all(required_values %in% names(values)),
  "Category-level Figure 1 input is missing required columns.")
stop_if(nrow(values) == 210L,
  "Category-level Figure 1 input must contain 210 source observations.")
stop_if(nrow(values[Dataset == "Lactylome (Kla)"]) == 92L,
  "Category-level Figure 1 input must contain 92 Kla observations.")
stop_if(nrow(values[Dataset == "Whole proteome"]) == 118L,
  "Category-level Figure 1 input must contain 118 whole-proteome observations.")

category_order <- c("normal_tissue", "cancer_tissue", "normal_cells", "cancer_cells")
dataset_order <- c("Whole proteome", "Lactylome (Kla)")
stop_if(setequal(unique(values$Category), category_order),
  "Category-level Figure 1 input does not contain exactly four categories.")
stop_if(setequal(unique(values$Dataset), dataset_order),
  "Category-level Figure 1 input does not contain exactly two datasets.")

expected_counts <- data.table(
  Category = rep(category_order, each = 2L),
  Dataset = rep(dataset_order, times = length(category_order)),
  ExpectedN = c(62L, 32L, 13L, 6L, 15L, 19L, 28L, 35L)
)
actual_counts <- values[, .(ActualN = .N), by = .(Category, Dataset)]
count_check <- merge(expected_counts, actual_counts, by = c("Category", "Dataset"), all = TRUE)
stop_if(nrow(count_check) == 8L,
  "Category-level Figure 1 input does not contain eight category-dataset groups.")
stop_if(all(count_check$ExpectedN == count_check$ActualN),
  "Category-level Figure 1 group counts do not match the confirmed sample counts.")

stop_if(!anyDuplicated(values[, .(PXD, SampleGroup, Dataset, SampleID)]),
  "Category-level Figure 1 source observations contain duplicate sample keys.")
stop_if(all(nzchar(trimws(values$SampleID))),
  "Category-level Figure 1 source observations contain an empty sample ID.")
stop_if(all(nzchar(trimws(values$SourceFile))),
  "Category-level Figure 1 source observations contain an empty source path.")
stop_if(all(!startsWith(values$SourceFile, "/")),
  "Category-level Figure 1 source paths must be relative.")
stop_if(all(is.finite(values$DdrFraction)),
  "Category-level Figure 1 fractions contain non-finite values.")
stop_if(all(is.finite(values$DdrFractionPercentage)),
  "Category-level Figure 1 percentages contain non-finite values.")
stop_if(all(values$DdrFraction >= 0 & values$DdrFraction <= 1),
  "Category-level Figure 1 fractions must be between 0 and 1.")
stop_if(all(values$DdrFractionPercentage >= 0 & values$DdrFractionPercentage <= 100),
  "Category-level Figure 1 percentages must be between 0 and 100.")

message(
  "PASS: category-level Figure 1 inputs contain eight boxplot groups with 92 Kla and 118 whole-proteome observations."
)
