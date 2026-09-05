#!/usr/bin/env Rscript

# Validate the combined sample-level Figure 1 input and its count record.
# This check reads committed candidate tables and never writes publication
# outputs.

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

candidate_dir <- file.path(project_root, "data", "candidate")
values <- fread(file.path(candidate_dir, "figure1_sample_boxplot_values.csv"), check.names = FALSE)
counts <- fread(file.path(candidate_dir, "biological_sample_count_record.csv"), check.names = FALSE)
registry <- fread(file.path(candidate_dir, "figure1_sample_boxplot_source_registry.csv"), check.names = FALSE)

required_values <- c(
  "RowOrder", "PXD", "SampleGroup", "Category", "DisplayGroup", "Dataset",
  "SampleID", "ConditionLabel", "SampleClass", "ObservationType", "SourceMode",
  "SourceFile", "ReferencePXD", "ProteinCount", "DdrProteinCount",
  "DdrFraction", "DdrFractionPercentage"
)
required_counts <- c(
  "RowOrder", "PXD", "SampleGroup", "KlaSampleCount", "KlaSampleIDs",
  "ReferencePXD", "ReferenceSampleCount", "ReferenceSampleIDs"
)
required_registry <- c("PXD", "SampleGroup", "SampleID", "SampleClass", "Parser", "SourceFile", "Dataset")
stop_if(all(required_values %in% names(values)), "Combined Figure 1 value table is missing required columns.")
stop_if(all(required_counts %in% names(counts)), "Sample-count record is missing required columns.")
stop_if(all(required_registry %in% names(registry)), "Combined Figure 1 source registry is missing required columns.")

stop_if(nrow(counts) == 30L, "The sample-count record must cover 30 publication groups.")
stop_if(nrow(values) == 210L, "The combined Figure 1 input must contain 210 observations.")
stop_if(nrow(values[Dataset == "Lactylome (Kla)"]) == 92L,
  "The combined Figure 1 input must contain 92 Kla observations.")
stop_if(nrow(values[Dataset == "Whole proteome"]) == 118L,
  "The combined Figure 1 input must contain 118 whole-proteome observations.")
stop_if(uniqueN(values[, .(PXD, SampleGroup)]) == 30L,
  "The combined Figure 1 input does not cover the 30-group publication scope.")
stop_if(uniqueN(values[, .(PXD, SampleGroup, Dataset)]) == 60L,
  "The combined Figure 1 input does not contain both datasets for all groups.")
stop_if(uniqueN(registry[, .(PXD, SampleGroup, Dataset)]) == 60L,
  "The source registry does not contain both datasets for all groups.")
stop_if(!anyDuplicated(values[, .(PXD, SampleGroup, Dataset, SampleID)]),
  "A combined Figure 1 observation is duplicated.")
stop_if(!anyDuplicated(registry[, .(PXD, SampleGroup, Dataset, SampleID)]),
  "A combined Figure 1 registry observation is duplicated.")
stop_if(all(nzchar(trimws(values$SampleID))), "A combined Figure 1 sample identifier is empty.")
stop_if(all(nzchar(trimws(values$SourceFile))), "A combined Figure 1 source path is empty.")
stop_if(all(!startsWith(values$SourceFile, "/")), "A combined Figure 1 source path is absolute.")
stop_if(all(!startsWith(registry$SourceFile, "/")), "A combined Figure 1 registry path is absolute.")
stop_if(all(is.finite(values$DdrFraction)), "A combined Figure 1 fraction is not finite.")
stop_if(all(is.finite(values$DdrFractionPercentage)), "A combined Figure 1 percentage is not finite.")
stop_if(all(values$DdrFraction >= 0 & values$DdrFraction <= 1),
  "A combined Figure 1 fraction is outside 0-1.")
stop_if(all(values$DdrFractionPercentage >= 0 & values$DdrFractionPercentage <= 100),
  "A combined Figure 1 percentage is outside 0-100.")

expected <- melt(
  counts,
  id.vars = c("RowOrder", "PXD", "SampleGroup"),
  measure.vars = c("KlaSampleCount", "ReferenceSampleCount"),
  variable.name = "CountType",
  value.name = "ExpectedN"
)
expected[, Dataset := fifelse(CountType == "KlaSampleCount", "Lactylome (Kla)", "Whole proteome")]
actual <- values[, .(ActualN = .N), by = .(RowOrder, PXD, SampleGroup, Dataset)]
aggregate_n <- values[, .(
  AggregateOnly = all(ObservationType == "aggregate")
), by = .(RowOrder, PXD, SampleGroup, Dataset)]
count_check <- merge(
  expected[, .(RowOrder, PXD, SampleGroup, Dataset, ExpectedN)],
  actual,
  by = c("RowOrder", "PXD", "SampleGroup", "Dataset"),
  all = TRUE
)
count_check <- merge(
  count_check,
  aggregate_n,
  by = c("RowOrder", "PXD", "SampleGroup", "Dataset"),
  all.x = TRUE
)
stop_if(nrow(count_check) == 60L, "The combined Figure 1 count check does not cover 60 dataset rows.")
count_ok <- !is.na(count_check$ExpectedN) & !is.na(count_check$ActualN) & (
  count_check$ExpectedN == count_check$ActualN |
    (count_check$ActualN == 1L & count_check$ExpectedN > 1L & count_check$AggregateOnly)
)
stop_if(all(count_ok), paste(
  "The combined Figure 1 observations do not match the sample-count record; only a single",
  "explicit aggregate observation may represent multiple source replicates."
))

message("PASS: combined Figure 1 boxplot inputs cover 92 Kla and 118 whole-proteome observations across 30 groups.")
