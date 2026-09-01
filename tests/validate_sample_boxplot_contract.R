#!/usr/bin/env Rscript

# Validate the candidate-only sample-level Figure 1 inputs.  This check reads
# only the committed candidate tables and never writes publication outputs.

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
values <- fread(file.path(candidate_dir, "sample_boxplot_values.csv"), check.names = FALSE)
reconciliation <- fread(file.path(candidate_dir, "sample_boxplot_reconciliation.csv"), check.names = FALSE)
registry <- fread(file.path(candidate_dir, "sample_boxplot_source_registry.csv"), check.names = FALSE)

required_values <- c(
  "RowOrder", "PXD", "SampleGroup", "Category", "DisplayGroup", "SampleID",
  "ConditionLabel", "SampleClass", "ObservationType", "SourceMode", "SourceFile",
  "KlaProteinCount", "KlaDdrProteinCount", "KlaDdrFractionPercentage", "ReferenceFraction"
)
required_reconciliation <- c(
  "RowOrder", "PXD", "SampleGroup", "FrozenKlaProteinCount", "FrozenKlaDdrProteinCount",
  "ObservedKlaProteinCount", "ObservedKlaDdrProteinCount", "ObservedSampleCount",
  "KlaProteinCountMatchesFrozen", "KlaDdrProteinCountMatchesFrozen", "GroupUnionStatus"
)
required_registry <- c("PXD", "SampleGroup", "SampleID", "SampleClass", "Parser", "SourceFile")
stop_if(all(required_values %in% names(values)), "Sample-level value table is missing required columns.")
stop_if(all(required_reconciliation %in% names(reconciliation)), "Reconciliation table is missing required columns.")
stop_if(all(required_registry %in% names(registry)), "Source registry is missing required columns.")

stop_if(nrow(values) == 88L, "The sample-level candidate must contain 88 observations.")
stop_if(nrow(registry) == nrow(values), "The source registry and value table have different observation counts.")
stop_if(nrow(reconciliation) == 30L, "The sample-level candidate must cover 30 publication groups.")
stop_if(!anyDuplicated(values[, .(PXD, SampleGroup, SampleID)]),
  "A sample identifier is duplicated within a publication group.")
stop_if(!anyDuplicated(registry[, .(PXD, SampleGroup, SampleID)]),
  "A source-registry sample identifier is duplicated within a publication group.")
stop_if(all(nzchar(trimws(values$SampleID))), "A sample identifier is empty.")
stop_if(all(is.finite(values$KlaDdrFractionPercentage)), "A sample fraction is not finite.")
stop_if(all(values$KlaDdrFractionPercentage >= 0 & values$KlaDdrFractionPercentage <= 100),
  "A sample fraction is outside 0-100 percent.")
stop_if(all(values$ReferenceFraction >= 0 & values$ReferenceFraction <= 100),
  "A reference fraction is outside 0-100 percent.")

expected <- data.table(
  PXD = c(
    "PXD033146", "PXD036307", "PXD046800", "PXD046800", "PXD050470", "PXD063047",
    "PXD064912", "PXD066054", "PXD075377", "PXD066054", "PXD075377", "PXD060185",
    "PXD028488", "PXD053474", "PXD066351", "PXD028488", "PXD050147", "PXD054919",
    "PXD060185", "PXD060185", "PXD063266", "PXD070007", "PXD078013", "PXD028488",
    "PXD028737", "PXD058534", "PXD078736", "PXD060185", "PXD070007", "PXD073311"
  ),
  SampleGroup = c(
    "pathological rotator cuff tendon", "normal human lung", "hypertrophic scar", "adjacent skin",
    "human hippocampus", "normal pregnancy placenta", "human sperm", "BPH", "adjacent liver",
    "prostate cancer", "HCC", "MCF7", "HCT116", "HCT116", "HCT116 control and Roseburia co-culture",
    "TALL-104", "HepG2 WT and SIRT1 or SIRT3 KO", "A549", "MDA-MB-468", "T-47D", "PC-3M",
    "glioblastoma stem cells", "RKO WT and GSK3B KO", "HEK293T", "HMC3", "pretreated HK-2",
    "HK-2 control and mannitol", "MCF10A", "neural stem cells", "HUVEC control and Pg infection"
  ),
  ExpectedN = c(1L, 6L, 4L, 4L, 3L, 3L, 3L, 5L, 1L, 5L, 1L, 1L, 1L, 1L, 2L, 1L, 9L, 3L, 1L, 1L, 3L, 6L, 4L, 1L, 2L, 1L, 6L, 1L, 2L, 6L)
)
actual <- values[, .(ActualN = .N), by = .(PXD, SampleGroup)]
counts <- merge(expected, actual, by = c("PXD", "SampleGroup"), all = TRUE)
stop_if(nrow(counts) == 30L, "Expected sample-count table does not match the 30-group scope.")
stop_if(all(counts$ExpectedN == counts$ActualN),
  "Observed sample counts do not match the source-defined sample design.")

category_counts <- unique(values[, .(PXD, SampleGroup, Category)])[, .N, by = Category]
expected_category_counts <- data.table(
  Category = c("normal_tissue", "cancer_tissue", "cancer_cells", "normal_cells"),
  N = c(9L, 2L, 12L, 7L)
)
category_check <- merge(expected_category_counts, category_counts, by = "Category", suffixes = c("_expected", "_actual"), all = TRUE)
stop_if(all(category_check$N_expected == category_check$N_actual),
  "Publication category counts are not 9/2/12/7.")

stop_if(all(reconciliation$KlaProteinCountMatchesFrozen),
  "A sample-level union does not match the frozen Kla protein count.")
stop_if(all(reconciliation$KlaDdrProteinCountMatchesFrozen),
  "A sample-level union does not match the frozen Kla-DDR protein count.")
stop_if(all(reconciliation$GroupUnionStatus == "PASS"),
  "The sample-level reconciliation contains a failed publication group.")

collapsed_types <- values[ObservationType %in% c("pool", "dataset_union"), uniqueN(paste(PXD, SampleGroup))]
stop_if(collapsed_types == 8L,
  "The candidate does not retain all eight transparent single-observation groups.")

stop_if(all(!startsWith(values$SourceFile, "/")), "A sample provenance path is absolute.")
stop_if(all(!startsWith(registry$SourceFile, "/")), "A registry provenance path is absolute.")

message("PASS: sample-level boxplot inputs cover 88 observations across 30 publication groups.")
