#!/usr/bin/env Rscript

# Validate the original-layout Figure 1 source-sample boxplot input.

suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args)) normalizePath(args[[1L]], mustWork = TRUE) else normalizePath(".", mustWork = TRUE)
stop_if <- function(condition, message) if (!isTRUE(condition)) stop(message, call. = FALSE)
candidate_dir <- normalizePath(Sys.getenv(
  "KLA_CANDIDATE_INPUT", unset = file.path(project_root, "data", "candidate")
), mustWork = TRUE)
values <- fread(file.path(candidate_dir, "figure1_sample_boxplot_values.csv"), check.names = FALSE)

required <- c("PXD", "SampleGroup", "Category", "Dataset", "SampleID", "ObservationType", "SourceMode", "SourceFile", "DdrFraction", "DdrFractionPercentage")
stop_if(all(required %in% names(values)), "Figure 1 sample input is missing required provenance columns.")
stop_if(setequal(unique(values$Dataset), c("Whole proteome", "Lactylome (Kla)")),
  "Figure 1 must contain whole-proteome and Kla source observations.")
stop_if(setequal(unique(values$Category), c("normal_tissue", "cancer_tissue", "normal_cells", "cancer_cells")),
  "Figure 1 must retain exactly four biological categories.")
stop_if(!anyDuplicated(values[, .(PXD, SampleGroup, Dataset, SampleID)]),
  "Figure 1 has duplicated source-sample points.")
stop_if(all(nzchar(trimws(values$SourceFile)) & !startsWith(values$SourceFile, "/")),
  "Figure 1 requires non-empty relative source paths.")
stop_if(all(is.finite(values$DdrFraction) & is.finite(values$DdrFractionPercentage) &
            values$DdrFraction >= 0 & values$DdrFraction <= 1 &
            abs(values$DdrFractionPercentage - values$DdrFraction * 100) < 1e-12),
  "Figure 1 fractions are inconsistent.")

group_counts <- values[, .N, by = .(Category, Dataset)]
stop_if(all(group_counts$N > 0L), "Every Figure 1 category/modality panel must have points.")

if (any(values$PXD == "PXD064038")) {
  escc_whole <- values[PXD == "PXD064038" & Dataset == "Whole proteome"]
  escc_kla <- values[PXD == "PXD064038" & Dataset == "Lactylome (Kla)"]
  stop_if(nrow(escc_whole) == 94L && uniqueN(escc_whole$SampleID) == 94L &&
            all(escc_whole$ObservationType == "sample") &&
            all(escc_whole$SourceMode == "external_tumor_reference_sample"),
    "PXD065830 must contribute its 94 ESCC-T source columns individually to Figure 1.")
  stop_if(nrow(escc_kla) == 6L && uniqueN(escc_kla$SampleID) == 6L,
    "PXD064038 must contribute its six Kla source samples to Figure 1.")
}

message("PASS: Figure 1 retains the four-panel source-sample boxplot contract.")
