#!/usr/bin/env Rscript

# Validate the seven-figure, four-category Kla pathway-summary input.

suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args)) normalizePath(args[[1L]], mustWork = TRUE) else normalizePath(".", mustWork = TRUE)
stop_if <- function(condition, message) if (!isTRUE(condition)) stop(message, call. = FALSE)
candidate_dir <- normalizePath(Sys.getenv(
  "KLA_CANDIDATE_INPUT", unset = file.path(project_root, "data", "candidate")
), mustWork = TRUE)
values <- fread(file.path(candidate_dir, "figure1_pathway_summary_sample_boxplot_values.csv"), check.names = FALSE)

pathway_order <- c("BER", "NER", "MMR", "FA", "HR", "AEJ", "NHEJ")
category_order <- c("normal_tissue", "cancer_tissue", "normal_cells", "cancer_cells")
required <- c("PXD", "SampleGroup", "Category", "Dataset", "SampleID", "SourceFile", "Pathway", "PositiveFraction", "NegativeFraction", "SignedFraction", "KlaDdrProteinCount")
stop_if(all(required %in% names(values)), "Pathway summary input is missing required columns.")
stop_if(all(values$Dataset == "Lactylome (Kla)"), "The pathway summary must be Kla-only.")
stop_if(setequal(unique(values$Pathway), pathway_order), "The pathway summary must contain seven DDR pathways.")
stop_if(setequal(unique(values$Category), category_order), "The pathway summary must contain four biological categories.")
stop_if(!anyDuplicated(values[, .(PXD, SampleGroup, SampleID, Pathway)]),
  "Each source sample must contribute once to every pathway.")
stop_if(all(nzchar(trimws(values$SourceFile)) & !startsWith(values$SourceFile, "/")),
  "Pathway-summary source paths must be non-empty and relative.")
stop_if(all(is.finite(values$PositiveFraction) & is.finite(values$NegativeFraction) & is.finite(values$SignedFraction) &
            values$PositiveFraction >= 0 & values$PositiveFraction <= 1 &
            values$NegativeFraction >= 0 & values$NegativeFraction <= 1 &
            abs(values$SignedFraction - values$PositiveFraction + values$NegativeFraction) < 1e-12),
  "Pathway summary fractions are inconsistent.")

sample_counts <- values[, .(SampleN = uniqueN(paste(PXD, SampleGroup, SampleID, sep = "|"))), by = .(Pathway, Category)]
stop_if(nrow(sample_counts) == length(pathway_order) * length(category_order) && all(sample_counts$SampleN > 0L),
  "Every pathway figure must retain four populated category rows.")

message("PASS: Kla pathway input supports seven pathway-specific figures, four categories, and Up/Down boxes on one positive axis.")
