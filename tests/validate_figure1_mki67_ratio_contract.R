#!/usr/bin/env Rscript

# Validate the isolated MKI67-normalized Figure 1 candidate inputs. This check
# does not touch the approved publication outputs.

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
key <- fread(file.path(candidate_dir, "figure1_mki67_ratio_protein_key.csv"), check.names = FALSE)
audit <- fread(file.path(candidate_dir, "figure1_mki67_ratio_source_audit.csv"), check.names = FALSE)
values <- fread(file.path(candidate_dir, "figure1_mki67_ratio_sample_values.csv"), check.names = FALSE)
coverage <- fread(file.path(candidate_dir, "figure1_mki67_ratio_coverage.csv"), check.names = FALSE)
significance <- fread(file.path(candidate_dir, "figure1_mki67_ratio_significance.csv"), check.names = FALSE)

target_labels <- c("MKI67", "ACTB", "TUBB", "H3C1")
target_accessions <- c(MKI67 = "P46013", ACTB = "P60709", TUBB = "P07437", H3C1 = "P68431")
denominators <- c("ACTB", "TUBB", "H3C1")
categories <- c("normal_tissue", "cancer_tissue", "normal_cells", "cancer_cells")

stop_if(nrow(key) == 4L, "The MKI67 ratio protein key must contain four entries.")
stop_if(identical(key$ProteinLabel, target_labels), "The MKI67 ratio protein-key order changed.")
stop_if(identical(key$BaseAccession, unname(target_accessions)),
  "The MKI67 ratio protein key does not use the confirmed stable accessions.")
stop_if(any(grepl("^BaseAccession$", names(key))), "The ratio protein key lacks the stable accession field.")

required_audit <- c(
  "ObsKey", "PXD", "SampleGroup", "Category", "Dataset", "SampleID", "SourceFile",
  "Parser", "QuantitationField", "QuantitationTransform", "SourceColumnStatus"
)
required_audit <- c(required_audit, unlist(lapply(target_labels, function(label) c(
  paste0(label, "_SourceValue"), paste0(label, "_Intensity"),
  paste0(label, "_RowCount"), paste0(label, "_ReportedGroup")
))))
required_audit <- c(required_audit, unlist(lapply(denominators, function(label) c(
  paste0("MKI67_over_", label), paste0("MKI67_over_", label, "_Status")
))))
stop_if(all(required_audit %in% names(audit)), "The MKI67 ratio source audit schema is incomplete.")
stop_if(nrow(audit) == 118L, "The source audit must contain the 118 whole-proteome observations.")
stop_if(uniqueN(audit$ObsKey) == 118L, "The source audit observation keys are not unique.")
stop_if(all(audit$Dataset == "Whole proteome"), "The source audit contains a non-whole-proteome record.")
stop_if(all(audit$Category %in% categories), "The source audit contains an unknown category.")
stop_if(all(nzchar(trimws(audit$SourceFile))), "The source audit contains an empty source path.")
stop_if(all(!startsWith(audit$SourceFile, "/")), "The source audit contains an absolute source path.")
stop_if(all(audit$QuantitationTransform %in% c("identity", "back_transform_2^x", "not_available")),
  "The source audit contains an undocumented quantitative transform.")

required_values <- c(
  "ObsKey", "PXD", "SampleGroup", "Category", "Dataset", "SampleID", "SourceFile",
  "Denominator", "RatioLabel", "Ratio"
)
stop_if(all(required_values %in% names(values)), "The MKI67 ratio sample-value schema is incomplete.")
stop_if(nrow(values) == 130L, "The MKI67 ratio sample-value table must contain 130 valid ratios.")
stop_if(all(values$Dataset == "Whole proteome"), "The ratio sample-value table contains a non-whole-proteome record.")
stop_if(all(values$Category %in% categories), "The ratio sample-value table contains an unknown category.")
stop_if(all(values$Denominator %in% denominators), "The ratio sample-value table contains an unknown denominator.")
stop_if(all(is.finite(values$Ratio) & values$Ratio > 0), "The ratio sample-value table contains a non-positive value.")
stop_if(uniqueN(values[, .(ObsKey, Denominator)]) == nrow(values),
  "The ratio sample-value table contains duplicate observation-denominator records.")
stop_if(all(values$RatioLabel == paste0("MKI67 / ", values$Denominator)),
  "The ratio labels do not match their denominator.")

expected_n <- c(ACTB = 44L, TUBB = 51L, H3C1 = 35L)
observed_n <- values[, .(N = .N), by = Denominator]
stop_if(all(expected_n[observed_n$Denominator] == observed_n$N),
  "The valid MKI67 ratio counts changed.")

for (denominator in denominators) {
  ratio_column <- paste0("MKI67_over_", denominator)
  denominator_column <- paste0(denominator, "_Intensity")
  status_column <- paste0(ratio_column, "_Status")
  valid <- audit[[status_column]] == "valid"
  valid[is.na(valid)] <- FALSE
  expected_ratio <- audit[[ratio_column]][valid]
  calculated_ratio <- audit$MKI67_Intensity[valid] / audit[[denominator_column]][valid]
  stop_if(all(is.finite(calculated_ratio) & calculated_ratio > 0),
    paste0("The ", denominator, " ratios contain an invalid calculated value."))
  stop_if(all(abs(expected_ratio - calculated_ratio) < 1e-12 * pmax(1, abs(calculated_ratio))),
    paste0("The ", denominator, " ratios do not match the stored intensities."))
  stored_keys <- audit$ObsKey[valid]
  value_keys <- values[Denominator == denominator, ObsKey]
  stop_if(setequal(stored_keys, value_keys),
    paste0("The ", denominator, " sample values do not match the valid audit records."))
}

stop_if(nrow(coverage) == length(denominators) * length(categories),
  "The MKI67 ratio coverage table must contain all denominator-category cells.")
stop_if(all(coverage$Denominator %in% denominators & coverage$Category %in% categories),
  "The MKI67 ratio coverage table contains an unknown key.")

expected_global_rows <- length(denominators)
stop_if(nrow(significance) == expected_global_rows,
  "The MKI67 ratio significance table has an unexpected number of rows.")
stop_if(all(significance$Denominator %in% denominators), "The significance table contains an unknown denominator.")
global_significance <- significance[Test == "one-way ANOVA"]
stop_if(nrow(global_significance) == expected_global_rows,
  "The MKI67 ratio significance table must contain one omnibus test per denominator.")
stop_if(identical(global_significance$Denominator, denominators),
  "The MKI67 ratio omnibus denominator order changed.")
stop_if(all(global_significance$AdjustmentFamily == "MKI67 ratio: three omnibus one-way ANOVA tests"),
  "The MKI67 ratio omnibus adjustment family changed.")
stop_if(all(is.finite(significance$PValue) & is.finite(significance$QValueBH)),
  "The MKI67 ratio significance table contains a missing p or q value.")
stop_if(all(significance$Significance %in% c("****", "***", "**", "*", "ns")),
  "The MKI67 ratio significance labels contain an unexpected level.")

expected_global_p <- values[, .(
  PValue = {
    fit <- stats::aov(log10(Ratio) ~ factor(Category, levels = categories), data = .SD)
    as.numeric(summary(fit)[[1]][1L, "Pr(>F)"])
  }
), by = Denominator]
setkey(expected_global_p, Denominator)
observed_global_p <- copy(global_significance[, .(Denominator, PValue)])
setkey(observed_global_p, Denominator)
stop_if(all(abs(expected_global_p$PValue - observed_global_p$PValue) < 1e-12),
  "The stored MKI67 ratio omnibus p values do not match the input values.")

message("PASS: MKI67/ACTB, MKI67/TUBB and MKI67/H3C1 candidate inputs are source-audited.")
