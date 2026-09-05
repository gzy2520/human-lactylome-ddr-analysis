#!/usr/bin/env Rscript

# Validate the 10 regulator sample-level whole-proteome expression percentiles contract.

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args)) normalizePath(args[[1L]], mustWork = TRUE) else normalizePath(".", mustWork = TRUE)
stop_if <- function(condition, message) if (!isTRUE(condition)) stop(message, call. = FALSE)

expected_genes <- c("AARS1", "ACAT2", "KRT18", "SIRT2", "PARK7", "HDAC1", "HDAC2", "BRD4", "SMARCA4", "TRIM33")
expected_accessions <- c(
  AARS1 = "P49588", ACAT2 = "Q9BWD1", KRT18 = "P05783", SIRT2 = "Q8IXJ6",
  PARK7 = "Q99497", HDAC1 = "Q13547", HDAC2 = "Q92769", BRD4 = "O60885",
  SMARCA4 = "P51532", TRIM33 = "Q9UPN9"
)

csv_path <- file.path(
  project_root,
  "data/candidate/escc_inclusion_20260903_pxd065830_tumor_reference/target_10_regulators_sample_percentiles.csv"
)
stop_if(file.exists(csv_path), paste0("Missing target regulators sample percentiles CSV: ", csv_path))

dt <- fread(csv_path)
required_cols <- c("PXD", "SampleGroup", "Category", "QuantSample", "GeneSymbol", "BaseAccession", "Role", "WholeProteomePercentile", "Detected", "Scope")
stop_if(all(required_cols %in% names(dt)), "Missing required columns in sample percentiles table.")
stop_if(setequal(unique(dt$GeneSymbol), expected_genes), "Target gene set mismatch.")
stop_if(setequal(unique(dt$BaseAccession), unname(expected_accessions)), "Target BaseAccession set mismatch.")
stop_if(all(dt$WholeProteomePercentile >= 0 & dt$WholeProteomePercentile <= 100), "Percentiles must be within [0, 100].")

# Verify sample counts per gene:
# Candidate scope must contain 634 unique sample observations per gene
counts_candidate <- dt[, uniqueN(SampleKey), by = GeneSymbol]
stop_if(all(counts_candidate$V1 == 634L), "Candidate scope must contain exactly 634 sample observations per gene.")

# Baseline scope must contain 540 unique sample observations per gene
counts_baseline <- dt[Scope == "baseline_30_datasets", uniqueN(SampleKey), by = GeneSymbol]
stop_if(all(counts_baseline$V1 == 540L), "Baseline scope must contain exactly 540 sample observations per gene.")

# Verify image outputs in candidate scope
out_dir_candidate <- file.path(
  project_root,
  "results/candidate/escc_inclusion_20260903_pxd065830_tumor_reference/regulator_sample_percentiles"
)
stop_if(dir.exists(out_dir_candidate), paste0("Candidate output directory missing: ", out_dir_candidate))

for (gene in expected_genes) {
  acc <- expected_accessions[[gene]]
  png_path <- file.path(out_dir_candidate, paste0("Figure_3_regulator_sample_percentile_", gene, "_", acc, ".png"))
  pdf_path <- file.path(out_dir_candidate, paste0("Figure_3_regulator_sample_percentile_", gene, "_", acc, ".pdf"))
  stop_if(file.exists(png_path) && file.info(png_path)$size > 10000L, paste("Missing or invalid PNG for", gene))
  stop_if(file.exists(pdf_path) && file.info(pdf_path)$size > 10000L, paste("Missing or invalid PDF for", gene))
}

# Verify image outputs in baseline scope
out_dir_baseline <- file.path(out_dir_candidate, "baseline_30_datasets")
stop_if(dir.exists(out_dir_baseline), paste0("Baseline output directory missing: ", out_dir_baseline))

for (gene in expected_genes) {
  acc <- expected_accessions[[gene]]
  png_path <- file.path(out_dir_baseline, paste0("Figure_3_regulator_sample_percentile_", gene, "_", acc, "_30datasets.png"))
  pdf_path <- file.path(out_dir_baseline, paste0("Figure_3_regulator_sample_percentile_", gene, "_", acc, "_30datasets.pdf"))
  stop_if(file.exists(png_path) && file.info(png_path)$size > 10000L, paste("Missing or invalid baseline PNG for", gene))
  stop_if(file.exists(pdf_path) && file.info(pdf_path)$size > 10000L, paste("Missing or invalid baseline PDF for", gene))
}

message("PASS: All 10 regulator sample-level scatter plots and percentiles validated successfully.")
