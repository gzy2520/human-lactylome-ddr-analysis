#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
out <- args[[1L]]
read <- function(name) fread(file.path(out, name))
eligibility <- read("sample_eligibility_108.csv")
qc <- read("sample_qc_28.csv")
samples <- read("HTS_sample_log2tpm_28.tsv.gz")
donors <- read("HTS_donor_log2tpm_26.tsv.gz")
profile <- read("HTS_material_profile.csv.gz")
comparison <- read("old_new_panel_comparison.csv")
summary <- read("processing_summary.csv")
hashes <- read("input_md5.csv")

stopifnot(nrow(eligibility) == 108L,
          sum(eligibility$Eligibility == "include_HTS") == 28L,
          all(eligibility[Eligibility == "include_HTS", Characteristic03] == "HTS"),
          all(eligibility[Eligibility == "include_HTS", Characteristic10] == "scar"),
          !any(eligibility[Eligibility == "include_HTS", Characteristic03] ==
                 "Normal scar"),
          nrow(qc) == 28L, uniqueN(qc$Subject) == 26L,
          ncol(samples) == 29L, ncol(donors) == 27L,
          nrow(samples) == nrow(donors),
          nrow(profile) == nrow(samples),
          !anyDuplicated(profile$Ensembl),
          setequal(profile$Ensembl, samples$Ensembl),
          nrow(comparison) == 128L,
          all(qc$RawLibrarySize > 1e6),
          all(qc$GenesDetected > 5000L),
          all(abs(qc$TPMTotal - 1e6) < 1e-5),
          all(file.exists(hashes$File)),
          identical(unname(tools::md5sum(hashes$File)), hashes$MD5))

# Reconstruct normalized TPM from the stored log scale for every sample, and
# verify the two-stage donor-balanced material median on selected genes.
sample_expr <- as.matrix(samples[, -1L])
stopifnot(all(abs(colSums(2^sample_expr - 0.5) - 1e6) < 1e-4))
setkey(profile, Ensembl)
setkey(donors, Ensembl)
for (gene in c("ENSG00000134333", "ENSG00000111716",
               "ENSG00000166796", "ENSG00000171989")) {
  value <- as.numeric(donors[gene, -1L, with = FALSE])
  stopifnot(length(value) == 26L,
            abs(profile[gene, Log2TPM] - median(value)) < 1e-10)
}
stopifnot(summary$MatchedOldPanelGenes == 128L,
          abs(summary$SpearmanWithOldProfile -
                cor(comparison$OldLog2TPM, comparison$NewLog2TPM,
                    method = "spearman")) < 1e-10)
cat("PASS: 28 true HTS libraries, 26 donors, stable IDs, TPM totals,",
    "donor-balanced profile, all 128 panel genes, and input fingerprints.\n")
