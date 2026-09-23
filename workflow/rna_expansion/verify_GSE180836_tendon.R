#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
out <- args[[1L]]
read <- function(name) fread(file.path(out, name))
meta <- read("sample_eligibility_31.csv")
qc <- read("sample_qc_31.csv")
samples <- read("tendon_sample_log2tpm_31.tsv.gz")
all_profile <- read("tendon_material_profile.csv.gz")
deg_profile <- read("tendon_degenerative_profile.csv.gz")
trauma_profile <- read("tendon_traumatic_profile.csv.gz")
comparison <- read("old_new_panel_comparison.csv")
summary <- read("processing_summary.csv")
hashes <- read("input_md5.csv")
stopifnot(nrow(meta) == 31L, uniqueN(meta$Subject) == 31L,
          all(meta$Characteristic02 == "Tendon"),
          all(meta$LibraryStrategy == "RNA-Seq"),
          all(meta$Condition %in% c("Degenerative", "Traumatic")),
          sum(meta$Condition == "Degenerative") == 15L,
          sum(meta$Condition == "Traumatic") == 16L,
          identical(meta$GSM, qc$GSM),
          nrow(samples) == nrow(all_profile),
          ncol(samples) == 32L,
          setequal(samples$Ensembl, all_profile$Ensembl),
          setequal(samples$Ensembl, deg_profile$Ensembl),
          setequal(samples$Ensembl, trauma_profile$Ensembl),
          !anyDuplicated(samples$Ensembl),
          !anyDuplicated(all_profile$Ensembl),
          all(qc$RawLibrarySize > 1e6),
          all(qc$GenesDetected > 5000L),
          all(abs(qc$TPMTotal - 1e6) < 1e-5),
          nrow(comparison) == 128L,
          summary$MatchedOldPanelGenes == 128L,
          all(file.exists(hashes$File)),
          identical(unname(tools::md5sum(hashes$File)), hashes$MD5))
mat <- as.matrix(samples[, -1L])
stopifnot(all(abs(colSums(2^mat - 0.5) - 1e6) < 1e-4))
setkey(all_profile, Ensembl)
setkey(deg_profile, Ensembl)
setkey(trauma_profile, Ensembl)
setkey(samples, Ensembl)
for (gene in c("ENSG00000134333", "ENSG00000111716",
               "ENSG00000166796", "ENSG00000171989")) {
  row <- as.numeric(samples[gene, -1L, with = FALSE])
  stopifnot(length(row) == 31L,
            abs(all_profile[gene, Log2TPM] - median(row)) < 1e-10,
            abs(deg_profile[gene, Log2TPM] -
                  median(row[meta$Condition == "Degenerative"])) < 1e-10,
            abs(trauma_profile[gene, Log2TPM] -
                  median(row[meta$Condition == "Traumatic"])) < 1e-10)
}
stopifnot(abs(summary$SpearmanWithOldProfile -
                cor(comparison$OldLog2TPM, comparison$NewLog2TPM,
                    method = "spearman")) < 1e-10)
cat("PASS: 31 independent torn-tendon libraries, stable Ensembl IDs,",
    "TPM totals, subgroup profiles, all 128 panel genes, and input fingerprints.\n")
