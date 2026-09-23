#!/usr/bin/env Rscript
# Untreated human torn shoulder tendon RNA from GSE180836 for KLA31_01.
# The author gene-count archive is read by stable Ensembl gene ID only.
suppressPackageStartupMessages(library(data.table))
set.seed(25)
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
out <- args[[1L]]
if (dir.exists(out) && length(list.files(out, all.files = TRUE, no.. = TRUE)))
  stop("Refusing to overwrite nonempty output")
root <- "data/candidate/rna_expansion"
tar_file <- file.path(root, "GSE180836/GSE180836_RAW.tar")
series_file <- file.path(root, "GSE180836/GSE180836_series_matrix.txt.gz")
length_file <- file.path(root, "annotation/human_gene_lengths_ensembl111.tsv")
inputs <- c(tar_file, series_file, length_file)
stopifnot(all(file.exists(inputs)))
source("workflow/lib_kla31_expression_20260916.R")

lines <- readLines(gzfile(series_file), warn = FALSE)
field <- function(prefix) {
  row <- lines[startsWith(lines, prefix)]
  stopifnot(length(row) == 1L)
  gsub('^"|"$', "", strsplit(row, "\t", fixed = TRUE)[[1L]][-1L])
}
meta <- data.table(GSM = field("!Sample_geo_accession"),
                   Title = field("!Sample_title"),
                   LibraryStrategy = field("!Sample_library_strategy"),
                   Organism = field("!Sample_organism_ch1"),
                   SourceName = field("!Sample_source_name_ch1"))
characteristics <- lines[startsWith(lines, "!Sample_characteristics_ch1")]
stopifnot(length(characteristics) == 5L, nrow(meta) == 31L,
          !anyDuplicated(meta$GSM))
for (i in seq_along(characteristics)) {
  value <- gsub('^"|"$', "",
                strsplit(characteristics[[i]], "\t", fixed = TRUE)[[1L]][-1L])
  stopifnot(length(value) == nrow(meta))
  meta[, (sprintf("Characteristic%02d", i)) := sub("^[^:]+: *", "", value)]
}
meta[, Condition := Characteristic05]
meta[, Subject := sub("^Tendon (RCT-[0-9]+) .*", "\\1", Title)]
stopifnot(all(meta$Characteristic02 == "Tendon"),
          all(meta$Organism == "Homo sapiens"),
          all(meta$SourceName == "Tendon"),
          all(meta$LibraryStrategy == "RNA-Seq"),
          all(meta$Condition %in% c("Traumatic", "Degenerative")),
          uniqueN(meta$Subject) == 31L,
          all(grepl("^RCT-[0-9]+$", meta$Subject)))

archive_names <- utils::untar(tar_file, list = TRUE)
stopifnot(length(archive_names) == 31L,
          all(grepl("gene_counts.txt.gz$", archive_names)),
          setequal(sub("_.*$", "", archive_names), meta$GSM))
tmp <- tempfile("GSE180836_extract_")
dir.create(tmp)
on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
utils::untar(tar_file, exdir = tmp)
files <- file.path(tmp, archive_names)
stopifnot(all(file.exists(files)))
names(files) <- sub("_.*$", "", archive_names)

read_one <- function(path) {
  tab <- fread(path, select = c("Feature", "Count"))
  stopifnot(nrow(tab) > 25000L, all(grepl("^ENSG[0-9]+", tab$Feature)),
            all(is.finite(tab$Count)), all(tab$Count >= 0),
            all(tab$Count == floor(tab$Count)))
  tab[, Ensembl := strip_ensembl_version(Feature)]
  tab[, .(Count = sum(Count)), by = Ensembl]
}
sample_tables <- lapply(files[meta$GSM], read_one)
names(sample_tables) <- meta$GSM
genes <- sort(sample_tables[[1L]]$Ensembl)
stopifnot(all(vapply(sample_tables, function(tab)
  setequal(tab$Ensembl, genes), logical(1))))
counts <- vapply(sample_tables, function(tab)
  tab$Count[match(genes, tab$Ensembl)], numeric(length(genes)))
rownames(counts) <- genes
stopifnot(!anyNA(counts), ncol(counts) == 31L)

lengths_tab <- fread(length_file, header = FALSE,
                     col.names = c("Ensembl", "LengthBp"))
lengths <- setNames(lengths_tab$LengthBp, lengths_tab$Ensembl)
tpm <- tpm_from_counts(counts, lengths)
logtpm <- log2_tpm(tpm)
stopifnot(nrow(logtpm) > 15000L, ncol(logtpm) == 31L,
          all(abs(colSums(tpm) - 1e6) < 1e-5))
qc <- data.table(GSM = meta$GSM, Subject = meta$Subject,
                 Condition = meta$Condition,
                 RawLibrarySize = colSums(counts),
                 GenesDetected = colSums(counts > 0),
                 TPMTotal = colSums(tpm))
stopifnot(all(qc$RawLibrarySize > 1e6), all(qc$GenesDetected > 5000L))

profile_for <- function(sel, key) {
  value <- apply(logtpm[, sel, drop = FALSE], 1L, median)
  data.table(Ensembl = names(value), Log2TPM = as.numeric(value),
             ReferenceKey = key)
}
profile_all <- profile_for(seq_len(nrow(meta)), "GSE180836_torn_tendon")
profile_deg <- profile_for(which(meta$Condition == "Degenerative"),
                           "GSE180836_degenerative_tendon")
profile_trauma <- profile_for(which(meta$Condition == "Traumatic"),
                              "GSE180836_traumatic_tendon")
old_file <- "outputs/20260920_lactate_metabolism/tables/material_gene_expression.csv"
old <- fread(old_file)[ReferenceKey == "tendon_nondiabetic",
                         .(Ensembl, OldLog2TPM = UnsmoothedLog2TPM)]
comparison <- merge(old, profile_all[, .(Ensembl, NewLog2TPM = Log2TPM)],
                    by = "Ensembl")
stopifnot(nrow(old) == 128L, nrow(comparison) > 110L)
summary <- data.table(
  Source = "GSE180836", GroupID = "KLA31_01",
  SelectedSamples = nrow(meta), UniqueSubjects = uniqueN(meta$Subject),
  Traumatic = sum(meta$Condition == "Traumatic"),
  Degenerative = sum(meta$Condition == "Degenerative"),
  SourceEnsemblRows = nrow(counts), TPMGeneRows = nrow(tpm),
  OldPanelGenes = nrow(old), MatchedOldPanelGenes = nrow(comparison),
  SpearmanWithOldProfile = cor(comparison$OldLog2TPM,
                               comparison$NewLog2TPM, method = "spearman"),
  MinimumRawLibrarySize = min(qc$RawLibrarySize),
  MinimumGenesDetected = min(qc$GenesDetected)
)

dir.create(out, recursive = TRUE, showWarnings = FALSE)
fwrite(meta, file.path(out, "sample_eligibility_31.csv"))
fwrite(qc, file.path(out, "sample_qc_31.csv"))
fwrite(data.table(Ensembl = rownames(logtpm), as.data.table(logtpm)),
       file.path(out, "tendon_sample_log2tpm_31.tsv.gz"), sep = "\t")
fwrite(profile_all, file.path(out, "tendon_material_profile.csv.gz"))
fwrite(profile_deg, file.path(out, "tendon_degenerative_profile.csv.gz"))
fwrite(profile_trauma, file.path(out, "tendon_traumatic_profile.csv.gz"))
fwrite(comparison, file.path(out, "old_new_panel_comparison.csv"))
fwrite(summary, file.path(out, "processing_summary.csv"))
fwrite(data.table(File = c(inputs, old_file),
                  MD5 = unname(tools::md5sum(c(inputs, old_file)))),
       file.path(out, "input_md5.csv"))
writeLines(trimws(capture.output(sessionInfo()), which = "right"),
           file.path(out, "sessionInfo.txt"))
print(summary)
