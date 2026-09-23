#!/usr/bin/env Rscript
# Independent, untreated hypertrophic-scar RNA profile for KLA31_03.
# NCBI Entrez gene IDs are mapped to Ensembl IDs; gene symbols are never keys.
suppressPackageStartupMessages(library(data.table))
set.seed(25)
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
out <- args[[1L]]
if (dir.exists(out) && length(list.files(out, all.files = TRUE, no.. = TRUE)))
  stop("Refusing to overwrite nonempty output")

root <- "data/candidate/rna_expansion"
counts_file <- file.path(root, "GSE178411/GSE178411_counts.txt.gz")
series_file <- file.path(root, "GSE178411/GSE178411_series_matrix.txt.gz")
map_file <- file.path(root, "annotation/human_gene2ensembl.tsv")
length_file <- file.path(root, "annotation/human_gene_lengths_ensembl111.tsv")
inputs <- c(counts_file, series_file, map_file, length_file)
stopifnot(all(file.exists(inputs)))
source("workflow/lib_kla31_expression_20260916.R")

# Parse each GEO sample field in its original order. The series uses two fields
# called "wound type"; field 3 is the precise HTS/normal-scar diagnosis, field
# 10 is only a broad scar/wound/uninjured classification.
lines <- readLines(gzfile(series_file), warn = FALSE)
field <- function(prefix) {
  row <- lines[startsWith(lines, prefix)]
  stopifnot(length(row) == 1L)
  gsub('^"|"$', "", strsplit(row, "\t", fixed = TRUE)[[1L]][-1L])
}
meta <- data.table(GSM = field("!Sample_geo_accession"),
                   Title = field("!Sample_title"),
                   LibraryStrategy = field("!Sample_library_strategy"),
                   Organism = field("!Sample_organism_ch1"))
characteristics <- lines[startsWith(lines, "!Sample_characteristics_ch1")]
stopifnot(length(characteristics) == 11L, nrow(meta) == 108L,
          !anyDuplicated(meta$GSM))
for (i in seq_along(characteristics)) {
  value <- gsub('^"|"$', "",
                strsplit(characteristics[[i]], "\t", fixed = TRUE)[[1L]][-1L])
  stopifnot(length(value) == nrow(meta))
  meta[, (sprintf("Characteristic%02d", i)) := sub("^[^:]+: *", "", value)]
}
meta[, SourceColumn := sub(":.*$", "", Title)]
meta[, Eligibility := fifelse(Characteristic03 == "HTS" &
                               Characteristic10 == "scar" &
                               LibraryStrategy == "RNA-Seq",
                             "include_HTS", "exclude_other_material")]
selected <- meta[Eligibility == "include_HTS"]
stopifnot(nrow(selected) == 28L, uniqueN(selected$Characteristic01) == 26L,
          !anyDuplicated(meta$SourceColumn),
          all(selected$Organism == "Homo sapiens"),
          all(selected$Characteristic02 == "Skin"),
          all(grepl(" HTS$", selected$Title)))

tab <- suppressWarnings(fread(counts_file))
setnames(tab, 1L, "Entrez") # Author file has 108 names for 109 columns.
stopifnot(nrow(tab) == 28395L, ncol(tab) == 109L,
          !anyDuplicated(tab$Entrez), all(grepl("^[0-9]+$", tab$Entrez)),
          setequal(names(tab)[-1L], meta$SourceColumn),
          all(selected$SourceColumn %in% names(tab)))
cols <- selected$SourceColumn
raw <- as.matrix(tab[, ..cols])
rownames(raw) <- as.character(tab$Entrez)
storage.mode(raw) <- "double"
stopifnot(all(is.finite(raw)), all(raw >= 0), all(raw == floor(raw)))

mapping <- load_id_mapping(map_file)
entrez <- build_entrez_map(mapping)
mapped <- collapse_by_stable_id(raw, entrez$map$GeneID,
                                entrez$map$Ensembl_gene_identifier, "counts")
counts <- mapped$matrix
lengths_tab <- fread(length_file, header = FALSE,
                     col.names = c("Ensembl", "LengthBp"))
stopifnot(!anyDuplicated(lengths_tab$Ensembl),
          all(is.finite(lengths_tab$LengthBp) & lengths_tab$LengthBp > 0))
lengths <- setNames(lengths_tab$LengthBp, lengths_tab$Ensembl)
tpm <- tpm_from_counts(counts, lengths)
logtpm <- log2_tpm(tpm)
stopifnot(nrow(logtpm) > 15000L, ncol(logtpm) == 28L,
          all(abs(colSums(tpm) - 1e6) < 1e-5),
          all(is.finite(logtpm)))

sample_qc <- data.table(
  GSM = selected$GSM, SourceColumn = selected$SourceColumn,
  Subject = selected$Characteristic01,
  RawLibrarySize = colSums(raw), MappedLibrarySize = colSums(counts),
  GenesDetected = colSums(counts > 0), TPMTotal = colSums(tpm),
  ZeroFraction = colMeans(counts == 0)
)
stopifnot(all(sample_qc$RawLibrarySize > 1e6),
          all(sample_qc$MappedLibrarySize > 1e6),
          all(sample_qc$GenesDetected > 5000L))

# Repeated sections from one patient do not receive extra weight in the
# material profile. First take each donor's gene median, then the donor median.
donors <- sort(unique(selected$Characteristic01))
donor_profile <- vapply(donors, function(subject) {
  ix <- which(selected$Characteristic01 == subject)
  if (length(ix) == 1L) logtpm[, ix]
  else apply(logtpm[, ix, drop = FALSE], 1L, median)
}, numeric(nrow(logtpm)))
rownames(donor_profile) <- rownames(logtpm)
colnames(donor_profile) <- donors
material_profile <- apply(donor_profile, 1L, median)
profile <- data.table(Ensembl = names(material_profile),
                      Log2TPM = as.numeric(material_profile),
                      ReferenceKey = "GSE178411_HTS")

old_file <- "outputs/20260920_lactate_metabolism/tables/material_gene_expression.csv"
old <- fread(old_file)[ReferenceKey == "GSE181540_HS_input",
                         .(Ensembl, OldLog2TPM = UnsmoothedLog2TPM)]
comparison <- merge(old, profile[, .(Ensembl, NewLog2TPM = Log2TPM)],
                    by = "Ensembl")
stopifnot(nrow(old) == 128L, nrow(comparison) > 110L)
summary <- data.table(
  Source = "GSE178411", GroupID = "KLA31_03",
  RNAReferenceKey = "GSE178411_HTS",
  SelectedSamples = nrow(selected), UniqueDonors = length(donors),
  SourceEntrezRows = nrow(tab), MappedEnsemblRows = nrow(counts),
  TPMGeneRows = nrow(tpm), OldPanelGenes = nrow(old),
  MatchedOldPanelGenes = nrow(comparison),
  SpearmanWithOldProfile = cor(comparison$OldLog2TPM,
                               comparison$NewLog2TPM, method = "spearman"),
  MedianRawLibrarySize = median(sample_qc$RawLibrarySize),
  MinimumRawLibrarySize = min(sample_qc$RawLibrarySize),
  MinimumGenesDetected = min(sample_qc$GenesDetected),
  AmbiguousEntrezDropped = entrez$dropped_ambiguous_entrez
)

dir.create(out, recursive = TRUE, showWarnings = FALSE)
fwrite(meta, file.path(out, "sample_eligibility_108.csv"))
fwrite(sample_qc, file.path(out, "sample_qc_28.csv"))
fwrite(data.table(Ensembl = rownames(logtpm), as.data.table(logtpm)),
       file.path(out, "HTS_sample_log2tpm_28.tsv.gz"), sep = "\t")
fwrite(data.table(Ensembl = rownames(donor_profile),
                  as.data.table(donor_profile)),
       file.path(out, "HTS_donor_log2tpm_26.tsv.gz"), sep = "\t")
fwrite(profile, file.path(out, "HTS_material_profile.csv.gz"))
fwrite(comparison, file.path(out, "old_new_panel_comparison.csv"))
fwrite(summary, file.path(out, "processing_summary.csv"))
fwrite(data.table(File = c(inputs, old_file),
                  MD5 = unname(tools::md5sum(c(inputs, old_file)))),
       file.path(out, "input_md5.csv"))
writeLines(trimws(capture.output(sessionInfo()), which = "right"),
           file.path(out, "sessionInfo.txt"))
print(summary)
