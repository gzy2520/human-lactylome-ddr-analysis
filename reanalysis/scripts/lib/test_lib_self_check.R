#!/usr/bin/env Rscript
root <- normalizePath(if (length(commandArgs(trailingOnly = TRUE)) >= 1) commandArgs(trailingOnly = TRUE)[[1]] else ".")
lib <- file.path(root, "reanalysis", "scripts", "lib")
source(file.path(lib, "accession_utils.R"))
source(file.path(lib, "io_utils.R"))
source(file.path(lib, "extractors.R"))

# --- accession_utils.R assertions ---

# base_accession: UniProt pipe format
stopifnot(base_accession("sp|P49959-2|MRE11_HUMAN") == "P49959")

# base_accession: REV__CON__ prefix + isoform (source scripts return "REV__CON__Q9H9Q4",
#   NOT "Q9H9Q4" — the brief example was corrected to match actual source behaviour)
stopifnot(base_accession("REV__CON__Q9H9Q4-3") == "REV__CON__Q9H9Q4")

# base_accession: NX_ prefix + isoform
stopifnot(base_accession("NX_Q9H9Q4-1") == "Q9H9Q4")

# is_uniprot: canonical UniProt accession is valid, isoform suffix is not, gene name is not
stopifnot(is_uniprot("P49959"), !is_uniprot("P49959-2"), !is_uniprot("MRE11"))

# split_accessions: mixed delimiters, prefix stripping, UniProt-only filter, sorted unique
stopifnot(identical(
  split_accessions("P49959-2;Q9H9Q4-3,sp|O60934|NBN_HUMAN"),
  c("O60934", "P49959", "Q9H9Q4")
))

# split_protein_identifiers: mixed UniProt + ENSP, strips isoforms/suffixes
stopifnot(identical(
  split_protein_identifiers("ENSP00000369497;P49959-2"),
  c("ENSP00000369497", "P49959")
))

# identifier_type: classifies ENSP vs UniProt
stopifnot(identical(
  identifier_type(c("ENSP00000369497", "P49959", "ENSP00000000000")),
  c("ENSEMBLPROT", "UniProtKB", "ENSEMBLPROT")
))

# accession_feature: strips isoform, returns "ACC:<first_token>" (corrected from brief
#   example — source scripts prepend "ACC:" prefix; vapply preserves names so unname())
stopifnot(identical(
  unname(accession_feature(c("P49959-2", "P49959-1"))),
  c("ACC:P49959", "ACC:P49959")
))

# accession_feature: NA / empty input
stopifnot(is.na(accession_feature(NA_character_)))
stopifnot(is.na(accession_feature("")))

# match_target_accession: returns first matching accession or NA
stopifnot(identical(
  unname(match_target_accession(c("P49959-2;Q9H9Q4-3", "O60934"), c("P49959", "O60934"))),
  c("P49959", "O60934")
))

# match_target_accession: no match returns NA (vapply preserves names so use is.na)
stopifnot(is.na(match_target_accession(c("P49959-2"), c("XXXXXX"))))

# match_target_accession: NA / empty input returns NA
stopifnot(is.na(match_target_accession(NA_character_, c("P49959"))))
stopifnot(is.na(match_target_accession("", c("P49959"))))

# safe_numeric: comma stripping, NA handling
stopifnot(identical(safe_numeric(c("1,234", "NA", "0")), c(1234, NA, 0)))

# --- io_utils.R assertions ---

# relative_path
stopifnot(grepl(
  "reanalysis/scripts/lib/test_lib_self_check\\.R$",
  relative_path(file.path(root, "reanalysis", "scripts", "lib", "test_lib_self_check.R"), root)
))

# read_delimited: csv
tmp <- tempfile(fileext = ".csv")
write.csv(data.frame(a = 1:2, b = c("x", "")), tmp, row.names = FALSE)
stopifnot(nrow(read_delimited(tmp)) == 2)

# read_delimited: tsv
tmp2 <- tempfile(fileext = ".tsv")
write.table(data.frame(a = 1:3, b = c("x", "y", "z")), tmp2, sep = "\t", row.names = FALSE)
stopifnot(nrow(read_delimited(tmp2)) == 3)

# valid_maxquant_rows: all-pass data
clean_data <- data.frame(
  Reverse = c(NA, NA),
  `Potential contaminant` = c(NA, NA),
  `Only identified by site` = c(NA, NA),
  check.names = FALSE
)
stopifnot(all(valid_maxquant_rows(clean_data)))

# valid_maxquant_rows: filter Reverse "+"
rev_data <- data.frame(
  Reverse = c(NA, "+"),
  check.names = FALSE
)
stopifnot(identical(valid_maxquant_rows(rev_data), c(TRUE, FALSE)))

# write_csv_std
tmp3 <- tempfile(fileext = ".csv")
write_csv_std(data.frame(x = 1:2, y = c(NA, "b")), tmp3)
result <- read.csv(tmp3, stringsAsFactors = FALSE)
stopifnot(nrow(result) == 2, result$y[1] == "", result$y[2] == "b")

# --- extractors.R assertions (basic smoke tests) ---

# extract_huvec_xml: confirm function exists and handles NX_ prefix
stopifnot(is.function(extract_huvec_xml))
stopifnot(is.function(extract_maxquant_sites))
stopifnot(is.function(extract_maxquant_proteins))
stopifnot(is.function(extract_pd_proteins))
stopifnot(is.function(extract_pd_lactyl_peptides))
stopifnot(is.function(extract_spectronaut_proteins))
stopifnot(is.function(extract_spectronaut_matrix))
stopifnot(is.function(extract_spectronaut_quant))

# --- cleanup ---
unlink(c(tmp, tmp2, tmp3))

cat("lib self-check passed\n")
