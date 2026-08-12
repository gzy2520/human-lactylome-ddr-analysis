# accession_utils.R - Shared accession handling utilities for kla reanalysis
lib_loaded <- TRUE

base_accession <- function(values) {
  values <- trimws(as.character(values))
  values <- sub("^.*\\|([^|]+)\\|.*$", "\\1", values)
  values <- sub("^NX_", "", values)
  values <- sub("-[0-9]+$", "", values)
  values
}

is_uniprot <- function(values) {
  grepl(
    paste0(
      "^(?:",
      "[OPQ][0-9][A-Z0-9]{3}[0-9]|",
      "[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9](?:[A-Z0-9]{3}[0-9])?",
      ")$"
    ),
    values
  )
}

split_accessions <- function(values) {
  values <- unlist(strsplit(as.character(values), "[;,]"))
  values <- base_accession(values)
  sort(unique(values[is_uniprot(values)]))
}

split_protein_identifiers <- function(values) {
  values <- unlist(strsplit(as.character(values), "[;,]"))
  values <- trimws(values)
  uniprot <- base_accession(values)
  uniprot <- uniprot[is_uniprot(uniprot)]
  ensembl <- sub("_RNA$", "", values)
  ensembl <- sub("\\.[0-9]+$", "", ensembl)
  ensembl <- ensembl[grepl("^ENSP[0-9]+$", ensembl)]
  sort(unique(c(uniprot, ensembl)))
}

identifier_type <- function(values) {
  ifelse(grepl("^ENSP[0-9]+$", values), "ENSEMBLPROT", "UniProtKB")
}

# match_target_accession: from heatmap scripts (analyze_kla_regulator_intensity.R 91-99,
#   analyze_kla_regulator_whole_proteome_intensity.R 98-105).
#   Adapted for lib: target_accessions passed as explicit parameter instead of closure.
#   Behavior: for each input, split by [;, ]+, base_accession, return first hit
#   matching target_accessions, or NA_character_.
match_target_accession <- function(values, target_accessions) {
  vapply(as.character(values), function(value) {
    if (is.na(value) || !nzchar(value)) return(NA_character_)
    tokens <- unique(trimws(unlist(strsplit(value, "[;, ]+"))))
    tokens <- base_accession(tokens)
    hits <- tokens[tokens %in% target_accessions]
    if (length(hits)) hits[[1]] else NA_character_
  }, character(1))
}

# accession_feature: from heatmap scripts (analyze_kla_regulator_intensity.R 101-108,
#   analyze_kla_regulator_whole_proteome_intensity.R 108-115).
#   Behavior: split by [;, ]+, base_accession, return "ACC:<first_token>"
#   or NA_character_. Does NOT depend on target_accessions.
accession_feature <- function(values) {
  vapply(as.character(values), function(value) {
    if (is.na(value) || !nzchar(value)) return(NA_character_)
    tokens <- unique(trimws(unlist(strsplit(value, "[;, ]+"))))
    tokens <- base_accession(tokens)
    tokens <- tokens[nzchar(tokens)]
    if (length(tokens)) paste0("ACC:", tokens[[1]]) else NA_character_
  }, character(1))
}

safe_numeric <- function(values) {
  suppressWarnings(as.numeric(gsub(",", "", as.character(values), fixed = TRUE)))
}
