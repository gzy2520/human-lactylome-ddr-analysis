#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) normalizePath(args[[1]]) else normalizePath(".")
config_path <- file.path(
  project_root,
  "reanalysis",
  "config",
  "healthy_tissue_reference_files.csv"
)
table_dir <- file.path(project_root, "reanalysis", "results", "tables")
config <- read.csv(config_path, stringsAsFactors = FALSE, check.names = FALSE)

sha256_file <- function(path) {
  result <- system2("shasum", c("-a", "256", shQuote(path)), stdout = TRUE)
  sub("[[:space:]].*$", "", result[[1]])
}

base_accession <- function(values) {
  values <- sub("^.*\\|([^|]+)\\|.*$", "\\1", values)
  sub("-[0-9]+$", "", values)
}

count_proteins <- function(path) {
  data <- read.delim(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = ""
  )
  keep <- rep(TRUE, nrow(data))
  if ("Reverse" %in% names(data)) {
    keep <- keep & data$Reverse != "+"
  }
  if ("Potential contaminant" %in% names(data)) {
    keep <- keep & data$`Potential contaminant` != "+"
  }
  if ("Only identified by site" %in% names(data)) {
    keep <- keep & data$`Only identified by site` != "+"
  }
  column <- intersect(c("Protein IDs", "Majority protein IDs"), names(data))[[1]]
  accessions <- unlist(strsplit(as.character(data[[column]][keep]), ";", fixed = TRUE))
  accessions <- trimws(accessions)
  accessions <- accessions[nzchar(accessions) & !is.na(accessions)]
  length(unique(base_accession(accessions)))
}

download_one <- function(i) {
  row <- config[i,]
  archive <- file.path(
    project_root,
    "data",
    row$PXD,
    "search_results",
    row$FileName
  )
  dir.create(dirname(archive), recursive = TRUE, showWarnings = FALSE)
  complete <- file.exists(archive) && file.info(archive)$size == row$SizeBytes
  if (!complete) {
    message(sprintf("[%d/%d] downloading %s", i, nrow(config), row$FileName))
    command <- paste(
      "curl --fail --location --continue-at - --retry 4 --retry-delay 2 --silent --show-error",
      "--output", shQuote(archive), shQuote(row$SourceURL)
    )
    status <- system(command)
    if (status != 0) {
      return(data.frame(
        TissueKey = row$TissueKey,
        DisplayName = row$DisplayName,
        PXD = row$PXD,
        FileName = row$FileName,
        SizeBytes = if (file.exists(archive)) file.info(archive)$size else 0,
        SHA256 = "",
        ProteinGroupsPath = "",
        ProteinCount = NA_integer_,
        Status = paste0("download_failed_exit_", status),
        SourceURL = row$SourceURL,
        stringsAsFactors = FALSE
      ))
    }
  }

  extract_dir <- file.path(
    project_root,
    "data",
    row$PXD,
    "search_results",
    "extracted_healthy_tissues",
    row$TissueKey
  )
  dir.create(extract_dir, recursive = TRUE, showWarnings = FALSE)
  existing <- list.files(
    extract_dir,
    pattern = "proteinGroups\\.txt$",
    recursive = TRUE,
    full.names = TRUE
  )
  if (!length(existing)) {
    status <- system2(
      "unzip",
      c("-j", "-n", shQuote(archive), "*proteinGroups.txt", "-d", shQuote(extract_dir))
    )
    if (status != 0) {
      return(data.frame(
        TissueKey = row$TissueKey,
        DisplayName = row$DisplayName,
        PXD = row$PXD,
        FileName = row$FileName,
        SizeBytes = file.info(archive)$size,
        SHA256 = sha256_file(archive),
        ProteinGroupsPath = "",
        ProteinCount = NA_integer_,
        Status = paste0("proteinGroups_extract_failed_exit_", status),
        SourceURL = row$SourceURL,
        stringsAsFactors = FALSE
      ))
    }
    existing <- list.files(
      extract_dir,
      pattern = "proteinGroups\\.txt$",
      recursive = TRUE,
      full.names = TRUE
    )
  }
  protein_path <- existing[[1]]
  data.frame(
    TissueKey = row$TissueKey,
    DisplayName = row$DisplayName,
    PXD = row$PXD,
    FileName = row$FileName,
    SizeBytes = file.info(archive)$size,
    SHA256 = sha256_file(archive),
    ProteinGroupsPath = substring(protein_path, nchar(project_root) + 2),
    ProteinCount = count_proteins(protein_path),
    Status = "downloaded_extracted_counted",
    SourceURL = row$SourceURL,
    stringsAsFactors = FALSE
  )
}

worker_count <- min(3L, nrow(config))
records <- parallel::mclapply(seq_len(nrow(config)), download_one, mc.cores = worker_count)
result <- bind_rows(records) |>
  arrange(TissueKey)
write.csv(
  result,
  file.path(table_dir, "healthy_tissue_reference_acquisition_manifest.csv"),
  row.names = FALSE,
  na = ""
)
message("Healthy tissue reference acquisition complete.")
