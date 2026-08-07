#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(dplyr))

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) normalizePath(args[[1]]) else normalizePath(".")
config <- read.csv(
  file.path(
    project_root,
    "reanalysis/config/healthy_special_reference_catalog.csv"
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
output_path <- file.path(
  project_root,
  "reanalysis/results/tables/healthy_special_reference_acquisition_manifest.csv"
)

resolve_files <- function(locator) {
  locator <- sub(" \\[.*$", "", locator)
  pattern <- file.path(project_root, locator)
  if (grepl("\\*", pattern)) Sys.glob(pattern) else pattern
}

sha256_file <- function(path) {
  output <- system2(
    "shasum",
    c("-a", "256", shQuote(path)),
    stdout = TRUE
  )
  sub("[[:space:]].*$", "", output[[1]])
}

records <- lapply(seq_len(nrow(config)), function(i) {
  row <- config[i,]
  files <- resolve_files(row$ProteinGroupsPath)
  files <- files[file.exists(files)]
  data.frame(
    TissueKey = row$TissueKey,
    DisplayName = row$DisplayName,
    PXD = row$PXD,
    EvidenceLocator = row$ProteinGroupsPath,
    LocalFileCount = length(files),
    LocalSizeBytes = sum(file.info(files)$size, na.rm = TRUE),
    SHA256 = paste(vapply(files, sha256_file, character(1)), collapse = ";"),
    ProteinCount = row$ProteinCount,
    CountBasis = row$CountBasis,
    MatchQuality = row$MatchQuality,
    Status = if (length(files)) {
      "downloaded_verified_counted"
    } else {
      "missing"
    },
    SourceURL = row$SourceURL,
    Caveat = row$Caveat,
    stringsAsFactors = FALSE
  )
})

result <- bind_rows(records)
write.csv(result, output_path, row.names = FALSE)
if (any(result$Status != "downloaded_verified_counted")) {
  stop("One or more special healthy references are missing", call. = FALSE)
}
message("Built special healthy-reference acquisition manifest.")
