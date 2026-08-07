#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(dplyr))

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) normalizePath(args[[1]]) else normalizePath(".")
config_path <- file.path(
  project_root,
  "reanalysis/config/additional_lactylome_pair_files.csv"
)
table_dir <- file.path(project_root, "reanalysis/results/tables")
manifest_path <- file.path(table_dir, "lactylome_pair_download_manifest.csv")
remote_path <- file.path(table_dir, "lactylome_pair_remote_file_sizes.csv")

config <- read.csv(config_path, stringsAsFactors = FALSE, check.names = FALSE)
remote <- read.csv(remote_path, stringsAsFactors = FALSE, check.names = FALSE)
existing <- read.csv(manifest_path, stringsAsFactors = FALSE, check.names = FALSE)

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
  path <- file.path(project_root, row$LocalRelativePath)
  remote_match <- remote[
    remote$PXD == row$PXD &
      remote$FileName == row$FileName,
  ]
  remote_size <- if (nrow(remote_match)) {
    remote_match$RemoteSizeBytes[[1]]
  } else {
    NA_real_
  }
  size <- if (file.exists(path)) file.info(path)$size else 0
  complete <- file.exists(path) &&
    (is.na(remote_size) || size == remote_size)
  data.frame(
    PXD = row$PXD,
    FileName = row$FileName,
    FileRole = row$FileRole,
    LocalRelativePath = row$LocalRelativePath,
    SizeBytes = size,
    SHA256 = if (complete) sha256_file(path) else "",
    SourceURL = row$SourceURL,
    RemoteSizeBytes = remote_size,
    SizeMatchesProbe = complete,
    DownloadStatus = if (complete) "downloaded" else "incomplete_or_missing",
    stringsAsFactors = FALSE
  )
})

additional <- bind_rows(records)
combined <- bind_rows(
  existing |>
    filter(
      !paste(PXD, FileName) %in%
        paste(additional$PXD, additional$FileName)
    ),
  additional
) |>
  arrange(PXD, FileName)

write.csv(combined, manifest_path, row.names = FALSE, na = "")
if (!all(additional$SizeMatchesProbe)) {
  stop("One or more additional files are incomplete or missing", call. = FALSE)
}
message("Registered ", nrow(additional), " additional lactylome pair files.")
