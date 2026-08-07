#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) normalizePath(args[[1]]) else normalizePath(".")
manifest_path <- file.path(
  project_root,
  "reanalysis",
  "results",
  "tables",
  "lactylome_pair_download_manifest.csv"
)
output_path <- file.path(
  project_root,
  "reanalysis",
  "results",
  "tables",
  "lactylome_pair_extraction_manifest.csv"
)

manifest <- read.csv(manifest_path, stringsAsFactors = FALSE, check.names = FALSE)
archives <- manifest[
  grepl("\\.(zip|rar)$", manifest$FileName, ignore.case = TRUE) &
    file.exists(file.path(project_root, manifest$LocalRelativePath)),
]

records <- list()
for (i in seq_len(nrow(archives))) {
  archive <- file.path(project_root, archives$LocalRelativePath[[i]])
  output_dir <- file.path(
    project_root,
    "data",
    archives$PXD[[i]],
    "search_results",
    "extracted_pairing",
    sub("\\.(zip|rar)$", "", archives$FileName[[i]], ignore.case = TRUE)
  )
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  message(sprintf("[%d/%d] extracting %s", i, nrow(archives), archive))
  status <- if (grepl("\\.zip$", archive, ignore.case = TRUE)) {
    system2("unzip", c("-q", "-n", shQuote(archive), "-d", shQuote(output_dir)))
  } else {
    system2("tar", c("-xf", shQuote(archive), "-C", shQuote(output_dir)))
  }
  files <- list.files(output_dir, recursive = TRUE, full.names = TRUE)
  records[[i]] <- data.frame(
    PXD = archives$PXD[[i]],
    Archive = archives$FileName[[i]],
    OutputRelativePath = substring(output_dir, nchar(project_root) + 2),
    ExtractionStatus = if (status == 0) "archive_extracted" else paste0("extract_exit_", status),
    ExtractedFileCount = length(files),
    ExtractedBytes = if (length(files)) sum(file.info(files)$size, na.rm = TRUE) else 0,
    stringsAsFactors = FALSE
  )
}

result <- if (length(records)) {
  do.call(rbind, records)
} else {
  data.frame(
    PXD = character(),
    Archive = character(),
    OutputRelativePath = character(),
    ExtractionStatus = character(),
    ExtractedFileCount = integer(),
    ExtractedBytes = numeric(),
    stringsAsFactors = FALSE
  )
}
write.csv(result, output_path, row.names = FALSE)
message("Wrote ", output_path)
