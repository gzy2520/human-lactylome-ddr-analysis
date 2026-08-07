#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) normalizePath(args[[1]]) else normalizePath(".")

manifest_path <- file.path(
  project_root,
  "reanalysis",
  "results",
  "tables",
  "human_lactylome_repository_file_manifest.csv"
)
output_path <- file.path(
  project_root,
  "reanalysis",
  "results",
  "tables",
  "lactylome_pair_remote_file_sizes.csv"
)

priority_pxd <- c(
  "PXD028737", "PXD033146", "PXD037371", "PXD037530", "PXD045967",
  "PXD046344", "PXD046800", "PXD047535", "PXD047673", "PXD048995",
  "PXD050147", "PXD052772", "PXD053029", "PXD055230", "PXD057709",
  "PXD058534", "PXD062720", "PXD063266", "PXD063945", "PXD064038",
  "PXD065104", "PXD065831", "PXD066351", "PXD068838", "PXD070007",
  "PXD070427", "PXD073311", "PXD075014", "PXD077426"
)

manifest <- read.csv(manifest_path, stringsAsFactors = FALSE, check.names = FALSE)
selected <- manifest[
  manifest$PXD %in% priority_pxd &
    manifest$FileCategory != "Associated raw file URI" &
    nzchar(manifest$FileURL),
]

selected$CanonicalURL <- sub("^ftp://", "https://", selected$FileURL)
selected$ProbeStatus <- ""
selected$RemoteSizeBytes <- NA_real_
selected$RemoteSizeMiB <- NA_real_
selected$LastModified <- ""

probe_one <- function(url) {
  command <- c(
    "-sL", "--max-time", "90", "-r", "0-0", "-D", "-", "-o", "/dev/null",
    url
  )
  command[[length(command)]] <- shQuote(command[[length(command)]])
  headers <- tryCatch(
    system2("curl", command, stdout = TRUE, stderr = TRUE),
    error = function(e) paste("ERROR", conditionMessage(e))
  )
  status <- attr(headers, "status")
  content_range <- grep("^Content-Range:", headers, ignore.case = TRUE, value = TRUE)
  content_length <- grep("^Content-Length:", headers, ignore.case = TRUE, value = TRUE)
  last_modified <- grep("^Last-Modified:", headers, ignore.case = TRUE, value = TRUE)
  http_status <- grep("^HTTP/", headers, value = TRUE)

  size <- NA_real_
  if (length(content_range)) {
    size <- suppressWarnings(as.numeric(sub(".*/", "", tail(content_range, 1))))
  } else if (length(content_length)) {
    size <- suppressWarnings(as.numeric(sub("^[^:]+:[[:space:]]*", "", tail(content_length, 1))))
  }

  list(
    status = if (is.null(status) || status == 0) {
      paste(tail(http_status, 1), collapse = "")
    } else {
      paste0("curl_exit_", status)
    },
    size = size,
    last_modified = if (length(last_modified)) {
      sub("^[^:]+:[[:space:]]*", "", tail(last_modified, 1))
    } else {
      ""
    }
  )
}

message("Probing ", nrow(selected), " processed/metadata files...")
for (i in seq_len(nrow(selected))) {
  result <- probe_one(selected$CanonicalURL[[i]])
  selected$ProbeStatus[[i]] <- result$status
  selected$RemoteSizeBytes[[i]] <- result$size
  selected$RemoteSizeMiB[[i]] <- if (is.na(result$size)) NA_real_ else result$size / 1024^2
  selected$LastModified[[i]] <- result$last_modified
  message(
    sprintf(
      "[%d/%d] %s %s %.1f MiB",
      i,
      nrow(selected),
      selected$PXD[[i]],
      selected$FileName[[i]],
      selected$RemoteSizeMiB[[i]]
    )
  )
}

selected <- selected[
  order(selected$PXD, selected$FileCategory, selected$FileName, selected$CanonicalURL),
]
write.csv(selected, output_path, row.names = FALSE, na = "")
message("Wrote ", output_path)
