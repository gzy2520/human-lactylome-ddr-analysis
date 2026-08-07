#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) normalizePath(args[[1]]) else normalizePath(".")
size_path <- file.path(
  project_root,
  "reanalysis",
  "results",
  "tables",
  "lactylome_pair_remote_file_sizes.csv"
)
output_path <- file.path(
  project_root,
  "reanalysis",
  "results",
  "tables",
  "lactylome_pair_download_manifest.csv"
)

remote <- read.csv(size_path, stringsAsFactors = FALSE, check.names = FALSE)

selected_files <- list(
  PXD028737 = c("txt.zip"),
  PXD033146 = c("search_result-HA119TPLa.zip", "search_result-HA119TQ.zip"),
  PXD037371 = c("txt.rar"),
  PXD037530 = c("checksum.txt"),
  PXD045967 = c("8.16-mqpar.xml"),
  PXD046800 = c(
    "HFX2_LFQ_QB001_Lacty_PeptideGroups.txt",
    "HFX2_LFQ_QB001_Lacty_Proteins.txt",
    "HFX2_LFQ_QB002_PeptideGroups.txt",
    "HFX2_LFQ_QB002_Proteins.txt"
  ),
  PXD047535 = c("checksum.txt"),
  PXD047673 = c("checksum.txt"),
  PXD048995 = c("mqpar.xml"),
  PXD050147 = c(
    "checksum.txt",
    "Kla_evidence.txt",
    "Lactyl_K_Sites.txt",
    "SIRT_proteinGroups.txt"
  ),
  PXD052772 = c("mqpar.xml", "txt.zip"),
  PXD055230 = c("checksum.txt"),
  PXD057709 = c("checksum.txt"),
  PXD058534 = c("checksum.txt", "txt.zip"),
  PXD062720 = c("checksum.txt", "txt.zip"),
  PXD063266 = c("checksum.txt", "LactylSites.xlsx", "mqpar.xml"),
  PXD063945 = c("lac_mqpar.xml", "pro_mqpar.xml"),
  PXD064038 = c("Clinical information of samples.docx", "txt.zip"),
  PXD065831 = c("YAS202408210011-1.rar"),
  PXD066351 = c(
    "XB01472B1DA-DIA_result.tsv",
    "XB01472B1DA-MSstats_Input.tsv",
    "XB01472B1DA-Protein_Quant.tsv",
    "XB01472B1DPAc-DIA_result.csv",
    "XB01472B1DPAc-MSstats_Input.csv",
    "XB01472B1DPLa-DIA_result.csv",
    "XB01472B1DPLa-MSstats_Input.csv"
  ),
  PXD068838 = c("checksum.txt"),
  PXD070007 = c("SA206LPLaB1_Annotation.xlsx"),
  PXD070427 = c("checksum.txt"),
  PXD073311 = c("Database_search_result.zip")
)

keep <- vapply(seq_len(nrow(remote)), function(i) {
  pxd <- remote$PXD[[i]]
  pxd %in% names(selected_files) && remote$FileName[[i]] %in% selected_files[[pxd]]
}, logical(1))
queue <- remote[keep,]

external <- data.frame(
  PXD = c("PXD066517", "PXD066517"),
  FileCategory = c("Search engine output file URI", "Other type file URI"),
  FileName = c("20240275.tsv", "Supplementary_Table.xlsx"),
  FileURL = c(
    "ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2025/07/PXD066517/20240275.tsv",
    "ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2025/07/PXD066517/Supplementary_Table.xlsx"
  ),
  LocalMatches = "",
  LocalStatus = "not_downloaded",
  DownloadPriority = "reference_proteome_first",
  RepositoryFileNote = "Normal human sperm conventional proteome reference",
  CanonicalURL = c(
    "https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/07/PXD066517/20240275.tsv",
    "https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/07/PXD066517/Supplementary_Table.xlsx"
  ),
  ProbeStatus = "repository_api_confirmed",
  RemoteSizeBytes = c(5560825, 12383625),
  RemoteSizeMiB = c(5560825, 12383625) / 1024^2,
  LastModified = "",
  stringsAsFactors = FALSE
)
queue <- rbind(queue, external[, names(queue)])

sha256_file <- function(path) {
  result <- system2("shasum", c("-a", "256", shQuote(path)), stdout = TRUE)
  sub("[[:space:]].*$", "", result[[1]])
}

download_one <- function(url, destination) {
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
  url <- gsub(" ", "%20", url, fixed = TRUE)
  command <- paste(
    "curl --fail --location --continue-at - --retry 3 --retry-delay 2 --silent --show-error",
    "--output", shQuote(destination),
    shQuote(url)
  )
  status <- system(command)
  if (status != 0) {
    stop("Download failed: ", url)
  }
}

records <- vector("list", nrow(queue))
for (i in seq_len(nrow(queue))) {
  pxd <- queue$PXD[[i]]
  role_dir <- if (grepl(
    "checksum|mqpar|Clinical information",
    queue$FileName[[i]],
    ignore.case = TRUE
  )) {
    "metadata"
  } else {
    "search_results"
  }
  destination_name <- queue$FileName[[i]]
  if (
    pxd == "PXD073311" &&
      sum(queue$PXD == pxd & queue$FileName == destination_name) > 1
  ) {
    component <- sub(".*/(IPX[0-9]+)/.*", "\\1", queue$CanonicalURL[[i]])
    destination_name <- paste0(component, "_", destination_name)
  }
  destination <- file.path(project_root, "data", pxd, role_dir, destination_name)
  preexisting <- file.exists(destination) &&
    file.info(destination)$size > 0 &&
    (
      is.na(queue$RemoteSizeBytes[[i]]) ||
        file.info(destination)$size == queue$RemoteSizeBytes[[i]]
    )
  if (!preexisting) {
    message(sprintf("[%d/%d] downloading %s %s", i, nrow(queue), pxd, destination_name))
    download_one(queue$CanonicalURL[[i]], destination)
  } else {
    message(sprintf("[%d/%d] present %s %s", i, nrow(queue), pxd, destination_name))
  }
  records[[i]] <- data.frame(
    PXD = pxd,
    FileName = destination_name,
    FileRole = role_dir,
    LocalRelativePath = substring(destination, nchar(project_root) + 2),
    SizeBytes = file.info(destination)$size,
    SHA256 = sha256_file(destination),
    SourceURL = queue$CanonicalURL[[i]],
    RemoteSizeBytes = queue$RemoteSizeBytes[[i]],
    SizeMatchesProbe = is.na(queue$RemoteSizeBytes[[i]]) ||
      file.info(destination)$size == queue$RemoteSizeBytes[[i]],
    DownloadStatus = if (preexisting) "already_present" else "downloaded",
    stringsAsFactors = FALSE
  )
}

result <- do.call(rbind, records)
write.csv(result, output_path, row.names = FALSE, na = "")

for (pxd in unique(result$PXD)) {
  per_pxd <- result[result$PXD == pxd,]
  per_pxd_path <- file.path(project_root, "data", pxd, "metadata", "pairing_download_manifest.csv")
  dir.create(dirname(per_pxd_path), recursive = TRUE, showWarnings = FALSE)
  write.csv(per_pxd, per_pxd_path, row.names = FALSE, na = "")
}

message("Wrote ", output_path)
