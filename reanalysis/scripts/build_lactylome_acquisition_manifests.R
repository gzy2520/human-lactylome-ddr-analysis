#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(jsonlite)
  library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) normalizePath(args[[1]]) else normalizePath(".")
table_dir <- file.path(project_root, "reanalysis", "results", "tables")
config_dir <- file.path(project_root, "reanalysis", "config")

selected_pxd <- c(
  "PXD036307",
  "PXD054919",
  "PXD063047",
  "PXD064912",
  "PXD066054",
  "PXD075377"
)

inventory <- read.csv(
  file.path(table_dir, "human_lactylome_mass_spectrometry_inventory.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
repository_files <- read.csv(
  file.path(table_dir, "human_lactylome_repository_file_manifest.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
validation <- read.csv(
  file.path(config_dir, "lactylome_acquisition_validation.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

relative_to_project <- function(path) {
  str_remove(normalizePath(path, mustWork = FALSE), paste0("^", fixed(project_root), "/"))
}

repository_alias <- function(file_name, file_url) {
  url_parent <- basename(dirname(sub("\\?.*$", "", file_url)))
  extension <- tools::file_ext(file_name)
  stem <- tools::file_path_sans_ext(file_name)
  if (nzchar(extension)) {
    paste0(stem, "_", url_parent, ".", extension)
  } else {
    paste0(stem, "_", url_parent)
  }
}

classify_role <- function(relative_path) {
  case_when(
    str_detect(relative_path, "/raw/") ~ "raw",
    str_detect(relative_path, "/search_results/") ~ "search_result",
    str_detect(relative_path, "/supplementary/") ~ "article_supplement",
    str_detect(relative_path, "/metadata/articles/") ~ "article_original",
    str_detect(relative_path, "/metadata/") ~ "metadata",
    TRUE ~ "other"
  )
}

source_files_for_pxd <- function(pxd) {
  pxd_dir <- file.path(project_root, "data", pxd)
  if (!dir.exists(pxd_dir)) {
    return(character())
  }
  files <- list.files(
    pxd_dir,
    recursive = TRUE,
    full.names = TRUE,
    all.files = TRUE,
    no.. = TRUE
  )
  files <- files[file.info(files)$isdir %in% FALSE]
  files <- files[!str_detect(files, "/search_results/extracted/")]
  files <- files[!str_detect(files, "\\.partial_download$")]
  files <- files[!basename(files) %in% c(
    "dataset_metadata.csv",
    "repository_file_manifest.csv",
    "download_manifest.csv",
    "extracted_file_inventory.csv"
  )]
  files
}

find_repository_record <- function(pxd, local_name) {
  candidates <- repository_files |>
    filter(PXD == pxd) |>
    mutate(LocalAlias = mapply(repository_alias, FileName, FileURL, USE.NAMES = FALSE)) |>
    filter(FileName == local_name | LocalAlias == local_name)
  if (nrow(candidates) == 0) {
    return(tibble(
      RepositoryFileCategory = "",
      SourceURL = "",
      DownloadPriority = ""
    ))
  }
  candidates |>
    slice(1) |>
    transmute(
      RepositoryFileCategory = FileCategory,
      SourceURL = FileURL,
      DownloadPriority
    )
}

all_downloads <- list()
for (pxd in selected_pxd) {
  metadata_dir <- file.path(project_root, "data", pxd, "metadata")
  dir.create(metadata_dir, recursive = TRUE, showWarnings = FALSE)

  dataset_row <- inventory |>
    filter(PXD == pxd)
  write.csv(
    dataset_row,
    file.path(metadata_dir, "dataset_metadata.csv"),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )

  repository_subset <- repository_files |>
    filter(PXD == pxd)
  write.csv(
    repository_subset,
    file.path(metadata_dir, "repository_file_manifest.csv"),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )

  local_files <- source_files_for_pxd(pxd)
  local_rows <- lapply(local_files, function(path) {
    local_name <- basename(path)
    source_row <- find_repository_record(pxd, local_name)
    validation_row <- validation |>
      filter(PXD == pxd, FileName == local_name) |>
      slice(1)
    if (nrow(validation_row) == 0) {
      validation_row <- tibble(
        SHA256 = "",
        RepositorySHA1 = "",
        ValidationStatus = "downloaded_not_independently_validated",
        ValidationNote = ""
      )
    }
    relative_path <- relative_to_project(path)
    tibble(
      PXD = pxd,
      LocalRelativePath = relative_path,
      FileName = local_name,
      FileRole = classify_role(relative_path),
      SizeBytes = unname(file.info(path)$size),
      SHA256 = validation_row$SHA256,
      RepositorySHA1 = validation_row$RepositorySHA1,
      ValidationStatus = validation_row$ValidationStatus,
      ValidationNote = validation_row$ValidationNote,
      RepositoryFileCategory = source_row$RepositoryFileCategory,
      SourceURL = source_row$SourceURL,
      DownloadPriority = source_row$DownloadPriority
    )
  })
  local_manifest <- bind_rows(local_rows) |>
    arrange(FileRole, FileName)
  write.csv(
    local_manifest,
    file.path(metadata_dir, "download_manifest.csv"),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )
  all_downloads[[pxd]] <- local_manifest

  extracted_root <- file.path(project_root, "data", pxd, "search_results", "extracted")
  extracted_files <- if (dir.exists(extracted_root)) {
    list.files(extracted_root, recursive = TRUE, full.names = TRUE)
  } else {
    character()
  }
  extracted_files <- extracted_files[file.info(extracted_files)$isdir %in% FALSE]
  extracted_manifest <- tibble(
    PXD = pxd,
    LocalRelativePath = vapply(extracted_files, relative_to_project, character(1)),
    FileName = basename(extracted_files),
    SizeBytes = unname(file.info(extracted_files)$size),
    EvidenceFile = str_detect(
      basename(extracted_files),
      regex("La \\(K\\)Sites|PTMSiteReport|MS_identified_information|evidence|proteinGroups|Identification", ignore_case = TRUE)
    )
  ) |>
    arrange(desc(EvidenceFile), LocalRelativePath)
  write.csv(
    extracted_manifest,
    file.path(metadata_dir, "extracted_file_inventory.csv"),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )

  proxi_path <- file.path(metadata_dir, "proxi_detail.json")
  if (file.exists(proxi_path)) {
    parsed <- tryCatch(fromJSON(proxi_path, simplifyVector = FALSE), error = function(error) NULL)
    if (is.null(parsed)) {
      warning("Invalid PROXI JSON for ", pxd)
    }
  }
}

combined_downloads <- bind_rows(all_downloads) |>
  arrange(PXD, FileRole, FileName)
write.csv(
  combined_downloads,
  file.path(table_dir, "priority_dataset_acquisition_manifest.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

message("Acquisition manifest rows: ", nrow(combined_downloads))
message("Per-PXD metadata written for: ", paste(selected_pxd, collapse = ", "))
