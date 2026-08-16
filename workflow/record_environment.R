#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args)) normalizePath(args[[1]]) else normalizePath(".")
output_dir <- file.path(project_root, "results", "provenance")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

packages <- c(
  "data.table", "digest", "dplyr", "eulerr", "ggplot2", "ggVennDiagram", "GO.db", "Matrix",
  "patchwork", "readr", "readxl", "Rtsne", "stringr", "tidyr", "uwot"
)
versions <- data.frame(
  Package = packages,
  Version = vapply(
    packages,
    function(package) {
      if (requireNamespace(package, quietly = TRUE)) {
        as.character(packageVersion(package))
      } else {
        NA_character_
      }
    },
    character(1)
  ),
  stringsAsFactors = FALSE
)
write.csv(
  versions,
  file.path(output_dir, "r_package_versions.csv"),
  row.names = FALSE,
  na = ""
)
writeLines(
  capture.output(sessionInfo()),
  file.path(output_dir, "session_info.txt"),
  useBytes = TRUE
)
message("Environment recorded in results/provenance.")
