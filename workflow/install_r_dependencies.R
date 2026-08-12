#!/usr/bin/env Rscript

# Install the R dependencies declared by the publication workflow.
# Exact versions used for a run are recorded by workflow/record_environment.R.

cran_packages <- c(
  "data.table", "digest", "dplyr", "eulerr", "ggplot2", "Matrix",
  "patchwork", "readr", "readxl", "Rtsne", "stringr", "tidyr", "uwot"
)
missing_cran <- cran_packages[
  !vapply(cran_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_cran)) {
  install.packages(missing_cran, repos = "https://cloud.r-project.org")
}

if (!requireNamespace("GO.db", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }
  BiocManager::install("GO.db", ask = FALSE, update = FALSE)
}

required <- c(cran_packages, "GO.db")
missing <- required[
  !vapply(required, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing)) {
  stop("R package installation incomplete: ", paste(missing, collapse = ", "))
}

message("All R dependencies are available.")
