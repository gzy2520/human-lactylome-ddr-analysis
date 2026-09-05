#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args)) normalizePath(args[[1L]], mustWork = TRUE) else normalizePath(".", mustWork = TRUE)

stop_if <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

groups <- fread(file.path(project_root, "data", "publication_input", "group_summary_30.csv"))
design <- fread(file.path(project_root, "data", "candidate", "sample_design_30.csv"), na.strings = c("", "NA"))
required <- c(
  "RowOrder", "PXD", "SampleGroup", "Category", "KlaN", "ReferenceN",
  "KlaSampleDesign", "ReferenceSampleDesign", "Aggregation", "MatchClass"
)
stop_if(identical(names(design), required), "Candidate sample design schema changed.")
stop_if(nrow(groups) == 30L && nrow(design) == 30L, "Candidate sample design must contain 30 rows.")
stop_if(!anyDuplicated(design[, .(PXD, SampleGroup)]), "Candidate sample design has duplicate group keys.")
stop_if(
  setequal(
    paste(groups$PXD, groups$SampleGroup, sep = "__"),
    paste(design$PXD, design$SampleGroup, sep = "__")
  ),
  "Candidate sample design does not match the frozen 30-group scope."
)
stop_if(
  !any(!is.na(design$KlaN) & (design$KlaN < 0 | design$KlaN != floor(design$KlaN))),
  "Invalid Kla sample count."
)
stop_if(
  !any(!is.na(design$ReferenceN) & (design$ReferenceN < 0 | design$ReferenceN != floor(design$ReferenceN))),
  "Invalid reference sample count."
)
stop_if(all(nchar(trimws(design$KlaSampleDesign)) > 0), "A Kla sample-design label is empty.")
stop_if(all(nchar(trimws(design$ReferenceSampleDesign)) > 0), "A reference sample-design label is empty.")
stop_if(all(nchar(trimws(design$Aggregation)) > 0), "An aggregation label is empty.")
stop_if(all(nchar(trimws(design$MatchClass)) > 0), "A match-class label is empty.")
message("PASS: candidate sample design covers the exact frozen 30-group scope.")
