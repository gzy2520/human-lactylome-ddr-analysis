#!/usr/bin/env Rscript

# The two supported paths converge at the same validated 30-group input:
# `publication` uses the versioned release tables; `source` first rebuilds
# their source-derived fields from downloaded processed proteomics files.

args <- commandArgs(trailingOnly = TRUE)
target <- if (length(args)) args[[1L]] else "publication"
if (!target %in% c("publication", "source", "validate")) {
  stop("Target must be 'publication', 'source', or 'validate'.", call. = FALSE)
}

command <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", command, value = TRUE)
script_path <- if (length(file_arg)) {
  normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE)
} else {
  normalizePath("workflow/run_pipeline.R", mustWork = TRUE)
}
project_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)

run <- function(executable, arguments, label, environment = character()) {
  message("[", label, "]")
  status <- system2(executable, arguments, env = environment)
  if (!identical(status, 0L)) stop("Workflow step failed: ", label, call. = FALSE)
}

input_dir <- file.path(project_root, "data", "publication_input")
if (target == "source") {
  source_arguments <- c(file.path(project_root, "workflow", "prepare_from_source.R"))
  if (length(args) >= 2L) source_arguments <- c(source_arguments, args[[2L]])
  run("Rscript", source_arguments, "source data to validated publication input")
  input_dir <- file.path(project_root, "work", "source_publication_input")
}
input_dir <- normalizePath(input_dir, mustWork = TRUE)
input_environment <- paste0("KLA_PUBLICATION_INPUT=", input_dir)

if (target %in% c("publication", "source")) {
  # Generated results are a closed publication set.  Clear only the two
  # declared output directories so a previous run cannot leave an extra panel
  # or supplementary workbook behind.
  unlink(file.path(project_root, "results", "figures"), recursive = TRUE, force = TRUE)
  unlink(file.path(project_root, "results", "supplementary"), recursive = TRUE, force = TRUE)
  run("Rscript", c(file.path(project_root, "R", "publication", "build_publication_outputs.R"), project_root), "manuscript figures", input_environment)
  run("Rscript", c(file.path(project_root, "workflow", "build_supplementary_workbooks.R"), project_root), "Supplementary Tables S1-S6", input_environment)
}

run("Rscript", c(file.path(project_root, "tests", "validate_publication_contract.R"), project_root), "publication contract", input_environment)
message("PASS: final 30-group publication workflow completed.")
