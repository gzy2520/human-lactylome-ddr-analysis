#!/usr/bin/env Rscript

# The publication workflow intentionally has one target. It starts from the
# frozen 30-group input and creates only figures and Supplementary Tables S1-S6
# described in the final manuscript.

args <- commandArgs(trailingOnly = TRUE)
target <- if (length(args)) args[[1L]] else "publication"
if (!target %in% c("publication", "validate")) {
  stop("Target must be 'publication' or 'validate'.", call. = FALSE)
}

command <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", command, value = TRUE)
script_path <- if (length(file_arg)) {
  normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE)
} else {
  normalizePath("workflow/run_pipeline.R", mustWork = TRUE)
}
project_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)

run <- function(executable, arguments, label) {
  message("[", label, "]")
  status <- system2(executable, arguments)
  if (!identical(status, 0L)) stop("Workflow step failed: ", label, call. = FALSE)
}

if (target == "publication") {
  # Generated results are a closed publication set.  Clear only the two
  # declared output directories so a previous run cannot leave an extra panel
  # or supplementary workbook behind.
  unlink(file.path(project_root, "results", "figures"), recursive = TRUE, force = TRUE)
  unlink(file.path(project_root, "results", "supplementary"), recursive = TRUE, force = TRUE)
  run("Rscript", c(file.path(project_root, "R", "publication", "build_publication_outputs.R"), project_root), "manuscript figures")
  run("Rscript", c(file.path(project_root, "workflow", "build_supplementary_workbooks.R"), project_root), "Supplementary Tables S1-S6")
}

run("Rscript", c(file.path(project_root, "tests", "validate_publication_contract.R"), project_root), "publication contract")
message("PASS: final 30-group publication workflow completed.")
