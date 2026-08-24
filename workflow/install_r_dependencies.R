#!/usr/bin/env Rscript

# Restore the exact R package environment used for this frozen publication
# release.  Installing current CRAN package versions is intentionally avoided:
# a package update may change a figure layout or workbook serialization.
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg)) {
  normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE)
} else {
  normalizePath("workflow/install_r_dependencies.R", mustWork = TRUE)
}
project_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)

required_r <- "4.4.3"
if (as.character(getRversion()) != required_r) {
  stop("This release requires R ", required_r, "; found R ", getRversion(), ".", call. = FALSE)
}

required_renv <- "1.1.8"
if (!requireNamespace("renv", quietly = TRUE) ||
    as.character(utils::packageVersion("renv")) != required_renv) {
  install.packages(
    "https://cran.r-project.org/src/contrib/Archive/renv/renv_1.1.8.tar.gz",
    repos = NULL,
    type = "source"
  )
}
if (!requireNamespace("renv", quietly = TRUE) ||
    as.character(utils::packageVersion("renv")) != required_renv) {
  stop("Unable to install renv ", required_renv, ".", call. = FALSE)
}

renv::restore(project = project_root, prompt = FALSE)
renv::activate(project = project_root)
status <- renv::status(project = project_root)
if (!isTRUE(status$synchronized)) {
  stop("The locked publication environment did not restore cleanly.", call. = FALSE)
}
message("Locked R 4.4.3 publication environment restored.")
