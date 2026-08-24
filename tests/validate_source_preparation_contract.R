#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args)) normalizePath(args[[1L]], mustWork = TRUE) else normalizePath(".", mustWork = TRUE)
path <- file.path(
  project_root, "workflow", "source_legacy", "python",
  "data_preparation", "build_reference_proteome_membership.py"
)
lines <- readLines(path, warn = FALSE)
analyzer_lines <- readLines(
  file.path(project_root, "workflow", "source_legacy", "R", "analyze_ddr_fraction.R"),
  warn = FALSE
)
preparation_lines <- readLines(
  file.path(project_root, "workflow", "prepare_from_source.R"),
  warn = FALSE
)
core_lines <- readLines(
  file.path(project_root, "workflow", "source_legacy", "python", "data_preparation", "build_core_kla_inputs.py"),
  warn = FALSE
)
assert <- function(condition, message) if (!isTRUE(condition)) stop(message, call. = FALSE)
assert(
  any(grepl("if KLA_STATS_PATH.exists():", lines, fixed = TRUE)),
  "Source preparation must not require a stale Kla statistics table before reference extraction."
)
assert(
  any(grepl("if len(statistics_frame) != 10", lines, fixed = TRUE)),
  "Reference-statistics completeness must be checked independently of an optional Kla comparison table."
)
assert(
  !any(grepl("if len(comparison) != 10", lines, fixed = TRUE)),
  "An absent optional Kla comparison table must not fail reference-proteome extraction."
)
assert(
  any(grepl("ensembl_mapping_identifiers", analyzer_lines, fixed = TRUE)),
  "Frozen Ensembl records must retain known unmapped identifiers as well as mapped identifiers."
)
assert(
  !any(grepl("KLA_STATS_PATH", lines[grep("^def build_source_manifest", lines):grep("^def main", lines)], fixed = TRUE)),
  "A stale Kla statistics table must not be recorded as a reference-source file."
)
assert(
  any(grepl("group_summary_30.csv", analyzer_lines, fixed = TRUE)),
  "Source membership construction must take its scope from the frozen final 30-group contract."
)
assert(
  !any(grepl("kla_regulator_intensity_availability_audit.csv", analyzer_lines, fixed = TRUE)),
  "Source membership construction must not require a regulator-analysis cache."
)
assert(
  !any(grepl("..names(groups)", preparation_lines, fixed = TRUE)),
  "Source-table column selection must use a materialized column-name vector."
)
assert(
  any(grepl("venn_analytical_columns", preparation_lines, fixed = TRUE)),
  "Source-derived Venn memberships must be compared separately from frozen display annotations."
)
assert(
  any(grepl("expand_reference_membership", preparation_lines, fixed = TRUE)),
  "Reference Venn memberships must expand stable BaseAccessions from their source-protein records."
)
assert(
  any(grepl("ddr_accession_set", preparation_lines, fixed = TRUE)),
  "Reference Venn DDR membership must be assigned after BaseAccession expansion from the frozen DDR annotation."
)
assert(
  any(grepl("GENE PRODUCT DB", preparation_lines, fixed = TRUE)) &&
    any(grepl("QUALIFIER", preparation_lines, fixed = TRUE)),
  "Source Venn DDR membership must preserve the frozen UniProtKB and non-NOT annotation filters."
)
assert(
  !any(grepl('extract_pxd014870(', core_lines, fixed = TRUE)),
  "A source parser outside the final 30-group publication scope must not run."
)
message("PASS: source preparation has no stale-reference-cache dependency.")
