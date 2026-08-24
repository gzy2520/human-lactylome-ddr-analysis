#!/usr/bin/env Rscript

# Prepare the exact 30-group publication input from downloaded processed
# ProteomeXchange/iProX source files.  The historical per-dataset parsers run
# in a disposable staging directory; only the candidate publication input is
# retained under work/.  Manual source assets (S4, S5, regulator percentiles
# and the frozen human DDR annotation) are deliberately not recalculated.

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
command <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", command, value = TRUE)
script_path <- if (length(file_arg)) {
  normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE)
} else {
  normalizePath("workflow/prepare_from_source.R", mustWork = TRUE)
}
project_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
source_dir <- if (length(args)) args[[1L]] else Sys.getenv(
  "KLA_SOURCE_DATA",
  unset = file.path(project_root, "data", "source_cache")
)
source_dir <- normalizePath(source_dir, mustWork = TRUE)
frozen_input <- normalizePath(
  file.path(project_root, "data", "publication_input"),
  mustWork = TRUE
)
legacy_dir <- normalizePath(
  file.path(project_root, "workflow", "source_legacy"),
  mustWork = TRUE
)
candidate_dir <- file.path(project_root, "work", "source_publication_input")
stage_dir <- file.path(project_root, "work", "source_preparation_stage")

assert <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}
is_true <- function(x) as.character(x) %in% c("TRUE", "True", "true", "T", "1")
key <- function(data) paste(data$PXD, data$SampleGroup, sep = "__")

link_path <- function(target, link) {
  assert(file.exists(target), paste("Missing source required for link:", target))
  assert(!file.exists(link) && !dir.exists(link), paste("Staging path already exists:", link))
  assert(isTRUE(file.symlink(target, link)), paste("Could not create staging link:", link))
}

run <- function(executable, arguments, label) {
  message("[", label, "]")
  status <- system2(executable, arguments)
  if (!identical(status, 0L)) stop("Source preparation step failed: ", label, call. = FALSE)
}

compare_table <- function(observed, expected, sort_columns, label) {
  assert(identical(names(observed), names(expected)), paste(label, "has an unexpected schema."))
  observed <- copy(as.data.table(observed))
  expected <- copy(as.data.table(expected))
  setorderv(observed, sort_columns)
  setorderv(expected, sort_columns)
  for (column in names(expected)) {
    set(observed, j = column, value = fifelse(is.na(observed[[column]]), "", as.character(observed[[column]])))
    set(expected, j = column, value = fifelse(is.na(expected[[column]]), "", as.character(expected[[column]])))
  }
  assert(identical(observed, expected), paste(label, "does not reproduce the frozen release table."))
}

region_name <- function(present) {
  if (length(present) == 4L) return("all_four")
  paste0(paste(present, collapse = "_and_"), "_only")
}

build_venn <- function(membership, groups, ddr_only) {
  categories <- c("normal_tissue", "cancer_tissue", "cancer_cells", "normal_cells")
  membership <- copy(as.data.table(membership))
  if (ddr_only) membership <- membership[is_true(IsDdr)]
  membership <- unique(membership[, .(PXD, SampleGroup, BaseAccession)])
  membership <- merge(
    membership,
    groups[, .(PXD, SampleGroup, Category)],
    by = c("PXD", "SampleGroup"),
    all.x = TRUE,
    sort = FALSE
  )
  assert(!anyNA(membership$Category), "A source membership row is outside the final 30 groups.")
  accessions <- sort(unique(membership$BaseAccession))
  out <- data.table(BaseAccession = accessions)
  for (category in categories) {
    present <- unique(membership[Category == category, BaseAccession])
    out[, (paste0("In_", category)) := BaseAccession %in% present]
  }
  out[, Region := vapply(seq_len(.N), function(index) {
    present <- categories[vapply(categories, function(category) {
      isTRUE(out[[paste0("In_", category)]][[index]])
    }, logical(1))]
    region_name(present)
  }, character(1))]
  setcolorder(out, c("BaseAccession", "Region", paste0("In_", categories)))
  setorderv(out, "BaseAccession")
  out
}

ddr_accession_set <- function(annotation_path) {
  annotations <- fread(annotation_path)
  assert(
    all(c("GENE PRODUCT DB", "GENE PRODUCT ID", "QUALIFIER", "TAXON ID") %in% names(annotations)),
    "The frozen DDR annotation is missing required provenance columns."
  )
  taxon <- suppressWarnings(as.integer(annotations[["TAXON ID"]]))
  qualifier <- as.character(annotations[["QUALIFIER"]])
  keep <- annotations[["GENE PRODUCT DB"]] == "UniProtKB" &
    (is.na(taxon) | taxon == 9606L) &
    !grepl("(^|\\|)NOT(\\||$)", qualifier)
  accessions <- unlist(strsplit(
    as.character(annotations[["GENE PRODUCT ID"]][keep]), "[;,]"
  ))
  accessions <- trimws(accessions)
  accessions <- sub("^.*\\|([^|]+)\\|.*$", "\\1", accessions)
  accessions <- sub("^([^|;]+)\\|.*$", "\\1", accessions)
  accessions <- sub("^NX_", "", accessions)
  accessions <- sub("-[0-9]+$", "", accessions)
  valid_uniprot <- grepl(
    "^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9](?:[A-Z0-9]{3}[0-9])?)$",
    accessions
  )
  accessions <- sort(unique(accessions[!is.na(accessions) & nzchar(accessions) & valid_uniprot]))
  assert(length(accessions) > 0L, "The frozen DDR annotation contains no BaseAccessions.")
  accessions
}

expand_reference_membership <- function(membership, ddr_accessions) {
  membership <- copy(as.data.table(membership))
  required <- c("PXD", "SampleGroup", "MappedBaseAccessions", "IsDdr")
  assert(
    all(required %in% names(membership)),
    "Reference membership is missing its source-to-BaseAccession mapping columns."
  )
  membership <- membership[
    !is.na(MappedBaseAccessions) & nzchar(trimws(MappedBaseAccessions))
  ]
  expanded <- membership[, .(
    BaseAccession = unlist(strsplit(MappedBaseAccessions, ";", fixed = TRUE))
  ), by = .(PXD, SampleGroup, IsDdr)]
  expanded[, BaseAccession := trimws(BaseAccession)]
  expanded <- unique(expanded[nzchar(BaseAccession)])
  expanded[, IsDdr := BaseAccession %in% ddr_accessions]
  expanded
}

unlink(stage_dir, recursive = TRUE, force = TRUE)
dir.create(file.path(stage_dir, "data"), recursive = TRUE, showWarnings = FALSE)
on.exit(unlink(stage_dir, recursive = TRUE, force = TRUE), add = TRUE)

link_path(file.path(legacy_dir, "config"), file.path(stage_dir, "config"))
link_path(file.path(legacy_dir, "python"), file.path(stage_dir, "python"))
link_path(file.path(legacy_dir, "R"), file.path(stage_dir, "R"))
link_path(frozen_input, file.path(stage_dir, "data", "publication_input"))
dir.create(file.path(stage_dir, "data", "annotations"), recursive = TRUE, showWarnings = FALSE)
link_path(
  file.path(frozen_input, "human_ddr_go_annotations.tsv"),
  file.path(stage_dir, "data", "annotations", "GO-repair+damage(human).tsv")
)

pxd_dirs <- list.dirs(source_dir, recursive = FALSE, full.names = TRUE)
pxd_dirs <- pxd_dirs[grepl("^PXD[0-9]+$", basename(pxd_dirs))]
assert(length(pxd_dirs) > 0L, "The supplied source directory contains no PXD folders.")
for (path in pxd_dirs) {
  link_path(path, file.path(stage_dir, "data", basename(path)))
}

run(
  "python3",
  c(file.path(stage_dir, "python", "data_preparation", "build_core_kla_inputs.py"),
    "--project-root", stage_dir),
  "core Kla source extraction"
)
run(
  "python3",
  c(file.path(stage_dir, "python", "data_preparation", "build_reference_proteome_membership.py"),
    "--project-root", stage_dir),
  "whole-proteome source extraction"
)
run(
  "Rscript",
  c(file.path(stage_dir, "R", "analyze_ddr_fraction.R"), stage_dir),
  "30-group BaseAccession membership"
)

groups <- fread(file.path(frozen_input, "group_summary_30.csv"))
assert(nrow(groups) == 30L && !anyDuplicated(groups[, .(PXD, SampleGroup)]),
  "The frozen release must define exactly 30 unique groups.")
target_keys <- key(groups)
publication_columns <- names(groups)
source_ddr_accessions <- ddr_accession_set(
  file.path(frozen_input, "human_ddr_go_annotations.tsv")
)

source_summary <- fread(file.path(
  stage_dir, "results", "tables",
  "cell_type_kla_vs_reference_ddr_statistics_accession_only.csv"
))
source_summary <- source_summary[match(target_keys, key(source_summary))]
assert(!anyNA(source_summary$PXD), "A final group was not recovered from downloaded source data.")
source_summary <- source_summary[, ..publication_columns]
compare_table(source_summary, groups, c("RowOrder", "PXD", "SampleGroup"), "group_summary_30.csv")

source_kla <- fread(file.path(
  stage_dir, "work", "intermediate", "expanded_ddr_by_accession",
  "kla_proteins_by_sample_group.csv"
))
source_kla <- source_kla[key(source_kla) %in% target_keys]
frozen_kla <- fread(file.path(frozen_input, "kla_protein_membership_30.csv"))
compare_table(source_kla, frozen_kla, c("PXD", "SampleGroup", "BaseAccession"),
  "kla_protein_membership_30.csv")

source_reference <- fread(file.path(
  stage_dir, "work", "intermediate", "expanded_ddr_by_accession",
  "reference_proteins_by_sample_group.csv"
))
source_reference <- source_reference[key(source_reference) %in% target_keys]
frozen_reference <- fread(file.path(frozen_input, "reference_protein_membership_30.csv"))
compare_table(source_reference, frozen_reference,
  c("PXD", "SampleGroup", "SourceProteinID"), "reference_protein_membership_30.csv")
source_reference_for_venn <- expand_reference_membership(
  source_reference, source_ddr_accessions
)

venn_specs <- list(
  list(file = "venn_all_kla.csv", membership = source_kla, ddr = FALSE),
  list(file = "venn_kla_ddr.csv", membership = source_kla, ddr = TRUE),
  list(file = "venn_reference.csv", membership = source_reference_for_venn, ddr = FALSE),
  list(file = "venn_reference_ddr.csv", membership = source_reference_for_venn, ddr = TRUE)
)
venn_analytical_columns <- c(
  "BaseAccession", paste0("In_", c(
    "normal_tissue", "cancer_tissue", "cancer_cells", "normal_cells"
  )), "Region"
)
for (spec in venn_specs) {
  observed <- build_venn(spec$membership, groups, spec$ddr)
  expected <- fread(file.path(frozen_input, spec$file))
  assert(
    all(venn_analytical_columns %in% names(expected)),
    paste(spec$file, "is missing its analytical Venn columns.")
  )
  observed <- observed[, ..venn_analytical_columns]
  expected <- expected[, ..venn_analytical_columns]
  compare_table(observed, expected, "BaseAccession", spec$file)
}

unlink(candidate_dir, recursive = TRUE, force = TRUE)
dir.create(candidate_dir, recursive = TRUE, showWarnings = FALSE)
for (filename in list.files(frozen_input, full.names = FALSE, no.. = TRUE)) {
  assert(
    isTRUE(file.copy(file.path(frozen_input, filename), file.path(candidate_dir, filename))),
    paste("Could not materialize validated release input:", filename)
  )
}
message("PASS: downloaded-source processing reproduced all source-derived 30-group tables.")
message("Validated publication input: ", candidate_dir)
