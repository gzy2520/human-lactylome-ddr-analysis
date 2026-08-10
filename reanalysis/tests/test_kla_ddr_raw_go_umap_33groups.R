#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
script_path <- if (length(file_arg) > 0L) {
  normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE)
} else {
  normalizePath("reanalysis/tests/test_kla_ddr_raw_go_umap_33groups.R", mustWork = TRUE)
}
project_root <- normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = TRUE)

table_dir <- file.path(
  project_root,
  "reanalysis/results/tables/kla_ddr_raw_go_umap_33groups"
)
figure_dir <- file.path(
  project_root,
  "reanalysis/results/figures/kla_ddr_raw_go_umap_33groups"
)

required_outputs <- c(
  file.path(table_dir, "project_sample_group_scope_33.csv"),
  file.path(table_dir, "removed_sample_groups_4.csv"),
  file.path(table_dir, "scope_summary.csv"),
  file.path(table_dir, "kla_ddr_protein_scope.csv"),
  file.path(table_dir, "protein_raw_go_term_long.csv"),
  file.path(table_dir, "raw_go_term_summary.csv"),
  file.path(table_dir, "raw_go_pattern_summary.csv"),
  file.path(table_dir, "protein_raw_go_term_binary_matrix.csv"),
  file.path(table_dir, "umap_coordinates_fixed.csv"),
  file.path(table_dir, "umap_parameters.csv"),
  file.path(table_dir, "input_file_audit.csv"),
  file.path(table_dir, "umap_raw_go_33groups_bundle.rds"),
  file.path(table_dir, "session_info.txt"),
  file.path(figure_dir, "kla_ddr_raw_go_umap_33groups_solid_blue.png"),
  file.path(figure_dir, "kla_ddr_raw_go_umap_33groups_solid_blue.pdf"),
  file.path(figure_dir, "kla_ddr_raw_go_umap_33groups_solid_blue.svg"),
  file.path(project_root, "reanalysis/reports/UMAP_RAW_GO_33GROUP_DATA_SCOPE.md")
)
stopifnot(all(file.exists(required_outputs)))
stopifnot(all(file.info(required_outputs)$size > 0))

scope <- fread(file.path(table_dir, "project_sample_group_scope_33.csv"))
removed <- fread(file.path(table_dir, "removed_sample_groups_4.csv"))
protein_scope <- fread(file.path(table_dir, "kla_ddr_protein_scope.csv"))
long <- fread(file.path(table_dir, "protein_raw_go_term_long.csv"))
term_summary <- fread(file.path(table_dir, "raw_go_term_summary.csv"))
pattern_summary <- fread(file.path(table_dir, "raw_go_pattern_summary.csv"))
binary <- fread(file.path(table_dir, "protein_raw_go_term_binary_matrix.csv"))
coordinates <- fread(file.path(table_dir, "umap_coordinates_fixed.csv"))
parameters <- fread(file.path(table_dir, "umap_parameters.csv"))

stopifnot(nrow(scope) == 33L)
stopifnot(uniqueN(scope, by = c("PXD", "SampleGroup")) == 33L)
stopifnot(identical(
  scope[, .N, by = Category][order(Category)],
  data.table(
    Category = c("cancer_cells", "cancer_tissue", "normal_cells", "normal_tissue"),
    N = c(13L, 2L, 9L, 9L)
  )
))

stopifnot(nrow(removed) == 4L)
stopifnot(all(removed$PresentIn37GroupSource))
removed_keys <- paste(removed$PXD, removed$SampleGroup, sep = "\r")
scope_keys <- paste(scope$PXD, scope$SampleGroup, sep = "\r")
stopifnot(!any(removed_keys %chin% scope_keys))

stopifnot(nrow(protein_scope) == 507L)
stopifnot(uniqueN(protein_scope$BaseAccession) == 507L)
stopifnot(all(protein_scope$NumberOfRawGOTerms >= 1L))

stopifnot(nrow(long) == 1029L)
stopifnot(uniqueN(long$BaseAccession) == 507L)
stopifnot(uniqueN(long$GO_TERM) == 66L)
stopifnot(uniqueN(long, by = c("BaseAccession", "GO_TERM")) == 1029L)
stopifnot(nrow(term_summary) == 66L)
stopifnot(nrow(pattern_summary) == 205L)
stopifnot(pattern_summary$PatternID[[1L]] == "PATTERN_001")
stopifnot(pattern_summary$RawGOTermPattern[[1L]] == "GO:0006974")
stopifnot(pattern_summary$NumberOfProteins[[1L]] == 102L)

stopifnot(nrow(binary) == 507L)
stopifnot(ncol(binary) == 67L)
stopifnot(identical(binary$BaseAccession, sort(binary$BaseAccession)))
binary_matrix <- as.matrix(binary[, -1L])
stopifnot(all(binary_matrix %in% c(0, 1)))
stopifnot(all(rowSums(binary_matrix) >= 1L))
stopifnot(sum(binary_matrix) == 1029L)
stopifnot(setequal(names(binary)[-1L], long$GO_TERM))

stopifnot(identical(names(coordinates), c("BaseAccession", "UMAP_1", "UMAP_2")))
stopifnot(nrow(coordinates) == 507L)
stopifnot(uniqueN(coordinates$BaseAccession) == 507L)
stopifnot(setequal(coordinates$BaseAccession, protein_scope$BaseAccession))
stopifnot(all(is.finite(coordinates$UMAP_1)))
stopifnot(all(is.finite(coordinates$UMAP_2)))

parameter_value <- setNames(parameters$Value, parameters$Parameter)
stopifnot(parameter_value[["ProteinAnalysisKey"]] == "isoform-stripped UniProt BaseAccession")
stopifnot(parameter_value[["UMAPInput"]] == "protein x raw GO term binary matrix only")
stopifnot(parameter_value[["GOCategoryAggregation"]] == "none")
stopifnot(parameter_value[["GeneSymbolUsedAsAnalysisKey"]] == "FALSE")
stopifnot(parameter_value[["DistanceMetric"]] == "cosine")
stopifnot(parameter_value[["Initialization"]] == "random")
stopifnot(parameter_value[["RandomSeed"]] == "25")

analysis_output_names <- c(
  names(binary),
  names(coordinates),
  names(long),
  names(protein_scope)
)
stopifnot(!any(grepl("GeneSymbol|Symbol", analysis_output_names, ignore.case = TRUE)))
stopifnot(!any(names(binary) %chin% c("HR", "NHEJ", "BER", "NER", "MMR", "TLS", "DRR", "CP", "Other")))

message("PASS: 33-group raw-GO Kla-DDR UMAP outputs are internally consistent.")
