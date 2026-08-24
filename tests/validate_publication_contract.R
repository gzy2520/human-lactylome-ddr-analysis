#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
})

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args)) normalizePath(args[[1]], mustWork = TRUE) else normalizePath(".", mustWork = TRUE)
input_dir <- file.path(project_root, "data", "publication_input")
figure_dir <- file.path(project_root, "results", "figures")
supplementary_dir <- file.path(project_root, "results", "supplementary")

assert <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}
is_true <- function(x) as.character(x) %in% c("TRUE", "True", "true", "T", "1")
md5_file <- function(path) digest::digest(file = path, algo = "md5", serialize = FALSE)

manifest <- fread(file.path(input_dir, "INPUT_MANIFEST.csv"))
assert(identical(names(manifest), c("File", "Bytes", "MD5")), "The frozen input manifest has an unexpected schema.")
actual_input_files <- list.files(input_dir, recursive = TRUE, full.names = FALSE, no.. = TRUE)
assert(
  setequal(actual_input_files, c("INPUT_MANIFEST.csv", manifest$File)),
  "The frozen publication input contains an unmanifested or missing file."
)
for (index in seq_len(nrow(manifest))) {
  filename <- manifest$File[[index]]
  path <- file.path(input_dir, filename)
  assert(file.exists(path), paste("Frozen input is missing:", filename))
  assert(file.info(path)$size == manifest$Bytes[[index]], paste("Frozen input byte count changed:", filename))
  assert(md5_file(path) == manifest$MD5[[index]], paste("Frozen input checksum changed:", filename))
}

groups <- fread(file.path(input_dir, "group_summary_30.csv"))
assert(nrow(groups) == 30L, "The project must begin with exactly 30 publication groups.")
assert(!anyDuplicated(groups[, .(PXD, SampleGroup)]), "Publication groups must be unique.")
expected_categories <- data.table(
  Category = c("cancer_cells", "cancer_tissue", "normal_cells", "normal_tissue"),
  N = c(12L, 2L, 7L, 9L)
)
observed_categories <- groups[, .N, by = Category][order(Category)]
assert(identical(observed_categories, expected_categories), "Publication category counts must remain 12/2/7/9.")

group_keys <- paste(groups$PXD, groups$SampleGroup, sep = "__")
for (filename in c("kla_protein_membership_30.csv", "reference_protein_membership_30.csv", "regulator_kla_percentiles_30.csv", "regulator_reference_percentiles_30.csv")) {
  data <- fread(file.path(input_dir, filename))
  keys <- unique(paste(data$PXD, data$SampleGroup, sep = "__"))
  assert(setequal(keys, group_keys), paste("Frozen input contains a non-final group:", filename))
}

kla_ddr <- fread(file.path(input_dir, "venn_kla_ddr.csv"))
assert(nrow(kla_ddr) == 399L && uniqueN(kla_ddr$BaseAccession) == 399L, "The final Kla-DDR union must contain 399 BaseAccessions.")
set_sizes <- c(
  normal_tissue = sum(is_true(kla_ddr$In_normal_tissue)),
  cancer_tissue = sum(is_true(kla_ddr$In_cancer_tissue)),
  cancer_cells = sum(is_true(kla_ddr$In_cancer_cells)),
  normal_cells = sum(is_true(kla_ddr$In_normal_cells))
)
assert(identical(as.integer(set_sizes), c(183L, 178L, 381L, 292L)), "Final pathway panels must contain 183/178/381/292 proteins.")

pathways <- c("BER", "NER", "MMR", "FA", "HR", "AEJ", "NHEJ")
s4_input <- file.path(input_dir, "Supplementary_Table_S4_Pathway_Protein_Ranking.xlsx")
s4_specs <- list(
  list(key = "normal_tissue", sheet = "NonTumorTissues", rows = 183L),
  list(key = "cancer_tissue", sheet = "TumorTissues", rows = 178L),
  list(key = "cancer_cells", sheet = "CancerCellLines", rows = 381L),
  list(key = "normal_cells", sheet = "NormalCellLines", rows = 292L)
)
for (spec in s4_specs) {
  panel <- as.data.table(read_excel(s4_input, sheet = spec$sheet))
  panel <- panel[!is.na(BaseAccession) & nzchar(trimws(BaseAccession))]
  panel[, BaseAccession := trimws(as.character(BaseAccession))]
  assert(nrow(panel) == spec$rows && !anyDuplicated(panel$BaseAccession), paste("Frozen S4 has an invalid panel:", spec$sheet))
  assert(all(c("BaseAccession", "SignedScore", pathways) %in% names(panel)), paste("Frozen S4 is missing columns:", spec$sheet))
  expected_ids <- kla_ddr$BaseAccession[is_true(kla_ddr[[paste0("In_", spec$key)]])]
  assert(setequal(panel$BaseAccession, expected_ids), paste("Frozen S4 membership changed:", spec$sheet))
  states <- as.matrix(panel[, ..pathways])
  storage.mode(states) <- "numeric"
  assert(all(states %in% c(-1, 0, 1)), paste("Frozen S4 states changed:", spec$sheet))
  expected_scores <- as.numeric(states %*% seq_along(pathways))
  assert(identical(as.numeric(panel$SignedScore), expected_scores), paste("Frozen S4 scores changed:", spec$sheet))
  assert(identical(order(panel$SignedScore, panel$BaseAccession), seq_len(nrow(panel))), paste("Frozen S4 order changed:", spec$sheet))
}
s5_input <- file.path(input_dir, "Supplementary_Table_S5_Lactylation_Regulators.xlsx")
s5_annotations <- as.data.table(read_excel(s5_input, sheet = "Regulator_Annotations"))
assert(all(c("Role", "GeneSymbol", "BaseAccession") %in% names(s5_annotations)), "Frozen S5 lacks accession-keyed regulator annotations.")
s5_role_map <- unique(s5_annotations[!is.na(BaseAccession) & nzchar(trimws(BaseAccession)), .(Role, BaseAccession)])
assert(nrow(s5_role_map) == 49L, "Frozen S5 must retain its 49 unique role/BaseAccession display mappings.")

figure_stems <- c(
  "Figure_1_DDR_fraction",
  "Figure_2a_whole_proteome_DDR_Venn",
  "Figure_2b_Kla_DDR_Venn",
  "Figure_2c_DDR_pathway_matrices_tissues",
  "Figure_2d_DDR_pathway_summary_tumor_tissue",
  "Figure_2e_DDR_pathway_summary_non_tumor_tissue",
  "Figure_3a_reference_regulator_percentiles",
  "Figure_3b_Kla_regulator_percentiles",
  "Supplementary_Figure_S1a_whole_proteome_Venn",
  "Supplementary_Figure_S1b_Kla_proteome_Venn",
  "Supplementary_Figure_S2a_DDR_pathway_matrices_cell_lines",
  "Supplementary_Figure_S2b_DDR_pathway_summary_cancer_cell_lines",
  "Supplementary_Figure_S2c_DDR_pathway_summary_normal_cell_lines"
)
expected_figures <- as.vector(outer(figure_stems, c(".png", ".pdf"), paste0))
actual_figures <- list.files(figure_dir, recursive = TRUE, full.names = FALSE)
assert(setequal(actual_figures, expected_figures), "Results contain a figure not described by the final manuscript, or a required figure is missing.")
assert(all(file.info(file.path(figure_dir, expected_figures))$size > 0L), "A required figure is empty.")

expected_workbooks <- c(
  "Supplementary_Table_S1_Kla_Data.xlsx",
  "Supplementary_Table_S2_Reference_Data.xlsx",
  "Supplementary_Table_S3_Human_DDR_GO_Annotations.xlsx",
  "Supplementary_Table_S4_Pathway_Protein_Ranking.xlsx",
  "Supplementary_Table_S5_Lactylation_Regulators.xlsx",
  "Supplementary_Table_S6_Venn_Membership.xlsx"
)
actual_workbooks <- list.files(supplementary_dir, recursive = TRUE, full.names = FALSE)
assert(setequal(actual_workbooks, expected_workbooks), "Results contain a supplementary workbook outside Tables S1-S6, or a required table is missing.")
assert(all(file.info(file.path(supplementary_dir, expected_workbooks))$size > 0L), "A supplementary workbook is empty.")

expected_sheets <- list(
  Supplementary_Table_S1_Kla_Data.xlsx = c("Group_Summary", "Kla_Protein_Membership", "Kla_DDR_Membership"),
  Supplementary_Table_S2_Reference_Data.xlsx = c("Reference_Group_Summary", "Reference_Protein_Membership", "Reference_DDR_Membership"),
  Supplementary_Table_S3_Human_DDR_GO_Annotations.xlsx = "Human_DDR_GO_Annotations",
  Supplementary_Table_S4_Pathway_Protein_Ranking.xlsx = c("NonTumorTissues", "TumorTissues", "CancerCellLines", "NormalCellLines"),
  Supplementary_Table_S5_Lactylation_Regulators.xlsx = c("Regulator_Annotations", "Regulator_ID_Mapping"),
  Supplementary_Table_S6_Venn_Membership.xlsx = c("AllKla_Members", "KlaDDR_Members", "Reference_Members", "ReferenceDDR_Members", "Set_Counts", "Region_Counts")
)
for (filename in names(expected_sheets)) {
  assert(identical(excel_sheets(file.path(supplementary_dir, filename)), expected_sheets[[filename]]), paste("Unexpected workbook sheet layout:", filename))
}

s1_path <- file.path(supplementary_dir, "Supplementary_Table_S1_Kla_Data.xlsx")
s2_path <- file.path(supplementary_dir, "Supplementary_Table_S2_Reference_Data.xlsx")
s3_path <- file.path(supplementary_dir, "Supplementary_Table_S3_Human_DDR_GO_Annotations.xlsx")
s4_path <- file.path(supplementary_dir, "Supplementary_Table_S4_Pathway_Protein_Ranking.xlsx")
s5_path <- file.path(supplementary_dir, "Supplementary_Table_S5_Lactylation_Regulators.xlsx")
kla_membership <- fread(file.path(input_dir, "kla_protein_membership_30.csv"))
reference_membership <- fread(file.path(input_dir, "reference_protein_membership_30.csv"))
assert(nrow(read_excel(s1_path, sheet = "Group_Summary")) == 30L, "S1 Group_Summary must have 30 rows.")
assert(nrow(read_excel(s1_path, sheet = "Kla_Protein_Membership")) == nrow(kla_membership), "S1 Kla membership row count changed.")
assert(nrow(read_excel(s1_path, sheet = "Kla_DDR_Membership")) == sum(is_true(kla_membership$IsDdr)), "S1 Kla-DDR membership row count changed.")
assert(nrow(read_excel(s2_path, sheet = "Reference_Group_Summary")) == 30L, "S2 Reference_Group_Summary must have 30 rows.")
assert(nrow(read_excel(s2_path, sheet = "Reference_Protein_Membership")) == nrow(reference_membership), "S2 reference membership row count changed.")
assert(nrow(read_excel(s2_path, sheet = "Reference_DDR_Membership")) == sum(is_true(reference_membership$IsDdr)), "S2 reference DDR membership row count changed.")
assert(nrow(read_excel(s3_path, sheet = "Human_DDR_GO_Annotations", col_types = "text")) == nrow(fread(file.path(input_dir, "human_ddr_go_annotations.tsv"), sep = "\t", quote = "")), "S3 GO annotation row count changed.")
assert(md5_file(s4_path) == md5_file(s4_input), "S4 must be copied unchanged from its frozen release asset.")
assert(md5_file(s5_path) == md5_file(s5_input), "S5 must be copied unchanged from its frozen release asset.")

s6_path <- file.path(supplementary_dir, "Supplementary_Table_S6_Venn_Membership.xlsx")
s6_regions <- as.data.table(read_excel(s6_path, sheet = "Region_Counts"))
assert(nrow(s6_regions) == 60L, "S6 must retain all 15 Venn regions for each of its four analyses.")
assert(all(s6_regions[, .N, by = Analysis]$N == 15L) && uniqueN(s6_regions[, .(Analysis, Region)]) == 60L, "Each S6 analysis must contain exactly 15 unique Venn regions, including zero-count regions.")
venn_sources <- c(AllKla = "venn_all_kla.csv", KlaDDR = "venn_kla_ddr.csv", Reference = "venn_reference.csv", ReferenceDDR = "venn_reference_ddr.csv")
for (analysis_name in names(venn_sources)) {
  membership <- fread(file.path(input_dir, venn_sources[[analysis_name]]))
  output_membership <- read_excel(s6_path, sheet = paste0(analysis_name, "_Members"))
  expected_regions <- membership[, .N, by = Region]
  observed_regions <- s6_regions[Analysis == analysis_name]
  reconstructed <- expected_regions$N[match(observed_regions$Region, expected_regions$Region)]
  reconstructed[is.na(reconstructed)] <- 0L
  assert(nrow(output_membership) == nrow(membership), paste("S6 membership row count changed for", analysis_name))
  assert(identical(as.integer(observed_regions$ProteinCount), as.integer(reconstructed)), paste("S6 region counts changed for", analysis_name))
}

message("PASS: exact 30-group scope, 399 BaseAccessions, four signed pathway panels, manuscript figures, and Tables S1-S6 only.")
