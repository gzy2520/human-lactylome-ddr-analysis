#!/usr/bin/env Rscript

# Build only Supplementary Tables S1-S6 named in DDR_Kla_manuscript_V3.docx.
# Inputs are frozen to the final 30-group publication boundary.  All set
# operations use BaseAccession; symbols and protein names are display fields.

suppressPackageStartupMessages({
  library(data.table)
  library(openxlsx2)
})

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args)) normalizePath(args[[1L]], mustWork = TRUE) else normalizePath(".", mustWork = TRUE)
input_dir <- file.path(project_root, "data", "publication_input")
output_dir <- file.path(project_root, "results", "supplementary")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

assert <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

input_path <- function(filename) file.path(input_dir, filename)
md5_file <- function(path) digest::digest(file = path, algo = "md5", serialize = FALSE)
is_true <- function(x) as.character(x) %in% c("TRUE", "True", "true", "T", "1")

require_columns <- function(data, columns, label) {
  missing <- setdiff(columns, names(data))
  assert(!length(missing), paste(label, "is missing column(s):", paste(missing, collapse = ", ")))
}

column_letter <- function(index) {
  out <- ""
  while (index > 0L) {
    index <- index - 1L
    out <- paste0(LETTERS[index %% 26L + 1L], out)
    index <- index %/% 26L
  }
  out
}

column_widths <- function(headers) {
  widths <- pmin(pmax(nchar(headers, type = "width") + 3L, 12L), 28L)
  wide <- grepl("Source|File|Material|Evidence|ProteinName|Reference|Note|Description", headers, ignore.case = TRUE)
  widths[wide] <- pmax(widths[wide], 36L)
  as.numeric(pmin(widths, 42L))
}

add_sheet <- function(workbook, sheet_name, data, table_name) {
  data <- as.data.table(data)
  assert(nrow(data) > 0L && ncol(data) > 0L, paste("No data available for", sheet_name))
  workbook <- wb_add_worksheet(workbook, sheet = sheet_name, grid_lines = FALSE)
  workbook <- wb_add_data_table(
    workbook,
    sheet = sheet_name,
    x = as.data.frame(data),
    table_name = table_name,
    table_style = "TableStyleLight1",
    banded_rows = FALSE,
    with_filter = TRUE,
    inline_strings = FALSE
  )
  last_column <- column_letter(ncol(data))
  workbook <- wb_add_font(workbook, sheet = sheet_name, dims = paste0("A1:", last_column, "1"), name = "Arial", size = 10, bold = TRUE)
  workbook <- wb_add_cell_style(workbook, sheet = sheet_name, dims = paste0("A1:", last_column, "1"), horizontal = "center", vertical = "center", wrap_text = TRUE)
  workbook <- wb_freeze_pane(workbook, sheet = sheet_name, first_row = TRUE)
  workbook <- wb_set_col_widths(workbook, sheet = sheet_name, cols = seq_len(ncol(data)), widths = column_widths(names(data)))
  workbook
}

write_workbook <- function(sheets, filename) {
  workbook <- wb_workbook(creator = "klaPublicationAnalysis")
  for (sheet in sheets) {
    workbook <- add_sheet(workbook, sheet$name, sheet$data, sheet$table_name)
  }
  wb_save(workbook, file.path(output_dir, filename), overwrite = TRUE)
}

copy_frozen_workbook <- function(filename) {
  source <- input_path(filename)
  target <- file.path(output_dir, filename)
  assert(file.exists(source), paste("Frozen supplementary workbook is missing:", filename))
  assert(isTRUE(file.copy(source, target, overwrite = TRUE)), paste("Could not copy frozen supplementary workbook:", filename))
  assert(md5_file(source) == md5_file(target), paste("Frozen supplementary workbook changed while copying:", filename))
}

groups <- fread(input_path("group_summary_30.csv"))
assert(nrow(groups) == 30L, "Supplementary Tables must begin with exactly 30 Kla groups.")
kla_membership <- fread(input_path("kla_protein_membership_30.csv"))
reference_membership <- fread(input_path("reference_protein_membership_30.csv"))
require_columns(kla_membership, c("PXD", "SampleGroup", "BaseAccession", "IsDdr"), "Kla membership")
require_columns(reference_membership, c("PXD", "SampleGroup", "SourceProteinID", "IsDdr"), "Reference membership")

# S1: final Kla groups and their exact protein memberships.
write_workbook(
  list(
    list(name = "Group_Summary", data = groups, table_name = "S1GroupSummary"),
    list(name = "Kla_Protein_Membership", data = kla_membership, table_name = "S1KlaMembership"),
    list(name = "Kla_DDR_Membership", data = kla_membership[is_true(IsDdr)], table_name = "S1KlaDDRMembership")
  ),
  "Supplementary_Table_S1_Kla_Data.xlsx"
)

# S2 mirrors S1's summary-first layout while retaining sample-level
# whole-proteome records and the exact DDR subset.
summary_columns <- c(
  "RowOrder", "Category", "CategoryEn", "PXD", "SampleGroup", "BiologicalMaterial",
  "ReferencePXD", "ReferenceEvidenceFile", "ReferenceProteinCount", "ReferenceDdrProteinCount",
  "ReferenceDdrFraction", "ReferenceMatchNote", "MatchMode"
)
require_columns(groups, summary_columns, "Group summary")
write_workbook(
  list(
    list(name = "Reference_Group_Summary", data = groups[, ..summary_columns], table_name = "S2ReferenceSummary"),
    list(name = "Reference_Protein_Membership", data = reference_membership, table_name = "S2ReferenceMembership"),
    list(name = "Reference_DDR_Membership", data = reference_membership[is_true(IsDdr)], table_name = "S2ReferenceDDRMembership")
  ),
  "Supplementary_Table_S2_Reference_Data.xlsx"
)

# S3: frozen human DDR GO annotation records.
go_annotations <- fread(input_path("human_ddr_go_annotations.tsv"), sep = "\t", quote = "", encoding = "UTF-8")
write_workbook(
  list(list(name = "Human_DDR_GO_Annotations", data = go_annotations, table_name = "S3HumanDDRGO")),
  "Supplementary_Table_S3_Human_DDR_GO_Annotations.xlsx"
)

# S4 and S5 are author-approved release assets. They are retained byte-for-byte
# and are not recalculated or reformatted by this workflow.
copy_frozen_workbook("Supplementary_Table_S4_Pathway_Protein_Ranking.xlsx")
copy_frozen_workbook("Supplementary_Table_S5_Lactylation_Regulators.xlsx")

# S6: exact four-category memberships, set sizes and all fifteen possible
# Venn regions (including regions whose frozen count is zero).
venn_categories <- c("normal_tissue", "cancer_tissue", "cancer_cells", "normal_cells")
all_venn_regions <- function() {
  vapply(seq_len(15L), function(mask) {
    members <- venn_categories[vapply(seq_along(venn_categories), function(index) bitwAnd(mask, bitwShiftL(1L, index - 1L)) != 0L, logical(1))]
    if (length(members) == length(venn_categories)) "all_four" else paste0(paste(members, collapse = "_and_"), "_only")
  }, character(1))
}
venn_specs <- list(
  list(label = "AllKla", file = "venn_all_kla.csv"),
  list(label = "KlaDDR", file = "venn_kla_ddr.csv"),
  list(label = "Reference", file = "venn_reference.csv"),
  list(label = "ReferenceDDR", file = "venn_reference_ddr.csv")
)
s6_sheets <- list()
set_counts <- list()
region_counts <- list()
for (spec in venn_specs) {
  membership <- fread(input_path(spec$file))
  require_columns(membership, c("BaseAccession", "Region", paste0("In_", venn_categories)), paste("S6", spec$label))
  s6_sheets[[length(s6_sheets) + 1L]] <- list(name = paste0(spec$label, "_Members"), data = membership, table_name = paste0("S6", spec$label, "Members"))
  set_counts[[length(set_counts) + 1L]] <- data.table(
    Analysis = spec$label,
    Category = venn_categories,
    ProteinCount = vapply(venn_categories, function(category) sum(is_true(membership[[paste0("In_", category)]])), integer(1))
  )
  observed <- membership[, .N, by = Region]
  regions <- all_venn_regions()
  region_counts[[length(region_counts) + 1L]] <- data.table(
    Analysis = spec$label,
    Region = regions,
    ProteinCount = fifelse(is.na(observed$N[match(regions, observed$Region)]), 0L, observed$N[match(regions, observed$Region)])
  )
}
s6_sheets[[length(s6_sheets) + 1L]] <- list(name = "Set_Counts", data = rbindlist(set_counts), table_name = "S6SetCounts")
s6_sheets[[length(s6_sheets) + 1L]] <- list(name = "Region_Counts", data = rbindlist(region_counts), table_name = "S6RegionCounts")
write_workbook(s6_sheets, "Supplementary_Table_S6_Venn_Membership.xlsx")

message("PASS: built Supplementary Tables S1-S3/S6 and copied frozen S4/S5 unchanged.")
