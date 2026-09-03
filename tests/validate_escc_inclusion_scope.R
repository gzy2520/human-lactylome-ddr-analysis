#!/usr/bin/env Rscript

# Validate the isolated PXD064038 inclusion scope and its rendered outputs.
# All biological membership checks use stable UniProt BaseAccessions.

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
})

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1L) {
  normalizePath(args[[1L]], mustWork = TRUE)
} else {
  normalizePath(".", mustWork = TRUE)
}

stop_if <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

scope_tag <- Sys.getenv(
  "KLA_ESCC_SCOPE_TAG",
  unset = "escc_inclusion_20260903_pxd065830_tumor_reference"
)
scope_dir <- file.path(project_root, "data", "candidate", scope_tag)
publication_dir <- file.path(scope_dir, "publication_input")
candidate_dir <- file.path(scope_dir, "candidate_input")
formal_output_dir <- file.path(project_root, "results", scope_tag, "formal_figures")
candidate_output_dir <- file.path(project_root, "results", "candidate", scope_tag)
legacy_output_dir <- file.path(candidate_output_dir, "legacy_layout")
pathway_output_dir <- file.path(legacy_output_dir, "pathway_summary_by_pathway")

stop_if(dir.exists(scope_dir), "The dated ESCC inclusion scope is missing.")

audit <- fread(file.path(scope_dir, "escc_inclusion_audit.csv"))
audit_values <- stats::setNames(audit$Value, audit$Item)
expected_audit <- c(
  selected_group = "MEC and NEC ESCC groups",
  excluded_escc_related_dataset_1 = "PXD048995 (KYSE30 histone-focused; excluded)",
  excluded_escc_related_dataset_2 = "PXD063945 (ESCC neoadjuvant; excluded)",
  escc_source_pxd = "PXD064038",
  escc_sample_count = "6",
  escc_kla_union_count = "1239",
  escc_kla_ddr_count = "92",
  strict_reference_status = "available; independent ESCC tumor ordinary-proteome union by BaseAccession",
  reference_pxd = "PXD065830",
  reference_tumor_sample_count = "94",
  reference_protein_union_count = "8083",
  reference_ddr_count = "420",
  reference_ddr_fraction = "0.0519609056",
  expanded_group_count = "31",
  expanded_tumor_S4_count = "192"
)
stop_if(all(unname(expected_audit) == unname(audit_values[names(expected_audit)])),
  "The ESCC inclusion audit does not match the recorded scope.")

groups <- fread(file.path(publication_dir, "group_summary_30.csv"))
stop_if(nrow(groups) == 31L, "Expanded publication input does not contain 31 groups.")
stop_if(!anyDuplicated(groups[, .(PXD, SampleGroup)]), "Expanded publication groups are not unique.")
stop_if(any(groups$PXD == "PXD064038" & groups$SampleGroup == "MEC and NEC ESCC groups"),
  "The selected PXD064038 group is missing.")
stop_if(!any(groups$PXD %in% c("PXD048995", "PXD063945")),
  "One of the two excluded ESCC-related datasets entered the expanded publication scope.")

category_counts <- groups[, .N, by = Category][match(
  c("normal_tissue", "cancer_tissue", "cancer_cells", "normal_cells"), Category
)]$N
stop_if(identical(as.integer(category_counts), c(9L, 3L, 12L, 7L)),
  "Expanded category counts are not 9/3/12/7.")

selected_group <- groups[PXD == "PXD064038" & SampleGroup == "MEC and NEC ESCC groups"]
stop_if(nrow(selected_group) == 1L, "The selected ESCC group is not unique.")
stop_if(selected_group$Category == "cancer_tissue", "The selected ESCC group is not classified as tumor tissue.")
stop_if(selected_group$KlaProteinCount == 1239L && selected_group$KlaDdrProteinCount == 92L,
  "PXD064038 group-level Kla counts changed.")
stop_if(
  selected_group$PairedAnalysisIncluded == TRUE &&
    selected_group$ReferencePXD == "PXD065830" &&
    selected_group$ReferenceProteinCount == 8083L &&
    selected_group$ReferenceDdrProteinCount == 420L &&
    abs(selected_group$ReferenceDdrFraction - 420 / 8083) < 1e-12 &&
    abs(selected_group$DdrFractionPercentagePointDifference - (92 / 1239 - 420 / 8083)) < 1e-12,
  "PXD064038 was not assigned the 94-sample ESCC tumor reference correctly.")
stop_if(selected_group$MatchMode == "BaseAccession_external_ESCC_tumor",
  "PXD064038 reference mode changed unexpectedly.")

membership <- fread(file.path(publication_dir, "kla_protein_membership_30.csv"))
selected_membership <- membership[PXD == "PXD064038" & SampleGroup == "MEC and NEC ESCC groups"]
stop_if(nrow(selected_membership) == 1239L && uniqueN(selected_membership$BaseAccession) == 1239L,
  "PXD064038 Kla membership does not contain 1239 unique BaseAccessions.")
stop_if(sum(selected_membership$IsDdr == TRUE) == 92L, "PXD064038 Kla-DDR membership does not contain 92 BaseAccessions.")

reference_membership <- fread(file.path(publication_dir, "reference_protein_membership_30.csv"))
selected_reference_membership <- reference_membership[
  PXD == "PXD064038" & SampleGroup == "MEC and NEC ESCC groups"
]
stop_if(
  nrow(selected_reference_membership) == 8083L &&
    uniqueN(selected_reference_membership$SourceProteinID) == 8083L &&
    uniqueN(selected_reference_membership$MappedBaseAccessions) == 8083L &&
    all(selected_reference_membership$ReferencePXD == "PXD065830") &&
    all(selected_reference_membership$IdentifierType == "UniProtKB") &&
    sum(selected_reference_membership$IsDdr == TRUE) == 420L,
  "PXD064038 reference membership does not contain the 8083-protein ESCC tumor union.")

all_venn <- fread(file.path(publication_dir, "venn_all_kla.csv"))
ddr_venn <- fread(file.path(publication_dir, "venn_kla_ddr.csv"))
stop_if(uniqueN(all_venn$BaseAccession) == nrow(all_venn) && nrow(all_venn) == 5814L,
  "Expanded all-Kla union does not contain 5814 unique BaseAccessions.")
stop_if(uniqueN(ddr_venn$BaseAccession) == nrow(ddr_venn) && nrow(ddr_venn) == 401L,
  "Expanded Kla-DDR union does not contain 401 unique BaseAccessions.")
stop_if(sum(ddr_venn$In_cancer_tissue == TRUE) == 192L, "Expanded tumor Kla-DDR membership is not 192.")
old_all_venn <- fread(file.path(project_root, "data", "publication_input", "venn_all_kla.csv"))
old_ddr_venn <- fread(file.path(project_root, "data", "publication_input", "venn_kla_ddr.csv"))
stop_if(length(setdiff(all_venn$BaseAccession, old_all_venn$BaseAccession)) == 72L,
  "The expanded all-Kla union did not add exactly 72 BaseAccessions.")
stop_if(length(setdiff(ddr_venn$BaseAccession, old_ddr_venn$BaseAccession)) == 2L,
  "The expanded Kla-DDR union did not add exactly 2 BaseAccessions.")
reference_venn <- fread(file.path(publication_dir, "venn_reference.csv"))
reference_ddr_venn <- fread(file.path(publication_dir, "venn_reference_ddr.csv"))
stop_if(uniqueN(reference_venn$BaseAccession) == nrow(reference_venn) && nrow(reference_venn) == 24397L,
  "Expanded reference union does not contain 24397 unique BaseAccessions.")
stop_if(uniqueN(reference_ddr_venn$BaseAccession) == nrow(reference_ddr_venn) && nrow(reference_ddr_venn) == 836L,
  "Expanded reference-DDR union does not contain 836 unique BaseAccessions.")

sample_values <- fread(file.path(candidate_dir, "sample_boxplot_values.csv"))
sample_summary <- fread(file.path(scope_dir, "escc_sample_summary.csv"))
expected_samples <- data.table(
  SampleID = c("MEC_1", "MEC_2", "MEC_3", "NEC_1", "NEC_2", "NEC_3"),
  KlaProteinCount = c(820L, 794L, 838L, 464L, 483L, 573L),
  KlaDdrProteinCount = c(65L, 61L, 79L, 26L, 27L, 33L)
)
selected_samples <- sample_values[PXD == "PXD064038", .(
  SampleID, KlaProteinCount, KlaDdrProteinCount, SampleClass, SourceMode, ReferenceFraction
)]
setorder(selected_samples, SampleID)
setorder(expected_samples, SampleID)
stop_if(nrow(selected_samples) == 6L, "PXD064038 does not contribute six sample observations.")
stop_if(all(selected_samples[, .(SampleID, KlaProteinCount, KlaDdrProteinCount)] == expected_samples),
  "PXD064038 sample-level counts changed.")
stop_if(all(selected_samples$SourceMode == "deposited_sample_table") &&
          all(abs(selected_samples$ReferenceFraction - 420 / 8083 * 100) < 1e-12),
  "PXD064038 sample provenance or reference status changed.")
stop_if(setequal(sample_summary$SampleID, expected_samples$SampleID), "ESCC sample summary is incomplete.")

reconciliation <- fread(file.path(candidate_dir, "sample_boxplot_reconciliation.csv"))
stop_if(nrow(reconciliation) == 31L && !anyDuplicated(reconciliation[, .(PXD, SampleGroup)]),
  "Expanded sample reconciliation does not cover 31 unique groups.")
stop_if(reconciliation[PXD == "PXD064038" & SampleGroup == "MEC and NEC ESCC groups", GroupUnionStatus] == "PASS",
  "PXD064038 sample observations do not reconcile to the expanded group union.")

figure1_values <- fread(file.path(candidate_dir, "figure1_sample_boxplot_values.csv"))
pathway_values <- fread(file.path(candidate_dir, "figure1_pathway_summary_sample_boxplot_values.csv"))
stop_if(nrow(figure1_values) == 310L, "Expanded Figure 1 sample input does not contain 310 observations.")
stop_if(nrow(figure1_values[Dataset == "Lactylome (Kla)"]) == 98L &&
          nrow(figure1_values[Dataset == "Whole proteome"]) == 212L,
  "Expanded Figure 1 dataset counts changed.")
selected_sample_figure1 <- figure1_values[PXD == "PXD064038"]
stop_if(nrow(selected_sample_figure1[Dataset == "Lactylome (Kla)"]) == 6L &&
          nrow(selected_sample_figure1[Dataset == "Whole proteome"]) == 94L &&
          uniqueN(selected_sample_figure1[Dataset == "Whole proteome", SampleID]) == 94L,
  "PXD065830 ESCC-T observations were not retained as 94 individual Figure 1 points.")
stop_if(nrow(pathway_values) == 686L && nrow(pathway_values[PXD == "PXD064038"]) == 42L,
  "Expanded pathway-summary input does not contain the six-by-seven ESCC panel.")

dataset_figure1_values <- fread(file.path(candidate_dir, "figure1_dataset_boxplot_values.csv"))
dataset_pathway_values <- fread(file.path(candidate_dir, "pathway_summary_dataset_boxplot_values.csv"))
stop_if(nrow(dataset_figure1_values) == 62L &&
          uniqueN(dataset_figure1_values$DatasetPointID) == 62L,
  "Expanded dataset-level Figure 1 input does not contain 31 points per modality.")
stop_if(nrow(dataset_pathway_values) == 434L &&
          uniqueN(dataset_pathway_values$DatasetPointID) == 62L &&
          uniqueN(dataset_pathway_values[, .(DatasetPointID, Pathway)]) == 434L,
  "Expanded dataset-level pathway input does not contain 31 points x 2 modalities x 7 pathways.")
selected_dataset_figure1 <- dataset_figure1_values[
  PXD == "PXD064038" & SampleGroup == "MEC and NEC ESCC groups"
]
stop_if(nrow(selected_dataset_figure1) == 2L &&
          selected_dataset_figure1[Dataset == "Lactylome (Kla)", DdrProteinCount] == 92L &&
          selected_dataset_figure1[Dataset == "Whole proteome", DdrProteinCount] == 420L &&
          abs(selected_dataset_figure1[Dataset == "Lactylome (Kla)", DdrFractionPercentage] - 92 / 1239 * 100) < 1e-12 &&
          abs(selected_dataset_figure1[Dataset == "Whole proteome", DdrFractionPercentage] - 420 / 8083 * 100) < 1e-12,
  "Expanded dataset-level Figure 1 values do not retain the selected ESCC source fractions.")

s4_path <- file.path(publication_dir, "Supplementary_Table_S4_Pathway_Protein_Ranking.xlsx")
s4_expected <- c(NonTumorTissues = 183L, TumorTissues = 192L, CancerCellLines = 381L, NormalCellLines = 292L)
sheet_category <- c(
  NonTumorTissues = "normal_tissue",
  TumorTissues = "cancer_tissue",
  CancerCellLines = "cancer_cells",
  NormalCellLines = "normal_cells"
)
for (sheet_name in names(s4_expected)) {
  panel <- as.data.table(read_excel(s4_path, sheet = sheet_name))
  panel <- panel[!is.na(BaseAccession) & nzchar(trimws(BaseAccession))]
  expected_ids <- ddr_venn[get(paste0("In_", sheet_category[[sheet_name]])) == TRUE, BaseAccession]
  stop_if(nrow(panel) == s4_expected[[sheet_name]] && setequal(panel$BaseAccession, expected_ids),
    paste("S4 panel does not match expanded Kla-DDR membership for", sheet_name))
}

manifest_path <- file.path(publication_dir, "INPUT_MANIFEST.csv")
manifest <- fread(manifest_path)
actual_files <- list.files(publication_dir, recursive = TRUE, full.names = FALSE, no.. = TRUE)
stop_if(setequal(actual_files, c("INPUT_MANIFEST.csv", manifest$File)),
  "Expanded publication input contains an unmanifested or missing file.")
for (index in seq_len(nrow(manifest))) {
  item_path <- file.path(publication_dir, manifest$File[[index]])
  stop_if(file.exists(item_path), paste("Manifested expanded input is missing:", manifest$File[[index]]))
  stop_if(file.info(item_path)$size == manifest$Bytes[[index]], paste("Expanded input byte count changed:", manifest$File[[index]]))
  stop_if(digest::digest(file = item_path, algo = "md5", serialize = FALSE) == manifest$MD5[[index]],
    paste("Expanded input checksum changed:", manifest$File[[index]]))
}

formal_stems <- c(
  "Figure_1_DDR_fraction", "Figure_2a_whole_proteome_DDR_UpSet", "Figure_2b_Kla_DDR_UpSet",
  "Figure_2c_DDR_pathway_matrices_tissues", "Figure_2d_DDR_pathway_summary_tumor_tissue",
  "Figure_2e_DDR_pathway_summary_non_tumor_tissue", "Figure_3a_reference_regulator_percentiles",
  "Figure_3b_Kla_regulator_percentiles", "Supplementary_Figure_S1a_whole_proteome_UpSet",
  "Supplementary_Figure_S1b_Kla_proteome_UpSet", "Supplementary_Figure_S2a_DDR_pathway_matrices_cell_lines",
  "Supplementary_Figure_S2b_DDR_pathway_summary_cancer_cell_lines",
  "Supplementary_Figure_S2c_DDR_pathway_summary_normal_cell_lines"
)
candidate_stems <- file.path(legacy_output_dir, "Figure_1_DDR_fraction_candidate_category_boxplot_refined")
formal_files <- file.path(formal_output_dir, paste0(rep(formal_stems, each = 2L), rep(c(".png", ".pdf"), length(formal_stems))))
stop_if(all(file.exists(formal_files)),
  "One or more expanded formal figure files are missing.")
candidate_files <- c(paste0(candidate_stems, ".png"), paste0(candidate_stems, ".pdf"))
stop_if(all(file.exists(candidate_files)),
  "One or more expanded candidate figure files are missing.")
stop_if(file.exists(file.path(legacy_output_dir, "figure1_category_one_way_anova.csv")) &&
          file.exists(file.path(legacy_output_dir, "figure1_category_boxplot_mean_median.csv")),
  "Expanded Figure 1 statistics sidecars are missing.")
pathway_manifest <- fread(file.path(pathway_output_dir, "pathway_summary_by_pathway_manifest.csv"))
pathway_anova <- fread(file.path(pathway_output_dir, "pathway_summary_two_way_anova.csv"))
stop_if(nrow(pathway_manifest) == 7L &&
          all(pathway_manifest$Dataset == "Lactylome (Kla)") &&
          uniqueN(pathway_manifest$Pathway) == 7L &&
          all(pathway_manifest$CategoryPanels == 4L) &&
          all(pathway_manifest$BoxesPerFigure == 8L),
  "Expanded pathway-summary manifest does not contain seven four-category Kla plots.")
stop_if(all(file.exists(file.path(pathway_output_dir, pathway_manifest$PNG))) &&
          all(file.exists(file.path(pathway_output_dir, pathway_manifest$PDF))),
  "One or more expanded pathway-summary boxplot files are missing.")
stop_if(nrow(pathway_anova) == 21L &&
          setequal(unique(pathway_anova$Term), c("CategoryFactor", "DirectionFactor", "CategoryFactor:DirectionFactor")) &&
          all(pathway_anova$Test == "two-way ANOVA"),
  "Expanded pathway-summary two-way ANOVA sidecar is incomplete.")

message("PASS: PXD064038-only ESCC inclusion scope, expanded S4, provenance, manifest, and figures are consistent.")
