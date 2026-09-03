#!/usr/bin/env Rscript

# Validate the expanded Supplementary Tables S1-S6 for the candidate 31-group scope
# with teacher-requested renumbering:
# - Table S1: Kla Data
# - Table S2: Reference Data
# - Table S3: Human DDR GO Annotations
# - Table S4: Venn Membership (formerly Table S6)
# - Table S5: Pathway Protein Ranking (formerly Table S4)
# - Table S6: Lactylation Regulators (formerly Table S5)

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
})

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args)) normalizePath(args[[1L]], mustWork = TRUE) else normalizePath(".", mustWork = TRUE)
stop_if <- function(condition, message) if (!isTRUE(condition)) stop(message, call. = FALSE)

scope_tag <- Sys.getenv("KLA_ESCC_SCOPE_TAG", unset = "escc_inclusion_20260903_pxd065830_tumor_reference")
candidate_pub_dir <- file.path(project_root, "data", "candidate", scope_tag, "publication_input")
supp_dir <- file.path(project_root, "results", scope_tag, "supplementary")

stop_if(dir.exists(candidate_pub_dir), paste("Candidate publication directory missing:", candidate_pub_dir))
stop_if(dir.exists(supp_dir), paste("Candidate supplementary directory missing:", supp_dir))

# Table S1: Kla Data
s1_path <- file.path(supp_dir, "Supplementary_Table_S1_Kla_Data.xlsx")
stop_if(file.exists(s1_path), "Missing Supplementary_Table_S1_Kla_Data.xlsx")
s1_groups <- as.data.table(read_excel(s1_path, sheet = "Group_Summary"))
s1_kla <- as.data.table(read_excel(s1_path, sheet = "Kla_Protein_Membership"))
s1_ddr <- as.data.table(read_excel(s1_path, sheet = "Kla_DDR_Membership"))

stop_if(nrow(s1_groups) == 31L, "S1 Group_Summary must contain exactly 31 groups.")
escc_row <- s1_groups[PXD == "PXD064038" & SampleGroup == "MEC and NEC ESCC groups"]
stop_if(nrow(escc_row) == 1L, "PXD064038 group missing in S1.")
stop_if(escc_row$KlaProteinCount == 1239L && escc_row$KlaDdrProteinCount == 92L, "PXD064038 Kla counts in S1 are incorrect.")
input_kla <- fread(file.path(candidate_pub_dir, "kla_protein_membership_30.csv"))
stop_if(nrow(s1_kla) == nrow(input_kla), "S1 Kla membership row count does not match input.")
stop_if(nrow(s1_ddr) == sum(input_kla$IsDdr == TRUE), "S1 Kla-DDR membership row count does not match input.")

# Table S2: Reference Data
s2_path <- file.path(supp_dir, "Supplementary_Table_S2_Reference_Data.xlsx")
stop_if(file.exists(s2_path), "Missing Supplementary_Table_S2_Reference_Data.xlsx")
s2_groups <- as.data.table(read_excel(s2_path, sheet = "Reference_Group_Summary"))
s2_ref <- as.data.table(read_excel(s2_path, sheet = "Reference_Protein_Membership"))
s2_ddr <- as.data.table(read_excel(s2_path, sheet = "Reference_DDR_Membership"))

stop_if(nrow(s2_groups) == 31L, "S2 Reference_Group_Summary must contain exactly 31 groups.")
escc_ref_row <- s2_groups[PXD == "PXD064038" & SampleGroup == "MEC and NEC ESCC groups"]
stop_if(escc_ref_row$ReferencePXD == "PXD065830", "PXD064038 must reference PXD065830 in S2.")
stop_if(escc_ref_row$ReferenceProteinCount == 8083L && escc_ref_row$ReferenceDdrProteinCount == 420L,
  "PXD064038 reference counts in S2 are incorrect.")
input_ref <- fread(file.path(candidate_pub_dir, "reference_protein_membership_30.csv"))
stop_if(nrow(s2_ref) == nrow(input_ref), "S2 reference membership row count does not match input.")
stop_if(nrow(s2_ddr) == sum(input_ref$IsDdr == TRUE), "S2 reference DDR membership row count does not match input.")

# Table S3: Human DDR GO Annotations
s3_path <- file.path(supp_dir, "Supplementary_Table_S3_Human_DDR_GO_Annotations.xlsx")
stop_if(file.exists(s3_path), "Missing Supplementary_Table_S3_Human_DDR_GO_Annotations.xlsx")
s3_data <- as.data.table(read_excel(s3_path, sheet = "Human_DDR_GO_Annotations", col_types = "text"))
input_go <- fread(file.path(candidate_pub_dir, "human_ddr_go_annotations.tsv"), sep = "\t", quote = "")
stop_if(nrow(s3_data) == nrow(input_go), "S3 GO annotations row count does not match input.")

# Table S4: Venn Membership (Renamed from former Table S6)
s4_path <- file.path(supp_dir, "Supplementary_Table_S4_Venn_Membership.xlsx")
stop_if(file.exists(s4_path), "Missing Supplementary_Table_S4_Venn_Membership.xlsx (formerly S6)")
s4_reg <- as.data.table(read_excel(s4_path, sheet = "Region_Counts"))
stop_if(nrow(s4_reg) == 60L, "S4 Region_Counts must contain 60 rows (4 analyses x 15 regions).")
venn_totals <- s4_reg[, .(Total = sum(ProteinCount)), by = Analysis]
stop_if(venn_totals[Analysis == "AllKla", Total] == 5814L, "S4 AllKla total mismatch.")
stop_if(venn_totals[Analysis == "KlaDDR", Total] == 401L, "S4 KlaDDR total mismatch.")
stop_if(venn_totals[Analysis == "Reference", Total] == 24397L, "S4 Reference total mismatch.")
stop_if(venn_totals[Analysis == "ReferenceDDR", Total] == 836L, "S4 ReferenceDDR total mismatch.")

for (spec in list(
  c("AllKla_Members", "venn_all_kla.csv"),
  c("KlaDDR_Members", "venn_kla_ddr.csv"),
  c("Reference_Members", "venn_reference.csv"),
  c("ReferenceDDR_Members", "venn_reference_ddr.csv")
)) {
  m_sheet <- as.data.table(read_excel(s4_path, sheet = spec[1]))
  m_file <- fread(file.path(candidate_pub_dir, spec[2]))
  stop_if(nrow(m_sheet) == nrow(m_file), paste("S4 membership mismatch for", spec[1]))
}

# Table S5: Pathway Protein Ranking (Renamed from former Table S4)
s5_path <- file.path(supp_dir, "Supplementary_Table_S5_Pathway_Protein_Ranking.xlsx")
stop_if(file.exists(s5_path), "Missing Supplementary_Table_S5_Pathway_Protein_Ranking.xlsx (formerly S4)")
s5_expected <- c(NonTumorTissues = 183L, TumorTissues = 192L, CancerCellLines = 381L, NormalCellLines = 292L)
for (sheet_name in names(s5_expected)) {
  panel <- as.data.table(read_excel(s5_path, sheet = sheet_name))
  panel <- panel[!is.na(BaseAccession) & nzchar(trimws(BaseAccession))]
  stop_if(nrow(panel) == s5_expected[[sheet_name]], paste("S5 row count mismatch for", sheet_name))
}

# Table S6: Lactylation Regulators (Renamed from former Table S5)
s6_path <- file.path(supp_dir, "Supplementary_Table_S6_Lactylation_Regulators.xlsx")
stop_if(file.exists(s6_path), "Missing Supplementary_Table_S6_Lactylation_Regulators.xlsx (formerly S5)")
s6_ann <- as.data.table(read_excel(s6_path, sheet = "Regulator_Annotations"))
stop_if(nrow(unique(s6_ann[!is.na(BaseAccession) & nzchar(trimws(BaseAccession)), .(Role, BaseAccession)])) == 49L,
  "S6 regulator count must be 49.")

# Verify obsolete mappings do not exist to ensure no filename ambiguity
stop_if(!file.exists(file.path(supp_dir, "Supplementary_Table_S4_Pathway_Protein_Ranking.xlsx")),
  "Obsolete file Supplementary_Table_S4_Pathway_Protein_Ranking.xlsx should not exist in 31-group supplementary dir.")
stop_if(!file.exists(file.path(supp_dir, "Supplementary_Table_S6_Venn_Membership.xlsx")),
  "Obsolete file Supplementary_Table_S6_Venn_Membership.xlsx should not exist in 31-group supplementary dir.")

message("PASS: Renamed Supplementary Tables S1-S6 (S4=Venn, S5=Pathways, S6=Regulators) verified successfully for 31-group scope.")
