#!/usr/bin/env Rscript
# Audit-only observation contract builder (no model training).
# Reads only local frozen inputs + repair release; never touches raw proteomics stores.
# Output: outputs/20260920_lactylation_ml/observation_contract.csv + audit hash log
suppressPackageStartupMessages({library(data.table); library(digest)})
set.seed(25)
root <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE))[1L]), "..", ".."))
setwd(root)
out_dir <- "outputs/20260920_lactylation_ml"
dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)

gs <- fread("data/publication_input/group_summary_31.csv")
reg <- fread("outputs/20260919_repair_release/core/regulator_rna_qsmooth_28materials.csv")
map28 <- unique(reg[, .(GroupID, PXD, SampleGroup, ReferenceKey)])
# expected N from contract
contract <- fread("config/rnaseq_expression_extraction_contract_20260919.csv")
expN <- contract[, .(GroupIDs, ReferenceKey, ExpectedN)]
# shared-reference flag
map28[, SharedReference := ReferenceKey %in% c("HCT116_control", "HK2_untreated")]
# linkage class (from sample_design_31 where available; else external label)
sdes <- fread("data/candidate/sample_design_31.csv")

obs <- merge(gs, map28, by=c("PXD", "SampleGroup"), all.x=TRUE)
# KLA31_15/16 (shared HCT116) and KLA31_28 (shared HK2) have no regulator GroupID row (28-material table);
# fill them explicitly as shared references with no independent RNA.
fill <- data.table(
  GroupID=c("KLA31_15", "KLA31_16", "KLA31_28"),
  PXD=c("PXD053474", "PXD066351", "PXD078736"),
  ReferenceKey=c("HCT116_control", "HCT116_control", "HK2_untreated"),
  SharedReference=c(TRUE, TRUE, TRUE))
for (i in seq_len(nrow(fill))) {
  idx <- which(obs$PXD==fill$PXD[i] & is.na(obs$GroupID))
  if (length(idx)>=1L) {
    # assign in RowOrder order
    idx <- idx[1L]
    obs$GroupID[idx] <- fill$GroupID[i]
    obs$ReferenceKey[idx] <- fill$ReferenceKey[i]
    obs$SharedReference[idx] <- TRUE
  }
}
obs[, ObservationID := sprintf("OBS_%02d_%s", seq_len(.N), gsub("[^A-Za-z0-9]+", "_", paste0(PXD, "_", SampleGroup)))]
obs[, RNASource := ReferenceKey]
obs[, ProteomeSource := KlaEvidenceFile]
obs[, StudyID := PXD]
obs[, MaterialType := Category]
obs[, LabelDefinition := "Kla detected protein count (deduplicated BaseAccession; detection feature only; NOT concentration/occupancy)"]
obs[, LabelUnit := "count"]
obs[, LabelLevel := "material_aggregate"]
obs[, Linkage := fcase(
  SharedReference==TRUE, "external_reference_shared",
  PXD %in% c("PXD046800","PXD050470","PXD066054"), "same_study_external_RNA",
  default="external_reference_match")]
obs[, SharedSource := fcase(
  ReferenceKey=="HCT116_control", "HCT116_control shared across KLA31_14/15/16",
  ReferenceKey=="HK2_untreated", "HK2_untreated shared across KLA31_27/28",
  default="none")]
obs[, BatchNote := paste0("Proteome study=", PXD, "; RNA source=", ReferenceKey, "; cross-MS intensities uncalibrated")]
obs[, ComparabilityLimit := "Cross-study Kla counts NOT comparable as unified continuous label; no total-protein denominator at same-sample level for most groups; external RNA is NOT paired measurement"]
setcolorder(obs, c("ObservationID","GroupID","ReferenceKey","RNASource","ProteomeSource","StudyID","MaterialType",
  "LabelDefinition","LabelUnit","LabelLevel","Linkage","SharedSource","BatchNote","ComparabilityLimit",
  "PXD","SampleGroup","Category","KlaProteinCount","ReferenceProteinCount","KlaDdrProteinCount","SharedReference"))
setorder(obs, GroupID)
fwrite(obs, file.path(out_dir, "observation_contract.csv"))
# connected-component grouping over shared RNA references (leakage graph)
leak <- obs[, .(GroupID, ReferenceKey, StudyID)]
leak[, ConnGroup := .GRP, by=ReferenceKey]
fwrite(leak[order(ConnGroup, GroupID)], file.path(out_dir, "connectivity_groups.csv"))
# hashes
files <- c("outputs/20260920_lactylation_ml/observation_contract.csv",
  "outputs/20260920_lactylation_ml/connectivity_groups.csv",
  "config/lactylation_ml/task_feasibility.csv",
  "config/lactylation_ml/feature_dictionary.csv",
  "config/lactylation_ml/frozen_params.csv")
hh <- data.table(File=files, SHA256=vapply(files, function(f) digest(file=f, algo="sha256"), character(1)))
fwrite(hh, "audit/20260920_lactylation_ml/input_hashes.csv")
cat(sprintf("OBS_CONTRACT_OK: n=%d shared_refs=%d indep_RNA_refs=%d\n", nrow(obs),
  sum(obs$SharedReference, na.rm=TRUE), uniqueN(obs$ReferenceKey, na.rm=TRUE)))
