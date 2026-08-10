#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args)) normalizePath(args[[1]]) else normalizePath(".")
table_dir <- file.path(project_root, "reanalysis", "results", "tables")

audit <- read.csv(
  file.path(table_dir, "strict_reference_material_identity_audit.csv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

stopifnot(nrow(audit) == 33)
stopifnot(!anyDuplicated(audit[c("KlaPXD", "SampleGroup")]))
stopifnot(sum(audit$MaterialIdentityMatch) == 33)
stopifnot(sum(audit$AnalysisDecision == "excluded_no_exact_material_reference") == 0)

hippocampus <- audit[
  audit$KlaPXD == "PXD050470" &
    audit$SampleGroup == "human hippocampus",
]
stopifnot(nrow(hippocampus) == 1)
stopifnot(hippocampus$ReferencePXD == "PXD050470")
stopifnot(
  hippocampus$ReferenceFile ==
    "data/PXD050470/supplementary/prca2331-sup-0006-tables4.xlsx"
)
stopifnot(hippocampus$ReferenceSubset == "H072;H081;H0187")
stopifnot(hippocampus$ExperimentalStateMatch == "exact_same_samples")

hcc <- audit[audit$KlaPXD == "PXD075377", ]
stopifnot(nrow(hcc) == 2)
stopifnot(
  hcc$DistinctSubsetVerified[hcc$SampleGroup == "HCC"] == "CISs_sheet_only"
)
stopifnot(
  hcc$DistinctSubsetVerified[hcc$SampleGroup == "adjacent liver"] ==
    "ANTs_sheet_only"
)

scar <- audit[audit$KlaPXD == "PXD046800", ]
stopifnot(
  scar$DistinctSubsetVerified[scar$SampleGroup == "hypertrophic scar"] ==
    "HSP1_HSP4_only"
)
stopifnot(
  scar$DistinctSubsetVerified[scar$SampleGroup == "adjacent skin"] ==
    "NSP1_NSP4_only"
)

hk2 <- audit[audit$ReferencePXD == "PXD072220", ]
stopifnot(nrow(hk2) == 2)
stopifnot(all(
  hk2$DistinctSubsetVerified == "same_HK2_baseline_intentionally_reused"
))
stopifnot(all(hk2$ExperimentalStateMatch == "baseline_reference_only"))

message("Strict reference material-identity tests passed.")
