#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) normalizePath(args[[1]]) else normalizePath(".")
table_dir <- file.path(project_root, "reanalysis", "results", "tables")
figure_dir <- file.path(project_root, "reanalysis", "results", "figures")

assert <- function(condition, message) {
  if (!isTRUE(condition)) {
    stop(message, call. = FALSE)
  }
}

long <- read.csv(
  file.path(table_dir, "kla_regulator_cell_type_long.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
summary <- read.csv(
  file.path(table_dir, "kla_regulator_detection_summary.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

assert(nrow(long) == 49 * 40, "Expected 49 role entries across 40 cell/tissue groups")
assert(
  length(unique(long$SampleGroupID)) == 40,
  "The expanded plot must retain all 40 pair-ready lactylome groups"
)
assert(
  sum(!duplicated(long$SampleGroupID) & long$GeneLevelAuditStatus == "逐蛋白可审计") >= 30,
  "At least 30 groups should have auditable protein-identity evidence"
)
assert(
  sum(!duplicated(long$SampleGroupID) & long$GeneLevelAuditStatus == "逐蛋白可审计") == 37,
  "Exactly 37 groups should have auditable protein-identity evidence"
)
assert(
  all(
    unique(long$PXD[long$GeneLevelAuditStatus != "逐蛋白可审计"]) ==
      "PXD037371"
  ),
  "Only the three PXD037371 groups should remain unavailable"
)
assert(length(unique(long$GeneSymbol)) == 48, "Expected 48 unique regulator genes")
assert(
  any(
    long$GeneSymbol == "HDAC8" &
      long$Role == "Eraser"
  ) &&
    any(
      long$GeneSymbol == "HDAC8" &
        long$Role == "Writer-Eraser"
    ),
  "A regulator assigned to multiple mechanisms must retain every distinct role"
)
assert(
  all(is.na(long$ExactSiteCount) | long$ExactSiteCount >= 0),
  "Exact Kla site counts must be non-negative when available"
)
assert(
  all(
    !long$Detected |
      !is.na(long$ExactSiteCount) |
      grepl(
        "site|protein|modified_peptide|modified_precursor",
        long$EvidenceResolution
      )
  ),
  "Every detected regulator must have site-, peptide-, or protein-level evidence"
)
assert(
  any(
    long$GeneSymbol == "HAT1" &
      long$PXD == "PXD054919" &
      long$SampleGroup == "A549" &
      long$ExactSiteCount == 1
  ),
  "PXD054919 must contribute the HAT1 K15 site in A549"
)
assert(
  sum(summary$DetectedGroupCount > 0) >= 19,
  "At least 19 unique regulator genes should be detected"
)
assert(
  all(grepl(
    "Lacty_PeptideGroups",
    long$SourceFile[long$PXD == "PXD046800" & long$Detected],
    fixed = TRUE
  )),
  "PXD046800 detections must come from explicit lactylated peptide groups"
)
assert(
  all(grepl(
    "DPLa-MSstats_Input",
    long$SourceFile[long$PXD == "PXD066351" & long$Detected],
    fixed = TRUE
  )),
  "PXD066351 detections must come from explicit Lac(K) PTM rows"
)

for (name in c(
  "kla_regulator_cell_type_detection_landscape.png",
  "kla_regulator_cell_type_detection_landscape.pdf"
)) {
  path <- file.path(figure_dir, name)
  assert(file.exists(path), paste("Missing figure:", name))
  assert(file.info(path)$size > 10000, paste("Figure is unexpectedly small:", name))
}

message("All Kla regulator landscape tests passed.")
