#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(dplyr))

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) normalizePath(args[[1]]) else normalizePath(".")
table_dir <- file.path(project_root, "reanalysis", "results", "tables")
figure_dir <- file.path(project_root, "reanalysis", "results", "figures")

required <- c(
  file.path(table_dir, "kla_regulator_intensity_availability_audit.csv"),
  file.path(table_dir, "kla_regulator_intensity_sample_level_long.csv"),
  file.path(table_dir, "kla_regulator_normalized_intensity_long.csv"),
  file.path(table_dir, "kla_regulator_within_pxd_zscore_long.csv"),
  file.path(table_dir, "kla_regulator_intensity_plot_exclusions.csv"),
  file.path(table_dir, "kla_regulator_intensity_id_mapping_audit.csv"),
  file.path(table_dir, "kla_regulator_intensity_pure_white_audit.csv"),
  file.path(table_dir, "kla_regulator_heatmap_axis_order.csv"),
  file.path(
    table_dir,
    "kla_regulator_whole_proteome_heatmap_display_long.csv"
  ),
  file.path(figure_dir, "kla_regulator_cross_study_relative_intensity_heatmap.png"),
  file.path(figure_dir, "kla_regulator_cross_study_relative_intensity_heatmap.pdf"),
  file.path(figure_dir, "kla_regulator_within_pxd_zscore_heatmaps.pdf")
)
stopifnot(all(file.exists(required)))

audit <- read.csv(required[[1]], check.names = FALSE, stringsAsFactors = FALSE)
normalized <- read.csv(required[[3]], check.names = FALSE, stringsAsFactors = FALSE)
sample_level <- read.csv(required[[2]], check.names = FALSE, stringsAsFactors = FALSE)
zscore <- read.csv(required[[4]], check.names = FALSE, stringsAsFactors = FALSE)
id_audit <- read.csv(required[[6]], check.names = FALSE, stringsAsFactors = FALSE)
pure_white_audit <- read.csv(required[[7]], check.names = FALSE, stringsAsFactors = FALSE)
axis_audit <- read.csv(required[[8]], check.names = FALSE, stringsAsFactors = FALSE)
whole_display <- read.csv(
  required[[9]],
  check.names = FALSE,
  stringsAsFactors = FALSE
)

stopifnot(nrow(audit) == 40)
stopifnot(sum(audit$定量可用) == 37)
stopifnot("严格配对分析纳入" %in% names(audit))
stopifnot("严格配对排除原因" %in% names(audit))
stopifnot(sum(audit$严格配对分析纳入) == 33)
stopifnot(sum(audit$定量可用 & !audit$严格配对分析纳入) == 4)
stopifnot(nrow(axis_audit) == 33)
stopifnot(all(axis_audit$QuantificationAvailable))
stopifnot(all(!is.na(axis_audit$WholeProteomeDisplayRowOrder)))
stopifnot(audit$定量可用[audit$PXD == "PXD050470"])
stopifnot(audit$定量可用[audit$PXD == "PXD028737"])
stopifnot(audit$定量可用[audit$PXD == "PXD073311"])
stopifnot(audit$定量可用[audit$PXD == "PXD075014"])
stopifnot(all(audit$PXD[!audit$定量可用] == "PXD037371"))
stopifnot(nrow(normalized) == 33 * nrow(unique(normalized[c("Role", "GeneSymbol")])))
sample_keys <- unique(sample_level[c("PXD", "SampleGroup", "QuantSample")])
regulator_keys <- unique(
  sample_level[c("GeneSymbol", "RegulatorBaseAccession")]
)
stopifnot(nrow(sample_level) == nrow(sample_keys) * nrow(regulator_keys))
excluded_keys <- c(
  "PXD062720::bladder cancer cells treated with EPI",
  "PXD063047::severe preeclampsia placenta",
  "PXD064038::MEC and NEC ESCC groups",
  "PXD075014::AC16 control and hypoxia"
)
normalized_keys <- paste(normalized$PXD, normalized$SampleGroup, sep = "::")
axis_keys <- paste(axis_audit$PXD, axis_audit$SampleGroup, sep = "::")
stopifnot(length(intersect(normalized_keys, excluded_keys)) == 0)
stopifnot(length(intersect(axis_keys, excluded_keys)) == 0)
stopifnot("PXD063047::normal pregnancy placenta" %in% axis_keys)
stopifnot(all(
  normalized$RelativeKlaPercentile[
    !is.na(normalized$RelativeKlaPercentile)
  ] >= 0 &
    normalized$RelativeKlaPercentile[
      !is.na(normalized$RelativeKlaPercentile)
    ] <= 100
))
stopifnot(all(sample_level$Signal >= 0))
stopifnot(all(sample_level$WithinSamplePercentile >= 0))
stopifnot(all(sample_level$WithinSamplePercentile <= 100))
stopifnot("RegulatorBaseAccession" %in% names(sample_level))
stopifnot("IdentityMatchMode" %in% names(sample_level))
row_order_check <- normalized |>
  group_by(PXD, SampleGroup) |>
  summarise(
    ComparisonRowOrderCount = n_distinct(ComparisonRowOrder),
    .groups = "drop"
  )
stopifnot(all(row_order_check$ComparisonRowOrderCount == 1L))
stopifnot("RegulatorDisplayName" %in% names(normalized))
stopifnot(all(
  normalized$RegulatorDisplayName[normalized$RegulatorBaseAccession == "Q92830"] ==
    "GCN5 (KAT2A)"
))
role_levels <- c("Writer", "Eraser", "Writer-Eraser", "Reader")
kla_columns <- normalized |>
  arrange(match(Role, role_levels), RoleEntryOrder) |>
  distinct(Role, RegulatorDisplayName, RoleEntryOrder)
whole_columns <- whole_display |>
  arrange(match(Role, role_levels), RoleEntryOrder) |>
  distinct(Role, RegulatorDisplayName, RoleEntryOrder)
stopifnot(identical(kla_columns, whole_columns))
stopifnot(all(
  sample_level$IdentityMatchMode[sample_level$Detected] == "BaseAccession"
))
stopifnot(all(
  !is.na(sample_level$RegulatorBaseAccession[sample_level$Detected]) &
    nzchar(sample_level$RegulatorBaseAccession[sample_level$Detected])
))
stopifnot(all(id_audit$IdentityMatchMode == "UniProt_BaseAccession_only"))
stopifnot(sum(id_audit$GeneSymbolFallbackCount) == 0)
stopifnot(all(id_audit$MappingSource == "UniProtKB reviewed Homo sapiens"))

expected_pure_white <- data.frame(
  PXD = c(
    "PXD014870",
    "PXD028488", "PXD028488", "PXD028488",
    "PXD064912"
  ),
  SampleGroup = c(
    "MCF7",
    "HEK293T", "HCT116", "TALL-104",
    "human sperm"
  ),
  stringsAsFactors = FALSE
)
observed_pure_white <- pure_white_audit[
  pure_white_audit$QuantificationAvailable & pure_white_audit$PureWhite,
  c("PXD", "SampleGroup")
]
expected_keys <- sort(paste(expected_pure_white$PXD, expected_pure_white$SampleGroup, sep = "::"))
observed_keys <- sort(paste(observed_pure_white$PXD, observed_pure_white$SampleGroup, sep = "::"))
stopifnot(identical(observed_keys, expected_keys))
stopifnot(
  pure_white_audit$Explanation[
    pure_white_audit$PXD == "PXD014870" &
      pure_white_audit$SampleGroup == "MCF7"
  ] == "sparse_detection_and_group_median_zero"
)
stopifnot(all(
  pure_white_audit$Explanation[
    pure_white_audit$PXD %in% c("PXD028488", "PXD064912")
  ] == "no_target_regulator_accession_hit_in_Kla_source"
))
stopifnot(all(
  pure_white_audit$DetectedSampleRows[
    pure_white_audit$PXD == "PXD014870" &
      pure_white_audit$SampleGroup == "MCF7"
  ] > 0
))
stopifnot(all(
  grepl(
    "Lacty_PeptideGroups",
    audit$使用的定量字段[audit$PXD == "PXD046800"],
    fixed = TRUE
  ) |
    grepl(
      "lactylated peptide",
      audit$使用的定量字段[audit$PXD == "PXD046800"],
      fixed = TRUE
    )
))
stopifnot(all(
  grepl(
    "MSstats Lac(K) PTM.Quantity",
    audit$使用的定量字段[audit$PXD == "PXD066351"],
    fixed = TRUE
  )
))
stopifnot(all(
  grepl(
    "Lacty_PeptideGroups",
    audit$SourceFile[audit$PXD == "PXD046800"],
    fixed = TRUE
  )
))
stopifnot(all(
  grepl(
    "DPLa-MSstats_Input",
    audit$SourceFile[audit$PXD == "PXD066351"],
    fixed = TRUE
  )
))
percentile_variation <- aggregate(
  WithinSamplePercentile ~ PXD + SampleGroup + QuantSample,
  sample_level[sample_level$Detected, ],
  function(values) length(unique(round(values, 8)))
)
stopifnot(any(percentile_variation$WithinSamplePercentile > 1))
z_variation <- aggregate(
  WithinPXDZ ~ PXD + GeneSymbol,
  zscore,
  function(values) length(unique(round(values, 8)))
)
stopifnot(any(z_variation$WithinPXDZ > 1))

cat("Kla regulator intensity tests passed.\n")
