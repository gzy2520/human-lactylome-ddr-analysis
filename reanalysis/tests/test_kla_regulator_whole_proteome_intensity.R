#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) normalizePath(args[[1]]) else normalizePath(".")
table_dir <- file.path(project_root, "reanalysis", "results", "tables")
figure_dir <- file.path(project_root, "reanalysis", "results", "figures")
script_path <- file.path(
  project_root,
  "reanalysis",
  "scripts",
  "analyze_kla_regulator_whole_proteome_intensity.R"
)

assert <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

audit <- read.csv(
  file.path(
    table_dir,
    "kla_regulator_whole_proteome_intensity_availability_audit.csv"
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
kla_scope <- read.csv(
  file.path(
    table_dir,
    "kla_regulator_intensity_availability_audit.csv"
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
normalized <- read.csv(
  file.path(
    table_dir,
    "kla_regulator_whole_proteome_normalized_long.csv"
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
algorithm_audit <- read.csv(
  file.path(
    table_dir,
    "kla_regulator_whole_proteome_algorithm_audit.csv"
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
ensembl_audit <- read.csv(
  file.path(
    table_dir,
    "kla_regulator_whole_proteome_ensembl_mapping_audit.csv"
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
regulator_sample <- read.csv(
  file.path(
    table_dir,
    "kla_regulator_whole_proteome_regulator_sample_level_long.csv"
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
scope_exclusions <- read.csv(
  file.path(
    table_dir,
    "kla_regulator_whole_proteome_strict_scope_exclusions.csv"
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
heatmap_rows <- read.csv(
  file.path(
    table_dir,
    "kla_regulator_whole_proteome_heatmap_rows.csv"
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
heatmap_long <- read.csv(
  file.path(
    table_dir,
    "kla_regulator_whole_proteome_heatmap_display_long.csv"
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_figures <- c(
  file.path(
    figure_dir,
    "kla_regulator_whole_proteome_relative_intensity_heatmap.png"
  ),
  file.path(
    figure_dir,
    "kla_regulator_whole_proteome_relative_intensity_heatmap.pdf"
  ),
  file.path(
    figure_dir,
    "kla_regulator_whole_proteome_relative_intensity_heatmap_zh.png"
  ),
  file.path(
    figure_dir,
    "kla_regulator_whole_proteome_relative_intensity_heatmap_zh.pdf"
  ),
  file.path(
    figure_dir,
    "kla_regulator_whole_proteome_relative_intensity_heatmap_en.png"
  ),
  file.path(
    figure_dir,
    "kla_regulator_whole_proteome_relative_intensity_heatmap_en.pdf"
  )
)
assert(all(file.exists(required_figures)), "Missing bilingual whole-proteome heatmap files")
assert(all(file.info(required_figures)$size > 10000), "Whole-proteome heatmap file is unexpectedly small")
assert(nrow(heatmap_rows) == 30, "Whole-proteome heatmap must contain 30 unique reference rows")
assert(
  identical(heatmap_rows$HeatmapDisplayRowOrder, seq_len(30)),
  "Whole-proteome row order is invalid"
)
assert(all(c("CategoryEn", "RowLabelZh", "RowLabelEn") %in% names(heatmap_rows)))
assert(all(file.exists(required_figures)), "Whole-proteome heatmap files are missing")
assert(all(file.info(required_figures)$size > 10000), "Whole-proteome figures are unexpectedly small")
assert(file.exists(script_path), "Whole-proteome analysis script is missing")
assert(
  nrow(algorithm_audit) == 15 &&
    algorithm_audit$Value[algorithm_audit$Field == "AlgorithmVersion"] ==
      "whole_proteome_regulator_rank_v6_unique_reference_display",
  "Whole-proteome algorithm audit is missing or outdated"
)
assert(
  algorithm_audit$Value[
    algorithm_audit$Field == "KlaSignalSubstitution"
  ] == "never; Kla-enriched intensity is not used as a fallback",
  "Whole-proteome heatmap must not fall back to Kla signal"
)

assert(nrow(audit) == 33, "Expected exactly 33 strict-reference sample groups in the audit")
assert(
  nrow(heatmap_rows) == 30 &&
    !anyDuplicated(heatmap_rows$ReferenceDisplayKey),
  "Whole-proteome heatmap must contain exactly 30 unique reference rows"
)
category_order <- c(
  "normal_tissue",
  "cancer_tissue",
  "normal_cells",
  "cancer_cells"
)
assert(
  all(c("Category", "CategoryZh", "HeatmapRowOrder") %in% names(audit)),
  "Whole-proteome heatmap audit must record four-class row grouping"
)
ordered_audit <- audit[order(audit$HeatmapRowOrder), ]
assert(
  identical(unique(ordered_audit$Category), category_order),
  "Heatmap rows must be ordered as normal tissue, cancer tissue, normal cells, cancer cells"
)
assert(
  identical(
    as.integer(table(factor(audit$Category, levels = category_order))),
    c(9L, 2L, 9L, 13L)
  ),
  "Strict-reference heatmap category counts must be 9, 2, 9, and 13"
)
assert(
  identical(
    as.integer(table(factor(heatmap_rows$Category, levels = category_order))),
    c(9L, 2L, 8L, 11L)
  ),
  "Unique-reference heatmap category counts must be 9, 2, 8, and 11"
)
shared_reference_expectations <- list(
  "PXD058534;PXD078736",
  "PXD014870;PXD060185",
  "PXD028488;PXD053474"
)
for (linked_pxd in shared_reference_expectations) {
  shared_row <- heatmap_rows[heatmap_rows$LinkedKlaPXD == linked_pxd, ]
  assert(
    nrow(shared_row) == 1 && shared_row$LinkedKlaStudyCount == 2,
    paste("Missing unique shared-reference heatmap row for", linked_pxd)
  )
}
assert(
  sum(heatmap_rows$LinkedKlaStudyCount == 2) == 3 &&
    sum(heatmap_rows$LinkedKlaStudyCount == 1) == 27,
  "Only HK-2, MCF7, and HCT116 may share a displayed reference row"
)
fibroblast_rows <- heatmap_rows[
  heatmap_rows$MaterialCluster == "human_fibroblasts",
]
assert(
  nrow(fibroblast_rows) == 2 &&
    setequal(
      fibroblast_rows$MaterialDisplayName,
      c(
        "human fibroblasts mock and HCMV or HSV-1",
        "human fibroblasts mock and HCMV"
      )
    ) &&
    diff(sort(fibroblast_rows$HeatmapDisplayRowOrder)) == 1,
  paste(
    "The two fibroblast studies must remain adjacent while retaining",
    "their distinct infection labels"
  )
)
display_labels_preserve_sample_groups <- mapply(
  function(linked_groups, row_label) {
    sample_groups <- unique(strsplit(
      linked_groups,
      ";",
      fixed = TRUE
    )[[1]])
    all(vapply(
      sample_groups,
      function(sample_group) grepl(sample_group, row_label, fixed = TRUE),
      logical(1)
    ))
  },
  heatmap_rows$LinkedKlaSampleGroup,
  heatmap_rows$RowLabel
)
assert(
  all(display_labels_preserve_sample_groups),
  paste(
    "Every original sample/treatment name must remain visible in its",
    "deduplicated heatmap row label"
  )
)
assert(
  all(audit$WholeProteomeQuantAvailable),
  "Every strict-reference sample group must have whole-proteome quantification"
)
kla_scope_keys <- with(
  kla_scope[kla_scope$定量可用 %in% c(TRUE, "TRUE", "True", 1, "1"), ],
  paste(PXD, 样本组, sep = "__")
)
whole_proteome_keys <- with(audit, paste(PXD, SampleGroup, sep = "__"))
assert(
  !length(setdiff(whole_proteome_keys, kla_scope_keys)) &&
    setequal(
      setdiff(kla_scope_keys, whole_proteome_keys),
      paste(
        scope_exclusions$PXD,
        scope_exclusions$SampleGroup,
        sep = "__"
      )
    ),
  "Whole-proteome heatmap must be the exact-reference subset of the 37 Kla groups"
)
assert(
  nrow(scope_exclusions) == 4 &&
    setequal(
      scope_exclusions$PXD,
      c("PXD062720", "PXD063047", "PXD064038", "PXD075014")
    ),
  "Strict-reference exclusions must contain the four unsupported sample groups"
)
assert(
  !any(audit$PXD == "PXD037371"),
  "The three PXD037371 groups must not be present in the final heatmap audit"
)
assert(
  any(
    audit$PXD == "PXD014870" &
      audit$ReferencePXD == "PXD030304" &
      grepl(
        "ProCan-DepMapSanger_protein_matrix_6692_averaged.txt$",
        audit$ReferenceEvidenceLocator
      )
  ),
  "PXD014870 must be paired to the PXD030304 ordinary-proteome reference"
)
assert(
  all(
    audit$ReferencePXD[
      audit$PXD %in% c("PXD058534", "PXD078736")
    ] == "PXD072220"
  ) &&
    all(
      audit$WholeProteomeSampleCount[
        audit$PXD %in% c("PXD058534", "PXD078736")
      ] == 3
    ),
  "Both HK-2 Kla groups must use the three untreated PXD072220 controls"
)
huvec <- audit[audit$PXD == "PXD073311", ]
assert(
  nrow(huvec) == 1 &&
    huvec$ReferencePXD == "PXD073311" &&
    huvec$ReferenceSampleSubset == "A0h_1;A0h_2;A0h_3" &&
    huvec$WholeProteomeSampleCount == 3 &&
    huvec$WholeProteomeFeatureCount == 7709 &&
    grepl("report.pg_matrix.tsv$", huvec$WholeProteomeSource),
  "HUVEC must use the same-study PXD073311 A0h ordinary-proteome matrix"
)
assert(
  all(
    audit$WholeProteomeFeatureCount[audit$WholeProteomeQuantAvailable] > 0
  ),
  "Every available whole-proteome group must contain quantified features"
)
assert(
  nrow(ensembl_audit) == 33 &&
    sum(ensembl_audit$EnsemblMappingApplied) == 2 &&
    all(
      ensembl_audit$EnsemblMappingFile ==
        "reanalysis/config/ensembl_protein_to_uniprot_biomart.tsv"
    ),
  "Ensembl-to-UniProt mapping audit must cover the two PXD010154 reference groups"
)
assert(
  all(audit$WholeProteomeMappedRegulatorCount > 0),
  "Every final whole-proteome heatmap row should contain at least one mapped regulator after ID conversion"
)
assert(
  all(
    normalized$WholeProteomeRelativePercentile[
      !is.na(normalized$WholeProteomeRelativePercentile)
    ] >= 0 &
      normalized$WholeProteomeRelativePercentile[
        !is.na(normalized$WholeProteomeRelativePercentile)
      ] <= 100
  ),
  "Whole-proteome percentiles must be between 0 and 100"
)
assert(
  "IdentityMatchMode" %in% names(regulator_sample),
  "Regulator sample table must record its identity matching mode"
)
assert(
  all(regulator_sample$IdentityMatchMode == "UniProt_BaseAccession_only"),
  "Whole-proteome regulator matching must use UniProt BaseAccession only"
)
assert(
  all(regulator_sample$Signal >= 0),
  "Whole-proteome regulator signals must be non-negative"
)
sample_key <- paste(
  normalized$PXD,
  normalized$SampleGroup,
  normalized$RegulatorBaseAccession,
  sep = "__"
)
assert(
  !anyDuplicated(sample_key),
  "Normalized whole-proteome output must have one row per sample group and regulator accession"
)
raw <- read.csv(
  file.path(
    table_dir,
    "kla_regulator_whole_proteome_sample_level_long.csv"
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
raw_key <- paste(
  raw$PXD,
  raw$SampleGroup,
  raw$QuantSample,
  raw$CanonicalFeature,
  sep = "__"
)
assert(
  !anyDuplicated(raw_key),
  "Whole-proteome rank universe contains duplicated feature/sample rows"
)
percentile_variation <- aggregate(
  WholeProteomePercentile ~ PXD + SampleGroup + QuantSample,
  raw,
  function(values) length(unique(round(values, 8)))
)
assert(
  all(percentile_variation$WholeProteomePercentile > 1),
  "At least one whole-proteome sample collapsed to a single percentile/color"
)
assert(
  !any(grepl("Kla|Lactyl|Lac\\(K\\)|乳酸化", regulator_sample$Measurement,
             ignore.case = TRUE)),
  "Regulator sample table contains a Kla-enriched measurement"
)

pxd028488 <- audit[audit$PXD == "PXD028488", ]
assert(
  all(pxd028488$WholeProteomeMappedRegulatorCount > 0),
  "PXD028488 must retain mapped regulator detections after PEAKS accession parsing"
)
assert(
  !any(is.na(normalized$WholeProteomeRelativePercentile)),
  "All 33 strict-reference rows must have plottable regulator percentiles"
)
gcn5 <- heatmap_long[
  heatmap_long$RegulatorBaseAccession == "Q92830" &
    heatmap_long$Role == "Writer",
]
assert(
  nrow(gcn5) == 30 &&
    all(gcn5$GeneSymbol == "KAT2A") &&
    all(gcn5$RegulatorDisplayName == "GCN5 (KAT2A)"),
  "GCN5 must be displayed once per unique reference row as Writer Q92830/KAT2A"
)
assert(
  !anyDuplicated(
    paste(
      heatmap_long$ReferenceDisplayKey,
      heatmap_long$Role,
      heatmap_long$RegulatorBaseAccession,
      sep = "__"
    )
  ),
  "Heatmap display table must not duplicate a regulator within a role/reference row"
)
assert(
  all(heatmap_long$IdentityMatchMode == "UniProt_BaseAccession_only"),
  "Heatmap display rows must retain accession-only matching"
)
hippocampus <- audit[audit$PXD == "PXD050470", ]
assert(
  nrow(hippocampus) == 1 &&
    hippocampus$ReferencePXD == "PXD050470" &&
    hippocampus$WholeProteomeSampleCount == 3 &&
    hippocampus$WholeProteomeFeatureCount == 6082 &&
    hippocampus$WholeProteomeMappedRegulatorCount > 0,
  "PXD050470 must use same-study Table S4 for H072/H081/H0187"
)
assert(
  hippocampus$WholeProteomeSource ==
    "data/PXD050470/supplementary/prca2331-sup-0006-tables4.xlsx",
  "PXD050470 must use direct UniProt IDs from the same-study Table S4"
)
assert(
  file.exists(file.path(
    table_dir,
    "kla_regulator_whole_proteome_ensembl_mapping_audit.csv"
  )),
  "Ensembl mapping audit is missing"
)

pxd062720_path <- file.path(
  project_root,
  "data/PXD062720/search_results/extracted_pairing/txt/proteinGroups.txt"
)
pxd062720 <- read.delim(
  pxd062720_path,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = ""
)
pxd062720_valid <- pxd062720$Reverse != "+" &
  pxd062720$`Potential contaminant` != "+" &
  pxd062720$`Only identified by site` != "+"
pxd062720_kla <- pxd062720_valid &
  !is.na(pxd062720$`La (K) site IDs`) &
  nzchar(pxd062720$`La (K) site IDs`)
assert(
  sum(pxd062720_valid) == 2223 &&
    sum(pxd062720_kla) == 2060 &&
    !any(audit$PXD == "PXD062720"),
  "PXD062720 Kla-enriched proteinGroups must never enter the whole-proteome heatmap"
)

script_text <- paste(readLines(script_path, warn = FALSE), collapse = "\n")
assert(
  grepl("#FFFFFF", script_text, fixed = TRUE) &&
    grepl("#B2182B", script_text, fixed = TRUE),
  "Whole-proteome heatmap must use a white-to-warm palette"
)
assert(
  !grepl("#2166AC", script_text, fixed = TRUE) &&
    !grepl("#2C7FB8", script_text, fixed = TRUE),
  "Whole-proteome heatmap must not use a cool low-value palette"
)
assert(
  grepl("CategoryZh ~ Role", script_text, fixed = TRUE),
  "Whole-proteome heatmap must visibly facet rows by the four sample classes"
)

message("Whole-proteome regulator intensity tests passed.")
