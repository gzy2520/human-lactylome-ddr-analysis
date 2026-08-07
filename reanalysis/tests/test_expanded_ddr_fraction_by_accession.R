#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args)) normalizePath(args[[1]]) else normalizePath(".")
table_dir <- file.path(project_root, "reanalysis", "results", "tables")
figure_dir <- file.path(project_root, "reanalysis", "results", "figures")

stats <- read.csv(
  file.path(
    table_dir,
    "cell_type_kla_vs_reference_ddr_statistics_accession_only.csv"
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
audit <- read.csv(
  file.path(
    table_dir,
    "cell_type_kla_vs_reference_ddr_accession_only_audit.csv"
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
reference_members <- read.csv(
  file.path(
    project_root,
    "reanalysis", "intermediate", "expanded_ddr_by_accession",
    "reference_proteins_by_sample_group.csv"
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
evidence_audit <- read.csv(
  file.path(
    table_dir,
    "cell_type_kla_ddr_lactylation_evidence_audit.csv"
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
stopifnot(file.exists(file.path(
  table_dir,
  "cell_type_kla_vs_reference_ddr_statistics_accession_only_zh.csv"
)))

stopifnot(nrow(stats) == 37)
stopifnot(sum(audit$Included) == 37)
stopifnot(sum(!audit$Included) == 3)
stopifnot(all(unique(audit$PXD[!audit$Included]) == "PXD037371"))
stopifnot(identical(
  unique(stats$Category),
  c("normal_tissue", "cancer_tissue", "normal_cells", "cancer_cells")
))
stopifnot(identical(
  as.integer(table(factor(
    stats$Category,
    levels = c("normal_tissue", "cancer_tissue", "normal_cells", "cancer_cells")
  ))),
  c(10L, 3L, 10L, 14L)
))
stopifnot(all(stats$MatchMode %in% c(
  "BaseAccession_only",
  "BaseAccession_after_reviewed_UniProt_symbol_conversion"
)))
stopifnot(all(stats$SymbolFallbackCount == 0))
stopifnot(nrow(evidence_audit) == nrow(stats))
stopifnot(all(evidence_audit$DirectLactylationEvidence))
stopifnot(sum(evidence_audit$GeneSymbolFallbackCount) == 0)
stopifnot(sum(evidence_audit$Pre2019ReviewRequired) == 1)
stopifnot(all(
  evidence_audit$PXD[evidence_audit$Pre2019ReviewRequired] == "PXD014870"
))
stopifnot(all(grepl(
  "pan anti-Kla",
  evidence_audit$EnrichmentOrValidationMethod[
    evidence_audit$Pre2019ReviewRequired
  ],
  fixed = TRUE
)))
stopifnot(all(grepl(
  "^include:",
  evidence_audit$Pre2019Eligibility[evidence_audit$Pre2019ReviewRequired]
)))
stopifnot(file.exists(file.path(
  table_dir,
  "pre_2019_lactylation_dataset_review.csv"
)))
all_source_files <- unique(unlist(strsplit(
  paste(evidence_audit$KlaEvidenceFile, collapse = ";"),
  ";",
  fixed = TRUE
)))
stopifnot(all(file.exists(file.path(project_root, all_source_files))))
stopifnot(all(stats$KlaDdrProteinCount <= stats$KlaProteinCount))
stopifnot(all(stats$ReferenceDdrProteinCount <= stats$ReferenceProteinCount))
stopifnot(all(stats$KlaDdrFraction >= 0 & stats$KlaDdrFraction <= 1))
stopifnot(all(stats$ReferenceDdrFraction >= 0 & stats$ReferenceDdrFraction <= 1))
stopifnot(any(stats$PXD == "PXD028737"))
stopifnot(any(stats$PXD == "PXD073311"))
stopifnot(any(stats$PXD == "PXD075014"))
stopifnot(
  stats$KlaProteinCount[
    stats$PXD == "PXD046800" &
      stats$SampleGroup == "hypertrophic scar"
  ] == 402
)
stopifnot(
  stats$KlaProteinCount[
    stats$PXD == "PXD046800" &
      stats$SampleGroup == "adjacent skin"
  ] == 527
)
stopifnot(stats$KlaProteinCount[stats$PXD == "PXD066351"] == 2398)
stopifnot(any(
  grepl(
    "^A0A",
    read.csv(
      file.path(
        project_root,
        "reanalysis", "intermediate", "expanded_ddr_by_accession",
        "kla_proteins_by_sample_group.csv"
      ),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )$BaseAccession
  )
))
stopifnot(any(stats$PXD == "PXD050470"))
stopifnot(
  stats$ReferenceProteinCount[stats$PXD == "PXD050470"] == 2105
)
stopifnot(
  stats$ReferenceDdrProteinCount[stats$PXD == "PXD050470"] == 83
)
stopifnot(
  stats$MatchMode[stats$PXD == "PXD050470"] ==
    "BaseAccession_after_reviewed_UniProt_symbol_conversion"
)
hippocampus_mapping_path <- file.path(
  table_dir,
  "PXD043880_hippocampus_symbol_to_reviewed_uniprot_mapping_audit.csv"
)
stopifnot(file.exists(hippocampus_mapping_path))
hippocampus_mapping <- read.csv(
  hippocampus_mapping_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
stopifnot(nrow(hippocampus_mapping) > 2000)
stopifnot(
  length(unique(
    hippocampus_mapping$BaseAccession[
      nzchar(hippocampus_mapping$BaseAccession)
    ]
  )) == 2105
)
ensembl_rows <- reference_members$IdentifierType == "ENSEMBLPROT"
stopifnot(sum(ensembl_rows) > 10000)
stopifnot(
  mean(nzchar(reference_members$MappedBaseAccessions[ensembl_rows])) >= 0.85
)

huvec_stats <- stats[
  stats$PXD == "PXD073311" &
    stats$SampleGroup == "HUVEC control and Pg infection",
]
stopifnot(nrow(huvec_stats) == 1)
stopifnot(huvec_stats$ReferencePXD == "PXD073311")
stopifnot(huvec_stats$ReferenceProteinCount == 7794)
stopifnot(huvec_stats$ReferenceDdrProteinCount == 512)
stopifnot(grepl(
  "report.pg_matrix.tsv$",
  huvec_stats$ReferenceEvidenceFile
))
huvec_members <- reference_members[
  reference_members$PXD == "PXD073311" &
    reference_members$SampleGroup == "HUVEC control and Pg infection",
]
stopifnot(nrow(huvec_members) == 7794)
stopifnot(all(huvec_members$ReferencePXD == "PXD073311"))
stopifnot(sum(huvec_members$IsDdr) == 512)

huvec_matrix_path <- file.path(project_root, huvec_stats$ReferenceEvidenceFile)
huvec_matrix <- read.delim(
  huvec_matrix_path,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = ""
)
huvec_a0h_columns <- grep("A0h_[123]\\.raw$", names(huvec_matrix), value = TRUE)
stopifnot(length(huvec_a0h_columns) == 3)
stopifnot(!any(grepl("A6h", huvec_a0h_columns, fixed = TRUE)))
huvec_detected <- rowSums(
  sapply(huvec_matrix[huvec_a0h_columns], function(values) {
    values <- suppressWarnings(as.numeric(values))
    !is.na(values) & values > 0
  }),
  na.rm = TRUE
) > 0
expected_huvec_ids <- unique(trimws(unlist(strsplit(
  as.character(huvec_matrix$Protein.Group[huvec_detected]),
  "[;,]"
))))
expected_huvec_ids <- sub("-[0-9]+$", "", expected_huvec_ids)
expected_huvec_ids <- sort(expected_huvec_ids[nzchar(expected_huvec_ids)])
stopifnot(setequal(expected_huvec_ids, huvec_members$SourceProteinID))

archive_stats_path <- file.path(
  project_root,
  "archive", "reanalysis_2026-08-07_pre_37group_reference_fix",
  "results", "tables",
  "cell_type_kla_vs_reference_ddr_statistics_accession_only.csv"
)
if (file.exists(archive_stats_path)) {
  archived <- read.csv(
    archive_stats_path,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  archived <- archived[
    match(
      paste(stats$PXD, stats$SampleGroup, sep = "__"),
      paste(archived$PXD, archived$SampleGroup, sep = "__")
    ),
  ]
  stopifnot(all(stats$KlaProteinCount == archived$KlaProteinCount))
  stopifnot(all(stats$KlaDdrProteinCount == archived$KlaDdrProteinCount))
}

for (name in c(
  "cell_type_kla_vs_reference_ddr_fraction_accession_only.png",
  "cell_type_kla_vs_reference_ddr_fraction_accession_only.pdf",
  "cell_type_kla_vs_reference_ddr_fraction.png",
  "cell_type_kla_vs_reference_ddr_fraction.pdf",
  "cell_type_kla_vs_reference_ddr_fraction_accession_only_en.png",
  "cell_type_kla_vs_reference_ddr_fraction_accession_only_en.pdf",
  "cell_type_kla_vs_reference_ddr_fraction_en.png",
  "cell_type_kla_vs_reference_ddr_fraction_en.pdf"
)) {
  path <- file.path(figure_dir, name)
  stopifnot(file.exists(path))
  stopifnot(file.info(path)$size > 10000)
}

cat("Expanded accession-only DDR fraction tests passed.\n")
