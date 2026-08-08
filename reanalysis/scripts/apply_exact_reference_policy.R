#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args)) normalizePath(args[[1]]) else normalizePath(".")
config_path <- file.path(
  project_root,
  "reanalysis",
  "config",
  "lactylome_reference_pairing.csv"
)
exclusion_path <- file.path(
  project_root,
  "reanalysis",
  "config",
  "strict_reference_exclusions.csv"
)
audit_path <- file.path(
  project_root,
  "reanalysis",
  "results",
  "tables",
  "exact_reference_selection_audit.csv"
)

pairing <- read.csv(
  config_path,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  na.strings = c("", "NA")
)
if (!"IncludeInStrictReferenceAnalysis" %in% names(pairing)) {
  pairing$IncludeInStrictReferenceAnalysis <- pairing$IncludeInPairedAnalysis
}
exclusions <- read.csv(
  exclusion_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

replacements <- data.frame(
  LactylomePXD = c(
    "PXD028488", "PXD028737", "PXD063266",
    "PXD070007", "PXD070007", "PXD075377", "PXD075377"
  ),
  SampleGroup = c(
    "TALL-104", "HMC3", "PC-3M",
    "glioblastoma stem cells", "neural stem cells",
    "HCC", "adjacent liver"
  ),
  ReferenceStrategy = c(
    "same_study_non_enrichment",
    "same_study_conventional_proteome",
    "external_exact_cell_line",
    "same_biospecimen_non_ptm_proteome",
    "same_biospecimen_non_ptm_proteome",
    "external_exact_disease_tissue",
    "external_exact_adjacent_tissue"
  ),
  ReferencePXD = c(
    "PXD028488", "PXD028737", "PXD022005",
    "PXD069969", "PXD069969", "PXD065775", "PXD065775"
  ),
  ReferenceSampleSubset = c(
    "TALL-104 non-enrichment Sample 1;Sample 2;Sample 3",
    "H0;H24",
    "PC-3M heavy SILAC channel",
    "G2907;G3028;G3264;GSC23;MES28;RKI",
    "ENSA;HMP1",
    "CISs sheet;Non-rec1-4;Rec1-4",
    "ANTs sheet;Non-rec1-4;Rec1-4"
  ),
  ReferenceEvidenceLocator = c(
    paste0(
      "data/PXD028488/search_results/Nonenrichment-Search files/",
      "TALL-Nonenrichment-Search files/proteins.csv"
    ),
    paste0(
      "data/PXD028737/search_results/extracted_reference/",
      "txt/proteinGroups.txt"
    ),
    "data/PXD022005/search_results/txt_proteomics.zip",
    "data/PXD069969/search_results/SA206LQB1_Annotation.xlsx",
    "data/PXD069969/search_results/SA206LQB1_Annotation.xlsx",
    "data/PXD065775/search_results/20170330_01-24_patients_iTRAQ.xlsx",
    "data/PXD065775/search_results/20170330_01-24_patients_iTRAQ.xlsx"
  ),
  ReferenceProteinCount = c(4023L, 6132L, 8049L, 8374L, 7246L, 5188L, 4251L),
  ReferenceAcquisitionStatus = "downloaded_and_counted",
  MatchQuality = c(
    "exact_same_study",
    "exact_same_study",
    "exact_cell_line",
    "exact_same_biospecimen",
    "exact_same_biospecimen",
    "exact_disease_tissue",
    "exact_adjacent_tissue"
  ),
  Caveat = c(
    "同一PXD的TALL-104非富集PEAKS蛋白表；使用三次普通蛋白Area",
    "同一PXD的HMC3普通MaxQuant蛋白组；H0和H24均有LFQ强度",
    "PC-3M在PXD022005中为heavy SILAC通道；只使用Intensity H",
    "与Kla研究相同的6个GSC模型，普通非PTM LFQ蛋白组",
    "与Kla研究相同的2个NSC模型，普通非PTM LFQ蛋白组",
    "独立HCC临床普通iTRAQ蛋白组；材料和疾病状态匹配但不是同一患者",
    "独立HCC邻近肝普通iTRAQ蛋白组；邻近组织状态匹配但不是同一患者"
  ),
  IncludeInPairedAnalysis = TRUE,
  IncludeInStrictReferenceAnalysis = TRUE,
  stringsAsFactors = FALSE
)
replacements <- bind_rows(
  replacements,
  data.frame(
    LactylomePXD = "PXD050470",
    SampleGroup = "human hippocampus",
    ReferenceStrategy = "same_study_same_biospecimen_non_ptm_proteome",
    ReferencePXD = "PXD050470",
    ReferenceSampleSubset = "H072;H081;H0187",
    ReferenceEvidenceLocator =
      "data/PXD050470/supplementary/prca2331-sup-0006-tables4.xlsx",
    ReferenceProteinCount = 6082L,
    ReferenceAcquisitionStatus = "downloaded_and_counted",
    MatchQuality = "exact_same_biospecimen",
    Caveat = paste0(
      "Same-study Table S4 ordinary whole-proteome intensities for the same ",
      "three hippocampus samples; no CA1 substitution and no GeneSymbol conversion"
    ),
    IncludeInPairedAnalysis = TRUE,
    IncludeInStrictReferenceAnalysis = TRUE,
    stringsAsFactors = FALSE
  )
)

for (index in seq_len(nrow(replacements))) {
  replacement <- replacements[index, ]
  row <- pairing$LactylomePXD == replacement$LactylomePXD &
    pairing$SampleGroup == replacement$SampleGroup
  if (sum(row) != 1) {
    stop(
      "Expected one pairing row for ",
      replacement$LactylomePXD,
      " / ",
      replacement$SampleGroup
    )
  }
  for (column in setdiff(names(replacements), c("LactylomePXD", "SampleGroup"))) {
    pairing[row, column] <- replacement[[column]]
  }
}

for (index in seq_len(nrow(exclusions))) {
  exclusion <- exclusions[index, ]
  row <- pairing$LactylomePXD == exclusion$LactylomePXD &
    pairing$SampleGroup == exclusion$SampleGroup
  if (sum(row) != 1) {
    stop(
      "Expected one pairing row for exclusion ",
      exclusion$LactylomePXD,
      " / ",
      exclusion$SampleGroup
    )
  }
  pairing$ReferenceStrategy[row] <- "unresolved_reference"
  pairing$ReferencePXD[row] <- NA_character_
  pairing$ReferenceSampleSubset[row] <- NA_character_
  pairing$ReferenceEvidenceLocator[row] <- NA_character_
  pairing$ReferenceProteinCount[row] <- NA_real_
  pairing$ReferenceAcquisitionStatus[row] <- "not_selected"
  pairing$MatchQuality[row] <- "no_exact_reference_found"
  pairing$Caveat[row] <- exclusion$Reason
  pairing$IncludeInPairedAnalysis[row] <- TRUE
  pairing$IncludeInStrictReferenceAnalysis[row] <- FALSE
}
pairing$IncludeInStrictReferenceAnalysis <-
  pairing$IncludeInStrictReferenceAnalysis & pairing$IncludeInPairedAnalysis

write.csv(
  pairing,
  config_path,
  row.names = FALSE,
  na = ""
)

selected_audit <- replacements |>
  transmute(
    LactylomePXD,
    SampleGroup,
    Decision = "included_exact_reference",
    ReferencePXD,
    ReferenceSampleSubset,
    ReferenceEvidenceLocator,
    ReferenceProteinCount,
    IntensityRequirement = case_when(
      ReferencePXD == "PXD028488" ~ "PEAKS Area Sample 1-3",
      ReferencePXD == "PXD028737" ~ "LFQ intensity H0/H24",
      ReferencePXD == "PXD022005" ~ "Intensity H",
      ReferencePXD == "PXD069969" ~ "LFQ intensity by named cell model",
      ReferencePXD == "PXD065775" ~ "iTRAQ Non-rec1-4 and Rec1-4",
      LactylomePXD == "PXD050470" ~
        "Table S4 relative intensity H072 H081 H0187",
      TRUE ~ ""
    ),
    MatchQuality,
    DecisionReason = Caveat
  )
excluded_audit <- exclusions |>
  transmute(
    LactylomePXD,
    SampleGroup,
    Decision = "excluded_no_exact_quantitative_reference",
    ReferencePXD = "",
    ReferenceSampleSubset = "",
    ReferenceEvidenceLocator = "",
    ReferenceProteinCount = NA_real_,
    IntensityRequirement = "exact biological material plus per-protein intensity",
    MatchQuality = "no_exact_reference_found",
    DecisionReason = Reason
  )

dir.create(dirname(audit_path), recursive = TRUE, showWarnings = FALSE)
write.csv(
  bind_rows(selected_audit, excluded_audit),
  audit_path,
  row.names = FALSE,
  na = ""
)

message(
  "Applied exact-reference policy: ",
  nrow(replacements),
  " newly resolved rows; ",
  nrow(exclusions),
  " excluded rows."
)
