#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) normalizePath(args[[1]]) else normalizePath(".")
config_dir <- file.path(project_root, "reanalysis", "config")
table_dir <- file.path(project_root, "reanalysis", "results", "tables")
report_dir <- file.path(project_root, "reanalysis", "reports")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)

pairing <- read.csv(
  file.path(config_dir, "lactylome_reference_pairing.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE,
  na.strings = c("", "NA")
)
if (!"IncludeInStrictReferenceAnalysis" %in% names(pairing)) {
  pairing$IncludeInStrictReferenceAnalysis <- pairing$IncludeInPairedAnalysis
}
pairing_config_columns <- names(pairing)
decisions <- read.csv(
  file.path(config_dir, "lactylome_dataset_decisions.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
inventory <- read.csv(
  file.path(table_dir, "human_lactylome_mass_spectrometry_inventory.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
remote <- read.csv(
  file.path(table_dir, "lactylome_pair_remote_file_sizes.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
download_path <- file.path(table_dir, "lactylome_pair_download_manifest.csv")
downloads <- if (file.exists(download_path)) {
  read.csv(download_path, stringsAsFactors = FALSE, check.names = FALSE)
} else {
  data.frame(PXD = character(), SizeBytes = numeric(), stringsAsFactors = FALSE)
}
healthy_config <- read.csv(
  file.path(config_dir, "healthy_tissue_reference_files.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
healthy_special_path <- file.path(
  config_dir,
  "healthy_special_reference_catalog.csv"
)
healthy_special <- if (file.exists(healthy_special_path)) {
  read.csv(
    healthy_special_path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
} else {
  data.frame(
    TissueKey = character(),
    DisplayName = character(),
    PXD = character(),
    FileName = character(),
    ProteinGroupsPath = character(),
    ProteinCount = integer(),
    Status = character(),
    SourceURL = character(),
    MatchQuality = character(),
    Caveat = character(),
    CountBasis = character(),
    stringsAsFactors = FALSE
  )
}
healthy_manifest_path <- file.path(
  table_dir,
  "healthy_tissue_reference_acquisition_manifest.csv"
)
healthy_manifest <- if (file.exists(healthy_manifest_path)) {
  read.csv(healthy_manifest_path, stringsAsFactors = FALSE, check.names = FALSE)
} else {
  data.frame(
    TissueKey = character(),
    ProteinGroupsPath = character(),
    ProteinCount = integer(),
    Status = character(),
    stringsAsFactors = FALSE
  )
}
strict_reference_path <- file.path(
  config_dir,
  "strict_reference_exclusions.csv"
)
strict_reference_exclusions <- read.csv(
  strict_reference_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

primary_path <- file.path(
  project_root,
  "reanalysis",
  "intermediate",
  "kla_by_dataset",
  "all_primary_sample_level_kla_sites.csv"
)
primary <- read.csv(primary_path, stringsAsFactors = FALSE, check.names = FALSE)

base_accession <- function(values) {
  values <- sub("^.*\\|([^|]+)\\|.*$", "\\1", values)
  sub("-[0-9]+$", "", values)
}

split_accessions <- function(values) {
  values <- unlist(strsplit(as.character(values), "[;,]"))
  values <- trimws(values)
  values <- values[nzchar(values) & !is.na(values)]
  unique(base_accession(values))
}

finite_value <- function(values) {
  result <- suppressWarnings(as.numeric(values))
  !is.na(result) & is.finite(result)
}

count_primary <- function(pxd, sample_group) {
  aliases <- c(
    "TALL-104" = "T-ALL",
    "human hippocampus" = "Human hippocampus",
    "RKO WT and GSK3B KO" = "RKO",
    "HK-2 control and mannitol" = "HK-2"
  )
  cell_type <- if (sample_group %in% names(aliases)) aliases[[sample_group]] else sample_group
  match <- primary[
    primary$PXD == pxd &
      tolower(primary$CellOrTissueType) == tolower(cell_type) &
      primary$PrimaryIncluded %in% c(TRUE, "TRUE", "True", 1, "1"),
  ]
  if (!nrow(match)) {
    return(NA_integer_)
  }
  length(unique(match$BaseAccession[nzchar(match$BaseAccession)]))
}

count_matrix_sample <- function(sample_subset) {
  matrix_path <- file.path(
    project_root,
    "data",
    "PXD030304",
    "search_results",
    "ProCan-DepMapSanger_protein_matrix_6692_averaged.txt"
  )
  if (!file.exists(matrix_path)) {
    return(NA_integer_)
  }
  sidm <- str_extract(sample_subset, "SIDM[0-9]+")
  if (is.na(sidm)) {
    return(NA_integer_)
  }
  data <- read.delim(matrix_path, check.names = FALSE, stringsAsFactors = FALSE)
  row <- data[grepl(sidm, data[[1]], fixed = TRUE), , drop = FALSE]
  if (nrow(row) != 1) {
    return(NA_integer_)
  }
  sum(!is.na(unlist(row[1, -1], use.names = FALSE)))
}

count_pxd066054 <- function(group, lactylome = FALSE) {
  if (lactylome) {
    path <- file.path(
      project_root,
      "data/PXD066054/search_results/extracted/PLa/XB08700B1DPLa_-PTMSiteReport.tsv"
    )
    if (!file.exists(path)) {
      return(NA_integer_)
    }
    data <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
    data <- data[data$PTM.ModificationTitle == "L-Lac(K)",]
    prefix <- if (group == "BPH") "^NAT" else "^PCa"
    return(length(unique(base_accession(
      data$PTM.ProteinId[grepl(prefix, data$R.Condition)]
    ))))
  }
  path <- file.path(
    project_root,
    "data/PXD066054/search_results/extracted/DA/Protein_Quant.tsv"
  )
  if (!file.exists(path)) {
    return(NA_integer_)
  }
  data <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  columns <- grep(if (group == "BPH") "NAT" else "PCa", names(data))
  keep <- rowSums(sapply(data[columns], finite_value)) > 0
  length(split_accessions(data[[1]][keep]))
}

count_pxd063047 <- function(group) {
  path <- file.path(
    project_root,
    "data/PXD063047/search_results/extracted/combined/txt/La (K)Sites.txt"
  )
  if (!file.exists(path)) {
    return(NA_integer_)
  }
  data <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  prefix <- if (group == "normal pregnancy placenta") "Localization prob Con_" else "Localization prob PE_"
  columns <- grep(paste0("^", prefix), names(data))
  keep <- rowSums(sapply(data[columns], function(x) suppressWarnings(as.numeric(x)) > 0), na.rm = TRUE) > 0
  keep <- keep & data$Reverse != "+" & data$`Potential contaminant` != "+" & !is.na(data$id)
  length(split_accessions(data$Proteins[keep]))
}

count_pxd075377 <- function(group) {
  path <- file.path(
    project_root,
    "data/PXD075377/search_results/extracted/2-Basic_analysis/MS_identified_information.txt"
  )
  if (!file.exists(path)) {
    return(NA_integer_)
  }
  data <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  column <- if (group == "HCC") "Intensity HCC" else "Intensity Control"
  keep <- finite_value(data[[column]]) & suppressWarnings(as.numeric(data[[column]])) > 0
  length(unique(base_accession(data$`Protein accession`[keep])))
}

count_simple_protein_table <- function(path) {
  if (!file.exists(path)) {
    return(NA_integer_)
  }
  delimiter <- if (grepl("\\.csv$", path, ignore.case = TRUE)) "," else "\t"
  data <- tryCatch(
    read.delim(
      path,
      sep = delimiter,
      check.names = FALSE,
      stringsAsFactors = FALSE,
      quote = "",
      comment.char = ""
    ),
    error = function(e) NULL
  )
  if (is.null(data) || !nrow(data)) {
    return(NA_integer_)
  }
  candidates <- c(
    "PG.ProteinGroups", "Protein IDs", "Protein.Group", "Protein Accessions",
    "Accession", "Protein accession", "Protein"
  )
  column <- intersect(candidates, names(data))
  if (!length(column)) {
    column <- names(data)[[1]]
  } else {
    column <- column[[1]]
  }
  length(split_accessions(data[[column]]))
}

count_maxquant_site_table <- function(path, group_pattern = NULL) {
  if (!file.exists(path)) {
    return(NA_integer_)
  }
  data <- tryCatch(
    read.delim(
      path,
      check.names = FALSE,
      stringsAsFactors = FALSE,
      quote = "",
      comment.char = ""
    ),
    error = function(e) NULL
  )
  if (is.null(data) || !nrow(data)) {
    return(NA_integer_)
  }
  keep <- rep(TRUE, nrow(data))
  if ("Reverse" %in% names(data)) {
    keep <- keep & data$Reverse != "+"
  }
  if ("Potential contaminant" %in% names(data)) {
    keep <- keep & data$`Potential contaminant` != "+"
  }
  if ("id" %in% names(data)) {
    keep <- keep & !is.na(data$id)
  }
  if (is.null(group_pattern)) {
    if ("Localization prob" %in% names(data)) {
      keep <- keep & suppressWarnings(as.numeric(data$`Localization prob`)) > 0
    }
  } else {
    columns <- grep(
      paste0("^Localization prob .*", group_pattern),
      names(data),
      ignore.case = TRUE
    )
    if (length(columns)) {
      localized <- rowSums(
        sapply(
          data[columns],
          function(values) suppressWarnings(as.numeric(values)) > 0
        ),
        na.rm = TRUE
      ) > 0
      keep <- keep & localized
    }
  }
  protein_column <- intersect(c("Proteins", "Protein", "Leading proteins"), names(data))
  if (!length(protein_column)) {
    return(NA_integer_)
  }
  length(split_accessions(data[[protein_column[[1]]]][keep]))
}

find_files <- function(pxd, pattern) {
  root <- file.path(project_root, "data", pxd, "search_results")
  if (!dir.exists(root)) {
    return(character())
  }
  list.files(root, pattern = pattern, recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
}

count_first_table <- function(pxd, pattern) {
  files <- find_files(pxd, pattern)
  if (!length(files)) {
    return(NA_integer_)
  }
  for (file in files) {
    value <- count_simple_protein_table(file)
    if (!is.na(value) && value > 0) {
      return(value)
    }
  }
  NA_integer_
}

count_pxd046800 <- function(group, lactylome = FALSE) {
  path <- file.path(
    project_root,
    "data",
    "PXD046800",
    "search_results",
    if (lactylome) {
      "HFX2_LFQ_QB001_Lacty_Proteins.txt"
    } else {
      "HFX2_LFQ_QB002_Proteins.txt"
    }
  )
  if (!file.exists(path)) {
    return(NA_integer_)
  }
  data <- read.delim(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    quote = "\"",
    comment.char = ""
  )
  if (!"Accession" %in% names(data)) {
    return(NA_integer_)
  }
  prefix <- if (group == "hypertrophic scar") "HSP" else "NSP"
  columns <- grep(paste0("Found in Sample:.*", prefix), names(data), ignore.case = TRUE)
  if (!length(columns)) {
    return(length(unique(base_accession(data$Accession))))
  }
  keep <- rowSums(sapply(data[columns], function(values) {
    !is.na(values) & nzchar(values) & values != "Not Found"
  })) > 0
  length(unique(base_accession(data$Accession[keep])))
}

count_spectronaut_report <- function(path, lactylome = FALSE) {
  if (!file.exists(path) || !requireNamespace("data.table", quietly = TRUE)) {
    return(NA_integer_)
  }
  header <- names(
    data.table::fread(
      path,
      nrows = 0,
      showProgress = FALSE,
      data.table = FALSE
    )
  )
  required <- c("Protein.Group", "PG.Q.Value")
  if (lactylome) {
    required <- c(
      required,
      "Modified.Sequence",
      "Q.Value"
    )
    if ("PTM.Site.Confidence" %in% header) {
      required <- c(required, "PTM.Site.Confidence")
    }
  }
  data <- tryCatch(
    data.table::fread(
      path,
      select = required,
      showProgress = FALSE,
      data.table = FALSE
    ),
    error = function(e) NULL
  )
  if (is.null(data) || !nrow(data)) {
    return(NA_integer_)
  }
  keep <- nzchar(data$Protein.Group) &
    suppressWarnings(as.numeric(data$PG.Q.Value)) <= 0.01
  if (lactylome) {
    keep <- keep &
      grepl("K\\(UniMod:378\\)", data$Modified.Sequence) &
      suppressWarnings(as.numeric(data$Q.Value)) <= 0.01
    if ("PTM.Site.Confidence" %in% names(data)) {
      keep <- keep &
        suppressWarnings(as.numeric(data$PTM.Site.Confidence)) > 0
    }
  }
  length(split_accessions(data$Protein.Group[keep]))
}

count_pxd065831_reference <- function(group) {
  path <- file.path(
    project_root,
    "data/PXD065831/search_results/extracted_pairing/",
    "YAS202408210011-1/YAS202408210011-1.pg_matrix.tsv"
  )
  if (!file.exists(path)) {
    return(NA_integer_)
  }
  data <- data.table::fread(
    path,
    showProgress = FALSE,
    data.table = FALSE
  )
  sample_pattern <- if (group == "normal endometrium") {
    "-NE-"
  } else {
    "-HEEC-|-L-MEEC-"
  }
  columns <- grep(sample_pattern, names(data))
  if (!length(columns)) {
    return(NA_integer_)
  }
  numeric_values <- sapply(
    data[columns],
    function(values) suppressWarnings(as.numeric(values))
  )
  keep <- grepl("_HUMAN", data$Protein.Names, fixed = TRUE) &
    rowSums(numeric_values > 0, na.rm = TRUE) > 0
  sum(keep)
}

count_pxd070007 <- function(group) {
  path <- file.path(
    project_root,
    "data/PXD070007/search_results/SA206LPLaB1_Annotation.xlsx"
  )
  if (!file.exists(path) || !requireNamespace("readxl", quietly = TRUE)) {
    return(NA_integer_)
  }
  data <- readxl::read_excel(path, sheet = "Annotation_Combine")
  sample_columns <- if (group == "glioblastoma stem cells") {
    c("G2907", "G3028", "G3264", "GSC23", "MES28", "RKI")
  } else {
    c("ENSA", "HMP1")
  }
  detected <- sapply(
    data[sample_columns],
    function(values) {
      numeric_values <- suppressWarnings(as.numeric(values))
      !is.na(numeric_values) & numeric_values > 0
    }
  )
  localized <- suppressWarnings(
    as.numeric(data$`Localization probability`)
  ) > 0
  keep <- localized & rowSums(detected, na.rm = TRUE) > 0
  length(unique(base_accession(
    data$`Protein accession`[keep & nzchar(data$`Protein accession`)]
  )))
}

count_pxd064912 <- function() {
  path <- file.path(
    project_root,
    "data/PXD064912/supplementary/europepmc/mmc1.xlsx"
  )
  if (!file.exists(path) || !requireNamespace("readxl", quietly = TRUE)) {
    return(NA_integer_)
  }
  data <- readxl::read_excel(path, skip = 1)
  probability_columns <- grep("^PTM.SiteProbability", names(data))
  probabilities <- sapply(
    data[probability_columns],
    function(values) suppressWarnings(as.numeric(values))
  )
  keep <- tolower(data$PTM.ModificationTitle) == "lactylation" &
    data$PTM.SiteAA == "K" &
    rowSums(probabilities > 0, na.rm = TRUE) > 0
  length(unique(base_accession(
    data$PTM.ProteinId[keep & nzchar(data$PTM.ProteinId)]
  )))
}

count_pxd073311_reference <- function() {
  path <- file.path(
    project_root,
    paste0(
      "data/PXD073311/search_results/extracted_pairing/",
      "IPX0015307001_Database_search_result/Database_search_result/",
      "report.pg_matrix.tsv"
    )
  )
  if (!file.exists(path)) return(NA_integer_)
  data <- read.delim(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = ""
  )
  sample_columns <- grep("A0h_[123]\\.raw$", names(data), value = TRUE)
  if (length(sample_columns) != 3 || !"Protein.Group" %in% names(data)) {
    return(NA_integer_)
  }
  detected <- rowSums(
    sapply(data[sample_columns], function(values) {
      values <- suppressWarnings(as.numeric(values))
      !is.na(values) & values > 0
    }),
    na.rm = TRUE
  ) > 0
  length(unique(split_accessions(data$Protein.Group[detected])))
}

known_lactylome_counts <- c(
  PXD028737 = 1270,
  PXD036307 = 476,
  PXD054919 = 1220,
  PXD063266 = 379,
  PXD066351 = 2386,
  PXD070007 = NA,
  PXD073311 = 1881,
  PXD075014 = 521
)

pairing$LactylomeProteinCount <- NA_integer_
pairing$LactylomeProteinCountBasis <- "not_yet_counted"

for (i in seq_len(nrow(pairing))) {
  pxd <- pairing$LactylomePXD[[i]]
  group <- pairing$SampleGroup[[i]]
  value <- count_primary(pxd, group)
  basis <- "sample-level primary Kla long table"
  if (pxd == "PXD063047") {
    value <- count_pxd063047(group)
    basis <- "group-specific positive localization columns in La (K)Sites.txt"
  } else if (pxd == "PXD066054") {
    value <- count_pxd066054(group, lactylome = TRUE)
    basis <- "unique L-Lac(K) PTM.ProteinId by Spectronaut condition"
  } else if (pxd == "PXD075377") {
    value <- count_pxd075377(group)
    basis <- "unique accession with positive group intensity"
  } else if (pxd == "PXD046800") {
    value <- count_pxd046800(group, lactylome = TRUE)
    basis <- "author Proteome Discoverer lactylome protein table by sample group"
  } else if (pxd == "PXD050147") {
    value <- count_maxquant_site_table(
      file.path(project_root, "data/PXD050147/search_results/Lactyl_K_Sites.txt")
    )
    basis <- "MaxQuant Lactyl_K_Sites unique protein union"
  } else if (pxd == "PXD055230") {
    path <- file.path(
      project_root,
      "data/PXD055230/search_results/LaIP_HSV1_DIA.tsv"
    )
    value <- count_spectronaut_report(path, lactylome = TRUE)
    basis <- "Spectronaut HSV-1 LaIP report with K(UniMod:378), q <= 0.01 and site confidence > 0"
    if (!is.na(value)) {
      pairing$LactylomeEvidenceLocator[[i]] <-
        "data/PXD055230/search_results/LaIP_HSV1_DIA.tsv"
      pairing$LactylomeAcquisitionStatus[[i]] <- "downloaded_and_qc_passed"
    }
  } else if (pxd == "PXD057709") {
    path <- file.path(
      project_root,
      "data/PXD057709/search_results/LaIP_report.tsv"
    )
    value <- count_spectronaut_report(path, lactylome = TRUE)
    basis <- "Spectronaut HCMV LaIP report with K(UniMod:378) and q <= 0.01; site confidence field unavailable"
    if (!is.na(value)) {
      pairing$LactylomeEvidenceLocator[[i]] <-
        "data/PXD057709/search_results/LaIP_report.tsv"
      pairing$LactylomeAcquisitionStatus[[i]] <- "downloaded_and_qc_passed"
    }
  } else if (pxd == "PXD064912") {
    value <- count_pxd064912()
    basis <- "author sperm supplementary Kla site table, localization probability > 0"
    if (!is.na(value)) {
      pairing$LactylomeEvidenceLocator[[i]] <-
        "data/PXD064912/supplementary/europepmc/mmc1.xlsx"
      pairing$LactylomeAcquisitionStatus[[i]] <- "downloaded_and_qc_passed"
    }
  } else if (pxd == "PXD070007") {
    value <- count_pxd070007(group)
    basis <- if (group == "glioblastoma stem cells") {
      "author annotation table, six GSC models, localization probability > 0"
    } else {
      "author annotation table, two NSC models, localization probability > 0"
    }
    if (!is.na(value)) {
      pairing$LactylomeEvidenceLocator[[i]] <-
        "data/PXD070007/search_results/SA206LPLaB1_Annotation.xlsx"
      pairing$LactylomeAcquisitionStatus[[i]] <- "downloaded_and_qc_passed"
    }
  } else if (is.na(value) && pxd %in% c("PXD033146", "PXD037371", "PXD058534", "PXD062720", "PXD064038")) {
    site_files <- find_files(pxd, "(La.*Sites|Lactyl.*Sites)\\.txt$")
    if (length(site_files)) {
      value <- count_maxquant_site_table(site_files[[1]])
      basis <- "downloaded search-result site table unique protein union"
    }
  } else if (is.na(value) && pxd %in% names(known_lactylome_counts)) {
    value <- known_lactylome_counts[[pxd]]
    basis <- "author or acquisition QC dataset-level count"
  }
  pairing$LactylomeProteinCount[[i]] <- value
  pairing$LactylomeProteinCountBasis[[i]] <- if (is.na(value)) "not_yet_counted" else basis

  if (
    !is.na(pairing$ReferencePXD[[i]]) &&
      pairing$ReferencePXD[[i]] == "PXD030304" &&
      is.na(pairing$ReferenceProteinCount[[i]])
  ) {
    pairing$ReferenceProteinCount[[i]] <- count_matrix_sample(pairing$ReferenceSampleSubset[[i]])
    if (!is.na(pairing$ReferenceProteinCount[[i]])) {
      pairing$ReferenceAcquisitionStatus[[i]] <- "downloaded_and_counted"
    }
  }
  if (pxd == "PXD066054") {
    pairing$ReferenceProteinCount[[i]] <- count_pxd066054(group, lactylome = FALSE)
    pairing$ReferenceAcquisitionStatus[[i]] <- "downloaded_and_counted"
  }
  if (pxd == "PXD046800") {
    pairing$ReferenceProteinCount[[i]] <- count_pxd046800(group, lactylome = FALSE)
    if (!is.na(pairing$ReferenceProteinCount[[i]])) {
      pairing$ReferenceAcquisitionStatus[[i]] <- "downloaded_and_counted"
    }
  }
  if (pxd == "PXD055230") {
    path <- file.path(
      project_root,
      "data/PXD055230/search_results/WP_HSV1_DIA.tsv"
    )
    pairing$ReferenceProteinCount[[i]] <-
      count_spectronaut_report(path, lactylome = FALSE)
    if (!is.na(pairing$ReferenceProteinCount[[i]])) {
      pairing$ReferenceEvidenceLocator[[i]] <-
        "data/PXD055230/search_results/WP_HSV1_DIA.tsv"
      pairing$ReferenceAcquisitionStatus[[i]] <- "downloaded_and_counted"
    }
  }
  if (pxd == "PXD057709") {
    path <- file.path(
      project_root,
      "data/PXD057709/search_results/WP_report.tsv"
    )
    pairing$ReferenceProteinCount[[i]] <-
      count_spectronaut_report(path, lactylome = FALSE)
    if (!is.na(pairing$ReferenceProteinCount[[i]])) {
      pairing$ReferenceEvidenceLocator[[i]] <-
        "data/PXD057709/search_results/WP_report.tsv"
      pairing$ReferenceAcquisitionStatus[[i]] <- "downloaded_and_counted"
    }
  }
  if (pxd == "PXD065831") {
    pairing$ReferenceProteinCount[[i]] <-
      count_pxd065831_reference(group)
    if (!is.na(pairing$ReferenceProteinCount[[i]])) {
      pairing$ReferenceEvidenceLocator[[i]] <-
        paste0(
          "data/PXD065831/search_results/extracted_pairing/",
          "YAS202408210011-1/YAS202408210011-1.pg_matrix.tsv"
        )
      pairing$ReferenceAcquisitionStatus[[i]] <- "downloaded_and_counted"
    }
  }
  if (pxd == "PXD073311") {
    pairing$ReferenceProteinCount[[i]] <- count_pxd073311_reference()
    if (!is.na(pairing$ReferenceProteinCount[[i]])) {
      pairing$ReferenceEvidenceLocator[[i]] <-
        paste0(
          "data/PXD073311/search_results/extracted_pairing/",
          "IPX0015307001_Database_search_result/Database_search_result/",
          "report.pg_matrix.tsv"
        )
      pairing$ReferenceSampleSubset[[i]] <- "A0h_1;A0h_2;A0h_3"
      pairing$ReferenceAcquisitionStatus[[i]] <- "downloaded_and_counted"
    }
  }
}

ordinary_paths <- c(
  PXD046800 = "data/PXD046800/search_results/HFX2_LFQ_QB002_Proteins.txt",
  PXD050147 = "data/PXD050147/search_results/SIRT_proteinGroups.txt",
  PXD066351 = "data/PXD066351/search_results/XB01472B1DA-Protein_Quant.tsv"
)
for (pxd in names(ordinary_paths)) {
  count <- count_simple_protein_table(file.path(project_root, ordinary_paths[[pxd]]))
  rows <- !is.na(pairing$ReferencePXD) &
    pairing$ReferencePXD == pxd &
    is.na(pairing$ReferenceProteinCount)
  if (!is.na(count) && any(rows)) {
    pairing$ReferenceProteinCount[rows] <- count
    pairing$ReferenceAcquisitionStatus[rows] <- "downloaded_and_counted"
  }
}

extracted_reference_patterns <- c(
  PXD033146 = "(proteinGroups|Proteins|Protein_Quant).*(txt|tsv)$",
  PXD052772 = "proteinGroups\\.txt$",
  PXD062720 = "proteinGroups\\.txt$"
)
for (pxd in names(extracted_reference_patterns)) {
  count <- count_first_table(pxd, extracted_reference_patterns[[pxd]])
  rows <- !is.na(pairing$ReferencePXD) &
    pairing$ReferencePXD == pxd &
    is.na(pairing$ReferenceProteinCount)
  if (!is.na(count) && any(rows)) {
    pairing$ReferenceProteinCount[rows] <- count
    pairing$ReferenceAcquisitionStatus[rows] <- "downloaded_and_counted"
  }
}

infer_tissue_key <- function(sample_group, biological_material) {
  value <- tolower(paste(sample_group, biological_material))
  rules <- list(
    lung = "lung|a549|luad",
    placenta = "placenta",
    liver = "liver|hcc|hepg2|cirrhos",
    stomach = "gastric|\\bags\\b",
    brain = "brain|microgl|hmc3|glioblastoma|neural stem",
    urinary_bladder = "bladder",
    esophagus = "esoph|escc|kyse",
    heart = "heart|cardiomy|ac16",
    endometrium = "endometri",
    colon = "colorectal|colon|hct116|rko",
    kidney = "kidney|renal|hk-2|achn",
    prostate = "prostate|pc-3|\\bbph\\b",
    breast = "breast|mammary|mcf7|mcf10a|mda-mb-468|t-47d",
    lymphoid = "tall|t-all|lymphoblastic|leukemia",
    tendon = "tendon",
    skin = "skin|scar|fibroblast",
    oral_mucosa = "oral|oscc",
    vascular = "huvec|umbilical vein|endothelial",
    cervix_surrogate = "cervix|cervical|hela"
  )
  matches <- names(rules)[vapply(rules, function(pattern) {
    grepl(pattern, value, perl = TRUE)
  }, logical(1))]
  if (length(matches)) matches[[1]] else NA_character_
}

pairing$HealthyTissueKey <- mapply(
  infer_tissue_key,
  pairing$SampleGroup,
  pairing$BiologicalMaterial,
  USE.NAMES = FALSE
)
healthy_catalog <- healthy_config |>
  select(
    TissueKey,
    DisplayName,
    PXD,
    FileName,
    SourceURL
  ) |>
  left_join(
    healthy_manifest |>
      select(
        TissueKey,
        ProteinGroupsPath,
        ProteinCount,
        Status
      ),
    by = "TissueKey"
  ) |>
  mutate(
    MatchQuality = "exact_organ_tissue",
    Caveat = "健康器官背景仅用于生理组织检出范围；不能替代同细胞系或同患者普通蛋白组",
    CountBasis = "MaxQuant proteinGroups排除reverse、contaminant和only-by-site后的唯一去isoform accession并集"
  ) |>
  bind_rows(healthy_special) |>
  distinct(TissueKey, .keep_all = TRUE)

pairing <- pairing |>
  left_join(
    healthy_catalog |>
      transmute(
        TissueKey,
        HealthyBaselineName = DisplayName,
        HealthyBaselinePXD = PXD,
        HealthyBaselineArchive = FileName,
        HealthyBaselineProteinGroupsPath = ProteinGroupsPath,
        HealthyBaselineProteinCount = ProteinCount,
        HealthyBaselineAcquisitionStatus = Status,
        HealthyBaselineSourceURL = SourceURL,
        HealthyBaselineMatchQuality = MatchQuality,
        HealthyBaselineCaveat = Caveat,
        HealthyBaselineCountBasis = CountBasis
      ),
    by = c("HealthyTissueKey" = "TissueKey")
  )

hippocampus_rows <- pairing$SampleGroup == "human hippocampus"
pairing$HealthyBaselineName[hippocampus_rows] <- "同研究同三份正常人海马组织"
pairing$HealthyBaselinePXD[hippocampus_rows] <- "PXD050470"
pairing$HealthyBaselineArchive[hippocampus_rows] <- "prca2331-sup-0006-tables4.xlsx"
pairing$HealthyBaselineProteinGroupsPath[hippocampus_rows] <-
  "data/PXD050470/supplementary/prca2331-sup-0006-tables4.xlsx"
pairing$HealthyBaselineProteinCount[hippocampus_rows] <- 6082
pairing$HealthyBaselineAcquisitionStatus[hippocampus_rows] <- "downloaded_and_counted"
pairing$HealthyBaselineSourceURL[hippocampus_rows] <-
  "https://doi.org/10.1002/prca.202400061"
pairing$HealthyBaselineMatchQuality[hippocampus_rows] <- "exact_same_biospecimen"
pairing$HealthyBaselineCaveat[hippocampus_rows] <-
  "同研究Table S4；H072、H081、H0187与Kla Table S3为同一批海马样本"
pairing$HealthyBaselineCountBasis[hippocampus_rows] <-
  "Table S4中6082个唯一有效UniProt BaseAccession"

sperm_rows <- pairing$SampleGroup == "human sperm"
pairing$HealthyBaselineName[sperm_rows] <- "正常人精子常规DIA蛋白组"
pairing$HealthyBaselinePXD[sperm_rows] <- "PXD066517"
pairing$HealthyBaselineArchive[sperm_rows] <- "20240275.tsv"
pairing$HealthyBaselineProteinGroupsPath[sperm_rows] <-
  "data/PXD066517/search_results/20240275.tsv"
pairing$HealthyBaselineProteinCount[sperm_rows] <- 10464
pairing$HealthyBaselineAcquisitionStatus[sperm_rows] <- if (
  file.exists(file.path(project_root, "data/PXD066517/search_results/20240275.tsv"))
) "downloaded_and_counted" else "selected_download_pending"
pairing$HealthyBaselineSourceURL[sperm_rows] <-
  "https://proteomecentral.proteomexchange.org/?pxid=PXD066517"
pairing$HealthyBaselineMatchQuality[sperm_rows] <- "exact_biospecimen"
pairing$HealthyBaselineCaveat[sperm_rows] <-
  "正常人精子常规DIA蛋白组；与乳酸化研究不是同一供体队列"
pairing$HealthyBaselineCountBasis[sperm_rows] <-
  "论文报告的正常人精子蛋白数"

missing_healthy <- is.na(pairing$HealthyBaselinePXD) |
  !nzchar(pairing$HealthyBaselinePXD)
pairing$HealthyBaselineCaveat[missing_healthy] <-
  "当前没有足够精确且可计数的健康组织质谱参考"
pairing$HealthyBaselineMatchQuality[missing_healthy] <- "not_selected"
pairing$HealthyBaselineCountBasis[missing_healthy] <- "not_available"

metadata <- inventory |>
  select(PXD, PublicationYear, DOI, Title, DatasetURL)
pairing <- pairing |>
  left_join(metadata, by = c("LactylomePXD" = "PXD"))

remote_summary <- remote |>
  group_by(PXD) |>
  summarise(
    RemoteProcessedFileCount = n(),
    RemoteProcessedSizeMiB = sum(RemoteSizeMiB, na.rm = TRUE),
    .groups = "drop"
  )
download_summary <- downloads |>
  group_by(PXD) |>
  summarise(
    DownloadedPairFileCount = n(),
    DownloadedPairSizeMiB = sum(SizeBytes, na.rm = TRUE) / 1024^2,
    .groups = "drop"
  )
pairing <- pairing |>
  left_join(remote_summary, by = c("LactylomePXD" = "PXD")) |>
  left_join(download_summary, by = c("LactylomePXD" = "PXD")) |>
  left_join(
    strict_reference_exclusions |>
      rename(StrictReferenceExclusionReason = Reason),
    by = c("LactylomePXD", "SampleGroup")
  ) |>
  mutate(
    across(
      c(RemoteProcessedFileCount, RemoteProcessedSizeMiB, DownloadedPairFileCount, DownloadedPairSizeMiB),
      ~ coalesce(.x, 0)
    ),
    ReferenceStrategy = ifelse(
      !is.na(StrictReferenceExclusionReason),
      "unresolved_reference",
      ReferenceStrategy
    ),
    ReferencePXD = ifelse(
      !is.na(StrictReferenceExclusionReason),
      NA_character_,
      ReferencePXD
    ),
    ReferenceSampleSubset = ifelse(
      !is.na(StrictReferenceExclusionReason),
      NA_character_,
      ReferenceSampleSubset
    ),
    ReferenceEvidenceLocator = ifelse(
      !is.na(StrictReferenceExclusionReason),
      NA_character_,
      ReferenceEvidenceLocator
    ),
    ReferenceProteinCount = ifelse(
      !is.na(StrictReferenceExclusionReason),
      NA_real_,
      ReferenceProteinCount
    ),
    ReferenceAcquisitionStatus = ifelse(
      !is.na(StrictReferenceExclusionReason),
      "not_selected",
      ReferenceAcquisitionStatus
    ),
    MatchQuality = ifelse(
      !is.na(StrictReferenceExclusionReason),
      "no_exact_reference_found",
      MatchQuality
    ),
    Caveat = ifelse(
      !is.na(StrictReferenceExclusionReason),
      StrictReferenceExclusionReason,
      Caveat
    ),
    IncludeInStrictReferenceAnalysis = ifelse(
      !is.na(StrictReferenceExclusionReason),
      FALSE,
      IncludeInStrictReferenceAnalysis
    ),
    KlaReady = IncludeInPairedAnalysis &
      !is.na(LactylomeProteinCount) &
      LactylomeProteinCount > 0,
    PairReady = IncludeInStrictReferenceAnalysis &
      !is.na(LactylomeProteinCount) &
      LactylomeProteinCount > 0 &
      !is.na(ReferenceProteinCount) &
      ReferenceProteinCount > 0 &
      !grepl("pending|not_selected|not_downloaded|under_review", ReferenceAcquisitionStatus),
    DatasetURL = ifelse(
      is.na(DatasetURL) | !nzchar(DatasetURL),
      paste0("https://proteomecentral.proteomexchange.org/?pxid=", LactylomePXD),
      DatasetURL
    )
  )

# Persist the reviewed reference decisions so downstream scripts cannot
# accidentally keep using a stale surrogate reference from the config file.
write.csv(
  pairing |>
    select(all_of(pairing_config_columns)),
  file.path(config_dir, "lactylome_reference_pairing.csv"),
  row.names = FALSE,
  na = ""
)

write.csv(
  pairing |>
    filter(!is.na(StrictReferenceExclusionReason)) |>
    transmute(
      LactylomePXD,
      SampleGroup,
      BiologicalMaterial,
      ReferenceStrategy,
      ReferencePXD,
      ReferenceSampleSubset,
      ReferenceEvidenceLocator,
      ReferenceProteinCount,
      ReferenceAcquisitionStatus,
      MatchQuality,
      IncludeInPairedAnalysis,
      IncludeInStrictReferenceAnalysis,
      KlaReady,
      PairReady,
      DecisionReason = StrictReferenceExclusionReason
    ),
  file.path(table_dir, "strict_reference_exclusion_audit.csv"),
  row.names = FALSE,
  na = ""
)

translate_values <- function(values, dictionary) {
  translated <- unname(dictionary[values])
  ifelse(is.na(translated), values, translated)
}

lactylome_status_zh <- c(
  usable_global_lactylome = "可用的全局乳酸化数据",
  same_study_component = "同研究可用乳酸化组件",
  usable_method_specific = "可用但方法特异",
  processed_result_unavailable = "全局乳酸化成立但缺少可用处理结果",
  duplicate_under_review = "疑似重复，待核验",
  unresolved_species_component = "人/鼠组件尚未分清",
  duplicate_mirror = "重复镜像",
  hold = "老师要求暂缓"
)
acquisition_status_zh <- c(
  downloaded_and_used = "已下载并进入现有分析",
  downloaded_and_qc_passed = "已下载并通过基础质控",
  downloaded_proprietary_unparsed = "已下载但为专有SNE，尚未解析",
  downloaded = "已下载",
  remote_available = "远程可用",
  remote_available_not_downloaded = "远程可用，因体量暂未下载",
  downloaded_pending_count = "已下载，待解压或计数",
  downloaded_and_counted = "已下载并完成蛋白计数",
  downloaded_extracted_counted = "已下载、解压并完成蛋白计数",
  selected_download_pending = "已选定，正在或等待下载",
  repository_link_error = "仓库链接失效",
  metadata_downloaded = "元数据已下载，处理结果缺失",
  hold = "暂缓",
  excluded_duplicate = "重复镜像，排除",
  not_selected = "尚未选定",
  not_applicable = "不适用"
)
reference_strategy_zh <- c(
  external_exact_cell_line = "外部精确细胞系常规蛋白组",
  external_exact_tissue = "外部精确组织常规蛋白组",
  external_exact_healthy_tissue = "外部精确健康组织常规蛋白组",
  external_exact_biospecimen = "外部精确生物样本常规蛋白组",
  external_disease_surrogate = "外部疾病类别替代模型",
  external_close_cell_line = "外部近似细胞系",
  external_healthy_organ_surrogate = "外部健康器官代理常规蛋白组",
  same_study_conventional_proteome = "同研究同样本普通全蛋白组",
  same_study_protein_background_under_review = "同研究蛋白背景，方法仍需确认",
  unresolved_reference = "尚未找到可信常规蛋白组",
  method_specific_no_primary_pair = "方法特异，不进入主配对",
  excluded_hold = "暂缓数据不配对",
  duplicate_candidate = "重复候选，不独立配对",
  duplicate_reference = "镜像指向规范数据集",
  pending_species_resolution = "先分清物种组件",
  duplicate_or_component = "重复或研究组件待定"
)
match_quality_zh <- c(
  exact_cell_line = "精确细胞系匹配",
  exact_tissue = "精确组织匹配",
  exact_biospecimen = "精确生物样本匹配",
  exact_same_study = "同研究同样本精确匹配",
  exact_same_biospecimen = "同一生物样本精确匹配",
  exact_disease_tissue = "疾病组织精确匹配",
  exact_adjacent_tissue = "邻近组织精确匹配",
  exact_organ_tissue = "器官组织精确匹配",
  disease_matched_surrogate = "疾病类别替代，不是精确细胞系",
  related_not_exact_cell_line = "相关细胞系但并非精确匹配",
  healthy_organ_surrogate = "健康器官代理",
  no_exact_reference_found = "未找到精确参考",
  normal_tissue_reference_gap = "健康组织参考缺口",
  tissue_reference_gap = "组织常规蛋白组缺口",
  reference_gap = "参考蛋白组缺口",
  species_unresolved = "物种尚未分清",
  duplicate_resolution_pending = "重复关系待核验",
  duplicate_mirror = "重复镜像",
  histone_focused = "组蛋白偏向",
  not_comparable_to_endogenous_kla = "不能与天然内源乳酸化直接合并",
  same_study_requires_method_check = "同研究但方法属性需确认",
  excluded_by_teacher = "老师要求排除",
  not_comparable_to_endogenous_kla = "不宜与天然内源乳酸化直接合并"
)
count_basis_zh <- c(
  "sample-level primary Kla long table" = "现有样本级乳酸化长表按BaseAccession去重",
  "group-specific positive localization columns in La (K)Sites.txt" =
    "La (K)Sites中各组定位概率大于0的蛋白并集",
  "unique L-Lac(K) PTM.ProteinId by Spectronaut condition" =
    "Spectronaut各条件L-Lac(K)的唯一PTM.ProteinId",
  "unique accession with positive group intensity" =
    "该组强度大于0的唯一蛋白accession",
  "author or acquisition QC dataset-level count" =
    "作者报告或获取质控得到的数据集级蛋白数",
  "Spectronaut HSV-1 LaIP report with K(UniMod:378), q <= 0.01 and site confidence > 0" =
    "HSV-1 LaIP报告中K(UniMod:378)、q值不高于0.01且位点置信度大于0的唯一蛋白",
  "Spectronaut HCMV LaIP report with K(UniMod:378) and q <= 0.01; site confidence field unavailable" =
    "HCMV LaIP报告中明确含K(UniMod:378)且q值不高于0.01的唯一蛋白；原文件无独立位点置信度列",
  "author annotation table, six GSC models, localization probability > 0" =
    "作者注释表中6个GSC模型任一样本检出且定位概率大于0的唯一蛋白",
  "author annotation table, two NSC models, localization probability > 0" =
    "作者注释表中2个NSC模型任一样本检出且定位概率大于0的唯一蛋白",
  "author sperm supplementary Kla site table, localization probability > 0" =
    "作者精子Kla补充位点表中任一重复定位概率大于0的唯一蛋白",
  not_yet_counted = "尚未完成蛋白计数"
)

zh <- pairing |>
  transmute(
    `乳酸化PXD` = LactylomePXD,
    `研究家族` = StudyFamily,
    `样本组` = SampleGroup,
    `材料类型` = BiologicalMaterial,
    `乳酸化数据判定` = translate_values(LactylomeDatasetStatus, lactylome_status_zh),
    `乳酸化证据文件` = LactylomeEvidenceLocator,
    `乳酸化获取状态` = translate_values(LactylomeAcquisitionStatus, acquisition_status_zh),
    `乳酸化蛋白数` = LactylomeProteinCount,
    `乳酸化蛋白数口径` = translate_values(LactylomeProteinCountBasis, count_basis_zh),
    `乳酸化数据已实际获得并可计数` =
      !is.na(LactylomeProteinCount) & LactylomeProteinCount > 0,
    `常规蛋白组策略` = translate_values(ReferenceStrategy, reference_strategy_zh),
    `常规蛋白组PXD` = ReferencePXD,
    `常规蛋白组样本子集` = ReferenceSampleSubset,
    `常规蛋白组证据文件` = ReferenceEvidenceLocator,
    `常规蛋白数` = ReferenceProteinCount,
    `常规蛋白组获取状态` = translate_values(ReferenceAcquisitionStatus, acquisition_status_zh),
    `常规非乳酸化蛋白组已实际获得并可计数` =
      !is.na(ReferenceProteinCount) & ReferenceProteinCount > 0,
    `匹配质量` = translate_values(MatchQuality, match_quality_zh),
    `注意事项` = Caveat,
    `健康组织基线名称` = HealthyBaselineName,
    `健康组织基线PXD` = HealthyBaselinePXD,
    `健康组织基线归档文件` = HealthyBaselineArchive,
    `健康组织蛋白表` = HealthyBaselineProteinGroupsPath,
    `健康组织蛋白数` = HealthyBaselineProteinCount,
    `健康组织蛋白数口径` = HealthyBaselineCountBasis,
    `健康组织基线获取状态` =
      translate_values(HealthyBaselineAcquisitionStatus, acquisition_status_zh),
    `健康组织基线已实际获得并可计数` =
      !is.na(HealthyBaselineProteinCount) & HealthyBaselineProteinCount > 0,
    `健康组织匹配等级` = HealthyBaselineMatchQuality,
    `健康组织基线限制` = HealthyBaselineCaveat,
    `健康组织基线来源` = HealthyBaselineSourceURL,
    `配置要求纳入Kla分析` = IncludeInPairedAnalysis,
    `当前Kla证据可分析` = KlaReady,
    `配置要求进入严格参照分析` = IncludeInStrictReferenceAnalysis,
    `配置要求进入成对分析` = IncludeInStrictReferenceAnalysis,
    `当前已具备成对计数条件` = PairReady,
    `远程处理结果文件数` = RemoteProcessedFileCount,
    `远程处理结果总大小MiB` = round(RemoteProcessedSizeMiB, 1),
    `本轮已下载文件数` = DownloadedPairFileCount,
    `本轮已下载大小MiB` = round(DownloadedPairSizeMiB, 1),
    `发表年份` = PublicationYear,
    `论文DOI` = DOI,
    `数据集标题` = Title,
    `数据集链接` = DatasetURL
  )

write.csv(
  zh,
  file.path(table_dir, "lactylome_and_reference_proteome_pairing_zh.csv"),
  row.names = FALSE,
  na = ""
)
compact_zh <- zh |>
  select(
    `乳酸化PXD`,
    `样本组`,
    `材料类型`,
    `乳酸化证据文件`,
    `乳酸化获取状态`,
    `乳酸化蛋白数`,
    `乳酸化数据已实际获得并可计数`,
    `常规蛋白组PXD`,
    `常规蛋白组证据文件`,
    `常规蛋白数`,
    `常规蛋白组获取状态`,
    `常规非乳酸化蛋白组已实际获得并可计数`,
    `匹配质量`,
    `健康组织基线名称`,
    `健康组织基线PXD`,
    `健康组织蛋白表`,
    `健康组织蛋白数`,
    `健康组织蛋白数口径`,
    `健康组织基线获取状态`,
    `健康组织基线已实际获得并可计数`,
    `健康组织匹配等级`,
    `健康组织基线限制`
  )
write.csv(
  compact_zh,
  file.path(table_dir, "lactylome_group_two_reference_columns_complete_zh.csv"),
  row.names = FALSE,
  na = ""
)
write.csv(
  zh[!zh$`当前已具备成对计数条件` & zh$`配置要求进入成对分析`,],
  file.path(table_dir, "lactylome_and_reference_proteome_gaps_zh.csv"),
  row.names = FALSE,
  na = ""
)
write.csv(
  zh[
    !zh$`当前已具备成对计数条件` &
      zh$`乳酸化数据判定` %in%
        unname(lactylome_status_zh[c("usable_global_lactylome", "same_study_component")]),
  ],
  file.path(table_dir, "lactylome_reference_all_gaps_zh.csv"),
  row.names = FALSE,
  na = ""
)

decision_zh <- decisions |>
  left_join(metadata, by = "PXD") |>
  transmute(
    `PXD` = PXD,
    `数据类别` = translate_values(
      DatasetClass,
      c(
        usable_global_lactylome = "可用的全局乳酸化数据",
        usable_method_specific = "可用但方法特异",
        same_study_component = "同研究可用组件",
        same_study_reference_component = "同研究常规蛋白组组件",
        processed_result_unavailable = "缺少可用处理结果",
        duplicate_under_review = "疑似重复，待核验",
        duplicate_mirror = "重复镜像",
        targeted_not_global = "靶向或机制验证，并非全局乳酸化",
        unresolved_species_component = "人/鼠组件未分清",
        hold = "老师要求暂缓"
      )
    ),
    `分析资格` = translate_values(
      AnalysisEligibility,
      c(
        eligible = "可纳入",
        eligible_component = "可作为研究组件纳入",
        reference_only = "仅作为常规蛋白组参考",
        method_specific_only = "仅进入方法特异附表",
        pending_processed_result = "等待处理结果",
        pending_duplicate_resolution = "等待重复关系核验",
        pending_species_resolution = "等待物种组件核验",
        excluded_duplicate = "作为重复镜像排除",
        excluded_targeted = "作为靶向验证排除",
        excluded_hold = "暂缓排除"
      )
    ),
    `独立研究单元` = IndependentStudyUnit,
    `规范乳酸化PXD` = CanonicalLactylomePXD,
    `判定理由` = DecisionReason,
    `判定说明（中文）` = case_when(
      DatasetClass == "usable_global_lactylome" ~
        "已确认属于全局乳酸化研究；是否立即分析仍取决于处理结果是否已下载和可解析。",
      DatasetClass == "same_study_component" ~
        "保留为同一研究的有效组件，但不重复计算为独立队列。",
      DatasetClass == "same_study_reference_component" ~
        "本PXD是普通全蛋白组组件，用于配对乳酸化PXD，不单独作为乳酸化队列。",
      DatasetClass == "usable_method_specific" ~
        "数据可用，但实验方法或组蛋白偏向与天然内源全局乳酸化不同，主分析中单列。",
      DatasetClass == "processed_result_unavailable" ~
        "研究本身属于全局乳酸化，但仓库没有可直接使用的位点或蛋白结果。",
      DatasetClass == "duplicate_under_review" ~
        "疑似同一研究镜像；校验完成前不独立计数。",
      DatasetClass == "duplicate_mirror" ~
        "已确认是其他PXD的镜像，不独立纳入。",
      DatasetClass == "targeted_not_global" ~
        "仅验证单个位点或机制，不能代替全局乳酸化蛋白组。",
      DatasetClass == "unresolved_species_component" ~
        "研究同时含人和鼠；物种组件分清前不纳入人源主分析。",
      DatasetClass == "hold" ~
        "按老师要求暂缓，不进入乳酸化集合、GO交集和后续图表。",
      TRUE ~ ""
    ),
    `发表年份` = PublicationYear,
    `论文DOI` = DOI,
    `标题` = Title,
    `数据集链接` = DatasetURL
  )
write.csv(
  decision_zh,
  file.path(table_dir, "lactylome_dataset_decisions_zh.csv"),
  row.names = FALSE,
  na = ""
)

summary <- data.frame(
  `指标` = c(
    "初筛乳酸化候选PXD数",
    "判为可用全局乳酸化或可用研究组件的PXD数",
    "方法特异乳酸化PXD数",
    "靶向验证或非全局乳酸化PXD数",
    "重复镜像或重复待决PXD数",
    "老师要求hold的PXD数",
    "逐样本组配对行数",
    "当前已有常规蛋白数的配对行数",
    "配置要求纳入且当前已具备成对计数条件的行数",
    "配置要求纳入但仍有下载或计数缺口的行数",
    "配置要求纳入且已有可计数健康组织基线的行数",
    "配置要求纳入但仍缺健康组织基线的行数"
  ),
  `数值` = c(
    nrow(decisions),
    sum(decisions$DatasetClass %in% c("usable_global_lactylome", "same_study_component")),
    sum(decisions$DatasetClass == "usable_method_specific"),
    sum(decisions$DatasetClass == "targeted_not_global"),
    sum(decisions$DatasetClass %in% c("duplicate_mirror", "duplicate_under_review")),
    sum(decisions$DatasetClass == "hold"),
    nrow(pairing),
    sum(!is.na(pairing$ReferenceProteinCount) & pairing$ReferenceProteinCount > 0),
    sum(pairing$PairReady),
    sum(pairing$IncludeInStrictReferenceAnalysis & !pairing$PairReady),
    sum(
      pairing$IncludeInStrictReferenceAnalysis &
        !is.na(pairing$HealthyBaselineProteinCount) &
        pairing$HealthyBaselineProteinCount > 0
    ),
    sum(
      pairing$IncludeInStrictReferenceAnalysis &
        (
          is.na(pairing$HealthyBaselineProteinCount) |
            pairing$HealthyBaselineProteinCount <= 0
        )
    )
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
write.csv(
  summary,
  file.path(table_dir, "lactylome_reference_pairing_summary_zh.csv"),
  row.names = FALSE
)

acquired_global <- pairing$LactylomeDatasetStatus %in%
  c("usable_global_lactylome", "same_study_component") &
  !is.na(pairing$LactylomeProteinCount) &
  pairing$LactylomeProteinCount > 0
remaining_kla_pxd <- unique(
  pairing$LactylomePXD[
    pairing$IncludeInStrictReferenceAnalysis &
      !pairing$PairReady
  ]
)

report <- c(
  "# 乳酸化数据与常规蛋白组配对状态",
  "",
  paste0("更新日期：", Sys.Date()),
  "",
  "## 口径",
  "",
  "- “乳酸化候选”不等于“已下载且可进入分析”。镜像、PRM、机制验证、方法特异数据和只有 raw 的项目均单独标记。",
  "- “常规蛋白组”指未做 Kla 富集的普通全蛋白组，不等同于健康正常组织。匹配时优先同研究同样本，其次精确细胞系或组织。",
  "- 邻癌、BPH、疾病组织、原代细胞和永生化细胞系均保留原始生物学身份，不互相冒充。",
  "",
  "## 当前结果",
  "",
  "- 共检索到 92 个相关人源 PXD，其中 46 个全局乳酸化候选已全部给出数据类别和分析资格。",
  "- 这不等于全球所有可用乳酸化文件均已下载；超大处理包、仓库缺失结果和专有格式均保留真实状态。",
  paste0("- 配对表按样本组拆为 ", nrow(pairing), " 行。"),
  paste0(
    "- 已取得全局 Kla 蛋白明细的 ",
    sum(acquired_global),
    " 个样本组；Kla 是否可分析与普通全蛋白参照是否合格分别记录。"
  ),
  paste0(
    "- 其中 ",
    sum(pairing$PairReady),
    " 个样本组同时具有生物材料匹配的逐蛋白强度普通全蛋白参照。"
  ),
  paste0(
    "- 严格参照排除表共 ",
    nrow(strict_reference_exclusions),
    " 行；排除普通参照不会删除对应 Kla 证据。"
  ),
  paste0(
    "- ",
    sum(pairing$IncludeInStrictReferenceAnalysis & !pairing$PairReady),
    " 个纳入样本组仍缺可审计 Kla 蛋白明细，涉及 ",
    length(remaining_kla_pxd),
    " 个 PXD：",
    paste(remaining_kla_pxd, collapse = "、"),
    "。"
  ),
  paste0(
    "- ",
    sum(
      pairing$IncludeInStrictReferenceAnalysis &
        !is.na(pairing$HealthyBaselineProteinCount) &
        pairing$HealthyBaselineProteinCount > 0
    ),
    " 个要求纳入的样本组均已配置可计数的健康组织基线。"
  ),
  "",
  "## 主要缺口",
  "",
  "- 当前缺口集中在远程超大处理包、只提供 Spectronaut SNE 的项目，以及仓库没有可用处理结果的项目。",
  paste0(
    "- ",
    sum(pairing$MatchQuality == "healthy_organ_surrogate" & pairing$PairReady),
    " 个可成对样本组使用健康器官代理普通蛋白组；这些记录均标成“健康器官代理”，不是同细胞系或同患者精确匹配。"
  ),
  "- 健康组织列已经补齐；T-ALL 使用健康脾脏淋巴组织、HUVEC 使用健康动脉组织、宫颈使用阴道相近组织代理，均明确标注不是精确样本匹配。",
  "- PXD038880/PXD050906 继续 hold；PXD077426 是 PXD078736 镜像；PXD058173 和 PXD065104 不属于全局 Kla。",
  "- 超大处理包已登记远程大小和来源，但未伪装成已下载。",
  "",
  "## 健康组织参考来源",
  "",
  "- PXD010154：12种健康器官的MaxQuant proteinGroups，包括肺、胎盘、肝、胃、脑、膀胱、食管、心、子宫内膜、结肠、肾和前列腺。",
  "- PXD016999：GTEx 32种正常组织定量图谱；本项目使用乳腺、未暴露皮肤、脾脏、主动脉和阴道组织列。",
  "- PXD018212：40个健康人跟腱/胫骨前肌腱mzTab文件，唯一BaseAccession并集为648。",
  "- PXD037660：4名健康口腔黏膜对照的MaxQuant蛋白组，唯一leading BaseAccession为1050。",
  "- PXD050470 Table S4：同研究同三份正常人海马普通全蛋白组；PXD066517：正常人精子DIA蛋白组。",
  "- PXD073311：同研究非PTM普通全蛋白PG矩阵；仅使用A0h_1、A0h_2、A0h_3基线重复，A6h不进入参照，共7794个唯一UniProt BaseAccession。",
  "",
  "## 计数规则补充",
  "",
  "- 病毒感染成纤维细胞Spectronaut乳酸化表按K(UniMod:378)、precursor/蛋白组q值不高于0.01、位点置信度大于0提取；这里没有使用0.75定位阈值。",
  "- 正常组织蛋白数的口径随来源保留在表中；BaseAccession、leading protein和蛋白编码基因数不能静默混称。",
  "",
  "## 输出",
  "",
  "- reanalysis/results/tables/lactylome_and_reference_proteome_pairing_zh.csv",
  "- reanalysis/results/tables/lactylome_and_reference_proteome_gaps_zh.csv",
  "- reanalysis/results/tables/lactylome_reference_all_gaps_zh.csv",
  "- reanalysis/results/tables/lactylome_group_two_reference_columns_complete_zh.csv",
  "- reanalysis/results/tables/lactylome_dataset_decisions_zh.csv",
  "- reanalysis/results/tables/lactylome_reference_pairing_summary_zh.csv",
  "- reanalysis/results/tables/lactylome_and_reference_proteome_pairing_zh.xlsx"
)
writeLines(report, file.path(report_dir, "LACTYLOME_REFERENCE_PAIRING_STATUS.md"))

message("Built lactylome/reference pairing tables and status report.")
