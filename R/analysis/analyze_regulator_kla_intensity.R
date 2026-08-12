#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readxl)
  library(stringr)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) normalizePath(args[[1]]) else normalizePath(".")
table_dir <- file.path(project_root, "results", "tables")
figure_dir <- file.path(project_root, "results", "figures")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

regulator_path <- file.path(
  project_root, "data", "identifier",
  "乳酸化调控因子_Writer-Eraser-Reader.xlsx"
)
detection_path <- file.path(table_dir, "kla_regulator_cell_type_long.csv")
mapping_path <- file.path(
  project_root, "config",
  "lactylation_regulator_uniprot_mapping.csv"
)

required_files <- c(regulator_path, detection_path, mapping_path)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) {
  stop("Missing required files: ", paste(missing_files, collapse = ", "))
}

regulators <- read_excel(regulator_path) |>
  transmute(
    Role = trimws(as.character(分组)),
    GeneSymbol = toupper(trimws(as.character(基因))),
    RoleEntryOrder = row_number()
  ) |>
  filter(
    Role %in% c("Writer", "Eraser", "Writer-Eraser", "Reader"),
    !is.na(GeneSymbol),
    nzchar(GeneSymbol)
  ) |>
  distinct(Role, GeneSymbol, .keep_all = TRUE)
target_genes <- unique(regulators$GeneSymbol)
plot_font <- "Arial Unicode MS"

accession_map <- read.csv(
  mapping_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
) |>
  transmute(
    BaseAccession = sub("-[0-9]+$", "", BaseAccession),
    GeneSymbol = toupper(GeneSymbol)
  ) |>
  filter(GeneSymbol %in% target_genes)
accession_to_gene <- setNames(accession_map$GeneSymbol, accession_map$BaseAccession)
regulators <- regulators |>
  left_join(accession_map, by = "GeneSymbol") |>
  mutate(
    RegulatorDisplayName = ifelse(
      BaseAccession == "Q92830",
      "GCN5 (KAT2A)",
      GeneSymbol
    )
  )
if (any(is.na(regulators$BaseAccession))) {
  stop(
    "Missing UniProt mapping for regulator genes: ",
    paste(regulators$GeneSymbol[is.na(regulators$BaseAccession)], collapse = ", ")
  )
}
target_accessions <- unique(regulators$BaseAccession)

detection <- read.csv(
  detection_path,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  na.strings = c("", "NA")
)
sample_catalog <- detection |>
  distinct(
    PXD, SampleGroup, SampleGroupID, RowLabel,
    BiologicalMaterial, LactylomeProteinCount,
    GeneLevelAuditStatus, Parser, AuditReason
  ) |>
  mutate(RowOrder = row_number())
if (nrow(sample_catalog) != 40) {
  stop("Expected 40 sample groups, found ", nrow(sample_catalog))
}

whole_heatmap_rows_path <- file.path(
  table_dir,
  "kla_regulator_whole_proteome_heatmap_rows.csv"
)
four_class_path <- file.path(
  project_root,
  "config",
  "four_class_sample_grouping.csv"
)
category_order <- c(
  "normal_tissue",
  "cancer_tissue",
  "normal_cells",
  "cancer_cells"
)
category_labels_zh <- c(
  normal_tissue = "正常/非肿瘤组织",
  cancer_tissue = "癌症组织",
  normal_cells = "正常/非肿瘤细胞",
  cancer_cells = "癌症细胞"
)
category_labels_en <- c(
  normal_tissue = "Normal/non-tumor tissues",
  cancer_tissue = "Cancer tissues",
  normal_cells = "Normal/non-tumor cells",
  cancer_cells = "Cancer cells"
)
missing_axis_files <- c(whole_heatmap_rows_path, four_class_path)[!file.exists(
  c(whole_heatmap_rows_path, four_class_path)
)]
if (length(missing_axis_files)) {
  stop(
    "Missing whole-proteome heatmap axis inputs; run the whole-proteome " ,
    "analysis first: ",
    paste(missing_axis_files, collapse = ", ")
  )
}
whole_heatmap_rows <- read.csv(
  whole_heatmap_rows_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
if (nrow(whole_heatmap_rows) != 30) {
  stop(
    "Expected 30 unique whole-proteome heatmap rows, found ",
    nrow(whole_heatmap_rows)
  )
}
whole_axis_members <- bind_rows(lapply(seq_len(nrow(whole_heatmap_rows)), function(i) {
  data.frame(
    PXD = strsplit(
      whole_heatmap_rows$LinkedKlaPXD[[i]],
      ";",
      fixed = TRUE
    )[[1]],
    SampleGroup = strsplit(
      whole_heatmap_rows$LinkedKlaSampleGroup[[i]],
      ";",
      fixed = TRUE
    )[[1]],
    WholeProteomeDisplayRowOrder = whole_heatmap_rows$HeatmapDisplayRowOrder[[i]],
    stringsAsFactors = FALSE
  )
}))
four_class <- read.csv(
  four_class_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
) |>
  transmute(
    PXD,
    SampleGroup,
    Category,
    CategoryZh = unname(category_labels_zh[Category]),
    CategoryEn = unname(category_labels_en[Category]),
    CategoryOrder = match(Category, category_order)
  )
category_max_order <- whole_heatmap_rows |>
  group_by(Category) |>
  summarise(
    MaxWholeProteomeDisplayRowOrder = max(HeatmapDisplayRowOrder),
    .groups = "drop"
  )
sample_catalog <- sample_catalog |>
  left_join(four_class, by = c("PXD", "SampleGroup")) |>
  left_join(
    whole_axis_members,
    by = c("PXD", "SampleGroup")
  ) |>
  left_join(category_max_order, by = "Category") |>
  group_by(Category) |>
  arrange(RowOrder, .by_group = TRUE) |>
  mutate(
    UnmatchedWithinCategory = cumsum(
      is.na(WholeProteomeDisplayRowOrder)
    ),
    ComparisonRowOrder = ifelse(
      !is.na(WholeProteomeDisplayRowOrder),
      WholeProteomeDisplayRowOrder,
      MaxWholeProteomeDisplayRowOrder +
        UnmatchedWithinCategory / 100
    )
  ) |>
  ungroup() |>
  arrange(CategoryOrder, ComparisonRowOrder, RowOrder) |>
  mutate(
    RowLabelZh = RowLabel,
    RowLabelEn = paste0(SampleGroup, " · Kla:", PXD)
  )

base_accession <- function(values) {
  values <- sub("^.*\\|([^|]+)\\|.*$", "\\1", as.character(values))
  sub("-[0-9]+$", "", values)
}

match_target_accession <- function(values) {
  vapply(as.character(values), function(value) {
    if (is.na(value) || !nzchar(value)) return(NA_character_)
    tokens <- unique(trimws(unlist(strsplit(value, "[;, ]+"))))
    tokens <- base_accession(tokens)
    hits <- tokens[tokens %in% target_accessions]
    if (length(hits)) hits[[1]] else NA_character_
  }, character(1))
}

accession_feature <- function(values) {
  vapply(as.character(values), function(value) {
    if (is.na(value) || !nzchar(value)) return(NA_character_)
    tokens <- unique(trimws(unlist(strsplit(value, "[;, ]+"))))
    tokens <- base_accession(tokens)
    tokens <- tokens[nzchar(tokens)]
    if (length(tokens)) paste0("ACC:", tokens[[1]]) else NA_character_
  }, character(1))
}

safe_numeric <- function(values) {
  suppressWarnings(as.numeric(gsub(",", "", as.character(values), fixed = TRUE)))
}

relative_path <- function(path) {
  sub(paste0("^", project_root, "/"), "", normalizePath(path))
}

quant_parts <- list()

add_quant <- function(
  data,
  pxd,
  sample_group,
  feature_values,
  target_values,
  sample_columns,
  sample_labels = sample_columns,
  measurement_level,
  quant_field,
  source_file
) {
  sample_columns <- intersect(sample_columns, names(data))
  if (!length(sample_columns) || !nrow(data)) return(invisible(NULL))
  if (length(sample_labels) != length(sample_columns)) {
    stop("sample_labels and sample_columns length mismatch for ", source_file)
  }
  base <- data.frame(
    FeatureID = as.character(feature_values),
    TargetAccession = as.character(target_values),
    stringsAsFactors = FALSE
  ) |>
    mutate(
      TargetGene = unname(accession_to_gene[TargetAccession])
    )
  values <- as.data.frame(lapply(data[sample_columns], safe_numeric))
  names(values) <- sample_labels
  long <- bind_cols(base, values) |>
    pivot_longer(
      cols = all_of(sample_labels),
      names_to = "QuantSample",
      values_to = "Signal"
    ) |>
    filter(
      !is.na(FeatureID),
      nzchar(FeatureID),
      is.finite(Signal),
      Signal > 0
    ) |>
    mutate(
      PXD = pxd,
      SampleGroup = sample_group,
      MeasurementLevel = measurement_level,
      QuantField = quant_field,
      SourceFile = relative_path(source_file),
      CanonicalFeature = ifelse(
        !is.na(TargetAccession) & nzchar(TargetAccession),
        paste0("ACC:", TargetAccession),
        FeatureID
      )
    ) |>
    group_by(
      PXD, SampleGroup, QuantSample, CanonicalFeature,
      TargetAccession, TargetGene,
      MeasurementLevel, QuantField, SourceFile
    ) |>
    summarise(Signal = sum(Signal, na.rm = TRUE), .groups = "drop")
  if (nrow(long)) quant_parts[[length(quant_parts) + 1]] <<- long
  invisible(NULL)
}

read_delimited <- function(path) {
  if (grepl("\\.csv$", path, ignore.case = TRUE)) {
    read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  } else {
    read.delim(
      path,
      check.names = FALSE,
      stringsAsFactors = FALSE,
      quote = "",
      comment.char = ""
    )
  }
}

extract_maxquant <- function(
  path,
  pxd,
  sample_group,
  sample_tokens = NULL,
  sample_labels = sample_tokens,
  sheet = NULL
) {
  if (!file.exists(path)) return(invisible(NULL))
  data <- if (grepl("\\.xlsx$", path, ignore.case = TRUE)) {
    read_excel(path, sheet = sheet)
  } else {
    read_delimited(path)
  }
  if (!"Gene names" %in% names(data)) return(invisible(NULL))
  keep <- rep(TRUE, nrow(data))
  if ("Reverse" %in% names(data)) {
    keep <- keep & (is.na(data$Reverse) | data$Reverse != "+")
  }
  if ("Potential contaminant" %in% names(data)) {
    keep <- keep &
      (is.na(data$`Potential contaminant`) | data$`Potential contaminant` != "+")
  }
  if ("id" %in% names(data)) {
    keep <- keep & !is.na(data$id)
  }
  data <- data[keep, , drop = FALSE]
  if (!nrow(data)) return(invisible(NULL))

  if (is.null(sample_tokens)) {
    sample_columns <- "Intensity"
    sample_labels <- basename(dirname(dirname(path)))
    localization_columns <- "Localization prob"
  } else {
    sample_columns <- paste("Intensity", sample_tokens)
    localization_columns <- paste("Localization prob", sample_tokens)
  }
  valid <- sample_columns %in% names(data)
  sample_columns <- sample_columns[valid]
  sample_labels <- sample_labels[valid]
  localization_columns <- localization_columns[valid]
  for (index in seq_along(sample_columns)) {
    local_keep <- rep(TRUE, nrow(data))
    if (localization_columns[[index]] %in% names(data)) {
      local_keep <- safe_numeric(data[[localization_columns[[index]]]]) > 0
      local_keep[is.na(local_keep)] <- FALSE
    }
    subset <- data[local_keep, , drop = FALSE]
    add_quant(
      subset, pxd, sample_group,
      feature_values = accession_feature(subset$Proteins),
      target_values = match_target_accession(subset$Proteins),
      sample_columns = sample_columns[[index]],
      sample_labels = sample_labels[[index]],
      measurement_level = "Kla_site",
      quant_field = "MaxQuant site Intensity",
      source_file = path
    )
  }
}

# MaxQuant site tables.
pxd014_files <- list.files(
  file.path(project_root, "data", "PXD014870", "search_results"),
  pattern = "^Lactyl \\(K\\)Sites\\.txt$",
  recursive = TRUE,
  full.names = TRUE
)
for (path in pxd014_files) {
  extract_maxquant(path, "PXD014870", "MCF7")
}

extract_maxquant(
  file.path(
    project_root,
    "data/PXD033146/search_results/extracted_pairing/search_result-HA119TPLa/La (K)Sites.txt"
  ),
  "PXD033146", "pathological rotator cuff tendon"
)
extract_maxquant(
  file.path(project_root, "data/PXD036307/search_results/extracted/txt/La (K)Sites.txt"),
  "PXD036307", "normal human lung",
  c("PTB340", "PTB342", "PTB344", "PTB346", "PTB364", "PTB372")
)
extract_maxquant(
  file.path(project_root, "data/PXD050147/search_results/Lactyl_K_Sites.txt"),
  "PXD050147", "HepG2 WT and SIRT1 or SIRT3 KO",
  c(
    "SIRT1KO_Kla_rep1", "SIRT1KO_Kla_rep2", "SIRT1KO_Kla_rep3",
    "SIRT3KO_Kla_rep1", "SIRT3KO_Kla_rep2", "SIRT3KO_Kla_rep3",
    "WT_Kla_rep1", "WT_Kla_rep2", "WT_Kla_rep3"
  )
)
extract_maxquant(
  file.path(project_root, "data/PXD058534/search_results/extracted_pairing/txt/La (K)Sites.txt"),
  "PXD058534", "pretreated HK-2", "HK2"
)
for (spec in list(
  list("A", "MCF7"),
  list("B", "MDA-MB-468"),
  list("C", "MCF10A"),
  list("D", "T-47D")
)) {
  extract_maxquant(
    file.path(project_root, "data/PXD060185/search_results/RESULT/combined/txt/La (K)Sites.txt"),
    "PXD060185", spec[[2]], spec[[1]], spec[[2]]
  )
}
extract_maxquant(
  file.path(project_root, "data/PXD062720/search_results/extracted_pairing/txt/La (K)Sites.txt"),
  "PXD062720", "bladder cancer cells treated with EPI",
  c("A_1", "A_3")
)
extract_maxquant(
  file.path(project_root, "data/PXD063047/search_results/extracted/combined/txt/La (K)Sites.txt"),
  "PXD063047", "normal pregnancy placenta",
  c("Con_1", "Con_2", "Con_3")
)
extract_maxquant(
  file.path(project_root, "data/PXD063047/search_results/extracted/combined/txt/La (K)Sites.txt"),
  "PXD063047", "severe preeclampsia placenta",
  c("PE_1", "PE_2", "PE_3")
)
extract_maxquant(
  file.path(project_root, "data/PXD063266/search_results/LactylSites.xlsx"),
  "PXD063266", "PC-3M", c("1", "2", "3"),
  sheet = "Lactyl (K)Sites"
)
extract_maxquant(
  file.path(project_root, "data/PXD064038/search_results/extracted_pairing/txt/txt/La (K)Sites.txt"),
  "PXD064038", "MEC and NEC ESCC groups",
  c("MEC_1", "MEC_2", "MEC_3", "NEC_1", "NEC_2", "NEC_3")
)
extract_maxquant(
  file.path(project_root, "data/PXD078736/search_results/txt/La(K)Sites.txt"),
  "PXD078736", "HK-2 control and mannitol",
  c("ctr_1", "ctr_2", "ctr_3", "man_1", "man_2", "man_3")
)
extract_maxquant(
  file.path(
    project_root,
    "data/PXD028737/search_results/extracted_pairing/combined/combined/txt/La(K)Sites.txt"
  ),
  "PXD028737", "HMC3",
  c("H0", "H24"),
  c("HMC3_normoxia", "HMC3_hypoxia")
)

# PEAKS PXD028488: sum lactylated peptide areas to protein accession.
peaks_specs <- list(
  list(
    "HEK293T",
    "data/PXD028488/search_results/Enrichment-Search files/HEK293T-Enrichment-all HCD-Search files/protein-peptides.csv"
  ),
  list(
    "HCT116",
    "data/PXD028488/search_results/Enrichment-Search files/HCT116-Enrichment-Search files/protein-peptides.csv"
  ),
  list(
    "TALL-104",
    "data/PXD028488/search_results/Enrichment-Search files/TALL-NALAC-Search files/protein-peptides.csv"
  )
)
for (spec in peaks_specs) {
  path <- file.path(project_root, spec[[2]])
  data <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  keep <- grepl("Lactyl|\\+72\\.02", data$PTM, ignore.case = TRUE) |
    grepl("K\\(\\+72\\.02\\)|K\\(Lactyl\\)", data$Peptide, ignore.case = TRUE)
  data <- data[keep, , drop = FALSE]
  area_columns <- grep("^Area", names(data), value = TRUE)
  add_quant(
    data, "PXD028488", spec[[1]],
    feature_values = accession_feature(data$`Protein Accession`),
    target_values = match_target_accession(data$`Protein Accession`),
    sample_columns = area_columns,
    sample_labels = paste(spec[[1]], sub("^Area ?", "", area_columns), sep = "__"),
    measurement_level = "modified_peptide",
    quant_field = "PEAKS peptide Area",
    source_file = path
  )
}

# Proteome Discoverer normalized abundance from explicitly lactylated peptides.
path <- file.path(
  project_root,
  "data/PXD046800/search_results/HFX2_LFQ_QB001_Lacty_PeptideGroups.txt"
)
data <- read.delim(
  path, check.names = FALSE, stringsAsFactors = FALSE,
  quote = "\"", comment.char = ""
)
data <- data[
  grepl("Lacty|Lactyl|La \\(K\\)", data$Modifications, ignore.case = TRUE),
  ,
  drop = FALSE
]
for (spec in list(
  list("hypertrophic scar", "HSP"),
  list("adjacent skin", "NSP")
)) {
  columns <- grep(
    paste0('^"Abundances \\(Normalized\\):.*', spec[[2]]),
    names(data), value = TRUE
  )
  if (!length(columns)) {
    columns <- grep(
      paste0("^Abundances \\(Normalized\\):.*", spec[[2]]),
      names(data), value = TRUE
    )
  }
  add_quant(
    data, "PXD046800", spec[[1]],
    feature_values = accession_feature(data$`Master Protein Accessions`),
    target_values = match_target_accession(data$`Master Protein Accessions`),
    sample_columns = columns,
    sample_labels = sub("^.*Sample, |\"$", "", columns),
    measurement_level = "modified_peptide",
    quant_field = "Proteome Discoverer lactylated peptide normalized abundance",
    source_file = path
  )
}
rm(data)

# PXD053474: use enriched DDA peptide areas and enriched DIA PTM quantities.
dda_files <- list.files(
  file.path(
    project_root, "data", "PXD053474", "search_results", "extracted",
    "Subcellular-Whole-cell-lysates-Enriched-DDA-Search-files"
  ),
  pattern = "^protein-peptides\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)
for (path in dda_files) {
  data <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  keep <- grepl("Lactyl|\\+72\\.02", data$PTM, ignore.case = TRUE) |
    grepl("K\\(\\+72\\.02\\)|K\\(Lactyl\\)", data$Peptide, ignore.case = TRUE)
  data <- data[keep, , drop = FALSE]
  columns <- grep("^Area", names(data), value = TRUE)
  labels <- paste("DDA", basename(dirname(path)), sub("^Area ?", "", columns), sep = "__")
  add_quant(
    data, "PXD053474", "HCT116",
    feature_values = accession_feature(data$`Protein Accession`),
    target_values = match_target_accession(data$`Protein Accession`),
    sample_columns = columns,
    sample_labels = labels,
    measurement_level = "modified_peptide",
    quant_field = "PEAKS peptide Area",
    source_file = path
  )
}

dia_files <- list.files(
  file.path(
    project_root, "data", "PXD053474", "search_results", "extracted",
    "Subcellular-Whole-cell-lysates-Enriched-DIA-Search-files"
  ),
  pattern = "ptm-site_+Report\\.tsv$",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)
for (path in dia_files) {
  data <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  keep <- grepl("Lac|Lactyl", data$PTM.ModificationTitle, ignore.case = TRUE) &
    data$PTM.SiteAA == "K"
  probability_columns <- grep("PTM\\.SiteProbability$", names(data), value = TRUE)
  if (length(probability_columns)) {
    keep <- keep & rowSums(
      sapply(data[probability_columns], function(x) safe_numeric(x) > 0),
      na.rm = TRUE
    ) > 0
  }
  data <- data[keep, , drop = FALSE]
  columns <- grep("PTM\\.Quantity$", names(data), value = TRUE)
  add_quant(
    data, "PXD053474", "HCT116",
    feature_values = accession_feature(data$PTM.ProteinId),
    target_values = match_target_accession(data$PTM.ProteinId),
    sample_columns = columns,
    sample_labels = paste(
      "DIA",
      basename(dirname(path)),
      seq_along(columns),
      sep = "__"
    ),
    measurement_level = "Kla_site",
    quant_field = "Spectronaut PTM.Quantity",
    source_file = path
  )
}

# Author supplementary tables with explicit site quantities.
path <- file.path(
  project_root,
  "data/PXD050470/supplementary/prca2331-sup-0005-tables3.xlsx"
)
data <- read_excel(path, sheet = 1, skip = 12)
keep <- safe_numeric(data$`Localization probability`) > 0
data <- data[keep, , drop = FALSE]
columns <- grep("^Intensity_H", names(data), value = TRUE)
add_quant(
  data, "PXD050470", "human hippocampus",
  feature_values = accession_feature(data$`Proteins accession`),
  target_values = match_target_accession(data$`Proteins accession`),
  sample_columns = columns,
  sample_labels = sub("^Intensity_", "", columns),
  measurement_level = "Kla_site",
  quant_field = "author supplementary site Intensity",
  source_file = path
)
rm(data)

path <- file.path(
  project_root,
  "data/PXD054919/supplementary/41419_2025_8113_MOESM2_ESM.xlsx"
)
data <- read_excel(path, skip = 1)
columns <- grep("^Intensity A549-", names(data), value = TRUE)
add_quant(
  data, "PXD054919", "A549",
  feature_values = accession_feature(data$`Protein accession`),
  target_values = match_target_accession(data$`Protein accession`),
  sample_columns = columns,
  sample_labels = sub("^Intensity ", "", columns),
  measurement_level = "Kla_site",
  quant_field = "author site Intensity",
  source_file = path
)
rm(data)

path <- file.path(
  project_root,
  "data/PXD064912/supplementary/europepmc/mmc1.xlsx"
)
data <- read_excel(path, skip = 1)
keep <- tolower(data$PTM.ModificationTitle) == "lactylation" &
  data$PTM.SiteAA == "K"
data <- data[keep, , drop = FALSE]
columns <- grep("^PTM\\.Quantity", names(data), value = TRUE)
add_quant(
  data, "PXD064912", "human sperm",
  feature_values = accession_feature(data$PTM.ProteinId),
  target_values = match_target_accession(data$PTM.ProteinId),
  sample_columns = columns,
  sample_labels = sub("^PTM\\.Quantity\\[|\\]$", "", columns),
  measurement_level = "Kla_site",
  quant_field = "author PTM.Quantity",
  source_file = path
)
rm(data)

path <- file.path(
  project_root,
  "data/PXD070007/search_results/SA206LPLaB1_Annotation.xlsx"
)
data <- read_excel(path, sheet = "Annotation_Combine")
for (spec in list(
  list("glioblastoma stem cells", c("G2907", "G3028", "G3264", "GSC23", "MES28", "RKI")),
  list("neural stem cells", c("ENSA", "HMP1"))
)) {
  keep <- safe_numeric(data$`Localization probability`) > 0
  subset <- data[keep, , drop = FALSE]
  add_quant(
    subset, "PXD070007", spec[[1]],
    feature_values = accession_feature(subset$`Protein accession`),
    target_values = match_target_accession(subset$`Protein accession`),
    sample_columns = spec[[2]],
    measurement_level = "Kla_site",
    quant_field = "author annotation sample intensity",
    source_file = path
  )
}
rm(data)

path <- file.path(
  project_root,
  "data/PXD075377/search_results/extracted/2-Basic_analysis/MS_identified_information.txt"
)
data <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
for (spec in list(
  list("HCC", "Intensity HCC", "HCC"),
  list("adjacent liver", "Intensity Control", "Control")
)) {
  keep <- safe_numeric(data$`Localization probability`) > 0
  subset <- data[keep, , drop = FALSE]
  add_quant(
    subset, "PXD075377", spec[[1]],
    feature_values = accession_feature(subset$`Protein accession`),
    target_values = match_target_accession(subset$`Protein accession`),
    sample_columns = spec[[2]],
    sample_labels = spec[[3]],
    measurement_level = "Kla_site",
    quant_field = "author site Intensity",
    source_file = path
  )
}
rm(data)

# Spectronaut modified-peptide reports. Read only required columns from large files.
extract_spectronaut_precursor <- function(path, pxd, sample_group, require_site_confidence) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("data.table is required for large Spectronaut reports")
  }
  header <- names(data.table::fread(path, nrows = 0, data.table = FALSE))
  columns <- intersect(
    c(
      "Protein.Group", "Genes", "Modified.Sequence", "Q.Value", "PG.Q.Value",
      "PTM.Site.Confidence", "Run", "Precursor.Quantity"
    ),
    header
  )
  data <- data.table::fread(
    path, select = columns, showProgress = FALSE, data.table = FALSE
  )
  keep <- grepl("K\\(UniMod:378\\)", data$Modified.Sequence) &
    safe_numeric(data$Q.Value) <= 0.01 &
    safe_numeric(data$PG.Q.Value) <= 0.01
  if (require_site_confidence && "PTM.Site.Confidence" %in% names(data)) {
    keep <- keep & safe_numeric(data$PTM.Site.Confidence) > 0
  }
  data <- data[keep, , drop = FALSE]
  runs <- unique(data$Run)
  for (run in runs) {
    subset <- data[data$Run == run, , drop = FALSE]
    add_quant(
      subset, pxd, sample_group,
      feature_values = accession_feature(subset$Protein.Group),
      target_values = match_target_accession(subset$Protein.Group),
      sample_columns = "Precursor.Quantity",
      sample_labels = run,
      measurement_level = "modified_precursor",
      quant_field = "Spectronaut Precursor.Quantity",
      source_file = path
    )
  }
  rm(data)
  gc(verbose = FALSE)
}

extract_spectronaut_precursor(
  file.path(project_root, "data/PXD055230/search_results/LaIP_HSV1_DIA.tsv"),
  "PXD055230", "human fibroblasts mock and HCMV or HSV-1", TRUE
)
extract_spectronaut_precursor(
  file.path(project_root, "data/PXD057709/search_results/LaIP_report.tsv"),
  "PXD057709", "human fibroblasts mock and HCMV", FALSE
)

# Spectronaut PTM site report with condition labels.
path <- file.path(
  project_root,
  "data/PXD066054/search_results/extracted/PLa/XB08700B1DPLa_-PTMSiteReport.tsv"
)
data <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
keep <- data$PTM.ModificationTitle == "L-Lac(K)" &
  data$PTM.SiteAA == "K" &
  safe_numeric(data$PTM.SiteProbability) > 0
data <- data[keep, , drop = FALSE]
for (spec in list(
  list("BPH", "^NAT"),
  list("prostate cancer", "^PCa")
)) {
  subset <- data[grepl(spec[[2]], data$R.Condition), , drop = FALSE]
  for (condition in unique(subset$R.Condition)) {
    condition_data <- subset[subset$R.Condition == condition, , drop = FALSE]
    add_quant(
      condition_data, "PXD066054", spec[[1]],
      feature_values = accession_feature(condition_data$PTM.ProteinId),
      target_values = match_target_accession(condition_data$PTM.ProteinId),
      sample_columns = "PTM.Quantity",
      sample_labels = condition,
      measurement_level = "Kla_site",
      quant_field = "Spectronaut PTM.Quantity",
      source_file = path
    )
  }
}
rm(data)

# Explicit Lac(K) PTM quantities with accession-based target mapping.
path <- file.path(
  project_root,
  "data/PXD066351/search_results/XB01472B1DPLa-MSstats_Input.csv"
)
data <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
keep <- data$PTM.ModificationTitle == "Lac (K)" &
  data$PTM.SiteAA == "K" &
  safe_numeric(data$PTM.SiteProbability) > 0 &
  safe_numeric(data$PTM.Quantity) > 0
data <- data[keep, , drop = FALSE]
for (condition in unique(data$R.Condition)) {
  condition_data <- data[data$R.Condition == condition, , drop = FALSE]
  add_quant(
    condition_data, "PXD066351",
    "HCT116 control and Roseburia co-culture",
    feature_values = accession_feature(condition_data$PTM.ProteinId),
    target_values = match_target_accession(condition_data$PTM.ProteinId),
    sample_columns = "PTM.Quantity",
    sample_labels = condition,
    measurement_level = "Kla_site",
    quant_field = "MSstats Lac(K) PTM.Quantity",
    source_file = path
  )
}
rm(data)

# HUVEC Kla precursor matrix. UniMod:2114 is L-lactyllysine.
path <- file.path(
  project_root,
  paste0(
    "data/PXD073311/search_results/extracted_pairing/",
    "IPX0015307003_Database_search_result/Database_search_result/",
    "report.pr_matrix.tsv"
  )
)
data <- read.delim(
  path,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = ""
)
data <- data[
  grepl("K\\(UniMod:2114\\)", data$Modified.Sequence),
  ,
  drop = FALSE
]
columns <- names(data)[grepl("\\.raw$", names(data), ignore.case = TRUE)]
labels <- basename(gsub("\\\\", "/", columns))
labels <- sub("^FWMS20253547_02_", "", labels)
labels <- sub("\\.raw$", "", labels, ignore.case = TRUE)
add_quant(
  data, "PXD073311", "HUVEC control and Pg infection",
  feature_values = accession_feature(data$Protein.Group),
  target_values = match_target_accession(data$Protein.Group),
  sample_columns = columns,
  sample_labels = labels,
  measurement_level = "modified_precursor",
  quant_field = "Spectronaut precursor matrix quantity",
  source_file = path
)
rm(data)

# AC16 TMT Kla peptide supplementary table.
path <- file.path(
  project_root,
  "data/PXD075014/supplementary/Table2.XLSX"
)
data <- read_excel(path, sheet = "Kla peptides")
data <- data[
  grepl("lactylation", data$Modifications, ignore.case = TRUE),
  ,
  drop = FALSE
]
columns <- grep("^Abundance: F1:", names(data), value = TRUE)
labels <- c(
  "hypoxia_126", "hypoxia_127", "hypoxia_128",
  "control_129", "control_130", "control_131"
)
add_quant(
  data, "PXD075014", "AC16 control and hypoxia",
  feature_values = accession_feature(data$`Master Protein Accessions`),
  target_values = match_target_accession(data$`Master Protein Accessions`),
  sample_columns = columns,
  sample_labels = labels,
  measurement_level = "modified_peptide",
  quant_field = "author TMT peptide abundance",
  source_file = path
)
rm(data)

# PXD078013 has peptide evidence intensity and explicit Kla site IDs.
path <- file.path(project_root, "data/PXD078013/search_results/evidence.txt")
data <- read.delim(
  path, check.names = FALSE, stringsAsFactors = FALSE,
  quote = "", comment.char = ""
)
keep <- safe_numeric(data$`La (K)`) > 0 &
  !is.na(data$`La (K) site IDs`) &
  nzchar(as.character(data$`La (K) site IDs`)) &
  (is.na(data$Reverse) | data$Reverse != "+") &
  (is.na(data$`Potential contaminant`) | data$`Potential contaminant` != "+")
data <- data[keep, , drop = FALSE]
for (raw_file in unique(data$`Raw file`)) {
  subset <- data[data$`Raw file` == raw_file, , drop = FALSE]
  add_quant(
    subset, "PXD078013", "RKO WT and GSK3B KO",
    feature_values = accession_feature(subset$Proteins),
    target_values = match_target_accession(subset$Proteins),
    sample_columns = "Intensity",
    sample_labels = raw_file,
    measurement_level = "modified_peptide",
    quant_field = "MaxQuant evidence Intensity",
    source_file = path
  )
}
rm(data)

if (!length(quant_parts)) stop("No quantitative data were extracted")
quant_features <- bind_rows(quant_parts) |>
  group_by(
    PXD, SampleGroup, QuantSample, CanonicalFeature,
    TargetAccession, TargetGene,
    MeasurementLevel, QuantField
  ) |>
  summarise(
    Signal = sum(Signal, na.rm = TRUE),
    SourceFile = paste(sort(unique(SourceFile)), collapse = ";"),
    .groups = "drop"
  ) |>
  group_by(PXD, SampleGroup, QuantSample) |>
  mutate(
    Log2Signal = log2(Signal + 1),
    SampleMedianLog2 = median(Log2Signal, na.rm = TRUE),
    SampleNormalizedLog2 = Log2Signal - SampleMedianLog2,
    WithinSamplePercentile = {
      values <- Log2Signal
      if (length(values) == 1) {
        100
      } else {
        100 * (rank(values, ties.method = "average") - 1) /
          (length(values) - 1)
      }
    }
  ) |>
  ungroup()

sample_registry <- quant_features |>
  group_by(PXD, SampleGroup, QuantSample) |>
  summarise(
    SampleMedianLog2 = first(SampleMedianLog2),
    .groups = "drop"
  )

target_sample_grid <- sample_registry |>
  crossing(
    regulators |>
      distinct(GeneSymbol, BaseAccession) |>
      rename(RegulatorBaseAccession = BaseAccession)
  ) |>
  left_join(
    quant_features |>
      filter(!is.na(TargetAccession), TargetAccession %in% target_accessions) |>
      group_by(
        PXD, SampleGroup, QuantSample, TargetAccession, TargetGene
      ) |>
      summarise(
        Signal = sum(Signal),
        Log2Signal = log2(Signal + 1),
        SampleNormalizedLog2 = median(SampleNormalizedLog2),
        WithinSamplePercentile = max(WithinSamplePercentile),
        MeasurementLevel = paste(sort(unique(MeasurementLevel)), collapse = ";"),
        QuantField = paste(sort(unique(QuantField)), collapse = ";"),
        SourceFile = paste(sort(unique(SourceFile)), collapse = ";"),
        .groups = "drop"
      ),
    by = c(
      "PXD", "SampleGroup", "QuantSample",
      "RegulatorBaseAccession" = "TargetAccession"
    )
  ) |>
  mutate(
    Detected = !is.na(Signal) & Signal > 0,
    Signal = replace_na(Signal, 0),
    Log2Signal = replace_na(Log2Signal, 0),
    SampleNormalizedLog2 = ifelse(
      is.na(SampleNormalizedLog2),
      -SampleMedianLog2,
      SampleNormalizedLog2
    ),
    WithinSamplePercentile = replace_na(WithinSamplePercentile, 0)
  ) |>
  mutate(
    IdentityMatchMode = ifelse(
      Detected,
      "BaseAccession",
      "not_detected"
    )
  )

paired_sample_keys <- sample_catalog |>
  filter(!is.na(WholeProteomeDisplayRowOrder)) |>
  select(PXD, SampleGroup)
if (nrow(paired_sample_keys) != 33) {
  stop(
    "Expected 33 Kla sample groups with exact whole-proteome references, found ",
    nrow(paired_sample_keys)
  )
}
paired_target_sample_grid <- target_sample_grid |>
  semi_join(paired_sample_keys, by = c("PXD", "SampleGroup"))

group_quant <- target_sample_grid |>
  group_by(PXD, SampleGroup, RegulatorBaseAccession, GeneSymbol) |>
  summarise(
    RelativeKlaPercentile = median(WithinSamplePercentile),
    DetectedSampleCount = sum(Detected),
    QuantSampleCount = n(),
    DetectedSampleFraction = mean(Detected),
    MedianLog2SignalDetected = ifelse(
      any(Detected),
      median(Log2Signal[Detected]),
      NA_real_
    ),
    MeasurementLevel = paste(sort(unique(na.omit(MeasurementLevel))), collapse = ";"),
    QuantField = paste(sort(unique(na.omit(QuantField))), collapse = ";"),
    SourceFile = paste(sort(unique(na.omit(SourceFile))), collapse = ";"),
    .groups = "drop"
  )

heatmap_data <- detection |>
  select(
    Role, GeneSymbol, PXD, SampleGroup, SampleGroupID, RowLabel,
    BiologicalMaterial, GeneLevelAuditStatus, Detected
  ) |>
  rename(IdentificationDetected = Detected) |>
  left_join(
    regulators |>
      select(
        Role,
        GeneSymbol,
        BaseAccession,
        RegulatorDisplayName,
        RoleEntryOrder
      ) |>
      rename(RegulatorBaseAccession = BaseAccession),
    by = c("Role", "GeneSymbol")
  ) |>
  left_join(
    sample_catalog |>
      select(
        PXD,
        SampleGroup,
        Category,
        CategoryZh,
        CategoryEn,
        ComparisonRowOrder,
        WholeProteomeDisplayRowOrder,
        RowLabelEn
      ),
    by = c("PXD", "SampleGroup")
  ) |>
  left_join(
    group_quant,
    by = c("PXD", "SampleGroup", "GeneSymbol", "RegulatorBaseAccession")
  ) |>
  mutate(
    QuantificationAvailable = !is.na(QuantSampleCount),
    QuantState = case_when(
      !QuantificationAvailable ~ "无可比定量值",
      DetectedSampleCount == 0 ~ "定量样本中未检出",
      TRUE ~ "有定量信号"
    ),
    RelativeKlaPercentile = ifelse(
      QuantificationAvailable,
      replace_na(RelativeKlaPercentile, 0),
      NA_real_
    )
  )

heatmap_data <- heatmap_data |>
  semi_join(paired_sample_keys, by = c("PXD", "SampleGroup")) |>
  arrange(
    match(
      Category,
      category_order
    ),
    ComparisonRowOrder,
    RoleEntryOrder
  )

write.csv(
  paired_target_sample_grid,
  file.path(table_dir, "kla_regulator_intensity_sample_level_long.csv"),
  row.names = FALSE,
  na = ""
)
write.csv(
  heatmap_data,
  file.path(table_dir, "kla_regulator_normalized_intensity_long.csv"),
  row.names = FALSE,
  na = ""
)
id_mapping_audit <- target_sample_grid |>
  semi_join(paired_sample_keys, by = c("PXD", "SampleGroup")) |>
  group_by(GeneSymbol, RegulatorBaseAccession) |>
  summarise(
    QuantitativeSampleRows = n(),
    DetectedSampleRows = sum(Detected),
    IdentityMatchMode = "UniProt_BaseAccession_only",
    GeneSymbolFallbackCount = 0L,
    MappingSource = "UniProtKB reviewed Homo sapiens",
    .groups = "drop"
  ) |>
  arrange(GeneSymbol)
write.csv(
  id_mapping_audit,
  file.path(table_dir, "kla_regulator_intensity_id_mapping_audit.csv"),
  row.names = FALSE,
  na = ""
)

quant_summary <- quant_features |>
  group_by(PXD, SampleGroup) |>
  summarise(
    QuantSampleCount = n_distinct(QuantSample),
    QuantifiedFeatureCount = n_distinct(CanonicalFeature),
    QuantifiedTargetGeneCount = n_distinct(TargetGene[!is.na(TargetGene)]),
    MeasurementLevel = paste(sort(unique(MeasurementLevel)), collapse = ";"),
    QuantField = paste(sort(unique(QuantField)), collapse = ";"),
    SourceFile = paste(sort(unique(SourceFile)), collapse = ";"),
    .groups = "drop"
  )

unavailable_reason <- function(pxd, status) {
  if (pxd == "PXD037371") {
    return("位点表含TMT reporter intensity，但当前缺少三类临床组与通道的可靠映射，不能拆分到对应组织组。")
  }
  if (status == "无逐蛋白明细") {
    return("当前文件不能把定量值可靠关联到Kla蛋白身份。")
  }
  "未从当前证据文件中提取到可解释的Kla定量字段，需进一步人工核对。"
}

audit <- sample_catalog |>
  left_join(quant_summary, by = c("PXD", "SampleGroup")) |>
  mutate(
    定量可用 = !is.na(QuantSampleCount),
    定量样本数 = replace_na(QuantSampleCount, 0L),
    定量层级 = ifelse(定量可用, MeasurementLevel, ""),
    使用的定量字段 = ifelse(定量可用, QuantField, ""),
    定量Kla特征数 = replace_na(QuantifiedFeatureCount, 0L),
    定量调控蛋白数 = replace_na(QuantifiedTargetGeneCount, 0L),
    调控蛋白身份匹配 = "UniProt BaseAccession only",
    GeneSymbol回退数 = 0L,
    跨研究热图数值 = ifelse(
      定量可用,
      "每个定量样本内log2(signal+1)后按全部Kla特征计算百分位，再对本样本组取中位数",
      ""
    ),
    可做PXD内Z分数 = 定量可用 & 定量样本数 >= 2,
    严格配对分析纳入 = !is.na(WholeProteomeDisplayRowOrder),
    严格配对排除原因 = case_when(
      !is.na(WholeProteomeDisplayRowOrder) ~ "",
      定量可用 ~ "没有完全匹配且可审计的普通全蛋白参照；保留Kla原始审计，但不进入配对热图和后续对照分析。",
      TRUE ~ "Kla本身没有可用的逐组定量值；保留原始审计，但不进入配对分析。"
    ),
    限制或不可用原因 = ifelse(
      定量可用,
      case_when(
        grepl(";", 定量层级, fixed = TRUE) ~
          "同一PXD包含多个测量层级；保留来源并仅用组内百分位汇总，不直接比较原始强度。",
        TRUE ~ "不同PXD的原始数值不可直接比较；跨研究只解释为组内相对Kla信号。"
      ),
      mapply(unavailable_reason, PXD, GeneLevelAuditStatus)
    ),
    定量来源文件 = ifelse(定量可用, SourceFile, "")
  ) |>
  transmute(
    PXD,
    样本组 = SampleGroup,
    行标签 = RowLabel,
    材料类型 = BiologicalMaterial,
    逐蛋白身份审计状态 = GeneLevelAuditStatus,
    定量可用,
    定量样本数,
    定量层级,
    使用的定量字段,
    定量Kla特征数,
    定量调控蛋白数,
    调控蛋白身份匹配,
    GeneSymbol回退数,
    跨研究热图数值,
    可做PXD内Z分数,
    严格配对分析纳入,
    严格配对排除原因,
    限制或不可用原因,
    定量来源文件
  )
write.csv(
  audit,
  file.path(table_dir, "kla_regulator_intensity_availability_audit.csv"),
  row.names = FALSE,
  na = ""
)
axis_audit <- sample_catalog |>
  filter(!is.na(WholeProteomeDisplayRowOrder)) |>
  left_join(
    audit |>
      transmute(
        PXD,
        SampleGroup = 样本组,
        QuantificationAvailable = 定量可用
    ),
    by = c("PXD", "SampleGroup")
  ) |>
  left_join(
    whole_heatmap_rows |>
      select(
        WholeProteomeDisplayRowOrder = HeatmapDisplayRowOrder,
        WholeProteomeRowLabel = RowLabel,
        WholeProteomeReferencePXD = ReferencePXD
      ),
    by = "WholeProteomeDisplayRowOrder"
  ) |>
  arrange(WholeProteomeDisplayRowOrder, RowOrder) |>
  mutate(KlaDisplayRowOrder = row_number()) |>
  transmute(
    PXD,
    SampleGroup,
    Category,
    CategoryZh,
    CategoryEn,
    KlaDisplayRowOrder,
    ComparisonRowOrder,
    WholeProteomeDisplayRowOrder,
    WholeProteomeReferencePXD,
    WholeProteomeRowLabel,
    QuantificationAvailable,
    RowLabel = RowLabelZh,
    RowLabelZh,
    RowLabelEn,
    OrderAlignedToWholeProteome = TRUE,
    OriginalRowOrder = RowOrder
  )
write.csv(
  axis_audit,
  file.path(table_dir, "kla_regulator_heatmap_axis_order.csv"),
  row.names = FALSE,
  na = ""
)
write.csv(
  axis_audit,
  file.path(
    table_dir,
    "kla_vs_whole_proteome_heatmap_axis_alignment.csv"
  ),
  row.names = FALSE,
  na = ""
)

pure_white_audit <- heatmap_data |>
  group_by(PXD, SampleGroup) |>
  summarise(
    QuantificationAvailable = all(QuantificationAvailable),
    PureWhite = QuantificationAvailable &&
      all(is.na(RelativeKlaPercentile) | RelativeKlaPercentile == 0),
    DetectedRegulatorCount = sum(DetectedSampleCount > 0, na.rm = TRUE),
    DetectedSampleRows = sum(DetectedSampleCount, na.rm = TRUE),
    MaxRelativeKlaPercentile = ifelse(
      all(is.na(RelativeKlaPercentile)),
      NA_real_,
      max(RelativeKlaPercentile, na.rm = TRUE)
    ),
    .groups = "drop"
  ) |>
  left_join(
    audit |>
      transmute(
        PXD,
        SampleGroup = 样本组,
        QuantifiedKlaFeatureCount = 定量Kla特征数,
        QuantificationSourceFile = 定量来源文件
      ),
    by = c("PXD", "SampleGroup")
  ) |>
  mutate(
    Explanation = case_when(
      !QuantificationAvailable ~ "excluded_or_no_usable_quantification",
      !PureWhite ~ "has_positive_regulator_signal",
      DetectedSampleRows == 0 ~
        "no_target_regulator_accession_hit_in_Kla_source",
      TRUE ~ "sparse_detection_and_group_median_zero"
    )
  ) |>
  arrange(match(PXD, sample_catalog$PXD), match(SampleGroup, sample_catalog$SampleGroup))
write.csv(
  pure_white_audit,
  file.path(table_dir, "kla_regulator_intensity_pure_white_audit.csv"),
  row.names = FALSE,
  na = ""
)

write.csv(
  audit |>
    filter(!严格配对分析纳入) |>
    select(
      PXD,
      样本组,
      定量可用,
      严格配对分析纳入,
      严格配对排除原因,
      限制或不可用原因,
      定量来源文件
    ),
  file.path(table_dir, "kla_regulator_intensity_plot_exclusions.csv"),
  row.names = FALSE,
  na = ""
)

# Main cross-study heatmap.
role_levels <- c("Writer", "Eraser", "Writer-Eraser", "Reader")
role_labels_zh <- c(
  Writer = "写入因子",
  Eraser = "去除因子",
  `Writer-Eraser` = "写入/去除因子",
  Reader = "读取因子"
)
role_labels_en <- c(
  Writer = "Writer",
  Eraser = "Eraser",
  `Writer-Eraser` = "Writer-Eraser",
  Reader = "Reader"
)
gene_levels <- regulators |>
  arrange(factor(Role, levels = role_levels), RoleEntryOrder) |>
  pull(RegulatorDisplayName) |>
  unique()
plot_sample_catalog <- sample_catalog |>
  inner_join(
    axis_audit |>
      filter(QuantificationAvailable) |>
      select(PXD, SampleGroup),
    by = c("PXD", "SampleGroup")
  ) |>
  left_join(
    axis_audit |>
      select(PXD, SampleGroup, KlaDisplayRowOrder),
    by = c("PXD", "SampleGroup")
  ) |>
  arrange(KlaDisplayRowOrder)
row_levels_zh <- rev(plot_sample_catalog$RowLabelZh)
row_levels_en <- rev(plot_sample_catalog$RowLabelEn)
plot_data_base <- heatmap_data |>
  filter(QuantificationAvailable) |>
  mutate(
    Role = factor(Role, levels = role_levels),
    RegulatorDisplayName = factor(
      RegulatorDisplayName,
      levels = gene_levels
    )
  )

make_main_plot <- function(language = c("zh", "en")) {
  language <- match.arg(language)
  is_zh <- language == "zh"
  category_labels <- if (is_zh) category_labels_zh else category_labels_en
  role_labels <- if (is_zh) role_labels_zh else role_labels_en
  row_levels <- if (is_zh) row_levels_zh else row_levels_en
  data <- plot_data_base |>
    mutate(
      CategoryLabel = factor(
        unname(category_labels[Category]),
        levels = unname(category_labels[category_order])
      ),
      RoleLabel = factor(
        unname(role_labels[as.character(Role)]),
        levels = unname(role_labels[role_levels])
      ),
      PlotRowLabel = factor(
        if (is_zh) RowLabel else RowLabelEn,
        levels = row_levels
      )
    )

  ggplot(
    data,
    aes(RegulatorDisplayName, PlotRowLabel, fill = RelativeKlaPercentile)
  ) +
    geom_tile(color = "white", linewidth = 0.22) +
    facet_grid(
      CategoryLabel ~ RoleLabel,
      scales = "free",
      space = "free"
    ) +
    scale_fill_gradientn(
      colours = c("#FFFFFF", "#FFF3E0", "#FDBB84", "#FC8D59", "#B2182B"),
      values = scales::rescale(c(0, 20, 50, 80, 100)),
      limits = c(0, 100),
      na.value = "#D9D9D9",
      name = if (is_zh) "组内Kla信号\n百分位" else "Within-sample Kla\npercentile"
    ) +
    labs(
      title = if (is_zh) {
        "乳酸化调控蛋白在细胞系与组织中的相对Kla信号"
      } else {
        "Relative Kla signals of lactylation regulators across tissues and cell lines"
      },
      subtitle = if (is_zh) {
        paste0(
          "严格配对范围为33个Kla样本组；行顺序完全跟随普通全蛋白热图，",
          "并按正常组织、癌症组织、正常细胞、癌症细胞分区。颜色由白色向暖色递增。"
        )
      } else {
        paste0(
          "The strict paired scope contains 33 Kla sample groups. Row order follows the whole-proteome heatmap and is divided into ",
          "normal tissues, cancer tissues, normal cells, and cancer cells. Warmer colors indicate higher within-sample percentiles."
        )
      },
      x = NULL,
      y = NULL,
      caption = if (is_zh) {
        paste0(
          "每个定量样本内按全部Kla特征计算百分位，再在样本组内取中位数；",
          "0表示未检出。身份仅按人源UniProt BaseAccession匹配，GeneSymbol只用于显示。"
        )
      } else {
        paste0(
          "Percentiles are calculated among all Kla features within each quantitative sample and aggregated by the median; 0 denotes not detected. ",
          "Protein matching uses human UniProt BaseAccession only; gene symbols are display-only."
        )
      }
    ) +
    theme_minimal(base_size = 8.5, base_family = plot_font) +
    theme(
      panel.grid = element_blank(),
      strip.text.x = element_text(face = "bold", size = 9),
      strip.text.y = element_text(face = "bold", size = 8),
      strip.background = element_rect(fill = "#F2F2F2", color = NA),
      axis.text.x = element_text(angle = 55, hjust = 1, vjust = 1, size = 7),
      axis.text.y = element_text(size = 7.2),
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 8.2, lineheight = 1.12),
      plot.caption = element_text(size = 7.5, hjust = 0, color = "#5C626A"),
      legend.position = "bottom",
      legend.key.width = grid::unit(35, "mm"),
      plot.margin = margin(8, 10, 8, 8)
    )
}

main_plot_zh <- make_main_plot("zh")
main_plot_en <- make_main_plot("en")
save_heatmap <- function(plot, stem) {
  ggsave(
    file.path(figure_dir, paste0(stem, ".png")),
    plot,
    width = 15.5,
    height = 11.5,
    dpi = 320,
    bg = "white"
  )
  ggsave(
    file.path(figure_dir, paste0(stem, ".pdf")),
    plot,
    width = 15.5,
    height = 11.5,
    device = cairo_pdf,
    bg = "white"
  )
}
save_heatmap(
  main_plot_zh,
  "kla_regulator_cross_study_relative_intensity_heatmap"
)
save_heatmap(
  main_plot_zh,
  "kla_regulator_cross_study_relative_intensity_heatmap_zh"
)
save_heatmap(
  main_plot_en,
  "kla_regulator_cross_study_relative_intensity_heatmap_en"
)

# Within-PXD gene-wise z-scores are valid only across samples from the same PXD.
zscore_data <- target_sample_grid |>
  semi_join(paired_sample_keys, by = c("PXD", "SampleGroup")) |>
  left_join(
    regulators |>
      distinct(
        GeneSymbol,
        RegulatorBaseAccession = BaseAccession,
        RegulatorDisplayName
      ),
    by = c("GeneSymbol", "RegulatorBaseAccession")
  ) |>
  group_by(PXD, GeneSymbol) |>
  mutate(
    WithinPXDZ = {
      values <- SampleNormalizedLog2
      value_sd <- sd(values, na.rm = TRUE)
      if (
        n_distinct(QuantSample) >= 2 &&
          sum(Detected, na.rm = TRUE) >= 2 &&
          is.finite(value_sd) &&
          value_sd > 0
      ) {
        as.numeric(scale(values))
      } else {
        rep(NA_real_, length(values))
      }
    }
  ) |>
  ungroup() |>
  filter(!is.na(WithinPXDZ))
write.csv(
  zscore_data,
  file.path(table_dir, "kla_regulator_within_pxd_zscore_long.csv"),
  row.names = FALSE,
  na = ""
)

z_pxd <- zscore_data |>
  group_by(PXD) |>
  summarise(
    Samples = n_distinct(QuantSample),
    Genes = n_distinct(GeneSymbol),
    .groups = "drop"
  ) |>
  filter(Samples >= 2, Genes >= 2) |>
  pull(PXD)

pdf_path <- file.path(figure_dir, "kla_regulator_within_pxd_zscore_heatmaps.pdf")
grDevices::cairo_pdf(pdf_path, width = 12, height = 8.5, onefile = TRUE)
for (pxd in z_pxd) {
  data <- zscore_data |>
    filter(PXD == pxd) |>
    mutate(
      RegulatorDisplayName = factor(
        RegulatorDisplayName,
        levels = gene_levels
      ),
      QuantSample = factor(QuantSample, levels = rev(unique(QuantSample)))
    )
  p <- ggplot(
    data,
    aes(
      RegulatorDisplayName,
      QuantSample,
      fill = pmax(-2.5, pmin(2.5, WithinPXDZ))
    )
  ) +
    geom_tile(color = "white", linewidth = 0.25) +
    scale_fill_gradient2(
      low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
      midpoint = 0, limits = c(-2.5, 2.5),
      name = "PXD内Z分数"
    ) +
    labs(
      title = paste0(pxd, "：乳酸化调控蛋白的PXD内相对变化"),
      subtitle = "log2(signal+1)后先做样本中位数中心化，再对每个蛋白跨本PXD样本计算Z分数；不跨PXD比较。",
      x = NULL,
      y = NULL
    ) +
    theme_minimal(base_size = 9, base_family = plot_font) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 55, hjust = 1, size = 7),
      axis.text.y = element_text(size = 7),
      plot.title = element_text(face = "bold"),
      legend.position = "bottom"
    )
  print(p)
}
grDevices::dev.off()

message(
  "Quantitative sample groups: ",
  sum(audit$定量可用),
  "/",
  nrow(audit),
  "; within-PXD heatmaps: ",
  length(z_pxd)
)
