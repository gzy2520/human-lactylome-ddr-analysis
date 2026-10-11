#!/usr/bin/env Rscript

# Build candidate-only sample-level Figure 1 inputs from the deposited processed
# tables. This script does not alter the publication workflow or its frozen
# tables. One observation is one biological sample, model, or experimental
# condition when that identity is recoverable from the source table.

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
})

set.seed(25)

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1L) {
  normalizePath(args[[1L]], mustWork = TRUE)
} else {
  normalizePath(".", mustWork = TRUE)
}
source_root <- if (length(args) >= 2L) {
  normalizePath(args[[2L]], mustWork = TRUE)
} else {
  stop(
    "Usage: Rscript R/candidate/prepare_sample_boxplot_inputs.R <project-root> <source-data-root>",
    call. = FALSE
  )
}

stop_if <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

require_file <- function(path) {
  stop_if(file.exists(path), paste0("Required source file is missing: ", path))
  path
}

base_accession <- function(values) {
  values <- trimws(as.character(values))
  values[is.na(values)] <- ""
  values <- sub("^.*\\|([^|]+)\\|.*$", "\\1", values)
  values <- sub("^([^|;]+)\\|.*$", "\\1", values)
  values <- sub("^NX_", "", values)
  sub("-[0-9]+$", "", values)
}

is_uniprot <- function(values) {
  grepl(
    "^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9](?:[A-Z0-9]{3}[0-9])?)$",
    values
  )
}

split_accessions <- function(values) {
  tokens <- unlist(strsplit(as.character(values), "[;,]"), use.names = FALSE)
  accessions <- base_accession(tokens)
  sort(unique(accessions[is_uniprot(accessions)]))
}

first_accession <- function(value) {
  accessions <- split_accessions(value)
  if (length(accessions)) accessions[[1L]] else ""
}

first_accession_current <- function(value) {
  values <- trimws(unlist(strsplit(as.character(value), "[;,]"), use.names = FALSE))
  values[is.na(values)] <- ""
  values <- sub("^(?:REV__|CON__)+", "", values, perl = TRUE)
  values <- sub("^(?:sp|tr)\\|", "", values, perl = TRUE)
  values <- sub("\\|.*$", "", values)
  values <- sub("^NX_", "", values)
  values <- sub("-[0-9]+$", "", values)
  values <- values[is_uniprot(values)]
  if (length(values)) values[[1L]] else ""
}

is_true_flag <- function(values) {
  values <- tolower(trimws(as.character(values)))
  values[is.na(values)] <- ""
  values %in% c("+", "1", "true", "yes", "y")
}

positive_numeric <- function(values) {
  values <- suppressWarnings(as.numeric(values))
  is.finite(values) & values > 0
}

detected_text <- function(values) {
  values <- trimws(as.character(values))
  !is.na(values) & nzchar(values) & values != "Not Found"
}

valid_maxquant_rows <- function(data) {
  keep <- rep(TRUE, nrow(data))
  for (column in c("Reverse", "Potential contaminant", "Contaminant", "Only identified by site")) {
    if (column %in% names(data)) keep <- keep & !is_true_flag(data[[column]])
  }
  keep[is.na(keep)] <- FALSE
  keep
}

records <- function(pxd, group, sample_id, accessions, source_mode, source_file, sample_class, observation_type = "sample") {
  accessions <- sort(unique(accessions[is_uniprot(accessions)]))
  if (!length(accessions)) {
    return(data.table(
      PXD = character(), SampleGroup = character(), SampleID = character(),
      BaseAccession = character(), SourceMode = character(), SourceFile = character(),
      SampleClass = character(), ObservationType = character()
    ))
  }
  data.table(
    PXD = pxd,
    SampleGroup = group,
    SampleID = sample_id,
    BaseAccession = accessions,
    SourceMode = source_mode,
    SourceFile = source_file,
    SampleClass = sample_class,
    ObservationType = observation_type
  )
}

read_maxquant_by_localization <- function(path, pxd, group, sample_map, accession_column = "Proteins") {
  data <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = "")
  stop_if(accession_column %in% names(data), paste0("Missing ", accession_column, " in ", path))
  base_keep <- valid_maxquant_rows(data)
  if ("id" %in% names(data)) base_keep <- base_keep & !is.na(data$id)
  output <- lapply(seq_len(nrow(sample_map)), function(index) {
    sample_id <- sample_map$SampleID[[index]]
    probability_column <- intersect(
      c(paste0("Localization prob ", sample_id), paste0("Localization prob", sample_id)),
      names(data)
    )
    if (!length(probability_column)) probability_column <- intersect("Localization prob", names(data))
    stop_if(length(probability_column) == 1L, paste0("Missing localization column for ", sample_id, " in ", path))
    keep <- base_keep & positive_numeric(data[[probability_column]])
    records(
      pxd, group, sample_id, split_accessions(data[[accession_column]][keep]),
      "deposited_sample_table", path, sample_map$SampleClass[[index]], "sample"
    )
  })
  rbindlist(output, fill = TRUE)
}

read_maxquant_by_identification <- function(path, pxd, group, sample_map, accession_column = "Protein", accession_fn = first_accession) {
  data <- fread(path, check.names = FALSE, showProgress = FALSE)
  stop_if(accession_column %in% names(data), paste0("Missing ", accession_column, " in ", path))
  base_keep <- valid_maxquant_rows(data)
  if ("Amino acid" %in% names(data)) base_keep <- base_keep & data[["Amino acid"]] == "K"
  output <- lapply(seq_len(nrow(sample_map)), function(index) {
    sample_id <- sample_map$SampleID[[index]]
    identification_column <- paste0("Identification type ", sample_id)
    stop_if(identification_column %in% names(data), paste0("Missing identification column for ", sample_id, " in ", path))
    keep <- base_keep & detected_text(data[[identification_column]])
    accessions <- vapply(data[[accession_column]][keep], accession_fn, character(1))
    records(
      pxd, group, sample_id, accessions[nzchar(accessions)],
      "deposited_sample_table", path, sample_map$SampleClass[[index]], "sample"
    )
  })
  rbindlist(output, fill = TRUE)
}

read_pd_lactyl_by_sample <- function(path, pxd, group, sample_map, accession_column = "Master Protein Accessions") {
  header <- names(read.delim(path, nrows = 0L, check.names = FALSE, stringsAsFactors = FALSE, quote = "\"", comment.char = ""))
  sample_columns <- vapply(sample_map$SampleID, function(sample_id) {
    found <- header[startsWith(header, "Found in Sample:")]
    hits <- found[endsWith(found, paste0(", ", sample_id)) | endsWith(found, paste0(", ", sample_id, "_1"))]
    stop_if(length(hits) == 1L, paste0("Expected one sample column for ", sample_id, " in ", path))
    hits[[1L]]
  }, character(1))
  data <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE, quote = "\"", comment.char = "")
  keep_lactyl <- grepl("Lacty|Lactyl|La \\(K\\)", data$Modifications, ignore.case = TRUE)
  output <- lapply(seq_len(nrow(sample_map)), function(index) {
    keep <- keep_lactyl & detected_text(data[[sample_columns[[index]]]])
    records(
      pxd, group, sample_map$SampleID[[index]], split_accessions(data[[accession_column]][keep]),
      "deposited_sample_table", path, sample_map$SampleClass[[index]], "sample"
    )
  })
  rbindlist(output, fill = TRUE)
}

read_hippocampus_by_sample <- function(path, pxd, group, sample_map) {
  data <- as.data.table(read_excel(path, sheet = "Sheet1", skip = 12L))
  output <- lapply(seq_len(nrow(sample_map)), function(index) {
    column <- paste0("Intensity_", sample_map$SampleID[[index]])
    stop_if(column %in% names(data), paste0("Missing ", column, " in ", path))
    keep <- positive_numeric(data[[column]]) & !is.na(data$`Positions within proteins`)
    records(
      pxd, group, sample_map$SampleID[[index]], split_accessions(data$`Proteins accession`[keep]),
      "deposited_supplementary_table", path, sample_map$SampleClass[[index]], "sample"
    )
  })
  rbindlist(output, fill = TRUE)
}

read_spectronaut_kla_by_sample <- function(path, pxd, group, sample_map) {
  data <- fread(
    path,
    select = c("R.Condition", "PTM.ProteinId", "PTM.ModificationTitle", "PTM.SiteAA", "PTM.SiteProbability"),
    check.names = FALSE,
    showProgress = FALSE
  )
  output <- lapply(seq_len(nrow(sample_map)), function(index) {
    sample_id <- sample_map$SampleID[[index]]
    keep <- data$R.Condition == sample_id &
      data$PTM.ModificationTitle == "L-Lac(K)" &
      data$PTM.SiteAA == "K" &
      positive_numeric(data$PTM.SiteProbability)
    records(
      pxd, group, sample_id, split_accessions(data$PTM.ProteinId[keep]),
      "deposited_spectronaut_table", path, sample_map$SampleClass[[index]], "sample"
    )
  })
  rbindlist(output, fill = TRUE)
}

read_itraq_pool <- function(path, pxd, group, sample_id, sample_class, intensity_column) {
  data <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = "")
  keep <- positive_numeric(data[[intensity_column]]) & positive_numeric(data$`Localization probability`)
  records(pxd, group, sample_id, split_accessions(data$`Protein accession`[keep]),
    "deposited_group_union", path, sample_class, "pool")
}

read_pxd064912_by_sample <- function(path, pxd, group, sample_map) {
  data <- as.data.table(read_excel(path, sheet = "Sheet1", skip = 1L))
  probability_columns <- grep("^PTM.SiteProbability", names(data), value = TRUE)
  sample_labels <- trimws(gsub("(^PTM.SiteProbability\\[|\\]$)", "", probability_columns))
  output <- lapply(seq_len(nrow(sample_map)), function(index) {
    sample_id <- sample_map$SampleID[[index]]
    probability_column <- probability_columns[sample_labels == sample_id]
    stop_if(length(probability_column) == 1L, paste0("Missing sperm probability column for ", sample_id))
    keep <- tolower(trimws(as.character(data$PTM.ModificationTitle))) == "lactylation" &
      data$PTM.SiteAA == "K" & positive_numeric(data[[probability_column]])
    records(
      pxd, group, sample_id, split_accessions(data$PTM.ProteinId[keep]),
      "deposited_supplementary_table", path, sample_map$SampleClass[[index]], "sample"
    )
  })
  rbindlist(output, fill = TRUE)
}

read_pxd066351_by_sample <- function(path, pxd, group, sample_map) {
  data <- fread(path, check.names = FALSE, showProgress = FALSE)
  output <- lapply(seq_len(nrow(sample_map)), function(index) {
    sample_id <- sample_map$SampleID[[index]]
    keep <- data$R.Condition == sample_id &
      data$PTM.ModificationTitle == "Lac (K)" &
      data$PTM.SiteAA == "K" &
      positive_numeric(data$PTM.SiteProbability) &
      positive_numeric(data$PTM.Quantity)
    records(
      pxd, group, sample_id, split_accessions(data$PTM.ProteinId[keep]),
      "deposited_spectronaut_table", path, sample_map$SampleClass[[index]], "condition"
    )
  })
  rbindlist(output, fill = TRUE)
}

read_author_intensity_by_sample <- function(path, pxd, group, sample_map) {
  data <- as.data.table(read_excel(path, sheet = 1L, skip = 1L))
  output <- lapply(seq_len(nrow(sample_map)), function(index) {
    sample_id <- sample_map$SampleID[[index]]
    column <- paste0("Intensity ", sample_id)
    stop_if(column %in% names(data), paste0("Missing ", column, " in ", path))
    keep <- positive_numeric(data[[column]])
    records(
      pxd, group, sample_id, split_accessions(data$`Protein accession`[keep]),
      "deposited_supplementary_table", path, sample_map$SampleClass[[index]], "sample"
    )
  })
  rbindlist(output, fill = TRUE)
}

read_pxd063266_by_sample <- function(path, pxd, group, sample_map) {
  data <- as.data.table(read_excel(path, sheet = "Lactyl (K)Sites"))
  base_keep <- valid_maxquant_rows(data)
  if ("id" %in% names(data)) base_keep <- base_keep & !is.na(data$id)
  output <- lapply(seq_len(nrow(sample_map)), function(index) {
    sample_id <- sample_map$SampleID[[index]]
    probability_column <- paste0("Localization prob ", sample_id)
    stop_if(probability_column %in% names(data), paste0("Missing ", probability_column, " in ", path))
    keep <- base_keep & positive_numeric(data[[probability_column]])
    records(
      pxd, group, sample_map$SampleID[[index]], split_accessions(data$Proteins[keep]),
      "deposited_supplementary_table", path, sample_map$SampleClass[[index]], "sample"
    )
  })
  rbindlist(output, fill = TRUE)
}

read_pxd070007_by_model <- function(path, pxd, group, sample_map) {
  data <- as.data.table(read_excel(path, sheet = "Annotation_Combine"))
  base_keep <- positive_numeric(data$`Localization probability`)
  output <- lapply(seq_len(nrow(sample_map)), function(index) {
    sample_id <- sample_map$SampleID[[index]]
    stop_if(sample_id %in% names(data), paste0("Missing model column ", sample_id, " in ", path))
    keep <- base_keep & positive_numeric(data[[sample_id]])
    records(
      pxd, group, sample_id, split_accessions(data$`Protein accession`[keep]),
      "deposited_supplementary_table", path, sample_map$SampleClass[[index]], "model"
    )
  })
  rbindlist(output, fill = TRUE)
}

read_pxd073311_by_raw_file <- function(path, pxd, group, sample_map) {
  data <- fread(path, check.names = FALSE, showProgress = FALSE)
  raw_columns <- names(data)[grepl("A[06]h_[123]\\.raw$", names(data))]
  output <- lapply(seq_len(nrow(sample_map)), function(index) {
    sample_id <- sample_map$SampleID[[index]]
    raw_column <- raw_columns[grepl(paste0("_", sample_id, "\\.raw$"), raw_columns)]
    stop_if(length(raw_column) == 1L, paste0("Missing HUVEC raw-file column for ", sample_id))
    keep <- grepl("K\\(UniMod:2114\\)", data$Modified.Sequence) & positive_numeric(data[[raw_column]])
    records(
      pxd, group, sample_id, split_accessions(data$Protein.Group[keep]),
      "deposited_spectronaut_matrix", path, sample_map$SampleClass[[index]], "condition"
    )
  })
  rbindlist(output, fill = TRUE)
}

read_maxquant_single_group <- function(path, pxd, group, sample_id, sample_class) {
  data <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = "")
  keep <- valid_maxquant_rows(data)
  if ("id" %in% names(data)) keep <- keep & !is.na(data$id)
  if ("Localization prob" %in% names(data)) keep <- keep & positive_numeric(data$`Localization prob`)
  records(pxd, group, sample_id, split_accessions(data$Proteins[keep]),
    "deposited_group_union", path, sample_class, "dataset_union")
}

sample_map <- function(pxd, group, ids, classes = rep("replicate", length(ids))) {
  data.table(PXD = pxd, SampleGroup = group, SampleID = ids, SampleClass = classes)
}

source_rel <- function(path) {
  normalized <- normalizePath(path, mustWork = FALSE)
  source_prefix <- paste0(source_root, "/")
  project_prefix <- paste0(project_root, "/")
  if (startsWith(normalized, source_prefix)) {
    sub(source_prefix, "", normalized, fixed = TRUE)
  } else if (startsWith(normalized, project_prefix)) {
    sub(project_prefix, "", normalized, fixed = TRUE)
  } else {
    normalized
  }
}

input_dir <- file.path(project_root, "data", "publication_input")
candidate_dir <- file.path(project_root, "data", "candidate")
dir.create(candidate_dir, recursive = TRUE, showWarnings = FALSE)

groups <- fread(file.path(input_dir, "group_summary_30.csv"), check.names = FALSE)
design <- fread(file.path(candidate_dir, "group_sample_design.csv"), check.names = FALSE, na.strings = c("", "NA"), fill = TRUE)
membership <- fread(file.path(input_dir, "kla_protein_membership_30.csv"), check.names = FALSE)
stop_if(nrow(groups) == 30L, "Frozen release must contain exactly 30 groups.")
stop_if(nrow(design) == 30L, "Sample design must contain exactly 30 groups.")
groups[, GroupKey := paste(PXD, SampleGroup, sep = "__")]
design[, GroupKey := paste(PXD, SampleGroup, sep = "__")]
stop_if(!anyDuplicated(groups$GroupKey), "Frozen groups are not unique.")
stop_if(!anyDuplicated(design$GroupKey), "Sample design rows are not unique.")
stop_if(setequal(groups$GroupKey, design$GroupKey), "Sample design does not cover the frozen 30-group scope.")

sample_records <- list()
registry <- list()
add_source <- function(pxd, group, sample_ids, classes, parser, path) {
  registry[[length(registry) + 1L]] <<- data.table(
    PXD = pxd, SampleGroup = group, SampleID = sample_ids,
    SampleClass = classes, Parser = parser, SourceFile = source_rel(path)
  )
}
add_records <- function(value) {
  sample_records[[length(sample_records) + 1L]] <<- value
}

# Samples with a direct per-sample signal in the deposited table.
path <- require_file(file.path(source_root, "PXD036307/search_results/extracted/txt/La (K)Sites.txt"))
map <- sample_map("PXD036307", "normal human lung", c("PTB340", "PTB342", "PTB344", "PTB346", "PTB364", "PTB372"))
add_records(read_maxquant_by_localization(path, map$PXD[[1L]], map$SampleGroup[[1L]], map, "Proteins"))
add_source(map$PXD[[1L]], map$SampleGroup[[1L]], map$SampleID, map$SampleClass, "MaxQuant site table / sample localization", path)

path <- require_file(file.path(source_root, "PXD046800/search_results/HFX2_LFQ_QB001_Lacty_PeptideGroups.txt"))
for (group in c("hypertrophic scar", "adjacent skin")) {
  ids <- if (group == "hypertrophic scar") paste0("HSP", 1:4) else paste0("NSP", 1:4)
  map <- sample_map("PXD046800", group, ids)
  add_records(read_pd_lactyl_by_sample(path, "PXD046800", group, map))
  add_source("PXD046800", group, ids, map$SampleClass, "Proteome Discoverer lactylated peptide groups", path)
}

path <- require_file(file.path(source_root, "PXD050470/supplementary/prca2331-sup-0005-tables3.xlsx"))
map <- sample_map("PXD050470", "human hippocampus", c("H072", "H081", "H187"))
add_records(read_hippocampus_by_sample(path, "PXD050470", "human hippocampus", map))
add_source("PXD050470", "human hippocampus", map$SampleID, map$SampleClass, "Author supplementary Table S3 intensity columns", path)

path <- require_file(file.path(source_root, "PXD063047/search_results/extracted/combined/txt/La (K)Sites.txt"))
map <- sample_map("PXD063047", "normal pregnancy placenta", paste0("Con_", 1:3))
add_records(read_maxquant_by_localization(path, "PXD063047", "normal pregnancy placenta", map, "Proteins"))
add_source("PXD063047", "normal pregnancy placenta", map$SampleID, map$SampleClass, "MaxQuant site table / sample localization", path)

path <- require_file(file.path(source_root, "PXD064912/supplementary/europepmc/mmc1.xlsx"))
map <- sample_map("PXD064912", "human sperm", c("Rep 1", "Rep 2", "Rep 3"))
add_records(read_pxd064912_by_sample(path, "PXD064912", "human sperm", map))
add_source("PXD064912", "human sperm", map$SampleID, map$SampleClass, "Author supplementary PTM table / sample probability", path)

path <- require_file(file.path(source_root, "PXD066054/search_results/extracted/PLa/XB08700B1DPLa_-PTMSiteReport.tsv"))
for (group in c("BPH", "prostate cancer")) {
  ids <- if (group == "BPH") paste0("NAT", 1:5) else paste0("PCa", 1:5)
  map <- sample_map("PXD066054", group, ids)
  add_records(read_spectronaut_kla_by_sample(path, "PXD066054", group, map))
  add_source("PXD066054", group, ids, map$SampleClass, "Spectronaut PTM site report / condition", path)
}

path <- require_file(file.path(source_root, "PXD066351/search_results/XB01472B1DPLa-MSstats_Input.csv"))
map <- sample_map("PXD066351", "HCT116 control and Roseburia co-culture", c("NC116", "R116"), c("control", "Roseburia co-culture"))
add_records(read_pxd066351_by_sample(path, "PXD066351", "HCT116 control and Roseburia co-culture", map))
add_source("PXD066351", "HCT116 control and Roseburia co-culture", map$SampleID, map$SampleClass, "Spectronaut PTM table / condition", path)

path <- require_file(file.path(source_root, "PXD050147/search_results/Lactyl_K_Sites.txt"))
source_ids <- c(
  "WT_Kla_rep1", "WT_Kla_rep2", "WT_Kla_rep3",
  "SIRT1KO_Kla_rep1", "SIRT1KO_Kla_rep2", "SIRT1KO_Kla_rep3",
  "SIRT3KO_Kla_rep1", "SIRT3KO_Kla_rep2", "SIRT3KO_Kla_rep3"
)
display_ids <- c("WT_rep1", "WT_rep2", "WT_rep3", "SIRT1KO_rep1", "SIRT1KO_rep2", "SIRT1KO_rep3", "SIRT3KO_rep1", "SIRT3KO_rep2", "SIRT3KO_rep3")
map <- sample_map("PXD050147", "HepG2 WT and SIRT1 or SIRT3 KO", source_ids, c("WT", "WT", "WT", "SIRT1 KO", "SIRT1 KO", "SIRT1 KO", "SIRT3 KO", "SIRT3 KO", "SIRT3 KO"))
parsed <- read_maxquant_by_localization(path, "PXD050147", "HepG2 WT and SIRT1 or SIRT3 KO", map, "Proteins")
parsed[, SampleID := display_ids[match(SampleID, source_ids)]]
add_records(parsed)
add_source("PXD050147", "HepG2 WT and SIRT1 or SIRT3 KO", display_ids, map$SampleClass, "MaxQuant site table / sample localization", path)

path <- require_file(file.path(source_root, "PXD054919/supplementary/41419_2025_8113_MOESM2_ESM.xlsx"))
map <- sample_map("PXD054919", "A549", paste0("A549-", 1:3))
add_records(read_author_intensity_by_sample(path, "PXD054919", "A549", map))
add_source("PXD054919", "A549", map$SampleID, map$SampleClass, "Author supplementary lactylation table / intensity", path)

path <- require_file(file.path(source_root, "PXD060185/search_results/RESULT/combined/txt/La (K)Sites.txt"))
for (item in list(
  list(group = "MCF7", id = "A"),
  list(group = "MDA-MB-468", id = "B"),
  list(group = "T-47D", id = "D"),
  list(group = "MCF10A", id = "C")
)) {
  map <- sample_map("PXD060185", item$group, item$id, "single")
  add_records(read_maxquant_by_identification(path, "PXD060185", item$group, map, "Protein"))
  add_source("PXD060185", item$group, item$id, "single", "MaxQuant site table / sample identification", path)
}

path <- require_file(file.path(source_root, "PXD063266/search_results/LactylSites.xlsx"))
map <- sample_map("PXD063266", "PC-3M", as.character(1:3))
add_records(read_pxd063266_by_sample(path, "PXD063266", "PC-3M", map))
add_source("PXD063266", "PC-3M", map$SampleID, map$SampleClass, "Author MaxQuant site workbook / sample localization", path)

path <- require_file(file.path(source_root, "PXD070007/search_results/SA206LPLaB1_Annotation.xlsx"))
for (item in list(
  list(group = "glioblastoma stem cells", ids = c("G2907", "G3028", "G3264", "GSC23", "MES28", "RKI")),
  list(group = "neural stem cells", ids = c("ENSA", "HMP1"))
)) {
  map <- sample_map("PXD070007", item$group, item$ids, "model")
  add_records(read_pxd070007_by_model(path, "PXD070007", item$group, map))
  add_source("PXD070007", item$group, map$SampleID, map$SampleClass, "Author annotation workbook / model intensity", path)
}

# MaxQuant evidence plus proteinGroups for the two RKO conditions.
evidence_path <- require_file(file.path(source_root, "PXD078013/search_results/evidence.txt"))
protein_path <- require_file(file.path(source_root, "PXD078013/search_results/proteinGroups.txt"))
evidence <- fread(evidence_path, check.names = FALSE, showProgress = FALSE)
proteins <- fread(protein_path, check.names = FALSE, showProgress = FALSE)
evidence <- evidence[!is_true_flag(evidence$Reverse) & !is_true_flag(evidence$`Potential contaminant`) & positive_numeric(evidence$`La (K)`) & nzchar(evidence$`La (K) site IDs`), , drop = FALSE]
evidence_sites <- lapply(seq_len(nrow(evidence)), function(index) {
  ids <- trimws(unlist(strsplit(as.character(evidence$`La (K) site IDs`[[index]]), ";")))
  data.table(SiteID = ids, SampleID = evidence$Experiment[[index]])
})
evidence_sites <- rbindlist(evidence_sites, fill = TRUE)
proteins <- proteins[!is_true_flag(proteins$Reverse) & !is_true_flag(proteins$`Potential contaminant`) & !is_true_flag(proteins$`Only identified by site`), , drop = FALSE]
protein_sites <- lapply(seq_len(nrow(proteins)), function(index) {
  ids <- trimws(unlist(strsplit(as.character(proteins$`La (K) site IDs`[[index]]), ";")))
  accession <- first_accession_current(ifelse(nzchar(proteins$`Majority protein IDs`[[index]]), proteins$`Majority protein IDs`[[index]], proteins$`Protein IDs`[[index]]))
  if (!nzchar(accession) || !length(ids)) return(data.table())
  data.table(BaseAccession = accession, SiteID = ids)
})
protein_sites <- rbindlist(protein_sites, fill = TRUE)
protein_sites <- unique(merge(protein_sites, evidence_sites, by = "SiteID", allow.cartesian = TRUE)[, .(BaseAccession, SampleID)])
protein_sites[, PXD := "PXD078013"]
protein_sites[, SampleGroup := "RKO WT and GSK3B KO"]
protein_sites[, SourceMode := "deposited_sample_table"]
protein_sites[, SourceFile := paste(source_rel(evidence_path), source_rel(protein_path), sep = ";")]
protein_sites[, SampleClass := ifelse(grepl("sg737", SampleID), "GSK3B KO", "GSK3B WT")]
protein_sites[, ObservationType := "condition"]
setcolorder(protein_sites, c("PXD", "SampleGroup", "SampleID", "BaseAccession", "SourceMode", "SourceFile", "SampleClass", "ObservationType"))
add_records(protein_sites)
add_source("PXD078013", "RKO WT and GSK3B KO", c("sg_1", "sg_2", "sg737_1", "sg737_2"), c("GSK3B WT", "GSK3B WT", "GSK3B KO", "GSK3B KO"), "MaxQuant evidence plus proteinGroups / condition", evidence_path)

path <- require_file(file.path(source_root, "PXD078736/search_results/txt/La(K)Sites.txt"))
map <- sample_map("PXD078736", "HK-2 control and mannitol", c(paste0("ctr_", 1:3), paste0("man_", 1:3)), c(rep("control", 3), rep("mannitol", 3)))
add_records(read_maxquant_by_identification(path, "PXD078736", "HK-2 control and mannitol", map, "Protein", first_accession_current))
add_source("PXD078736", "HK-2 control and mannitol", map$SampleID, map$SampleClass, "MaxQuant site table / sample identification", path)

path <- require_file(file.path(source_root, "PXD028737/search_results/extracted_pairing/combined/combined/txt/La(K)Sites.txt"))
map <- sample_map("PXD028737", "HMC3", c("H0", "H24"), c("normoxia", "hypoxia"))
add_records(read_maxquant_by_localization(path, "PXD028737", "HMC3", map, "Proteins"))
add_source("PXD028737", "HMC3", map$SampleID, map$SampleClass, "MaxQuant site table / sample localization", path)

path <- require_file(file.path(source_root, "PXD073311/search_results/extracted_pairing/IPX0015307003_Database_search_result/Database_search_result/report.pr_matrix.tsv"))
map <- sample_map("PXD073311", "HUVEC control and Pg infection", c("A0h_1", "A0h_2", "A0h_3", "A6h_1", "A6h_2", "A6h_3"), c(rep("A0h control", 3), rep("A6h Pg infection", 3)))
add_records(read_pxd073311_by_raw_file(path, "PXD073311", "HUVEC control and Pg infection", map))
add_source("PXD073311", "HUVEC control and Pg infection", map$SampleID, map$SampleClass, "Spectronaut peptide matrix / positive raw-file intensity", path)

# Dataset-level observations: the source contains technical/structural runs or
# pooled channels, but it does not support independent biological sample IDs.
single_groups <- list(
  list(pxd = "PXD033146", group = "pathological rotator cuff tendon", id = "TMT_plex", class = "plex_channel", file = "PXD033146/search_results/extracted_pairing/search_result-HA119TPLa/La (K)Sites.txt"),
  list(pxd = "PXD075377", group = "adjacent liver", id = "Control_pool", class = "pool", file = "PXD075377/search_results/extracted/2-Basic_analysis/MS_identified_information.txt", intensity = "Intensity Control"),
  list(pxd = "PXD075377", group = "HCC", id = "HCC_pool", class = "pool", file = "PXD075377/search_results/extracted/2-Basic_analysis/MS_identified_information.txt", intensity = "Intensity HCC"),
  list(pxd = "PXD028488", group = "HCT116", id = "HCT116_dataset_union", class = "single"),
  list(pxd = "PXD053474", group = "HCT116", id = "HCT116_dataset_union", class = "single"),
  list(pxd = "PXD028488", group = "TALL-104", id = "TALL-104_dataset_union", class = "single"),
  list(pxd = "PXD028488", group = "HEK293T", id = "HEK293T_dataset_union", class = "single"),
  list(pxd = "PXD058534", group = "pretreated HK-2", id = "HK2", class = "single", file = "PXD058534/search_results/extracted_pairing/txt/La (K)Sites.txt")
)
for (item in single_groups) {
  if (item$pxd == "PXD075377") {
    path <- require_file(file.path(source_root, item$file))
    add_records(read_itraq_pool(path, item$pxd, item$group, item$id, item$class, item$intensity))
    add_source(item$pxd, item$group, item$id, item$class, "iTRAQ source table / pooled channel", path)
  } else if (item$pxd == "PXD058534") {
    path <- require_file(file.path(source_root, item$file))
    add_records(read_maxquant_single_group(path, item$pxd, item$group, item$id, item$class))
    add_source(item$pxd, item$group, item$id, item$class, "MaxQuant site table / single sample", path)
  } else {
    source_path <- file.path(input_dir, "kla_protein_membership_30.csv")
    accessions <- membership[PXD == item$pxd & SampleGroup == item$group, BaseAccession]
    add_records(records(item$pxd, item$group, item$id, accessions,
      "validated_publication_membership", source_path, item$class, "dataset_union"))
    add_source(item$pxd, item$group, item$id, item$class, "Validated frozen group union / structural runs collapsed", source_path)
  }
}

all_records <- unique(rbindlist(sample_records, fill = TRUE))
all_records[, GroupKey := paste(PXD, SampleGroup, sep = "__")]
all_records[, SourceFile := vapply(SourceFile, source_rel, character(1))]
stop_if(setequal(unique(all_records$GroupKey), groups$GroupKey), paste0(
  "Sample preparation did not cover exactly the frozen groups. Missing: ",
  paste(setdiff(groups$GroupKey, unique(all_records$GroupKey)), collapse = ", ")
))
stop_if(!anyDuplicated(all_records[, .(PXD, SampleGroup, SampleID, BaseAccession)]), "Duplicate sample-level membership rows detected.")

all_records <- merge(
  all_records,
  unique(membership[, .(PXD, SampleGroup, BaseAccession, IsDdr)]),
  by = c("PXD", "SampleGroup", "BaseAccession"),
  all.x = TRUE,
  sort = FALSE
)
all_records[is.na(IsDdr), IsDdr := FALSE]

sample_values <- all_records[, .(
  KlaProteinCount = uniqueN(BaseAccession),
  KlaDdrProteinCount = uniqueN(BaseAccession[IsDdr == TRUE]),
  SourceMode = first(SourceMode),
  SourceFile = first(SourceFile),
  SampleClass = first(SampleClass),
  ObservationType = first(ObservationType)
), by = .(PXD, SampleGroup, SampleID)]
sample_values[, KlaDdrFraction := fifelse(KlaProteinCount > 0L, KlaDdrProteinCount / KlaProteinCount, NA_real_)]
sample_values[, KlaDdrFractionPercentage := KlaDdrFraction * 100]

sample_values <- merge(
  sample_values,
  groups[, .(PXD, SampleGroup, RowOrder, Category, DisplayGroup = KlaLabelEn, ReferenceFraction = ReferenceDdrFraction * 100, ReferenceDdr = ReferenceDdrProteinCount, ReferenceTotal = ReferenceProteinCount)],
  by = c("PXD", "SampleGroup"),
  all.x = TRUE,
  sort = FALSE
)
sample_values[, ConditionLabel := fifelse(
  SampleClass %in% c("GSK3B WT", "GSK3B KO", "control", "Roseburia co-culture", "mannitol", "A0h control", "A6h Pg infection"),
  SampleClass,
  SampleID
)]

frozen_group_values <- groups[, .(
  FrozenKlaProteinCount = KlaProteinCount,
  FrozenKlaDdrProteinCount = KlaDdrProteinCount,
  FrozenKlaDdrFraction = KlaDdrFraction * 100
), by = .(PXD, SampleGroup)]
observed_group_values <- all_records[, .(
  ObservedKlaProteinCount = uniqueN(BaseAccession),
  ObservedKlaDdrProteinCount = uniqueN(BaseAccession[IsDdr == TRUE]),
  ObservedSampleCount = uniqueN(SampleID),
  SourceModes = paste(sort(unique(SourceMode)), collapse = ";")
), by = .(PXD, SampleGroup)]
reconciliation <- merge(frozen_group_values, observed_group_values, by = c("PXD", "SampleGroup"), all = TRUE, sort = FALSE)
reconciliation[, KlaProteinCountMatchesFrozen := ObservedKlaProteinCount == FrozenKlaProteinCount]
reconciliation[, KlaDdrProteinCountMatchesFrozen := ObservedKlaDdrProteinCount == FrozenKlaDdrProteinCount]
reconciliation[, GroupUnionStatus := fifelse(
  KlaProteinCountMatchesFrozen & KlaDdrProteinCountMatchesFrozen, "PASS", "SAMPLE_SCOPE_DIFFERENCE"
)]
reconciliation <- merge(
  reconciliation,
  groups[, .(PXD, SampleGroup, RowOrder)],
  by = c("PXD", "SampleGroup"),
  all.x = TRUE,
  sort = FALSE
)

stop_if(!anyNA(sample_values$KlaDdrFraction), "A sample-level fraction is missing.")
stop_if(all(sample_values$KlaDdrFraction >= 0 & sample_values$KlaDdrFraction <= 1), "A sample-level fraction is outside 0-1.")
stop_if(all(reconciliation$ObservedSampleCount >= 1L), "A frozen group has no sample-level observation.")

registry_table <- rbindlist(registry, fill = TRUE)
setorder(registry_table, PXD, SampleGroup, SampleID)
setorder(sample_values, RowOrder, SampleID)
setorder(reconciliation, RowOrder)

sample_values <- sample_values[, .(
  RowOrder, PXD, SampleGroup, Category, DisplayGroup, SampleID, ConditionLabel,
  SampleClass, ObservationType, SourceMode, SourceFile,
  KlaProteinCount, KlaDdrProteinCount, KlaDdrFraction, KlaDdrFractionPercentage,
  ReferenceFraction, ReferenceDdr, ReferenceTotal
)]
reconciliation <- reconciliation[, .(
  RowOrder, PXD, SampleGroup, FrozenKlaProteinCount, FrozenKlaDdrProteinCount, FrozenKlaDdrFraction,
  ObservedKlaProteinCount, ObservedKlaDdrProteinCount, ObservedSampleCount,
  SourceModes, KlaProteinCountMatchesFrozen, KlaDdrProteinCountMatchesFrozen, GroupUnionStatus
)]

fwrite(sample_values, file.path(candidate_dir, "sample_boxplot_values.csv"), na = "")
fwrite(reconciliation, file.path(candidate_dir, "sample_boxplot_reconciliation.csv"), na = "")
fwrite(registry_table, file.path(candidate_dir, "sample_boxplot_source_registry.csv"), na = "")

message(
  "Wrote sample-level boxplot inputs for ", nrow(sample_values),
  " observations across ", uniqueN(sample_values[, .(PXD, SampleGroup)]),
  " publication groups."
)
