#!/usr/bin/env Rscript

# Prepare source-resolved whole-proteome values for three exploratory
# MKI67-normalized Figure 1 plots. The approved publication inputs and
# renderers are not changed.

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
    "Usage: Rscript R/candidate/prepare_figure1_mki67_ratio_inputs.R <project-root> <source-data-root>",
    call. = FALSE
  )
}

stop_if <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

require_file <- function(file_path) {
  stop_if(file.exists(file_path), paste0("Required source file is missing: ", file_path))
  file_path
}

target_labels <- c("MKI67", "ACTB", "TUBB", "H3C1")
target_accessions <- c(MKI67 = "P46013", ACTB = "P60709", TUBB = "P07437", H3C1 = "P68431")
denominator_labels <- c("ACTB", "TUBB", "H3C1")

base_accession <- function(value) {
  value <- trimws(as.character(value))
  value[is.na(value)] <- ""
  value <- sub("^(?:REV__|CON__)+", "", value, perl = TRUE)
  value <- sub("^(sp|tr)[|]", "", value, perl = TRUE)
  value <- sub("[|].*$", "", value, perl = TRUE)
  value <- sub("^NX_", "", value)
  sub("-[0-9]+$", "", value)
}

split_accessions <- function(value) {
  tokens <- unlist(strsplit(as.character(value), "[;,]"), use.names = FALSE)
  tokens <- base_accession(tokens)
  sort(unique(tokens[nzchar(tokens)]))
}

has_accession <- function(value, target) {
  target %in% split_accessions(value)
}

as_numeric_value <- function(value) {
  value <- trimws(as.character(value))
  value[is.na(value) | value %in% c("", "NA", "NaN", "nan", "Not Found")] <- NA_character_
  value <- gsub(",", "", value, fixed = TRUE)
  suppressWarnings(as.numeric(value))
}

read_text_table <- function(file_path, separator = "\t", skip = 0L) {
  if (separator == ",") {
    return(fread(
      file_path,
      sep = separator,
      skip = skip,
      check.names = FALSE,
      colClasses = "character",
      showProgress = FALSE
    ))
  }
  as.data.table(read.delim(
    file_path,
    header = TRUE,
    sep = separator,
    skip = skip,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    quote = "\"",
    comment.char = "",
    colClasses = "character"
  ))
}

choose_accession_column <- function(data) {
  candidates <- c(
    "Majority protein IDs", "Protein IDs", "Accession", "Protein accession",
    "PG.ProteinGroups", "PG.ProteinAccessions", "Protein.Group", "Protein.Group IDs"
  )
  present <- candidates[candidates %in% names(data)]
  if (length(present)) present[[1L]] else NA_character_
}

choose_first_column <- function(data, candidates) {
  present <- candidates[candidates %in% names(data)]
  if (length(present)) present[[1L]] else NA_character_
}

target_row_indices <- function(data, accession_column) {
  stop_if(!is.na(accession_column) && accession_column %in% names(data),
    paste0("No stable accession column was found in source table. Columns: ", paste(names(data), collapse = ", ")))
  values <- as.character(data[[accession_column]])
  lapply(target_accessions, function(target) {
    which(vapply(values, has_accession, logical(1), target = target))
  }) |> setNames(target_labels)
}

make_suffix_column_map <- function(sample_ids, columns, suffix_builder) {
  values <- setNames(rep(NA_character_, length(sample_ids)), sample_ids)
  for (sample_id in sample_ids) {
    suffix <- suffix_builder(sample_id)
    hits <- columns[endsWith(columns, suffix)]
    if (length(hits) == 1L) values[[sample_id]] <- hits[[1L]]
  }
  values
}

apply_transform <- function(values, transform) {
  if (transform == "identity") return(values)
  if (transform == "back_transform_2^x") {
    return(ifelse(is.finite(values), 2^values, NA_real_))
  }
  stop(paste0("Unknown quantitative transform: ", transform), call. = FALSE)
}

extract_table_values <- function(
  data, accession_column, observations, sample_columns,
  parser, quantitation_field, transform = "identity", source_file
) {
  target_rows <- target_row_indices(data, accession_column)
  output <- lapply(seq_len(nrow(observations)), function(index) {
    observation <- observations[index]
    sample_id <- observation$SampleID[[1L]]
    quant_column <- unname(sample_columns[[sample_id]])
    source_values <- setNames(as.list(rep(NA_real_, length(target_labels))), target_labels)
    row_counts <- setNames(integer(length(target_labels)), target_labels)
    reported_groups <- setNames(as.list(rep(NA_character_, length(target_labels))), target_labels)

    if (length(quant_column) == 1L && !is.na(quant_column) && nzchar(quant_column)) {
      for (label in target_labels) {
        indices <- target_rows[[label]]
        row_counts[[label]] <- length(indices)
        if (length(indices)) {
          reported_groups[[label]] <- paste(unique(as.character(data[[accession_column]][indices])), collapse = ";")
          numeric_values <- as_numeric_value(data[[quant_column]][indices])
          numeric_values <- numeric_values[is.finite(numeric_values)]
          if (length(numeric_values)) source_values[[label]] <- sum(numeric_values)
        }
      }
    }

    result <- data.table(
      ObsKey = observation$ObsKey[[1L]],
      SourceFile = source_file,
      Parser = parser,
      QuantitationField = quantitation_field,
      QuantitationTransform = transform,
      QuantitationColumn = if (length(quant_column) == 1L && !is.na(quant_column)) quant_column else NA_character_,
      SourceColumnStatus = if (length(quant_column) == 1L && !is.na(quant_column)) "sample_column_found" else "sample_column_missing"
    )
    for (label in target_labels) {
      result[[paste0(label, "_SourceValue")]] <- source_values[[label]]
      result[[paste0(label, "_Intensity")]] <- apply_transform(source_values[[label]], transform)
      result[[paste0(label, "_RowCount")]] <- row_counts[[label]]
      result[[paste0(label, "_ReportedGroup")]] <- reported_groups[[label]]
    }
    result
  })
  rbindlist(output, fill = TRUE)
}

extract_matrix_values <- function(data, observations, matrix_cell_lines, source_file) {
  target_columns <- setNames(
    c("P46013;KI67_HUMAN", "P60709;ACTB_HUMAN", "P07437;TBB5_HUMAN", "P68431;H3C1_HUMAN"),
    target_labels
  )
  output <- lapply(seq_len(nrow(observations)), function(index) {
    observation <- observations[index]
    cell_line <- unname(matrix_cell_lines[[observation$SampleGroup[[1L]]]])
    row <- data[data$Project_Identifier == cell_line, , drop = FALSE]
    result <- data.table(
      ObsKey = observation$ObsKey[[1L]],
      SourceFile = source_file,
      Parser = "ProCan averaged protein matrix",
      QuantitationField = "matrix protein abundance",
      QuantitationTransform = "back_transform_2^x",
      QuantitationColumn = "Project_Identifier",
      SourceColumnStatus = ifelse(nrow(row) == 1L, "sample_profile_found", "sample_profile_missing")
    )
    for (label in target_labels) {
      column <- target_columns[[label]]
      raw_value <- if (length(column) == 1L && column %in% names(row) && nrow(row) == 1L) {
        as_numeric_value(row[[column]][[1L]])
      } else {
        NA_real_
      }
      result[[paste0(label, "_SourceValue")]] <- raw_value
      result[[paste0(label, "_Intensity")]] <- apply_transform(raw_value, "back_transform_2^x")
      result[[paste0(label, "_RowCount")]] <- if (length(column) == 1L && column %in% names(row)) nrow(row) else 0L
      result[[paste0(label, "_ReportedGroup")]] <- if (length(column) == 1L && column %in% names(row)) column else NA_character_
    }
    result
  })
  rbindlist(output, fill = TRUE)
}

unavailable_values <- function(observations, source_file, parser, quantitation_field, reason) {
  output <- copy(observations[, .(ObsKey)])
  output[, c(
    "SourceFile", "Parser", "QuantitationField", "QuantitationTransform",
    "QuantitationColumn", "SourceColumnStatus"
  ) := list(
    source_file, parser, quantitation_field, "not_available", NA_character_, reason
  )]
  for (label in target_labels) {
    output[, (paste0(label, "_SourceValue")) := NA_real_]
    output[, (paste0(label, "_Intensity")) := NA_real_]
    output[, (paste0(label, "_RowCount")) := 0L]
    output[, (paste0(label, "_ReportedGroup")) := NA_character_]
  }
  output
}

candidate_dir <- file.path(project_root, "data", "candidate")
dir.create(candidate_dir, recursive = TRUE, showWarnings = FALSE)

sample_values <- fread(file.path(candidate_dir, "figure1_sample_boxplot_values.csv"), check.names = FALSE)
registry <- fread(file.path(candidate_dir, "figure1_sample_boxplot_source_registry.csv"), check.names = FALSE)
whole_values <- sample_values[Dataset == "Whole proteome"]
whole_registry <- registry[Dataset == "Whole proteome"]
stop_if(nrow(whole_values) == 118L, "Whole-proteome Figure 1 source input must contain 118 observations.")
stop_if(uniqueN(whole_values[, .(PXD, SampleGroup, SampleID)]) == 118L,
  "Whole-proteome Figure 1 observations are not unique.")

observation_columns <- c(
  "RowOrder", "PXD", "SampleGroup", "Category", "DisplayGroup", "Dataset", "SampleID",
  "ConditionLabel", "SampleClass", "ObservationType", "ReferencePXD", "SourceFile"
)
observations <- merge(
  whole_values[, ..observation_columns],
  whole_registry[, .(PXD, SampleGroup, SampleID, RegistrySourceFile = SourceFile)],
  by = c("PXD", "SampleGroup", "SampleID"),
  all.x = TRUE,
  sort = FALSE
)
observations[, ObsKey := paste(PXD, SampleGroup, SampleID, sep = "__")]
observations[, SourceFile := fifelse(
  nzchar(trimws(as.character(RegistrySourceFile))), RegistrySourceFile, SourceFile
)]
stop_if(!anyNA(observations$RegistrySourceFile), "A whole-proteome observation is missing from the source registry.")
stop_if(!anyDuplicated(observations$ObsKey), "Duplicate whole-proteome observation keys detected.")

parsed_values <- list()
add_parsed <- function(result) {
  parsed_values[[length(parsed_values) + 1L]] <<- result
}

add_source_table <- function(source_file, observation_filter, data, accession_column, sample_columns,
                             parser, quantitation_field, transform = "identity") {
  selected <- observations[eval(substitute(observation_filter), envir = parent.frame())]
  if (!nrow(selected)) return(invisible(NULL))
  add_parsed(extract_table_values(
    data, accession_column, selected, sample_columns, parser,
    quantitation_field, transform, source_file
  ))
  invisible(NULL)
}

source_file <- "PXD028488/search_results/Nonenrichment-Search files/TALL-Nonenrichment-Search files/proteins.csv"
if (any(observations$SourceFile == source_file)) {
  file_path <- require_file(file.path(source_root, source_file))
  data <- read_text_table(file_path, separator = ",")
  columns <- setNames(paste0("Area Sample ", 1:3), paste0("Sample ", 1:3))
  add_source_table(source_file, SourceFile == source_file, data, choose_accession_column(data), columns,
    "TALL protein table", "Area", "identity")
}

source_file <- "PXD028737/search_results/extracted_reference/txt/proteinGroups.txt"
if (any(observations$SourceFile == source_file)) {
  file_path <- require_file(file.path(source_root, source_file))
  data <- read_text_table(file_path)
  columns <- setNames(paste0("Intensity ", c("H0", "H24")), c("H0", "H24"))
  add_source_table(source_file, SourceFile == source_file, data, choose_accession_column(data), columns,
    "MaxQuant proteinGroups", "Intensity", "identity")
}

source_file <- "PXD033146/search_results/extracted_pairing/search_result-HA119TQ/proteinGroups.txt"
if (any(observations$SourceFile == source_file)) {
  file_path <- require_file(file.path(source_root, source_file))
  data <- read_text_table(file_path)
  quant_columns <- grep("^Reporter intensity corrected ", names(data), value = TRUE)
  channel_map <- c(RCT1 = "1", RCT2 = "2", RCT3 = "3")
  columns <- setNames(
    quant_columns[match(channel_map, sub("^Reporter intensity corrected ", "", quant_columns))],
    names(channel_map)
  )
  add_source_table(source_file, SourceFile == source_file, data, choose_accession_column(data), columns,
    "MaxQuant reporter proteinGroups", "Reporter intensity corrected", "identity")
}

source_file <- "PXD010154/search_results/extracted_healthy_tissues/lung/P013163_lung_proteinGroups.txt"
if (any(observations$SourceFile == source_file)) {
  file_path <- require_file(file.path(source_root, source_file))
  data <- read_text_table(file_path)
  quant_column <- choose_first_column(data, c("Intensity", "Intensity P013163"))
  selected <- observations[SourceFile == source_file]
  columns <- setNames(rep(quant_column, nrow(selected)), selected$SampleID)
  add_source_table(source_file, SourceFile == source_file, data, choose_accession_column(data), columns,
    "MaxQuant aggregate proteinGroups", quant_column, "identity")
}

source_file <- "PXD010154/search_results/extracted_healthy_tissues/placenta/P013680_placenta_proteinGroups.txt"
if (any(observations$SourceFile == source_file)) {
  file_path <- require_file(file.path(source_root, source_file))
  data <- read_text_table(file_path)
  quant_column <- choose_first_column(data, c("Intensity", "Intensity P013680"))
  selected <- observations[SourceFile == source_file]
  columns <- setNames(rep(quant_column, nrow(selected)), selected$SampleID)
  add_source_table(source_file, SourceFile == source_file, data, sample_columns = columns,
    accession_column = choose_accession_column(data), parser = "MaxQuant aggregate proteinGroups",
    quantitation_field = quant_column, transform = "identity")
}

source_file <- "PXD046800/search_results/HFX2_LFQ_QB002_Proteins.txt"
if (any(observations$SourceFile == source_file)) {
  file_path <- require_file(file.path(source_root, source_file))
  data <- read_text_table(file_path)
  abundance_columns <- grep("^Abundance:", names(data), value = TRUE)
  columns <- make_suffix_column_map(
    observations[SourceFile == source_file, SampleID],
    abundance_columns,
    function(sample_id) paste0(", ", sample_id)
  )
  add_source_table(source_file, SourceFile == source_file, data, choose_accession_column(data), columns,
    "Proteome Discoverer protein table", "Abundance:", "identity")
}

source_file <- "PXD050147/search_results/SIRT_proteinGroups.txt"
if (any(observations$SourceFile == source_file)) {
  file_path <- require_file(file.path(source_root, source_file))
  data <- read_text_table(file_path)
  selected <- observations[SourceFile == source_file]
  columns <- setNames(paste0("Intensity ", selected$SampleID), selected$SampleID)
  add_source_table(source_file, SourceFile == source_file, data, choose_accession_column(data), columns,
    "MaxQuant proteinGroups", "Intensity", "identity")
}

source_file <- "PXD050470/supplementary/prca2331-sup-0006-tables4.xlsx"
if (any(observations$SourceFile == source_file)) {
  file_path <- require_file(file.path(source_root, source_file))
  data <- as.data.table(read_excel(file_path, sheet = "Sheet1", skip = 5L, col_types = "text"))
  selected <- observations[SourceFile == source_file]
  columns <- setNames(paste0("Intensity_", selected$SampleID), selected$SampleID)
  add_source_table(source_file, SourceFile == source_file, data, choose_accession_column(data), columns,
    "Author supplementary protein table", "Intensity_H072/H081/H187", "identity")
}

source_file <- "PXD002400/search_results/msms.zip"
if (any(observations$SourceFile == source_file)) {
  selected <- observations[SourceFile == source_file]
  add_parsed(unavailable_values(
    selected, source_file, "PXD002400 msms archive", "protein intensity",
    "no processed protein intensity table in archive"
  ))
}

source_file <- "PXD022005/search_results/txt_proteomics.zip"
if (any(observations$SourceFile == source_file)) {
  file_path <- require_file(file.path(source_root, source_file))
  temporary_dir <- tempfile("pxd022005_proteins_")
  dir.create(temporary_dir)
  unzip(file_path, files = "proteinGroups.txt", exdir = temporary_dir)
  protein_path <- require_file(file.path(temporary_dir, "proteinGroups.txt"))
  data <- read_text_table(protein_path)
  selected <- observations[SourceFile == source_file]
  columns <- setNames(rep("Intensity", nrow(selected)), selected$SampleID)
  add_source_table(source_file, SourceFile == source_file, data, choose_accession_column(data), columns,
    "MaxQuant proteinGroups from archive", "Intensity", "identity")
}

source_file <- "PXD030304/search_results/ProCan-DepMapSanger_protein_matrix_6692_averaged.txt"
if (any(observations$SourceFile == source_file)) {
  file_path <- require_file(file.path(source_root, source_file))
  target_columns <- c(
    "Project_Identifier", "P46013;KI67_HUMAN", "P60709;ACTB_HUMAN",
    "P07437;TBB5_HUMAN", "P68431;H3C1_HUMAN"
  )
  header <- names(fread(file_path, nrows = 0L, check.names = FALSE, showProgress = FALSE))
  available_columns <- target_columns[target_columns %in% header]
  data <- fread(file_path, select = available_columns, check.names = FALSE, showProgress = FALSE)
  matrix_cell_lines <- c(
    HCT116 = "SIDM00783;HCT-116",
    MCF7 = "SIDM00148;MCF7",
    "MDA-MB-468" = "SIDM00628;MDA-MB-468",
    "T-47D" = "SIDM00097;T47D",
    A549 = "SIDM00903;A549",
    "RKO WT and GSK3B KO" = "SIDM01090;RKO"
  )
  selected <- observations[SourceFile == source_file]
  add_parsed(extract_matrix_values(data, selected, matrix_cell_lines, source_file))
}

source_file <- "PXD066517/search_results/20240275.tsv"
if (any(observations$SourceFile == source_file)) {
  file_path <- require_file(file.path(source_root, source_file))
  data <- read_text_table(file_path)
  columns <- make_suffix_column_map(
    observations[SourceFile == source_file, SampleID],
    names(data),
    function(sample_id) paste0("_", sample_id, ".raw.PG.Quantity")
  )
  add_source_table(source_file, SourceFile == source_file, data, choose_accession_column(data), columns,
    "Spectronaut protein quantity table", ".PG.Quantity", "identity")
}

for (source_file in c(
  "PXD066054/search_results/extracted/DA/Protein_Quant.tsv",
  "PXD066351/search_results/XB01472B1DA-Protein_Quant.tsv"
)) {
  if (!any(observations$SourceFile == source_file)) next
  file_path <- require_file(file.path(source_root, source_file))
  data <- read_text_table(file_path)
  selected <- observations[SourceFile == source_file]
  quantity_columns <- names(data)[endsWith(names(data), ".PG.Quantity")]
  columns <- make_suffix_column_map(
    selected$SampleID,
    quantity_columns,
    function(sample_id) paste0("_", sample_id, "_")
  )
  for (sample_id in selected$SampleID) {
    if (!is.na(columns[[sample_id]])) next
    hits <- quantity_columns[grepl(paste0("_", sample_id, "_"), quantity_columns, fixed = TRUE)]
    if (length(hits) == 1L) columns[[sample_id]] <- hits[[1L]]
  }
  add_source_table(source_file, SourceFile == source_file, data, choose_accession_column(data), columns,
    "Spectronaut protein quantity table", ".PG.Quantity", "identity")
}

source_file <- "PXD065775/search_results/20170330_01-24_patients_iTRAQ.xlsx"
if (any(observations$SourceFile == source_file)) {
  file_path <- require_file(file.path(source_root, source_file))
  selected <- observations[SourceFile == source_file]
  outputs <- lapply(unique(selected$SampleGroup), function(group) {
    group_observations <- selected[SampleGroup == group]
    sheet <- if (group == "HCC") "CISs" else "ANTs"
    data <- as.data.table(read_excel(file_path, sheet = sheet, col_types = "text"))
    columns <- setNames(rep(NA_character_, nrow(group_observations)), group_observations$SampleID)
    for (sample_id in group_observations$SampleID) {
      if (sample_id %in% names(data)) columns[[sample_id]] <- sample_id
    }
    extract_table_values(
      data, choose_accession_column(data), group_observations, columns,
      "iTRAQ protein table", "sample iTRAQ ratio", "identity", source_file
    )
  })
  add_parsed(rbindlist(outputs, fill = TRUE))
}

source_file <- "PXD069969/search_results/SA206LQB1_Annotation.xlsx"
if (any(observations$SourceFile == source_file)) {
  file_path <- require_file(file.path(source_root, source_file))
  data <- as.data.table(read_excel(file_path, sheet = "Annotation_Combine", col_types = "text"))
  selected <- observations[SourceFile == source_file]
  columns <- setNames(paste0("LFQ intensity ", selected$SampleID), selected$SampleID)
  add_source_table(source_file, SourceFile == source_file, data, choose_accession_column(data), columns,
    "Author annotation protein table", "LFQ intensity", "identity")
}

source_file <- "PXD073311/search_results/extracted_pairing/IPX0015307001_Database_search_result/Database_search_result/report.pg_matrix.tsv"
if (any(observations$SourceFile == source_file)) {
  file_path <- require_file(file.path(source_root, source_file))
  data <- read_text_table(file_path)
  selected <- observations[SourceFile == source_file]
  columns <- make_suffix_column_map(
    selected$SampleID,
    names(data),
    function(sample_id) paste0("_", sample_id, ".raw")
  )
  add_source_table(source_file, SourceFile == source_file, data, choose_accession_column(data), columns,
    "Spectronaut protein matrix", "positive raw-file intensity", "identity")
}

source_file <- "PXD072220/search_results/HK-2_Spectronaut-report_PG_Quantity.txt"
if (any(observations$SourceFile == source_file)) {
  file_path <- require_file(file.path(source_root, source_file))
  data <- read_text_table(file_path, skip = 2L)
  selected <- observations[SourceFile == source_file]
  columns <- make_suffix_column_map(
    selected$SampleID,
    names(data),
    function(sample_id) paste0("_", sample_id, ".raw.PG.MS2Quantity")
  )
  add_source_table(source_file, SourceFile == source_file, data, choose_accession_column(data), columns,
    "Spectronaut protein-group quantity table", "PG.MS2Quantity", "identity")
}

source_file <- "PXD030304/search_results/ProCan-DepMapSanger_peptide_counts_per_protein_per_sample.txt"
if (any(observations$SourceFile == source_file)) {
  selected <- observations[SourceFile == source_file]
  add_parsed(unavailable_values(
    selected, source_file, "ProCan peptide-count matrix",
    "peptide count, not protein expression intensity",
    "source contains peptide counts only"
  ))
}

parsed <- rbindlist(parsed_values, fill = TRUE)
stop_if(uniqueN(parsed$ObsKey) == nrow(observations),
  paste0("Source parsing did not return exactly one record per observation. Got ",
    uniqueN(parsed$ObsKey), " of ", nrow(observations), "."))
stop_if(!anyDuplicated(parsed$ObsKey), "Source parsing returned duplicate observation records.")

parsed[, SourceFile := NULL]
output <- merge(observations, parsed, by = "ObsKey", all.x = TRUE, sort = FALSE)
stop_if(!anyNA(output$PXD), "A parsed ratio record could not be matched to Figure 1 metadata.")

for (label in target_labels) {
  output[, (paste0(label, "_Intensity")) := as.numeric(get(paste0(label, "_Intensity")))]
}

for (denominator in denominator_labels) {
  ratio_column <- paste0("MKI67_over_", denominator)
  intensity_column <- paste0(denominator, "_Intensity")
  output[, (ratio_column) := fifelse(
    is.finite(MKI67_Intensity) & MKI67_Intensity > 0 &
      is.finite(get(intensity_column)) & get(intensity_column) > 0,
    MKI67_Intensity / get(intensity_column),
    NA_real_
  )]
  output[, (paste0("MKI67_over_", denominator, "_Status")) := fifelse(
    !is.finite(MKI67_Intensity) | MKI67_Intensity <= 0,
    "MKI67 missing or non-positive",
    fifelse(
      !is.finite(get(intensity_column)) | get(intensity_column) <= 0,
      paste0(denominator, " missing or non-positive"),
      "valid"
    )
  )]
}

target_key <- data.table(
  ProteinLabel = target_labels,
  BaseAccession = unname(target_accessions),
  Role = c("numerator", "denominator", "denominator", "denominator"),
  DisplayName = c("Ki-67", "beta-actin", "beta-tubulin", "histone H3.1")
)
fwrite(target_key, file.path(candidate_dir, "figure1_mki67_ratio_protein_key.csv"))

audit_columns <- c(
  "RowOrder", "PXD", "SampleGroup", "Category", "DisplayGroup", "Dataset", "SampleID",
  "ConditionLabel", "SampleClass", "ObservationType", "ReferencePXD", "ObsKey",
  "SourceFile", "Parser", "QuantitationField", "QuantitationTransform", "QuantitationColumn",
  "SourceColumnStatus"
)
audit_columns <- c(audit_columns, unlist(lapply(target_labels, function(label) c(
  paste0(label, "_SourceValue"), paste0(label, "_Intensity"),
  paste0(label, "_RowCount"), paste0(label, "_ReportedGroup")
))))
audit_columns <- c(audit_columns, unlist(lapply(denominator_labels, function(label) c(
  paste0("MKI67_over_", label), paste0("MKI67_over_", label, "_Status")
))))
audit <- output[, ..audit_columns]
setorder(audit, RowOrder, PXD, SampleGroup, SampleID)
fwrite(audit, file.path(candidate_dir, "figure1_mki67_ratio_source_audit.csv"))

ratio_long <- melt(
  output,
  id.vars = c(
    "RowOrder", "PXD", "SampleGroup", "Category", "DisplayGroup", "Dataset", "SampleID",
    "ConditionLabel", "SampleClass", "ObservationType", "ReferencePXD", "ObsKey", "SourceFile",
    "Parser", "QuantitationField", "QuantitationTransform", "QuantitationColumn", "SourceColumnStatus"
  ),
  measure.vars = paste0("MKI67_over_", denominator_labels),
  variable.name = "RatioName",
  value.name = "Ratio"
)
ratio_long[, Denominator := sub("^MKI67_over_", "", RatioName)]
ratio_long[, RatioLabel := paste0("MKI67 / ", Denominator)]
ratio_long <- ratio_long[is.finite(Ratio) & Ratio > 0]
ratio_long[, RatioName := NULL]
setcolorder(ratio_long, c(
  "RowOrder", "PXD", "SampleGroup", "Category", "DisplayGroup", "Dataset", "SampleID",
  "ConditionLabel", "SampleClass", "ObservationType", "ReferencePXD", "ObsKey", "SourceFile",
  "Parser", "QuantitationField", "QuantitationTransform", "QuantitationColumn", "SourceColumnStatus",
  "Denominator", "RatioLabel", "Ratio"
))
setorder(ratio_long, Denominator, RowOrder, PXD, SampleGroup, SampleID)
fwrite(ratio_long, file.path(candidate_dir, "figure1_mki67_ratio_sample_values.csv"))

coverage <- ratio_long[, .(
  N = .N,
  Min = min(Ratio),
  Q1 = quantile(Ratio, 0.25, names = FALSE),
  Median = median(Ratio),
  Q3 = quantile(Ratio, 0.75, names = FALSE),
  Max = max(Ratio)
), by = .(Denominator, Category)]
all_coverage <- CJ(Denominator = denominator_labels, Category = unique(observations$Category), unique = TRUE)
coverage <- merge(all_coverage, coverage, by = c("Denominator", "Category"), all.x = TRUE)
setorder(coverage, Denominator, Category)
fwrite(coverage, file.path(candidate_dir, "figure1_mki67_ratio_coverage.csv"))

message(
  "Prepared MKI67 ratio inputs: ", nrow(ratio_long),
  " valid sample ratios across ", uniqueN(ratio_long$ObsKey),
  " whole-proteome observations."
)
print(ratio_long[, .(N = .N, Observations = uniqueN(ObsKey)), by = Denominator])
print(coverage)
