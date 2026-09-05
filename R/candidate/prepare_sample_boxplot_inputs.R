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

pathway_order <- c("BER", "NER", "MMR", "FA", "HR", "AEJ", "NHEJ")

read_frozen_pathway_scores <- function(path) {
  score_tables <- lapply(excel_sheets(path), function(sheet) {
    data <- as.data.table(read_excel(path, sheet = sheet))
    stop_if(all(c("BaseAccession", pathway_order) %in% names(data)),
      paste0("Missing pathway score columns in ", sheet))
    data[, c("BaseAccession", pathway_order), with = FALSE]
  })
  scores <- rbindlist(score_tables, fill = TRUE)
  scores[, BaseAccession := base_accession(BaseAccession)]
  scores <- scores[nzchar(BaseAccession)]
  for (pathway in pathway_order) {
    conflicts <- scores[, .(DistinctScoreCount = uniqueN(get(pathway))), by = BaseAccession][DistinctScoreCount > 1L]
    stop_if(!nrow(conflicts), paste0("Frozen pathway scores disagree for ", pathway))
  }
  score_matrix <- as.matrix(scores[, ..pathway_order])
  storage.mode(score_matrix) <- "numeric"
  stop_if(all(score_matrix %in% c(-1, 0, 1)), "Frozen pathway scores must be -1, 0 or +1.")
  unique(scores[, lapply(.SD, function(values) values[[1L]]), by = BaseAccession, .SDcols = pathway_order])
}

build_pathway_summary_profile <- function(all_records, groups, pathway_scores) {
  sample_keys <- unique(all_records[, .(PXD, SampleGroup, SampleID)])
  ddr_totals <- all_records[IsDdr == TRUE, .(
    KlaDdrProteinCount = uniqueN(BaseAccession)
  ), by = .(PXD, SampleGroup, SampleID)]
  stop_if(nrow(ddr_totals) == nrow(sample_keys),
    "A source sample is missing a Kla-DDR denominator.")

  sample_ddr <- unique(all_records[
    IsDdr == TRUE & BaseAccession %in% pathway_scores$BaseAccession,
    .(PXD, SampleGroup, SampleID, BaseAccession)
  ])
  joined <- merge(sample_ddr, pathway_scores, by = "BaseAccession", all = FALSE)
  long <- melt(
    joined,
    id.vars = c("PXD", "SampleGroup", "SampleID", "BaseAccession"),
    measure.vars = pathway_order,
    variable.name = "Pathway",
    value.name = "SignedState"
  )
  counts <- long[, .(
    PositiveProteinCount = uniqueN(BaseAccession[SignedState == 1]),
    NegativeProteinCount = uniqueN(BaseAccession[SignedState == -1]),
    AnyPathwayProteinCount = uniqueN(BaseAccession[SignedState != 0])
  ), by = .(PXD, SampleGroup, SampleID, Pathway)]

  grid <- sample_keys[, .(Pathway = pathway_order), by = .(PXD, SampleGroup, SampleID)]
  output <- merge(grid, counts, by = c("PXD", "SampleGroup", "SampleID", "Pathway"), all.x = TRUE)
  output[is.na(PositiveProteinCount), PositiveProteinCount := 0L]
  output[is.na(NegativeProteinCount), NegativeProteinCount := 0L]
  output[is.na(AnyPathwayProteinCount), AnyPathwayProteinCount := 0L]
  output <- merge(output, ddr_totals, by = c("PXD", "SampleGroup", "SampleID"), all.x = TRUE)
  stop_if(!anyNA(output$KlaDdrProteinCount), "A pathway profile denominator is missing.")

  source_meta <- all_records[, .(
    SampleClass = first(SampleClass),
    ObservationType = first(ObservationType),
    SourceMode = paste(sort(unique(SourceMode)), collapse = ";"),
    SourceFile = paste(sort(unique(SourceFile)), collapse = ";")
  ), by = .(PXD, SampleGroup, SampleID)]
  group_meta <- groups[, .(
    RowOrder, PXD, SampleGroup, Category, DisplayGroup = KlaLabelEn
  )]
  output <- merge(output, source_meta, by = c("PXD", "SampleGroup", "SampleID"), all.x = TRUE)
  output <- merge(output, group_meta, by = c("PXD", "SampleGroup"), all.x = TRUE)
  output[, ConditionLabel := fifelse(
    SampleClass %in% c("GSK3B WT", "GSK3B KO", "control", "Roseburia co-culture", "mannitol", "A0h control", "A6h Pg infection"),
    SampleClass,
    SampleID
  )]
  output[, Dataset := "Lactylome (Kla)"]
  output[, `:=`(
    PositiveFraction = PositiveProteinCount / KlaDdrProteinCount,
    NegativeFraction = NegativeProteinCount / KlaDdrProteinCount
  )]
  output[, SignedFraction := PositiveFraction - NegativeFraction]
  stop_if(nrow(output) == nrow(sample_keys) * length(pathway_order),
    "The sample-level pathway profile is incomplete.")
  setcolorder(output, c(
    "RowOrder", "PXD", "SampleGroup", "Category", "DisplayGroup", "Dataset",
    "SampleID", "ConditionLabel", "SampleClass", "ObservationType", "SourceMode", "SourceFile",
    "Pathway", "PositiveProteinCount", "NegativeProteinCount", "AnyPathwayProteinCount",
    "KlaDdrProteinCount", "PositiveFraction", "NegativeFraction", "SignedFraction"
  ))
  output
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

read_maxquant_by_numeric_columns <- function(
  path, pxd, group, sample_map, column_map,
  accession_columns = c("Majority protein IDs", "Protein IDs"),
  source_mode = "deposited_sample_table",
  observation_type = "sample"
) {
  data <- fread(path, check.names = FALSE, showProgress = FALSE)
  accession_column <- intersect(accession_columns, names(data))[[1L]]
  stop_if(!is.na(accession_column), paste0("Missing protein accession column in ", path))
  base_keep <- valid_maxquant_rows(data)
  output <- lapply(seq_len(nrow(sample_map)), function(index) {
    sample_id <- sample_map$SampleID[[index]]
    columns <- column_map[[sample_id]]
    stop_if(length(columns) > 0L && all(columns %in% names(data)),
      paste0("Missing quantitative column for ", sample_id, " in ", path))
    present <- Reduce(`|`, lapply(columns, function(column) positive_numeric(data[[column]])))
    keep <- base_keep & present
    records(
      pxd, group, sample_id, split_accessions(data[[accession_column]][keep]),
      source_mode, path, sample_map$SampleClass[[index]], observation_type
    )
  })
  rbindlist(output, fill = TRUE)
}

read_pxd033146_kla_by_channel <- function(path, pxd, group, sample_map) {
  data <- fread(path, check.names = FALSE, showProgress = FALSE)
  base_keep <- valid_maxquant_rows(data)
  if ("id" %in% names(data)) base_keep <- base_keep & !is.na(data$id)
  if ("Localization prob" %in% names(data)) {
    base_keep <- base_keep & positive_numeric(data$`Localization prob`)
  }
  output <- lapply(seq_len(nrow(sample_map)), function(index) {
    channel <- sample_map$Channel[[index]]
    columns <- paste0("Reporter intensity corrected ", channel, "___", 1:3)
    stop_if(all(columns %in% names(data)),
      paste0("Missing PXD033146 Kla reporter columns for channel ", channel))
    present <- Reduce(`|`, lapply(columns, function(column) positive_numeric(data[[column]])))
    keep <- base_keep & present
    records(
      pxd, group, sample_map$SampleID[[index]], split_accessions(data$Proteins[keep]),
      "deposited_sample_table", path, sample_map$SampleClass[[index]], "sample"
    )
  })
  rbindlist(output, fill = TRUE)
}

read_pd_proteins_by_sample <- function(path, pxd, group, sample_map, accession_column = "Accession") {
  header <- names(read.delim(path, nrows = 0L, check.names = FALSE, stringsAsFactors = FALSE, quote = "\"", comment.char = ""))
  data <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE, quote = "\"", comment.char = "")
  stop_if(accession_column %in% names(data), paste0("Missing ", accession_column, " in ", path))
  found_columns <- header[startsWith(header, "Found in Sample:")]
  output <- lapply(seq_len(nrow(sample_map)), function(index) {
    sample_id <- sample_map$SampleID[[index]]
    hits <- found_columns[
      endsWith(found_columns, paste0(", ", sample_id)) |
        endsWith(found_columns, paste0(", ", sample_id, "_1"))
    ]
    stop_if(length(hits) == 1L, paste0("Missing Proteome Discoverer sample column for ", sample_id))
    keep <- detected_text(data[[hits[[1L]]]])
    records(
      pxd, group, sample_id, split_accessions(data[[accession_column]][keep]),
      "deposited_sample_table", path, sample_map$SampleClass[[index]], "sample"
    )
  })
  rbindlist(output, fill = TRUE)
}

read_tall_by_sample <- function(path, pxd, group, sample_map) {
  data <- fread(path, check.names = FALSE, showProgress = FALSE)
  output <- lapply(seq_len(nrow(sample_map)), function(index) {
    sample_id <- sample_map$SampleID[[index]]
    area_column <- paste0("Area ", sample_id)
    stop_if(all(c("Accession", area_column) %in% names(data)),
      paste0("Missing TALL sample columns for ", sample_id))
    keep <- positive_numeric(data[[area_column]])
    records(
      pxd, group, sample_id, split_accessions(data$Accession[keep]),
      "deposited_sample_table", path, sample_map$SampleClass[[index]], "sample"
    )
  })
  rbindlist(output, fill = TRUE)
}

read_pg_quantity_by_sample <- function(path, pxd, group, sample_map, accession_column = "PG.ProteinGroups") {
  data <- fread(path, check.names = FALSE, showProgress = FALSE)
  stop_if(accession_column %in% names(data), paste0("Missing ", accession_column, " in ", path))
  header <- names(data)
  output <- lapply(seq_len(nrow(sample_map)), function(index) {
    sample_id <- sample_map$SampleID[[index]]
    hits <- header[
      endsWith(header, ".PG.Quantity") & grepl(paste0("_", sample_id, "_"), header, fixed = TRUE)
    ]
    if (!length(hits)) {
      hits <- header[
        endsWith(header, ".PG.Quantity") & grepl(paste0("_", sample_id, "\\."), header, fixed = FALSE)
      ]
    }
    stop_if(length(hits) == 1L, paste0("Expected one protein-quantity column for ", sample_id, " in ", path))
    keep <- positive_numeric(data[[hits[[1L]]]])
    records(
      pxd, group, sample_id, split_accessions(data[[accession_column]][keep]),
      "deposited_sample_table", path, sample_map$SampleClass[[index]], "sample"
    )
  })
  rbindlist(output, fill = TRUE)
}

read_itraq_by_sample <- function(path, pxd, group, sample_map, sheet) {
  data <- as.data.table(read_excel(path, sheet = sheet))
  stop_if("Accession" %in% names(data), paste0("Missing Accession in ", path, " / ", sheet))
  output <- lapply(seq_len(nrow(sample_map)), function(index) {
    sample_id <- sample_map$SampleID[[index]]
    stop_if(sample_id %in% names(data), paste0("Missing iTRAQ sample column for ", sample_id))
    keep <- positive_numeric(data[[sample_id]])
    records(
      pxd, group, sample_id, split_accessions(data$Accession[keep]),
      "deposited_supplementary_table", path, sample_map$SampleClass[[index]], "sample"
    )
  })
  rbindlist(output, fill = TRUE)
}

read_excel_numeric_by_sample <- function(
  path, pxd, group, sample_map, sheet = 1L, skip = 0L,
  accession_column, column_map, source_mode = "deposited_supplementary_table",
  observation_type = "sample"
) {
  data <- as.data.table(read_excel(path, sheet = sheet, skip = skip))
  stop_if(accession_column %in% names(data), paste0("Missing ", accession_column, " in ", path))
  output <- lapply(seq_len(nrow(sample_map)), function(index) {
    sample_id <- sample_map$SampleID[[index]]
    column <- column_map[[sample_id]]
    stop_if(length(column) == 1L && column %in% names(data),
      paste0("Missing quantitative column for ", sample_id, " in ", path))
    keep <- positive_numeric(data[[column]])
    records(
      pxd, group, sample_id, split_accessions(data[[accession_column]][keep]),
      source_mode, path, sample_map$SampleClass[[index]], observation_type
    )
  })
  rbindlist(output, fill = TRUE)
}

read_pxd066517_by_sample <- function(path, pxd, group, sample_map) {
  data <- fread(path, check.names = FALSE, showProgress = FALSE)
  stop_if("PG.ProteinAccessions" %in% names(data), paste0("Missing PG.ProteinAccessions in ", path))
  output <- lapply(seq_len(nrow(sample_map)), function(index) {
    sample_id <- sample_map$SampleID[[index]]
    hits <- names(data)[grepl(paste0("_", sample_id, "\\.raw\\.PG\\.Quantity$"), names(data))]
    stop_if(length(hits) == 1L, paste0("Missing sperm reference sample column for ", sample_id))
    keep <- positive_numeric(data[[hits[[1L]]]])
    records(
      pxd, group, sample_id, split_accessions(data$PG.ProteinAccessions[keep]),
      "deposited_sample_table", path, sample_map$SampleClass[[index]], "sample"
    )
  })
  rbindlist(output, fill = TRUE)
}

read_hk2_by_sample <- function(path, pxd, group, sample_map) {
  data <- fread(path, skip = 2L, check.names = FALSE, showProgress = FALSE)
  stop_if("PG.ProteinGroups" %in% names(data), paste0("Missing PG.ProteinGroups in ", path))
  output <- lapply(seq_len(nrow(sample_map)), function(index) {
    sample_id <- sample_map$SampleID[[index]]
    hits <- names(data)[grepl(paste0("_", sample_id, "\\.raw\\.PG\\.MS2Quantity$"), names(data))]
    stop_if(length(hits) == 1L, paste0("Missing HK-2 reference sample column for ", sample_id))
    keep <- positive_numeric(data[[hits[[1L]]]])
    records(
      pxd, group, sample_id, split_accessions(data$PG.ProteinGroups[keep]),
      "deposited_sample_table", path, sample_map$SampleClass[[index]], "sample"
    )
  })
  rbindlist(output, fill = TRUE)
}

read_pxd073311_proteome_by_sample <- function(path, pxd, group, sample_map) {
  data <- fread(path, check.names = FALSE, showProgress = FALSE)
  stop_if("Protein.Group" %in% names(data), paste0("Missing Protein.Group in ", path))
  raw_columns <- names(data)[grepl("A[06]h_[123]\\.raw$", names(data))]
  output <- lapply(seq_len(nrow(sample_map)), function(index) {
    sample_id <- sample_map$SampleID[[index]]
    hits <- raw_columns[grepl(paste0("_", sample_id, "\\.raw$"), raw_columns)]
    stop_if(length(hits) == 1L, paste0("Missing HUVEC reference sample column for ", sample_id))
    keep <- positive_numeric(data[[hits[[1L]]]])
    records(
      pxd, group, sample_id, split_accessions(data$Protein.Group[keep]),
      "deposited_sample_table", path, sample_map$SampleClass[[index]], "sample"
    )
  })
  rbindlist(output, fill = TRUE)
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
reference_membership <- fread(file.path(input_dir, "reference_protein_membership_30.csv"), check.names = FALSE)
pathway_scores <- read_frozen_pathway_scores(file.path(
  input_dir, "Supplementary_Table_S4_Pathway_Protein_Ranking.xlsx"
))
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

# PXD033146 contains three pathological rotator-cuff tendon samples (RCT1-RCT3)
# and a separate normal-control arm (NC1-NC3).  The Figure 1 tendon row uses
# only the RCT arm.  Each channel has three acquisition columns in the
# deposited site table; those are technical acquisitions of the same sample.
path <- require_file(file.path(source_root, "PXD033146/search_results/extracted_pairing/search_result-HA119TPLa/La (K)Sites.txt"))
map <- sample_map("PXD033146", "pathological rotator cuff tendon", paste0("RCT", 1:3))
map[, Channel := 1:3]
add_records(read_pxd033146_kla_by_channel(path, "PXD033146", "pathological rotator cuff tendon", map))
add_source("PXD033146", "pathological rotator cuff tendon", map$SampleID, map$SampleClass, "MaxQuant Kla site table / RCT channel", path)

# The TALL-104 deposited protein tables expose Sample 1-Sample 3 in both the
# enriched Kla and non-enriched reference arms.  Keep those as three source
# observations rather than collapsing the table to one dataset-level point.
path <- require_file(file.path(source_root, "PXD028488/search_results/Enrichment-Search files/TALL-NALAC-Search files/proteins.csv"))
map <- sample_map("PXD028488", "TALL-104", paste0("Sample ", 1:3))
add_records(read_tall_by_sample(path, "PXD028488", "TALL-104", map))
add_source("PXD028488", "TALL-104", map$SampleID, map$SampleClass, "Deposited protein table / sample area", path)

# Dataset-level observations: the source contains technical/structural runs or
# pooled channels, but it does not support independent biological sample IDs.
single_groups <- list(
  list(pxd = "PXD075377", group = "adjacent liver", id = "Control_pool", class = "pool", file = "PXD075377/search_results/extracted/2-Basic_analysis/MS_identified_information.txt", intensity = "Intensity Control"),
  list(pxd = "PXD075377", group = "HCC", id = "HCC_pool", class = "pool", file = "PXD075377/search_results/extracted/2-Basic_analysis/MS_identified_information.txt", intensity = "Intensity HCC"),
  list(pxd = "PXD028488", group = "HCT116", id = "HCT116_dataset_union", class = "single"),
  list(pxd = "PXD053474", group = "HCT116", id = "HCT116_dataset_union", class = "single"),
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
missing_pathway_scores <- setdiff(
  unique(all_records$BaseAccession[all_records$IsDdr == TRUE]),
  pathway_scores$BaseAccession
)
stop_if(!length(missing_pathway_scores),
  "A source-derived Kla-DDR accession is absent from the frozen pathway-score table.")
pathway_profile <- build_pathway_summary_profile(all_records, groups, pathway_scores)

dataset_boxplot_helper <- file.path(
  project_root, "R", "candidate", "build_dataset_level_boxplot_inputs.R"
)
stop_if(file.exists(dataset_boxplot_helper),
  paste0("Missing dataset-level boxplot helper: ", dataset_boxplot_helper))
source(dataset_boxplot_helper, local = TRUE)
dataset_level_figure1_values <- build_dataset_level_figure1_values(groups)
dataset_level_pathway_values <- build_dataset_level_pathway_summary(
  groups, membership, reference_membership, pathway_scores
)

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

# Build the matched whole-proteome side at the same observation level.  Where
# the reference file contains named source observations, use those columns;
# where it is an averaged cell-line or single-material profile, retain one
# transparent aggregate observation rather than manufacturing replicates.
reference_records_list <- list()
add_reference_records <- function(value, reference_pxd) {
  if (!nrow(value)) return(invisible(NULL))
  value[, ReferencePXD := reference_pxd]
  reference_records_list[[length(reference_records_list) + 1L]] <<- value
  invisible(NULL)
}

reference_membership_long <- rbindlist(lapply(seq_len(nrow(reference_membership)), function(index) {
  accessions <- split_accessions(reference_membership$MappedBaseAccessions[[index]])
  if (!length(accessions)) return(data.table())
  data.table(
    PXD = reference_membership$PXD[[index]],
    SampleGroup = reference_membership$SampleGroup[[index]],
    BaseAccession = accessions,
    IsDdr = is_true_flag(reference_membership$IsDdr[[index]])
  )
}), fill = TRUE)

add_reference_aggregate <- function(row) {
  source_path <- require_file(file.path(
    source_root,
    sub("^data/", "", row$ReferenceEvidenceFile[[1L]])
  ))
  accessions <- reference_membership_long[
    PXD == row$PXD[[1L]] & SampleGroup == row$SampleGroup[[1L]],
    BaseAccession
  ]
  value <- records(
    row$PXD[[1L]], row$SampleGroup[[1L]],
    paste0(row$ReferencePXD[[1L]], "_aggregate"), accessions,
    "validated_reference_membership", source_path, "aggregate", "aggregate"
  )
  add_reference_records(value, row$ReferencePXD[[1L]])
}

# Same-study, sample-resolved references.
path <- require_file(file.path(source_root, "PXD033146/search_results/extracted_pairing/search_result-HA119TQ/proteinGroups.txt"))
map <- sample_map("PXD033146", "pathological rotator cuff tendon", paste0("RCT", 1:3))
map[, Channel := 1:3]
channel_columns <- setNames(
  lapply(map$Channel, function(channel) paste0("Reporter intensity corrected ", channel)),
  map$SampleID
)
add_reference_records(
  read_maxquant_by_numeric_columns(path, "PXD033146", "pathological rotator cuff tendon", map, channel_columns),
  "PXD033146"
)

path <- require_file(file.path(source_root, "PXD046800/search_results/HFX2_LFQ_QB002_Proteins.txt"))
for (group in c("hypertrophic scar", "adjacent skin")) {
  ids <- if (group == "hypertrophic scar") paste0("HSP", 1:4) else paste0("NSP", 1:4)
  map <- sample_map("PXD046800", group, ids)
  add_reference_records(read_pd_proteins_by_sample(path, "PXD046800", group, map), "PXD046800")
}

path <- require_file(file.path(source_root, "PXD050470/supplementary/prca2331-sup-0006-tables4.xlsx"))
map <- sample_map("PXD050470", "human hippocampus", c("H072", "H081", "H187"))
intensity_map <- setNames(paste0("Intensity_", map$SampleID), map$SampleID)
add_reference_records(
  read_excel_numeric_by_sample(
    path, "PXD050470", "human hippocampus", map, sheet = "Sheet1", skip = 5L,
    accession_column = "Protein accession", column_map = intensity_map
  ),
  "PXD050470"
)

path <- require_file(file.path(source_root, "PXD066517/search_results/20240275.tsv"))
ids <- c(paste0("L", 1:17), paste0("N", 1:16))
map <- sample_map("PXD064912", "human sperm", ids, c(rep("L sperm", 17), rep("N sperm", 16)))
add_reference_records(read_pxd066517_by_sample(path, "PXD064912", "human sperm", map), "PXD066517")

path <- require_file(file.path(source_root, "PXD066054/search_results/extracted/DA/Protein_Quant.tsv"))
for (group in c("BPH", "prostate cancer")) {
  ids <- if (group == "BPH") paste0("NAT", 1:5) else paste0("PCa", 1:5)
  map <- sample_map("PXD066054", group, ids)
  add_reference_records(read_pg_quantity_by_sample(path, "PXD066054", group, map), "PXD066054")
}

path <- require_file(file.path(source_root, "PXD065775/search_results/20170330_01-24_patients_iTRAQ.xlsx"))
ids <- c(paste0("Non-rec", 1:4), paste0("Rec", 1:4))
classes <- c(rep("Non-rec", 4), rep("Rec", 4))
for (item in list(
  list(group = "adjacent liver", sheet = "ANTs"),
  list(group = "HCC", sheet = "CISs")
)) {
  map <- sample_map("PXD075377", item$group, ids, classes)
  add_reference_records(read_itraq_by_sample(path, "PXD075377", item$group, map, item$sheet), "PXD065775")
}

path <- require_file(file.path(source_root, "PXD066351/search_results/XB01472B1DA-Protein_Quant.tsv"))
map <- sample_map("PXD066351", "HCT116 control and Roseburia co-culture", c("NC116", "R116"), c("control", "Roseburia co-culture"))
add_reference_records(read_pg_quantity_by_sample(path, "PXD066351", "HCT116 control and Roseburia co-culture", map), "PXD066351")

path <- require_file(file.path(source_root, "PXD028488/search_results/Nonenrichment-Search files/TALL-Nonenrichment-Search files/proteins.csv"))
map <- sample_map("PXD028488", "TALL-104", paste0("Sample ", 1:3))
add_reference_records(read_tall_by_sample(path, "PXD028488", "TALL-104", map), "PXD028488")

path <- require_file(file.path(source_root, "PXD050147/search_results/SIRT_proteinGroups.txt"))
ids <- c(
  "WT_pro_rep1", "WT_pro_rep2", "WT_pro_rep3",
  "SIRT1KO_pro_rep1", "SIRT1KO_pro_rep2", "SIRT1KO_pro_rep3",
  "SIRT3KO_pro_rep1", "SIRT3KO_pro_rep2", "SIRT3KO_pro_rep3"
)
classes <- c(rep("WT", 3), rep("SIRT1 KO", 3), rep("SIRT3 KO", 3))
map <- sample_map("PXD050147", "HepG2 WT and SIRT1 or SIRT3 KO", ids, classes)
intensity_map <- setNames(paste0("Intensity ", ids), ids)
add_reference_records(
  read_maxquant_by_numeric_columns(path, "PXD050147", "HepG2 WT and SIRT1 or SIRT3 KO", map, intensity_map),
  "PXD050147"
)

path <- require_file(file.path(source_root, "PXD028737/search_results/extracted_reference/txt/proteinGroups.txt"))
map <- sample_map("PXD028737", "HMC3", c("H0", "H24"), c("normoxia", "hypoxia"))
intensity_map <- setNames(paste0("LFQ intensity ", map$SampleID), map$SampleID)
add_reference_records(
  read_maxquant_by_numeric_columns(path, "PXD028737", "HMC3", map, intensity_map),
  "PXD028737"
)

path <- require_file(file.path(source_root, "PXD072220/search_results/HK-2_Spectronaut-report_PG_Quantity.txt"))
map <- sample_map("PXD058534", "pretreated HK-2", c("amostra1", "amostra3", "amostra4"), "untreated control")
add_reference_records(read_hk2_by_sample(path, "PXD058534", "pretreated HK-2", map), "PXD072220")
map <- sample_map("PXD078736", "HK-2 control and mannitol", c("amostra1", "amostra3", "amostra4"), "untreated control")
add_reference_records(read_hk2_by_sample(path, "PXD078736", "HK-2 control and mannitol", map), "PXD072220")

path <- require_file(file.path(source_root, "PXD069969/search_results/SA206LQB1_Annotation.xlsx"))
for (item in list(
  list(group = "glioblastoma stem cells", ids = c("G2907", "G3028", "G3264", "GSC23", "MES28", "RKI")),
  list(group = "neural stem cells", ids = c("ENSA", "HMP1"))
)) {
  map <- sample_map("PXD070007", item$group, item$ids, "model")
  intensity_map <- setNames(paste0("LFQ intensity ", item$ids), item$ids)
  add_reference_records(
    read_excel_numeric_by_sample(
      path, "PXD070007", item$group, map, sheet = "Annotation_Combine",
      accession_column = "Protein accession", column_map = intensity_map
    ),
    "PXD069969"
  )
}

path <- require_file(file.path(source_root, "PXD073311/search_results/extracted_pairing/IPX0015307001_Database_search_result/Database_search_result/report.pg_matrix.tsv"))
map <- sample_map("PXD073311", "HUVEC control and Pg infection", paste0("A0h_", 1:3), "A0h control")
add_reference_records(read_pxd073311_proteome_by_sample(path, "PXD073311", "HUVEC control and Pg infection", map), "PXD073311")

direct_reference_keys <- c(
  "PXD033146__pathological rotator cuff tendon",
  "PXD046800__hypertrophic scar", "PXD046800__adjacent skin",
  "PXD050470__human hippocampus", "PXD064912__human sperm",
  "PXD066054__BPH", "PXD066054__prostate cancer",
  "PXD075377__adjacent liver", "PXD075377__HCC",
  "PXD066351__HCT116 control and Roseburia co-culture",
  "PXD028488__TALL-104", "PXD050147__HepG2 WT and SIRT1 or SIRT3 KO",
  "PXD028737__HMC3", "PXD058534__pretreated HK-2",
  "PXD078736__HK-2 control and mannitol",
  "PXD070007__glioblastoma stem cells", "PXD070007__neural stem cells",
  "PXD073311__HUVEC control and Pg infection"
)
for (index in seq_len(nrow(groups))) {
  row <- groups[index]
  key <- paste(row$PXD, row$SampleGroup, sep = "__")
  if (!(key %in% direct_reference_keys)) add_reference_aggregate(row)
}

reference_records <- unique(rbindlist(reference_records_list, fill = TRUE))
reference_records[, SourceFile := vapply(SourceFile, source_rel, character(1))]
reference_records[, GroupKey := paste(PXD, SampleGroup, sep = "__")]
stop_if(setequal(unique(reference_records$GroupKey), groups$GroupKey),
  "Reference preparation did not cover exactly the frozen groups.")
stop_if(!anyDuplicated(reference_records[, .(PXD, SampleGroup, SampleID, BaseAccession)]),
  "Duplicate reference sample-level membership rows detected.")

reference_records <- merge(
  reference_records,
  unique(reference_membership_long[, .(PXD, SampleGroup, BaseAccession, IsDdr)]),
  by = c("PXD", "SampleGroup", "BaseAccession"),
  all.x = TRUE,
  sort = FALSE
)
reference_records[is.na(IsDdr), IsDdr := FALSE]

reference_values <- reference_records[, .(
  ReferencePXD = first(ReferencePXD),
  WholeProteomeProteinCount = uniqueN(BaseAccession),
  WholeProteomeDdrProteinCount = uniqueN(BaseAccession[IsDdr == TRUE]),
  SourceMode = first(SourceMode),
  SourceFile = first(SourceFile),
  SampleClass = first(SampleClass),
  ObservationType = first(ObservationType)
), by = .(PXD, SampleGroup, SampleID)]
reference_values[, WholeProteomeDdrFraction := fifelse(
  WholeProteomeProteinCount > 0L,
  WholeProteomeDdrProteinCount / WholeProteomeProteinCount,
  NA_real_
)]
reference_values[, WholeProteomeDdrFractionPercentage := WholeProteomeDdrFraction * 100]
reference_values[, ConditionLabel := SampleID]

plot_group_meta <- groups[, .(
  RowOrder, PXD, SampleGroup, Category, DisplayGroup = KlaLabelEn,
  ReferencePXD, FrozenKlaProteinCount = KlaProteinCount,
  FrozenKlaDdrProteinCount = KlaDdrProteinCount,
  FrozenKlaDdrFraction = KlaDdrFraction * 100,
  FrozenReferenceProteinCount = ReferenceProteinCount,
  FrozenReferenceDdrProteinCount = ReferenceDdrProteinCount,
  FrozenReferenceDdrFraction = ReferenceDdrFraction * 100
)]

kla_plot_values <- merge(
  sample_values[, .(
    RowOrder, PXD, SampleGroup, Category, DisplayGroup, SampleID, ConditionLabel,
    SampleClass, ObservationType, SourceMode, SourceFile,
    ProteinCount = KlaProteinCount, DdrProteinCount = KlaDdrProteinCount,
    DdrFraction = KlaDdrFraction, DdrFractionPercentage = KlaDdrFractionPercentage
  )],
  plot_group_meta,
  by = c("RowOrder", "PXD", "SampleGroup", "Category", "DisplayGroup"),
  all.x = TRUE,
  sort = FALSE
)
kla_plot_values[, Dataset := "Lactylome (Kla)"]

reference_plot_values <- merge(
  reference_values,
  plot_group_meta,
  by = c("PXD", "SampleGroup"),
  all.x = TRUE,
  sort = FALSE
)
reference_plot_values[, ReferencePXD := ReferencePXD.x]
reference_plot_values[, c("ReferencePXD.x", "ReferencePXD.y") := NULL]
reference_plot_values[, `:=`(
  ProteinCount = WholeProteomeProteinCount,
  DdrProteinCount = WholeProteomeDdrProteinCount,
  DdrFraction = WholeProteomeDdrFraction,
  DdrFractionPercentage = WholeProteomeDdrFractionPercentage,
  Dataset = "Whole proteome"
)]

figure1_values <- rbindlist(list(
  kla_plot_values[, .(
    RowOrder, PXD, SampleGroup, Category, DisplayGroup, Dataset,
    SampleID, ConditionLabel, SampleClass, ObservationType, SourceMode, SourceFile,
    ReferencePXD, ProteinCount, DdrProteinCount, DdrFraction, DdrFractionPercentage,
    FrozenKlaProteinCount, FrozenKlaDdrProteinCount, FrozenKlaDdrFraction,
    FrozenReferenceProteinCount, FrozenReferenceDdrProteinCount, FrozenReferenceDdrFraction
  )],
  reference_plot_values[, .(
    RowOrder, PXD, SampleGroup, Category, DisplayGroup, Dataset,
    SampleID, ConditionLabel, SampleClass, ObservationType, SourceMode, SourceFile,
    ReferencePXD, ProteinCount, DdrProteinCount, DdrFraction, DdrFractionPercentage,
    FrozenKlaProteinCount, FrozenKlaDdrProteinCount, FrozenKlaDdrFraction,
    FrozenReferenceProteinCount, FrozenReferenceDdrProteinCount, FrozenReferenceDdrFraction
  )]
), fill = TRUE)
setorder(figure1_values, RowOrder, Dataset, SampleID)

stop_if(nrow(figure1_values[Dataset == "Lactylome (Kla)"]) == 92L,
  "The source-defined Kla sample count must be 92.")
stop_if(nrow(figure1_values[Dataset == "Whole proteome"]) == 118L,
  "The source-defined whole-proteome sample count must be 118.")
stop_if(all(is.finite(figure1_values$DdrFractionPercentage)),
  "A Figure 1 sample fraction is not finite.")
stop_if(all(figure1_values$DdrFractionPercentage >= 0 & figure1_values$DdrFractionPercentage <= 100),
  "A Figure 1 sample fraction is outside 0-100 percent.")

reference_registry <- reference_records[, .(
  SampleClass = first(SampleClass),
  Parser = paste0("whole-proteome / ", first(SourceMode)),
  SourceFile = first(SourceFile)
), by = .(PXD, SampleGroup, SampleID)]
reference_registry[, Dataset := "Whole proteome"]
registry_table <- rbindlist(registry, fill = TRUE)
registry_table[, Dataset := "Lactylome (Kla)"]
figure1_registry <- rbindlist(list(registry_table, reference_registry), fill = TRUE)
setorder(figure1_registry, Dataset, PXD, SampleGroup, SampleID)

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
fwrite(figure1_values, file.path(candidate_dir, "figure1_sample_boxplot_values.csv"), na = "")
fwrite(figure1_registry, file.path(candidate_dir, "figure1_sample_boxplot_source_registry.csv"), na = "")
fwrite(pathway_profile, file.path(candidate_dir, "figure1_pathway_summary_sample_boxplot_values.csv"), na = "")
fwrite(dataset_level_figure1_values, file.path(candidate_dir, "figure1_dataset_boxplot_values.csv"), na = "")
fwrite(dataset_level_pathway_values, file.path(candidate_dir, "pathway_summary_dataset_boxplot_values.csv"), na = "")

message(
  "Wrote sample-level boxplot inputs for ", nrow(sample_values),
  " observations across ", uniqueN(sample_values[, .(PXD, SampleGroup)]),
  " publication groups and ", nrow(pathway_profile),
  " sample-level pathway records."
)
