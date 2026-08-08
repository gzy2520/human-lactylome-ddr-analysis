#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(readxl)
  library(stringr)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) normalizePath(args[[1]]) else normalizePath(".")
table_dir <- file.path(project_root, "reanalysis", "results", "tables")
figure_dir <- file.path(project_root, "reanalysis", "results", "figures")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

regulator_path <- file.path(
  project_root, "data", "identifier",
  "乳酸化调控因子_Writer-Eraser-Reader.xlsx"
)
mapping_path <- file.path(
  project_root, "reanalysis", "config",
  "lactylation_regulator_uniprot_mapping.csv"
)
ensembl_mapping_path <- file.path(
  project_root, "reanalysis", "config",
  "ensembl_protein_to_uniprot_biomart.tsv"
)
detection_path <- file.path(
  table_dir, "kla_regulator_cell_type_long.csv"
)
kla_scope_path <- file.path(
  table_dir, "kla_regulator_intensity_availability_audit.csv"
)
reviewed_uniprot_path <- file.path(
  project_root, "reanalysis", "config",
  "uniprot_human_reviewed_2026-08-05.tsv"
)

required_files <- c(
  regulator_path,
  mapping_path,
  ensembl_mapping_path,
  detection_path,
  kla_scope_path,
  reviewed_uniprot_path
)
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

accession_map <- read.csv(
  mapping_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
) |>
  transmute(
    BaseAccession = sub("-[0-9]+$", "", BaseAccession),
    GeneSymbol = toupper(GeneSymbol)
  )
accession_to_gene <- setNames(accession_map$GeneSymbol, accession_map$BaseAccession)
regulators <- regulators |>
  left_join(accession_map, by = "GeneSymbol")
if (any(is.na(regulators$BaseAccession))) {
  stop(
    "Missing UniProt mapping for: ",
    paste(regulators$GeneSymbol[is.na(regulators$BaseAccession)], collapse = ", ")
  )
}
target_accessions <- unique(regulators$BaseAccession)

is_uniprot <- function(values) {
  grepl(
    paste0(
      "^(?:",
      "[OPQ][0-9][A-Z0-9]{3}[0-9]|",
      "[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9](?:[A-Z0-9]{3}[0-9])?",
      ")$"
    ),
    values
  )
}

normalize_ensembl_protein <- function(values) {
  values <- trimws(as.character(values))
  values <- sub("^ACC:", "", values)
  values <- sub("_RNA$", "", values)
  values <- sub("\\.[0-9]+$", "", values)
  values
}

ensembl_mapping <- read.delim(
  ensembl_mapping_path,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = ""
) |>
  transmute(
    EnsemblProteinID = normalize_ensembl_protein(EnsemblProteinID),
    BaseAccession = sub("-[0-9]+$", "", as.character(BaseAccession))
  ) |>
  filter(
    grepl("^ENSP[0-9]+$", EnsemblProteinID),
    is_uniprot(BaseAccession)
  ) |>
  distinct()
ensembl_to_uniprot <- split(
  ensembl_mapping$BaseAccession,
  ensembl_mapping$EnsemblProteinID
)

detection <- read.csv(
  detection_path,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  na.strings = c("", "NA")
)
kla_scope <- read.csv(
  kla_scope_path,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  na.strings = c("", "NA")
) |>
  filter(`定量可用` %in% c(TRUE, "TRUE", "True", 1, "1")) |>
  transmute(
    PXD,
    SampleGroup = `样本组`,
    KlaRowLabel = `行标签`
  )
if (nrow(kla_scope) != 37) {
  stop("Kla quantitative heatmap scope must contain exactly 37 sample groups")
}
pairing_path <- file.path(
  project_root,
  "reanalysis",
  "config",
  "lactylome_reference_pairing.csv"
)
if (!file.exists(pairing_path)) {
  stop("Missing lactylome/reference pairing table: ", pairing_path)
}
pairing <- read.csv(
  pairing_path,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  na.strings = c("", "NA")
) |>
  filter(
    !is.na(LactylomePXD),
    !is.na(SampleGroup)
  )
sample_catalog_all <- detection |>
  distinct(
    PXD, SampleGroup, SampleGroupID, RowLabel,
    BiologicalMaterial, GeneLevelAuditStatus
  ) |>
  inner_join(kla_scope, by = c("PXD", "SampleGroup")) |>
  left_join(
    pairing |>
      transmute(
        PXD = LactylomePXD,
        SampleGroup,
        ReferencePXD,
        PairingInclude = IncludeInStrictReferenceAnalysis,
        ReferenceMatchQuality = MatchQuality,
        ReferenceCaveat = Caveat
      ),
    by = c("PXD", "SampleGroup")
  ) |>
  mutate(RowOrder = row_number())
if (nrow(sample_catalog_all) != 37) {
  stop("Kla-derived candidate scope must contain exactly 37 sample groups")
}
write.csv(
  sample_catalog_all |>
    filter(
      PairingInclude %in% c(FALSE, "FALSE", "False", 0, "0") |
        is.na(ReferencePXD) |
        !nzchar(ReferencePXD)
    ) |>
    transmute(
      PXD,
      SampleGroup,
      Decision = "excluded_no_exact_quantitative_reference",
      Reason = ReferenceCaveat
    ),
  file.path(
    table_dir,
    "kla_regulator_whole_proteome_strict_scope_exclusions.csv"
  ),
  row.names = FALSE,
  na = ""
)
sample_catalog <- sample_catalog_all |>
  filter(
    PairingInclude %in% c(TRUE, "TRUE", "True", 1, "1"),
    !is.na(ReferencePXD),
    nzchar(ReferencePXD)
  )
if (nrow(sample_catalog) != 33) {
  stop(
    "Strict whole-proteome heatmap scope must contain 33 sample groups, found ",
    nrow(sample_catalog)
  )
}
sample_catalog$RowLabel <- ifelse(
  !is.na(sample_catalog$ReferencePXD) &
    nzchar(sample_catalog$ReferencePXD),
  paste0(
    sample_catalog$SampleGroup,
    " · Ref:", sample_catalog$ReferencePXD,
    " · Kla:", sample_catalog$PXD
  ),
  paste0(
    sample_catalog$SampleGroup,
    " · Ref:未找到严格参照",
    " · Kla:", sample_catalog$PXD
  )
)

base_accession <- function(values) {
  values <- as.character(values)
  values <- sub("^.*\\|([^|]+)\\|.*$", "\\1", values)
  values <- sub("^([^|;]+)\\|.*$", "\\1", values)
  sub("-[0-9]+$", "", values)
}

match_target_accession <- function(values) {
  vapply(as.character(values), function(value) {
    if (is.na(value) || !nzchar(value)) return(NA_character_)
    tokens <- unique(trimws(unlist(strsplit(value, "[;, ]+"))))
    direct_tokens <- base_accession(tokens)
    direct_hits <- direct_tokens[direct_tokens %in% target_accessions]
    ensembl_tokens <- normalize_ensembl_protein(tokens)
    ensembl_tokens <- ensembl_tokens[
      grepl("^ENSP[0-9]+$", ensembl_tokens)
    ]
    mapped_hits <- unlist(
      ensembl_to_uniprot[intersect(
        ensembl_tokens,
        names(ensembl_to_uniprot)
      )],
      use.names = FALSE
    )
    hits <- unique(c(
      direct_hits,
      mapped_hits[mapped_hits %in% target_accessions]
    ))
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

safe_max <- function(values) {
  values <- values[is.finite(values)]
  if (length(values)) max(values) else NA_real_
}

relative_path <- function(path) {
  sub(paste0("^", project_root, "/"), "", normalizePath(path))
}

reviewed_uniprot <- data.table::fread(
  reviewed_uniprot_path,
  sep = "\t",
  quote = "",
  data.table = FALSE,
  check.names = FALSE,
  showProgress = FALSE
)
reviewed_primary_summary <- reviewed_uniprot |>
  transmute(
    BaseAccession = base_accession(Entry),
    Symbols = `Gene Names (primary)`
  ) |>
  separate_rows(Symbols, sep = "[;[:space:]]+") |>
  transmute(
    SourceSymbol = toupper(trimws(Symbols)),
    BaseAccession
  ) |>
  filter(nzchar(SourceSymbol)) |>
  distinct() |>
  group_by(SourceSymbol) |>
  summarise(
    CandidateCount = n_distinct(BaseAccession),
    CandidateAccessions = paste(
      sort(unique(BaseAccession)),
      collapse = ";"
    ),
    .groups = "drop"
  )
reviewed_alias_summary <- reviewed_uniprot |>
  transmute(
    BaseAccession = base_accession(Entry),
    Symbols = `Gene Names`
  ) |>
  separate_rows(Symbols, sep = "[;[:space:]]+") |>
  transmute(
    SourceSymbol = toupper(trimws(Symbols)),
    BaseAccession
  ) |>
  filter(nzchar(SourceSymbol)) |>
  distinct() |>
  group_by(SourceSymbol) |>
  summarise(
    CandidateCount = n_distinct(BaseAccession),
    CandidateAccessions = paste(
      sort(unique(BaseAccession)),
      collapse = ";"
    ),
    .groups = "drop"
  )

map_reviewed_symbol_features <- function(source_features) {
  data.frame(
    SourceFeature = as.character(source_features),
    stringsAsFactors = FALSE
  ) |>
    mutate(SourceSymbol = strsplit(SourceFeature, "_", fixed = TRUE)) |>
    unnest(SourceSymbol) |>
    mutate(SourceSymbol = toupper(trimws(SourceSymbol))) |>
    filter(nzchar(SourceSymbol)) |>
    distinct() |>
    left_join(
      reviewed_primary_summary |>
        rename(
          PrimaryCandidateCount = CandidateCount,
          PrimaryCandidates = CandidateAccessions
        ),
      by = "SourceSymbol"
    ) |>
    left_join(
      reviewed_alias_summary |>
        rename(
          AliasCandidateCount = CandidateCount,
          AliasCandidates = CandidateAccessions
        ),
      by = "SourceSymbol"
    ) |>
    mutate(
      MappingMode = case_when(
        !is.na(PrimaryCandidateCount) & PrimaryCandidateCount == 1 ~
          "reviewed_primary_symbol",
        is.na(PrimaryCandidateCount) & AliasCandidateCount == 1 ~
          "reviewed_unique_alias",
        !is.na(PrimaryCandidateCount) & PrimaryCandidateCount > 1 ~
          "ambiguous_reviewed_primary",
        AliasCandidateCount > 1 ~ "ambiguous_reviewed_alias",
        TRUE ~ "unmapped"
      ),
      BaseAccession = case_when(
        MappingMode == "reviewed_primary_symbol" ~ PrimaryCandidates,
        MappingMode == "reviewed_unique_alias" ~ AliasCandidates,
        TRUE ~ ""
      ),
      TargetAccession = ifelse(
        BaseAccession %in% target_accessions,
        BaseAccession,
        NA_character_
      ),
      MappingSource =
        "UniProtKB reviewed Homo sapiens snapshot 2026-08-05"
    )
}

quant_parts <- list()
quant_audit <- sample_catalog |>
  transmute(
    PXD,
    SampleGroup,
    SampleGroupID,
    WholeProteomeQuantAvailable = FALSE,
    WholeProteomeSource = "",
    WholeProteomeMeasurement = "",
    WholeProteomeStatus = "unavailable",
    WholeProteomeReason =
      "未配置可解析的普通全蛋白定量文件；不使用Kla信号替代"
  ) |>
  left_join(
    pairing |>
      transmute(
        PXD = LactylomePXD,
        SampleGroup,
        ReferencePXD = ReferencePXD,
        ReferenceSampleSubset,
        ReferenceEvidenceLocator,
        ReferenceProteinCount,
        ReferenceMatchQuality = MatchQuality,
        ReferenceCaveat = Caveat,
        PairingInclude = IncludeInStrictReferenceAnalysis
      ),
    by = c("PXD", "SampleGroup")
  )
quant_audit$WholeProteomeReason[
  quant_audit$PairingInclude %in% c(FALSE, "FALSE", "False", 0, "0") &
    quant_audit$ReferenceMatchQuality == "no_exact_reference_found"
] <- "未找到同组织、同细胞系或同生物状态的严格普通全蛋白参照；按老师要求留空"

update_audit <- function(
  pxd,
  sample_group,
  source_file,
  measurement,
  reason = "同研究或同样本普通全蛋白定量可用"
) {
  rows <- quant_audit$PXD == pxd & quant_audit$SampleGroup == sample_group
  if (!any(rows)) return(invisible(NULL))
  quant_audit$WholeProteomeQuantAvailable[rows] <<- TRUE
  old_sources <- quant_audit$WholeProteomeSource[rows]
  old_measurements <- quant_audit$WholeProteomeMeasurement[rows]
  quant_audit$WholeProteomeSource[rows] <<- vapply(
    old_sources,
    function(old) paste(
      sort(unique(c(
        strsplit(old, ";", fixed = TRUE)[[1]][nzchar(old)],
        source_file
      ))),
      collapse = ";"
    ),
    character(1)
  )
  quant_audit$WholeProteomeMeasurement[rows] <<- vapply(
    old_measurements,
    function(old) paste(
      sort(unique(c(
        strsplit(old, ";", fixed = TRUE)[[1]][nzchar(old)],
        measurement
      ))),
      collapse = ";"
    ),
    character(1)
  )
  quant_audit$WholeProteomeStatus[rows] <<- "available"
  quant_audit$WholeProteomeReason[rows] <<- reason
  invisible(NULL)
}

add_total_quant <- function(
  data,
  pxd,
  sample_group,
  feature_values,
  target_values,
  sample_columns,
  sample_labels = sample_columns,
  measurement,
  source_file,
  signal_is_log2 = FALSE
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
      Measurement = measurement,
      SignalIsLog2 = signal_is_log2,
      SourceFile = relative_path(source_file),
      CanonicalFeature = ifelse(
        !is.na(TargetAccession) & nzchar(TargetAccession),
        paste0("ACC:", TargetAccession),
        FeatureID
      )
    ) |>
    group_by(
    PXD, SampleGroup, QuantSample, CanonicalFeature,
    TargetAccession, TargetGene, Measurement, SignalIsLog2, SourceFile
  ) |>
    summarise(Signal = sum(Signal, na.rm = TRUE), .groups = "drop")
  if (nrow(long)) {
    quant_parts[[length(quant_parts) + 1]] <<- long
    update_audit(
      pxd,
      sample_group,
      relative_path(source_file),
      measurement
    )
  }
  invisible(NULL)
}

add_maxquant_proteome <- function(
  path,
  pxd,
  sample_group,
  token_map,
  measurement = "MaxQuant proteinGroups"
) {
  if (!file.exists(path)) return(invisible(NULL))
  data <- read.delim(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = ""
  )
  keep <- rep(TRUE, nrow(data))
  if ("Reverse" %in% names(data)) keep <- keep & data$Reverse != "+"
  if ("Potential contaminant" %in% names(data)) {
    keep <- keep & data$`Potential contaminant` != "+"
  }
  if ("Only identified by site" %in% names(data)) {
    keep <- keep & data$`Only identified by site` != "+"
  }
  if ("id" %in% names(data)) keep <- keep & !is.na(data$id)
  data <- data[keep, , drop = FALSE]
  accession_col <- intersect(
    c("Majority protein IDs", "Protein IDs"),
    names(data)
  )
  if (!length(accession_col)) return(invisible(NULL))
  accession_col <- accession_col[[1]]
  for (spec in token_map) {
    token <- spec[[1]]
    label <- spec[[2]]
    candidate_columns <- c(
      token,
      paste("LFQ intensity", token),
      paste("Intensity", token),
      paste("iBAQ", token)
    )
    column <- candidate_columns[candidate_columns %in% names(data)]
    if (!length(column)) next
    add_total_quant(
      data,
      pxd,
      sample_group,
      accession_feature(data[[accession_col]]),
      match_target_accession(data[[accession_col]]),
      column[[1]],
      label,
      measurement,
      path
    )
  }
}

add_pd_proteome <- function(path, pxd, group_map) {
  if (!file.exists(path)) return(invisible(NULL))
  data <- read.delim(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    quote = "\"",
    comment.char = ""
  )
  for (spec in group_map) {
    group <- spec[[1]]
    token <- spec[[2]]
    columns <- grep(
      paste0("^Abundances \\(Normalized\\):.*", token),
      names(data),
      value = TRUE
    )
    add_total_quant(
      data,
      pxd,
      group,
      accession_feature(data$Accession),
      match_target_accession(data$Accession),
      columns,
      sub("^.*Sample, ", "", columns),
      "Proteome Discoverer normalized total-protein abundance",
      path
    )
  }
}

add_spectronaut_report <- function(path, pxd, sample_group, measurement) {
  if (!file.exists(path)) return(invisible(NULL))
  data <- data.table::fread(path, data.table = FALSE, showProgress = FALSE)
  if (!all(c("Protein.Group", "PG.Quantity") %in% names(data))) {
    return(invisible(NULL))
  }
  if ("PG.Q.Value" %in% names(data)) {
    data <- data[safe_numeric(data$PG.Q.Value) <= 0.01, , drop = FALSE]
  }
  run_column <- if ("Run" %in% names(data)) "Run" else "File.Name"
  for (run in unique(data[[run_column]])) {
    subset <- data[data[[run_column]] == run, , drop = FALSE] |>
      group_by(Protein.Group) |>
      summarise(PG.Quantity = max(safe_numeric(PG.Quantity)), .groups = "drop")
    add_total_quant(
      subset,
      pxd,
      sample_group,
      accession_feature(subset$Protein.Group),
      match_target_accession(subset$Protein.Group),
      "PG.Quantity",
      as.character(run),
      measurement,
      path
    )
  }
}

add_spectronaut_standard_report <- function(
  path,
  pxd,
  sample_group,
  measurement,
  column_pattern = NULL,
  skip_lines = 0,
  signal_is_log2 = FALSE
) {
  if (!file.exists(path)) return(invisible(NULL))
  data <- data.table::fread(
    path,
    skip = skip_lines,
    data.table = FALSE,
    showProgress = FALSE
  )
  accession_column <- intersect(
    c("PG.ProteinAccessions", "PG.ProteinGroups", "Protein.Group"),
    names(data)
  )
  quantity_columns <- grep(
    "PG\\.(Quantity|Log2Quantity|MS2Quantity)$",
    names(data),
    value = TRUE
  )
  if (!is.null(column_pattern)) {
    quantity_columns <- quantity_columns[
      grepl(column_pattern, quantity_columns, ignore.case = TRUE)
    ]
  }
  if (!length(accession_column) || !length(quantity_columns)) {
    return(invisible(NULL))
  }
  accession_column <- accession_column[[1]]
  for (column in quantity_columns) {
    sample_label <- sub("^\\[[0-9]+\\] ", "", column)
    sample_label <- sub(
      "\\.PG\\.(Quantity|Log2Quantity|MS2Quantity)$",
      "",
      sample_label
    )
    subset <- data |>
      transmute(
        ProteinAccessions = .data[[accession_column]],
        Quantity = safe_numeric(.data[[column]])
      ) |>
      group_by(ProteinAccessions) |>
      summarise(Quantity = safe_max(Quantity), .groups = "drop")
    add_total_quant(
      subset,
      pxd,
      sample_group,
      accession_feature(subset$ProteinAccessions),
      match_target_accession(subset$ProteinAccessions),
      "Quantity",
      sample_label,
      measurement,
      path,
      signal_is_log2
    )
  }
}

add_spectronaut_matrix <- function(
  path,
  pxd,
  sample_group,
  measurement,
  column_pattern = NULL
) {
  if (!file.exists(path)) return(invisible(NULL))
  data <- read.delim(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = ""
  )
  columns <- grep("PG\\.Quantity$", names(data), value = TRUE)
  if (!is.null(column_pattern)) {
    columns <- columns[grepl(column_pattern, columns)]
  }
  if (!length(columns)) return(invisible(NULL))
  add_total_quant(
    data,
    pxd,
    sample_group,
    accession_feature(data$PG.ProteinGroups),
    match_target_accession(data$PG.ProteinGroups),
    columns,
    sub("^.*\\] |\\.PG\\.Quantity$", "", columns),
    measurement,
    path
  )
}

add_peaks_proteins <- function(paths, pxd, sample_group) {
  for (path in paths) {
    data <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
    area_columns <- grep("^Area", names(data), value = TRUE)
    if (!length(area_columns) || !"Accession" %in% names(data)) next
    add_total_quant(
      data,
      pxd,
      sample_group,
      accession_feature(data$Accession),
      match_target_accession(data$Accession),
      area_columns,
      paste0(sample_group, "__", basename(dirname(path)), "__", sub("^Area ", "", area_columns)),
      "PEAKS non-enriched protein Area",
      path
    )
  }
}

add_pxd043880_hippocampus <- function(path) {
  if (!file.exists(path)) return(invisible(NULL))
  data <- read_excel(
    path,
    sheet = "Source Data Proteins",
    col_names = FALSE,
    col_types = "text",
    .name_repair = "minimal"
  )
  if (nrow(data) < 3 || ncol(data) < 9) {
    stop("Unexpected PXD043880 protein matrix layout: ", path)
  }

  donor_rows <- 3:nrow(data)
  donor_ids <- trimws(as.character(data[[1]][donor_rows]))
  donor_keep <- !is.na(donor_ids) &
    nzchar(donor_ids) &
    !grepl("\\((removed|outlier)\\)", donor_ids, ignore.case = TRUE)
  donor_rows <- donor_rows[donor_keep]
  donor_ids <- donor_ids[donor_keep]
  if (length(donor_ids) != 74) {
    stop("Expected 74 retained PXD043880 donors, found ", length(donor_ids))
  }

  source_features <- trimws(as.character(
    unlist(data[2, 9:ncol(data)], use.names = FALSE)
  ))
  feature_keep <- !is.na(source_features) & nzchar(source_features)
  source_features <- source_features[feature_keep]
  intensity <- as.matrix(
    data[donor_rows, 9:ncol(data), drop = FALSE]
  )[, feature_keep, drop = FALSE]
  quant_data <- as.data.frame(
    t(intensity),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  sample_labels <- make.unique(paste0("CA1_", donor_ids))
  names(quant_data) <- sample_labels

  mapping <- map_reviewed_symbol_features(source_features)
  target_by_feature <- mapping |>
    group_by(SourceFeature) |>
    summarise(
      TargetHitCount = n_distinct(TargetAccession[!is.na(TargetAccession)]),
      TargetAccession = paste(
        sort(unique(TargetAccession[!is.na(TargetAccession)])),
        collapse = ";"
      ),
      .groups = "drop"
    )
  if (any(target_by_feature$TargetHitCount > 1)) {
    stop(
      "PXD043880 source feature maps to multiple regulator accessions: ",
      paste(
        target_by_feature$SourceFeature[
          target_by_feature$TargetHitCount > 1
        ],
        collapse = ", "
      )
    )
  }
  feature_targets <- target_by_feature$TargetAccession[
    match(source_features, target_by_feature$SourceFeature)
  ]
  feature_targets[!nzchar(feature_targets)] <- NA_character_

  mapping |>
    mutate(
      SourceFile = relative_path(path),
      RetainedDonorCount = length(donor_ids),
      TargetRegulatorMatch = !is.na(TargetAccession)
    ) |>
    write.csv(
      file.path(
        table_dir,
        "kla_regulator_whole_proteome_hippocampus_id_mapping_audit.csv"
      ),
      row.names = FALSE,
      na = ""
    )

  add_total_quant(
    quant_data,
    "PXD050470",
    "human hippocampus",
    paste0("PXD043880_SYMBOL_FEATURE:", source_features),
    feature_targets,
    sample_labels,
    sample_labels,
    paste0(
      "PXD043880 CA1 LFQ intensity after reviewed UniProt ",
      "symbol conversion"
    ),
    path
  )
  update_audit(
    "PXD050470",
    "human hippocampus",
    relative_path(path),
    paste0(
      "PXD043880 CA1 LFQ intensity after reviewed UniProt ",
      "symbol conversion"
    ),
    paste0(
      "独立正常CA1海马全蛋白参照；74名供体；",
      "symbol先转换为人源reviewed UniProt BaseAccession"
    )
  )
}

# Read the ordinary-proteome reference selected for each lactylome row.
# The lactylome PXD defines the row; the reference PXD supplies the signal.
add_procan_average <- function(path, pxd, sample_group, project_identifier) {
  if (!file.exists(path)) return(invisible(NULL))
  data <- data.table::fread(path, data.table = FALSE, showProgress = FALSE)
  id_col <- names(data)[1]
  hit <- data[grepl(project_identifier, data[[id_col]], fixed = TRUE), , drop = FALSE]
  if (!nrow(hit)) return(invisible(NULL))
  feature_columns <- setdiff(names(hit), id_col)
  values <- as.numeric(hit[1, feature_columns, drop = TRUE])
  quant_data <- data.frame(
    FeatureID = feature_columns,
    ReferenceSignal = values,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  add_total_quant(
    quant_data,
    pxd,
    sample_group,
    feature_columns,
    match_target_accession(feature_columns),
    "ReferenceSignal",
    paste0("reference_", project_identifier),
    "ordinary proteome log2 protein abundance; rank-preserving",
    path,
    TRUE
  )
}

add_pxd073311_huvec_reference <- function(path, pxd, sample_group) {
  if (!file.exists(path)) return(invisible(NULL))
  data <- read.delim(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = ""
  )
  sample_columns <- grep("A0h_[123]\\.raw$", names(data), value = TRUE)
  if (!length(sample_columns) || !"Protein.Group" %in% names(data)) {
    return(invisible(NULL))
  }
  add_total_quant(
    data,
    pxd,
    sample_group,
    accession_feature(data$Protein.Group),
    match_target_accession(data$Protein.Group),
    sample_columns,
    sub("\\.raw$", "", basename(sample_columns)),
    "PXD073311 ordinary whole-proteome A0h baseline PG matrix",
    path,
    FALSE
  )
}

add_procan_hek293t_counts <- function(path, pxd, sample_group) {
  if (!file.exists(path)) return(invisible(NULL))
  counts <- data.table::fread(path, data.table = FALSE, showProgress = FALSE)
  mapping <- data.table::fread(
    file.path(
      project_root,
      "data/PXD030304/search_results",
      "ProCan-DepMapSanger_mapping_file_replicates.txt"
    ),
    data.table = FALSE,
    showProgress = FALSE
  )
  runs <- mapping$Automatic_MS_filename[
    mapping$Cell_line == "Control_HEK293T_lys"
  ]
  runs <- intersect(runs, counts[[1]])
  if (!length(runs)) return(invisible(NULL))
  counts <- counts[counts[[1]] %in% runs, , drop = FALSE]
  feature_columns <- names(counts)[-1]
  quant_data <- as.data.frame(
    t(as.matrix(counts[, feature_columns, drop = FALSE])),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  names(quant_data) <- paste0("HEK293T_", runs)
  add_total_quant(
    quant_data,
    pxd,
    sample_group,
    feature_columns,
    match_target_accession(feature_columns),
    names(quant_data),
    names(quant_data),
    "ordinary proteome peptide-count abundance proxy; no PTM enrichment",
    path
  )
}

add_reference_maxquant <- function(
  path,
  pxd,
  sample_group,
  sample_columns = "Intensity",
  sample_label = "reference"
) {
  if (!file.exists(path)) return(invisible(NULL))
  data <- read.delim(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = ""
  )
  accession_col <- intersect(
    c("Majority protein IDs", "Protein IDs", "Protein IDs"),
    names(data)
  )
  if (!length(accession_col)) return(invisible(NULL))
  accession_col <- accession_col[[1]]
  keep <- rep(TRUE, nrow(data))
  for (flag in c("Reverse", "Potential contaminant", "Only identified by site")) {
    if (flag %in% names(data)) keep <- keep & data[[flag]] != "+"
  }
  data <- data[keep, , drop = FALSE]
  for (token in sample_columns) {
    candidates <- c(
      token,
      paste("Intensity", token),
      paste("LFQ intensity", token),
      paste("iBAQ", token)
    )
    column <- candidates[candidates %in% names(data)]
    if (!length(column)) next
    label <- if (length(sample_label) == length(sample_columns)) {
      sample_label[match(token, sample_columns)]
    } else {
      sample_label
    }
    add_total_quant(
      data,
      pxd,
      sample_group,
      accession_feature(data[[accession_col]]),
      match_target_accession(data[[accession_col]]),
      column[[1]],
      label,
      "ordinary MaxQuant protein intensity; no PTM enrichment",
      path
    )
  }
}

add_mcf10a_reference <- function(path, pxd, sample_group) {
  if (!file.exists(path)) return(invisible(NULL))
  evidence <- read.delim(
    unz(path, "evidence.txt"),
    check.names = FALSE,
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = ""
  )
  keep <- grepl("_10A_", evidence$`Raw file`, fixed = TRUE) &
    safe_numeric(evidence$Intensity) > 0 &
    (is.na(evidence$Reverse) | evidence$Reverse != "+") &
    (is.na(evidence$Contaminant) | evidence$Contaminant != "+")
  evidence <- evidence[keep, , drop = FALSE]
  if (!nrow(evidence)) return(invisible(NULL))
  evidence$FeatureID <- accession_feature(evidence$`Leading Razor Protein`)
  evidence <- evidence |>
    filter(!is.na(FeatureID), nzchar(FeatureID)) |>
    group_by(FeatureID, `Raw file`) |>
    summarise(Signal = sum(safe_numeric(Intensity), na.rm = TRUE), .groups = "drop") |>
    pivot_wider(names_from = `Raw file`, values_from = Signal, values_fill = 0)
  feature_columns <- setdiff(names(evidence), "FeatureID")
  target_values <- match_target_accession(
    sub("^ACC:", "", evidence$FeatureID)
  )
  add_total_quant(
    evidence,
    pxd,
    sample_group,
    evidence$FeatureID,
    target_values,
    feature_columns,
    feature_columns,
    "ordinary MaxQuant evidence intensity; no PTM enrichment",
    path
  )
}

add_maxquant_zip_proteome <- function(
  zip_path,
  pxd,
  sample_group,
  sample_columns,
  sample_labels = sample_columns,
  measurement = "ordinary MaxQuant proteinGroups; no PTM enrichment"
) {
  if (!file.exists(zip_path)) return(invisible(NULL))
  data <- read.delim(
    unz(zip_path, "proteinGroups.txt"),
    check.names = FALSE,
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = ""
  )
  accession_col <- intersect(
    c("Majority protein IDs", "Protein IDs"),
    names(data)
  )
  if (!length(accession_col)) return(invisible(NULL))
  accession_col <- accession_col[[1]]
  keep <- rep(TRUE, nrow(data))
  for (flag in c("Reverse", "Potential contaminant", "Only identified by site")) {
    if (flag %in% names(data)) keep <- keep & data[[flag]] != "+"
  }
  if ("id" %in% names(data)) keep <- keep & !is.na(data$id)
  data <- data[keep, , drop = FALSE]
  sample_columns <- intersect(sample_columns, names(data))
  if (!length(sample_columns)) return(invisible(NULL))
  if (length(sample_labels) != length(sample_columns)) {
    stop("sample_labels and sample_columns length mismatch for ", zip_path)
  }
  for (index in seq_along(sample_columns)) {
    add_total_quant(
      data,
      pxd,
      sample_group,
      accession_feature(data[[accession_col]]),
      match_target_accession(data[[accession_col]]),
      sample_columns[[index]],
      sample_labels[[index]],
      measurement,
      zip_path
    )
  }
}

add_pxd069969_proteome <- function(path, pxd, sample_group, sample_names) {
  if (!file.exists(path)) return(invisible(NULL))
  data <- read_excel(path, sheet = "Annotation_Combine")
  columns <- paste0("LFQ intensity ", sample_names)
  columns <- intersect(columns, names(data))
  if (!length(columns) || !"Protein accession" %in% names(data)) {
    return(invisible(NULL))
  }
  add_total_quant(
    data,
    pxd,
    sample_group,
    accession_feature(data$`Protein accession`),
    match_target_accession(data$`Protein accession`),
    columns,
    sub("^LFQ intensity ", "", columns),
    "ordinary protein LFQ; no PTM enrichment",
    path
  )
}

add_pxd055025_proteome <- function(path, pxd, sample_group) {
  if (!file.exists(path)) return(invisible(NULL))
  data <- read_excel(path, sheet = "B VS A")
  columns <- paste0("Abundances B", 1:3)
  columns <- intersect(columns, names(data))
  if (!length(columns) || !"Accession" %in% names(data)) {
    return(invisible(NULL))
  }
  add_total_quant(
    data,
    pxd,
    sample_group,
    accession_feature(data$Accession),
    match_target_accession(data$Accession),
    columns,
    paste0("EOPE_", seq_along(columns)),
    "ordinary TMT proteome; early-onset preeclampsia placenta",
    path
  )
}

add_pxd065775_proteome <- function(path, pxd, sample_group, sheet_name) {
  if (!file.exists(path)) return(invisible(NULL))
  data <- read_excel(path, sheet = sheet_name)
  sample_columns <- c(
    "Non-rec1", "Non-rec2", "Non-rec3", "Non-rec4",
    "Rec1", "Rec2", "Rec3", "Rec4"
  )
  sample_columns <- intersect(sample_columns, names(data))
  if (!length(sample_columns) || !"Accession" %in% names(data)) {
    return(invisible(NULL))
  }
  add_total_quant(
    data,
    pxd,
    sample_group,
    accession_feature(data$Accession),
    match_target_accession(data$Accession),
    sample_columns,
    paste0(sheet_name, "_", sample_columns),
    "ordinary iTRAQ clinical liver proteome",
    path
  )
}

add_pxd059985_ac16 <- function(path, pxd, sample_group) {
  if (!file.exists(path)) return(invisible(NULL))
  data <- read.delim(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = ""
  )
  quantity_columns <- grep(
    "\\.PG\\.Quantity$",
    names(data),
    value = TRUE
  )
  quantity_columns <- quantity_columns[
    grepl("AC16|V8", quantity_columns)
  ]
  if (!length(quantity_columns) || !"PG.ProteinGroups" %in% names(data)) {
    return(invisible(NULL))
  }
  labels <- sub("^\\[[0-9]+\\] ", "", quantity_columns)
  labels <- sub("\\.htrms\\.PG\\.Quantity$", "", labels)
  add_total_quant(
    data,
    pxd,
    sample_group,
    accession_feature(data$PG.ProteinGroups),
    match_target_accession(data$PG.ProteinGroups),
    quantity_columns,
    labels,
    "ordinary Spectronaut DIA protein quantity; AC16 study",
    path
  )
}

reference_rows <- pairing |>
  inner_join(
    sample_catalog |>
      select(PXD, SampleGroup),
    by = c("LactylomePXD" = "PXD", "SampleGroup")
  )
for (i in seq_len(nrow(reference_rows))) {
  ref <- reference_rows[i, ]
  pxd <- ref$LactylomePXD
  group <- ref$SampleGroup
  if (
    is.na(ref$ReferencePXD) ||
      !nzchar(ref$ReferencePXD) ||
      ref$IncludeInStrictReferenceAnalysis %in%
        c(FALSE, "FALSE", "False", 0, "0")
  ) {
    next
  }
  path <- file.path(project_root, ref$ReferenceEvidenceLocator)
  if (ref$ReferencePXD == "PXD028488" && group == "TALL-104") {
    add_peaks_proteins(path, pxd, group)
  } else if (ref$ReferencePXD == "PXD028737" && group == "HMC3") {
    add_maxquant_proteome(
      path,
      pxd,
      group,
      list(c("H0", "H0"), c("H24", "H24")),
      "same-study ordinary MaxQuant protein LFQ; HMC3 H0/H24"
    )
  } else if (ref$ReferencePXD == "PXD030304" && group == "HEK293T") {
    add_procan_hek293t_counts(path, pxd, group)
  } else if (ref$ReferencePXD == "PXD030304") {
    project_id <- c(
      "MCF7" = "SIDM00148;MCF7",
      "HCT116" = "SIDM00783;HCT-116",
      "TALL-104" = "SIDM00370;TALL-1",
      "A549" = "SIDM00903;A549",
      "MCF10A" = "not-used",
      "MCF7" = "SIDM00148;MCF7",
      "MDA-MB-468" = "SIDM00628;MDA-MB-468",
      "T-47D" = "SIDM00097;T47D",
      "RKO WT and GSK3B KO" = "SIDM01090;RKO",
      "PC-3M" = "SIDM00088;PC-3"
    )[[group]]
    if (!is.null(project_id)) {
      add_procan_average(path, pxd, group, project_id)
    }
  } else if (ref$ReferencePXD == "PXD043880") {
    add_pxd043880_hippocampus(path)
  } else if (ref$ReferencePXD == "PXD002400") {
    add_mcf10a_reference(path, pxd, group)
  } else if (ref$ReferencePXD == "PXD010154") {
    add_reference_maxquant(path, pxd, group, "Intensity", "healthy_reference")
  } else if (ref$ReferencePXD == "PXD033146") {
    add_reference_maxquant(
      file.path(
        project_root,
        "data/PXD033146/search_results/extracted_pairing",
        "search_result-HA119TQ/proteinGroups.txt"
      ),
      pxd,
      group,
      paste0("Reporter intensity corrected ", 1:6),
      paste0("channel_", 1:6)
    )
  } else if (ref$ReferencePXD == "PXD046800") {
    add_pd_proteome(
      path,
      pxd,
      list(list(group, ifelse(group == "hypertrophic scar", "HSP", "NSP")))
    )
  } else if (ref$ReferencePXD == "PXD050147") {
    add_reference_maxquant(
      path,
      pxd,
      group,
      c(
        "SIRT1KO_pro_rep1", "SIRT1KO_pro_rep2", "SIRT1KO_pro_rep3",
        "SIRT3KO_pro_rep1", "SIRT3KO_pro_rep2", "SIRT3KO_pro_rep3",
        "WT_pro_rep1", "WT_pro_rep2", "WT_pro_rep3"
      ),
      c(
        "SIRT1KO_pro_rep1", "SIRT1KO_pro_rep2", "SIRT1KO_pro_rep3",
        "SIRT3KO_pro_rep1", "SIRT3KO_pro_rep2", "SIRT3KO_pro_rep3",
        "WT_pro_rep1", "WT_pro_rep2", "WT_pro_rep3"
      )
    )
  } else if (ref$ReferencePXD == "PXD055230") {
    add_spectronaut_report(
      if (file.exists(path)) path else file.path(
        project_root, "data/PXD055230/search_results/WP_HSV1_DIA.tsv"
      ),
      pxd,
      group,
      "ordinary Spectronaut PG.Quantity"
    )
  } else if (ref$ReferencePXD == "PXD057709") {
    add_spectronaut_report(
      if (file.exists(path)) path else file.path(
        project_root, "data/PXD057709/search_results/WP_report.tsv"
      ),
      pxd,
      group,
      "ordinary Spectronaut PG.Quantity"
    )
  } else if (ref$ReferencePXD == "PXD072220") {
    add_spectronaut_standard_report(
      path,
      pxd,
      group,
      "ordinary Spectronaut PG.Log2Quantity; untreated HK-2 reference",
      "amostra(1|3|4)\\.raw\\.PG\\.Log2Quantity$",
      2,
      TRUE
    )
  } else if (ref$ReferencePXD == "PXD066517") {
    add_spectronaut_standard_report(
      path,
      pxd,
      group,
      "ordinary Spectronaut PG.Quantity; human sperm reference",
      "PG.Quantity"
    )
  } else if (ref$ReferencePXD == "PXD073311") {
    add_pxd073311_huvec_reference(path, pxd, group)
  } else if (ref$ReferencePXD == "PXD062720") {
    add_reference_maxquant(
      file.path(
        project_root,
        "data/PXD062720/search_results/extracted_pairing/txt/proteinGroups.txt"
      ),
      pxd,
      group,
      c("A_1", "A_3"),
      c("A_1", "A_3")
    )
  } else if (ref$ReferencePXD == "PXD066054") {
    add_spectronaut_matrix(
      path,
      pxd,
      group,
      "ordinary Spectronaut PG.Quantity",
      ifelse(group == "BPH", "NAT", "PCa")
    )
  } else if (ref$ReferencePXD == "PXD066351") {
    add_spectronaut_matrix(
      path,
      pxd,
      group,
      "ordinary Spectronaut PG.Quantity"
    )
  } else if (ref$ReferencePXD == "PXD022005" && group == "PC-3M") {
    add_maxquant_zip_proteome(
      path,
      pxd,
      group,
      "Intensity H",
      "PC-3M_SILAC_proteome",
      "ordinary PC-3M heavy-channel SILAC proteome protein intensity"
    )
  } else if (
    ref$ReferencePXD == "PXD055025" &&
      group == "severe preeclampsia placenta"
  ) {
    add_pxd055025_proteome(path, pxd, group)
  } else if (
    ref$ReferencePXD == "PXD069969" &&
      group %in% c("glioblastoma stem cells", "neural stem cells")
  ) {
    samples <- if (group == "glioblastoma stem cells") {
      c("G2907", "G3028", "G3264", "GSC23", "MES28", "RKI")
    } else {
      c("ENSA", "HMP1")
    }
    add_pxd069969_proteome(path, pxd, group, samples)
  } else if (
    ref$ReferencePXD == "PXD059985" &&
      group == "AC16 control and hypoxia"
  ) {
    add_pxd059985_ac16(path, pxd, group)
  } else if (
    ref$ReferencePXD == "PXD065775" &&
      group %in% c("HCC", "adjacent liver")
  ) {
    add_pxd065775_proteome(
      path,
      pxd,
      group,
      ifelse(group == "HCC", "CISs", "ANTs")
    )
  }
}

if (!length(quant_parts)) stop("No whole-proteome quantitative data extracted")

quant_features <- bind_rows(quant_parts) |>
  group_by(
    PXD, SampleGroup, QuantSample, CanonicalFeature,
    TargetAccession, TargetGene, Measurement, SignalIsLog2
  ) |>
  summarise(
    Signal = sum(Signal, na.rm = TRUE),
    SourceFile = paste(sort(unique(SourceFile)), collapse = ";"),
    .groups = "drop"
  ) |>
  group_by(PXD, SampleGroup, QuantSample) |>
  mutate(
    Log2Signal = ifelse(
      SignalIsLog2,
      Signal,
      log2(Signal + 1)
    ),
    WholeProteomePercentile = {
      feature_count <- n()
      if (feature_count <= 1) {
        100
      } else {
        100 * (rank(Log2Signal, ties.method = "average") - 1) /
          (feature_count - 1)
      }
    }
  ) |>
  ungroup()

# This is the only ranking universe for the ordinary-proteome heatmap:
# all positive, finite whole-proteome features from the same quantitative sample.
# Kla-enriched tables are deliberately never joined into this table.
target_quant <- quant_features |>
  filter(
    !is.na(TargetAccession),
    TargetAccession %in% target_accessions
  ) |>
  group_by(
    PXD, SampleGroup, QuantSample,
    TargetAccession, TargetGene
  ) |>
  summarise(
    Signal = sum(Signal),
    Log2Signal = log2(Signal + 1),
    WholeProteomePercentile = max(WholeProteomePercentile),
    Measurement = paste(sort(unique(Measurement)), collapse = ";"),
    SourceFile = paste(sort(unique(SourceFile)), collapse = ";"),
    .groups = "drop"
  )

sample_registry <- quant_features |>
  distinct(PXD, SampleGroup, QuantSample)
target_sample_grid <- sample_registry |>
  crossing(
    regulators |>
      distinct(GeneSymbol, BaseAccession) |>
      rename(RegulatorBaseAccession = BaseAccession)
  ) |>
  left_join(
    target_quant,
    by = c(
      "PXD", "SampleGroup", "QuantSample",
      "RegulatorBaseAccession" = "TargetAccession"
    )
  ) |>
  mutate(
    GeneSymbol = unname(accession_to_gene[RegulatorBaseAccession]),
    IdentityMatchMode = "UniProt_BaseAccession_only",
    Detected = !is.na(Signal) & Signal > 0,
    Signal = replace_na(Signal, 0),
    Log2Signal = replace_na(Log2Signal, 0),
    WholeProteomePercentile = replace_na(WholeProteomePercentile, 0)
  )

group_summary <- target_sample_grid |>
  group_by(PXD, SampleGroup, RegulatorBaseAccession, GeneSymbol) |>
  summarise(
    WholeProteomeRelativePercentile = median(WholeProteomePercentile),
    DetectedSampleCount = sum(Detected),
    QuantSampleCount = n(),
    DetectedSampleFraction = mean(Detected),
    MedianLog2SignalDetected = ifelse(
      any(Detected),
      median(Log2Signal[Detected]),
      NA_real_
    ),
    Measurement = paste(sort(unique(na.omit(Measurement))), collapse = ";"),
    SourceFile = paste(sort(unique(na.omit(SourceFile))), collapse = ";"),
    .groups = "drop"
  )

quant_summary <- target_sample_grid |>
  group_by(PXD, SampleGroup) |>
  summarise(
    WholeProteomeSampleCount = n_distinct(QuantSample),
    WholeProteomeMappedRegulatorCount = n_distinct(
      RegulatorBaseAccession[Detected]
    ),
    .groups = "drop"
  )
whole_feature_summary <- quant_features |>
  group_by(PXD, SampleGroup) |>
  summarise(
    WholeProteomeFeatureCount = n_distinct(CanonicalFeature),
    .groups = "drop"
  )
quant_summary <- quant_summary |>
  left_join(whole_feature_summary, by = c("PXD", "SampleGroup"))

quant_audit <- quant_audit |>
  left_join(quant_summary, by = c("PXD", "SampleGroup")) |>
  mutate(
    WholeProteomeSampleCount = replace_na(WholeProteomeSampleCount, 0L),
    WholeProteomeFeatureCount = replace_na(
      WholeProteomeFeatureCount, 0L
    ),
    WholeProteomeMappedRegulatorCount = replace_na(
      WholeProteomeMappedRegulatorCount, 0L
    ),
    WholeProteomeQuantAvailable =
      WholeProteomeSampleCount > 0,
    WholeProteomeStatus = ifelse(
      WholeProteomeQuantAvailable, "available", "unavailable"
    )
  )

write.csv(
  quant_audit,
  file.path(
    table_dir,
    "kla_regulator_whole_proteome_intensity_availability_audit.csv"
  ),
  row.names = FALSE,
  na = ""
)
write.csv(
  quant_audit |>
    transmute(
      PXD,
      SampleGroup,
      ReferencePXD,
      WholeProteomeSource,
      EnsemblMappingApplied = ReferencePXD == "PXD010154",
      EnsemblMappingFile = relative_path(ensembl_mapping_path),
      WholeProteomeFeatureCount,
      WholeProteomeMappedRegulatorCount
    ),
  file.path(
    table_dir,
    "kla_regulator_whole_proteome_ensembl_mapping_audit.csv"
  ),
  row.names = FALSE,
  na = ""
)
write.csv(
  quant_features,
  file.path(
    table_dir,
    "kla_regulator_whole_proteome_sample_level_long.csv"
  ),
  row.names = FALSE,
  na = ""
)
write.csv(
  target_sample_grid,
  file.path(
    table_dir,
    "kla_regulator_whole_proteome_regulator_sample_level_long.csv"
  ),
  row.names = FALSE,
  na = ""
)
write.csv(
  group_summary,
  file.path(
    table_dir,
    "kla_regulator_whole_proteome_normalized_long.csv"
  ),
  row.names = FALSE,
  na = ""
)

algorithm_audit <- data.frame(
  Field = c(
    "AlgorithmVersion",
    "InputData",
    "SignalFilter",
    "Transform",
    "RankingUniverse",
    "PercentileFormula",
    "TieHandling",
    "SampleGroupAggregation",
    "UndetectedRegulatorHandling",
    "NoWholeProteomeHandling",
    "IdentityKey",
    "EnsemblProteinMapping",
    "GeneSymbolRole",
    "KlaSignalSubstitution"
  ),
  Value = c(
    "whole_proteome_regulator_rank_v4_exact_reference",
    "ordinary whole-proteome quantitative files or selected normal whole-proteome references",
    "finite Signal > 0; reverse/contaminant/only-identified-by-site rows excluded when available",
    "log2(Signal + 1)",
    "all retained whole-proteome features within each QuantSample",
    "100 * (average_rank - 1) / (n_features - 1)",
    "average rank for ties",
    "median percentile across QuantSample replicates/conditions within PXD and SampleGroup",
    "0 percentile for a regulator not detected in an otherwise usable whole-proteome sample",
    "sample group excluded from the ordinary-proteome heatmap when no exact per-protein quantitative reference exists",
    "UniProt human reviewed BaseAccession after isoform suffix removal; Ensembl protein IDs are converted to UniProt by the project mapping table",
    "ENSP(_RNA/version) -> UniProt BaseAccession via ensembl_protein_to_uniprot_biomart.tsv; GeneSymbol is never used as a fallback",
    "display/audit only; never used for matching or aggregation",
    "never; Kla-enriched intensity is not used as a fallback"
  ),
  stringsAsFactors = FALSE
)
write.csv(
  algorithm_audit,
  file.path(
    table_dir,
    "kla_regulator_whole_proteome_algorithm_audit.csv"
  ),
  row.names = FALSE,
  na = ""
)

available_sample_groups <- paste(
  quant_audit$PXD[quant_audit$WholeProteomeQuantAvailable],
  quant_audit$SampleGroup[quant_audit$WholeProteomeQuantAvailable],
  sep = "__"
)
plot_data <- sample_catalog |>
  select(PXD, SampleGroup, SampleGroupID, RowLabel, RowOrder) |>
  crossing(
    regulators |>
      select(Role, GeneSymbol, BaseAccession, RoleEntryOrder) |>
      rename(RegulatorBaseAccession = BaseAccession)
  ) |>
  left_join(
    group_summary,
    by = c("PXD", "SampleGroup", "GeneSymbol", "RegulatorBaseAccession")
  ) |>
  mutate(
    QuantificationAvailable = paste(PXD, SampleGroup, sep = "__") %in%
      available_sample_groups,
    WholeProteomeRelativePercentile = ifelse(
      QuantificationAvailable,
      replace_na(WholeProteomeRelativePercentile, 0),
      NA_real_
    ),
    Role = factor(
      Role,
      levels = c("Writer", "Eraser", "Writer-Eraser", "Reader")
    ),
    GeneSymbol = factor(
      GeneSymbol,
      levels = unique(regulators$GeneSymbol[order(regulators$RoleEntryOrder)])
    ),
    RowLabel = factor(RowLabel, levels = rev(sample_catalog$RowLabel))
  )

plot_font <- "Arial Unicode MS"
main_plot <- ggplot(
  plot_data,
  aes(GeneSymbol, RowLabel, fill = WholeProteomeRelativePercentile)
) +
  geom_tile(color = "white", linewidth = 0.22) +
  geom_text(
    data = plot_data |> filter(!QuantificationAvailable),
    aes(label = "?"),
    color = "#4D4D4D",
    size = 2.3,
    fontface = "bold",
    family = plot_font
  ) +
  facet_grid(. ~ Role, scales = "free_x", space = "free_x") +
  scale_fill_gradientn(
    colours = c("#FFFFFF", "#FFF3E0", "#FDBB84", "#FC8D59", "#B2182B"),
    values = scales::rescale(c(0, 20, 50, 80, 100)),
    limits = c(0, 100),
    na.value = "#D9D9D9",
    name = "全蛋白组信号\n百分位"
  ) +
  labs(
    title = "乳酸化调控蛋白在全蛋白组中的相对信号",
    subtitle = paste(
      "仅纳入生物材料匹配且具有逐蛋白强度的普通全蛋白定量文件，不使用Kla富集信号；",
      "颜色由白色向暖色递增。"
    ),
    x = NULL,
    y = NULL,
    caption = paste0(
      "身份按人源UniProt reviewed BaseAccession匹配；",
      "全蛋白组百分位只在各自样本内计算，不是Log FC。"
    )
  ) +
  theme_minimal(base_size = 8.5, base_family = plot_font) +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold", size = 9),
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

ggsave(
  file.path(
    figure_dir,
    "kla_regulator_whole_proteome_relative_intensity_heatmap.png"
  ),
  main_plot,
  width = 15.5,
  height = 11.5,
  dpi = 320,
  bg = "white"
)
ggsave(
  file.path(
    figure_dir,
    "kla_regulator_whole_proteome_relative_intensity_heatmap.pdf"
  ),
  main_plot,
  width = 15.5,
  height = 11.5,
  device = cairo_pdf,
  bg = "white"
)

message(
  "Whole-proteome regulator heatmap: ",
  sum(quant_audit$WholeProteomeQuantAvailable),
  "/",
  nrow(quant_audit),
  " sample groups with usable total-proteome intensity."
)
