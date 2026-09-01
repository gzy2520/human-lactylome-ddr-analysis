#!/usr/bin/env Rscript

# Build candidate-only, sample-resolved Kla–DDR inputs from processed public
# source files. The resulting small CSVs are committed so the companion figure
# can be rendered from the repository without the local source cache.

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
  stop("Usage: Rscript R/candidate/prepare_sample_resolved_inputs.R <project-root> <source-data-root>", call. = FALSE)
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
  values <- sub("^(?:REV__|CON__)+", "", values, perl = TRUE)
  values <- sub("^(?:sp|tr)\\|", "", values, perl = TRUE)
  values <- sub("\\|.*$", "", values)
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
  values <- base_accession(unlist(strsplit(as.character(value), "[;,]"), use.names = FALSE))
  values <- values[is_uniprot(values)]
  if (length(values)) values[[1L]] else ""
}

is_true_flag <- function(values) {
  values <- tolower(trimws(as.character(values)))
  values[is.na(values)] <- ""
  values %in% c("+", "1", "true", "yes", "y")
}

is_detected <- function(values) {
  values <- trimws(as.character(values))
  !is.na(values) & nzchar(values) & values != "Not Found"
}

positive_numeric <- function(values) {
  values <- suppressWarnings(as.numeric(values))
  is.finite(values) & values > 0
}

records_from_accession_rows <- function(rows, sample_key, accession_column) {
  if (!nrow(rows)) {
    return(data.table(SampleKey = character(), BaseAccession = character()))
  }
  accession_lists <- lapply(rows[[accession_column]], split_accessions)
  lengths <- lengths(accession_lists)
  if (!sum(lengths)) {
    return(data.table(SampleKey = character(), BaseAccession = character()))
  }
  unique(data.table(
    SampleKey = rep(sample_key, sum(lengths)),
    BaseAccession = unlist(accession_lists, use.names = FALSE)
  ))
}

records_from_first_accession_rows <- function(rows, sample_key, accession_column) {
  if (!nrow(rows)) {
    return(data.table(SampleKey = character(), BaseAccession = character()))
  }
  accessions <- vapply(rows[[accession_column]], first_accession, character(1))
  unique(data.table(
    SampleKey = sample_key,
    BaseAccession = accessions[nzchar(accessions)]
  ))
}

pd_sample_column <- function(columns, sample_id) {
  found_columns <- columns[startsWith(columns, "Found in Sample:")]
  hits <- found_columns[
    endsWith(found_columns, paste0(", ", sample_id)) |
      endsWith(found_columns, paste0(", ", sample_id, "_1"))
  ]
  stop_if(
    length(hits) == 1L,
    paste0("Expected one Proteome Discoverer sample column for ", sample_id, "; found ", length(hits), ".")
  )
  hits[[1L]]
}

read_pd_by_sample <- function(path, registry_rows, accession_column, require_lactylation = FALSE) {
  header <- names(fread(path, nrows = 0L, check.names = FALSE, showProgress = FALSE))
  sample_columns <- vapply(registry_rows$SampleID, function(sample_id) {
    pd_sample_column(header, sample_id)
  }, character(1))
  selected <- unique(c(accession_column, sample_columns, if (require_lactylation) "Modifications" else character()))
  data <- fread(path, select = selected, check.names = FALSE, showProgress = FALSE)
  if (require_lactylation) {
    data <- data[grepl("Lacty|Lactyl|La \\(K\\)", Modifications, ignore.case = TRUE)]
  }
  rbindlist(lapply(seq_len(nrow(registry_rows)), function(index) {
    sample_column <- sample_columns[[index]]
    rows <- data[is_detected(data[[sample_column]])]
    records_from_accession_rows(rows, registry_rows$SampleKey[[index]], accession_column)
  }))
}

read_intensity_table_by_sample <- function(path, skip, registry_rows, accession_column, position_column = NULL) {
  data <- as.data.table(read_excel(path, sheet = "Sheet1", skip = skip))
  expected_columns <- c(accession_column, paste0("Intensity_", registry_rows$SampleID))
  stop_if(
    all(expected_columns %in% names(data)),
    paste0("Unexpected intensity-table layout in ", path)
  )
  if (!is.null(position_column)) {
    stop_if(position_column %in% names(data), paste0("Missing ", position_column, " in ", path))
    data <- data[!is.na(get(position_column))]
  }
  rbindlist(lapply(seq_len(nrow(registry_rows)), function(index) {
    intensity_column <- paste0("Intensity_", registry_rows$SampleID[[index]])
    rows <- data[positive_numeric(data[[intensity_column]])]
    records_from_accession_rows(rows, registry_rows$SampleKey[[index]], accession_column)
  }))
}

read_pxd066054_kla <- function(path, registry_rows) {
  columns <- c(
    "R.Condition", "PTM.ProteinId", "PTM.ModificationTitle", "PTM.SiteAA", "PTM.SiteProbability"
  )
  data <- fread(path, select = columns, check.names = FALSE, showProgress = FALSE)
  rbindlist(lapply(seq_len(nrow(registry_rows)), function(index) {
    sample_id <- registry_rows$SampleID[[index]]
    rows <- data[
      R.Condition == sample_id &
        PTM.ModificationTitle == "L-Lac(K)" &
        PTM.SiteAA == "K" &
        positive_numeric(PTM.SiteProbability)
    ]
    records_from_accession_rows(rows, registry_rows$SampleKey[[index]], "PTM.ProteinId")
  }))
}

read_pxd066054_proteome <- function(path, registry_rows) {
  header <- names(fread(path, nrows = 0L, check.names = FALSE, showProgress = FALSE))
  quantity_columns <- vapply(registry_rows$SampleID, function(sample_id) {
    hits <- header[grepl(paste0("_", sample_id, "_"), header, fixed = TRUE)]
    stop_if(length(hits) == 1L, paste0("Expected one PXD066054 protein quantity column for ", sample_id, "."))
    hits[[1L]]
  }, character(1))
  data <- fread(
    path,
    select = unique(c("PG.ProteinGroups", quantity_columns)),
    check.names = FALSE,
    showProgress = FALSE
  )
  rbindlist(lapply(seq_len(nrow(registry_rows)), function(index) {
    rows <- data[positive_numeric(data[[quantity_columns[[index]]]])]
    records_from_accession_rows(rows, registry_rows$SampleKey[[index]], "PG.ProteinGroups")
  }))
}

split_tokens <- function(value) {
  tokens <- trimws(unlist(strsplit(as.character(value), ";"), use.names = FALSE))
  tokens[!is.na(tokens) & nzchar(tokens)]
}

expand_site_tokens <- function(site_ids, sample_keys = NULL, accessions = NULL) {
  tokens <- lapply(site_ids, split_tokens)
  sizes <- lengths(tokens)
  if (!sum(sizes)) {
    return(data.table(SiteID = character(), SampleKey = character(), BaseAccession = character()))
  }
  result <- data.table(SiteID = unlist(tokens, use.names = FALSE))
  if (!is.null(sample_keys)) result[, SampleKey := rep(sample_keys, sizes)]
  if (!is.null(accessions)) result[, BaseAccession := rep(accessions, sizes)]
  result
}

read_pxd078013_kla <- function(evidence_path, proteins_path, registry_rows) {
  sample_to_key <- setNames(registry_rows$SampleKey, registry_rows$SampleID)
  evidence <- fread(
    evidence_path,
    select = c("Experiment", "La (K)", "La (K) site IDs", "Reverse", "Potential contaminant"),
    check.names = FALSE,
    showProgress = FALSE
  )
  evidence[, SampleKey := unname(sample_to_key[Experiment])]
  evidence <- evidence[
    !is.na(SampleKey) &
      !is_true_flag(Reverse) &
      !is_true_flag(`Potential contaminant`) &
      positive_numeric(`La (K)`) &
      !is.na(`La (K) site IDs`) & nzchar(trimws(`La (K) site IDs`))
  ]
  evidence_sites <- unique(expand_site_tokens(evidence$`La (K) site IDs`, sample_keys = evidence$SampleKey))

  proteins <- fread(
    proteins_path,
    select = c(
      "Majority protein IDs", "Protein IDs", "La (K) site IDs", "La (K) site positions",
      "Reverse", "Potential contaminant", "Only identified by site"
    ),
    check.names = FALSE,
    showProgress = FALSE
  )
  proteins <- proteins[
    !is_true_flag(Reverse) &
      !is_true_flag(`Potential contaminant`) &
      !is_true_flag(`Only identified by site`)
  ]
  proteins[, BaseAccession := vapply(
    fifelse(!is.na(`Majority protein IDs`) & nzchar(trimws(`Majority protein IDs`)), `Majority protein IDs`, `Protein IDs`),
    first_accession,
    character(1)
  )]
  proteins[, SiteIds := lapply(`La (K) site IDs`, split_tokens)]
  proteins[, SitePositions := lapply(`La (K) site positions`, split_tokens)]
  proteins <- proteins[
    nzchar(BaseAccession) &
      lengths(SiteIds) > 0L &
      lengths(SiteIds) == lengths(SitePositions)
  ]
  protein_sites <- unique(expand_site_tokens(
    proteins$`La (K) site IDs`,
    accessions = proteins$BaseAccession
  )[, .(SiteID, BaseAccession)])
  unique(merge(protein_sites, evidence_sites, by = "SiteID", allow.cartesian = TRUE)[, .(SampleKey, BaseAccession)])
}

read_pxd078736_kla <- function(path, registry_rows) {
  header <- names(fread(path, nrows = 0L, check.names = FALSE, showProgress = FALSE))
  id_columns <- paste0("Identification type ", registry_rows$SampleID)
  required <- c("Protein", "Leading proteins", "Amino acid", "Reverse", "Potential contaminant", id_columns)
  stop_if(all(required %in% header), paste0("Unexpected PXD078736 site-table layout in ", path))
  data <- fread(path, select = required, check.names = FALSE, showProgress = FALSE)
  data <- data[
    !is_true_flag(Reverse) &
      !is_true_flag(`Potential contaminant`) &
      `Amino acid` == "K"
  ]
  data[, BaseAccession := vapply(
    fifelse(!is.na(Protein) & nzchar(trimws(Protein)), Protein, `Leading proteins`),
    first_accession,
    character(1)
  )]
  data <- data[nzchar(BaseAccession)]
  rbindlist(lapply(seq_len(nrow(registry_rows)), function(index) {
    sample_column <- id_columns[[index]]
    rows <- data[is_detected(data[[sample_column]])]
    unique(rows[, .(SampleKey = registry_rows$SampleKey[[index]], BaseAccession)])
  }))
}

annotate_ddr_membership <- function(records, registry, ddr_lookup) {
  context <- merge(
    records,
    registry[, .(SampleKey, PXD, PublicationGroup)],
    by = "SampleKey",
    all.x = TRUE,
    sort = FALSE
  )
  context <- merge(
    context,
    unique(ddr_lookup[, .(PXD, PublicationGroup, BaseAccession, IsDdr = TRUE)]),
    by = c("PXD", "PublicationGroup", "BaseAccession"),
    all.x = TRUE,
    sort = FALSE
  )
  context[is.na(IsDdr), IsDdr := FALSE]
  context
}

summarize_detection <- function(records, registry, modality) {
  records <- unique(records)
  total_counts <- records[, .(ProteinCount = uniqueN(BaseAccession)), by = SampleKey]
  ddr_counts <- records[IsDdr == TRUE, .(DdrProteinCount = uniqueN(BaseAccession)), by = SampleKey]
  output <- merge(registry, total_counts, by = "SampleKey", all.x = TRUE, sort = FALSE)
  output <- merge(output, ddr_counts, by = "SampleKey", all.x = TRUE, sort = FALSE)
  output[is.na(ProteinCount), ProteinCount := 0L]
  output[is.na(DdrProteinCount), DdrProteinCount := 0L]
  output[, DdrFraction := fifelse(ProteinCount > 0L, DdrProteinCount / ProteinCount, NA_real_)]
  setnames(output, c("ProteinCount", "DdrProteinCount", "DdrFraction"), paste0(modality, c("ProteinCount", "DdrProteinCount", "DdrFraction")))
  output
}

read_frozen_pathway_scores <- function(path) {
  pathways <- c("BER", "NER", "MMR", "FA", "HR", "AEJ", "NHEJ")
  tables <- lapply(excel_sheets(path), function(sheet) {
    data <- as.data.table(read_excel(path, sheet = sheet))
    stop_if(all(c("BaseAccession", pathways) %in% names(data)), paste0("Unexpected pathway score layout in ", sheet))
    data[, c("BaseAccession", pathways), with = FALSE]
  })
  scores <- rbindlist(tables, fill = TRUE)
  scores[, BaseAccession := base_accession(BaseAccession)]
  scores <- scores[nzchar(BaseAccession)]
  for (pathway in pathways) {
    inconsistent <- scores[, .(DistinctScoreCount = uniqueN(get(pathway))), by = BaseAccession][DistinctScoreCount > 1L]
    stop_if(!nrow(inconsistent), paste0("Frozen pathway scores disagree for ", pathway, "."))
  }
  unique(scores[, lapply(.SD, function(values) values[[1L]]), by = BaseAccession, .SDcols = pathways])
}

build_pathway_profile <- function(kla_records, sample_summary, pathway_scores) {
  pathways <- c("BER", "NER", "MMR", "FA", "HR", "AEJ", "NHEJ")
  sample_ddr <- unique(kla_records[IsDdr == TRUE & BaseAccession %in% pathway_scores$BaseAccession])
  joined <- merge(sample_ddr, pathway_scores, by = "BaseAccession", all = FALSE)
  long <- melt(
    joined,
    id.vars = c("SampleKey", "BaseAccession"),
    measure.vars = pathways,
    variable.name = "Pathway",
    value.name = "SignedState"
  )
  counts <- long[, .(
    PositiveProteinCount = uniqueN(BaseAccession[SignedState > 0]),
    NegativeProteinCount = uniqueN(BaseAccession[SignedState < 0]),
    AnyPathwayProteinCount = uniqueN(BaseAccession[SignedState != 0])
  ), by = .(SampleKey, Pathway)]
  grid <- CJ(SampleKey = sample_summary$SampleKey, Pathway = pathways, unique = TRUE)
  output <- merge(grid, counts, by = c("SampleKey", "Pathway"), all.x = TRUE)
  output[is.na(PositiveProteinCount), PositiveProteinCount := 0L]
  output[is.na(NegativeProteinCount), NegativeProteinCount := 0L]
  output[is.na(AnyPathwayProteinCount), AnyPathwayProteinCount := 0L]
  output <- merge(
    output,
    sample_summary[, .(SampleKey, PXD, SampleID, PublicationGroup, Category, ConditionLabel, EvidenceTier, KlaDdrProteinCount)],
    by = "SampleKey",
    all.x = TRUE,
    sort = FALSE
  )
  output[, FractionOfKlaDdrProteins := fifelse(
    KlaDdrProteinCount > 0L,
    AnyPathwayProteinCount / KlaDdrProteinCount,
    NA_real_
  )]
  setcolorder(output, c(
    "SampleKey", "PXD", "SampleID", "PublicationGroup", "Category", "ConditionLabel", "EvidenceTier",
    "Pathway", "PositiveProteinCount", "NegativeProteinCount", "AnyPathwayProteinCount",
    "KlaDdrProteinCount", "FractionOfKlaDdrProteins"
  ))
  output
}

candidate_dir <- file.path(project_root, "data", "candidate")
input_dir <- file.path(project_root, "data", "publication_input")
registry_path <- file.path(candidate_dir, "sample_resolved_source_registry.csv")
registry <- fread(registry_path, na.strings = c("", "NA"))
required_registry_columns <- c(
  "SampleKey", "PXD", "SampleID", "PublicationGroup", "Category", "ConditionLabel", "EvidenceTier", "UnitType",
  "PairingStatus", "KlaParser", "WholeProteomeParser", "KlaSourceFile", "WholeProteomeSourceFile"
)
stop_if(identical(names(registry), required_registry_columns), "Sample-resolved source registry schema changed.")
stop_if(nrow(registry) == 31L, "Sample-resolved source registry must contain 31 rows.")
stop_if(!anyDuplicated(registry$SampleKey), "Sample-resolved source registry has duplicate sample keys.")
stop_if(all(registry$EvidenceTier %in% c("A_same_source_same_sample", "B_kla_sample_only")), "Unknown source evidence tier.")
stop_if(
  registry[EvidenceTier == "B_kla_sample_only", all(is.na(WholeProteomeSourceFile))],
  "Tier B samples must not claim a paired whole-proteome source file."
)

frozen_kla_membership <- fread(file.path(input_dir, "kla_protein_membership_30.csv"), check.names = FALSE)
kla_ddr_lookup <- unique(frozen_kla_membership[
  IsDdr %in% c(TRUE, "TRUE", "True", 1, "1"),
  .(PXD, PublicationGroup = SampleGroup, BaseAccession)
])
stop_if(
  uniqueN(kla_ddr_lookup$BaseAccession) == 399L,
  "Frozen publication Kla-DDR membership must contain 399 BaseAccessions."
)

frozen_reference_membership <- fread(file.path(input_dir, "reference_protein_membership_30.csv"), check.names = FALSE)
reference_ddr_rows <- frozen_reference_membership[IsDdr %in% c(TRUE, "TRUE", "True", 1, "1")]
reference_ddr_lookup <- unique(rbindlist(lapply(seq_len(nrow(reference_ddr_rows)), function(index) {
  mapped_accessions <- split_accessions(reference_ddr_rows$MappedBaseAccessions[[index]])
  if (!length(mapped_accessions)) {
    return(data.table(PXD = character(), PublicationGroup = character(), BaseAccession = character()))
  }
  data.table(
    PXD = reference_ddr_rows$PXD[[index]],
    PublicationGroup = reference_ddr_rows$SampleGroup[[index]],
    BaseAccession = mapped_accessions
  )
})))
stop_if(nrow(reference_ddr_lookup) > 0L, "Frozen reference DDR membership is empty.")

registry_046800 <- registry[PXD == "PXD046800"]
registry_050470 <- registry[PXD == "PXD050470"]
registry_066054 <- registry[PXD == "PXD066054"]
registry_078013 <- registry[PXD == "PXD078013"]
registry_078736 <- registry[PXD == "PXD078736"]

kla_046800 <- read_pd_by_sample(
  require_file(file.path(source_root, "PXD046800/search_results/HFX2_LFQ_QB001_Lacty_PeptideGroups.txt")),
  registry_046800,
  "Master Protein Accessions",
  require_lactylation = TRUE
)
reference_046800 <- read_pd_by_sample(
  require_file(file.path(source_root, "PXD046800/search_results/HFX2_LFQ_QB002_Proteins.txt")),
  registry_046800,
  "Accession"
)

kla_050470 <- read_intensity_table_by_sample(
  require_file(file.path(source_root, "PXD050470/supplementary/prca2331-sup-0005-tables3.xlsx")),
  skip = 12L,
  registry_050470,
  accession_column = "Proteins accession",
  position_column = "Positions within proteins"
)
reference_050470 <- read_intensity_table_by_sample(
  require_file(file.path(source_root, "PXD050470/supplementary/prca2331-sup-0006-tables4.xlsx")),
  skip = 5L,
  registry_050470,
  accession_column = "Protein accession"
)

kla_066054 <- read_pxd066054_kla(
  require_file(file.path(source_root, "PXD066054/search_results/extracted/PLa/XB08700B1DPLa_-PTMSiteReport.tsv")),
  registry_066054
)
reference_066054 <- read_pxd066054_proteome(
  require_file(file.path(source_root, "PXD066054/search_results/extracted/DA/Protein_Quant.tsv")),
  registry_066054
)

kla_078013 <- read_pxd078013_kla(
  require_file(file.path(source_root, "PXD078013/search_results/evidence.txt")),
  require_file(file.path(source_root, "PXD078013/search_results/proteinGroups.txt")),
  registry_078013
)
kla_078736 <- read_pxd078736_kla(
  require_file(file.path(source_root, "PXD078736/search_results/txt/La(K)Sites.txt")),
  registry_078736
)

kla_records <- unique(rbindlist(list(kla_046800, kla_050470, kla_066054, kla_078013, kla_078736)))
reference_records <- unique(rbindlist(list(reference_046800, reference_050470, reference_066054)))
stop_if(all(is_uniprot(kla_records$BaseAccession)), "Kla records contain a non-UniProt identifier.")
stop_if(all(is_uniprot(reference_records$BaseAccession)), "Reference records contain a non-UniProt identifier.")
stop_if(setequal(unique(kla_records$SampleKey), registry$SampleKey), "A registered Kla sample has no source-derived detection record.")
stop_if(
  setequal(unique(reference_records$SampleKey), registry[EvidenceTier == "A_same_source_same_sample", SampleKey]),
  "A Tier A sample has no matched whole-proteome detection record."
)

kla_records <- annotate_ddr_membership(kla_records, registry, kla_ddr_lookup)
reference_records <- annotate_ddr_membership(
  reference_records,
  registry[EvidenceTier == "A_same_source_same_sample"],
  reference_ddr_lookup
)
kla_summary <- summarize_detection(kla_records, registry, "Kla")
reference_summary <- summarize_detection(
  reference_records,
  registry[EvidenceTier == "A_same_source_same_sample"],
  "WholeProteome"
)
paired_summary <- merge(
  kla_summary,
  reference_summary[, .(SampleKey, WholeProteomeProteinCount, WholeProteomeDdrProteinCount, WholeProteomeDdrFraction)],
  by = "SampleKey",
  all = FALSE,
  sort = FALSE
)
paired_summary[, DeltaPercentagePoints := (KlaDdrFraction - WholeProteomeDdrFraction) * 100]

frozen_groups <- fread(file.path(input_dir, "group_summary_30.csv"))
observed_kla_group <- kla_records[, .(
  SourceKlaProteinCount = uniqueN(BaseAccession),
  SourceKlaDdrProteinCount = uniqueN(BaseAccession[IsDdr == TRUE])
), by = .(PXD, PublicationGroup)]
observed_reference_group <- reference_records[, .(
  SourceWholeProteomeProteinCount = uniqueN(BaseAccession),
  SourceWholeProteomeDdrProteinCount = uniqueN(BaseAccession[IsDdr == TRUE])
), by = .(PXD, PublicationGroup)]
reconciliation <- merge(
  frozen_groups[PXD %in% registry$PXD, .(
    PXD, PublicationGroup = SampleGroup,
    FrozenKlaProteinCount = KlaProteinCount,
    FrozenKlaDdrProteinCount = KlaDdrProteinCount,
    FrozenWholeProteomeProteinCount = ReferenceProteinCount,
    FrozenWholeProteomeDdrProteinCount = ReferenceDdrProteinCount
  )],
  observed_kla_group,
  by = c("PXD", "PublicationGroup"),
  all.x = TRUE,
  sort = FALSE
)
reconciliation <- merge(reconciliation, observed_reference_group, by = c("PXD", "PublicationGroup"), all.x = TRUE, sort = FALSE)
reconciliation[, KlaProteinCountMatchesFrozen := SourceKlaProteinCount == FrozenKlaProteinCount]
reconciliation[, KlaDdrProteinCountMatchesFrozen := SourceKlaDdrProteinCount == FrozenKlaDdrProteinCount]
reconciliation[, WholeProteomeProteinCountMatchesFrozen := is.na(SourceWholeProteomeProteinCount) | SourceWholeProteomeProteinCount == FrozenWholeProteomeProteinCount]
reconciliation[, WholeProteomeDdrProteinCountMatchesFrozen := is.na(SourceWholeProteomeDdrProteinCount) | SourceWholeProteomeDdrProteinCount == FrozenWholeProteomeDdrProteinCount]
stop_if(
  all(reconciliation$KlaProteinCountMatchesFrozen & reconciliation$KlaDdrProteinCountMatchesFrozen),
  paste0(
    "Sample-derived Kla unions do not reconcile to frozen inputs: ",
    paste(
      reconciliation[
        !(KlaProteinCountMatchesFrozen & KlaDdrProteinCountMatchesFrozen),
        paste0(
          PXD, ":", PublicationGroup,
          " [source ", SourceKlaProteinCount, "/", SourceKlaDdrProteinCount,
          "; frozen ", FrozenKlaProteinCount, "/", FrozenKlaDdrProteinCount, "]"
        )
      ],
      collapse = ", "
    )
  )
)
stop_if(
  all(reconciliation$WholeProteomeProteinCountMatchesFrozen & reconciliation$WholeProteomeDdrProteinCountMatchesFrozen),
  paste0(
    "Sample-derived whole-proteome unions do not reconcile to frozen inputs: ",
    paste(
      reconciliation[
        !(WholeProteomeProteinCountMatchesFrozen & WholeProteomeDdrProteinCountMatchesFrozen),
        paste0(
          PXD, ":", PublicationGroup,
          " [source ", SourceWholeProteomeProteinCount, "/", SourceWholeProteomeDdrProteinCount,
          "; frozen ", FrozenWholeProteomeProteinCount, "/", FrozenWholeProteomeDdrProteinCount, "]"
        )
      ],
      collapse = ", "
    )
  )
)

pathway_scores <- read_frozen_pathway_scores(file.path(input_dir, "Supplementary_Table_S4_Pathway_Protein_Ranking.xlsx"))
missing_pathway_scores <- setdiff(unique(kla_records$BaseAccession[kla_records$IsDdr]), pathway_scores$BaseAccession)
stop_if(!length(missing_pathway_scores), "A source-derived Kla-DDR accession is absent from the frozen pathway-score table.")
pathway_profile <- build_pathway_profile(kla_records, kla_summary, pathway_scores)

setorder(kla_summary, PXD, ConditionLabel, SampleID)
setorder(paired_summary, PXD, ConditionLabel, SampleID)
setorder(pathway_profile, PXD, ConditionLabel, SampleID, Pathway)
setorder(reconciliation, PXD, PublicationGroup)

fwrite(kla_summary, file.path(candidate_dir, "sample_resolved_kla_summary.csv"), na = "")
fwrite(paired_summary, file.path(candidate_dir, "sample_resolved_matched_modalities.csv"), na = "")
fwrite(pathway_profile, file.path(candidate_dir, "sample_resolved_pathway_profile.csv"), na = "")
fwrite(reconciliation, file.path(candidate_dir, "sample_resolved_source_reconciliation.csv"), na = "")

message("Wrote candidate sample-resolved inputs for ", nrow(kla_summary), " source samples and ", nrow(paired_summary), " same-sample modality pairs.")
