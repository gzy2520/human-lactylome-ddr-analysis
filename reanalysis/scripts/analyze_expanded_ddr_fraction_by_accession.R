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
project_root <- if (length(args)) normalizePath(args[[1]]) else normalizePath(".")
table_dir <- file.path(project_root, "reanalysis", "results", "tables")
figure_dir <- file.path(project_root, "reanalysis", "results", "figures")
intermediate_dir <- file.path(
  project_root, "reanalysis", "intermediate", "expanded_ddr_by_accession"
)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(intermediate_dir, recursive = TRUE, showWarnings = FALSE)

pairing_path <- file.path(
  table_dir, "lactylome_and_reference_proteome_pairing_zh.csv"
)
kla_scope_path <- file.path(
  table_dir, "kla_regulator_intensity_availability_audit.csv"
)
primary_path <- file.path(
  project_root, "reanalysis", "intermediate", "kla_by_dataset",
  "all_primary_sample_level_kla_sites.csv"
)
existing_reference_path <- file.path(
  table_dir, "reference_proteome_all_proteins.csv"
)
go_path <- file.path(
  project_root, "data", "annotations", "GO-repair+damage(human).tsv"
)
base_accession <- function(values) {
  values <- trimws(as.character(values))
  values <- sub("^.*\\|([^|]+)\\|.*$", "\\1", values)
  values <- sub("^([^|;]+)\\|.*$", "\\1", values)
  values <- sub("^NX_", "", values)
  values <- sub("-[0-9]+$", "", values)
  values
}

relative_path <- function(path) {
  sub(
    paste0("^", project_root, "/"),
    "",
    normalizePath(path)
  )
}

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

split_accessions <- function(values) {
  values <- unlist(strsplit(as.character(values), "[;,]"))
  values <- base_accession(values)
  sort(unique(values[is_uniprot(values)]))
}

split_protein_identifiers <- function(values) {
  values <- unlist(strsplit(as.character(values), "[;,]"))
  values <- trimws(values)
  uniprot <- base_accession(values)
  uniprot <- uniprot[is_uniprot(uniprot)]
  ensembl <- sub("_RNA$", "", values)
  ensembl <- sub("\\.[0-9]+$", "", ensembl)
  ensembl <- ensembl[grepl("^ENSP[0-9]+$", ensembl)]
  sort(unique(c(uniprot, ensembl)))
}

biomart_mapping_path <- file.path(
  project_root,
  "reanalysis", "config", "ensembl_protein_to_uniprot_biomart.tsv"
)
ensembl_mapping_records <- if (file.exists(biomart_mapping_path)) {
  read.delim(
    biomart_mapping_path,
    check.names = FALSE,
    stringsAsFactors = FALSE
  ) |>
    transmute(
      EnsemblProteinID,
      BaseAccession = base_accession(BaseAccession)
    ) |>
    filter(is_uniprot(BaseAccession)) |>
    distinct()
} else {
  data.frame(
    EnsemblProteinID = character(),
    BaseAccession = character(),
    stringsAsFactors = FALSE
  )
}

map_ensembl_proteins <- function(values) {
  values <- split_protein_identifiers(values)
  values <- sort(unique(values[grepl("^ENSP[0-9]+$", values)]))
  if (!length(values)) return(character())
  missing_values <- setdiff(
    values,
    unique(ensembl_mapping_records$EnsemblProteinID)
  )
  if (!length(missing_values)) {
    return(sort(unique(
      ensembl_mapping_records$BaseAccession[
        ensembl_mapping_records$EnsemblProteinID %in% values
      ]
    )))
  }
  if (!requireNamespace("AnnotationDbi", quietly = TRUE) ||
      !requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    stop("AnnotationDbi and org.Hs.eg.db are required for ENSEMBLPROT mapping")
  }
  valid_keys <- AnnotationDbi::keys(
    org.Hs.eg.db::org.Hs.eg.db,
    keytype = "ENSEMBLPROT"
  )
  query_values <- intersect(missing_values, valid_keys)
  if (!length(query_values)) {
    return(sort(unique(
      ensembl_mapping_records$BaseAccession[
        ensembl_mapping_records$EnsemblProteinID %in% values
      ]
    )))
  }
  mapped <- suppressMessages(AnnotationDbi::select(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = query_values,
    keytype = "ENSEMBLPROT",
    columns = "UNIPROT"
  )) |>
    transmute(
      EnsemblProteinID = ENSEMBLPROT,
      BaseAccession = base_accession(UNIPROT)
    ) |>
    filter(is_uniprot(BaseAccession)) |>
    distinct()
  ensembl_mapping_records <<- bind_rows(
    ensembl_mapping_records,
    mapped
  ) |>
    distinct()
  sort(unique(mapped$BaseAccession))
}

ddr_matching_identifiers <- function(values) {
  values <- sort(unique(as.character(values)))
  direct <- values[is_uniprot(values) & values %in% ddr_accessions]
  ensembl <- values[grepl("^ENSP[0-9]+$", values)]
  if (length(ensembl)) {
    map_ensembl_proteins(ensembl)
  }
  mapped_ensembl <- ensembl_mapping_records$EnsemblProteinID[
    ensembl_mapping_records$EnsemblProteinID %in% ensembl &
      ensembl_mapping_records$BaseAccession %in% ddr_accessions
  ]
  sort(unique(c(direct, mapped_ensembl)))
}

identifier_type <- function(values) {
  ifelse(grepl("^ENSP[0-9]+$", values), "ENSEMBLPROT", "UniProtKB")
}

mapped_accessions_for <- function(values) {
  direct <- ifelse(is_uniprot(values), values, NA_character_)
  ensembl <- values[grepl("^ENSP[0-9]+$", values)]
  if (length(ensembl)) {
    map_ensembl_proteins(ensembl)
  }
  mapped_lookup <- split(
    ensembl_mapping_records$BaseAccession,
    ensembl_mapping_records$EnsemblProteinID
  )
  vapply(seq_along(values), function(index) {
    value <- values[[index]]
    if (is_uniprot(value)) return(value)
    matches <- sort(unique(mapped_lookup[[value]]))
    if (!length(matches)) "" else paste(matches, collapse = ";")
  }, character(1))
}

extract_pxd050470_reference <- function(path) {
  data <- read_excel(
    path,
    sheet = "Sheet1",
    skip = 5
  )
  required_columns <- c(
    "Protein accession",
    "Intensity_H072",
    "Intensity_H081",
    "Intensity_H187"
  )
  if (!all(required_columns %in% names(data))) {
    stop("Unexpected PXD050470 Table S4 layout: ", path)
  }
  intensity_columns <- required_columns[-1]
  keep <- rowSums(
    sapply(data[intensity_columns], function(values) {
      values <- suppressWarnings(as.numeric(values))
      is.finite(values) & values > 0
    }),
    na.rm = TRUE
  ) > 0
  accessions <- split_accessions(data$`Protein accession`[keep])
  if (length(accessions) != 6082) {
    stop("Expected 6082 unique UniProt accessions in PXD050470 Table S4")
  }
  accessions
}

valid_maxquant_rows <- function(data) {
  keep <- rep(TRUE, nrow(data))
  if ("Reverse" %in% names(data)) {
    keep <- keep & (is.na(data$Reverse) | data$Reverse != "+")
  }
  if ("Potential contaminant" %in% names(data)) {
    keep <- keep &
      (is.na(data$`Potential contaminant`) |
         data$`Potential contaminant` != "+")
  }
  if ("Only identified by site" %in% names(data)) {
    keep <- keep &
      (is.na(data$`Only identified by site`) |
         data$`Only identified by site` != "+")
  }
  keep[is.na(keep)] <- FALSE
  keep
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

extract_maxquant_sites <- function(
  path,
  sample_tokens = NULL,
  sheet = NULL
) {
  data <- if (grepl("\\.xlsx$", path, ignore.case = TRUE)) {
    read_excel(path, sheet = sheet)
  } else {
    read_delimited(path)
  }
  keep <- valid_maxquant_rows(data)
  if ("id" %in% names(data)) {
    keep <- keep & !is.na(data$id)
  }
  localization_columns <- if (is.null(sample_tokens)) {
    intersect("Localization prob", names(data))
  } else {
    intersect(paste("Localization prob", sample_tokens), names(data))
  }
  if (length(localization_columns)) {
    localized <- rowSums(
      sapply(
        data[localization_columns],
        function(x) suppressWarnings(as.numeric(x)) > 0
      ),
      na.rm = TRUE
    ) > 0
    keep <- keep & localized
  }
  split_accessions(data$Proteins[keep])
}

extract_maxquant_proteins <- function(path, abundance_pattern = NULL) {
  data <- read_delimited(path)
  keep <- valid_maxquant_rows(data)
  if (!is.null(abundance_pattern)) {
    columns <- grep(abundance_pattern, names(data), value = TRUE)
    if (!length(columns)) return(character())
    present <- rowSums(
      sapply(
        data[columns],
        function(x) {
          values <- suppressWarnings(as.numeric(gsub(",", "", x, fixed = TRUE)))
          !is.na(values) & values > 0
        }
      ),
      na.rm = TRUE
    ) > 0
    keep <- keep & present
  }
  accession_column <- intersect(
    c("Majority protein IDs", "Protein IDs"), names(data)
  )[[1]]
  raw_ids <- data[[accession_column]][keep]
  identifiers <- split_protein_identifiers(raw_ids)
  map_ensembl_proteins(identifiers)
  identifiers
}

extract_pd_proteins <- function(path, sample_token, lactylome = FALSE) {
  data <- read.delim(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    quote = "\"",
    comment.char = ""
  )
  columns <- grep(
    paste0("Found in Sample:.*", sample_token),
    names(data),
    value = TRUE
  )
  keep <- rowSums(
    sapply(data[columns], function(x) {
      !is.na(x) & nzchar(x) & x != "Not Found"
    }),
    na.rm = TRUE
  ) > 0
  split_accessions(data$Accession[keep])
}

extract_pd_lactyl_peptides <- function(path, sample_token) {
  data <- read.delim(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    quote = "\"",
    comment.char = ""
  )
  columns <- grep(
    paste0("Found in Sample:.*", sample_token),
    names(data),
    value = TRUE
  )
  keep <- grepl(
    "Lacty|Lactyl|La \\(K\\)",
    data$Modifications,
    ignore.case = TRUE
  ) &
    rowSums(
      sapply(data[columns], function(x) {
        !is.na(x) & nzchar(x) & x != "Not Found"
      }),
      na.rm = TRUE
    ) > 0
  split_accessions(data$`Master Protein Accessions`[keep])
}

extract_spectronaut_proteins <- function(
  path,
  lactyl_pattern = NULL,
  group_pattern = NULL,
  accession_column = "Protein.Group"
) {
  header <- names(fread(path, nrows = 0, data.table = FALSE))
  needed <- unique(intersect(
    c(
      accession_column, "Modified.Sequence", "Q.Value", "PG.Q.Value",
      "PTM.Site.Confidence", "R.Condition"
    ),
    header
  ))
  data <- fread(
    path,
    select = needed,
    showProgress = FALSE,
    data.table = FALSE
  )
  keep <- rep(TRUE, nrow(data))
  if (!is.null(lactyl_pattern)) {
    keep <- keep & grepl(lactyl_pattern, data$Modified.Sequence)
  }
  if ("Q.Value" %in% names(data)) {
    keep <- keep & suppressWarnings(as.numeric(data$Q.Value)) <= 0.01
  }
  if ("PG.Q.Value" %in% names(data)) {
    keep <- keep & suppressWarnings(as.numeric(data$PG.Q.Value)) <= 0.01
  }
  if ("PTM.Site.Confidence" %in% names(data) && !is.null(lactyl_pattern)) {
    keep <- keep &
      suppressWarnings(as.numeric(data$PTM.Site.Confidence)) > 0
  }
  if (!is.null(group_pattern) && "R.Condition" %in% names(data)) {
    keep <- keep & grepl(group_pattern, data$R.Condition)
  }
  split_accessions(data[[accession_column]][keep])
}

extract_spectronaut_matrix <- function(path, lactyl_pattern) {
  data <- read_delimited(path)
  keep <- grepl(lactyl_pattern, data$Modified.Sequence)
  split_accessions(data$Protein.Group[keep])
}

extract_spectronaut_quant <- function(path, group_pattern = NULL) {
  data <- read_delimited(path)
  columns <- grep("\\.PG\\.Quantity$", names(data), value = TRUE)
  if (!is.null(group_pattern)) {
    columns <- columns[grepl(group_pattern, columns)]
  }
  keep <- rowSums(
    sapply(data[columns], function(x) {
      values <- suppressWarnings(as.numeric(x))
      !is.na(values) & values > 0
    }),
    na.rm = TRUE
  ) > 0
  split_accessions(data$PG.ProteinGroups[keep])
}

extract_pxd030304_row <- function(row_key) {
  path <- file.path(
    project_root, "data", "PXD030304", "search_results",
    "ProCan-DepMapSanger_protein_matrix_6692_averaged.txt"
  )
  data <- read.delim(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = ""
  )
  row <- data[data[[1]] == row_key, , drop = FALSE]
  if (nrow(row) != 1) stop("Cannot find PXD030304 row: ", row_key)
  accessions <- base_accession(sub(";.*$", "", names(data)[-1]))
  present <- !is.na(row[1, -1]) & trimws(as.character(row[1, -1])) != ""
  sort(unique(accessions[present & is_uniprot(accessions)]))
}

extract_pxd073311_huvec_reference <- function(path) {
  data <- read_delimited(path)
  sample_columns <- grep("A0h_[123]\\.raw$", names(data), value = TRUE)
  if (length(sample_columns) != 3 || !"Protein.Group" %in% names(data)) {
    stop(
      "PXD073311 ordinary-proteome matrix must contain Protein.Group and ",
      "exactly three A0h baseline columns"
    )
  }
  detected <- rowSums(
    sapply(data[sample_columns], function(values) {
      values <- suppressWarnings(as.numeric(values))
      !is.na(values) & values > 0
    }),
    na.rm = TRUE
  ) > 0
  split_accessions(data$Protein.Group[detected])
}

extract_peaks_reference <- function(path) {
  data <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (!"Accession" %in% names(data)) return(character())
  area_columns <- grep("^Area", names(data), value = TRUE)
  keep <- if (length(area_columns)) {
    rowSums(
      sapply(data[area_columns], function(values) {
        values <- suppressWarnings(as.numeric(values))
        !is.na(values) & values > 0
      }),
      na.rm = TRUE
    ) > 0
  } else {
    rep(TRUE, nrow(data))
  }
  split_accessions(data$Accession[keep])
}

extract_maxquant_zip_proteins <- function(
  zip_path,
  abundance_pattern = NULL
) {
  data <- read.delim(
    unz(zip_path, "proteinGroups.txt"),
    check.names = FALSE,
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = ""
  )
  keep <- valid_maxquant_rows(data)
  if (!is.null(abundance_pattern)) {
    columns <- grep(abundance_pattern, names(data), value = TRUE)
    if (!length(columns)) return(character())
    present <- rowSums(
      sapply(data[columns], function(values) {
        values <- suppressWarnings(as.numeric(values))
        !is.na(values) & values > 0
      }),
      na.rm = TRUE
    ) > 0
    keep <- keep & present
  }
  accession_column <- intersect(
    c("Majority protein IDs", "Protein IDs"),
    names(data)
  )[[1]]
  split_accessions(data[[accession_column]][keep])
}

extract_pxd069969_reference <- function(path, sample_names) {
  data <- read_excel(path, sheet = "Annotation_Combine")
  columns <- intersect(
    paste0("LFQ intensity ", sample_names),
    names(data)
  )
  if (!length(columns)) return(character())
  keep <- rowSums(
    sapply(data[columns], function(values) {
      values <- suppressWarnings(as.numeric(values))
      !is.na(values) & values > 0
    }),
    na.rm = TRUE
  ) > 0
  split_accessions(data$`Protein accession`[keep])
}

extract_pxd055025_reference <- function(path) {
  data <- read_excel(path, sheet = "B VS A")
  columns <- intersect(paste0("Abundances B", 1:3), names(data))
  keep <- rowSums(
    sapply(data[columns], function(values) {
      values <- suppressWarnings(as.numeric(values))
      !is.na(values) & values > 0
    }),
    na.rm = TRUE
  ) > 0
  split_accessions(data$Accession[keep])
}

extract_pxd065775_reference <- function(path, sheet_name) {
  data <- read_excel(path, sheet = sheet_name)
  columns <- intersect(
    c(
      "Non-rec1", "Non-rec2", "Non-rec3", "Non-rec4",
      "Rec1", "Rec2", "Rec3", "Rec4"
    ),
    names(data)
  )
  keep <- rowSums(
    sapply(data[columns], function(values) {
      values <- suppressWarnings(as.numeric(values))
      !is.na(values) & values > 0
    }),
    na.rm = TRUE
  ) > 0
  split_accessions(data$Accession[keep])
}

extract_pxd059985_reference <- function(path) {
  data <- read.delim(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = ""
  )
  columns <- grep("\\.PG\\.Quantity$", names(data), value = TRUE)
  columns <- columns[grepl("AC16|V8", columns)]
  keep <- rowSums(
    sapply(data[columns], function(values) {
      values <- suppressWarnings(as.numeric(values))
      !is.na(values) & values > 0
    }),
    na.rm = TRUE
  ) > 0
  split_accessions(data$PG.ProteinGroups[keep])
}

go <- read.delim(
  go_path,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = ""
)
go_keep <- (is.na(go$`TAXON ID`) | go$`TAXON ID` == 9606) &
  !grepl("(^|\\|)NOT($|\\|)", go$QUALIFIER) &
  go$`GENE PRODUCT DB` == "UniProtKB"
ddr_accessions <- split_accessions(go$`GENE PRODUCT ID`[go_keep])

kla_scope <- read.csv(
  kla_scope_path,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  na.strings = c("", "NA")
) |>
  filter(`定量可用` %in% c(TRUE, "TRUE", "True", 1, "1")) |>
  transmute(
    乳酸化PXD = PXD,
    样本组,
    KlaRowOrder = row_number()
  )
if (nrow(kla_scope) != 37) {
  stop("Kla quantitative scope must contain exactly 37 sample groups")
}

pairing <- read.csv(
  pairing_path,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  na.strings = c("", "NA")
) |>
  inner_join(kla_scope, by = c("乳酸化PXD", "样本组")) |>
  mutate(
    SampleGroupID = paste(乳酸化PXD, 样本组, sep = "__"),
    RowOrder = KlaRowOrder
  )
if (nrow(pairing) != 37) {
  stop("Reference comparison scope must retain all 37 Kla sample groups")
}

primary <- read.csv(
  primary_path,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  na.strings = c("", "NA")
)
existing_reference <- read.csv(
  existing_reference_path,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  na.strings = c("", "NA")
)

core_alias <- c(
  "PXD014870__MCF7" = "MCF7",
  "PXD028488__HEK293T" = "HEK293T",
  "PXD028488__HCT116" = "HCT116",
  "PXD028488__TALL-104" = "T-ALL",
  "PXD050470__human hippocampus" = "Human hippocampus",
  "PXD053474__HCT116" = "HCT116",
  "PXD060185__MCF10A" = "MCF10A",
  "PXD060185__MCF7" = "MCF7",
  "PXD060185__MDA-MB-468" = "MDA-MB-468",
  "PXD060185__T-47D" = "T-47D",
  "PXD078013__RKO WT and GSK3B KO" = "RKO",
  "PXD078736__HK-2 control and mannitol" = "HK-2"
)

kla_set <- function(pxd, group) {
  key <- paste(pxd, group, sep = "__")
  if (key %in% names(core_alias)) {
    values <- primary$BaseAccession[
      primary$PXD == pxd &
        primary$CellOrTissueType == core_alias[[key]] &
        primary$PrimaryIncluded %in% c(TRUE, "TRUE", "True", 1, "1")
    ]
    return(split_accessions(values))
  }

  if (pxd == "PXD028737") {
    return(extract_maxquant_sites(
      file.path(
        project_root,
        "data/PXD028737/search_results/extracted_pairing/combined/combined/txt/La(K)Sites.txt"
      ),
      c("H0", "H24")
    ))
  }
  if (pxd == "PXD033146") {
    return(extract_maxquant_sites(file.path(
      project_root,
      "data/PXD033146/search_results/extracted_pairing/search_result-HA119TPLa/La (K)Sites.txt"
    )))
  }
  if (pxd == "PXD036307") {
    return(extract_maxquant_sites(
      file.path(
        project_root,
        "data/PXD036307/search_results/extracted/txt/La (K)Sites.txt"
      ),
      c("PTB340", "PTB342", "PTB344", "PTB346", "PTB364", "PTB372")
    ))
  }
  if (pxd == "PXD037371") return(character())
  if (pxd == "PXD046800") {
    token <- if (group == "hypertrophic scar") "HSP" else "NSP"
    return(extract_pd_lactyl_peptides(
      file.path(
        project_root,
        "data/PXD046800/search_results/HFX2_LFQ_QB001_Lacty_PeptideGroups.txt"
      ),
      token
    ))
  }
  if (pxd == "PXD050147") {
    return(extract_maxquant_sites(file.path(
      project_root, "data/PXD050147/search_results/Lactyl_K_Sites.txt"
    )))
  }
  if (pxd == "PXD054919") {
    data <- read_excel(
      file.path(
        project_root,
        "data/PXD054919/supplementary/41419_2025_8113_MOESM2_ESM.xlsx"
      ),
      skip = 1
    )
    return(split_accessions(data$`Protein accession`))
  }
  if (pxd == "PXD055230") {
    return(extract_spectronaut_proteins(
      file.path(project_root, "data/PXD055230/search_results/LaIP_HSV1_DIA.tsv"),
      "K\\(UniMod:378\\)"
    ))
  }
  if (pxd == "PXD057709") {
    return(extract_spectronaut_proteins(
      file.path(project_root, "data/PXD057709/search_results/LaIP_report.tsv"),
      "K\\(UniMod:378\\)"
    ))
  }
  if (pxd == "PXD058534") {
    return(extract_maxquant_sites(
      file.path(
        project_root,
        "data/PXD058534/search_results/extracted_pairing/txt/La (K)Sites.txt"
      ),
      "HK2"
    ))
  }
  if (pxd == "PXD062720") {
    return(extract_maxquant_sites(
      file.path(
        project_root,
        "data/PXD062720/search_results/extracted_pairing/txt/La (K)Sites.txt"
      ),
      c("A_1", "A_3")
    ))
  }
  if (pxd == "PXD063047") {
    tokens <- if (group == "normal pregnancy placenta") {
      c("Con_1", "Con_2", "Con_3")
    } else {
      c("PE_1", "PE_2", "PE_3")
    }
    return(extract_maxquant_sites(
      file.path(
        project_root,
        "data/PXD063047/search_results/extracted/combined/txt/La (K)Sites.txt"
      ),
      tokens
    ))
  }
  if (pxd == "PXD063266") {
    return(extract_maxquant_sites(
      file.path(project_root, "data/PXD063266/search_results/LactylSites.xlsx"),
      c("1", "2", "3"),
      sheet = "Lactyl (K)Sites"
    ))
  }
  if (pxd == "PXD064038") {
    return(extract_maxquant_sites(
      file.path(
        project_root,
        "data/PXD064038/search_results/extracted_pairing/txt/txt/La (K)Sites.txt"
      ),
      c("MEC_1", "MEC_2", "MEC_3", "NEC_1", "NEC_2", "NEC_3")
    ))
  }
  if (pxd == "PXD064912") {
    data <- read_excel(
      file.path(
        project_root, "data/PXD064912/supplementary/europepmc/mmc1.xlsx"
      ),
      skip = 1
    )
    probability_columns <- grep("^PTM.SiteProbability", names(data), value = TRUE)
    keep <- tolower(data$PTM.ModificationTitle) == "lactylation" &
      data$PTM.SiteAA == "K" &
      rowSums(
        sapply(
          data[probability_columns],
          function(x) suppressWarnings(as.numeric(x)) > 0
        ),
        na.rm = TRUE
      ) > 0
    return(split_accessions(data$PTM.ProteinId[keep]))
  }
  if (pxd == "PXD066054") {
    data <- read_delimited(file.path(
      project_root,
      "data/PXD066054/search_results/extracted/PLa/XB08700B1DPLa_-PTMSiteReport.tsv"
    ))
    prefix <- if (group == "BPH") "^NAT" else "^PCa"
    keep <- data$PTM.ModificationTitle == "L-Lac(K)" &
      data$PTM.SiteAA == "K" &
      grepl(prefix, data$R.Condition) &
      suppressWarnings(as.numeric(data$PTM.SiteProbability)) > 0
    return(split_accessions(data$PTM.ProteinId[keep]))
  }
  if (pxd == "PXD066351") {
    data <- read.csv(
      file.path(
        project_root,
        "data/PXD066351/search_results/XB01472B1DPLa-MSstats_Input.csv"
      ),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    keep <- data$PTM.ModificationTitle == "Lac (K)" &
      data$PTM.SiteAA == "K" &
      suppressWarnings(as.numeric(data$PTM.SiteProbability)) > 0 &
      suppressWarnings(as.numeric(data$PTM.Quantity)) > 0
    return(split_accessions(data$PTM.ProteinId[keep]))
  }
  if (pxd == "PXD070007") {
    data <- read_excel(
      file.path(
        project_root,
        "data/PXD070007/search_results/SA206LPLaB1_Annotation.xlsx"
      ),
      sheet = "Annotation_Combine"
    )
    columns <- if (group == "glioblastoma stem cells") {
      c("G2907", "G3028", "G3264", "GSC23", "MES28", "RKI")
    } else {
      c("ENSA", "HMP1")
    }
    keep <- suppressWarnings(as.numeric(data$`Localization probability`)) > 0 &
      rowSums(
        sapply(data[columns], function(x) suppressWarnings(as.numeric(x)) > 0),
        na.rm = TRUE
      ) > 0
    return(split_accessions(data$`Protein accession`[keep]))
  }
  if (pxd == "PXD073311") {
    return(extract_spectronaut_matrix(
      file.path(
        project_root,
        paste0(
          "data/PXD073311/search_results/extracted_pairing/",
          "IPX0015307003_Database_search_result/Database_search_result/",
          "report.pr_matrix.tsv"
        )
      ),
      "K\\(UniMod:2114\\)"
    ))
  }
  if (pxd == "PXD075014") {
    data <- read_excel(
      file.path(project_root, "data/PXD075014/supplementary/Table2.XLSX"),
      sheet = "Kla peptides"
    )
    keep <- grepl("lactylation", data$Modifications, ignore.case = TRUE)
    return(split_accessions(data$`Master Protein Accessions`[keep]))
  }
  if (pxd == "PXD075377") {
    data <- read_delimited(file.path(
      project_root,
      "data/PXD075377/search_results/extracted/2-Basic_analysis/MS_identified_information.txt"
    ))
    intensity_column <- if (group == "HCC") "Intensity HCC" else "Intensity Control"
    keep <- suppressWarnings(as.numeric(data[[intensity_column]])) > 0 &
      suppressWarnings(as.numeric(data$`Localization probability`)) > 0
    return(split_accessions(data$`Protein accession`[keep]))
  }
  stop("No Kla accession parser for ", pxd, " / ", group)
}

reference_set <- function(row) {
  pxd <- row$乳酸化PXD
  group <- row$样本组
  reference_pxd <- row$常规蛋白组PXD

  existing_alias <- c(
    "PXD028488__HEK293T" = "HEK293T",
    "PXD028488__HCT116" = "HCT116",
    "PXD053474__HCT116" = "HCT116",
    "PXD058534__pretreated HK-2" = "HK-2",
    "PXD060185__MCF10A" = "MCF10A",
    "PXD060185__MCF7" = "MCF7",
    "PXD060185__MDA-MB-468" = "MDA-MB-468",
    "PXD060185__T-47D" = "T-47D",
    "PXD078013__RKO WT and GSK3B KO" = "RKO",
    "PXD078736__HK-2 control and mannitol" = "HK-2"
  )
  key <- paste(pxd, group, sep = "__")
  if (pxd == "PXD028488" && group == "TALL-104") {
    return(extract_peaks_reference(file.path(
      project_root,
      row$常规蛋白组证据文件
    )))
  }
  if (pxd == "PXD028737" && group == "HMC3") {
    return(extract_maxquant_proteins(
      file.path(project_root, row$常规蛋白组证据文件),
      "LFQ intensity H0|LFQ intensity H24"
    ))
  }
  if (
    reference_pxd == "PXD055025" &&
      group == "severe preeclampsia placenta"
  ) {
    return(extract_pxd055025_reference(file.path(
      project_root,
      row$常规蛋白组证据文件
    )))
  }
  if (reference_pxd == "PXD022005" && group == "PC-3M") {
    return(extract_maxquant_zip_proteins(
      file.path(project_root, row$常规蛋白组证据文件),
      "^Intensity H$"
    ))
  }
  if (
    reference_pxd == "PXD069969" &&
      group %in% c("glioblastoma stem cells", "neural stem cells")
  ) {
    samples <- if (group == "glioblastoma stem cells") {
      c("G2907", "G3028", "G3264", "GSC23", "MES28", "RKI")
    } else {
      c("ENSA", "HMP1")
    }
    return(extract_pxd069969_reference(
      file.path(project_root, row$常规蛋白组证据文件),
      samples
    ))
  }
  if (
    reference_pxd == "PXD059985" &&
      group == "AC16 control and hypoxia"
  ) {
    return(extract_pxd059985_reference(file.path(
      project_root,
      row$常规蛋白组证据文件
    )))
  }
  if (
    reference_pxd == "PXD065775" &&
      group %in% c("HCC", "adjacent liver")
  ) {
    return(extract_pxd065775_reference(
      file.path(project_root, row$常规蛋白组证据文件),
      ifelse(group == "HCC", "CISs", "ANTs")
    ))
  }
  if (key %in% names(existing_alias)) {
    return(split_accessions(
      existing_reference$BaseAccession[
        existing_reference$CellOrTissueType == existing_alias[[key]]
      ]
    ))
  }

  if (pxd == "PXD014870") {
    return(extract_pxd030304_row("SIDM00148;MCF7"))
  }
  if (
    reference_pxd == "PXD050470" &&
      group == "human hippocampus"
  ) {
    return(extract_pxd050470_reference(file.path(
      project_root,
      row$常规蛋白组证据文件
    )))
  }
  if (reference_pxd == "PXD010154") {
    return(extract_maxquant_proteins(
      file.path(project_root, row$常规蛋白组证据文件)
    ))
  }
  if (reference_pxd == "PXD033146") {
    return(extract_maxquant_proteins(file.path(
      project_root,
      "data/PXD033146/search_results/extracted_pairing/search_result-HA119TQ/proteinGroups.txt"
    )))
  }
  if (reference_pxd == "PXD046800") {
    token <- if (group == "hypertrophic scar") "HSP" else "NSP"
    return(extract_pd_proteins(
      file.path(
        project_root,
        "data/PXD046800/search_results/HFX2_LFQ_QB002_Proteins.txt"
      ),
      token
    ))
  }
  if (reference_pxd == "PXD050147") {
    return(extract_maxquant_proteins(file.path(
      project_root, "data/PXD050147/search_results/SIRT_proteinGroups.txt"
    )))
  }
  if (reference_pxd == "PXD055230") {
    return(extract_spectronaut_proteins(file.path(
      project_root, "data/PXD055230/search_results/WP_HSV1_DIA.tsv"
    )))
  }
  if (reference_pxd == "PXD057709") {
    return(extract_spectronaut_proteins(file.path(
      project_root, "data/PXD057709/search_results/WP_report.tsv"
    )))
  }
  if (reference_pxd == "PXD062720") {
    return(extract_maxquant_proteins(file.path(
      project_root,
      "data/PXD062720/search_results/extracted_pairing/txt/proteinGroups.txt"
    )))
  }
  if (reference_pxd == "PXD066517") {
    data <- read_delimited(file.path(
      project_root, "data/PXD066517/search_results/20240275.tsv"
    ))
    return(split_accessions(data$PG.ProteinAccessions))
  }
  if (reference_pxd == "PXD066054") {
    prefix <- if (group == "BPH") "NAT" else "PCa"
    return(extract_spectronaut_quant(
      file.path(
        project_root,
        "data/PXD066054/search_results/extracted/DA/Protein_Quant.tsv"
      ),
      prefix
    ))
  }
  if (reference_pxd == "PXD066351") {
    return(extract_spectronaut_quant(file.path(
      project_root,
      "data/PXD066351/search_results/XB01472B1DA-Protein_Quant.tsv"
    )))
  }
  if (reference_pxd == "PXD073311") {
    return(extract_pxd073311_huvec_reference(file.path(
      project_root,
      row$常规蛋白组证据文件
    )))
  }
  if (reference_pxd == "PXD030304" && group == "A549") {
    return(extract_pxd030304_row("SIDM00903;A549"))
  }
  if (reference_pxd == "PXD030304" && group == "PC-3M") {
    return(extract_pxd030304_row("SIDM00088;PC-3"))
  }
  stop("No reference accession parser for ", pxd, " / ", group)
}

display_names <- c(
  "human hippocampus" = "人海马组织",
  "pathological rotator cuff tendon" = "病理性肩袖肌腱",
  "normal human lung" = "正常人肺组织",
  "hypertrophic scar" = "增生性瘢痕",
  "adjacent skin" = "邻近皮肤",
  "HepG2 WT and SIRT1 or SIRT3 KO" = "HepG2 WT/SIRT1-KO/SIRT3-KO",
  "human fibroblasts mock and HCMV or HSV-1" = "人成纤维细胞 mock/HCMV/HSV-1",
  "human fibroblasts mock and HCMV" = "人成纤维细胞 mock/HCMV",
  "pretreated HK-2" = "预处理HK-2",
  "bladder cancer cells treated with EPI" = "EPI处理膀胱癌细胞",
  "normal pregnancy placenta" = "正常妊娠胎盘",
  "severe preeclampsia placenta" = "重度子痫前期胎盘",
  "MEC and NEC ESCC groups" = "ESCC MEC/NEC组",
  "human sperm" = "人精子",
  "BPH" = "良性前列腺增生",
  "prostate cancer" = "前列腺癌",
  "HCT116 control and Roseburia co-culture" = "HCT116对照/Roseburia共培养",
  "glioblastoma stem cells" = "胶质母细胞瘤干细胞",
  "neural stem cells" = "神经干细胞",
  "HUVEC control and Pg infection" = "HUVEC对照/Pg感染",
  "AC16 control and hypoxia" = "AC16对照/低氧",
  "HCC" = "肝细胞癌",
  "adjacent liver" = "邻近肝组织",
  "RKO WT and GSK3B KO" = "RKO WT/GSK3B-KO",
  "HK-2 control and mannitol" = "HK-2对照/甘露醇",
  "TALL-104" = "TALL-104"
)

statistics <- list()
audit <- list()
kla_members <- list()
reference_members <- list()

for (i in seq_len(nrow(pairing))) {
  row <- pairing[i, ]
  pxd <- row$乳酸化PXD
  group <- row$样本组
  key <- row$SampleGroupID

  kla <- kla_set(pxd, group)
  if (!length(kla)) {
    stop("Unable to extract Kla accessions for ", pxd, " / ", group)
  }

  reference_configured <-
    row$配置要求进入严格参照分析 %in%
      c(TRUE, "TRUE", "True", 1, "1") &&
    !is.na(row$常规蛋白组PXD) &&
    nzchar(row$常规蛋白组PXD) &&
    !is.na(row$常规蛋白组证据文件) &&
    nzchar(row$常规蛋白组证据文件)
  reference <- if (reference_configured) reference_set(row) else character()
  paired_included <- reference_configured && length(reference) > 0
  kla_ddr <- ddr_matching_identifiers(kla)
  reference_ddr <- if (paired_included) {
    ddr_matching_identifiers(reference)
  } else {
    character()
  }
  display <- if (group %in% names(display_names)) {
    unname(display_names[[group]])
  } else {
    group
  }
  reference_pxd <- if (paired_included) row$常规蛋白组PXD else NA_character_
  reference_note <- if (paired_included) row$匹配质量 else row$注意事项
  match_mode <- "BaseAccession_only"
  reference_protein_count <- if (paired_included) length(reference) else NA_integer_
  reference_ddr_count <- if (paired_included) length(reference_ddr) else NA_integer_
  reference_fraction <- if (paired_included) {
    length(reference_ddr) / length(reference)
  } else {
    NA_real_
  }

  statistics[[length(statistics) + 1]] <- data.frame(
    PXD = pxd,
    SampleGroup = group,
    DisplayLabel = paste0(display, " · ", pxd),
    BiologicalMaterial = row$材料类型,
    KlaEvidenceFile = row$乳酸化证据文件,
    KlaProteinCount = length(kla),
    KlaDdrProteinCount = length(kla_ddr),
    KlaDdrFraction = length(kla_ddr) / length(kla),
    PairedAnalysisIncluded = paired_included,
    ReferencePXD = reference_pxd,
    ReferenceEvidenceFile = if (paired_included) {
      row$常规蛋白组证据文件
    } else {
      NA_character_
    },
    ReferenceProteinCount = reference_protein_count,
    ReferenceDdrProteinCount = reference_ddr_count,
    ReferenceDdrFraction = reference_fraction,
    DdrFractionPercentagePointDifference =
      if (paired_included) {
        length(kla_ddr) / length(kla) - reference_fraction
      } else {
        NA_real_
      },
    ReferenceMatchNote = reference_note,
    MatchMode = match_mode,
    SymbolFallbackCount = 0L,
    RowOrder = row$RowOrder,
    stringsAsFactors = FALSE
  )
  audit[[length(audit) + 1]] <- data.frame(
    PXD = pxd,
    SampleGroup = group,
    Included = TRUE,
    ReferenceIncluded = paired_included,
    Reason = if (!paired_included) {
      row$注意事项
    } else {
      ""
    },
    KlaAccessionCount = length(kla),
    ReferenceAccessionCount = ifelse(
      paired_included,
      length(reference),
      0L
    ),
    MatchMode = match_mode,
    stringsAsFactors = FALSE
  )
  kla_members[[length(kla_members) + 1]] <- data.frame(
    PXD = pxd,
    SampleGroup = group,
    BaseAccession = kla,
    IsDdr = kla %in% ddr_accessions,
    stringsAsFactors = FALSE
  )
  if (paired_included) {
    reference_members[[length(reference_members) + 1]] <- data.frame(
      PXD = pxd,
      SampleGroup = group,
      ReferencePXD = reference_pxd,
      SourceProteinID = reference,
      IdentifierType = identifier_type(reference),
      MappedBaseAccessions = mapped_accessions_for(reference),
      IsDdr = reference %in% reference_ddr,
      stringsAsFactors = FALSE
    )
  }
}

core_source_alias <- data.frame(
  PXD = sub("__.*$", "", names(core_alias)),
  SampleGroup = sub("^[^_]+__", "", names(core_alias)),
  CellOrTissueType = unname(core_alias),
  stringsAsFactors = FALSE
)

primary_source_map <- primary |>
  filter(
    PrimaryIncluded %in% c(TRUE, "TRUE", "True", 1, "1")
  ) |>
  inner_join(core_source_alias, by = c("PXD", "CellOrTissueType")) |>
  group_by(PXD, SampleGroup) |>
  summarise(
    PrimaryEvidenceFiles = paste(sort(unique(SourceFile)), collapse = ";"),
    .groups = "drop"
  )

category_map <- c(
  "PXD033146__pathological rotator cuff tendon" = "normal_tissue",
  "PXD036307__normal human lung" = "normal_tissue",
  "PXD046800__hypertrophic scar" = "normal_tissue",
  "PXD046800__adjacent skin" = "normal_tissue",
  "PXD050470__human hippocampus" = "normal_tissue",
  "PXD063047__normal pregnancy placenta" = "normal_tissue",
  "PXD063047__severe preeclampsia placenta" = "normal_tissue",
  "PXD064912__human sperm" = "normal_tissue",
  "PXD066054__BPH" = "normal_tissue",
  "PXD075377__adjacent liver" = "normal_tissue",
  "PXD064038__MEC and NEC ESCC groups" = "cancer_tissue",
  "PXD066054__prostate cancer" = "cancer_tissue",
  "PXD075377__HCC" = "cancer_tissue",
  "PXD028488__HEK293T" = "normal_cells",
  "PXD028737__HMC3" = "normal_cells",
  "PXD055230__human fibroblasts mock and HCMV or HSV-1" = "normal_cells",
  "PXD057709__human fibroblasts mock and HCMV" = "normal_cells",
  "PXD058534__pretreated HK-2" = "normal_cells",
  "PXD060185__MCF10A" = "normal_cells",
  "PXD070007__neural stem cells" = "normal_cells",
  "PXD073311__HUVEC control and Pg infection" = "normal_cells",
  "PXD075014__AC16 control and hypoxia" = "normal_cells",
  "PXD078736__HK-2 control and mannitol" = "normal_cells",
  "PXD014870__MCF7" = "cancer_cells",
  "PXD028488__HCT116" = "cancer_cells",
  "PXD028488__TALL-104" = "cancer_cells",
  "PXD050147__HepG2 WT and SIRT1 or SIRT3 KO" = "cancer_cells",
  "PXD053474__HCT116" = "cancer_cells",
  "PXD054919__A549" = "cancer_cells",
  "PXD060185__MCF7" = "cancer_cells",
  "PXD060185__MDA-MB-468" = "cancer_cells",
  "PXD060185__T-47D" = "cancer_cells",
  "PXD062720__bladder cancer cells treated with EPI" = "cancer_cells",
  "PXD063266__PC-3M" = "cancer_cells",
  "PXD066351__HCT116 control and Roseburia co-culture" = "cancer_cells",
  "PXD070007__glioblastoma stem cells" = "cancer_cells",
  "PXD078013__RKO WT and GSK3B KO" = "cancer_cells"
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
english_display_names <- c(
  "pathological rotator cuff tendon" = "Pathological rotator cuff tendon",
  "normal human lung" = "Normal human lung",
  "hypertrophic scar" = "Hypertrophic scar",
  "adjacent skin" = "Adjacent skin",
  "human hippocampus" = "Human hippocampus",
  "normal pregnancy placenta" = "Normal pregnancy placenta",
  "severe preeclampsia placenta" = "Severe preeclampsia placenta",
  "human sperm" = "Human sperm",
  "BPH" = "Benign prostatic hyperplasia",
  "adjacent liver" = "Adjacent liver",
  "MEC and NEC ESCC groups" = "ESCC MEC/NEC groups",
  "prostate cancer" = "Prostate cancer",
  "HCC" = "HCC",
  "HEK293T" = "HEK293T",
  "HMC3" = "HMC3",
  "human fibroblasts mock and HCMV or HSV-1" =
    "Human fibroblasts mock/HCMV/HSV-1",
  "human fibroblasts mock and HCMV" = "Human fibroblasts mock/HCMV",
  "pretreated HK-2" = "Pretreated HK-2",
  "MCF10A" = "MCF10A",
  "neural stem cells" = "Neural stem cells",
  "HUVEC control and Pg infection" = "HUVEC control/Pg infection",
  "AC16 control and hypoxia" = "AC16 control/hypoxia",
  "HK-2 control and mannitol" = "HK-2 control/mannitol",
  "MCF7" = "MCF7",
  "HCT116" = "HCT116",
  "TALL-104" = "TALL-104",
  "HepG2 WT and SIRT1 or SIRT3 KO" = "HepG2 WT/SIRT1-KO/SIRT3-KO",
  "A549" = "A549",
  "MDA-MB-468" = "MDA-MB-468",
  "T-47D" = "T-47D",
  "bladder cancer cells treated with EPI" = "EPI-treated bladder cancer cells",
  "PC-3M" = "PC-3M",
  "HCT116 control and Roseburia co-culture" =
    "HCT116 control/Roseburia co-culture",
  "glioblastoma stem cells" = "Glioblastoma stem cells",
  "RKO WT and GSK3B KO" = "RKO WT/GSK3B-KO"
)

statistics <- bind_rows(statistics) |>
  left_join(primary_source_map, by = c("PXD", "SampleGroup")) |>
  mutate(
    KlaEvidenceFile = case_when(
      !is.na(PrimaryEvidenceFiles) & nzchar(PrimaryEvidenceFiles) ~
        PrimaryEvidenceFiles,
      PXD == "PXD046800" ~
        "data/PXD046800/search_results/HFX2_LFQ_QB001_Lacty_PeptideGroups.txt",
      PXD == "PXD066351" ~
        "data/PXD066351/search_results/XB01472B1DPLa-MSstats_Input.csv",
      TRUE ~ KlaEvidenceFile
    ),
    Category = unname(category_map[paste(PXD, SampleGroup, sep = "__")]),
    CategoryZh = unname(category_labels_zh[Category]),
    CategoryEn = unname(category_labels_en[Category]),
    ReferenceDisplay = ifelse(
      PairedAnalysisIncluded,
      ReferencePXD,
      "未纳入：无完全匹配强度参照"
    ),
    DisplayLabel = paste0(
      DisplayLabel,
      " · Kla:",
      PXD,
      " / Ref:",
      ReferenceDisplay
    ),
    EnglishDisplayLabel = paste0(
      unname(english_display_names[SampleGroup]),
      " · Kla: ",
      PXD,
      " / Ref: ",
      ifelse(
        PairedAnalysisIncluded,
        ReferencePXD,
        "excluded: no exact quantitative reference"
      )
    )
  ) |>
  select(-PrimaryEvidenceFiles, -ReferenceDisplay) |>
  arrange(
    match(Category, c("normal_tissue", "cancer_tissue", "normal_cells", "cancer_cells")),
    RowOrder
  )
if (any(is.na(statistics$Category))) {
  stop(
    "Missing four-category classification for: ",
    paste(
      paste(statistics$PXD, statistics$SampleGroup, sep = "__")[
        is.na(statistics$Category)
      ],
      collapse = ", "
    )
  )
}
audit <- bind_rows(audit)
kla_members <- bind_rows(kla_members)
reference_members <- bind_rows(reference_members)

write.csv(
  statistics,
  file.path(table_dir, "cell_type_kla_vs_reference_ddr_statistics_accession_only.csv"),
  row.names = FALSE,
  na = ""
)
statistics_zh <- statistics |>
  transmute(
    `乳酸化PXD` = PXD,
    `样本组` = SampleGroup,
    `图中标签` = DisplayLabel,
    `四分类` = CategoryZh,
    `材料类型` = BiologicalMaterial,
    `乳酸化证据文件` = KlaEvidenceFile,
    `乳酸化蛋白ID数` = KlaProteinCount,
    `乳酸化DDR蛋白ID数` = KlaDdrProteinCount,
    `乳酸化DDR占比` = KlaDdrFraction,
    `纳入严格配对分析` = PairedAnalysisIncluded,
    `常规蛋白组PXD` = ReferencePXD,
    `常规蛋白组证据文件` = ReferenceEvidenceFile,
    `常规蛋白ID数` = ReferenceProteinCount,
    `常规DDR蛋白ID数` = ReferenceDdrProteinCount,
    `常规DDR占比` = ReferenceDdrFraction,
    `乳酸化减常规占比差` = DdrFractionPercentagePointDifference,
    `参照匹配说明` = ReferenceMatchNote,
    `交集判定方式` = MatchMode,
    `GeneSymbol回退数` = SymbolFallbackCount
  )
write.csv(
  statistics_zh,
  file.path(
    table_dir,
    "cell_type_kla_vs_reference_ddr_statistics_accession_only_zh.csv"
  ),
  row.names = FALSE,
  na = ""
)
write.csv(
  audit,
  file.path(table_dir, "cell_type_kla_vs_reference_ddr_accession_only_audit.csv"),
  row.names = FALSE,
  na = ""
)
evidence_rule <- function(pxd) {
  case_when(
    pxd %in% c(
      "PXD014870", "PXD028737", "PXD033146", "PXD036307",
      "PXD050147", "PXD058534", "PXD060185", "PXD062720",
      "PXD063047", "PXD063266", "PXD064038", "PXD078736"
    ) ~ "乳酸化位点表；定位概率>0；排除reverse和contaminant",
    pxd %in% c("PXD028488", "PXD053474") ~
      "修饰肽/PTM表明确含Lactyl或K(+72.02)，并保留蛋白accession",
    pxd == "PXD046800" ~
      "Lacty_PeptideGroups中Modifications明确含Lacty且对应样本检出",
    pxd %in% c("PXD050470", "PXD054919", "PXD070007", "PXD075377") ~
      "作者位点表；Kla位点及定位概率/样本强度明确",
    pxd %in% c("PXD055230", "PXD057709") ~
      "Spectronaut修饰肽明确含K(UniMod:378)，q值通过",
    pxd == "PXD064912" ~
      "PTM.ModificationTitle=lactylation、SiteAA=K且定位概率>0",
    pxd == "PXD066054" ~
      "PTM.ModificationTitle=L-Lac(K)、SiteAA=K且定位概率>0",
    pxd == "PXD066351" ~
      "MSstats PTM表中ModificationTitle=Lac(K)、SiteAA=K、定位概率>0且PTM.Quantity>0",
    pxd == "PXD073311" ~
      "Spectronaut修饰前体明确含K(UniMod:2114)",
    pxd == "PXD075014" ~
      "作者修饰肽表中Modifications明确为lactylation",
    pxd == "PXD078013" ~
      "evidence中La(K)>0且存在La(K) site IDs",
    TRUE ~ "直接Kla证据"
  )
}
evidence_audit <- statistics |>
  transmute(
    PXD,
    SampleGroup,
    KlaEvidenceFile,
    DirectLactylationEvidence = TRUE,
    LactylationEvidenceRule = evidence_rule(PXD),
    Pre2019ReviewRequired = PXD == "PXD014870",
    AcquisitionPeriod = case_when(
      PXD == "PXD014870" ~
        "2018; DCA/rotenone raw names contain 20181127; same study published in 2019",
      TRUE ~ "not flagged as pre-2019 lactylome in current audit"
    ),
    EnrichmentOrValidationMethod = case_when(
      PXD == "PXD014870" ~ paste0(
        "pan anti-Kla (PTM-1401) peptide immunoprecipitation after tryptic digestion; ",
        "SILAC/isotope tracing; revalidated by lactyllysine cyclic immonium ion in 2022"
      ),
      TRUE ~ evidence_rule(PXD)
    ),
    Pre2019Eligibility = case_when(
      PXD == "PXD014870" ~
        "include: satisfies pan anti-Kla enrichment requirement",
      TRUE ~ "not applicable"
    ),
    AntibodyTerminologyNote = case_when(
      PXD == "PXD014870" ~
        "pan anti-Kla is the lactylation antibody; pan anti-Kac is the acetylation control",
      TRUE ~ ""
    ),
    ProteinIdentityMode = "protein ID/accession",
    GeneSymbolFallbackCount = 0L,
    KlaProteinCount,
    KlaDdrProteinCount
  )
write.csv(
  evidence_audit,
  file.path(
    table_dir,
    "cell_type_kla_ddr_lactylation_evidence_audit.csv"
  ),
  row.names = FALSE,
  na = ""
)
write.csv(
  evidence_audit |>
    filter(Pre2019ReviewRequired) |>
    select(
      PXD,
      SampleGroup,
      AcquisitionPeriod,
      EnrichmentOrValidationMethod,
      Pre2019Eligibility,
      AntibodyTerminologyNote,
      KlaEvidenceFile,
      KlaProteinCount,
      KlaDdrProteinCount
    ),
  file.path(
    table_dir,
    "pre_2019_lactylation_dataset_review.csv"
  ),
  row.names = FALSE,
  na = ""
)
write.csv(
  kla_members,
  file.path(intermediate_dir, "kla_proteins_by_sample_group.csv"),
  row.names = FALSE
)
write.csv(
  reference_members,
  file.path(intermediate_dir, "reference_proteins_by_sample_group.csv"),
  row.names = FALSE
)
write.csv(
  ensembl_mapping_records,
  file.path(
    intermediate_dir,
    "ensembl_protein_to_uniprot_mapping_used.csv"
  ),
  row.names = FALSE
)
paired_statistics <- statistics |>
  filter(PairedAnalysisIncluded) |>
  mutate(
    SampleLabelZh = sub(" ·.*$", "", DisplayLabel),
    SampleLabelEn = sub(" ·.*$", "", EnglishDisplayLabel),
    MaterialCluster = case_when(
      grepl("HK-2", SampleGroup, fixed = TRUE) ~ "HK-2",
      SampleGroup == "MCF7" ~ "MCF7",
      grepl("^HCT116", SampleGroup) ~ "HCT116",
      grepl("^human fibroblasts", SampleGroup) ~ "human_fibroblasts",
      TRUE ~ paste(PXD, SampleGroup, sep = "__")
    ),
    MaterialClusterZh = case_when(
      MaterialCluster == "HK-2" ~ "HK-2",
      MaterialCluster == "MCF7" ~ "MCF7",
      MaterialCluster == "HCT116" ~ "HCT116",
      MaterialCluster == "human_fibroblasts" ~ "人成纤维细胞",
      TRUE ~ SampleLabelZh
    ),
    MaterialClusterEn = case_when(
      MaterialCluster == "human_fibroblasts" ~ "Human fibroblasts",
      TRUE ~ ifelse(
        MaterialCluster %in% c("HK-2", "MCF7", "HCT116"),
        MaterialCluster,
        SampleLabelEn
      )
    ),
    ReferenceDisplayKey = paste(
      ReferencePXD,
      ReferenceEvidenceFile,
      ReferenceProteinCount,
      ReferenceDdrProteinCount,
      sep = "||"
    )
  )

font_family <- "Arial Unicode MS"
category_order <- c(
  "normal_tissue",
  "cancer_tissue",
  "normal_cells",
  "cancer_cells"
)

cluster_order <- paired_statistics |>
  group_by(Category, MaterialCluster) |>
  summarise(ClusterOrder = min(RowOrder), .groups = "drop")
reference_order <- paired_statistics |>
  group_by(Category, MaterialCluster, ReferenceDisplayKey) |>
  summarise(ReferenceOrder = min(RowOrder), .groups = "drop")
paired_statistics <- paired_statistics |>
  left_join(cluster_order, by = c("Category", "MaterialCluster")) |>
  left_join(
    reference_order,
    by = c("Category", "MaterialCluster", "ReferenceDisplayKey")
  )

reference_plot_rows <- paired_statistics |>
  group_by(
    Category,
    CategoryZh,
    CategoryEn,
    MaterialCluster,
    MaterialClusterZh,
    MaterialClusterEn,
    ClusterOrder,
    ReferenceDisplayKey,
    ReferenceOrder,
    ReferencePXD,
    ReferenceEvidenceFile,
    ReferenceProteinCount,
    ReferenceDdrProteinCount,
    ReferenceDdrFraction
  ) |>
  summarise(
    LinkedKlaPXD = paste(unique(PXD), collapse = ";"),
    LinkedKlaStudyCount = n(),
    .groups = "drop"
  ) |>
  transmute(
    Category,
    CategoryZh,
    CategoryEn,
    MaterialCluster,
    MaterialClusterZh,
    MaterialClusterEn,
    ClusterOrder,
    ReferenceOrder,
    BarPriority = 0L,
    StudyOrder = ReferenceOrder,
    BarType = "reference",
    PXD = NA_character_,
    SampleGroup = "ordinary whole-proteome reference",
    ReferencePXD,
    ReferenceEvidenceFile,
    ReferenceDisplayKey,
    LinkedKlaPXD,
    LinkedKlaStudyCount,
    DdrFraction = ReferenceDdrFraction,
    Ddr = ReferenceDdrProteinCount,
    Total = ReferenceProteinCount,
    DisplayLabelZh = paste0(
      MaterialClusterZh,
      " · 普通全蛋白 Ref:",
      ReferencePXD
    ),
    DisplayLabelEn = paste0(
      MaterialClusterEn,
      " · whole proteome Ref:",
      ReferencePXD
    )
  )

kla_plot_rows <- paired_statistics |>
  transmute(
    Category,
    CategoryZh,
    CategoryEn,
    MaterialCluster,
    MaterialClusterZh,
    MaterialClusterEn,
    ClusterOrder,
    ReferenceOrder,
    BarPriority = 1L,
    StudyOrder = RowOrder,
    BarType = "kla",
    PXD,
    SampleGroup,
    ReferencePXD,
    ReferenceEvidenceFile,
    ReferenceDisplayKey,
    LinkedKlaPXD = PXD,
    LinkedKlaStudyCount = 1L,
    DdrFraction = KlaDdrFraction,
    Ddr = KlaDdrProteinCount,
    Total = KlaProteinCount,
    DisplayLabelZh = paste0(SampleLabelZh, " · Kla:", PXD),
    DisplayLabelEn = paste0(SampleLabelEn, " · Kla:", PXD)
  )

plot_data_base <- bind_rows(reference_plot_rows, kla_plot_rows) |>
  arrange(
    match(Category, category_order),
    ClusterOrder,
    ReferenceOrder,
    BarPriority,
    StudyOrder
  ) |>
  mutate(
    BarOrder = row_number(),
    LabelZh = sprintf("%s/%s（%.1f%%）", Ddr, Total, DdrFraction * 100),
    LabelEn = sprintf("%s/%s (%.1f%%)", Ddr, Total, DdrFraction * 100)
  )

write.csv(
  plot_data_base |>
    select(
      BarOrder,
      Category,
      CategoryZh,
      CategoryEn,
      MaterialCluster,
      MaterialClusterZh,
      MaterialClusterEn,
      BarType,
      PXD,
      SampleGroup,
      ReferencePXD,
      ReferenceEvidenceFile,
      ReferenceDisplayKey,
      LinkedKlaPXD,
      LinkedKlaStudyCount,
      Ddr,
      Total,
      DdrFraction,
      DisplayLabelZh,
      DisplayLabelEn
    ),
  file.path(
    table_dir,
    "cell_type_kla_vs_reference_ddr_plot_rows.csv"
  ),
  row.names = FALSE,
  na = ""
)

make_ddr_plot <- function(language = c("zh", "en")) {
  language <- match.arg(language)
  is_zh <- language == "zh"
  dataset_labels <- if (is_zh) {
    c(
      reference = "常规全蛋白组参照",
      kla = "乳酸化蛋白组（Kla）"
    )
  } else {
    c(
      reference = "Reference whole proteome",
      kla = "Lactylome (Kla)"
    )
  }
  plot_data <- plot_data_base |>
    mutate(
      PlotLabel = if (is_zh) DisplayLabelZh else DisplayLabelEn,
      CategoryLabel = if (is_zh) CategoryZh else CategoryEn,
      BarLabel = if (is_zh) LabelZh else LabelEn,
      Dataset = factor(
        BarType,
        levels = c("reference", "kla"),
        labels = unname(dataset_labels[c("reference", "kla")])
      ),
      PlotLabel = factor(
        PlotLabel,
        levels = rev(
          if (is_zh) plot_data_base$DisplayLabelZh else
            plot_data_base$DisplayLabelEn
        )
      ),
      CategoryLabel = factor(
        Category,
        levels = category_order,
        labels = unname(
          if (is_zh) category_labels_zh[category_order] else
            category_labels_en[category_order]
        )
      )
    )
  max_fraction <- max(plot_data$DdrFraction, na.rm = TRUE)
  figure_height <- max(12, nrow(paired_statistics) * 0.42 + 4.2)
  ggplot(
    plot_data,
    aes(x = DdrFraction * 100, y = PlotLabel, fill = Dataset)
  ) +
    geom_col(
      width = 0.68,
      color = "white",
      linewidth = 0.2
    ) +
    geom_text(
      aes(label = BarLabel),
      hjust = -0.04,
      size = if (is_zh) 2.65 else 2.45,
      family = font_family,
      color = "#30343B"
    ) +
    facet_grid(
      CategoryLabel ~ .,
      scales = "free_y",
      space = "free_y",
      switch = "y"
    ) +
    scale_fill_manual(
      values = setNames(
        c("#4E79A7", "#F28E2B"),
        c(dataset_labels[["reference"]], dataset_labels[["kla"]])
      )
    ) +
    scale_x_continuous(
      limits = c(0, max_fraction * 100 + 5.2),
      expand = expansion(mult = c(0, 0))
    ) +
    labs(
      title = if (is_zh) {
        "正常组织、癌症组织、正常细胞与癌症细胞中的 DDR 蛋白占比"
      } else {
        "DDR protein fraction in normal tissues, cancer tissues, normal cells, and cancer cells"
      },
      subtitle = if (is_zh) {
        paste0(
          "纳入", nrow(paired_statistics), "个Kla研究组和",
          nrow(reference_plot_rows), "个唯一普通全蛋白参照；",
          "同一参照只显示一次，蓝色参照位于其橙色Kla研究上方；",
          "GO-DDR交集按去isoform后的BaseAccession计算"
        )
      } else {
        paste0(
          nrow(paired_statistics), " Kla study groups and ",
          nrow(reference_plot_rows), " unique whole-proteome references; ",
          "each reference is displayed once above its linked Kla studies; ",
          "GO-DDR overlap uses isoform-stripped BaseAccession"
        )
      },
      x = if (is_zh) "DDR 注释蛋白占比（%）" else
        "GO-DDR annotated protein fraction (%)",
      y = NULL,
      fill = NULL,
      caption = if (is_zh) {
        paste0(
          "柱端为DDR蛋白数/总蛋白数。PXD037371的3组因TMT通道映射不明而排除；",
          "HK-2、MCF7和HCT116的共享参照蓝柱不重复；",
          "相同细胞系的不同Kla研究在组内相邻显示。"
        )
      } else {
        paste0(
          "Bar labels show DDR proteins/total proteins. Three PXD037371 groups are excluded ",
          "because TMT channels cannot be reliably mapped; shared HK-2, MCF7, and HCT116 ",
          "reference bars are not duplicated, and Kla studies of the same cell line are adjacent."
        )
      }
    ) +
    theme_minimal(base_size = 10, base_family = font_family) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(color = "#D9DDE3", linewidth = 0.45),
      axis.text.y = element_text(size = if (is_zh) 7.2 else 6.8, color = "#30343B"),
      axis.text.x = element_text(color = "#30343B"),
      strip.placement = "outside",
      strip.text.y.left = element_text(
        size = 9,
        face = "bold",
        color = "#30343B",
        angle = 90
      ),
      strip.background = element_rect(fill = "#F2F2F2", color = NA),
      panel.spacing.y = grid::unit(0.28, "lines"),
      plot.title = element_text(
        size = if (is_zh) 17 else 15, face = "bold",
        color = "#20252B", hjust = 0
      ),
      plot.subtitle = element_text(size = 9.5, color = "#525A64"),
      plot.caption = element_text(
        size = 7.8, color = "#5C626A", hjust = 0, margin = margin(t = 10)
      ),
      legend.position = "top",
      legend.justification = "right",
      legend.text = element_text(size = 9),
      plot.margin = margin(14, 28, 12, 12)
    )
}

plot_zh <- make_ddr_plot("zh")
plot_en <- make_ddr_plot("en")
figure_height <- max(12, nrow(plot_data_base) * 0.28 + 4.2)

for (path in c(
  "cell_type_kla_vs_reference_ddr_fraction_accession_only.png",
  "cell_type_kla_vs_reference_ddr_fraction.png"
)) {
  ggsave(
    file.path(figure_dir, path),
    plot_zh,
    width = 13.5,
    height = figure_height,
    dpi = 350,
    bg = "white"
  )
}
for (path in c(
  "cell_type_kla_vs_reference_ddr_fraction_accession_only.pdf",
  "cell_type_kla_vs_reference_ddr_fraction.pdf"
)) {
  ggsave(
    file.path(figure_dir, path),
    plot_zh,
    width = 13.5,
    height = figure_height,
    device = cairo_pdf,
    bg = "white"
  )
}
for (path in c(
  "cell_type_kla_vs_reference_ddr_fraction_accession_only_en.png",
  "cell_type_kla_vs_reference_ddr_fraction_en.png"
)) {
  ggsave(
    file.path(figure_dir, path),
    plot_en,
    width = 13.5,
    height = figure_height,
    dpi = 350,
    bg = "white"
  )
}
for (path in c(
  "cell_type_kla_vs_reference_ddr_fraction_accession_only_en.pdf",
  "cell_type_kla_vs_reference_ddr_fraction_en.pdf"
)) {
  ggsave(
    file.path(figure_dir, path),
    plot_en,
    width = 13.5,
    height = figure_height,
    device = cairo_pdf,
    bg = "white"
  )
}

message(
  "Expanded accession-only DDR comparison: ",
  nrow(statistics),
  " included; ",
  sum(!audit$Included),
  " excluded."
)
