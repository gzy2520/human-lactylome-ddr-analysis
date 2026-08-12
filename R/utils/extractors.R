# extractors.R - Shared dataset-specific protein/peptide extractors for kla reanalysis
# Depends on: accession_utils.R (base_accession, is_uniprot, split_accessions,
#   split_protein_identifiers) and io_utils.R (read_delimited, valid_maxquant_rows).
# Load order: accession_utils.R -> io_utils.R -> extractors.R
extractors_loaded <- TRUE

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

extract_huvec_xml <- function(path) {
  lines <- readLines(path, warn = FALSE)
  hits <- unlist(str_extract_all(lines, 'protein_name="NX_[A-Z0-9]+(?:-[0-9]+)?"'))
  hits <- sub('^protein_name="', "", hits)
  hits <- sub('"$', "", hits)
  split_accessions(hits)
}
