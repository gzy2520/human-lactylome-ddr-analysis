#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(biomaRt)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args)) normalizePath(args[[1]]) else normalizePath(".")
output_path <- file.path(
  project_root,
  "reanalysis", "config", "ensembl_protein_to_uniprot_biomart.tsv"
)
unmapped_path <- file.path(
  project_root,
  "reanalysis", "config", "ensembl_protein_unmapped_biomart.tsv"
)

files <- list.files(
  file.path(
    project_root,
    "data", "PXD010154", "search_results", "extracted_healthy_tissues"
  ),
  pattern = "proteinGroups\\.txt$",
  recursive = TRUE,
  full.names = TRUE
)

ensembl_ids <- character()
for (path in files) {
  data <- read.delim(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = ""
  )
  keep <- (is.na(data$Reverse) | data$Reverse != "+") &
    (is.na(data$`Potential contaminant`) |
       data$`Potential contaminant` != "+") &
    (is.na(data$`Only identified by site`) |
       data$`Only identified by site` != "+")
  values <- unlist(strsplit(
    as.character(data$`Majority protein IDs`[keep]),
    "[;,]"
  ))
  values <- sub("_RNA$", "", trimws(values))
  values <- sub("\\.[0-9]+$", "", values)
  ensembl_ids <- c(
    ensembl_ids,
    values[grepl("^ENSP[0-9]+$", values)]
  )
}
ensembl_ids <- sort(unique(ensembl_ids))
if (!length(ensembl_ids)) stop("No Ensembl protein IDs found")

connect_mart <- function() {
  mirrors <- c("www", "useast", "asia")
  for (mirror in mirrors) {
    mart <- tryCatch(
      useEnsembl(
        biomart = "genes",
        dataset = "hsapiens_gene_ensembl",
        mirror = mirror
      ),
      error = function(e) NULL
    )
    if (!is.null(mart)) return(mart)
  }
  stop("Unable to connect to an Ensembl BioMart mirror")
}

mart <- connect_mart()
chunks <- split(ensembl_ids, ceiling(seq_along(ensembl_ids) / 2000))
mapped_parts <- vector("list", length(chunks))

for (i in seq_along(chunks)) {
  values <- chunks[[i]]
  result <- NULL
  for (attempt in seq_len(3)) {
    result <- tryCatch(
      getBM(
        attributes = c(
          "ensembl_peptide_id",
          "uniprotswissprot",
          "uniprotsptrembl"
        ),
        filters = "ensembl_peptide_id",
        values = values,
        mart = mart
      ),
      error = function(e) NULL
    )
    if (!is.null(result)) break
    Sys.sleep(attempt * 2)
    mart <- connect_mart()
  }
  if (is.null(result)) {
    stop("BioMart mapping failed for chunk ", i, " of ", length(chunks))
  }
  mapped_parts[[i]] <- result
  message(
    "Mapped BioMart chunk ", i, "/", length(chunks),
    " (", length(values), " Ensembl protein IDs)"
  )
}

raw_mapping <- bind_rows(mapped_parts)
mapping <- bind_rows(
  raw_mapping |>
    transmute(
      EnsemblProteinID = ensembl_peptide_id,
      BaseAccession = uniprotswissprot,
      UniProtSource = "Swiss-Prot"
    ),
  raw_mapping |>
    transmute(
      EnsemblProteinID = ensembl_peptide_id,
      BaseAccession = uniprotsptrembl,
      UniProtSource = "TrEMBL"
    )
) |>
  mutate(
    BaseAccession = sub("-[0-9]+$", "", trimws(BaseAccession)),
    MappingSource = "Ensembl BioMart hsapiens_gene_ensembl",
    MappingDate = format(Sys.Date(), "%Y-%m-%d")
  ) |>
  filter(
    nzchar(EnsemblProteinID),
    nzchar(BaseAccession)
  ) |>
  distinct() |>
  arrange(EnsemblProteinID, UniProtSource, BaseAccession)

unmapped <- data.frame(
  EnsemblProteinID = setdiff(ensembl_ids, unique(mapping$EnsemblProteinID)),
  MappingSource = "Ensembl BioMart hsapiens_gene_ensembl",
  MappingDate = format(Sys.Date(), "%Y-%m-%d"),
  stringsAsFactors = FALSE
)

write.table(
  mapping,
  output_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  unmapped,
  unmapped_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message(
  "BioMart mapping complete: ",
  length(unique(mapping$EnsemblProteinID)), "/", length(ensembl_ids),
  " Ensembl protein IDs mapped; ", nrow(unmapped), " unmapped."
)
