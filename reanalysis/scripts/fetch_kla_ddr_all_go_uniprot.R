#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(httr2)
})

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
script_path <- if (length(file_arg) > 0L) {
  normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE)
} else {
  normalizePath(
    "reanalysis/scripts/fetch_kla_ddr_all_go_uniprot.R",
    mustWork = TRUE
  )
}
project_root <- normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = TRUE)

accession_path <- file.path(
  project_root,
  paste0(
    "reanalysis/results/tables/kla_ddr_raw_go_umap_33groups/",
    "protein_raw_go_term_binary_matrix.csv"
  )
)
output_path <- file.path(
  project_root,
  "reanalysis/config/uniprot_kla_ddr_507_all_go_2026-08-10.tsv"
)
metadata_path <- file.path(
  project_root,
  "reanalysis/config/uniprot_kla_ddr_507_all_go_2026-08-10_metadata.tsv"
)

stop_if_false <- function(condition, message) {
  if (!isTRUE(condition)) {
    stop(message, call. = FALSE)
  }
}

base_accession <- function(x) {
  x <- trimws(as.character(x))
  x <- sub("^(sp|tr)\\|", "", x)
  x <- sub("\\|.*$", "", x)
  x <- sub("^.*:", "", x)
  sub("-[0-9]+$", "", x)
}

accessions <- fread(accession_path, select = "BaseAccession")$BaseAccession
accessions <- sort(unique(base_accession(accessions)))
stop_if_false(length(accessions) == 507L, "Expected 507 unique BaseAccessions.")

api_url <- "https://rest.uniprot.org/uniprotkb/search"
field_spec <- paste(
  c(
    "accession",
    "reviewed",
    "id",
    "protein_name",
    "gene_primary",
    "go_id",
    "go_p",
    "go_c",
    "go_f",
    "protein_families",
    "ft_domain"
  ),
  collapse = ","
)

chunks <- split(accessions, ceiling(seq_along(accessions) / 40L))
responses <- vector("list", length(chunks))
release_headers <- vector("list", length(chunks))

for (chunk_index in seq_along(chunks)) {
  chunk <- chunks[[chunk_index]]
  query <- paste0(
    "(",
    paste(sprintf("accession:%s", chunk), collapse = " OR "),
    ") AND organism_id:9606"
  )
  response <- request(api_url) |>
    req_url_query(
      query = query,
      format = "tsv",
      fields = field_spec,
      size = 500
    ) |>
    req_user_agent("kla-ddr-all-go-layout/1.0") |>
    req_retry(max_tries = 4L) |>
    req_perform()

  stop_if_false(
    resp_status(response) == 200L,
    sprintf("UniProt request failed for chunk %d.", chunk_index)
  )
  responses[[chunk_index]] <- fread(
    text = resp_body_string(response),
    sep = "\t",
    quote = "",
    na.strings = character(),
    check.names = FALSE
  )
  release_headers[[chunk_index]] <- data.table(
    Chunk = chunk_index,
    UniProtRelease = resp_header(response, "x-uniprot-release"),
    UniProtReleaseDate = resp_header(response, "x-uniprot-release-date"),
    APIDeploymentDate = resp_header(response, "x-api-deployment-date")
  )
}

annotation <- rbindlist(responses, use.names = TRUE, fill = TRUE)
setnames(annotation, "Entry", "BaseAccession")
annotation[, BaseAccession := base_accession(BaseAccession)]
setorder(annotation, BaseAccession)

stop_if_false(
  nrow(annotation) == 507L && uniqueN(annotation$BaseAccession) == 507L,
  "UniProt did not return exactly 507 unique human protein entries."
)
stop_if_false(
  setequal(annotation$BaseAccession, accessions),
  paste0(
    "Returned UniProt accessions do not exactly match the 507 requested ",
    "BaseAccessions."
  )
)

release_table <- unique(rbindlist(release_headers))
stop_if_false(
  uniqueN(release_table$UniProtRelease) == 1L,
  "UniProt release changed between annotation chunks."
)
stop_if_false(
  uniqueN(release_table$UniProtReleaseDate) == 1L,
  "UniProt release date changed between annotation chunks."
)

fwrite(annotation, output_path, sep = "\t", quote = TRUE, na = "")
metadata <- data.table(
  Key = c(
    "RetrievalDate",
    "RetrievalTimezone",
    "UniProtRelease",
    "UniProtReleaseDate",
    "APIDeploymentDate",
    "OrganismTaxonID",
    "RequestedProteinCount",
    "ReturnedProteinCount",
    "QueryEndpoint",
    "RequestedFields",
    "AccessionSource",
    "OutputMD5"
  ),
  Value = c(
    format(Sys.Date(), "%Y-%m-%d"),
    "Asia/Shanghai",
    unique(release_table$UniProtRelease),
    unique(release_table$UniProtReleaseDate),
    paste(sort(unique(release_table$APIDeploymentDate)), collapse = ";"),
    "9606",
    length(accessions),
    nrow(annotation),
    api_url,
    field_spec,
    sub(paste0("^", project_root, "/"), "", accession_path),
    unname(tools::md5sum(output_path))
  )
)
fwrite(metadata, metadata_path, sep = "\t", quote = TRUE)

message("UniProt entries retrieved: ", nrow(annotation))
message("UniProt release: ", unique(release_table$UniProtRelease))
message("Annotation output: ", output_path)
message("Metadata output: ", metadata_path)
