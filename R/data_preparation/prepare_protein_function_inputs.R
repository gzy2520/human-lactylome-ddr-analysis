#!/usr/bin/env Rscript

# Build the stable 507-protein annotation and pathway-assignment inputs used by
# all publication figures. This step does not fit an embedding or draw a figure.

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(readxl)
})

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args)) normalizePath(args[[1]]) else normalizePath(".")
output_dir <- file.path(project_root, "results", "tables", "protein_function_inputs")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

score_path <- file.path(
  project_root, "data", "identifier", "260810乳酸化DDR基因评分表.xlsx"
)
uniprot_path <- file.path(
  project_root, "config", "uniprot_kla_ddr_507_all_go_2026-08-10.tsv"
)
uniprot_metadata_path <- file.path(
  project_root, "config", "uniprot_kla_ddr_507_all_go_2026-08-10_metadata.tsv"
)
required_inputs <- c(score_path, uniprot_path, uniprot_metadata_path)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs)) {
  stop("Missing required input(s): ", paste(missing_inputs, collapse = "; "))
}

assert <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

base_accession <- function(values) {
  values <- trimws(as.character(values))
  values <- sub("^(sp|tr)\\|", "", values)
  values <- sub("\\|.*$", "", values)
  values <- sub("^.*:", "", values)
  sub("-[0-9]+$", "", values)
}

relative_path <- function(path) {
  sub(paste0("^", project_root, "/?"), "", path)
}

pathway_info <- data.table(
  Pathway = c(
    "HR", "NHEJ", "AEJ", "BER", "NER", "MMR", "FA",
    "Chromatin Interaction",
    "Others (Transcription, RNA processing and proteostasis)"
  ),
  DisplayLabel = c(
    "HR", "NHEJ", "AEJ", "BER", "NER", "MMR", "FA",
    "Chromatin interaction", "Other support"
  ),
  Color = c(
    "#3C5488", "#E64B35", "#00A087", "#4DBBD5", "#F39B7F",
    "#8491B4", "#91D1C2", "#7E6148", "#6F6F6F"
  ),
  PathwayOrder = seq_len(9L)
)
pathway_info[, ColorName := c(
  "Chambray", "Cinnabar", "PersianGreen", "Shakespeare", "Apricot",
  "WildBlueYonder", "MonteCarlo", "RomanCoffee", "NeutralCharcoal"
)]
setcolorder(
  pathway_info,
  c("Pathway", "DisplayLabel", "ColorName", "Color", "PathwayOrder")
)
pathway_columns <- pathway_info$Pathway

scores <- as.data.table(
  read_excel(score_path, sheet = "评分表", .name_repair = "minimal")
)
scores <- scores[!is.na(BaseAccession) & nzchar(trimws(BaseAccession))]
scores[, BaseAccession := base_accession(BaseAccession)]
scores[, ID := as.integer(ID)]
for (column in pathway_columns) {
  assert(column %in% names(scores), paste("Missing pathway column:", column))
  set(scores, j = column, value = as.integer(scores[[column]]))
}
assert(
  nrow(scores) == 507L && uniqueN(scores$BaseAccession) == 507L,
  "The score workbook must contain 507 unique BaseAccessions."
)
score_matrix <- as.matrix(scores[, ..pathway_columns])
assert(
  !anyNA(score_matrix) && all(score_matrix %in% c(-1L, 0L, 1L)),
  "Pathway scores must be complete and limited to -1, 0, and +1."
)

uniprot <- fread(
  uniprot_path,
  sep = "\t",
  quote = "\"",
  na.strings = character(),
  check.names = FALSE
)
uniprot[, BaseAccession := base_accession(BaseAccession)]
assert(
  nrow(uniprot) == 507L &&
    uniqueN(uniprot$BaseAccession) == 507L &&
    setequal(uniprot$BaseAccession, scores$BaseAccession),
  "The UniProt cache and score workbook must contain the same 507 proteins."
)

aspect_spec <- data.table(
  Ontology = c("BP", "CC", "MF"),
  SourceColumn = c(
    "Gene Ontology (biological process)",
    "Gene Ontology (cellular component)",
    "Gene Ontology (molecular function)"
  )
)
assert(
  all(aspect_spec$SourceColumn %in% names(uniprot)),
  "The UniProt cache lacks one or more GO aspect columns."
)

parse_go_cell <- function(accession, ontology, value) {
  if (is.na(value) || !nzchar(value)) return(NULL)
  pieces <- strsplit(value, "; ", fixed = TRUE)[[1L]]
  pieces <- pieces[grepl("\\[GO:[0-9]{7}\\]$", pieces)]
  if (!length(pieces)) return(NULL)
  data.table(
    BaseAccession = accession,
    Ontology = ontology,
    GO_ID = sub("^.*\\[(GO:[0-9]{7})\\]$", "\\1", pieces),
    GO_Term = sub("\\s*\\[GO:[0-9]{7}\\]$", "", pieces)
  )
}

go_long <- rbindlist(lapply(seq_len(nrow(aspect_spec)), function(aspect_index) {
  ontology <- aspect_spec$Ontology[[aspect_index]]
  source_column <- aspect_spec$SourceColumn[[aspect_index]]
  rbindlist(lapply(seq_len(nrow(uniprot)), function(protein_index) {
    parse_go_cell(
      uniprot$BaseAccession[[protein_index]],
      ontology,
      uniprot[[source_column]][[protein_index]]
    )
  }))
}))
go_long <- unique(go_long, by = c("BaseAccession", "Ontology", "GO_ID"))
setorder(go_long, BaseAccession, Ontology, GO_ID)
assert(
  nrow(go_long) == 13738L && uniqueN(go_long$GO_ID) == 3461L,
  "Expected 13,738 direct protein-GO annotations and 3,461 GO terms."
)
assert(
  uniqueN(go_long[Ontology == "BP", BaseAccession]) == 507L,
  "At least one protein lacks a direct biological-process annotation."
)

protein_levels <- sort(scores$BaseAccession)
term_levels <- sort(unique(go_long$GO_ID))
go_binary <- sparseMatrix(
  i = match(go_long$BaseAccession, protein_levels),
  j = match(go_long$GO_ID, term_levels),
  x = 1L,
  dims = c(length(protein_levels), length(term_levels)),
  dimnames = list(protein_levels, term_levels)
)
go_binary@x[] <- 1
assert(
  identical(dim(go_binary), c(507L, 3461L)) && sum(go_binary) == 13738L,
  "The direct-GO binary matrix contract was not recovered."
)

protein_metadata <- go_long[, .(
  DirectGOCount = .N,
  BPTermCount = sum(Ontology == "BP"),
  CCTermCount = sum(Ontology == "CC"),
  MFTermCount = sum(Ontology == "MF")
), by = BaseAccession]
scores <- scores[match(protein_levels, BaseAccession)]
protein_metadata <- merge(
  scores,
  protein_metadata,
  by = "BaseAccession",
  all.x = TRUE,
  sort = FALSE
)
protein_metadata[, PositiveAssignmentCount := rowSums(.SD == 1L), .SDcols = pathway_columns]
protein_metadata[, NegativeAssignmentCount := rowSums(.SD == -1L), .SDcols = pathway_columns]
protein_metadata[, TotalAssignmentCount :=
  PositiveAssignmentCount + NegativeAssignmentCount]
setorder(protein_metadata, BaseAccession)

assignment_long <- melt(
  scores,
  id.vars = c("ID", "BaseAccession", "GeneSymbol", "ProteinName", "Note"),
  measure.vars = pathway_columns,
  variable.name = "Pathway",
  value.name = "Score",
  variable.factor = FALSE
)[Score != 0L]
assignment_long <- merge(
  assignment_long,
  pathway_info,
  by = "Pathway",
  all.x = TRUE,
  sort = FALSE
)
assignment_long[, Direction :=
  fifelse(Score == 1L, "Promoting (+1)", "Suppressing (-1)")]
setorder(assignment_long, BaseAccession, PathwayOrder)
assert(
  nrow(assignment_long) == 1175L &&
    sum(assignment_long$Score == 1L) == 1108L &&
    sum(assignment_long$Score == -1L) == 67L &&
    protein_metadata[TotalAssignmentCount == 0L, .N] == 22L,
  "The signed pathway-assignment contract changed."
)

binary_dense <- as.data.table(as.matrix(go_binary))
binary_dense[, BaseAccession := protein_levels]
setcolorder(binary_dense, c("BaseAccession", term_levels))

input_audit <- data.table(
  InputRole = c(
    "seven-pathway score workbook",
    "UniProt full direct GO annotation cache",
    "UniProt retrieval metadata"
  ),
  Path = relative_path(required_inputs),
  SHA256 = vapply(
    required_inputs,
    function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE),
    character(1)
  )
)

fwrite(go_long, file.path(output_dir, "direct_go_annotations.csv"))
fwrite(binary_dense, file.path(output_dir, "direct_go_binary_matrix.csv"))
fwrite(protein_metadata, file.path(output_dir, "protein_metadata_and_scores.csv"))
fwrite(assignment_long, file.path(output_dir, "pathway_assignments.csv"))
fwrite(pathway_info, file.path(output_dir, "pathway_colors.csv"))
fwrite(input_audit, file.path(output_dir, "input_audit.csv"))
writeLines(
  capture.output(sessionInfo()),
  file.path(output_dir, "session_info.txt"),
  useBytes = TRUE
)

message("Prepared 507 proteins, 13,738 direct GO annotations, and 1,175 assignments.")
