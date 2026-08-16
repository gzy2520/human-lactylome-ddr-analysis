#!/usr/bin/env Rscript

# Build the current 30-group / 399-protein pathway scores from direct GO terms.
# Gene symbols and protein names are retained for display/audit only; the stable
# isoform-stripped UniProt BaseAccession is the analytical key.

suppressPackageStartupMessages({
  library(data.table)
  library(digest)
  library(GO.db)
})

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args)) normalizePath(args[[1]]) else normalizePath(".")
analysis_name <- "go_term_pathway_scoring_30groups"
output_dir <- file.path(project_root, "results", "tables", analysis_name)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

membership_path <- file.path(
  project_root,
  "results/tables/four_class_venn/kla_ddr_four_class_venn/membership.csv"
)
uniprot_path <- file.path(
  project_root,
  "config/uniprot_kla_ddr_507_all_go_2026-08-10.tsv"
)
uniprot_metadata_path <- file.path(
  project_root,
  "config/uniprot_kla_ddr_507_all_go_2026-08-10_metadata.tsv"
)
rules_path <- file.path(project_root, "config/go_term_pathway_seed_rules.csv")
required_inputs <- c(
  membership_path,
  uniprot_path,
  uniprot_metadata_path,
  rules_path
)

assert <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

assert(
  all(file.exists(required_inputs)),
  paste(
    "Missing required input(s):",
    paste(required_inputs[!file.exists(required_inputs)], collapse = "; ")
  )
)

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

membership <- fread(membership_path)
membership[, BaseAccession := base_accession(BaseAccession)]
assert(
  nrow(membership) == 399L &&
    uniqueN(membership$BaseAccession) == 399L,
  "Expected 399 unique current Kla-DDR BaseAccessions."
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
    all(membership$BaseAccession %in% uniprot$BaseAccession),
  "The UniProt cache must contain the current 399 proteins within 507 unique proteins."
)
uniprot <- uniprot[BaseAccession %in% membership$BaseAccession]

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
  "The UniProt cache lacks one or more direct GO aspect columns."
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
    UniProtGOName = sub("\\s*\\[GO:[0-9]{7}\\]$", "", pieces)
  )
}

direct_go <- rbindlist(lapply(seq_len(nrow(aspect_spec)), function(index) {
  ontology <- aspect_spec$Ontology[[index]]
  source_column <- aspect_spec$SourceColumn[[index]]
  rbindlist(lapply(seq_len(nrow(uniprot)), function(protein_index) {
    parse_go_cell(
      uniprot$BaseAccession[[protein_index]],
      ontology,
      uniprot[[source_column]][[protein_index]]
    )
  }))
}))
direct_go <- unique(
  direct_go,
  by = c("BaseAccession", "Ontology", "GO_ID")
)
setorder(direct_go, BaseAccession, Ontology, GO_ID)
assert(
  nrow(direct_go) == 10605L &&
    uniqueN(direct_go$GO_ID) == 2785L,
  "Expected 10,605 direct protein-GO pairs and 2,785 unique direct GO terms."
)
assert(
  uniqueN(direct_go[Ontology == "BP", BaseAccession]) == 399L,
  "At least one current protein lacks a direct BP annotation."
)

term_catalog <- unique(
  direct_go[, .(GO_ID, Ontology, UniProtGOName)],
  by = "GO_ID"
)
term_catalog[, `:=`(
  GODbName = vapply(
    GO_ID,
    function(id) {
      value <- tryCatch(Term(GOTERM[[id]]), error = function(e) NA_character_)
      if (is.null(value) || !length(value)) NA_character_ else as.character(value)
    },
    character(1)
  ),
  GODbOntology = vapply(
    GO_ID,
    function(id) {
      value <- tryCatch(Ontology(GOTERM[[id]]), error = function(e) NA_character_)
      if (is.null(value) || !length(value)) NA_character_ else as.character(value)
    },
    character(1)
  ),
  GODbDefinition = vapply(
    GO_ID,
    function(id) {
      value <- tryCatch(Definition(GOTERM[[id]]), error = function(e) NA_character_)
      if (is.null(value) || !length(value)) NA_character_ else as.character(value)
    },
    character(1)
  )
)]
term_catalog[, MissingFromGODb := is.na(GODbName)]
term_catalog[, NameMismatch :=
  !MissingFromGODb & tolower(UniProtGOName) != tolower(GODbName)]
assert(
  term_catalog[MissingFromGODb == TRUE, .N] == 20L,
  "Expected 20 recent UniProt GO terms to be absent from GO.db 3.20.0."
)

rules <- fread(rules_path)
rules[, IncludeDescendants :=
  IncludeDescendants %in% c(TRUE, "TRUE", "True", 1, "1")]
required_rule_columns <- c(
  "Pathway",
  "PathwayOrder",
  "SeedGOID",
  "IncludeDescendants",
  "RuleLabel",
  "Rationale"
)
assert(
  all(required_rule_columns %in% names(rules)) &&
    setequal(unique(rules$Pathway), c("BER", "NER", "MMR", "FA", "HR", "NHEJ", "AEJ")) &&
    !anyDuplicated(rules[, .(Pathway, SeedGOID)]),
  "The seed-rule table must define seven pathways with unique pathway-seed pairs."
)

go_descendants <- function(go_id) {
  ontology <- tryCatch(Ontology(GOTERM[[go_id]]), error = function(e) NA_character_)
  assert(!is.na(ontology), paste("Seed term is absent from GO.db:", go_id))
  offspring_environment <- switch(
    ontology,
    BP = GOBPOFFSPRING,
    CC = GOCCOFFSPRING,
    MF = GOMFOFFSPRING,
    stop("Unsupported ontology for seed term: ", go_id)
  )
  offspring <- tryCatch(offspring_environment[[go_id]], error = function(e) character())
  unique(c(go_id, offspring))
}

rule_hits <- rbindlist(lapply(seq_len(nrow(rules)), function(index) {
  rule <- rules[index]
  candidate_ids <- if (rule$IncludeDescendants) {
    go_descendants(rule$SeedGOID)
  } else {
    rule$SeedGOID
  }
  matched_ids <- intersect(term_catalog$GO_ID, candidate_ids)
  if (!length(matched_ids)) return(NULL)
  data.table(
    GO_ID = matched_ids,
    Pathway = rule$Pathway,
    PathwayOrder = rule$PathwayOrder,
    SeedGOID = rule$SeedGOID,
    IncludeDescendants = rule$IncludeDescendants,
    MatchType = ifelse(matched_ids == rule$SeedGOID, "seed_exact", "GO_descendant"),
    RuleLabel = rule$RuleLabel,
    Rationale = rule$Rationale
  )
}))
rule_hits <- unique(
  rule_hits,
  by = c("GO_ID", "Pathway", "SeedGOID", "MatchType")
)
setorder(rule_hits, PathwayOrder, GO_ID, SeedGOID)

pathway_mapping <- rule_hits[
  ,
  .(
    SeedGOIDs = paste(sort(unique(SeedGOID)), collapse = ";"),
    MatchTypes = paste(sort(unique(MatchType)), collapse = ";"),
    RuleLabels = paste(sort(unique(RuleLabel)), collapse = "; "),
    MappingRationale = paste(sort(unique(Rationale)), collapse = " | ")
  ),
  by = .(GO_ID, Pathway, PathwayOrder)
]
pathway_mapping <- merge(
  pathway_mapping,
  term_catalog,
  by = "GO_ID",
  all.x = TRUE,
  sort = FALSE
)
setorder(pathway_mapping, PathwayOrder, GO_ID)

mapped_term_ids <- unique(pathway_mapping$GO_ID)
others_mapping <- copy(term_catalog[!GO_ID %in% mapped_term_ids])
others_mapping[, `:=`(
  Pathway = "Others",
  PathwayOrder = 8L,
  SeedGOIDs = NA_character_,
  MatchTypes = "not_specific_to_seven_pathways",
  RuleLabels = "Others",
  MappingRationale = paste(
    "No curated seven-pathway seed or GO descendant matched;",
    "broad DDR/repair/binding terms remain Others."
  )
)]
setcolorder(others_mapping, names(pathway_mapping))
term_pathway_long <- rbindlist(
  list(pathway_mapping, others_mapping),
  use.names = TRUE
)
setorder(term_pathway_long, PathwayOrder, GO_ID)

term_decisions <- term_pathway_long[
  ,
  .(
    UniProtGOName = data.table::first(UniProtGOName),
    Ontology = data.table::first(Ontology),
    GODbName = data.table::first(GODbName),
    GODbDefinition = data.table::first(GODbDefinition),
    MissingFromGODb = data.table::first(MissingFromGODb),
    NameMismatch = data.table::first(NameMismatch),
    AssignedPathways = paste(Pathway[Pathway != "Others"], collapse = ";"),
    SevenPathwayCount = sum(Pathway != "Others"),
    Decision = ifelse(all(Pathway == "Others"), "Others", "Seven-pathway assignment")
  ),
  by = GO_ID
]
term_decisions[SevenPathwayCount == 0L, AssignedPathways := "Others"]
setorder(term_decisions, Ontology, GO_ID)
assert(
  nrow(term_decisions) == 2785L &&
    all(term_decisions$SevenPathwayCount >= 0L),
  "Every direct GO term must receive a seven-pathway assignment or Others."
)
assert(
  !any(term_pathway_long$Pathway == "Others" &
    term_pathway_long$GO_ID %in% mapped_term_ids),
  "Others must be mutually exclusive with the seven-pathway mappings."
)

multi_pathway_terms <- term_decisions[SevenPathwayCount > 1L]

protein_term_pathway <- merge(
  direct_go,
  term_pathway_long[, .(GO_ID, Pathway, PathwayOrder)],
  by = "GO_ID",
  all.x = TRUE,
  sort = FALSE
)
assert(
  !anyNA(protein_term_pathway$Pathway),
  "At least one direct protein-GO pair lacks a pathway/Others decision."
)
protein_term_pathway <- unique(
  protein_term_pathway,
  by = c("BaseAccession", "GO_ID", "Pathway")
)
setorder(protein_term_pathway, BaseAccession, PathwayOrder, GO_ID)

pathway_order <- c("BER", "NER", "MMR", "FA", "HR", "NHEJ", "AEJ", "Others")
protein_pathway_counts_long <- protein_term_pathway[
  ,
  .(DirectTermCount = uniqueN(GO_ID)),
  by = .(BaseAccession, Pathway)
]
count_grid <- CJ(
  BaseAccession = sort(membership$BaseAccession),
  Pathway = pathway_order,
  unique = TRUE
)
protein_pathway_counts_long <- merge(
  count_grid,
  protein_pathway_counts_long,
  by = c("BaseAccession", "Pathway"),
  all.x = TRUE,
  sort = FALSE
)
protein_pathway_counts_long[is.na(DirectTermCount), DirectTermCount := 0L]
protein_pathway_counts_long[
  ,
  PathwayOrder := match(Pathway, pathway_order)
]
setorder(protein_pathway_counts_long, BaseAccession, PathwayOrder)

protein_pathway_counts <- dcast(
  protein_pathway_counts_long,
  BaseAccession ~ Pathway,
  value.var = "DirectTermCount",
  fill = 0L
)
setcolorder(protein_pathway_counts, c("BaseAccession", pathway_order))
seven_pathways <- pathway_order[seq_len(7L)]
protein_pathway_counts[, SevenPathwayTermScore :=
  rowSums(.SD), .SDcols = seven_pathways]
protein_pathway_counts[, DirectGOCount :=
  direct_go[, .N, by = BaseAccession]$N[
    match(BaseAccession, direct_go[, .N, by = BaseAccession]$BaseAccession)
  ]]

protein_pathway_binary <- copy(protein_pathway_counts)
for (pathway in seven_pathways) {
  set(
    protein_pathway_binary,
    j = pathway,
    value = as.integer(protein_pathway_binary[[pathway]] > 0L)
  )
}
protein_pathway_binary[, c("Others", "SevenPathwayTermScore", "DirectGOCount") := NULL]

display_metadata <- membership[
  ,
  .(
    BaseAccession,
    GeneSymbolAudit,
    ProteinNameAudit,
    ReviewedStatus
  )
]
protein_pathway_counts <- merge(
  display_metadata,
  protein_pathway_counts,
  by = "BaseAccession",
  all.y = TRUE,
  sort = FALSE
)
protein_pathway_binary <- merge(
  display_metadata,
  protein_pathway_binary,
  by = "BaseAccession",
  all.y = TRUE,
  sort = FALSE
)
setorder(protein_pathway_counts, BaseAccession)
setorder(protein_pathway_binary, BaseAccession)

missing_go_terms <- term_catalog[MissingFromGODb == TRUE][order(Ontology, GO_ID)]
assert(
  all(missing_go_terms$GO_ID %in% others_mapping$GO_ID),
  "A GO.db-missing term was assigned without an exact curated seed."
)

input_audit <- data.table(
  InputRole = c(
    "current 30-group Kla-DDR membership",
    "UniProt full direct GO annotation cache",
    "UniProt retrieval metadata",
    "curated GO-term pathway seed rules"
  ),
  Path = vapply(required_inputs, relative_path, character(1)),
  SHA256 = vapply(
    required_inputs,
    function(path) digest(file = path, algo = "sha256", serialize = FALSE),
    character(1)
  )
)

ontology_summary <- direct_go[
  ,
  .(
    DirectProteinTermPairs = .N,
    UniqueTerms = uniqueN(GO_ID),
    Proteins = uniqueN(BaseAccession)
  ),
  by = Ontology
][order(match(Ontology, c("BP", "CC", "MF")))]

fwrite(direct_go, file.path(output_dir, "direct_go_annotations_399proteins.csv"))
fwrite(term_catalog[order(Ontology, GO_ID)], file.path(
  output_dir,
  "direct_go_term_catalog_2785.csv"
))
fwrite(rule_hits, file.path(output_dir, "term_pathway_rule_hit_evidence.csv"))
fwrite(term_pathway_long, file.path(output_dir, "go_term_to_pathway_long.csv"))
fwrite(term_decisions, file.path(output_dir, "go_term_decision_audit_2785.csv"))
fwrite(multi_pathway_terms, file.path(output_dir, "multi_pathway_go_terms.csv"))
fwrite(missing_go_terms, file.path(output_dir, "go_db_missing_terms.csv"))
fwrite(protein_term_pathway, file.path(output_dir, "protein_go_term_pathway_long.csv"))
fwrite(protein_pathway_counts_long, file.path(
  output_dir,
  "protein_pathway_direct_term_counts_long.csv"
))
fwrite(protein_pathway_counts, file.path(
  output_dir,
  "protein_pathway_direct_term_count_matrix.csv"
))
fwrite(protein_pathway_binary, file.path(
  output_dir,
  "protein_seven_pathway_binary_matrix.csv"
))
fwrite(ontology_summary, file.path(output_dir, "direct_go_ontology_summary.csv"))
fwrite(input_audit, file.path(output_dir, "input_file_audit.csv"))
writeLines(
  c(
    paste0("GO.db_version=", as.character(packageVersion("GO.db"))),
    paste0("GO_db_schema=", GO_dbInfo()[["DBSCHEMA"]]),
    paste0("GO_db_source_date=", GO_dbInfo()[["GOSOURCEDATE"]]),
    "",
    capture.output(sessionInfo())
  ),
  file.path(output_dir, "session_info.txt"),
  useBytes = TRUE
)

message(
  "Built GO-term pathway scores for 399 proteins: ",
  format(nrow(direct_go), big.mark = ","),
  " direct protein-GO pairs, ",
  format(nrow(term_catalog), big.mark = ","),
  " terms, ",
  format(length(mapped_term_ids), big.mark = ","),
  " seven-pathway terms, and ",
  format(nrow(multi_pathway_terms), big.mark = ","),
  " multi-pathway terms."
)
