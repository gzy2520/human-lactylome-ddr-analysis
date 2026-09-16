# Shared helpers for the 31-group expression extraction (2026-09-16).
#
# Conventions enforced here:
#   * every matrix reaching the final output is keyed by a stable gene identifier
#     (Ensembl ENSG for the merged output); gene symbols are never a join key;
#   * TPM is computed from counts against one shared merged-exon length table, so
#     count-based groups are normalised identically; FPKM/RPKM sources are converted
#     with the same length table and are flagged as converted;
#   * ambiguous stable-ID collapses are resolved by an explicit, recorded rule rather
#     than a silent join.
#
# Base R only: the server runs R 4.1.2 with no data.table / Bioconductor packages.

suppressPackageStartupMessages(library(parallel))

## ---- reading -------------------------------------------------------------

# R detects gzip compression from the file signature, so a .gz path is read directly:
# handing read.delim an explicit gzfile() connection leaves it closed on return, and
# isOpen() on that closed connection then errors rather than reporting FALSE.
# quote defaults to "" because the expression matrices are unquoted and a stray quote
# character inside a field must not swallow the rest of the line; files that really are
# quoted (the GDC manifest, author CSVs) pass quote = "\"" explicitly.
read_table_gz <- function(path, sep = "\t", header = TRUE, quote = "", ...) {
  read.delim(path, sep = sep, header = header, check.names = FALSE,
             quote = quote, comment.char = "", ...)
}

## ---- identifiers ---------------------------------------------------------

# "ENSG00000000003.15" -> "ENSG00000000003"; "ENSG00000182378.14_PAR_Y" -> "ENSG00000182378_PAR_Y"
strip_ensembl_version <- function(x) sub("\\.[0-9]+(_PAR_Y)?$", "\\1", x)

is_ensembl <- function(x) grepl("^ENSG[0-9]+(_PAR_Y)?$", x)
is_entrez  <- function(x) grepl("^[0-9]+$", x)
is_ensTranscript <- function(x) grepl("^ENST[0-9]+$", x)

load_id_mapping <- function(path) {
  m <- read_table_gz(path)
  names(m)[1] <- "tax_id"
  m <- m[m$tax_id == 9606, c("GeneID", "Ensembl_gene_identifier", "Ensembl_rna_identifier")]
  m$GeneID <- as.character(m$GeneID)
  m
}

# Entrez -> ENSG through the official NCBI gene2ensembl table.
# one Entrez -> one ENSG : kept
# Entrez -> many ENSG      : dropped (ambiguous), counted
# many Entrez -> one ENSG  : kept as a group, collapsed by `collapse`
build_entrez_map <- function(mapping) {
  m <- unique(mapping[, c("GeneID", "Ensembl_gene_identifier")])
  m <- m[!is.na(m$Ensembl_gene_identifier) & m$Ensembl_gene_identifier != "", ]
  n_per_entrez <- table(m$GeneID)
  ambiguous <- names(n_per_entrez)[n_per_entrez > 1L]
  kept <- m[!m$GeneID %in% ambiguous, ]
  list(map = kept, dropped_ambiguous_entrez = length(ambiguous))
}

# ENST -> ENSG through the same table. The table stores versioned transcript IDs
# ("ENST00000263100.8") while source matrices key on the bare accession, so versions are
# stripped on both sides before matching.
build_transcript_map <- function(mapping) {
  m <- mapping[, c("Ensembl_rna_identifier", "Ensembl_gene_identifier")]
  m$Ensembl_rna_identifier <- sub("\\.[0-9]+$", "", m$Ensembl_rna_identifier)
  m <- m[!is.na(m$Ensembl_rna_identifier) & m$Ensembl_rna_identifier != "" & !is.na(m$Ensembl_gene_identifier), ]
  unique(m[!duplicated(m$Ensembl_rna_identifier), ])
}

collapse_by_stable_id <- function(mat, map_from, map_to, value_type = c("counts", "fraction")) {
  value_type <- match.arg(value_type)
  keep <- map_from %in% rownames(mat)
  map_from <- map_from[keep]; map_to <- map_to[keep]
  sub <- mat[map_from, , drop = FALSE]
  if (value_type == "counts") {
    # counts are additive: many source IDs collapsing onto one ENSG are summed
    out <- rowsum(sub, group = map_to, reorder = TRUE)
    list(matrix = out, rule = "summed where several source IDs map to one Ensembl gene",
         n_source = nrow(sub), n_out = nrow(out))
  } else {
    # fraction values (TPM/FPKM/RPKM) are not additive or averageable across
    # distinct genes: any Ensembl gene receiving more than one source ID is dropped
    tab <- table(map_to)
    ambiguous <- names(tab)[tab > 1L]
    sel <- !map_to %in% ambiguous
    out <- sub[sel, , drop = FALSE]
    rownames(out) <- map_to[sel]
    list(matrix = out, rule = "ambiguous Ensembl genes dropped",
         n_source = nrow(sub), n_out = nrow(out), n_dropped_ambiguous = length(ambiguous))
  }
}

dedupe_ensembl_rows <- function(mat, value_type = c("counts", "fraction")) {
  value_type <- match.arg(value_type)
  ids <- strip_ensembl_version(rownames(mat))
  if (!anyDuplicated(ids)) { rownames(mat) <- ids; return(list(matrix = mat, rule = "none needed")) }
  if (value_type == "counts") {
    out <- rowsum(mat, group = ids, reorder = TRUE)
    list(matrix = out, rule = "rows with identical version-stripped Ensembl ID were summed")
  } else {
    dup <- unique(ids[duplicated(ids)])
    sel <- !ids %in% dup
    out <- mat[sel, , drop = FALSE]; rownames(out) <- ids[sel]
    list(matrix = out, rule = "rows with identical version-stripped Ensembl ID were dropped")
  }
}

## ---- expression values ---------------------------------------------------

tpm_from_counts <- function(counts, lengths_bp) {
  common <- intersect(rownames(counts), names(lengths_bp))
  stopifnot(length(common) > 1000L)
  len_kb <- lengths_bp[common] / 1000
  rate <- sweep(counts[common, , drop = FALSE], 1L, len_kb, "/")
  rate / rep(colSums(rate), each = nrow(rate)) * 1e6
}

# FPKM and RPKM are the same quantity on different scales; both convert to TPM by
# dividing by gene length and re-scaling to a 1e6 sum.
fpkm_to_tpm <- function(fpkm, lengths_bp) {
  common <- intersect(rownames(fpkm), names(lengths_bp))
  stopifnot(length(common) > 1000L)
  len_kb <- lengths_bp[common] / 1000
  rate <- sweep(fpkm[common, , drop = FALSE], 1L, len_kb, "/")
  rate / rep(colSums(rate), each = nrow(rate)) * 1e6
}

log2_tpm <- function(tpm) log2(tpm + 0.5)

## ---- QC ------------------------------------------------------------------

group_qc <- function(counts, tpm, sample_meta) {
  l2 <- log2_tpm(tpm)
  filt <- l2[rowSums(tpm >= 1) >= max(2L, ceiling(ncol(tpm) / 2L)), , drop = FALSE]
  if (nrow(filt) < 2L) filt <- l2
  cm <- suppressWarnings(cor(filt, method = "spearman"))
  off <- cm[upper.tri(cm)]
  data.frame(
    SampleID = colnames(tpm),
    GenesInMatrix = nrow(tpm),
    GenesDetected = if (is.null(counts)) NA_integer_ else colSums(counts > 0),
    LibrarySize = if (is.null(counts)) NA_real_ else colSums(counts),
    TPMTotal = colSums(tpm),
    ZeroFraction = colMeans(tpm == 0),
    PairwiseSpearmanMin = if (length(off)) min(off) else NA_real_,
    PairwiseSpearmanMax = if (length(off)) max(off) else NA_real_,
    stringsAsFactors = FALSE
  )
}

write_matrix_tsv <- function(mat, path) {
  out <- data.frame(stable_gene_id = rownames(mat), mat, check.names = FALSE)
  con <- gzfile(path, "wt")
  on.exit(close(con))
  write.table(out, con, sep = "\t", quote = FALSE, row.names = FALSE)
}

finish_group <- function(group_ids, reference_key, label, provenance, counts, tpm, sample_meta) {
  qc <- group_qc(counts, tpm, sample_meta)
  list(
    GroupIDs = group_ids, ReferenceKey = reference_key, GroupLabel = label,
    Provenance = provenance, SampleMeta = sample_meta,
    counts = counts, tpm = tpm, log2_tpm = log2_tpm(tpm), qc = qc
  )
}
