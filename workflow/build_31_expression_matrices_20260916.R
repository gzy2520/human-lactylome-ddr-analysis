#!/usr/bin/env Rscript
# Stage 2 of the 31-group expression extraction (2026-09-16).
#
# For every one of the 31 proteome/Kla group rows this reads the frozen source files on the
# server, extracts the selected unperturbed samples, keys genes by a stable identifier, and
# emits per-group expression matrices. It performs no comparison between groups and no
# differential-expression test: each group is reduced to its own expression profile, which is
# what the downstream cross-tissue qsmooth step expects.
#
# Main metric per group: log2(TPM + 0.5). Cleaned counts are archived alongside wherever the
# source actually provides counts. TPM is computed from counts with one shared merged-exon
# length table (metadata/annotation/human_gene_lengths_ensembl111.tsv) so that count-based
# groups are normalised identically; sources that ship only FPKM/RPKM are converted with the
# same table and are flagged as converted in the manifest.
#
# Usage: build_31_expression_matrices_20260916.R <server_root> <contract_csv> <out_dir> [--only=KLA31_01,...]
args <- commandArgs(TRUE)
stopifnot(length(args) >= 3L)
root <- normalizePath(args[[1L]], mustWork = TRUE)
contract_path <- normalizePath(args[[2L]], mustWork = TRUE)
out_dir <- args[[3L]]
only <- NULL
if (length(args) >= 4L && grepl("^--only=", args[[4L]])) {
  only <- strsplit(sub("^--only=", "", args[[4L]]), ",", fixed = TRUE)[[1L]]
}

script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
source(file.path(script_dir, "lib_kla31_expression_20260916.R"))

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "matrices"), showWarnings = FALSE)
dir.create(file.path(out_dir, "objects"), showWarnings = FALSE)

ann_dir <- file.path(root, "metadata", "annotation")
lengths_tab <- read_table_gz(file.path(ann_dir, "human_gene_lengths_ensembl111.tsv"), header = FALSE)
names(lengths_tab) <- c("gene_id", "length")
gene_lengths <- setNames(as.numeric(lengths_tab$length), lengths_tab$gene_id)
gene_lengths <- gene_lengths[!is.na(gene_lengths) & gene_lengths > 0]

id_map <- load_id_mapping(file.path(root, "metadata", "id_mapping", "human_gene2ensembl.tsv"))
entrez_map <- build_entrez_map(id_map)
transcript_map <- build_transcript_map(id_map)
# The placenta series was quantified against a GRCh37-era transcriptome, so its ENST IDs are
# resolved preferentially through the Ensembl GRCh37 transcript definitions, with the
# gene2ensembl table kept as a fallback for transcripts the older annotation does not carry.
grch37_map_path <- file.path(ann_dir, "human_grch37_transcript_to_gene.tsv")
if (file.exists(grch37_map_path)) {
  g37 <- read_table_gz(grch37_map_path, header = FALSE)
  names(g37) <- c("Ensembl_rna_identifier", "Ensembl_gene_identifier")
  transcript_map <- rbind(g37, transcript_map)
  transcript_map <- transcript_map[!duplicated(transcript_map$Ensembl_rna_identifier), ]
}
symbol_map_tab <- read_table_gz(file.path(ann_dir, "human_symbol_to_entrez.tsv"), header = FALSE)
names(symbol_map_tab) <- c("symbol", "entrez")
symbol_map_tab$entrez <- as.character(symbol_map_tab$entrez)
# one symbol -> one GeneID only; symbols with several official GeneIDs are dropped
sym_tab <- table(symbol_map_tab$symbol)
symbol_map <- symbol_map_tab[symbol_map_tab$symbol %in% names(sym_tab)[sym_tab == 1L], ]
symbol_map <- symbol_map[!duplicated(symbol_map$symbol), ]
symbol_lookup <- setNames(symbol_map$entrez, symbol_map$symbol)

gdc_manifest <- read_table_gz(file.path(root, "config", "gdc_star_counts_manifest_20260914.tsv"))
# the manifest is written fully quoted; R's quote handling is bypassed here so the
# delimiters are removed explicitly rather than depended upon
names(gdc_manifest) <- gsub('"', "", names(gdc_manifest), fixed = TRUE)
gdc_manifest[] <- lapply(gdc_manifest, function(x) {
  if (is.character(x)) trimws(gsub('"', "", x, fixed = TRUE)) else x
})

tmp_dir <- file.path(out_dir, "tmp")
dir.create(tmp_dir, showWarnings = FALSE)

## ---- source handlers -----------------------------------------------------

h_ncbi_counts <- function(spec) {
  tab <- read_table_gz(file.path(root, spec$SourcePath))
  ids <- as.character(tab[[1L]])
  stopifnot(!anyDuplicated(ids))
  sel <- spec$selector
  missing <- setdiff(sel, names(tab))
  stopifnot(length(missing) == 0L)
  m <- as.matrix(tab[, sel, drop = FALSE])
  rownames(m) <- ids
  probe <- m[seq_len(min(500L, nrow(m))), , drop = FALSE]
  integer_native <- all(probe == floor(probe))
  mapped <- collapse_by_stable_id(m, entrez_map$map$GeneID, entrez_map$map$Ensembl_gene_identifier, "counts")
  list(counts = mapped$matrix, tpm_native = NULL, id_rule = mapped$rule,
       id_note = sprintf("%d/%d Entrez genes mapped to Ensembl", mapped$n_out, mapped$n_source),
       integer_native = integer_native)
}

h_single_ensembl_counts <- function(spec) {
  tab <- read_table_gz(file.path(root, spec$SourcePath))
  ids <- as.character(tab[[1L]])
  sel <- spec$selector
  stopifnot(all(sel %in% names(tab)))
  m <- as.matrix(tab[, sel, drop = FALSE]); rownames(m) <- ids
  d <- dedupe_ensembl_rows(m, "counts")
  list(counts = d$matrix, tpm_native = NULL, id_rule = d$rule,
       id_note = sprintf("%d Ensembl rows", nrow(d$matrix)),
       integer_native = all(m == floor(m)))
}

h_single_symbol_counts <- function(spec) {
  path <- file.path(root, spec$SourcePath)
  # Some author matrices are space-delimited with a blank leading field rather than
  # tab-delimited, so the separator is taken from the file itself.
  first_line <- readLines(path, n = 1L, warn = FALSE)
  ws <- !grepl("\t", first_line)
  tab <- read_table_gz(path, sep = if (ws) "" else "\t", row.names = if (ws) 1L else NULL)
  syms <- if (ws) rownames(tab) else as.character(tab[[1L]])
  sel <- spec$selector
  stopifnot(all(sel %in% names(tab)))
  m <- as.matrix(tab[, sel, drop = FALSE]); rownames(m) <- syms
  entrez <- unname(symbol_lookup[syms])
  keep <- !is.na(entrez)
  m <- m[keep, , drop = FALSE]; entrez <- entrez[keep]
  map_from <- entrez; map_to <- unname(setNames(entrez_map$map$Ensembl_gene_identifier,
                                                entrez_map$map$GeneID)[entrez])
  keep2 <- !is.na(map_to)
  m <- m[keep2, , drop = FALSE]; map_from <- map_from[keep2]; map_to <- map_to[keep2]
  out <- rowsum(m, group = map_to, reorder = TRUE)
  list(counts = out, tpm_native = NULL, id_rule = "symbol -> official NCBI GeneID -> Ensembl; counts summed",
       id_note = sprintf("%d/%d symbol rows mapped to Ensembl", nrow(out), length(syms)),
       integer_native = all(m == floor(m)))
}

h_single_ensembl_fpkm <- function(spec) {
  tab <- read_table_gz(file.path(root, spec$SourcePath), sep = ",", quote = "\"")
  # this author table prefixes its human rows with the assembly ("hg19_ENSG...")
  ids <- sub("^.*?(ENSG[0-9]+.*)$", "\\1", as.character(tab[[1L]]))
  sel <- spec$selector
  stopifnot(all(sel %in% names(tab)))
  m <- as.matrix(tab[, sel, drop = FALSE]); rownames(m) <- ids
  d <- dedupe_ensembl_rows(m, "fraction")
  list(counts = NULL, tpm_native = NULL, fpkm = d$matrix, id_rule = d$rule,
       id_note = sprintf("%d Ensembl rows", nrow(d$matrix)), integer_native = NA)
}

h_xlsx_fpkm <- function(spec) {
  tsv <- file.path(tmp_dir, paste0(spec$ReferenceKey, "_GSE181540_sheet.tsv"))
  if (!file.exists(tsv)) {
    st <- system2("python3", c(file.path(root, "workflow", "xlsx_sheet_to_tsv_20260916.py"),
                              file.path(root, spec$SourcePath), tsv))
    stopifnot(st == 0L)
  }
  raw <- read.delim(tsv, sep = "\t", header = FALSE, check.names = FALSE,
                    quote = "", comment.char = "", colClasses = "character")
  hdr_row <- which(apply(raw, 1L, function(r) any(r == "gene_id")))[1L]
  hdr <- as.character(unlist(raw[hdr_row, ]))
  body <- raw[(hdr_row + 1L):nrow(raw), , drop = FALSE]
  names(body) <- make.unique(hdr)
  sel <- spec$selector
  stopifnot(all(sel %in% names(body)))
  ids <- body$gene_id
  m <- as.matrix(body[, sel, drop = FALSE])
  suppressWarnings(storage.mode(m) <- "double")
  rownames(m) <- ids
  m <- m[!is.na(ids) & ids != "", , drop = FALSE]
  d <- dedupe_ensembl_rows(m, "fraction")
  counts_cols <- sub("_FPKM$", "_count", sel)
  cnt <- NULL
  if (all(counts_cols %in% names(body))) {
    cm <- as.matrix(body[, counts_cols, drop = FALSE])
    suppressWarnings(storage.mode(cm) <- "double")
    rownames(cm) <- ids
    cm <- cm[!is.na(ids) & ids != "", , drop = FALSE]
    cd <- dedupe_ensembl_rows(cm, "fraction")
    cnt <- cd$matrix[rownames(d$matrix), , drop = FALSE]
    colnames(cnt) <- colnames(d$matrix)
  }
  list(counts = cnt, tpm_native = NULL, fpkm = d$matrix, id_rule = d$rule,
       id_note = sprintf("GSE181540 workbook, %d Ensembl rows (unversioned in source)", nrow(d$matrix)),
       integer_native = NA, counts_note = "source 'counts' block is cufflinks-style fractional, not integer")
}

h_per_sample_rpkm <- function(spec) {
  tar_path <- file.path(root, spec$SourcePath)
  members <- system2("tar", c("-tf", tar_path), stdout = TRUE)
  spec_members <- members[grepl("BPH_", members)]
  stopifnot(length(spec_members) == spec$ExpectedN)
  gsms <- sub("_.*$", "", spec_members)
  mats <- lapply(seq_along(spec_members), function(i) {
    v <- read.delim(pipe(sprintf("tar -xOf %s %s | zcat", shQuote(tar_path), shQuote(spec_members[[i]]))),
                    sep = "\t", header = FALSE, colClasses = c("character", "double"))
    setNames(v[[2L]], v[[1L]])
  })
  genes <- Reduce(intersect, lapply(mats, names))
  m <- sapply(mats, function(v) v[genes])
  rownames(m) <- genes
  colnames(m) <- gsms
  d <- dedupe_ensembl_rows(m, "fraction")
  list(counts = NULL, tpm_native = NULL, fpkm = d$matrix, id_rule = d$rule,
       id_note = sprintf("%d genes common to all %d samples", nrow(d$matrix), length(mats)),
       integer_native = NA, samples = gsms)
}

h_rsem_per_sample <- function(spec) {
  tar_path <- file.path(root, spec$SourcePath)
  members <- system2("tar", c("-tf", tar_path), stdout = TRUE)
  sel <- spec$selector
  chosen <- vapply(sel, function(g) grep(paste0("^", g, "_"), members, value = TRUE)[1L], character(1L))
  stopifnot(!any(is.na(chosen)))
  res <- lapply(seq_along(chosen), function(i) {
    read.delim(pipe(sprintf("tar -xOf %s %s | zcat", shQuote(tar_path), shQuote(chosen[[i]]))),
               sep = "\t", header = TRUE, check.names = FALSE, quote = "")
  })
  genes <- Reduce(intersect, lapply(res, function(d) d$gene_id))
  cnt <- sapply(res, function(d) d$expected_count[match(genes, d$gene_id)])
  tpmv <- sapply(res, function(d) d$TPM[match(genes, d$gene_id)])
  rownames(cnt) <- genes; rownames(tpmv) <- genes
  colnames(cnt) <- sel; colnames(tpmv) <- sel
  dc <- dedupe_ensembl_rows(cnt, "counts")
  dt <- dedupe_ensembl_rows(tpmv, "fraction")
  list(counts = dc$matrix, tpm_native = dt$matrix[rownames(dc$matrix), , drop = FALSE],
       id_rule = dc$rule, id_note = sprintf("%d RSEM gene rows", nrow(dc$matrix)),
       integer_native = FALSE,
       counts_note = "RSEM expected_count is fractional by construction; TPM is computed from it with the shared length table")
}

h_annotated_star_counts <- function(spec) {
  tab <- read_table_gz(file.path(root, spec$SourcePath))
  ids <- as.character(tab[[1L]])
  drop <- grepl("^N_", ids)
  tab <- tab[!drop, , drop = FALSE]; ids <- ids[!drop]
  sel <- spec$selector
  stopifnot(all(sel %in% names(tab)))
  m <- as.matrix(tab[, sel, drop = FALSE]); rownames(m) <- ids
  d <- dedupe_ensembl_rows(m, "counts")
  list(counts = d$matrix, tpm_native = NULL, id_rule = d$rule,
       id_note = sprintf("%d Ensembl gene rows", nrow(d$matrix)),
       integer_native = all(m == floor(m)))
}

h_recount3_runs <- function(spec) {
  path <- file.path(root, spec$SourcePath)
  # the file opens with "##annotation=" / "##date.generated=" metadata lines
  lead <- readLines(path, n = 5L, warn = FALSE)
  n_skip <- sum(grepl("^##", lead))
  tab <- read_table_gz(path, skip = n_skip)
  ids <- as.character(tab[[1L]])
  # run_groups holds one pipe-separated run list per sample, e.g. "SRR1|SRR2;SRR3|SRR4"
  stopifnot(all(unlist(strsplit(spec$run_groups, "[;|]")) %in% names(tab)))
  groups <- strsplit(spec$run_groups, ";", fixed = TRUE)[[1L]]
  labels <- spec$selector
  stopifnot(length(groups) == length(labels))
  m <- sapply(groups, function(g) {
    rr <- strsplit(g, "|", fixed = TRUE)[[1L]]
    rowSums(as.matrix(tab[, rr, drop = FALSE]))
  })
  rownames(m) <- ids; colnames(m) <- labels
  d <- dedupe_ensembl_rows(m, "counts")
  list(counts = d$matrix, tpm_native = NULL, id_rule = d$rule,
       id_note = sprintf("%d recount3 G026 gene rows", nrow(d$matrix)),
       integer_native = TRUE,
       counts_note = paste("recount3 sra.gene_sums stores base-level coverage sums rather than",
                           "read counts, so the archived counts block is on a coverage scale",
                           "(about 1e9 per sample here); it is proportional to read count, so the",
                           "derived TPM is unaffected. Four sequencing runs were summed per sample."))
}

h_transcript_counts <- function(spec) {
  tab <- read_table_gz(file.path(root, spec$SourcePath), header = FALSE)
  cols <- as.character(unlist(tab[1L, ]))[-1L]
  m <- as.matrix(tab[-1L, -1L, drop = FALSE])
  suppressWarnings(storage.mode(m) <- "double")
  rownames(m) <- as.character(tab[[1L]])[-1L]; colnames(m) <- cols
  want <- spec$selector
  stopifnot(all(want %in% colnames(m)))
  m <- m[, want, drop = FALSE]
  cmap <- transcript_map
  keep <- cmap$Ensembl_rna_identifier %in% rownames(m)
  cmap <- cmap[keep, ]
  sub <- m[cmap$Ensembl_rna_identifier, , drop = FALSE]
  out <- rowsum(sub, group = cmap$Ensembl_gene_identifier, reorder = TRUE)
  mapped_frac <- length(unique(cmap$Ensembl_rna_identifier)) / nrow(m)
  lost_frac <- 1 - sum(sub) / sum(m)
  list(counts = out, tpm_native = NULL,
       id_rule = "transcript counts summed to gene level (Ensembl GRCh37 annotation, NCBI gene2ensembl fallback)",
       id_note = sprintf("%d/%d source transcripts mapped (%.1f%%), reaching %d Ensembl genes; %.2f%% of total counts fall on unmapped transcripts",
                         length(unique(cmap$Ensembl_rna_identifier)), nrow(m), 100 * mapped_frac, nrow(out), 100 * lost_frac),
       integer_native = all(m == floor(m)))
}

h_gtex_v8 <- function(spec) {
  att <- read_table_gz(file.path(root, "raw", "gtex_v8",
                                 "GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"), header = TRUE)
  sel_all <- as.character(att$SAMPID[att$SMTSD == spec$selector])
  stopifnot(length(sel_all) > 20L)
  gct <- file.path(root, "raw", "gtex_v8", "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz")
  con <- gzfile(gct); hdr_raw <- readLines(con, n = 3L); close(con)
  all_samples <- strsplit(hdr_raw[[3L]], "\t")[[1L]][-(1:2)]
  # the sample-annotation table lists some additional sub-dissections that carry no
  # expression row in this release; only annotated samples with expression are used, and
  # the number dropped for that reason is recorded
  sel <- intersect(sel_all, all_samples)
  n_absent <- length(setdiff(sel_all, all_samples))
  stopifnot(length(sel) > 20L)
  idx <- match(sel, all_samples)
  stopifnot(!any(is.na(idx)))
  cols <- idx + 2L
  prog <- file.path(tmp_dir, paste0("gtex_", spec$ReferenceKey, ".awk"))
  # GCT layout: line 1 "#1.2", line 2 dimensions, line 3 header, line 4+ data.
  writeLines(sprintf(paste0("BEGIN { n = split(\"%s\", c, \",\") }",
                            " NR > 3 { printf \"%%s\", $1;",
                            " for (i = 1; i <= n; i++) printf \"\\t%%s\", $c[i]; printf \"\\n\" }"),
                     paste(cols, collapse = ",")), prog)
  tsv <- file.path(tmp_dir, paste0("gtex_", spec$ReferenceKey, ".tsv"))
  if (!file.exists(tsv)) {
    st <- system(sprintf("gzip -dc %s | awk -f %s > %s",
                         shQuote(gct), shQuote(prog), shQuote(tsv)))
    stopifnot(st == 0L)
  }
  tab <- read.delim(tsv, sep = "\t", header = FALSE, check.names = FALSE, quote = "")
  ids <- strip_ensembl_version(as.character(tab[[1L]]))
  m <- as.matrix(tab[, -1L, drop = FALSE])
  suppressWarnings(storage.mode(m) <- "double")
  rownames(m) <- ids; colnames(m) <- sel
  d <- dedupe_ensembl_rows(m, "counts")
  list(counts = d$matrix, tpm_native = NULL, id_rule = d$rule,
       id_note = sprintf("GTEx v8 RNASeQC gene reads, %s; %d samples with expression of %d annotated (%d annotated sub-dissections carry no expression row in this release)",
                         spec$selector, length(sel), length(sel_all), n_absent),
       integer_native = TRUE)
}

h_gdc_star <- function(spec) {
  man <- gdc_manifest[gdc_manifest$CohortLabel == spec$selector, ]
  stopifnot(nrow(man) > 0L)
  dir_ <- file.path(root, "raw", "gdc_star_counts", spec$selector)
  # files were saved as "<FileUUID>_<FileName>"; fall back to suffix matching if that
  # convention ever changes
  paths <- file.path(dir_, paste0(man$FileUUID, "_", man$FileName))
  if (!all(file.exists(paths))) {
    present <- list.files(dir_)
    paths <- vapply(man$FileName, function(fn) {
      hit <- present[endsWith(present, fn)]
      if (length(hit)) file.path(dir_, hit[[1L]]) else NA_character_
    }, character(1L))
  }
  stopifnot(!any(is.na(paths)), all(file.exists(paths)))
  one <- function(p) {
    d <- read.delim(p, sep = "\t", skip = 1L, header = TRUE, check.names = FALSE, quote = "",
                    colClasses = c("character", "NULL", "NULL", "integer",
                                   "NULL", "NULL", "NULL", "NULL", "NULL"))
    # columns 2, 3 and 5-9 of the file are dropped by colClasses, so the surviving
    # unstranded counts column is the second one of the frame, not the fourth
    keep <- grepl("^ENSG", d[[1L]])
    setNames(d[["unstranded"]][keep], d[[1L]][keep])
  }
  res <- mclapply(paths, one, mc.cores = min(16L, detectCores()))
  # intersect on gene identifiers, not on the count values they are attached to
  genes <- Reduce(intersect, lapply(res, names))
  m <- sapply(res, function(v) v[genes])
  rownames(m) <- genes
  colnames(m) <- man$SampleSubmitterID
  d <- dedupe_ensembl_rows(m, "counts")
  list(counts = d$matrix, tpm_native = NULL, id_rule = d$rule,
       id_note = sprintf("GDC augmented STAR gene counts (GENCODE v36), %s", spec$selector),
       integer_native = TRUE, samples = man$SampleSubmitterID)
}

handlers <- list(
  ncbi_counts = h_ncbi_counts,
  single_ensembl_counts = h_single_ensembl_counts,
  single_symbol_counts = h_single_symbol_counts,
  single_ensembl_fpkm = h_single_ensembl_fpkm,
  xlsx_fpkm = h_xlsx_fpkm,
  per_sample_rpkm = h_per_sample_rpkm,
  rsem_per_sample = h_rsem_per_sample,
  annotated_star_counts = h_annotated_star_counts,
  recount3_runs = h_recount3_runs,
  transcript_counts = h_transcript_counts,
  gtex_v8 = h_gtex_v8,
  gdc_star = h_gdc_star
)

## ---- main loop -----------------------------------------------------------

contract <- read.csv(contract_path, stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(all(c("GroupIDs", "ReferenceKey", "GroupLabel", "SourceKind", "SourcePath",
                "SelectorType", "SelectorValue", "ExpectedN") %in% names(contract)))

manifest_rows <- list(); qc_rows <- list(); fail_rows <- list()
for (i in seq_len(nrow(contract))) {
  spec <- as.list(contract[i, ])
  group_id <- strsplit(spec$GroupIDs, ";", fixed = TRUE)[[1L]][1L]
  if (!is.null(only) && !group_id %in% only) next
  obj_path <- file.path(out_dir, "objects", paste0(spec$ReferenceKey, ".rds"))
  cat(sprintf("[%02d/%d] %s %s (%s)\n", i, nrow(contract), spec$GroupIDs, spec$ReferenceKey, spec$SourceKind))
  if (file.exists(obj_path)) { cat("    cached\n"); next }

  res <- tryCatch({
    spec$selector <- if (spec$SelectorType == "run_sum") strsplit(spec$SelectorValue, ";", fixed = TRUE)[[1L]]
                     else strsplit(spec$SelectorValue, ";", fixed = TRUE)[[1L]]
    spec$run_groups <- spec$RunGroups
    spec$ExpectedN <- as.integer(spec$ExpectedN)
    h <- handlers[[spec$SourceKind]]
    stopifnot(!is.null(h))
    h(spec)
  }, error = function(e) e)

  if (inherits(res, "error")) {
    cat("    FAILED: ", conditionMessage(res), "\n", sep = "")
    fail_rows[[group_id]] <- data.frame(GroupIDs = spec$GroupIDs, ReferenceKey = spec$ReferenceKey,
                                        SourceKind = spec$SourceKind, Error = conditionMessage(res))
    next
  }

  counts <- res$counts
  if (!is.null(counts)) {
    # TPM is computed from the values as supplied (RSEM expected_count is fractional);
    # the archived counts block is the rounded integer matrix.
    tpm <- tpm_from_counts(counts, gene_lengths)
    counts <- round(counts)
    storage.mode(counts) <- "double"
  } else if (!is.null(res$fpkm)) {
    tpm <- fpkm_to_tpm(res$fpkm, gene_lengths)
  } else if (!is.null(res$tpm_native)) {
    tpm <- res$tpm_native
  } else stop("no expression values produced")
  if (is.null(colnames(tpm)) || !length(colnames(tpm))) {
    stop("handler produced a matrix without sample identifiers")
  }

  sample_meta <- data.frame(
    SampleID = colnames(tpm),
    GroupIDs = spec$GroupIDs,
    ReferenceKey = spec$ReferenceKey,
    SourceKind = spec$SourceKind,
    SelectionNote = spec$SelectionNote,
    stringsAsFactors = FALSE
  )
  obj <- finish_group(spec$GroupIDs, spec$ReferenceKey, spec$GroupLabel,
                      list(SourceKind = spec$SourceKind, SourcePath = spec$SourcePath,
                           SelectorType = spec$SelectorType, SelectorValue = spec$SelectorValue,
                           NativeValueRoute = if (!is.null(res$fpkm)) "FPKM/RPKM -> TPM (shared Ensembl 111 length table)"
                                              else if (!is.null(res$counts)) "counts -> TPM (shared Ensembl 111 length table)"
                                              else "source-native TPM",
                           IDRule = res$id_rule, IDNote = res$id_note,
                           NativeValuesInteger = res$integer_native,
                           SourceFileMD5 = if (file.exists(file.path(root, spec$SourcePath))) unname(tools::md5sum(file.path(root, spec$SourcePath))) else NA_character_,
                           Note = if (is.null(res$counts_note)) "" else res$counts_note),
                      counts, tpm, sample_meta)
  saveRDS(obj, obj_path, compress = "xz")

  write_matrix_tsv(obj$log2_tpm, file.path(out_dir, "matrices", paste0(spec$ReferenceKey, "_log2tpm.tsv.gz")))
  if (!is.null(counts)) write_matrix_tsv(counts, file.path(out_dir, "matrices", paste0(spec$ReferenceKey, "_counts.tsv.gz")))

  man <- data.frame(GroupIDs = spec$GroupIDs, ReferenceKey = spec$ReferenceKey,
                    GroupLabel = spec$GroupLabel, SourceKind = spec$SourceKind,
                    SourcePath = spec$SourcePath, SelectorValue = spec$SelectorValue,
                    Samples = ncol(obj$tpm), Genes = nrow(obj$tpm),
                    GenesWithCounts = if (is.null(counts)) NA_integer_ else nrow(counts),
                    ValueRoute = obj$Provenance$NativeValueRoute,
                    IDRule = obj$Provenance$IDRule, IDNote = obj$Provenance$IDNote,
                    NativeValuesInteger = obj$Provenance$NativeValuesInteger,
                    SourceFileMD5 = obj$Provenance$SourceFileMD5,
                    stringsAsFactors = FALSE)
  manifest_rows[[spec$ReferenceKey]] <- man
  q <- obj$qc; q$ReferenceKey <- spec$ReferenceKey; q$GroupIDs <- spec$GroupIDs
  qc_rows[[spec$ReferenceKey]] <- q
  cat(sprintf("    OK  %d samples x %d genes\n", ncol(obj$tpm), nrow(obj$tpm)))
}

if (length(manifest_rows)) {
  write.csv(do.call(rbind, manifest_rows), file.path(out_dir, "group_manifest.csv"), row.names = FALSE)
  write.csv(do.call(rbind, qc_rows), file.path(out_dir, "group_sample_qc.csv"), row.names = FALSE)
}
if (length(fail_rows)) {
  write.csv(do.call(rbind, fail_rows), file.path(out_dir, "group_failures.csv"), row.names = FALSE)
  cat("GROUPS_FAILED=", length(fail_rows), "\n", sep = "")
}
cat("EXTRACTION_DONE groups_written=", length(manifest_rows), " failures=", length(fail_rows), "\n", sep = "")
