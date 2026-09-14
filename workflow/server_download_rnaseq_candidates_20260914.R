#!/usr/bin/env Rscript
# Download candidate RNA source material to an isolated server directory.
# This script intentionally does not alter the 31-group registry or mark any group analysis-ready.
set.seed(25)

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 2L || length(args) > 3L) {
  stop('Usage: Rscript --vanilla server_download_rnaseq_candidates_20260914.R <server_root> <specification_tsv> [TaskID[,TaskID...]]')
}
root <- normalizePath(args[[1]], mustWork=FALSE)
spec_path <- normalizePath(args[[2]], mustWork=TRUE)
if (!dir.exists(root)) dir.create(root, recursive=TRUE, showWarnings=FALSE)
for (d in c('raw', 'metadata', 'logs')) dir.create(file.path(root, d), recursive=TRUE, showWarnings=FALSE)

spec <- read.delim(spec_path, check.names=FALSE, stringsAsFactors=FALSE, na.strings='')
need <- c('TaskID','Source','Accession','GroupIDs','SelectedSamples','URL','RelativeTarget','Action','SelectionContract')
stopifnot(identical(names(spec), need), nrow(spec) > 0L, !anyDuplicated(spec$TaskID))
file.copy(spec_path, file.path(root, 'metadata', 'download_specification.tsv'), overwrite=TRUE)
all_spec <- spec
if (length(args) == 3L) {
  selected_tasks <- strsplit(args[[3]], ',', fixed=TRUE)[[1]]
  spec <- spec[spec$TaskID %in% selected_tasks, , drop=FALSE]
  stopifnot(nrow(spec) > 0L)
}

status_path <- file.path(root, 'metadata', 'download_status.tsv')
status_head <- c('TimestampUTC','TaskID','Source','Accession','GroupIDs','SelectedSamples','SourceURL','RelativeTarget','State','Bytes','SHA256','Detail')
append_status <- function(task, url, relative, state, file='', detail='') {
  bytes <- if (nzchar(file) && file.exists(file)) as.character(file.info(file)$size) else ''
  sha <- ''
  if (nzchar(file) && file.exists(file) && identical(state, 'complete')) {
    sha <- sub(' .*$', '', system2('sha256sum', file, stdout=TRUE, stderr=FALSE)[1])
  }
  row <- data.frame(
    TimestampUTC=format(Sys.time(), tz='UTC', usetz=TRUE), TaskID=task$TaskID, Source=task$Source,
    Accession=task$Accession, GroupIDs=task$GroupIDs, SelectedSamples=task$SelectedSamples,
    SourceURL=url, RelativeTarget=relative, State=state, Bytes=bytes, SHA256=sha, Detail=detail,
    check.names=FALSE, stringsAsFactors=FALSE
  )
  write.table(row, status_path, sep='\t', row.names=FALSE, quote=TRUE,
              col.names=!file.exists(status_path), append=file.exists(status_path), na='')
}

fetch <- function(task, url, relative) {
  dest <- file.path(root, 'raw', relative)
  dir.create(dirname(dest), recursive=TRUE, showWarnings=FALSE)
  valid_file <- function(path) {
    !grepl('\\.gz$', path, ignore.case=TRUE) || identical(suppressWarnings(system2('gzip', c('-t', path))), 0L)
  }
  if (file.exists(dest) && isTRUE(file.info(dest)$size > 0L)) {
    if (valid_file(dest)) {
      append_status(task, url, relative, 'complete', dest, 'already present; gzip integrity confirmed where applicable; SHA256 recalculated')
      return(invisible(TRUE))
    }
    unlink(dest)
  }
  part <- paste0(dest, '.part')
  code <- suppressWarnings(system2('wget', c('--continue', '--timeout=60', '--read-timeout=60', '--tries=5',
    '--waitretry=5', '--retry-connrefused', '--no-verbose', '--show-progress', '-O', part, url)))
  if (identical(code, 0L) && file.exists(part) && isTRUE(file.info(part)$size > 0L) && valid_file(part)) {
    if (!file.rename(part, dest)) stop('Could not finalize: ', dest)
    append_status(task, url, relative, 'complete', dest, 'downloaded with wget continuation support')
    return(invisible(TRUE))
  }
  append_status(task, url, relative, 'failed', if (file.exists(part)) part else '', paste('wget exit', code, '; retained .part for resumption'))
  invisible(FALSE)
}

geo_bucket <- function(gse) paste0('GSE', floor(as.integer(sub('^GSE', '', gse)) / 1000), 'nnn')
geo_url <- function(gse, section, name='') {
  paste0('https://ftp.ncbi.nlm.nih.gov/geo/series/', geo_bucket(gse), '/', gse, '/', section, '/', name)
}
geo_supplements <- function(gse) {
  html <- suppressWarnings(system2('wget', c('-qO-', geo_url(gse, 'suppl')), stdout=TRUE, stderr=FALSE))
  hits <- unlist(regmatches(html, gregexpr('href="[^"]+"', html, perl=TRUE)), use.names=FALSE)
  names <- gsub('^href="|"$', '', hits)
  names <- names[!names %in% c('../', './') & !grepl('/$', names)]
  names <- unique(basename(names))
  ext <- grepl('\\.(txt|tsv|csv|gct|xlsx|xls|rds)(\\.gz)?$', names, ignore.case=TRUE)
  archive_name <- grepl('\\.zip$', names, ignore.case=TRUE) & grepl('(matrix|count|expression|rna|fpkm|tpm)', names, ignore.case=TRUE)
  non_expression <- grepl('(methylated.*site|m6a.*site|\\bpeak\\b)', names, ignore.case=TRUE)
  excluded_arm <- identical(gse, 'GSE235595') & grepl('KI169(Ctrl)?_counts', names, ignore.case=TRUE)
  sort(names[(ext | archive_name) & !non_expression & !excluded_arm])
}
geo_matrix_files <- function(gse) {
  html <- suppressWarnings(system2('wget', c('-qO-', geo_url(gse, 'matrix')), stdout=TRUE, stderr=FALSE))
  hits <- unlist(regmatches(html, gregexpr('href="[^"]+_series_matrix\\.txt\\.gz"', html, perl=TRUE)), use.names=FALSE)
  names <- gsub('^href="|"$', '', hits)
  names <- unique(basename(names[nzchar(names)]))
  names <- sort(names[grepl('_series_matrix\\.txt\\.gz$', names, ignore.case=TRUE)])
  if (!length(names)) names <- paste0(gse, '_series_matrix.txt.gz')
  names
}
run_geo <- function(task) {
  gse <- task$Accession
  matrix_names <- geo_matrix_files(gse)
  fixed <- c(setNames(geo_url(gse, 'soft', paste0(gse, '_family.soft.gz')), paste0('metadata/', gse, '_family.soft.gz')),
             setNames(geo_url(gse, 'matrix', matrix_names), paste0('metadata/', matrix_names)))
  for (relative in names(fixed)) fetch(task, fixed[[relative]], file.path(task$RelativeTarget, relative))
  supp <- geo_supplements(gse)
  if (!length(supp)) {
    append_status(task, geo_url(gse, 'suppl'), file.path(task$RelativeTarget, 'supplement_listing'), 'listed_no_whitelisted_matrix', '',
                  'No txt/tsv/csv/gct/xlsx/xls/rds matrix-like GEO supplement was listed; raw-read resolution remains a separate queue.')
  } else {
    for (name in supp) fetch(task, geo_url(gse, 'suppl', utils::URLencode(name, reserved=TRUE)), file.path(task$RelativeTarget, 'supplement', name))
  }
}

for (i in seq_len(nrow(spec))) {
  task <- spec[i, , drop=FALSE]
  message('[', format(Sys.time(), tz='UTC', usetz=TRUE), '] ', task$TaskID)
  if (identical(task$Action, 'geo_bundle')) run_geo(task) else if (identical(task$Action, 'direct')) fetch(task, task$URL, task$RelativeTarget) else stop('Unknown action: ', task$Action)
}

raw_queue <- all_spec[all_spec$Source == 'GEO', c('TaskID','Accession','GroupIDs','SelectedSamples','SelectionContract')]
raw_queue$RawReadState <- 'not_downloaded: resolve selected GSM to exact runs after processed-material QC; use Ensembl/Entrez for later analysis'
write.table(raw_queue, file.path(root, 'metadata', 'raw_read_resolution_queue.tsv'), sep='\t', row.names=FALSE, quote=TRUE)
writeLines(c(
  'This server directory contains source downloads only; it is not an analysis result.',
  'All candidates remain AnalysisReady=FALSE until sample identity, library protocol, matrix/read integrity, stable-ID mapping and units are checked.',
  'No gene symbol is used as an analytical key. The fixed seed for downstream stochastic procedures is 25.'
), file.path(root, 'README.txt'))
message('DOWNLOAD_STAGE_COMPLETE: ', nrow(spec), ' tasks inspected. See metadata/download_status.tsv.')
