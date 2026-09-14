#!/usr/bin/env Rscript
# Retrieve the fixed DepMap release rows with the MD5 values published by DepMap.
# The files are external candidate material; downloading them does not make a group analysis-ready.
set.seed(25)
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 2L) stop('Usage: Rscript --vanilla server_download_depmap_20260914.R <server_root> <depmap_manifest_tsv>')
root <- normalizePath(args[[1]], mustWork=TRUE)
manifest_path <- normalizePath(args[[2]], mustWork=TRUE)
manifest <- read.delim(manifest_path, check.names=FALSE, stringsAsFactors=FALSE, quote='')
need <- c('TaskID','Release','ReleaseDate','FileName','URL','MD5','RelativeTarget','GroupIDs','Purpose','SourceCatalogURL')
stopifnot(identical(names(manifest), need), nrow(manifest) > 0L,
          !anyDuplicated(manifest$TaskID), all(grepl('^[0-9a-f]{32}$', manifest$MD5)),
          all(nzchar(manifest$URL)), all(nzchar(manifest$RelativeTarget)))
dir.create(file.path(root, 'raw'), recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(root, 'metadata'), recursive=TRUE, showWarnings=FALSE)
file.copy(manifest_path, file.path(root, 'metadata', 'depmap_expression_20260914.tsv'), overwrite=TRUE)

md5 <- function(path) tolower(sub(' .*$', '', system2('md5sum', path, stdout=TRUE, stderr=FALSE)[1]))
sha256 <- function(path) tolower(sub(' .*$', '', system2('sha256sum', path, stdout=TRUE, stderr=FALSE)[1]))
status_path <- file.path(root, 'metadata', 'depmap_download_status.tsv')
append_status <- function(row, state, file='', detail='') {
  observed <- if (nzchar(file) && file.exists(file)) md5(file) else ''
  item <- data.frame(TimestampUTC=format(Sys.time(), tz='UTC', usetz=TRUE), TaskID=row$TaskID,
    Release=row$Release, FileName=row$FileName, RelativeTarget=row$RelativeTarget, State=state,
    Bytes=if (nzchar(file) && file.exists(file)) as.character(file.info(file)$size) else '',
    ExpectedMD5=row$MD5, ObservedMD5=observed, Detail=detail,
    stringsAsFactors=FALSE, check.names=FALSE)
  write.table(item, status_path, sep='\t', quote=TRUE, row.names=FALSE,
    col.names=!file.exists(status_path), append=file.exists(status_path), na='')
}

catalog_url <- manifest$SourceCatalogURL[[1]]
catalog_dest <- file.path(root, 'metadata', 'depmap_no_captcha_download_catalog_20260914.csv')
catalog_part <- paste0(catalog_dest, '.part')
if (!file.exists(catalog_dest)) {
  code <- suppressWarnings(system2('wget', c('--timeout=60','--read-timeout=60','--tries=5',
    '--waitretry=5','--retry-connrefused','--no-verbose','-O',catalog_part,catalog_url)))
  if (identical(code, 0L) && file.exists(catalog_part) && file.info(catalog_part)$size > 100L) {
    file.rename(catalog_part, catalog_dest)
  } else {
    warning('Could not save the DepMap no-captcha release catalog; candidate files will still be attempted.')
  }
}

fetch <- function(row) {
  dest <- file.path(root, 'raw', row$RelativeTarget)
  dir.create(dirname(dest), recursive=TRUE, showWarnings=FALSE)
  if (file.exists(dest)) {
    if (identical(md5(dest), tolower(row$MD5))) {
      append_status(row, 'complete', dest, 'already present; official DepMap MD5 verified')
      return(invisible(TRUE))
    }
    unlink(dest)
  }
  part <- paste0(dest, '.part')
  code <- suppressWarnings(system2('wget', c('--continue','--timeout=60','--read-timeout=60','--tries=5',
    '--waitretry=5','--retry-connrefused','--no-verbose','--show-progress','-O',part,row$URL)))
  if (identical(code, 0L) && file.exists(part) && identical(md5(part), tolower(row$MD5))) {
    if (!file.rename(part, dest)) stop('Could not finalize ', dest)
    append_status(row, 'complete', dest, paste('downloaded; official DepMap MD5 verified; SHA256=', sha256(dest)))
    return(invisible(TRUE))
  }
  append_status(row, 'failed', if (file.exists(part)) part else '',
    paste('wget exit', code, '; MD5 mismatch or incomplete file retained as .part'))
  invisible(FALSE)
}

for (i in seq_len(nrow(manifest))) {
  message('[', format(Sys.time(), tz='UTC', usetz=TRUE), '] ', i, '/', nrow(manifest), ' ', manifest$TaskID[i])
  fetch(manifest[i, , drop=FALSE])
}
readme <- c(
  'DepMap candidate data only.',
  'Release fixed to DepMap Public 24Q4 (2024-12-16), because the current 26Q1 no-captcha catalog exposes metadata but no usable URL for these files.',
  'Every completed file has been compared with the source MD5 recorded in metadata/depmap_expression_20260914.tsv.',
  'Use ACH model IDs and stable gene identifiers for analysis; gene symbols are display-only.',
  'No row is analysis-ready until release, model identity, stable-ID columns, scale and observation-unit QC are completed.'
)
writeLines(readme, file.path(root, 'metadata', 'depmap_README.txt'))
message('DEPMAP_DOWNLOAD_STAGE_COMPLETE: ', nrow(manifest), ' manifest rows inspected.')
