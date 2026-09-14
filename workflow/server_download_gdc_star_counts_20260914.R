#!/usr/bin/env Rscript
# Retrieve only the UUIDs fixed in the GDC manifest, with source MD5 verification.
# The files are external candidate material and do not make a group analysis-ready.
set.seed(25)
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 2L) stop('Usage: Rscript --vanilla server_download_gdc_star_counts_20260914.R <server_root> <gdc_manifest_tsv>')
root <- normalizePath(args[[1]], mustWork=TRUE)
manifest_path <- normalizePath(args[[2]], mustWork=TRUE)
manifest <- read.delim(manifest_path, check.names=FALSE, stringsAsFactors=FALSE)
need <- c('CohortLabel','GroupIDs','SelectionDefinition','FileUUID','FileName','MD5','Bytes','ProjectID','CaseUUID','SampleType','SampleSubmitterID','PrimaryDiagnosis','DownloadURL')
stopifnot(identical(names(manifest), need), nrow(manifest) > 0L, !anyDuplicated(manifest$FileUUID),
          all(grepl('^[0-9a-f-]{36}$', manifest$FileUUID)), all(grepl('^[0-9a-f]{32}$', manifest$MD5)))
dir.create(file.path(root, 'raw', 'gdc_star_counts'), recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(root, 'metadata'), recursive=TRUE, showWarnings=FALSE)
file.copy(manifest_path, file.path(root, 'metadata', 'gdc_star_counts_manifest_20260914.tsv'), overwrite=TRUE)
status_path <- file.path(root, 'metadata', 'gdc_star_counts_download_status.tsv')

md5 <- function(path) tolower(sub(' .*$', '', system2('md5sum', path, stdout=TRUE, stderr=FALSE)[1]))
append_status <- function(row, relative, state, file='', detail='') {
  observed <- if (nzchar(file) && file.exists(file)) md5(file) else ''
  item <- data.frame(TimestampUTC=format(Sys.time(), tz='UTC', usetz=TRUE), CohortLabel=row$CohortLabel,
    GroupIDs=row$GroupIDs, FileUUID=row$FileUUID, RelativeTarget=relative, State=state,
    Bytes=if (nzchar(file) && file.exists(file)) as.character(file.info(file)$size) else '',
    ExpectedMD5=row$MD5, ObservedMD5=observed, Detail=detail, stringsAsFactors=FALSE, check.names=FALSE)
  write.table(item, status_path, sep='\t', quote=TRUE, row.names=FALSE,
    col.names=!file.exists(status_path), append=file.exists(status_path), na='')
}
fetch <- function(row) {
  safe_name <- basename(row$FileName)
  relative <- file.path('gdc_star_counts', row$CohortLabel, paste0(row$FileUUID, '_', safe_name))
  dest <- file.path(root, 'raw', relative)
  dir.create(dirname(dest), recursive=TRUE, showWarnings=FALSE)
  if (file.exists(dest)) {
    if (identical(md5(dest), tolower(row$MD5))) {
      append_status(row, relative, 'complete', dest, 'already present; GDC MD5 verified')
      return(invisible(TRUE))
    }
    unlink(dest)
  }
  part <- paste0(dest, '.part')
  code <- suppressWarnings(system2('wget', c('--continue','--timeout=60','--read-timeout=60','--tries=5','--waitretry=5','--retry-connrefused','--no-verbose','--show-progress','-O',part,row$DownloadURL)))
  if (identical(code, 0L) && file.exists(part) && identical(md5(part), tolower(row$MD5))) {
    if (!file.rename(part, dest)) stop('Could not finalize ', dest)
    append_status(row, relative, 'complete', dest, 'downloaded; official GDC MD5 verified')
    return(invisible(TRUE))
  }
  append_status(row, relative, 'failed', if (file.exists(part)) part else '', paste('wget exit', code, '; MD5 mismatch or incomplete file retained as .part'))
  invisible(FALSE)
}
for (i in seq_len(nrow(manifest))) {
  message('[', format(Sys.time(), tz='UTC', usetz=TRUE), '] ', i, '/', nrow(manifest), ' ', manifest$CohortLabel[i], ' ', manifest$FileUUID[i])
  fetch(manifest[i, , drop=FALSE])
}
writeLines(c('GDC STAR-count candidate data only.', 'Every completed file has been compared with the source MD5 recorded in metadata/gdc_star_counts_manifest_20260914.tsv.', 'No row is analysis-ready until case/sample, stable-ID and count-unit QC is completed.'), file.path(root, 'metadata', 'gdc_star_counts_README.txt'))
message('GDC_DOWNLOAD_STAGE_COMPLETE: ', nrow(manifest), ' manifest rows inspected.')
