#!/usr/bin/env Rscript
# Build a reproducible GDC UUID/MD5 manifest for the four registered cancer-tissue candidates.
# This only fixes source files and never changes AnalysisReady in the RNA registry.
suppressPackageStartupMessages(library(jsonlite))
set.seed(25)

args <- commandArgs(trailingOnly=TRUE)
out_dir <- if (length(args)) args[[1]] else 'config'
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
endpoint <- 'https://api.gdc.cancer.gov/files'

`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x
value <- function(x, field) {
  y <- x[[field]] %||% ''
  if (length(y)) as.character(y[[1]]) else ''
}
fetch_hits <- function(project, sample_type) {
  filters <- list(op='and', content=list(
    list(op='in', content=list(field='cases.project.project_id', value=list(project))),
    list(op='in', content=list(field='data_type', value=list('Gene Expression Quantification'))),
    list(op='in', content=list(field='analysis.workflow_type', value=list('STAR - Counts'))),
    list(op='in', content=list(field='cases.samples.sample_type', value=list(sample_type)))
  ))
  request <- list(
    filters=filters, format='JSON', size=5000,
    fields='file_id,file_name,md5sum,file_size,data_type,analysis.workflow_type,cases.case_id,cases.project.project_id,cases.samples.sample_type,cases.samples.submitter_id,cases.diagnoses.primary_diagnosis'
  )
  json_path <- tempfile(fileext='.json')
  writeLines(toJSON(request, auto_unbox=TRUE), json_path)
  on.exit(unlink(json_path), add=TRUE)
  raw <- system2('curl', c('--fail','--silent','--show-error','--request','POST', '--header',shQuote('Content-Type: application/json'), '--data-binary',paste0('@',json_path),endpoint), stdout=TRUE, stderr=TRUE)
  answer <- fromJSON(paste(raw, collapse='\n'), simplifyVector=FALSE)
  hits <- answer$data$hits %||% list()
  total <- as.integer(answer$data$pagination$total %||% 0L)
  if (!identical(length(hits), total)) stop('GDC pagination incomplete for ', project, ' ', sample_type, ': ', length(hits), ' of ', total)
  hits
}
flatten_hit <- function(hit, label, group_ids, definition) {
  case <- (hit$cases %||% list(list()))[[1]]
  samples <- case$samples %||% list()
  diagnoses <- case$diagnoses %||% list()
  sample_type <- paste(unique(vapply(samples, value, '', field='sample_type')), collapse=';')
  submitter <- paste(unique(vapply(samples, value, '', field='submitter_id')), collapse=';')
  diagnosis <- paste(unique(vapply(diagnoses, value, '', field='primary_diagnosis')), collapse=';')
  data.frame(
    CohortLabel=label, GroupIDs=group_ids, SelectionDefinition=definition,
    FileUUID=value(hit,'file_id'), FileName=value(hit,'file_name'), MD5=value(hit,'md5sum'), Bytes=as.character(value(hit,'file_size')),
    ProjectID=value((case$project %||% list()),'project_id'), CaseUUID=value(case,'case_id'),
    SampleType=sample_type, SampleSubmitterID=submitter, PrimaryDiagnosis=diagnosis,
    DownloadURL=paste0('https://api.gdc.cancer.gov/data/', value(hit,'file_id')),
    stringsAsFactors=FALSE, check.names=FALSE
  )
}

queries <- list(
  list(label='TCGA_LIHC_adjacent_liver', groups='KLA31_09', project='TCGA-LIHC', sample_type='Solid Tissue Normal', definition='TCGA-LIHC Solid Tissue Normal: cancer-adjacent material; do not replace with healthy liver.'),
  list(label='TCGA_ESCA_ESCC_primary', groups='KLA31_10', project='TCGA-ESCA', sample_type='Primary Tumor', definition='TCGA-ESCA Primary Tumor with primary diagnosis containing squamous; excludes non-squamous ESCA.'),
  list(label='TCGA_PRAD_primary', groups='KLA31_11', project='TCGA-PRAD', sample_type='Primary Tumor', definition='TCGA-PRAD Primary Tumor.'),
  list(label='TCGA_LIHC_HCC_primary', groups='KLA31_12', project='TCGA-LIHC', sample_type='Primary Tumor', definition='TCGA-LIHC Primary Tumor.')
)
rows <- list()
for (q in queries) {
  x <- fetch_hits(q$project, q$sample_type)
  y <- do.call(rbind, lapply(x, flatten_hit, label=q$label, group_ids=q$groups, definition=q$definition))
  if (identical(q$label, 'TCGA_ESCA_ESCC_primary')) {
    y <- y[grepl('squamous', y$PrimaryDiagnosis, ignore.case=TRUE), , drop=FALSE]
  }
  if (!nrow(y)) stop('No eligible GDC files for ', q$label)
  if (anyDuplicated(y$FileUUID) || any(!nzchar(y$MD5)) || any(!nzchar(y$CaseUUID))) stop('Incomplete UUID, MD5, or case metadata for ', q$label)
  rows[[length(rows)+1L]] <- y
}
manifest <- do.call(rbind, rows)
stopifnot(!anyDuplicated(manifest$FileUUID), all(grepl('^[0-9a-f-]{36}$', manifest$FileUUID)), all(grepl('^[0-9a-f]{32}$', manifest$MD5)))
write.table(manifest, file.path(out_dir, 'gdc_star_counts_manifest_20260914.tsv'), sep='\t', row.names=FALSE, quote=TRUE)
summary <- aggregate(FileUUID ~ CohortLabel + GroupIDs + ProjectID + SampleType, data=manifest, FUN=length)
names(summary)[ncol(summary)] <- 'FileCount'
write.table(summary, file.path(out_dir, 'gdc_star_counts_manifest_20260914_summary.tsv'), sep='\t', row.names=FALSE, quote=TRUE)
cat('GDC_MANIFEST_PASS:', nrow(manifest), 'files across', nrow(summary), 'cohorts\n')
print(summary, row.names=FALSE)
