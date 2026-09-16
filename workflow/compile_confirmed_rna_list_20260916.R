#!/usr/bin/env Rscript
# Expand the user-criteria source-local scope to a reviewable KLA31 group list.
args <- commandArgs(TRUE)
stopifnot(length(args) == 3L)
suppressPackageStartupMessages(library(data.table))
scope <- fread(args[[1L]], logical01=TRUE)
analysis_dir <- normalizePath(args[[2L]], mustWork=TRUE)
out_dir <- args[[3L]]
if (dir.exists(out_dir) && length(list.files(out_dir, all.files=TRUE, no..=TRUE))) {
  stop('Output directory already exists and is nonempty: ', out_dir, call.=FALSE)
}
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

catalog <- fread('outputs/20260914_reference_material_audit/rnaseq_reference_candidate_31.csv')
samples <- rbindlist(list(
  fread('audit/20260915_replicate_resolution/resolved_matrices/selected_samples.csv'),
  fread('audit/20260915_replicate_resolution/resolved_matrices/native_selected_samples.csv')[, .(ReferenceKey, GSE, GSM, Condition)]
))
manifest <- fread(file.path(analysis_dir, 'analysis_manifest.csv'))
stopifnot(!anyDuplicated(scope$ReferenceKey), !anyDuplicated(manifest$ReferenceKey),
  all(scope$ReferenceKey %in% manifest$ReferenceKey), all(scope$ReferenceKey %in% samples$ReferenceKey))

expanded <- rbindlist(lapply(seq_len(nrow(scope)), function(i) {
  x <- scope[i]
  data.table(ReferenceKey=x$ReferenceKey, GroupID=unlist(strsplit(x$GroupIDs, ';', fixed=TRUE)),
    AnatomicalOrCellIdentityMatch=x$AnatomicalOrCellIdentityMatch,
    SelectedSamplePerturbation=x$SelectedSamplePerturbation,
    VehicleExposure=x$VehicleExposure,
    NoDrugOrGeneticPerturbation=x$NoDrugOrGeneticPerturbation,
    IncludeSourceLocalAnalysis=x$IncludeSourceLocalAnalysis,
    InclusionRationale=x$Rationale)
}))
setkey(catalog, GroupID); setkey(expanded, GroupID); setkey(manifest, ReferenceKey)
out <- catalog[expanded]
setnames(out, c('ReferenceKey', 'i.ReferenceKey'), c('CatalogReferenceKey', 'ReferenceKey'))
out <- manifest[out, on='ReferenceKey']
out[, GSE := vapply(ReferenceKey, function(key) unique(samples[ReferenceKey == key, GSE]), character(1))]
out[, SelectedGSM := vapply(ReferenceKey, function(key) paste(samples[ReferenceKey == key, GSM], collapse=';'), character(1))]
out[, SourceLocalAnalysisStatus := 'CONFIRMED: material/cell identity match; selected samples have no genetic perturbation or experimental drug']
out[VehicleExposure == TRUE, SourceLocalAnalysisStatus := paste0(SourceLocalAnalysisStatus, '; DMSO vehicle exposure recorded')]
setcolorder(out, c('GroupID','PXD','SampleGroup','CatalogReferenceKey','ReferenceKey','GSE','SelectedGSM','Samples',
  'AnatomicalOrCellIdentityMatch','SelectedSamplePerturbation','VehicleExposure',
  'NoDrugOrGeneticPerturbation','StableID','GenesAfterCPMFilter','PairwiseSpearmanMin',
  'PairwiseSpearmanMax','SourceLocalAnalysisStatus','InclusionRationale','SourceLocalResult'))
stopifnot(nrow(out) == 13L, !anyDuplicated(out$GroupID), all(out$AnatomicalOrCellIdentityMatch),
  all(out$NoDrugOrGeneticPerturbation), all(out$Samples >= 3L), all(out$IncludeSourceLocalAnalysis))
fwrite(out, file.path(out_dir, 'confirmed_source_local_rna_list.csv'))
fwrite(out[, .(GroupID, PXD, SampleGroup, ReferenceKey, GSE, SelectedGSM, Samples,
  SelectedSamplePerturbation, VehicleExposure, StableID, SourceLocalAnalysisStatus)],
  file.path(out_dir, 'confirmed_source_local_rna_list.tsv'), sep='\t')
writeLines(c(
  '# Confirmed source-local RNA list (2026-09-16)', '',
  'Inclusion rule confirmed with the user: tissue/cell identity must match; selected RNA samples cannot have knockdown, overexpression or experimental drug treatment.',
  'Vehicle-only DMSO samples are retained because they are vehicle controls, not the experimental drug arm. Their vehicle exposure is explicit in the list and they remain separate from untreated samples.',
  'All 13 group rows have a stable-ID gene-count matrix and at least three selected samples. Source-local QC and log2(CPM + 0.5) summaries have started in ../20260916_source_local_rna_analysis/.',
  'The list does not license raw-count merging across studies or claim matched-patient multi-omics.'
), file.path(out_dir, 'README.md'))
cat('CONFIRMED_RNA_LIST_PASS: ', nrow(out), ' group rows from ', uniqueN(out$ReferenceKey), ' sources.\n', sep='')
