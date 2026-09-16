#!/usr/bin/env Rscript
# A full 31-row ledger: distinguish analysis-ready sources from queued and blocked candidates.
args <- commandArgs(TRUE)
stopifnot(length(args) == 1L)
suppressPackageStartupMessages(library(data.table))
out_dir <- args[[1L]]
if (dir.exists(out_dir) && length(list.files(out_dir, all.files=TRUE, no..=TRUE))) {
  stop('Output directory already exists and is nonempty: ', out_dir, call.=FALSE)
}
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

catalog <- fread('outputs/20260914_reference_material_audit/rnaseq_reference_candidate_31.csv')
gate <- fread('config/rnaseq_replicate_resolution_gate_20260915.csv', logical01=TRUE)[, .(
  GroupID, GateReferenceKey=ReferenceKey, ProposedReplacementSource,
  ProposedReplacementSamples, MatrixSampleN, NumericQCPass, GateStableID=StableID
)]
confirmed <- fread('audit/20260916_confirmed_rna_list/confirmed_source_local_rna_list.csv', logical01=TRUE)[, .(
  GroupID, ConfirmedReferenceKey=ReferenceKey, ConfirmedGSE=GSE,
  ConfirmedGSM=SelectedGSM, ConfirmedSamples=Samples,
  ConfirmedPerturbation=SelectedSamplePerturbation,
  ConfirmedVehicleExposure=VehicleExposure, ConfirmedStableID=StableID
)]

status <- data.table(
  GroupID=sprintf('KLA31_%02d', 1:31),
  ScopeStatus=c(
    'STARTED_SOURCE_LOCAL',
    'QUEUED_SELECTION_AND_MATRIX_QC','QUEUED_SELECTION_AND_MATRIX_QC','QUEUED_SELECTION_AND_MATRIX_QC',
    'QUEUED_SELECTION_AND_MATRIX_QC','QUEUED_SELECTION_AND_MATRIX_QC','STARTED_SOURCE_LOCAL','QUEUED_SELECTION_AND_MATRIX_QC',
    'QUEUED_SELECTION_AND_MATRIX_QC','QUEUED_SELECTION_AND_MATRIX_QC','QUEUED_SELECTION_AND_MATRIX_QC','QUEUED_SELECTION_AND_MATRIX_QC',
    'STARTED_SOURCE_LOCAL','STARTED_SOURCE_LOCAL','STARTED_SOURCE_LOCAL','STARTED_SOURCE_LOCAL','BLOCKED_CURRENTLY',
    'STARTED_SOURCE_LOCAL','BLOCKED_CURRENTLY','STARTED_SOURCE_LOCAL','BLOCKED_CURRENTLY','STARTED_SOURCE_LOCAL','BLOCKED_CURRENTLY',
    'STARTED_SOURCE_LOCAL','QUEUED_SELECTION_AND_MATRIX_QC','QUEUED_SELECTION_AND_MATRIX_QC','STARTED_SOURCE_LOCAL','STARTED_SOURCE_LOCAL',
    'QUEUED_SELECTION_AND_MATRIX_QC','BLOCKED_CURRENTLY','STARTED_SOURCE_LOCAL'
  ),
  DataSituation=c(
    'Stable-ID count matrix validated.',
    'GTEx v8 counts/TPM and donor metadata are available; exact lung samples not selected.',
    'Three scar input candidates exist; pair identity and source-matrix stable-ID QC remain.',
    'Three adjacent-skin input candidates exist; pair identity and source-matrix stable-ID QC remain.',
    'GTEx v8 counts/TPM and donor metadata are available; exact hippocampus samples not selected.',
    'Published placenta count matrix exists; gestational-age and sampling-layer filter remains.',
    'Stable-ID count matrix validated.',
    'Eighteen BPH transition-zone candidates have source matrices; clinical selector remains.',
    'GDC STAR-count files are downloaded; exact Solid Tissue Normal donor/file set is not frozen.',
    'GDC STAR-count files are downloaded; ESCC histology and Primary Tumor selector is not frozen.',
    'GDC STAR-count files are downloaded; Primary Tumor sample selector is not frozen.',
    'GDC STAR-count files are downloaded; Primary Tumor sample selector is not frozen.',
    'Stable-ID count matrix validated.', 'Stable-ID count matrix validated.', 'Stable-ID count matrix validated.', 'Stable-ID count matrix validated.',
    'Exact TALL-104 bulk source has one library only; the three-replicate lead is microarray.',
    'Stable-ID count matrix validated.',
    'Three untreated A549 samples found, but the available author matrix uses gene symbols rather than stable IDs.',
    'Stable-ID count matrix validated.',
    'T47D control replicates found, but available complete count tables use gene symbols rather than stable IDs.',
    'Stable-ID count matrix validated.',
    'Exact MES28 shMCT1-control source has only two libraries; separate shKU70 controls cannot be pooled.',
    'Stable-ID count matrix validated.',
    'WT untreated HEK293T RSEM gene files are identified; source-local matrix QC has not run.',
    'Three HMC3 vehicle baseline candidates and a STAR count matrix are identified; source-local matrix QC has not run.',
    'Stable-ID count matrix validated.', 'Stable-ID count matrix validated.',
    'recount3 route and untreated MCF10A candidate libraries are identified; exact independent sample map/QC remains.',
    'Nine NSC libraries exist, but the exact neural-stem-cell model/developmental identity is not established.',
    'Stable-ID count matrix validated.'
  ),
  RequiredBeforeStart=c(
    'None; source-local analysis started.',
    'Freeze donor/sample IDs and run source-local QC.',
    'Confirm paired donor identity and validate gene matrix.',
    'Confirm paired donor identity and validate gene matrix.',
    'Freeze donor/sample IDs and run source-local QC.',
    'Lock gestational-age/sampling filters and run source-local QC.',
    'None; source-local analysis started.',
    'Lock clinical/transition-zone inclusion and run source-local QC.',
    'Freeze adjacent-liver UUIDs/donors and run source-local QC.',
    'Freeze ESCC histology/UUIDs and run source-local QC.',
    'Freeze primary-tumor UUIDs/donors and run source-local QC.',
    'Freeze primary-tumor UUIDs/donors and run source-local QC.',
    'None; source-local analysis started.', 'None; source-local analysis started.', 'None; source-local analysis started.', 'None; source-local analysis started.',
    'Find >=3 matched bulk RNA-seq samples or keep descriptive only.',
    'None; source-local analysis started.',
    'Obtain source annotation mapping or reprocess selected reads to stable IDs.',
    'None; source-local analysis started.',
    'Obtain source annotation mapping or reprocess selected reads to stable IDs.',
    'None; source-local analysis started.',
    'Find >=1 additional matched MES28 control library; do not pool designs.',
    'None; source-local analysis started.',
    'Validate three selected WT untreated RSEM files.',
    'Validate three selected HMC3 vehicle files and exclude SLO arm.',
    'None; source-local analysis started.', 'None; source-local analysis started.',
    'Validate independent MCF10A samples and source-local matrix.',
    'Find a model-identity-matched NSC reference; current panel remains approximate.',
    'None; source-local analysis started.'
  )
)

out <- merge(catalog, status, by='GroupID', all.x=TRUE, sort=FALSE)
out <- merge(out, gate, by='GroupID', all.x=TRUE, sort=FALSE)
out <- merge(out, confirmed, by='GroupID', all.x=TRUE, sort=FALSE)
setorder(out, RowOrder)
out[, CurrentReferenceKey := fifelse(!is.na(ConfirmedReferenceKey), ConfirmedReferenceKey, GateReferenceKey)]
out[, CurrentGSE := fifelse(!is.na(ConfirmedGSE), ConfirmedGSE,
  fifelse(nzchar(ProposedReplacementSource), ProposedReplacementSource, Source))]
out[, CurrentSamples := fifelse(!is.na(ConfirmedGSM), ConfirmedGSM,
  fifelse(nzchar(ProposedReplacementSamples), ProposedReplacementSamples, SampleSelection))]
out[, MatrixQCComplete := !is.na(ConfirmedSamples)]
out[, SelectedSamplePerturbation := ConfirmedPerturbation]
out[, VehicleExposure := ConfirmedVehicleExposure]
out[, StableID := fifelse(!is.na(ConfirmedStableID), ConfirmedStableID, GateStableID)]
out[, InAnalysisNow := ScopeStatus == 'STARTED_SOURCE_LOCAL']
setcolorder(out, c('GroupID','PXD','SampleGroup','Category','ScopeStatus','InAnalysisNow',
  'Source','CurrentGSE','CurrentReferenceKey','CurrentSamples','MatrixQCComplete',
  'DataSituation','RequiredBeforeStart','SelectedSamplePerturbation','VehicleExposure','StableID',
  'MatrixSampleN','NumericQCPass','EvidenceStatus','AssayAndUnits','Limitations'))
out <- out[, .SD, .SDcols=c('GroupID','PXD','SampleGroup','Category','ScopeStatus','InAnalysisNow',
  'Source','CurrentGSE','CurrentReferenceKey','CurrentSamples','MatrixQCComplete',
  'DataSituation','RequiredBeforeStart','SelectedSamplePerturbation','VehicleExposure','StableID',
  'MatrixSampleN','NumericQCPass','EvidenceStatus','AssayAndUnits','Limitations')]
stopifnot(nrow(out) == 31L, !anyDuplicated(out$GroupID),
  out[ScopeStatus == 'STARTED_SOURCE_LOCAL', .N] == 13L,
  out[ScopeStatus == 'QUEUED_SELECTION_AND_MATRIX_QC', .N] == 13L,
  out[ScopeStatus == 'BLOCKED_CURRENTLY', .N] == 5L,
  all(out[InAnalysisNow == TRUE, MatrixQCComplete]),
  all(out[InAnalysisNow == TRUE, NumericQCPass]))
fwrite(out, file.path(out_dir, 'rna_31_group_status.csv'))
fwrite(out, file.path(out_dir, 'rna_31_group_status.tsv'), sep='\t')
writeLines(c(
  '# Full 31-group RNA status (2026-09-16)', '',
  'This is the master ledger for all 31 proteome/Kla group rows. It distinguishes the 13 rows already in source-local analysis from rows with data that still need selection or matrix QC, and rows with an actual current eligibility gap.',
  'Queued means a relevant source or file exists, not that it has been accepted. Blocked means the current exact candidate fails a required condition: insufficient matched bulk replication, stable-ID absence, or material/model mismatch.',
  'The 13 started rows use material/cell identity matching and selected samples with no knockdown, overexpression or experimental drug. DMSO vehicle controls remain explicitly marked and separated from untreated sources.'
), file.path(out_dir, 'README.md'))
cat('FULL_31_STATUS_PASS: started=13 queued=13 blocked=5\n')
