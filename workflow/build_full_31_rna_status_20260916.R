#!/usr/bin/env Rscript
# A full 31-row ledger: distinguish analysis-ready sources from queued and blocked candidates.
args <- commandArgs(TRUE)
stopifnot(length(args) == 1L)
suppressPackageStartupMessages(library(data.table))
out_dir <- args[[1L]]
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
    'STARTED_SOURCE_LOCAL','QUEUED_SELECTION_AND_MATRIX_QC','STARTED_SOURCE_LOCAL','QUEUED_SELECTION_AND_MATRIX_QC','STARTED_SOURCE_LOCAL','QUEUED_SELECTION_AND_MATRIX_QC',
    'STARTED_SOURCE_LOCAL','QUEUED_SELECTION_AND_MATRIX_QC','QUEUED_SELECTION_AND_MATRIX_QC','STARTED_SOURCE_LOCAL','STARTED_SOURCE_LOCAL',
    'QUEUED_SELECTION_AND_MATRIX_QC','QUEUED_SELECTION_AND_MATRIX_QC','STARTED_SOURCE_LOCAL'
  ),
  CandidateSeries=c(
    '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '',
    'GSE163787', '', 'GSE171750', '', 'GSE283812', '', 'GSE266884',
    '', '', '', '', '', '', 'GSE119834', ''
  ),
  CandidateSamples=c(
    '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '',
    'GSM4987488; SRX9725233', '', 'GSM5233375;GSM5233376;GSM5233377', '', 'GSM8672552;GSM8672553;GSM8672554;GSM8672555', '', 'GSM8255530;GSM8255531;GSM8255534;GSM8255535',
    '', '', '', '', '', '', 'GSM3384847-GSM3384855', ''
  ),
  CandidateReferenceKey=c(
    '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '',
    'GSE163787_TALL104', '', 'GSE171750_A549_untreated', '', 'GSE283812_T47D_vehicle', '', 'GSE266884_MES28_vector_control',
    '', '', '', '', '', '', 'GSE119834_NSC', ''
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
    'Exact TALL-104 bulk source has one library only (GSM4987488); the three-replicate lead in GEO is microarray.',
    'Stable-ID count matrix validated.',
    'Three untreated A549 biological replicates confirmed in GSE171750; author RSEM table uses gene symbols.',
    'Stable-ID count matrix validated.',
    'Four T-47D vehicle control replicates confirmed in GSE283812; author count table uses gene symbols.',
    'Stable-ID count matrix validated.',
    'Four unperturbed MES28 vector control replicates (shMCT1-Ctrl + shKU70-Ctrl) confirmed; NCBI Entrez counts matrix is downloaded.',
    'Stable-ID count matrix validated.',
    'WT untreated HEK293T RSEM gene files are identified; source-local matrix QC has not run.',
    'Three HMC3 vehicle baseline candidates and a STAR count matrix are identified; source-local matrix QC has not run.',
    'Stable-ID count matrix validated.', 'Stable-ID count matrix validated.',
    'recount3 route and untreated MCF10A candidate libraries are identified; exact independent sample map/QC remains.',
    'Nine untreated, category-matched NSC libraries are selected by user decision. The downloaded author FPKM matrix is keyed only by gene symbol.',
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
    'Decision pending: include GSM4987488 as n=1 descriptive background or substitute matched T-ALL line with n>=3 replicates.',
    'None; source-local analysis started.',
    'Map candidate count matrix from Symbol to Ensembl/Entrez and run source-local QC.',
    'None; source-local analysis started.',
    'Map candidate count matrix from Symbol to Ensembl/Entrez and run source-local QC.',
    'None; source-local analysis started.',
    'Extract the four control samples from NCBI Entrez matrix and run source-local QC.',
    'None; source-local analysis started.',
    'Validate three selected WT untreated RSEM files.',
    'Validate three selected HMC3 vehicle files and exclude SLO arm.',
    'None; source-local analysis started.', 'None; source-local analysis started.',
    'Validate independent MCF10A samples and source-local matrix.',
    'Requantify the nine selected raw RNA-seq libraries against Ensembl; do not analyse the author gene-symbol FPKM table.',
    'None; source-local analysis started.'
  )
)

out <- merge(catalog, status, by='GroupID', all.x=TRUE, sort=FALSE)
out <- merge(out, gate, by='GroupID', all.x=TRUE, sort=FALSE)
out <- merge(out, confirmed, by='GroupID', all.x=TRUE, sort=FALSE)
setorder(out, RowOrder)
out[, CurrentReferenceKey := fifelse(!is.na(ConfirmedReferenceKey), ConfirmedReferenceKey,
  fifelse(nzchar(CandidateReferenceKey), CandidateReferenceKey, GateReferenceKey))]
out[, CurrentGSE := fifelse(!is.na(ConfirmedGSE), ConfirmedGSE,
  fifelse(nzchar(CandidateSeries), CandidateSeries,
    fifelse(nzchar(ProposedReplacementSource), ProposedReplacementSource, Source)))]
out[, CurrentSamples := fifelse(!is.na(ConfirmedGSM), ConfirmedGSM,
  fifelse(nzchar(CandidateSamples), CandidateSamples,
    fifelse(nzchar(ProposedReplacementSamples), ProposedReplacementSamples, SampleSelection)))]
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

char_cols <- names(out)[sapply(out, is.character)]
out[, (char_cols) := lapply(.SD, trimws), .SDcols=char_cols]

stopifnot(nrow(out) == 31L, !anyDuplicated(out$GroupID),
  out[ScopeStatus == 'STARTED_SOURCE_LOCAL', .N] == 13L,
  out[ScopeStatus == 'QUEUED_SELECTION_AND_MATRIX_QC', .N] == 17L,
  out[ScopeStatus == 'BLOCKED_CURRENTLY', .N] == 1L,
  all(out[InAnalysisNow == TRUE, MatrixQCComplete]),
  all(out[InAnalysisNow == TRUE, NumericQCPass]))
fwrite(out, file.path(out_dir, 'rna_31_group_status.csv'))
fwrite(out, file.path(out_dir, 'rna_31_group_status.tsv'), sep='\t')
writeLines(c(
  '# Full 31-group RNA status (2026-09-16)', '',
  'This is the master ledger for all 31 proteome/Kla group rows. It distinguishes the 13 rows already in source-local analysis from 17 rows with confirmed candidates queued for selection or matrix QC, and the 1 row with an actual remaining biological/replicate constraint (KLA31_17 TALL-104).',
  'Queued means a verified eligible candidate source exists with matching unperturbed biological replicates, pending stable-ID mapping, sample freezing, or local matrix QC. Blocked means the exact candidate fails a required condition (e.g. single bulk library without within-group biological replication).',
  'The 13 started rows use material/cell identity matching and selected samples with no knockdown, overexpression or experimental drug. DMSO vehicle controls remain explicitly marked and separated from untreated sources.'
), file.path(out_dir, 'README.md'))
cat('FULL_31_STATUS_PASS: started=13 queued=17 blocked=1\n')
