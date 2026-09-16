#!/usr/bin/env Rscript
# Source-local reference RNA analysis under the 2026-09-16 material/perturbation rule.
# Input IDs are retained as supplied (Entrez or Ensembl); no gene-symbol joins occur.
args <- commandArgs(TRUE)
stopifnot(length(args) == 3L)
matrix_dir <- normalizePath(args[[1L]], mustWork=TRUE)
scope_file <- normalizePath(args[[2L]], mustWork=TRUE)
out_dir <- args[[3L]]
if (dir.exists(out_dir) && length(list.files(out_dir, all.files=TRUE, no..=TRUE))) {
  stop('Output directory already exists and is nonempty: ', out_dir, call.=FALSE)
}
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(out_dir, 'source_summaries'), showWarnings=FALSE)
dir.create(file.path(out_dir, 'pca'), showWarnings=FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
set.seed(25)

scope <- fread(scope_file, logical01=TRUE)
required_scope <- c('ReferenceKey','GroupIDs','AnatomicalOrCellIdentityMatch',
  'SelectedSamplePerturbation','VehicleExposure','NoDrugOrGeneticPerturbation',
  'MatrixN','StableID','IncludeSourceLocalAnalysis','Rationale')
stopifnot(identical(names(scope), required_scope), !anyDuplicated(scope$ReferenceKey))
scope <- scope[IncludeSourceLocalAnalysis == TRUE]
stopifnot(nrow(scope) > 0L, all(scope$AnatomicalOrCellIdentityMatch),
  all(scope$NoDrugOrGeneticPerturbation), all(scope$MatrixN >= 3L))

inventory <- fread(file.path(matrix_dir, 'combined_reference_inventory.csv'), logical01=TRUE)
stopifnot(all(scope$ReferenceKey %in% inventory$ReferenceKey))

safe_name <- function(x) gsub('[^A-Za-z0-9_.-]', '_', x)
qc_rows <- list()
pca_rows <- list()
manifest_rows <- list()
expression_rds <- list()

for (i in seq_len(nrow(scope))) {
  s <- scope[i]
  key <- s$ReferenceKey
  input <- file.path(matrix_dir, paste0(key, '.rds'))
  object <- readRDS(input)
  counts <- object$counts
  id_type <- object$IDType
  stopifnot(is.matrix(counts), ncol(counts) == s$MatrixN, !anyDuplicated(rownames(counts)),
    !anyDuplicated(colnames(counts)), all(is.finite(counts)), all(counts >= 0),
    all(counts == floor(counts)), all(colSums(counts) > 0))
  if (s$StableID == 'Entrez') stopifnot(identical(id_type, 'NCBI GeneID (Entrez)'), all(grepl('^[0-9]+$', rownames(counts))))
  if (s$StableID == 'Ensembl') stopifnot(identical(id_type, 'Ensembl'), all(grepl('^ENSG[0-9]+(\\.[0-9]+)?$', rownames(counts))))

  library_size <- colSums(counts)
  cpm <- sweep(counts, 2L, library_size / 1e6, '/')
  log2_cpm <- log2(cpm + 0.5)
  min_samples <- ceiling(ncol(counts) / 2)
  expressed <- rowSums(cpm >= 1) >= min_samples
  filtered <- log2_cpm[expressed, , drop=FALSE]
  stopifnot(nrow(filtered) > 1L)
  pairwise <- suppressWarnings(cor(filtered, method='spearman'))
  correlations <- pairwise[upper.tri(pairwise)]

  sample_qc <- data.table(
    ReferenceKey=key, GroupIDs=s$GroupIDs, StableID=s$StableID, SampleID=colnames(counts),
    LibrarySize=as.numeric(library_size), GenesInput=nrow(counts), GenesDetected=colSums(counts > 0),
    ZeroFraction=colMeans(counts == 0), GenesCPM1=colSums(cpm >= 1),
    SelectedSamplePerturbation=s$SelectedSamplePerturbation, VehicleExposure=s$VehicleExposure
  )
  qc_rows[[key]] <- sample_qc

  summary <- data.table(
    StableGeneID=rownames(counts), SamplesWithCPM1=rowSums(cpm >= 1),
    MeanLog2CPM=rowMeans(log2_cpm), SDLog2CPM=apply(log2_cpm, 1L, sd),
    ExpressedBySourceFilter=expressed
  )
  saveRDS(list(
    ReferenceKey=key, GroupIDs=s$GroupIDs, StableID=s$StableID, Condition=object$Condition,
    SelectedSamplePerturbation=s$SelectedSamplePerturbation, VehicleExposure=s$VehicleExposure,
    Normalization='log2(CPM + 0.5), library-size CPM; source-local only',
    ExpressionFilter=sprintf('CPM >= 1 in >= %d of %d samples', min_samples, ncol(counts)),
    sample_qc=sample_qc, gene_summary=summary, log2_cpm=log2_cpm
  ), file.path(out_dir, 'source_summaries', paste0(safe_name(key), '_source_local.rds')), compress=TRUE)
  expression_rds[[key]] <- file.path('source_summaries', paste0(safe_name(key), '_source_local.rds'))

  variance <- apply(filtered, 1L, var)
  top_n <- min(500L, length(variance))
  top <- names(sort(variance, decreasing=TRUE))[seq_len(top_n)]
  pca <- prcomp(t(filtered[top, , drop=FALSE]), center=TRUE, scale.=FALSE, rank.=2L)
  scores <- data.table(ReferenceKey=key, GroupIDs=s$GroupIDs, SampleID=rownames(pca$x),
    PC1=pca$x[, 1L], PC2=if (ncol(pca$x) >= 2L) pca$x[, 2L] else NA_real_,
    VehicleExposure=s$VehicleExposure)
  pca_rows[[key]] <- scores
  variance_percent <- (pca$sdev^2 / sum(pca$sdev^2)) * 100
  p <- ggplot(scores, aes(PC1, PC2, label=SampleID)) +
    geom_point(size=3, colour='#2F6690') +
    geom_text(vjust=-0.7, size=3) +
    scale_x_continuous(expand=expansion(mult=c(0.12, 0.25))) +
    scale_y_continuous(expand=expansion(mult=c(0.12, 0.18))) +
    theme_classic(base_size=11) +
    labs(title=paste0(key, ' source-local PCA'),
      subtitle=paste0(ncol(counts), ' samples; ', s$SelectedSamplePerturbation),
      x=sprintf('PC1 (%.1f%%)', variance_percent[[1L]]),
      y=sprintf('PC2 (%.1f%%)', variance_percent[[2L]]))
  ggsave(file.path(out_dir, 'pca', paste0(safe_name(key), '_PCA.png')), p, width=5.2, height=4.2, dpi=200)

  manifest_rows[[key]] <- data.table(
    ReferenceKey=key, GroupIDs=s$GroupIDs, InputFile=normalizePath(input), InputMD5=unname(tools::md5sum(input)),
    StableID=s$StableID, Samples=ncol(counts), GenesInput=nrow(counts), GenesAfterCPMFilter=sum(expressed),
    PairwiseSpearmanMin=min(correlations), PairwiseSpearmanMax=max(correlations),
    VehicleExposure=s$VehicleExposure, SelectedSamplePerturbation=s$SelectedSamplePerturbation,
    SourceLocalResult=expression_rds[[key]]
  )
}

fwrite(scope, file.path(out_dir, 'analysis_scope.csv'))
fwrite(rbindlist(manifest_rows), file.path(out_dir, 'analysis_manifest.csv'))
fwrite(rbindlist(qc_rows), file.path(out_dir, 'sample_qc.csv'))
fwrite(rbindlist(pca_rows), file.path(out_dir, 'pca_scores.csv'))
writeLines(trimws(capture.output(sessionInfo()), which='right'), file.path(out_dir, 'sessionInfo.txt'))
writeLines(c(
  '# Source-local RNA analysis (2026-09-16)', '',
  'Scope: material/cell identity matched, no selected-sample genetic perturbation or experimental drug.',
  'Vehicle-only DMSO controls are included under the user-defined rule and explicitly labelled; they are not merged with untreated samples.',
  'All expression values are log2(CPM + 0.5) calculated within their source. No raw counts, TPM or FPKM from distinct sources are merged.',
  'Stable IDs are retained as Entrez or Ensembl. Gene symbols were not used for joins or filtering.',
  'This is source-local expression/QC preparation. It is not a differential-expression comparison and does not change the historical AnalysisReady flag.'
), file.path(out_dir, 'README.md'))
cat('SOURCE_LOCAL_RNA_ANALYSIS_PASS: ', nrow(scope), ' sources; ', sum(lengths(strsplit(scope$GroupIDs, ';', fixed=TRUE))), ' group rows.\n', sep='')
