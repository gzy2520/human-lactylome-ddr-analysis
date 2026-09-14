#!/usr/bin/env Rscript
# Read-only validation. No data regeneration, downloads or output writes.
suppressPackageStartupMessages({library(data.table); library(digest)})
set.seed(25)
arg <- grep('^--file=', commandArgs(FALSE), value=TRUE)
root <- normalizePath(file.path(dirname(sub('^--file=', '', arg[[1L]])), '..'))
setwd(root)
check <- function(ok, msg) if (!isTRUE(ok)) stop(msg, call.=FALSE)
verify_hashes <- function(table, base='.', path_col='Path') {
  paths <- file.path(base, table[[path_col]])
  check(all(file.exists(paths)), 'Missing files in hash manifest')
  hashes <- vapply(paths, function(p) digest(file=p, algo='sha256'), character(1))
  check(all(hashes==table$SHA256), paste('Hash mismatch:', base))
  nrow(table)
}
for (p in c('results','renv/library')) {
  check(dir.exists(p) && !nzchar(Sys.readlink(p)), paste('Expected independent directory:', p))
}
frozen <- 'corrected_final_result_20260905_sample_only'
check(normalizePath('final_result')==normalizePath(frozen), 'Wrong final_result target')
n_frozen <- verify_hashes(fread(file.path(frozen,'Corrected_Data/corrected_package_sha256.csv')), frozen, 'File')
n_local <- verify_hashes(fread('config/kla31_local_inputs_sha256.csv'))
audit <- 'audit/20260912_workspace_rnaseq'
g <- fread(file.path(audit,'groups_31.csv'))
r <- fread(file.path(audit,'rnaseq_reference_registry_31.csv'))
check(nrow(g)==31L && !anyDuplicated(g[,.(PXD,SampleGroup)]), 'Invalid group scope')
check(identical(g$GroupID,r$GroupID) && identical(g$PXD,r$PXD) && identical(g$SampleGroup,r$SampleGroup), 'RNA registry mismatch')
check(all(!r$AnalysisReady), 'RNA readiness changed; review eligibility evidence')
check(sum(g$PXD=='PXD064038')==1, 'ESCC missing/duplicated')
counts <- g[, .N, by=Category]
expected <- c(normal_tissue=9L,cancer_tissue=3L,cancer_cells=12L,normal_cells=7L)
check(setequal(counts$Category,names(expected)) && all(counts$N==expected[counts$Category]), 'Category scope mismatch')
validation <- fread(file.path(frozen,'Corrected_Data/correction_validation.csv'))
check(all(validation$Pass), 'Frozen correction validation includes failure')
sample_dir <- file.path(frozen,'Corrected_Data/sample_only_inputs')
expected_rows <- c(figure1_sample_boxplot_values.csv=272L,
                  figure1_mki67_ratio_sample_values.csv=264L,
                  figure1_pathway_summary_sample_boxplot_values.csv=504L)
for (f in names(expected_rows)) {
  x <- fread(file.path(sample_dir,f))
  check(nrow(x)==expected_rows[[f]] && all(x$ObservationType=='sample'), paste('Sample grain failure:',f))
}
cat(sprintf('KLA31_PREFLIGHT_PASS: 31 groups; %d frozen files and %d local input files hash-verified; sample rows 272/264/504.\n',n_frozen,n_local))
cat('Scope: workspace integrity and existing sample-only contract, not fresh statistical inference or RNA matrix validation.\n')
