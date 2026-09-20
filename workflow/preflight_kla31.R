#!/usr/bin/env Rscript
# Current release check. The pre-20260916 workspace gate is archived separately.
suppressPackageStartupMessages({library(data.table); library(digest)})
set.seed(25)
arg <- grep('^--file=', commandArgs(FALSE), value=TRUE)
root <- normalizePath(file.path(dirname(sub('^--file=', '', arg[[1L]])), '..'))
setwd(root)
verify <- function(tab, base='.', column='Path') {
  p <- file.path(base,tab[[column]])
  if (!all(file.exists(p))) stop('Missing manifest files: ',paste(p[!file.exists(p)],collapse=', '))
  actual <- vapply(p,function(f)digest(file=f,algo='sha256'),character(1))
  if (!identical(unname(actual),unname(tab$SHA256))) stop('SHA256 mismatch: ',paste(p[actual!=tab$SHA256],collapse=', '))
  nrow(tab)
}
frozen <- 'corrected_final_result_20260905_sample_only'
n_frozen <- verify(fread(file.path(frozen,'Corrected_Data/corrected_package_sha256.csv')), frozen,'File')
release <- fread('config/rna_current_release.csv')
stopifnot(!anyDuplicated(release$Role),all(c('expression','qsmooth','panel','assisted','delivery')%in%release$Role))
rpath <- function(role) release[Role==role,Path][[1L]]
n_input <- verify(fread('config/kla31_current_inputs_sha256.csv'))
n_release <- verify(fread(file.path(rpath('delivery'),'release_sha256.csv')))
q <- fread(file.path(rpath('expression'),'group_sample_qc.csv'))
g <- fread(file.path(rpath('qsmooth'),'group_expansion_31.csv'))
v <- fread(file.path(rpath('delivery'),'validation.csv'))
p <- fread('data/publication_input/group_summary_31.csv')
expected <- c(normal_tissue=9L,cancer_tissue=3L,cancer_cells=12L,normal_cells=7L)
category_counts <- p[,.N,by=Category]
stopifnot(nrow(p)==31L,uniqueN(p[,.(PXD,SampleGroup)])==31L,
          sum(p$PXD=='PXD064038')==1L,
          setequal(category_counts$Category,names(expected)),
          all(category_counts$N==expected[category_counts$Category]))
stopifnot(nrow(q)==1898L,uniqueN(q$SampleID)==1898L,uniqueN(q$ReferenceKey)==28L,
          nrow(g)==31L,uniqueN(g$GroupID)==31L,uniqueN(g$ReferenceKey)==28L,all(v$Pass))
for (p in c('results','renv/library')) stopifnot(dir.exists(p),!nzchar(Sys.readlink(p)))
cat(sprintf('KLA31_PREFLIGHT_PASS: %d frozen files; %d current inputs; %d release files; 28 RNA materials / 1898 samples / 31 mappings.\n',n_frozen,n_input,n_release))
cat('Scope: release integrity and numerical pipeline checks. Statistical method unchanged at user request; not an inference endorsement.\n')
