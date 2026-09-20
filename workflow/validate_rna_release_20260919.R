#!/usr/bin/env Rscript
# Numerical and scope checks, deliberately not a statistical-method change.
suppressPackageStartupMessages({library(data.table); library(digest)})
args <- commandArgs(TRUE); stopifnot(length(args)==5L)
expr <- args[[1]]; qs <- args[[2]]; panel <- args[[3]]; assisted <- args[[4]]; delivery <- args[[5]]
checks <- list()
check <- function(name, ok) {
  checks[[name]] <<- data.table(Check=name,Pass=isTRUE(ok))
}
readmat <- function(p) { x<-fread(p); m<-as.matrix(x[,-1]); rownames(m)<-x[[1]];m }
q<-fread(file.path(expr,'group_sample_qc.csv')); man<-fread(file.path(expr,'group_manifest.csv'))
g<-fread(file.path(qs,'group_expansion_31.csv'))
check('28 reference matrices / 1898 unique samples',nrow(man)==28L && nrow(q)==1898L && uniqueN(q$SampleID)==1898L)
check('31 mappings / 28 reference keys',nrow(g)==31L && uniqueN(g$GroupID)==31L && uniqueN(g$ReferenceKey)==28L)
check('TPM totals and valid IDs in server verification',nrow(fread(file.path(expr,'verification_summary.csv')))==28L && all(abs(q$TPMTotal-1e6)<1))
u<-readmat(file.path(qs,'matrices/unified_full_log2tpm.tsv.gz'))
b<-readmat(file.path(qs,'matrices/qsmooth_B_full_log2tpm.tsv.gz'))
a<-readmat(file.path(qs,'matrices/qsmooth_A_collapsed_log2tpm.tsv.gz'))
check('matrix shape and finite values',identical(dim(b),c(18332L,1898L)) && identical(dim(a),c(18332L,28L)) && all(is.finite(b)) && all(is.finite(a)))
check('stable unique Ensembl IDs',!anyDuplicated(rownames(b)) && all(grepl('^ENSG[0-9]+(_PAR_Y)?$',rownames(b))))
check('sample metadata complete',!anyDuplicated(colnames(b)) && setequal(colnames(b),q$SampleID))
changes<-list()
for(k in man$ReferenceKey) {
  m<-readmat(file.path(expr,'matrices',paste0(k,'_log2tpm.tsv.gz')))
  prev<-readmat(file.path('outputs/20260916_expression_extraction/matrices',paste0(k,'_log2tpm.tsv.gz')))
  check(paste('source-to-qsmooth input',k),max(abs(m[rownames(u),,drop=FALSE]-u[,colnames(m),drop=FALSE]))<1e-10)
  same_keys<-identical(dimnames(m),dimnames(prev))
  delta<-if(same_keys)max(abs(m-prev)) else NA_real_
  changes[[k]]<-data.table(ReferenceKey=k,SameKeys=same_keys,MaxAbsLog2Change=delta)
  expected_change<-k%in%c('GSE132714_BPH','GSE163787_TALL104')
  check(paste('expected source correction scope',k),same_keys && is.finite(delta) && if(expected_change)delta>0.1 else delta<1e-10)
}
fwrite(rbindlist(changes),file.path(delivery,'source_changes.csv'))
core<-file.path(delivery,'core'); ref<-file.path(delivery,'reference')
for (nm in c('fixed_ddr','lactylated_gene')) {
  x<-fread(file.path(core,paste0('figure1a_',nm,'_28materials.csv')))
  s<-fread(file.path(core,paste0('figure1a_',nm,'_sample_values.csv')))
  check(paste(nm,'28 materials'),nrow(x)==28L&&uniqueN(x$ReferenceKey)==28L)
  check(paste(nm,'1898 samples'),nrow(s)==1898L&&uniqueN(s$SampleID)==1898L)
}
x<-fread(file.path(core,'figure1b_hallmark_proliferation_28materials.csv'))
s<-fread(file.path(core,'figure1b_hallmark_proliferation_sample_values.csv'))
check('proliferation scope',nrow(x)==28L&&uniqueN(x$ReferenceKey)==28L&&nrow(s)==1898L&&uniqueN(s$SampleID)==1898L)
r<-fread(file.path(core,'regulator_rna_qsmooth_28materials.csv'))
check('regulators 49 role entries x 28 references',nrow(r)==49L*28L&&uniqueN(r$ReferenceKey)==28L&&all(is.finite(r$QsmoothLog2TPM))&&all(is.finite(r$QsmoothZScore)))
dual<-fread(file.path(core,'figure1a_rna_ddr_and_lactylated_28materials.csv'))
check('dual table 56 rows, two scores per material',nrow(dual)==56L&&all(dual[,.N,by=GroupID]$N==2L))
for(d in c(core,ref)) {
  pdfs<-list.files(d,pattern='\\.pdf$',full.names=TRUE)
  check(paste(basename(d),'PDF inventory'),length(pdfs)==if(d==core)14L else 4L)
  check(paste(basename(d),'PNG pairs'),all(file.exists(sub('\\.pdf$','.png',pdfs))))
  for(p in pdfs) {
    txt<-system2('pdftotext',c(shQuote(p),'-'),stdout=TRUE)
    check(paste('no missing statistics label',basename(p)),!any(grepl('(p|F) = NA',txt)))
  }
}
check('no historical 31-group core deliverables',!any(grepl('31groups|percentiles',list.files(core))))
out<-rbindlist(checks);fwrite(out,file.path(delivery,'validation.csv'))
if(!all(out$Pass)){ print(out[Pass==FALSE]); stop('Release validation failed') }
cat(sprintf('RNA_RELEASE_VALIDATION_PASS: %d checks\n',nrow(out)))
