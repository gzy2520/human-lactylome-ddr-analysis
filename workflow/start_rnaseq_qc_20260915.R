#!/usr/bin/env Rscript
# Source-local RNA preparation. No cross-source abundance comparison.
args <- commandArgs(TRUE)
stopifnot(length(args) == 2L)
root <- normalizePath(args[1], mustWork=TRUE)
out <- args[2]
if (dir.exists(out)) stop('Output already exists; choose a new run directory')
dir.create(out, recursive=TRUE)
set.seed(25)
write_tsv <- function(x, name) write.table(x, file.path(out,name), sep='\t', row.names=FALSE, quote=TRUE, na='NA')
qc <- list()
save_matrix <- function(m, key, unit) {
  stopifnot(is.numeric(m), !anyDuplicated(rownames(m)), !anyDuplicated(colnames(m)))
  saveRDS(list(expression=m, unit=unit, ReferenceKey=key), file.path(out,paste0(key,'.rds')))
  qc[[key]] <<- data.frame(ReferenceKey=key, SampleID=colnames(m), Unit=unit,
    Genes=nrow(m), Missing=colSums(is.na(m)), Nonfinite=colSums(!is.finite(m)),
    Negative=colSums(m<0,na.rm=TRUE), Detected=colSums(m>0,na.rm=TRUE))
  if(ncol(m)>1) write_tsv(as.data.frame(as.table(cor(m, method='spearman', use='pairwise.complete.obs'))),paste0(key,'_spearman.tsv'))
}
base <- file.path(root,'raw/depmap/DepMap_Public_24Q4')
targets <- c('ACH-000019','ACH-000971','ACH-000739','ACH-000681','ACH-000849','ACH-000147','ACH-000943')
con <- file(file.path(base,'OmicsExpressionProteinCodingGenesTPMLogp1.csv'),'r')
header <- readLines(con,n=1)
selected <- character()
repeat {
  lines <- readLines(con,n=100,warn=FALSE)
  if(!length(lines)) break
  ids <- gsub('"','',sub(',.*','',lines))
  selected <- c(selected,lines[ids %in% targets])
}
close(con)
d <- read.csv(text=paste(c(header,selected),collapse='\n'),check.names=FALSE)
stopifnot(nrow(d)==length(targets),setequal(d[[1]],targets))
ids <- sub('^.*\\(([0-9]+)\\)$','\\1',names(d)[-1])
stopifnot(all(grepl('^[0-9]+$',ids)))
dup <- duplicated(ids)|duplicated(ids,fromLast=TRUE)
write_tsv(data.frame(SourceColumn=names(d)[-1],EntrezID=ids,ExcludedAmbiguousDuplicate=dup),'depmap_gene_id_audit.tsv')
# Exclude every ambiguous repeated Entrez ID; never select or average silently.
m <- t(as.matrix(d[,-1,drop=FALSE]))
rownames(m) <- ids; colnames(m) <- d[[1]]
save_matrix(m[!dup,,drop=FALSE],'DepMap_Public_24Q4_unique_Entrez','TPMLogp1_as_released')
model <- read.csv(file.path(base,'Model.csv'),check.names=FALSE)
profile <- read.csv(file.path(base,'OmicsProfiles.csv'),check.names=FALSE)
write_tsv(model[model$ModelID %in% targets,,drop=FALSE],'depmap_selected_models.tsv')
write_tsv(profile[profile$ModelID %in% targets & profile$Datatype=='rna',,drop=FALSE],'depmap_selected_rna_profiles.tsv')
enc <- list.files(file.path(root,'raw/encode/ENCSR000COZ'),pattern='tsv$',full.names=TRUE)
e <- lapply(enc,read.delim,check.names=FALSE)
stopifnot(length(e)==2,identical(e[[1]]$gene_id,e[[2]]$gene_id))
keep <- grepl('^ENSG[0-9]+(\\.[0-9]+)?$',e[[1]]$gene_id)
write_tsv(data.frame(GeneID=e[[1]]$gene_id,IncludedHumanEnsembl=keep),'encode_gene_id_audit.tsv')
for(unit in c('expected_count','TPM','FPKM')) {
  mat <- do.call(cbind,lapply(e,function(x)x[[unit]][keep]))
  rownames(mat) <- e[[1]]$gene_id[keep]
  colnames(mat) <- sub('_GRCh38.*','',basename(enc))
  save_matrix(mat,paste0('ENCSR000COZ_',unit),unit)
}
write_tsv(do.call(rbind,qc),'matrix_sample_qc.tsv')
files <- c(list.files(base,full.names=TRUE),enc)
hash <- vapply(files,function(f)strsplit(system2('sha256sum',shQuote(f),stdout=TRUE),' ')[[1]][1],character(1))
write_tsv(data.frame(File=files,SHA256=hash),'input_sha256.tsv')
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
writeLines('MATRIX_PREPARATION_COMPLETE; biological eligibility and downstream integration remain pending',file.path(out,'COMPLETE.txt'))
cat('MATRIX_PREPARATION_COMPLETE\n')
