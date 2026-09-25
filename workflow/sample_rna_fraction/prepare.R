#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
set.seed(25)
args <- commandArgs(TRUE)
stopifnot(length(args) == 2L)
root <- normalizePath(args[1]); out <- args[2]
if (dir.exists(out) && length(list.files(out, all.files=TRUE, no..=TRUE))) stop('Output exists')
dir.create(out, recursive=TRUE)
exprdir <- file.path(root, 'outputs/20260919_expression_corrected')
inputs <- c(script='workflow/sample_rna_fraction/prepare.R',
  qc=file.path(exprdir,'group_sample_qc.csv'),
  source_manifest=file.path(exprdir,'group_manifest.csv'),
  previous_rna='outputs/20260920_lactate_metabolism/tables/material_gene_expression.csv',
  labels='outputs/20260923_rna_only_kla_fraction/material_labels_28.csv')
qc <- fread(inputs['qc']); sm <- fread(inputs['source_manifest'])
old <- fread(inputs['previous_rna']); ids <- sort(unique(old$Ensembl))
lab <- fread(inputs['labels']); setorder(lab,ReferenceKey)
stopifnot(length(ids)==128, nrow(qc)==1898, !anyDuplicated(qc$SampleID), nrow(lab)==28)
rows <- list(); checks <- list()
for (k in lab$ReferenceKey) {
  f <- file.path(exprdir,'matrices',paste0(k,'_log2tpm.tsv.gz'))
  inputs[paste0('sample_matrix_',k)] <- f
  x <- fread(f); setnames(x,1,'Ensembl')
  samples <- qc[ReferenceKey==k,SampleID]
  stopifnot(all(samples %in% names(x)),all(ids %in% x$Ensembl),!anyDuplicated(x$Ensembl))
  a <- as.matrix(x[match(ids,Ensembl),..samples]); storage.mode(a)<-'double'
  stopifnot(all(is.finite(a)))
  w <- as.data.table(t(a)); setnames(w,ids)
  w[, `:=`(SampleID=samples,ReferenceKey=k)]
  setcolorder(w,c('ReferenceKey','SampleID',ids)); rows[[k]]<-w
  target <- old[ReferenceKey==k][match(ids,Ensembl),UnsmoothedLog2TPM]
  checks[[k]] <- data.table(ReferenceKey=k,Samples=length(samples),
    MaxMedianDifference=max(abs(apply(a,1,median)-target)))
  rm(x,a); gc(FALSE)
}
wide <- rbindlist(rows); check <- rbindlist(checks)
stopifnot(nrow(wide)==1898,!anyDuplicated(wide$SampleID),max(check$MaxMedianDifference)<1e-8)
fwrite(wide,file.path(out,'sample_rna128.csv.gz'))
fwrite(check,file.path(out,'sample_median_crosscheck.csv'))
metadata <- merge(qc[,.(ReferenceKey,SampleID)],lab,by='ReferenceKey')
metadata <- merge(metadata,sm[,.(ReferenceKey,SourceKind,SourceNote)],by='ReferenceKey')
metadata[,OriginalConnGroup:=ConnGroup]
metadata[,DonorKey:=NA_character_]
metadata[grepl('^TCGA-',SampleID),DonorKey:=substr(SampleID,1,12)]
metadata[grepl('^GTEX-',SampleID),DonorKey:=vapply(strsplit(SampleID,'-',fixed=TRUE),
  function(x)paste(x[1:2],collapse='-'),character(1))]
metadata[,DonorKnown:=!is.na(DonorKey)]
shared <- metadata[DonorKnown==TRUE,.(Materials=uniqueN(ReferenceKey),Blocks=uniqueN(ConnGroup),
  BlockIDs=paste(sort(unique(ConnGroup)),collapse=';')),by=DonorKey]
fwrite(shared[Materials>1],file.path(out,'known_shared_donors.csv'))
# Union the existing source blocks when an identifiable donor crosses them.
parents <- seq_len(max(metadata$ConnGroup))
rootof <- function(i){while(parents[i]!=i)i<-parents[i];i}
for (ids0 in strsplit(shared[Blocks>1,BlockIDs],';',fixed=TRUE)) {
  ids0 <- as.integer(ids0); rr <- vapply(ids0,rootof,integer(1))
  parents[rr] <- min(rr)
}
oldblocks<-sort(unique(metadata$ConnGroup));roots<-vapply(oldblocks,rootof,integer(1))
lookup<-data.table(OriginalConnGroup=oldblocks,ConnGroup=match(roots,sort(unique(roots))))
metadata[,ConnGroup:=NULL];metadata<-merge(metadata,lookup,by='OriginalConnGroup')
fwrite(lookup,file.path(out,'source_block_contract.csv'))
metadata[DonorKnown==FALSE,DonorKey:=paste0('Unresolved:',SampleID)]
metadata[, MaterialSampleN:=.N,by=ReferenceKey]
metadata[, DonorRecordN:=.N,by=.(ReferenceKey,DonorKey)]
metadata[, MaterialDonorUnits:=uniqueN(DonorKey),by=ReferenceKey]
metadata[, SourceMaterials:=uniqueN(ReferenceKey),by=ConnGroup]
metadata[, MaterialBalancedWeight:=1/(DonorRecordN*MaterialDonorUnits)]
metadata[, SourceBalancedWeight:=MaterialBalancedWeight/SourceMaterials]
stopifnot(max(abs(metadata[,.(W=sum(SourceBalancedWeight)),by=ConnGroup]$W-1))<1e-10)
fwrite(metadata,file.path(out,'sample_metadata.csv'))
fwrite(sm,file.path(out,'rna_source_inventory.csv'))
lab[,OriginalConnGroup:=ConnGroup];lab[,ConnGroup:=NULL]
lab<-merge(lab,lookup,by='OriginalConnGroup');setorder(lab,ReferenceKey)
fwrite(lab,file.path(out,'material_labels_28.csv'))
fwrite(data.table(Role=names(inputs),Path=unname(inputs),
  SHA256=vapply(inputs,digest::digest,character(1),algo='sha256',file=TRUE,USE.NAMES=FALSE)),
  file.path(out,'input_sha256.csv'))
cat('PASS: sample RNA',nrow(wide),'x',length(ids),'; donor-aware source blocks',
    uniqueN(metadata$ConnGroup),'; unresolved donor identities',sum(!metadata$DonorKnown),'\n')

writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
