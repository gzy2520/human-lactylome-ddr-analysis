#!/usr/bin/env Rscript
args <- commandArgs(TRUE); stopifnot(length(args)==2)
root <- normalizePath(args[1]); out <- args[2]
dir.create(out, recursive=TRUE, showWarnings=FALSE)
specs <- list(
 list(key='HCT116_control', groups='KLA31_14;KLA31_15;KLA31_16', gse='GSE253699', file='GSE253699_raw_counts.txt.gz', cols=paste0('HCT_NC_',1:3), gsm=paste0('GSM',8026675:8026677), condition='WT control; IMEM 10% FBS; not coculture; culture independence documentation pending'),
 list(key='RKO_DMSO', groups='KLA31_24', gse='GSE318640', file='GSE318640_Expression_matrix.txt.gz', cols=paste0('C',1:3,'_count'), gsm=paste0('GSM',9499692:9499694), condition='DMSO control; duration and concentration unreported; parental not GSK3B knockout'))
inv <- qc <- samples <- list()
for(s in specs) {
 f <- file.path(root,'raw/replicate_replacements',s$gse,s$file)
 d <- read.delim(gzfile(f), check.names=FALSE, quote='', comment.char='')
 stopifnot(all(grepl('^ENSG[0-9]+(\\.[0-9]+)?$',d$gene_id)),!anyDuplicated(d$gene_id),all(s$cols %in% names(d)))
 m <- as.matrix(d[,s$cols,drop=FALSE]); storage.mode(m)<-'numeric'
 rownames(m)<-d$gene_id; colnames(m)<-s$gsm
 stopifnot(all(is.finite(m)),all(m>=0),all(m==floor(m)),all(colSums(m)>0))
 saveRDS(list(counts=m,IDType='Ensembl',source=f,ReferenceKey=s$key,Condition=s$condition,AnalysisReady=FALSE),file.path(out,paste0(s$key,'.rds')))
 inv[[s$key]]<-data.frame(ReferenceKey=s$key,GroupIDs=s$groups,GSE=s$gse,RequestedN=3,MatrixN=3,MissingSamples='',Genes=nrow(m),StableID='Ensembl',NumericQCPass=TRUE,ReplicateCountGate=TRUE,Condition=s$condition,AnalysisReady=FALSE,InputMD5=unname(tools::md5sum(f)))
 qc[[s$key]]<-data.frame(ReferenceKey=s$key,SampleID=colnames(m),Genes=nrow(m),Counts=colSums(m),DetectedGenes=colSums(m>0),ZeroFraction=colMeans(m==0))
 samples[[s$key]]<-data.frame(ReferenceKey=s$key,GSE=s$gse,GSM=s$gsm,SourceColumn=s$cols,Condition=s$condition)
 write.csv(cor(m,method='spearman'),file.path(out,paste0(s$key,'_spearman.csv')))
}
write.csv(do.call(rbind,inv),file.path(out,'native_reference_inventory.csv'),row.names=FALSE)
write.csv(do.call(rbind,qc),file.path(out,'native_sample_qc.csv'),row.names=FALSE)
write.csv(do.call(rbind,samples),file.path(out,'native_selected_samples.csv'),row.names=FALSE)
cat('NATIVE_REFERENCE_MATRIX_QC_COMPLETE\n')
