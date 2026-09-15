#!/usr/bin/env Rscript
args<-commandArgs(TRUE);stopifnot(length(args)==2)
root<-normalizePath(args[1]);out<-args[2];dir.create(out,recursive=TRUE,showWarnings=FALSE);set.seed(25)
sources<-list(
 list(key='HK2_untreated',groups='KLA31_27;KLA31_28',gse='GSE240748',gsm=paste0('GSM',7708662:7708664),condition='WT untreated 24h; independent culture documentation pending'),
 list(key='MCF7_DMSO7d',groups='KLA31_13',gse='GSE157383',gsm=paste0('GSM',4763907:4763909),condition='DMSO 7 days; external baseline; concentration not specified in GEO'),
 list(key='MDAMB468_DMSO7d',groups='KLA31_20',gse='GSE157383',gsm=paste0('GSM',4763917:4763919),condition='DMSO 7 days; external baseline; concentration not specified in GEO'),
 list(key='HepG2_DMSO24h',groups='KLA31_18',gse='GSE158552',gsm=paste0('GSM',4802741:4802743),condition='DMSO 24h; parental reference, not SIRT1/3 knockout'),
 list(key='HUVEC_untreated',groups='KLA31_31',gse='GSE203551',gsm=c('GSM6176178','GSM6176183','GSM6176188','GSM6176193'),condition='unstimulated untreated; donors 32/36/46/47; passage 3; not Pg infected'),
 list(key='sperm_GroupI',groups='KLA31_07',gse='GSE65683',gsm=paste0('GSM',1602977:1602983),condition='Group I birth-outcome reference; clinical eligibility and protocol review pending'),
 list(key='tendon_nondiabetic',groups='KLA31_01',gse='GSE236746',gsm=paste0('GSM',7574837:7574839),condition='nondiabetic torn tendon, not healthy tendon'),
 list(key='PC3M_DMSO',groups='KLA31_22',gse='GSE235595',gsm=paste0('GSM',7506018:7506020),condition='PC-3M DMSO 4 days; not KMI169/KMI169Ctrl'))
inventory<-list();qc<-list();samp<-list()
for(s in sources) {
 files<-list.files(file.path(root,'raw/ncbi_counts_20260915'),pattern=paste0('^',s$gse,'_raw_counts_.*tsv.gz$'),full.names=TRUE)
 stopifnot(length(files)==1)
 d<-read.delim(gzfile(files),check.names=FALSE,colClasses='character')
 gene<-d[[1]];stopifnot(all(grepl('^[0-9]+$',gene)),!anyDuplicated(gene))
 available<-intersect(s$gsm,names(d));missing<-setdiff(s$gsm,names(d))
 n<-length(available);if(n==0)stop('No selected samples: ',s$key)
 m<-as.matrix(d[,available,drop=FALSE]);storage.mode(m)<-'numeric';rownames(m)<-gene
 stopifnot(all(is.finite(m)),all(m>=0),all(m==floor(m)),!anyDuplicated(colnames(m)),all(colSums(m)>0))
 saveRDS(list(counts=m,IDType='NCBI GeneID (Entrez)',source=files,ReferenceKey=s$key,
  Condition=s$condition,Assembly='GRCh38.p13; NCBI generated counts',AnalysisReady=FALSE),file.path(out,paste0(s$key,'.rds')))
 for(j in seq_len(n))qc[[paste(s$key,j)]]<-data.frame(ReferenceKey=s$key,SampleID=colnames(m)[j],Genes=nrow(m),
  Counts=sum(m[,j]),DetectedGenes=sum(m[,j]>0),ZeroFraction=mean(m[,j]==0))
 inventory[[s$key]]<-data.frame(ReferenceKey=s$key,GroupIDs=s$groups,GSE=s$gse,RequestedN=length(s$gsm),MatrixN=n,
  MissingSamples=paste(missing,collapse=';'),Genes=nrow(m),StableID='Entrez',NumericQCPass=TRUE,
  ReplicateCountGate=n>=3,Condition=s$condition,AnalysisReady=FALSE,InputMD5=unname(tools::md5sum(files)))
 samp[[s$key]]<-data.frame(ReferenceKey=s$key,GSE=s$gse,GSM=available,Condition=s$condition)
 if(n>=2)write.csv(cor(m,method='spearman'),file.path(out,paste0(s$key,'_spearman.csv')))
}
write.csv(do.call(rbind,inventory),file.path(out,'resolved_reference_inventory.csv'),row.names=FALSE)
write.csv(do.call(rbind,samp),file.path(out,'selected_samples.csv'),row.names=FALSE)
write.csv(do.call(rbind,qc),file.path(out,'sample_qc.csv'),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
cat('NCBI_REFERENCE_MATRIX_QC_COMPLETE\n');print(do.call(rbind,inventory)[,1:7])
