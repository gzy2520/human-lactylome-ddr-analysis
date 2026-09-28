#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table);library(IRanges)})
a<-commandArgs(TRUE);stopifnot(length(a)==1);out<-a[1];if(dir.exists(out))stop('Fresh output required');dir.create(out,recursive=TRUE)
gtf<-'audit/20260927_external_material/raw/gencode.v19.annotation.gtf.gz'
g<-fread(gtf,sep='\t',header=FALSE,skip='chr1',select=c(3,4,5,9),quote='');g<-g[V3=='exon'];g[,Gene:=sub('[.].*','',sub('.*gene_id "([^"]+)".*','\\1',V9))]
len<-g[,.(ExonicLength=sum(width(reduce(IRanges(start=V4,end=V5))))),by=Gene];stopifnot(!anyDuplicated(len$Gene),all(len$ExonicLength>0));fwrite(len,file.path(out,'gencode19_exonic_lengths.csv'));rm(g);gc(FALSE)
fs<-sort(list.files('audit/20260927_external_material/rna',pattern='txt.gz$',full.names=TRUE));stopifnot(length(fs)==46)
need<-fread('outputs/20260926_full_rna_inputs/gene_ids.csv')$Ensembl;x<-matrix(NA_real_,46,length(need),dimnames=list(sub('_Patient.*','',basename(fs)),need));meta<-list();qc<-list()
for(i in seq_along(fs)){
 d<-fread(fs[i]);d[,Gene:=sub('[.].*','',Gene_id)];stopifnot(!anyDuplicated(d$Gene),all(need%in%d$Gene),all(d$Assigned_count>=0))
 d[,Length:=len$ExonicLength[match(Gene,len$Gene)]];lost<-sum(d$Assigned_count[is.na(d$Length)])/sum(d$Assigned_count);stopifnot(lost<.01,all(is.finite(d$Length[match(need,d$Gene)])))
 d[,Rate:=Assigned_count/Length];total<-sum(d$Rate,na.rm=TRUE);d[,TPM:=1e6*Rate/total];stopifnot(abs(sum(d$TPM,na.rm=TRUE)-1e6)<1e-5)
 x[i,]<-log2(d$TPM[match(need,d$Gene)]+1)
 material<-if(grepl('_cancer',fs[i]))'ESCC_untreated'else'ESCC_adjacent_untreated'
 meta[[i]]<-data.table(SampleID=rownames(x)[i],ReferenceKey='GSE130078',Material=material,DonorKey=sub('.*(Patient_[0-9]+).*','\\1',basename(fs[i])),SourceFile=fs[i],Library='rRNA-depleted total RNA',ExperimentalPerturbation='No KO/KD/OE or drug perturbation documented for tissue RNA',PretreatmentHistory='Not reported in inspected clinical table; not verified as treatment-naive')
 qc[[i]]<-data.table(SampleID=rownames(x)[i],InputGenes=nrow(d),ModelGenes=length(need),MissingModelGenes=0,GenesWithoutLength=sum(is.na(d$Length)),FractionReadsWithoutLength=lost,CountSum=sum(d$Assigned_count),TPMSum=sum(d$TPM,na.rm=TRUE))
}
s<-rbindlist(meta);old<-readRDS('outputs/20260926_full_rna_inputs/full_rna.rds')$metadata;stopifnot(!any(s$SampleID%in%old$SampleID),all(is.finite(x)))
saveRDS(list(x=x,metadata=s),file.path(out,'full_rna.rds'));fwrite(s,file.path(out,'metadata.csv'));fwrite(rbindlist(qc),file.path(out,'qc.csv'))
paths<-c(script='workflow/external_material/prepare_rna.R',gtf=gtf,fs);fwrite(data.table(Path=unname(paths),SHA256=vapply(paths,digest::digest,character(1),algo='sha256',file=TRUE)),file.path(out,'input_sha256.csv'))
cat('PASS: 46 new RNA records; all 18332 model genes present; count-to-TPM with source-era GENCODE19 exonic union lengths.\n');print(rbindlist(qc)[,.(MaxMissingReadFraction=max(FractionReadsWithoutLength),GenesWithoutLength=unique(GenesWithoutLength))])
