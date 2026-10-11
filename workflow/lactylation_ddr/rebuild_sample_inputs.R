#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
set.seed(25)
audit <- 'audit/20261011_ddr_sample_recovery'
out <- 'deliverables/20261011_lactylation_ddr_sample_rebuilt'
read <- function(p) {f<-file.path(audit,p); if(grepl('gz$',f)) fread(cmd=paste('gzip -dc',shQuote(f))) else fread(f)}
keys <- c('PXD','SampleGroup','Dataset','SampleID')
cols <- c(keys,'BaseAccession','SourceFile','UnitType','Condition','SemEligible','Evidence')
b <- read('base_rebuild/all_membership.csv.gz')[ObservationType=='sample' & !(PXD=='PXD028488' & SampleGroup=='TALL-104' & Dataset=='Lactylome (Kla)')]
b[,`:=`(UnitType=ifelse(SampleClass=='model','model','source_sample'),Condition=SampleClass,SemEligible=SampleClass!='model',Evidence='Deposited source columns; see base_rebuild/rebuild.R')]
b <- b[,..cols]
e <- read('base_rebuild/escc_kla_membership.csv.gz')
e[,`:=`(SemEligible=TRUE,Evidence='Positive sample-specific localization evidence in each MEC/NEC source column')]
e <- e[,..cols]
x <- read('existing_samples/recomputed_membership.csv.gz')
t <- x[Case=='TALL_Kla_legacy_mapping']
t[,`:=`(PXD='PXD028488',SampleGroup='TALL-104',Dataset='Lactylome (Kla)',UnitType='source_sample',Condition='NaLAC',SemEligible=TRUE,Evidence=Method)]
t <- t[,..cols]
s <- x[Case=='ESCC_proteome']
s[,`:=`(PXD='PXD064038',SampleGroup='MEC and NEC ESCC groups',Dataset='Whole proteome',UnitType='source_sample',Condition='ESCC tumor tissue',SemEligible=TRUE,Evidence=Method)]
s <- s[,..cols]
k <- read('kla/recovered_kla_membership.csv.gz')
k[,`:=`(Dataset='Lactylome (Kla)',Evidence=UnitEvidence)]
k <- k[,..cols]
w <- read('whole_proteome/recovered_protein_membership.csv.gz')[UseRole=='primary']
wm <- read('whole_proteome/recovered_unit_values.csv')[UseRole=='primary']
w <- merge(w,wm[,.(PXD,SampleGroup,SampleID,Dataset,SourceFile,UnitType=ObservationType,Evidence=UnitDefinition)],by=c('PXD','SampleGroup','SampleID'))
w[,`:=`(Condition='source material',SemEligible=UnitType %in% c('technical_preparation','sample'))]
w <- w[,..cols]
a <- rbindlist(list(b,e,t,s,k,w),use.names=TRUE)
stopifnot(!anyDuplicated(a[,c(keys,'BaseAccession'),with=FALSE]))
g <- fread('data/publication_input/human_ddr_go_annotations.tsv')
g <- g[(is.na(`TAXON ID`)|`TAXON ID`==9606) & `GENE PRODUCT DB`=='UniProtKB' & !grepl('(^|\\|)NOT($|\\|)',QUALIFIER)]
ddr <- unique(sub('-[0-9]+$','',g[['GENE PRODUCT ID']]))
ddr <- ddr[grepl('^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})$',ddr)]
stopifnot(length(ddr)==2660)
a[,IsDdr:=BaseAccession %in% ddr]
u <- a[,.(ProteinCount=.N,DdrProteinCount=sum(IsDdr),SourceFile=first(SourceFile),UnitType=first(UnitType),Condition=first(Condition),SemEligible=first(SemEligible),Evidence=first(Evidence)),by=keys]
u[,Fraction:=100*DdrProteinCount/ProteinCount]
# Two singleton conditions supply no within-condition replication.
cnt <- u[,.(ConditionN=.N),by=.(PXD,SampleGroup,Dataset,Condition)]
no_rep <- cnt[,.(NoReplication=all(ConditionN==1),Conditions=.N),by=.(PXD,SampleGroup,Dataset)][NoReplication & Conditions>1]
u[no_rep,on=.(PXD,SampleGroup,Dataset),SemEligible:=FALSE]
m <- fread(file.path(out,'tables/display_mapping_31.csv'))
u <- merge(u,m[,.(PXD,SampleGroup,GroupID,X,PlotLabel)],by=c('PXD','SampleGroup'),all.x=TRUE)
stopifnot(nrow(u)==325,!anyNA(u$X),uniqueN(u[,.(GroupID,Dataset)])==62)
ps <- u[,.(N=.N,Mean=mean(Fraction),SD=if(.N>1)sd(Fraction) else NA_real_,SEM=if(.N>1)sd(Fraction)/sqrt(.N) else NA_real_,SemEligible=all(SemEligible),UnitTypes=paste(sort(unique(UnitType)),collapse=';'),Conditions=paste(sort(unique(Condition)),collapse=';')),by=.(PXD,SampleGroup,GroupID,X,Dataset)]
ps[,PlotSEM:=ifelse(SemEligible & N>1,SEM,NA_real_)]
ps[,UnitMark:=fcase(grepl('technical_preparation',UnitTypes),'t',grepl('model',UnitTypes),'m',grepl('qc_source_pool',UnitTypes),'q',default='')]
setorder(u,X,Dataset,SampleID);setorder(ps,X,Dataset)
fwrite(a,file.path(out,'inputs/ddr_protein_membership.csv.gz'))
fwrite(u,file.path(out,'tables/ddr_source_unit_values.csv'))
fwrite(ps,file.path(out,'tables/ddr_mean_sem_62.csv'))
fwrite(cnt,file.path(out,'tables/ddr_condition_counts.csv'))
condition_summary <- u[,.(N=.N,Mean=mean(Fraction),SD=if(.N>1)sd(Fraction) else NA_real_,SEM=if(.N>1 && all(SemEligible))sd(Fraction)/sqrt(.N) else NA_real_,UnitTypes=paste(sort(unique(UnitType)),collapse=';')),by=.(PXD,SampleGroup,Dataset,Condition)]
fwrite(condition_summary,file.path(out,'tables/ddr_condition_summary.csv'))
fwrite(data.table(BaseAccession=sort(ddr)),file.path(out,'inputs/frozen_ddr_accessions.csv'))
old <- fread('deliverables/20261011_lactylation_ddr/inputs/proteome_samples.csv')
# Original TALL sample IDs are raw-file labels in the source; replace in comparison only.
old[,OriginalSampleID:=SampleID]
old[PXD=='PXD028488' & SampleGroup=='TALL-104' & Dataset=='Lactylome (Kla)',SampleID:=paste0('TALL-NaLAC-enrichment-',sub('Sample ','',SampleID),'.raw')]
v <- merge(old[,c(keys,'OriginalSampleID','ProteinCount','DdrProteinCount','DdrFractionPercentage'),with=FALSE],u[,c(keys,'ProteinCount','DdrProteinCount','Fraction'),with=FALSE],by=keys,all.x=TRUE,suffixes=c('_Original','_Rebuilt'))
v[,ExactCounts:=ProteinCount_Original==ProteinCount_Rebuilt & DdrProteinCount_Original==DdrProteinCount_Rebuilt]
stopifnot(nrow(v)==272,!anyNA(v$ExactCounts),sum(v$ExactCounts)==269)
fwrite(v,file.path(out,'tables/original_sample_comparison.csv'))
oldsum <- old[,.(OriginalMean=mean(DdrFractionPercentage),OriginalN=.N),by=.(PXD,SampleGroup,Dataset)]
comparison <- merge(oldsum,ps,by=c('PXD','SampleGroup','Dataset'),all.y=TRUE)
comparison[,Difference:=Mean-OriginalMean]
fwrite(comparison,file.path(out,'tables/original_figure_comparison.csv'))
print(ps[,.(Positions=.N,Units=sum(N)),by=.(UnitTypes,SemEligible)])
cat('325 source-unit display records; 62 bars; ',sum(!is.na(ps$PlotSEM)),' error bars.\n',sep='')
