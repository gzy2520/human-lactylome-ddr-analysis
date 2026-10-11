suppressPackageStartupMessages(library(data.table))
set.seed(25)
old<-'deliverables/20261011_lactylation_ddr_tall_verified'
out<-'deliverables/20261011_lactylation_ddr_all_groups_audited'
aud<-'audit/20261011_all_groups_algorithm_audit'
readgz<-function(p)fread(cmd=paste('gzip -dc',shQuote(p)))
k<-c('PXD','SampleGroup','Dataset','SampleID')
a<-readgz(file.path(old,'inputs/ddr_protein_membership.csv.gz'))
u0<-fread(file.path(old,'tables/ddr_source_unit_values.csv'))
meta<-u0[,setdiff(names(u0),c('ProteinCount','DdrProteinCount','Fraction')),with=FALSE]
meta[,`:=`(SourceObservationIDs=SampleID,SourceObservationN=1L)]
fix<-readgz(file.path(aud,'kla_cells/corrected_affected_membership.csv.gz'))
fixkeys<-unique(fix[,..k])
# Replace source sets after row-level decoy filtering; preserve the audit baseline.
repl<-merge(fix[,c(k,'BaseAccession'),with=FALSE],a[,lapply(.SD,first),by=k,.SDcols=setdiff(names(a),c(k,'BaseAccession','IsDdr'))],by=k)
repl[,IsDdr:=BaseAccession %in% fread(file.path(out,'inputs/frozen_ddr_accessions.csv'))$BaseAccession]
a<-rbindlist(list(a[!fixkeys,on=k],repl),use.names=TRUE)
meta[fixkeys,on=k,Evidence:='Independent source reparse: reject raw reverse representative before normalizing IDs; other forward support checked']
# Complete author-processed matrices provide one protein panel, not independent detection masks.
panels<-data.table(GroupID=c('KLA31_05','KLA31_05','KLA31_09','KLA31_12'),Dataset=c('Lactylome (Kla)','Whole proteome','Whole proteome','Whole proteome'))
for(i in seq_len(nrow(panels))){
 mm<-meta[panels[i],on=.(GroupID,Dataset)]
 aa<-a[mm[,..k],on=k]
 sig<-aa[,.(Signature=paste(sort(BaseAccession),collapse=';')),by=SampleID]
 stopifnot(uniqueN(sig$Signature)==1)
 new_id<-paste0(panels$GroupID[i],'_author_processed_matrix')
 aa<-aa[SampleID==mm$SampleID[1]];aa[,`:=`(SampleID=new_id,UnitType='processed_matrix',Condition='author processed matrix',SemEligible=FALSE,Evidence='One author-processed complete protein panel; source sample columns retained in audit; original per-sample detection mask unavailable')]
 a<-rbindlist(list(a[!mm[,..k],on=k],aa),use.names=TRUE)
 nm<-copy(mm[1]);nm[,`:=`(SampleID=new_id,UnitType='processed_matrix',Condition='author processed matrix',SemEligible=FALSE,Evidence=aa$Evidence[1],SourceObservationIDs=paste(mm$SampleID,collapse=';'),SourceObservationN=nrow(mm))]
 meta<-rbindlist(list(meta[!mm[,..k],on=k],nm),use.names=TRUE)
}
meta[GroupID=='KLA31_28'&Dataset=='Lactylome (Kla)',Condition:=ifelse(grepl('^ctr',SampleID),'ctr','man')]
a[PXD=='PXD078736'&Dataset=='Lactylome (Kla)',Condition:=ifelse(grepl('^ctr',SampleID),'ctr','man')]
u<-a[,.(ProteinCount=uniqueN(BaseAccession),DdrProteinCount=uniqueN(BaseAccession[IsDdr])),by=k]
u<-merge(u,meta,by=k);u[,Fraction:=100*DdrProteinCount/ProteinCount]
m<-fread(file.path(out,'tables/display_mapping_31.csv'))
m[GroupID=='KLA31_28',PlotLabel:='HK-2 (ctr/man)*']
u[GroupID=='KLA31_28',PlotLabel:='HK-2 (ctr/man)*']
fwrite(m,file.path(out,'tables/display_mapping_31.csv'))
s<-u[,.(N=.N,Mean=mean(Fraction),SD=if(.N>1)sd(Fraction)else NA_real_,SEM=if(.N>1)sd(Fraction)/sqrt(.N)else NA_real_,SemEligible=all(SemEligible),UnitTypes=paste(sort(unique(UnitType)),collapse=';'),Conditions=paste(sort(unique(Condition)),collapse=';')),by=.(PXD,SampleGroup,GroupID,X,Dataset)]
s[,PlotSEM:=ifelse(N>1&SemEligible,SEM,NA_real_)]
s[,UnitMark:=fcase(UnitTypes=='processed_matrix','p',grepl('technical_preparation',UnitTypes),'t',grepl('model',UnitTypes),'m',grepl('qc_source_pool',UnitTypes),'q',default='')]
setorder(u,X,Dataset,SampleID);setorder(s,X,Dataset)
stopifnot(nrow(u)==307,nrow(s)==62,sum(!is.na(s$PlotSEM))==36,!anyDuplicated(a[,c(k,'BaseAccession'),with=FALSE]))
fwrite(a,file.path(out,'inputs/ddr_protein_membership.csv.gz'))
fwrite(u,file.path(out,'tables/ddr_source_unit_values.csv'))
fwrite(s,file.path(out,'tables/ddr_mean_sem_62.csv'))
c<-u[,.(N=.N,Mean=mean(Fraction),SD=if(.N>1)sd(Fraction)else NA_real_,SEM=if(.N>1&&all(SemEligible))sd(Fraction)/sqrt(.N)else NA_real_),by=.(GroupID,PXD,SampleGroup,Dataset,Condition)]
fwrite(c,file.path(out,'tables/ddr_condition_summary.csv'))
v<-merge(fread(file.path(old,'tables/ddr_mean_sem_62.csv')),s,by=c('GroupID','Dataset'),suffixes=c('_Before','_After'))
v[,MeanChange:=Mean_After-Mean_Before]
fwrite(v,file.path(out,'tables/all_groups_before_after.csv'))
fwrite(u[UnitType=='processed_matrix'],file.path(out,'tables/processed_matrix_panels.csv'))
cat('307 analysis units; 62 bars; 36 SEM bars; 4 processed-matrix panels.\n')
