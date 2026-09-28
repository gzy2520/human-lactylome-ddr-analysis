#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table);library(readxl)})
args<-commandArgs(TRUE);stopifnot(length(args)==1);out<-args[1];if(dir.exists(out))stop('Fresh directory required');dir.create(out,recursive=TRUE)
root<-'audit/20260927_external_material/tables'
p<-as.data.table(read_excel(file.path(root,'mmc3.xlsx'),skip=1));k<-as.data.table(read_excel(file.path(root,'mmc4.xlsx'),skip=1))
ann<-fread('outputs/20260927_ddr_fraction_labels/protein_target_membership.csv')
pcols<-grep('^LFQ.intensity.',names(p),value=TRUE);kcols<-grep('^Intensity[.]',names(k),value=TRUE)
stopifnot(sum(rowSums(as.matrix(p[,..pcols])>0,na.rm=TRUE)>0)==8281,sum(rowSums(as.matrix(k[,..kcols])>0,na.rm=TRUE)>0)==1836)
# Preserve cohort aggregation; H/L numbering is NOT a patient crosswalk between assays.
groups<-data.table(Material=c('ESCC_untreated','ESCC_adjacent_untreated','ESCC_moderate','ESCC_poor'),Pattern=c('[HL][0-9]+T$','[HL][0-9]+N$','H[0-9]+T$','L[0-9]+T$'))
expand<-function(v){v<-unique(sub('-[0-9]+$','',trimws(unlist(strsplit(v,';',fixed=TRUE)))));stopifnot(all(grepl('^[A-Z0-9]+$',v)));v}
rows<-list();members<-list();selected<-list()
for(i in seq_len(nrow(groups))){
 pc<-grep(groups$Pattern[i],pcols,value=TRUE);kc<-grep(groups$Pattern[i],kcols,value=TRUE)
 stopifnot(length(pc)>0,length(kc)>0,!any(grepl('NAC',c(pc,kc))))
 pp<-expand(p[rowSums(as.matrix(p[,..pc])>0,na.rm=TRUE)>0,Protein]);kk<-expand(k[rowSums(as.matrix(k[,..kc])>0,na.rm=TRUE)>0,Proteins])
 selected[[i]]<-rbind(data.table(Material=groups$Material[i],Assay='Proteome',Column=pc),data.table(Material=groups$Material[i],Assay='Kla',Column=kc))
 for(t in c('DDR','DNA_repair')){
  a<-ann[Target==t,unique(BaseAccession)];den<-intersect(pp,a);num<-intersect(den,kk)
  rows[[length(rows)+1]]<-data.table(Material=groups$Material[i],Target=t,ProteomeSource='PXD063945',ProteinColumns=length(pc),KlaColumns=length(kc),AllProteinAccessions=length(pp),AllKlaAccessions=length(kk),Denominator=length(den),Numerator=length(num),KlaOutsideReference=length(setdiff(intersect(kk,a),pp)),ObservedPercent=100*length(num)/length(den))
  members[[length(members)+1]]<-data.table(Material=groups$Material[i],Target=t,BaseAccession=den,InKla=den%in%kk)
 }
}
fwrite(rbindlist(rows),file.path(out,'external_labels.csv'));fwrite(rbindlist(members),file.path(out,'protein_members.csv'));fwrite(rbindlist(selected),file.path(out,'selected_columns.csv'))
paths<-c(script='workflow/external_material/escc_labels.R',proteome=file.path(root,'mmc3.xlsx'),kla=file.path(root,'mmc4.xlsx'),go='outputs/20260927_ddr_fraction_labels/protein_target_membership.csv')
fwrite(data.table(Role=names(paths),Path=paths,SHA256=vapply(paths,digest::digest,character(1),algo='sha256',file=TRUE)),file.path(out,'input_sha256.csv'))
print(rbindlist(rows))
