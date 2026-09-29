#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
out<-'outputs/20260929_assay_composition/labels';stopifnot(!dir.exists(out));dir.create(out,recursive=TRUE)
ann<-fread('outputs/20260927_ddr_fraction_labels/protein_target_membership.csv')
r<-unique(fread('data/publication_input/reference_protein_membership_31.csv'))
k<-unique(fread('data/publication_input/kla_protein_membership_31.csv'))
map<-fread('outputs/20260920_lactate_metabolism/tables/rna_proteome_mapping_31.csv')[,.(GroupID,ReferenceKey,PXD,SampleGroup,Category,RNALabel)]
g<-fread('data/publication_input/group_summary_31.csv')
stopifnot(nrow(map)==31,!anyDuplicated(r,by=c('PXD','SampleGroup','SourceProteinID')),!anyDuplicated(k,by=c('PXD','SampleGroup','BaseAccession')))
# Reproduce Figure 1a counting units. Preserve Ensembl protein IDs and their mappings.
rows<-list();members<-list();sensitivity<-list()
for(pathway in c('DDR','DNA_repair')){
 a<-ann[Target==pathway,unique(BaseAccession)]
 for(assay in c('Proteome','Kla'))for(i in seq_len(nrow(map))){
  px<-map$PXD[i];sg<-map$SampleGroup[i];target<-paste(assay,pathway,sep='_')
  if(assay=='Proteome'){
   d<-r[PXD==px & SampleGroup==sg];ids<-d$SourceProteinID;maps<-strsplit(d$MappedBaseAccessions,';',fixed=TRUE)
   hit<-vapply(maps,function(z)any(z%in%a),logical(1));mapped<-d$MappedBaseAccessions
   if(pathway=='DDR')stopifnot(identical(hit,d$IsDdr))
   u<-unique(unlist(maps));sensitivity[[length(sensitivity)+1]]<-data.table(GroupID=map$GroupID[i],Target=target,SourceIDs=length(ids),UniqueMappedAccessions=length(u),MultiMappedSourceIDs=sum(lengths(maps)>1),SourceIDPercent=100*mean(hit),AccessionPercent=100*mean(u%in%a))
  }else{d<-k[PXD==px & SampleGroup==sg];ids<-d$BaseAccession;mapped<-ids;hit<-ids%in%a;if(pathway=='DDR')stopifnot(identical(hit,d$IsDdr))}
  stopifnot(length(ids)>0,!anyNA(hit))
  rows[[length(rows)+1]]<-cbind(map[i],data.table(Target=target,Assay=assay,Pathway=pathway,Numerator=sum(hit),Denominator=length(ids),ObservedPercent=100*mean(hit)))
  members[[length(members)+1]]<-data.table(GroupID=map$GroupID[i],Target=target,ProteinID=ids,MappedBaseAccessions=mapped,InPathway=hit)
  if(pathway=='DDR'){
   old<-g[PXD==px & SampleGroup==sg];stopifnot(nrow(old)==1)
   stopifnot(length(ids)==if(assay=='Proteome')old$ReferenceProteinCount else old$KlaProteinCount,
    sum(hit)==if(assay=='Proteome')old$ReferenceDdrProteinCount else old$KlaDdrProteinCount)
  }
 }
}
z<-rbindlist(rows);stopifnot(nrow(z)==124,all(z$Numerator<=z$Denominator))
fwrite(z,file.path(out,'protein_labels_31.csv'));fwrite(rbindlist(members),file.path(out,'protein_members.csv.gz'));fwrite(rbindlist(sensitivity),file.path(out,'identifier_sensitivity.csv'))
l<-z[,.(ObservedPercent=mean(ObservedPercent),MinimumGroupPercent=min(ObservedPercent),MaximumGroupPercent=max(ObservedPercent),ProteomeGroups=.N,GroupIDs=paste(GroupID,collapse=';')),by=.(Target,Assay,Pathway,ReferenceKey)]
stopifnot(nrow(l)==112);fwrite(l,file.path(out,'material_labels_28.csv'))
paths<-c('workflow/assay_composition/labels.R','outputs/20260927_ddr_fraction_labels/protein_target_membership.csv','data/publication_input/reference_protein_membership_31.csv','data/publication_input/kla_protein_membership_31.csv','data/publication_input/group_summary_31.csv','outputs/20260920_lactate_metabolism/tables/rna_proteome_mapping_31.csv')
fwrite(data.table(Path=paths,SHA256=vapply(paths,digest::digest,character(1),file=TRUE,algo='sha256')),file.path(out,'input_sha256.csv'))
writeLines('PASS: all 31 DDR numerators and denominators match frozen Figure 1a group summary, in each assay. Four targets; no cross-assay intersection. Stable source protein IDs retained for proteome; BaseAccessions for Kla.',file.path(out,'validation.txt'))
print(z[,.(Groups=.N,Min=min(ObservedPercent),Max=max(ObservedPercent)),by=Target])
