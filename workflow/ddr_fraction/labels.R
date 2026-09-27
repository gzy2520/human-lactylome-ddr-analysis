#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table);library(jsonlite)})
a<-commandArgs(TRUE);stopifnot(length(a)==1);out<-a[1];if(dir.exists(out))stop('Fresh output required');dir.create(out,recursive=TRUE)
paths<-c(script='workflow/ddr_fraction/labels.R',go='data/publication_input/human_ddr_go_annotations.tsv',
 ddr='audit/20260927_ddr_fraction/raw/ddr_descendants.json',repair='audit/20260927_ddr_fraction/raw/repair_descendants.json',
 kla='data/publication_input/kla_protein_membership_31.csv',reference='data/publication_input/reference_protein_membership_31.csv',
 mapping='outputs/20260920_lactate_metabolism/tables/rna_proteome_mapping_31.csv')
g<-fread(paths['go']);allann<-nrow(g)
g<-g[get('GENE PRODUCT DB')=='UniProtKB' & as.character(get('TAXON ID'))=='9606' & !grepl('(^|\\|)NOT($|\\|)',QUALIFIER)]
g[,BaseAccession:=sub('-[0-9]+$','',get('GENE PRODUCT ID'))]
sets<-list(DDR=unique(c('GO:0006974',unlist(fromJSON(paths['ddr'])$results$descendants))),DNA_repair=unique(c('GO:0006281',unlist(fromJSON(paths['repair'])$results$descendants))))
membership<-rbindlist(lapply(names(sets),function(n)unique(g[get('GO TERM')%chin%sets[[n]],.(Target=n,BaseAccession,GOterm=get('GO TERM'),Qualifier=QUALIFIER)])))
stopifnot(all(membership[Target=='DNA_repair',unique(BaseAccession)]%in%membership[Target=='DDR',unique(BaseAccession)]))
fwrite(membership,file.path(out,'protein_target_membership.csv'));fwrite(rbindlist(lapply(names(sets),function(n)data.table(Target=n,GOterm=sets[[n]]))),file.path(out,'ontology_terms.csv'))
k<-unique(fread(paths['kla'])[,.(PXD,SampleGroup,BaseAccession)])
r0<-fread(paths['reference']);r<-unique(r0[,.(BaseAccession=trimws(unlist(strsplit(MappedBaseAccessions,';',fixed=TRUE)))),by=.(PXD,SampleGroup)])
stopifnot(!anyNA(r$BaseAccession),all(nzchar(r$BaseAccession)))
map<-fread(paths['mapping']);map<-map[,.(GroupID,ReferenceKey,PXD,SampleGroup,Category,RNALabel)]
rows<-list();memberrows<-list()
for(n in names(sets))for(i in seq_len(nrow(map))){
 u<-unique(membership[Target==n,BaseAccession]);px<-map$PXD[i];sg<-map$SampleGroup[i]
 kk<-k[PXD==px & SampleGroup==sg,BaseAccession];rr<-r[PXD==px & SampleGroup==sg,BaseAccession]
 denominator<-intersect(rr,u);numerator<-intersect(denominator,kk);outs<-setdiff(intersect(kk,u),rr)
 rows[[length(rows)+1]]<-cbind(map[i],data.table(Target=n,Denominator=length(denominator),Numerator=length(numerator),KlaOutsideReference=length(outs),ObservedPercent=if(length(denominator))100*length(numerator)/length(denominator)else NA_real_))
 memberrows[[length(memberrows)+1]]<-data.table(Target=n,GroupID=map$GroupID[i],BaseAccession=denominator,InKla=denominator%in%kk)
}
z<-rbindlist(rows);stopifnot(nrow(z)==62,all(z$Denominator>0),all(is.finite(z$ObservedPercent)),all(z$Numerator<=z$Denominator))
fwrite(z,file.path(out,'protein_labels_31.csv'));fwrite(rbindlist(memberrows),file.path(out,'label_protein_members.csv'))
l<-z[,.(ObservedPercent=mean(ObservedPercent),MinimumGroupPercent=min(ObservedPercent),MaximumGroupPercent=max(ObservedPercent),ProteomeGroups=.N,GroupIDs=paste(GroupID,collapse=';')),by=.(Target,ReferenceKey)]
stopifnot(nrow(l)==56);fwrite(l,file.path(out,'material_labels_28.csv'))
fwrite(membership[,.(ProteinAccessions=uniqueN(BaseAccession)),by=Target],file.path(out,'panel_summary.csv'))
fwrite(data.table(Role=names(paths),Path=unname(paths),SHA256=vapply(paths,digest::digest,character(1),algo='sha256',file=TRUE)),file.path(out,'input_sha256.csv'))
print(z[,.(Groups=.N,MinDenominator=min(Denominator),MaxDenominator=max(Denominator),MinPercent=min(ObservedPercent),MaxPercent=max(ObservedPercent)),by=Target])
