suppressPackageStartupMessages(library(data.table))
r<-'deliverables/20261011_lactylation_ddr_tall_verified'
out<-'audit/20261011_all_groups_algorithm_audit/cross_group'
u<-fread(file.path(r,'tables/ddr_source_unit_values.csv'))
s<-fread(file.path(r,'tables/ddr_mean_sem_62.csv'))
g<-fread('data/publication_input/group_summary_31.csv')
a<-fread(cmd=paste('gzip -dc',file.path(r,'inputs/ddr_protein_membership.csv.gz')))
f<-rbindlist(list(g[,.(PXD,SampleGroup,Dataset='Lactylome (Kla)',HistoricalN=KlaProteinCount,HistoricalDdr=KlaDdrProteinCount,HistoricalFraction=100*KlaDdrFraction)],g[,.(PXD,SampleGroup,Dataset='Whole proteome',HistoricalN=ReferenceProteinCount,HistoricalDdr=ReferenceDdrProteinCount,HistoricalFraction=100*ReferenceDdrFraction)]))
unions<-a[,.(CurrentUnionN=uniqueN(BaseAccession),CurrentUnionDdr=uniqueN(BaseAccession[IsDdr])),by=.(PXD,SampleGroup,Dataset)]
unions[,CurrentUnionFraction:=100*CurrentUnionDdr/CurrentUnionN]
z<-merge(merge(s,f,by=c('PXD','SampleGroup','Dataset')),unions,by=c('PXD','SampleGroup','Dataset'))
z[,`:=`(MeanVsHistorical=Mean-HistoricalFraction,CurrentUnionVsHistorical=CurrentUnionFraction-HistoricalFraction,CurrentMeanVsUnion=Mean-CurrentUnionFraction)]
fwrite(z,file.path(out,'historical_union_vs_current_mean_62.csv'))
reuse<-u[,.(DisplayPositions=uniqueN(GroupID),Groups=paste(sort(unique(GroupID)),collapse=';'),FractionVariants=uniqueN(Fraction)),by=.(Dataset,SourceFile,SampleID)][DisplayPositions>1]
fwrite(reuse,file.path(out,'source_unit_reuse.csv'))
conditions<-u[,.(N=.N,Mean=mean(Fraction),UnitTypes=paste(sort(unique(UnitType)),collapse=';')),by=.(GroupID,PXD,SampleGroup,Dataset,Condition)]
fwrite(conditions,file.path(out,'conditions_and_units.csv'))
# Identity of source memberships; repeated complete sets highlight imputation or shared profile tables for source inspection.
sets<-a[,.(ProteinSignature=paste(sort(BaseAccession),collapse=';')),by=.(PXD,SampleGroup,Dataset,SampleID)]
identicalsets<-sets[,.(N=.N,SampleIDs=paste(sort(SampleID),collapse=';')),by=.(PXD,SampleGroup,Dataset,ProteinSignature)][N>1]
identicalsets[,ProteinSignature:=NULL]
fwrite(identicalsets,file.path(out,'identical_source_membership_sets.csv'))
stopifnot(nrow(z)==62,all(z$CurrentUnionDdr<=z$CurrentUnionN))
print(reuse);print(identicalsets)

unit_sets<-merge(u[,.(PXD,SampleGroup,Dataset,SampleID,SourceFile,GroupID)],sets,by=c('PXD','SampleGroup','Dataset','SampleID'))
sharing<-unit_sets[,.(DisplayPositions=uniqueN(GroupID),Groups=paste(sort(unique(GroupID)),collapse=';'),DistinctProteinSets=uniqueN(ProteinSignature)),by=.(Dataset,SourceFile,SampleID)][DisplayPositions>1]
sharing[,Finding:=ifelse(DistinctProteinSets==1,'same source profile reused across display groups','sample label reused for different material-specific source columns')]
fwrite(sharing,file.path(out,'source_identity_review.csv'))
