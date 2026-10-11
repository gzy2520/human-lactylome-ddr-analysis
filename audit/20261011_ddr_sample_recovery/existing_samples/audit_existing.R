suppressPackageStartupMessages({library(data.table);library(readxl)})
root <- '/Users/gzy2520/Desktop/Research/kla/.worktrees/kla31-rnaseq-20260912'
rawroot <- '/Users/gzy2520/Desktop/Research/kla/data'
out <- file.path(root,'audit/20261011_ddr_sample_recovery/existing_samples')
base <- function(x) {x<-trimws(x);x<-sub('^(sp|tr)\\|','',x);x<-sub('\\|.*$','',x);x<-sub('^NX_','',x);sub('-[0-9]+$','',x)}
ids <- function(x) {x<-base(unlist(strsplit(as.character(x),'[;,]'))); sort(unique(x[!is.na(x)&grepl('^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9](?:[A-Z0-9]{3}[0-9])?)$',x)]))}
pos <- function(x) {x<-suppressWarnings(as.numeric(x));is.finite(x)&x>0}
g <- fread(file.path(root,'data/publication_input/human_ddr_go_annotations.tsv'))
ddr <- ids(g[(is.na(`TAXON ID`)|`TAXON ID`==9606)&!grepl('(^|\\|)NOT($|\\|)',QUALIFIER)&`GENE PRODUCT DB`=='UniProtKB',`GENE PRODUCT ID`])
frozen <- fread(file.path(root,'data/publication_input/kla_protein_membership_31.csv'))
old <- fread(file.path(root,'deliverables/20261011_lactylation_ddr/inputs/proteome_samples.csv'))
paths <- unique(old[,.(SourceFile)])
paths[, ResolvedPath := vapply(SourceFile,function(p) {candidates<-c(p,file.path(rawroot,p),file.path(dirname(rawroot),p),file.path(root,p)); hit<-candidates[file.exists(candidates)];if(length(hit))normalizePath(hit[1])else ''},character(1))]
paths[,Exists:=nzchar(ResolvedPath)]
fwrite(paths,file.path(out,'source_path_audit.csv'))
all <- list()
add <- function(case,sample,acc,method,src) {all[[length(all)+1L]]<<-data.table(Case=case,SampleID=sample,BaseAccession=acc,IsDdr=acc%in%ddr,Method=method,SourceFile=src)}
tall <- file.path(rawroot,'PXD028488/search_results/Enrichment-Search files/TALL-NALAC-Search files/DB search psm.csv')
t <- fread(tall)
fwrite(t[,.(N=.N,PTM=paste(sort(unique(PTM)),collapse=';')),by='Source File'],file.path(out,'tall_source_files.csv'))
for(s in sort(unique(t$`Source File`))) {k<-t$`Source File`==s & grepl('(^|; *)Lac($|;)',t$PTM) & grepl('K\\(\\+72\\.02\\)',t$Peptide);add('TALL_Kla_corrected',s,ids(t$Accession[k]),'K(+72.02) peptide and Lac PTM; union over PSMs',tall)}
tp <- fread(file.path(dirname(tall),'proteins.csv'))
for(s in 1:3) {acc<-ids(tp$Accession[pos(tp[[paste('Area Sample',s)]])]);add('TALL_enrichment_all_proteins',paste('Sample',s),acc,'positive area without Kla restriction',file.path(dirname(tall),'proteins.csv'))}
pcpath<-file.path(rawroot,'PXD063266/search_results/LactylSites.xlsx');p<-as.data.table(read_excel(pcpath,sheet='Lactyl (K)Sites'))
for(s in 1:3) {k<-!is.na(p$id)&pos(p[[paste('Localization prob',s)]]);add('PC3M_Kla_localization',as.character(s),ids(p$Proteins[k]),'localized site probability >0',pcpath);add('PC3M_Kla_intensity',as.character(s),ids(p$Proteins[k&pos(p[[paste('Intensity',s)]])]),'localized site probability >0 and intensity >0',pcpath)}
espath<-file.path(rawroot,'PXD065830/supplementary/Dataset1.xlsx');e<-read_excel(espath,sheet='2.a protein raw information',col_names=FALSE,.name_repair='minimal');h<-trimws(as.character(unlist(e[2,])));ix<-which(grepl('^ESCC-[0-9]+T$',h));for(j in ix) add('ESCC_proteome',h[j],ids(e[[1]][-c(1,2)][pos(e[[j]][-c(1,2)])]),'ESCC T sample raw intensity >0',espath)
m<-unique(rbindlist(all));fwrite(m,file.path(out,'recomputed_membership.csv.gz'))
sum<-m[,.(ProteinCount=.N,DdrProteinCount=sum(IsDdr),FractionPercentage=100*mean(IsDdr)),by=.(Case,SampleID,Method,SourceFile)];fwrite(sum,file.path(out,'recomputed_counts.csv'))
agg<-sum[,.(N=.N,Mean=mean(FractionPercentage),SD=sd(FractionPercentage)),by=Case];print(agg);fwrite(agg,file.path(out,'recomputed_summary.csv'))
# Frozen join reproduces historical undercount where source contains new proteins.
for(case in c('PC3M_Kla_localization','TALL_enrichment_all_proteins')) {pxd<-if(case=='PC3M_Kla_localization')'PXD063266' else 'PXD028488';group<-if(pxd=='PXD063266')'PC-3M' else 'TALL-104';f<-frozen[PXD==pxd&SampleGroup==group&IsDdr==TRUE,BaseAccession];print(m[Case==case,.(N=.N,GlobalDDR=sum(IsDdr),FrozenDDR=sum(BaseAccession%in%f),NotInFrozen=sum(!BaseAccession%in%frozen[PXD==pxd&SampleGroup==group,BaseAccession])),by=SampleID])}
comp<-merge(sum[Case=='ESCC_proteome'],old[Dataset=='Whole proteome'&RowOrder==10,.(SampleID,OldCount=ProteinCount,OldDdr=DdrProteinCount)],by='SampleID');print(comp[,.(N=.N,CountMismatch=sum(ProteinCount!=OldCount),DdrMismatch=sum(DdrProteinCount!=OldDdr))]);fwrite(comp,file.path(out,'escc_old_vs_recomputed.csv'))
print(paths[,.(Paths=.N,Missing=sum(!Exists))])
pppath<-file.path(dirname(tall),'protein-peptides.csv'); pp<-fread(pppath)
for(s in 1:3) {k<-grepl('K\\(\\+72\\.02\\)',pp$Peptide)&grepl('(^|; *)Lac($|;)',pp$PTM)&pos(pp[[paste('Area Sample',s)]]);add('TALL_Kla_peptide_area',paste('Sample',s),ids(pp$`Protein Accession`[k]),'K(+72.02) peptide; Lac PTM; sample peptide Area >0',pppath)}
m<-unique(rbindlist(all));fwrite(m,file.path(out,'recomputed_membership.csv.gz'));sum<-m[,.(ProteinCount=.N,DdrProteinCount=sum(IsDdr),FractionPercentage=100*mean(IsDdr)),by=.(Case,SampleID,Method,SourceFile)];fwrite(sum,file.path(out,'recomputed_counts.csv'));fwrite(sum[,.(N=.N,Mean=mean(FractionPercentage),SD=sd(FractionPercentage)),by=Case],file.path(out,'recomputed_summary.csv'));print(sum[Case=='TALL_Kla_peptide_area']);print(m[grepl('TALL_Kla',Case),.(Union=.N,Ddr=sum(IsDdr)),by=.(Case,BaseAccession)][,.(Union=.N,Ddr=sum(Ddr>0)),by=Case])
for(s in 1:3) print(p[!is.na(id)&pos(get(paste('Localization prob',s)))&!pos(get(paste('Intensity',s))),.(Sites=.N,FiniteLoc=sum(pos(get(paste('Localization prob',s)))),MedianScore=median(get(paste('Score',s)),na.rm=TRUE))])
strippep<-function(x)gsub('\\([^)]*\\)|\\[[^]]*\\]','',sub('^[A-Z]\\.([A-Z].*)\\.[A-Z]$','\\1',x))
map<-unique(data.table(PeptideClean=strippep(pp$Peptide),BaseAccession=base(pp$`Protein Accession`),Start=pp$Start))[!is.na(Start)]
for(s in sort(unique(t$`Source File`))) {
 rows<-t[`Source File`==s & `-10lgP`>=20 & !grepl('DECOY|REV|REVERSED|CON_|CONTAM|CRAP',Accession,ignore.case=TRUE)]
 acc<-lapply(seq_len(nrow(rows)),function(i) {pep<-strippep(rows$Peptide[i]); if(!grepl('K[0-9]+:Lac',rows$AScore[i])&&!grepl('K\\(\\+72\\.02\\)',rows$Peptide[i]))return(character());intersect(map[PeptideClean==pep,BaseAccession],ids(rows$Accession[i]))})
 add('TALL_Kla_legacy_mapping',s,sort(unique(unlist(acc))),'PSM -10lgP>=20; Lac; matching peptide-protein mapping; frozen legacy parser criteria',tall)
 astr<-regmatches(rows$AScore,gregexpr('K[0-9]+:Lac:[0-9.]+',rows$AScore)); pass<-vapply(astr,function(x) any(as.numeric(sub('.*:','',x))>=20),logical(1))
 add('TALL_Kla_AScore20',s,sort(unique(unlist(acc[pass]))),'legacy mapping plus at least one Lac AScore>=20',tall)
}
m<-unique(rbindlist(all));fwrite(m,file.path(out,'recomputed_membership.csv.gz'));sum<-m[,.(ProteinCount=.N,DdrProteinCount=sum(IsDdr),FractionPercentage=100*mean(IsDdr)),by=.(Case,SampleID,Method,SourceFile)];fwrite(sum,file.path(out,'recomputed_counts.csv'));fwrite(sum[,.(N=.N,Mean=mean(FractionPercentage),SD=sd(FractionPercentage)),by=Case],file.path(out,'recomputed_summary.csv'));print(sum[grepl('TALL_Kla_(legacy|AScore)',Case),.(Case,SampleID,ProteinCount,DdrProteinCount,FractionPercentage)])
