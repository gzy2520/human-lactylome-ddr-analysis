suppressPackageStartupMessages(library(data.table))
root <- '/Users/gzy2520/Desktop/Research/kla/data'
p <- file.path(root,'PXD064038/search_results/extracted_pairing/txt/txt/La (K)Sites.txt')
x <- fread(p)
uniprot <- function(v) {
 a<-unlist(strsplit(v,'[;,]'));a<-sub('^.*\\|([^|]+)\\|.*$','\\1',a);a<-sub('-[0-9]+$','',a)
 unique(a[grepl('^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9](?:[A-Z0-9]{3}[0-9])?)$',a)])
}
flag<-function(z)tolower(as.character(z)) %in% c('+','1','true','yes','y')
keep<-!flag(x$Reverse)&!flag(x[['Potential contaminant']])&!is.na(x$id)
mem<-rbindlist(lapply(c(paste0('MEC_',1:3),paste0('NEC_',1:3)),function(s){
 v<-as.numeric(x[[paste('Localization prob',s)]]);data.table(SampleID=s,BaseAccession=uniprot(x$Proteins[keep & is.finite(v)&v>0]))
}))
a<-fread('data/publication_input/human_ddr_go_annotations.tsv')
g<-a[(is.na(get('TAXON ID'))|get('TAXON ID')==9606)&!grepl('(^|\\|)NOT($|\\|)',QUALIFIER)&get('GENE PRODUCT DB')=='UniProtKB']
ddr<-uniprot(g[['GENE PRODUCT ID']]);mem[,IsDdr:=BaseAccession %in% ddr]
mem[,`:=`(PXD='PXD064038',SampleGroup='MEC and NEC ESCC groups',Dataset='Lactylome (Kla)',SourceFile=p,UnitType='source_sample',Condition=sub('_[0-9]+$','',SampleID))]
score<-mem[,.(ProteinCount=uniqueN(BaseAccession),DdrProteinCount=uniqueN(BaseAccession[IsDdr])),by=.(PXD,SampleGroup,Dataset,SampleID,SourceFile,UnitType,Condition)]
score[,DdrFractionPercentage:=100*DdrProteinCount/ProteinCount]
old<-fread('deliverables/20261011_lactylation_ddr/inputs/proteome_samples.csv')[PXD=='PXD064038' & Dataset=='Lactylome (Kla)']
v<-merge(score,old,by=c('PXD','SampleGroup','Dataset','SampleID'));stopifnot(nrow(v)==6,all(v$ProteinCount.x==v$ProteinCount.y),all(v$DdrProteinCount.x==v$DdrProteinCount.y))
fwrite(mem,'audit/20261011_ddr_sample_recovery/base_rebuild/escc_kla_membership.csv.gz')
fwrite(score,'audit/20261011_ddr_sample_recovery/base_rebuild/escc_kla_scores.csv')
cat('All 6 ESCC Kla source samples reproduced exactly.\n')
