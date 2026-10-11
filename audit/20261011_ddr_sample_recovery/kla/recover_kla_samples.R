suppressPackageStartupMessages({library(data.table);library(readxl)})
set.seed(25)
project_root <- normalizePath('.')
source_root <- '/Users/gzy2520/Desktop/Research/kla/data'
out <- file.path(project_root,'audit/20261011_ddr_sample_recovery/kla')
# Reuse read-only parsers; stop before their shared-output writing section.
code <- readLines('R/candidate/prepare_sample_boxplot_inputs.R')
eval(parse(text=code[30:632]))
res <- list(); ledger <- list()
add <- function(d,unit,condition=NULL,evidence,sem){
 d[,UnitType:=unit]; d[,Condition:=if(is.null(condition)) SampleClass else condition]
 d[,ObservationType:=unit]; d[,SemEligible:=sem]; d[,UnitEvidence:=evidence]
 res[[length(res)+1L]] <<- d
}
p <- file.path(source_root,'PXD075377/search_results/extracted/2-Basic_analysis/MS_identified_information.txt')
for(g in c('adjacent liver','HCC')){
 cond<-if(g=='HCC') 'HCC' else 'Control'
 d<-read_itraq_pool(p,'PXD075377',g,paste0('patient1_',cond),cond,paste('Intensity',cond))
 add(d,'tissue_sample',cond,'dataset_metadata.csv: paired tissues from one clinical patient; each tissue channel is one sample',FALSE)
}
p <- file.path(source_root,'PXD058534/search_results/extracted_pairing/txt/La (K)Sites.txt')
d<-read_maxquant_single_group(p,'PXD058534','pretreated HK-2','HK2','pretreated HK-2')
add(d,'single_sample',NULL,'summary.txt: one raw file FE458LPLa_HK2_Slot1-35_1_14667, Experiment HK2',FALSE)
p<-file.path(source_root,'PXD066351/search_results/XB01472B1DPLa-MSstats_Input.csv')
map<-sample_map('PXD066351','HCT116 control and Roseburia co-culture',c('NC116','R116'),c('control','Roseburia co-culture'))
d<-read_pxd066351_by_sample(p,'PXD066351',map$SampleGroup[1],map)
add(d,'single_sample_per_condition',NULL,'R.FileName gives NC116 and R116 one raw file each; R.Replicate=0; R.Fraction=NA',FALSE)
p<-file.path(source_root,'PXD070007/search_results/SA206LPLaB1_Annotation.xlsx')
for(g in c('glioblastoma stem cells','neural stem cells')){
 ids<-if(g=='neural stem cells')c('ENSA','HMP1')else c('G2907','G3028','G3264','GSC23','MES28','RKI')
 map<-sample_map('PXD070007',g,ids,'model')
 d<-read_pxd070007_by_model(p,'PXD070007',g,map)
 add(d,'model',g,'Annotation_Combine contains named model columns; no within-model replicate columns or raw-file identifiers',FALSE)
}
p<-file.path(source_root,'PXD073311/search_results/extracted_pairing/IPX0015307003_Database_search_result/Database_search_result/report.pr_matrix.tsv')
map<-sample_map('PXD073311','HUVEC control and Pg infection',c(paste0('A0h_',1:3),paste0('A6h_',1:3)),rep(c('A0h control','A6h Pg infection'),each=3))
d<-read_pxd073311_by_raw_file(p,'PXD073311',map$SampleGroup[1],map)
add(d,'sample',NULL,'pr_matrix: six distinct sample raw-file columns A0h_1/2/3 and A6h_1/2/3',TRUE)
# Existing evidence-to-protein join preserves the historical protein-group representative rule.
eval(parse(text=code[752:776]))
protein_sites[,SampleClass:=ifelse(grepl('sg737',SampleID),'GSK3B KO','GSK3B WT')]
add(protein_sites,'sample',NULL,'Experiment sg_1/sg_2 and sg737_1/sg737_2; metadata states two replicates each; SDRF biological replicate blank',TRUE)
for(px in c('PXD028488','PXD053474')){
 d<-fread(file.path(out,paste0(tolower(px),'_source_sites.csv')))
 for(g in unique(d$CellOrTissueType)){
  x<-d[CellOrTissueType==g]
  unit<-if(px=='PXD028488'&&g=='HCT116')'single_sample' else 'technical_sourceunit'
  desc<-if(px=='PXD053474')'S3 plus site-level DDA/DIA corroboration; all compartments, fractions and acquisition modes merged'else if(g=='HEK293T')'HCD24/28/32/36/40 collision-energy runs of one mixed material merged'else 'one raw file 20211127-HCT116-LAC-Enrichment-1.raw'
  z<-records(px,g,paste0(g,if(unit=='technical_sourceunit')'_technical_merged' else '_1'),sort(unique(x$BaseAccession)),'reparsed_primary_source',paste(sort(unique(x$SourceFile)),collapse=';'),g,unit)
  add(z,unit,g,desc,FALSE)
 }
}
a<-unique(rbindlist(res,fill=TRUE),by=c('PXD','SampleGroup','SampleID','BaseAccession'))
go<-fread('data/publication_input/human_ddr_go_annotations.tsv',quote='',check.names=FALSE)
keep<-(is.na(go[['TAXON ID']])|go[['TAXON ID']]==9606)&!grepl('(^|\\|)NOT($|\\|)',go$QUALIFIER)&go[['GENE PRODUCT DB']]=='UniProtKB'
ddr<-sort(unique(base_accession(go[['GENE PRODUCT ID']][keep])))
a[,IsDdr:=BaseAccession %in% ddr]
a[,SourceFile:=gsub(paste0(source_root,'/'),'',SourceFile,fixed=TRUE)]
frozen<-fread('data/publication_input/kla_protein_membership_31.csv')
a<-merge(a,frozen[,.(PXD,SampleGroup,BaseAccession,FrozenIsDdr=IsDdr)],by=c('PXD','SampleGroup','BaseAccession'),all.x=TRUE,sort=FALSE)
a[,FrozenMembershipPresent:=!is.na(FrozenIsDdr)]
fwrite(a,file.path(out,'recovered_kla_membership.csv.gz'))
fwrite(a[!FrozenMembershipPresent|IsDdr!=FrozenIsDdr],file.path(out,'membership_annotation_differences.csv'))
s<-a[,.(KlaProteinCount=.N,KlaDdrProteinCount=sum(IsDdr),KlaDdrFraction=mean(IsDdr),KlaDdrFractionPercentage=100*mean(IsDdr),SourceFile=SourceFile[1],UnitEvidence=UnitEvidence[1]),by=.(PXD,SampleGroup,SampleID,Condition,UnitType,SemEligible)]
fwrite(s,file.path(out,'recovered_kla_sample_values.csv'))
z<-s[,.(N=.N,Mean=mean(KlaDdrFractionPercentage),SD=if(.N>1)sd(KlaDdrFractionPercentage)else NA_real_,SEM=if(.N>1&&all(SemEligible))sd(KlaDdrFractionPercentage)/sqrt(.N)else NA_real_),by=.(PXD,SampleGroup,Condition,UnitType)]
fwrite(z,file.path(out,'condition_summary.csv'))
print(s[,.(PXD,SampleGroup,SampleID,UnitType,KlaProteinCount,KlaDdrProteinCount,KlaDdrFractionPercentage)])
cat('Annotation differences among frozen members:',nrow(a[FrozenMembershipPresent&IsDdr!=FrozenIsDdr]),'\nOutside frozen membership:',nrow(a[FrozenMembershipPresent==FALSE]),'\n')
stopifnot(!anyDuplicated(a[,.(PXD,SampleGroup,SampleID,BaseAccession)]),all(s$KlaDdrProteinCount<=s$KlaProteinCount),all(s$KlaProteinCount>0))
