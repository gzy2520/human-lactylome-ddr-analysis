suppressPackageStartupMessages({library(data.table);library(digest)})
h<-fread('config/kla31_local_inputs_sha256.csv');print(h[!file.exists(Path),.(Path)])
e<-h[file.exists(Path)];print(e[vapply(Path,function(p)digest(file=p,algo='sha256'),'')!=SHA256,.(Path)])
f<-'corrected_final_result_20260905_sample_only';h<-fread(file.path(f,'Corrected_Data/corrected_package_sha256.csv'));cat('Frozen hash mismatches',sum(vapply(file.path(f,h$File),function(p)digest(file=p,algo='sha256'),'')!=h$SHA256),'\n')
a<-fread('outputs/20260918_qsmooth_31group_hgnc/AB_comparison.csv');cat('AB genes max delta >0.5',sum(a$MaxAbsDelta>0.5),'max',max(a$MaxAbsDelta),'genes mean delta >0.5',sum(a$MeanAbsDelta>0.5),'\n')
for(d in c('results/rna_teacher_questions','outputs/20260917_delivery_31group/rna')) {
 fs<-list.files(d,pattern='^figure1a_fixed_ddr_(31groups|28materials).csv$',full.names=TRUE)
 for(f in fs){x<-fread(f);cat(f,'\n'); print(summary(aov(DDRExpressionMedian~Category,x)));print(summary(aov(DdrFractionMedian~Category,x)))}
}
suppressPackageStartupMessages(library(data.table))
for(f in list.files('outputs/20260917_corrected_final_31group/Corrected_Data/sample_only_inputs',pattern='values.csv$',full.names=TRUE)){x<-fread(f);cat('PROTEOME SAMPLE',basename(f),nrow(x),'all_sample',all(x$ObservationType=='sample'),'\n')}
for(key in c('GSE132714_BPH','GSE163787_TALL104'))cat(key,'local_matrix_md5',unname(tools::md5sum(file.path('outputs/20260916_expression_extraction/matrices',paste0(key,'_log2tpm.tsv.gz')))),'\n')
d<-'outputs/20260917_delivery_31group/rna'; a<-fread(file.path(d,'figure1a_fixed_ddr_28materials.csv'));b<-fread(file.path(d,'figure1a_lactylated_gene_28materials.csv'));cat('DUAL columns',names(b),'\n')
x<-merge(a[,.(GroupID,Category,DDRExpressionMedian)],b[,.(GroupID,KlaExpressionMedian)],by='GroupID');z<-melt(x,id.vars=c('GroupID','Category'),variable.name='GeneSet',value.name='y');cat('Within material score correlation',cor(x$DDRExpressionMedian,x$KlaExpressionMedian),'\n');print(summary(aov(y~Category*GeneSet,z)));cat('Paired sensitivity only, not replacement official test\n');print(summary(aov(y~Category*GeneSet+Error(GroupID/GeneSet),z)))
