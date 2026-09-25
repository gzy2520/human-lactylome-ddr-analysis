#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(ragg)})
args<-commandArgs(TRUE);stopifnot(length(args)==2);inp<-args[1];out<-args[2]
m<-fread(file.path(inp,'sample_metadata.csv'))
p<-fread(file.path(out,'tables/heldout_material_predictions.csv'))
s<-fread(file.path(out,'tables/heldout_sample_weak_label_predictions.csv'))
metrics<-fread(file.path(out,'tables/performance.csv'))
stopifnot(nrow(m)==1898,!anyDuplicated(m$SampleID),uniqueN(m$ReferenceKey)==28,
 uniqueN(m$ConnGroup)==14,
 all(m[DonorKnown==TRUE,.(N=uniqueN(ConnGroup)),by=DonorKey]$N==1),
 max(abs(m[,sum(MaterialBalancedWeight),by=ReferenceKey]$V1-1))<1e-12,
 max(abs(m[,sum(SourceBalancedWeight),by=ConnGroup]$V1-1))<1e-12)
s<-merge(s,m[,.(SampleID,MaterialBalancedWeight)],by='SampleID')
agg<-s[,.(Recomputed=pmin(100,pmax(0,sum(RawPrediction*MaterialBalancedWeight)))),by=.(Model,ReferenceKey)]
cmp<-merge(agg,p,by=c('Model','ReferenceKey'))
stopifnot(nrow(cmp)==112,max(abs(cmp$Recomputed-cmp$Predicted))<1e-10)
sc<-p[,.(MAE=mean(abs(Observed-Predicted))),by=.(Model,ConnGroup)]
re<-merge(sc[,.(Recomputed=mean(MAE)),by=Model],metrics,by='Model')
stopifnot(max(abs(re$Recomputed-re$SourceMAE))<1e-10)
for(root in c(inp,out)){
 manifest<-file.path(root,if(root==out)'tables' else '', 'input_sha256.csv')
 h<-fread(manifest)
 stopifnot(all(file.exists(h$Path)),all(h$SHA256==vapply(h$Path,digest::digest,character(1),algo='sha256',file=TRUE)))
}
fit<-fread(file.path(out,'tables/descriptive_full_fit.csv'))
fit<-merge(fit,unique(m[,.(ReferenceKey,ConnGroup)]),by='ReferenceKey')
fmet<-fit[,.(MAE=mean(abs(Observed-Fitted))),by=.(Model,ConnGroup)][,.(TrainingSourceMAE=mean(MAE)),by=Model]
fwrite(merge(metrics,fmet,by='Model',all.x=TRUE),file.path(out,'tables/fit_and_heldout_metrics.csv'))
labels<-c(NoRNA='No RNA baseline',MaterialMedian='Material median: 16 modules',
 SampleUnweighted='Sample: equal record weight',SampleMaterialBalanced='Sample: equal material weight',
 SampleBalanced='Sample: equal source weight',SampleBalanced128='Sample: 128 genes, source weight',
 BagMean='Bag: mean of sample features',BagDistribution='Bag: sample distribution')
order<-names(labels)
dd<-melt(metrics,id.vars='Model',measure.vars=c('MaterialMAE','SourceMAE'),variable.name='Evaluation',value.name='MAE')
dd[,Method:=factor(labels[Model],levels=rev(labels[order]))]
dd[,Evaluation:=factor(Evaluation,levels=c('MaterialMAE','SourceMAE'),labels=c('Each material has equal weight','Each source block has equal weight'))]
baseline<-dd[Model=='NoRNA']
g<-ggplot(dd,aes(MAE,Method))+geom_vline(data=baseline,aes(xintercept=MAE),linetype=2,colour='#666666')+
 geom_point(aes(colour=Model=='NoRNA'),size=3)+geom_text(aes(label=sprintf('%.2f',MAE)),hjust=-.4,size=3.3)+
 facet_wrap(~Evaluation,ncol=2)+scale_colour_manual(values=c('FALSE'='#2166AC','TRUE'='#666666'),guide='none')+
 scale_x_continuous(limits=c(0,12.5),breaks=c(0,4,8,12))+
 labs(title='Sample-level RNA: weighting helps, but does not beat the source baseline',
 subtitle='1,898 RNA records; 28 material labels; 14 donor-aware source blocks',
 x='Held-out mean absolute error (percentage points; lower is better)',y=NULL,
 caption='Nested source-held-out predictions. Dashed lines: no-RNA baseline.\nRNA records share material-level labels; this is not patient-level prediction accuracy.')+
 theme_minimal(base_size=11,base_family='Arial')+theme(panel.grid.major.y=element_blank(),panel.grid.minor=element_blank(),plot.caption=element_text(hjust=0))
ggsave(file.path(out,'figures/teacher_comparison.pdf'),g,width=12,height=5.5,device=cairo_pdf)
ggsave(file.path(out,'figures/teacher_comparison.png'),g,width=12,height=5.5,device=agg_png,dpi=180)
writeLines(c('PASS: unique sample IDs and material counts',
 'PASS: known donors do not cross outer source blocks',
 'PASS: material and source weights each sum to one',
 'PASS: held-out sample scores aggregate to saved material predictions',
 'PASS: source-level MAE independently recomputed',
 'PASS: all preparation and model input hashes match'),file.path(out,'validation.txt'))
fwrite(data.table(Path='workflow/sample_rna_fraction/validate_and_render.R',
 SHA256=digest::digest(file='workflow/sample_rna_fraction/validate_and_render.R',algo='sha256')),
 file.path(out,'tables/renderer_sha256.csv'))
for(root in c(inp,out)){
 fs<-sort(list.files(root,recursive=TRUE,full.names=TRUE))
 fs<-fs[basename(fs)!='release_sha256.csv']
 fwrite(data.table(Path=substring(fs,nchar(root)+2),SHA256=vapply(fs,digest::digest,character(1),algo='sha256',file=TRUE)),
 file.path(root,'release_sha256.csv'))
}
cat('PASS: independent validation and teacher comparison rendered.\n')
