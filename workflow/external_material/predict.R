#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table);library(ranger);library(ggplot2);library(ragg)})
a<-commandArgs(TRUE);stopifnot(length(a)==3);input<-a[1];out<-a[2];scope<-a[3];if(dir.exists(out))stop('Fresh output required');dir.create(file.path(out,'figures'),recursive=TRUE)
b<-readRDS(input);x<-b$x;s<-as.data.table(b$metadata)
if(scope=='reused_TCGA'){ii<-which(s$ReferenceKey=='TCGA_ESCA_ESCC');x<-x[ii,,drop=FALSE];s<-s[ii];s[,Material:='ESCC_untreated']}
stopifnot(all(s$Material%in%c('ESCC_untreated','ESCC_adjacent_untreated')),identical(rownames(x),s$SampleID))
labels<-fread('outputs/20260928_external_escc_labels/external_labels.csv')
rows<-list()
for(t in c('DDR','DNA_repair')){
 mod<-readRDS(file.path('outputs/20260927_ddr_fraction_models',t,'RF_mtry135.rds'));stopifnot(all(mod$genes%in%colnames(x)),all(is.finite(x)))
 pr<-predict(mod$model,as.data.frame(x[,mod$genes,drop=FALSE]))$predictions
 lab<-labels[Target==t];d<-data.table(SampleID=s$SampleID,Material=s$Material,Target=t,Scope=scope,ObservedPercent=lab$ObservedPercent[match(s$Material,lab$Material)],PredictedPercent=pr)
 d[,ErrorPP:=PredictedPercent-ObservedPercent];rows[[t]]<-d
}
p<-rbindlist(rows);stopifnot(all(is.finite(p$ObservedPercent)));fwrite(p,file.path(out,'predictions.csv'))
m<-p[,.(N=.N,Label=ObservedPercent[1],MeanPrediction=mean(PredictedPercent),MedianPrediction=median(PredictedPercent),MinimumPrediction=min(PredictedPercent),MaximumPrediction=max(PredictedPercent),MAE=mean(abs(ErrorPP)),MaterialMeanError=abs(mean(PredictedPercent)-ObservedPercent[1]),Within5=100*mean(abs(ErrorPP)<=5)),by=.(Scope,Target,Material)]
fwrite(m,file.path(out,'summary.csv'))
g<-ggplot(p,aes(Material,PredictedPercent))+geom_boxplot(outlier.shape=NA,fill='#D6EAF8',width=.5)+geom_point(position=position_jitter(width=.12,seed=25),alpha=.5,size=1.5,colour='#2166AC')+geom_point(data=m,aes(Material,Label),inherit.aes=FALSE,shape=18,size=4,colour='#B2182B')+facet_wrap(~Target)+scale_x_discrete(labels=c(ESCC_untreated='ESCC tumor',ESCC_adjacent_untreated='Adjacent esophagus'))+coord_cartesian(ylim=c(0,70))+labs(title=paste('Frozen RNA models:',scope),subtitle='Red diamond: new proteome-derived material fraction; blue dots: RNA predictions',x=NULL,y='Kla-positive / detected target proteins (%)',caption=if(scope=='reused_TCGA')'Reused training RNA; new proteome reference. Tissue-level comparison, not an independent RNA test.'else'New tissue-matched RNA; no patient pairing. RNA clinical pretreatment history is not documented.\nRed labels use untreated proteome groups. Frozen models; no fitting on external data.')+theme_minimal(base_size=11,base_family='Arial')
ggsave(file.path(out,'figures/predictions.png'),g,width=10,height=5,device=agg_png,dpi=180);ggsave(file.path(out,'figures/predictions.pdf'),g,width=10,height=5,device=cairo_pdf)
paths<-c(script='workflow/external_material/predict.R',input=input,labels='outputs/20260928_external_escc_labels/external_labels.csv',DDR='outputs/20260927_ddr_fraction_models/DDR/RF_mtry135.rds',DNA_repair='outputs/20260927_ddr_fraction_models/DNA_repair/RF_mtry135.rds')
fwrite(data.table(Role=names(paths),Path=paths,SHA256=vapply(paths,digest::digest,character(1),algo='sha256',file=TRUE)),file.path(out,'input_sha256.csv'));print(m)
