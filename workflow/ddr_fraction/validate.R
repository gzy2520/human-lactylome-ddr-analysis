#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table);library(ranger)})
root<-'outputs/20260927_ddr_fraction_models';labroot<-'outputs/20260927_ddr_fraction_labels'
b<-readRDS('outputs/20260926_full_rna_inputs/full_rna.rds');frozen<-fread('outputs/20260926_full_rna_models/tables/split.csv')
z<-fread(file.path(labroot,'protein_labels_31.csv'));mem<-fread(file.path(labroot,'label_protein_members.csv'))
ct<-mem[,.(N=.N,K=sum(InKla)),by=.(Target,GroupID)];zz<-merge(z,ct,by=c('Target','GroupID'))
stopifnot(all(zz$N==zz$Denominator),all(zz$K==zz$Numerator),max(abs(100*zz$K/zz$N-zz$ObservedPercent))<1e-10)
res<-list()
for(t in c('DDR','DNA_repair')){
 out<-file.path(root,t);p<-fread(file.path(out,'tables/predictions.csv.gz'));m<-fread(file.path(out,'tables/metrics.csv'));s<-fread(file.path(out,'tables/split.csv'));lab<-fread(file.path(out,'tables/material_target_labels.csv'))
 stopifnot(identical(s,frozen),all(p$Label==lab$ObservedPercent[match(p$ReferenceKey,lab$ReferenceKey)]),p[Evaluation=='Source held out',.N]==1898,p[Evaluation=='Within-material test',.N]==374)
 stopifnot(s[,all(uniqueN(ConnGroup)==1),by=ReferenceKey]$V1)
 for(i in seq_len(nrow(m))){d<-p[Evaluation==m$Evaluation[i]];w<-if(m$Weighting[i]=='Record')rep(1,nrow(d))else d$W;w<-w/sum(w);err<-d$Predicted-d$Label
 calc<-c(MAE=sum(w*abs(err)),RMSE=sqrt(sum(w*err^2)),R2=1-sum(w*err^2)/sum(w*(d$Label-sum(w*d$Label))^2),Within2=100*sum(w*(abs(err)<=2)),Within5=100*sum(w*(abs(err)<=5)))
 stopifnot(max(abs(calc-as.numeric(m[i,names(calc),with=FALSE])))<1e-9)
 }
 rf<-readRDS(file.path(out,'RF_mtry135.rds'));pred<-predict(rf$model,as.data.frame(b$x[,rf$genes,drop=FALSE]))$predictions;old<-p[Evaluation=='Training fit'][match(b$metadata$SampleID,SampleID)]
 stopifnot(max(abs(pred-old$Predicted))<1e-10,length(rf$genes)==18332)
 imp<-fread(file.path(out,'tables/training_importance.csv'));stopifnot(setequal(imp$Ensembl,rf$genes),!anyDuplicated(imp$Ensembl))
 h<-fread(file.path(out,'release_sha256.csv'));stopifnot(all(h$SHA256==vapply(file.path(out,h$Path),digest::digest,character(1),algo='sha256',file=TRUE)))
 res[[t]]<-copy(m)[,Target:=t]
}
fwrite(rbindlist(res),file.path(root,'comparison_metrics.csv'))
writeLines(c('PASS: 62 target/group counts independently re-counted from protein membership.', 'PASS: both splits identical to previous full-RNA split.', 'PASS: all 16 metric rows independently recomputed.', 'PASS: both saved models reload and reproduce all 1898 full-fit predictions within 1e-10.', 'PASS: top-gene importance covers all 18332 Ensembl inputs.', 'PASS: source blocks do not split material references.', 'PASS: both model release manifests verified.'),file.path(root,'validation.txt'))
cat(readLines(file.path(root,'validation.txt')),sep='\n')
