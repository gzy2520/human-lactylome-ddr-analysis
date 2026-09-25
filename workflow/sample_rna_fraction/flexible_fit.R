#!/usr/bin/env Rscript
# Descriptive nonlinear bag regression; capacity selected to meet a TRAINING R2 target.
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(ragg)})
set.seed(25)
args<-commandArgs(TRUE);stopifnot(length(args)==2);inp<-args[1];out<-args[2]
if(dir.exists(out))stop('Use a new output directory')
dir.create(file.path(out,'tables'),recursive=TRUE);dir.create(file.path(out,'figures'))
paths<-c(script='workflow/sample_rna_fraction/flexible_fit.R',
 rna=file.path(inp,'sample_rna128.csv.gz'),metadata=file.path(inp,'sample_metadata.csv'),labels=file.path(inp,'material_labels_28.csv'))
d<-fread(paths['labels']);setorder(d,ReferenceKey)
s<-merge(fread(paths['rna']),fread(paths['metadata'])[,.(SampleID,MaterialBalancedWeight,SourceBalancedWeight)],by='SampleID')
genes<-sort(grep('^ENSG',names(s),value=TRUE));x<-as.matrix(s[,..genes])
bi<-match(s$ReferenceKey,d$ReferenceKey);y<-d$ObservedPercent;w<-s$MaterialBalancedWeight
features<-function(train){
 ii<-which(bi %in% train);ww<-s$SourceBalancedWeight[ii];ww<-ww/sum(ww)
 mu<-as.vector(crossprod(ww,x[ii,,drop=FALSE]));z<-sweep(x,2,mu,'-')
 sd<-sqrt(as.vector(crossprod(ww,z[ii,,drop=FALSE]^2)));sd[sd<1e-8]<-1
 z<-sweep(z,2,sd,'/');a<-rowsum(z*w,bi);b<-sqrt(pmax(rowsum(z*z*w,bi)-a*a,0))
 stopifnot(identical(dim(a),dim(b)),nrow(a)==28,ncol(a)==128)
 cbind(a,b) # Every sample contributes; singleton bag SD is zero.
}
dist2<-function(a,b)pmax(outer(rowSums(a*a),rowSums(b*b),'+')-2*tcrossprod(a,b),0)
grid<-10^seq(2,-6,by=-1)
fit<-function(train){
 f<-features(train);dd<-dist2(f[train,,drop=FALSE],f[train,,drop=FALSE]);bw<-median(dd[upper.tri(dd)])
 stopifnot(is.finite(bw),bw>0)
 k<-exp(-dd/(2*bw));ka<-exp(-dist2(f,f[train,,drop=FALSE])/(2*bw));mu<-mean(y[train])
 curve<-rbindlist(lapply(grid,function(lam){
   pred<-mu+as.vector(k%*%solve(k+diag(lam,length(train)),y[train]-mu))
   data.table(Lambda=lam,TrainingR2=1-sum((pred-y[train])^2)/sum((y[train]-mu)^2),TrainingMAE=mean(abs(pred-y[train])))
 }))
 # Deliberately optimize descriptive fit, not validation performance. No test label enters selection.
 j<-which(curve$TrainingR2>=.90)[1];if(is.na(j))j<-nrow(curve)
 lam<-curve$Lambda[j];beta<-solve(k+diag(lam,length(train)),y[train]-mu)
 list(pred=pmin(100,pmax(0,mu+as.vector(ka%*%beta))),curve=curve,lambda=lam,
      saved=list(TrainingIndices=train,FeatureMatrix=f,Bandwidth=bw,Beta=beta,Intercept=mu,Lambda=lam))
}
full<-fit(seq_len(nrow(d)));pred<-rep(NA_real_,nrow(d));folds<-list()
for(b in sort(unique(d$ConnGroup))){
 tr<-which(d$ConnGroup!=b);te<-which(d$ConnGroup==b);v<-fit(tr);pred[te]<-v$pred[te]
 folds[[as.character(b)]]<-data.table(ConnGroup=b,SelectedLambda=v$lambda)
}
p<-d[,.(ReferenceKey,ConnGroup,Observed=ObservedPercent)]
p[,`:=`(TrainingFitted=full$pred,HeldoutPrediction=pred)]
fwrite(p,file.path(out,'tables/predictions.csv'));fwrite(full$curve,file.path(out,'tables/training_capacity_curve.csv'))
fwrite(rbindlist(folds),file.path(out,'tables/fold_penalties.csv'));saveRDS(full$saved,file.path(out,'descriptive_model.rds'))
long<-melt(p,id.vars=c('ReferenceKey','ConnGroup','Observed'),variable.name='Evaluation',value.name='Predicted')
metrics<-long[,.(MAE=mean(abs(Predicted-Observed)),RMSE=sqrt(mean((Predicted-Observed)^2)),
 R2=1-sum((Predicted-Observed)^2)/sum((Observed-mean(Observed))^2)),by=Evaluation]
sm<-long[,.(E=mean(abs(Predicted-Observed))),by=.(Evaluation,ConnGroup)][,.(SourceMAE=mean(E)),by=Evaluation]
metrics<-merge(metrics,sm,by='Evaluation');fwrite(metrics,file.path(out,'tables/metrics.csv'))
lab<-c(TrainingFitted='Training fit: existing materials',HeldoutPrediction='Source held out: unseen materials')
long<-merge(long,metrics[,.(Evaluation,R2,MAE)],by='Evaluation')
long[,Panel:=factor(sprintf('%s\nR2 = %.3f; MAE = %.2f pp',lab[Evaluation],R2,MAE),
 levels=sprintf('%s\nR2 = %.3f; MAE = %.2f pp',lab[metrics$Evaluation],metrics$R2,metrics$MAE))]
g<-ggplot(long,aes(Observed,Predicted))+geom_abline(slope=1,intercept=0,colour='grey60',linetype=2)+
 geom_point(colour='#2166AC',size=2.5)+facet_wrap(~Panel,nrow=1)+coord_equal(xlim=c(0,40),ylim=c(0,40))+
 labs(title='A flexible RNA model can fit the existing group labels',
 subtitle='128 genes; sample mean + variability per material; nonlinear kernel regression',
 x='Observed group detection percentage',y='Model output (%)',
 caption='Capacity chosen to reach training R2 >= 0.90; this is an explicitly descriptive objective.\n1,898 RNA records form 28 labeled bags. Right: same rule refitted within each of 14 source-held-out folds.')+
 theme_minimal(base_size=11,base_family='Arial')+theme(panel.grid.minor=element_blank(),plot.caption=element_text(hjust=0))
ggsave(file.path(out,'figures/flexible_fit_and_heldout.pdf'),g,width=10,height=6,device=cairo_pdf)
ggsave(file.path(out,'figures/flexible_fit_and_heldout.png'),g,width=10,height=6,device=agg_png,dpi=180)
stopifnot(all(is.finite(pred)),nrow(p)==28,full$curve[Lambda==full$lambda,TrainingR2]>=.90)
fwrite(data.table(Role=names(paths),Path=unname(paths),SHA256=vapply(paths,digest::digest,character(1),algo='sha256',file=TRUE)),file.path(out,'input_sha256.csv'))
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
fs<-sort(list.files(out,recursive=TRUE,full.names=TRUE))
fwrite(data.table(Path=substring(fs,nchar(out)+2),SHA256=vapply(fs,digest::digest,character(1),algo='sha256',file=TRUE)),file.path(out,'release_sha256.csv'))
print(metrics);cat('Selected descriptive penalty:',full$lambda,'\n')
