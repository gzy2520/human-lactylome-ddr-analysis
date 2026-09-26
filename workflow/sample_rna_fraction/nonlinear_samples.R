#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table);library(ranger);library(ggplot2);library(ragg)})
set.seed(25)
a<-commandArgs(TRUE);stopifnot(length(a)==2);inp<-a[1];out<-a[2];if(dir.exists(out))stop('Fresh output required')
dir.create(file.path(out,'tables'),recursive=TRUE);dir.create(file.path(out,'figures'))
paths<-c(script='workflow/sample_rna_fraction/nonlinear_samples.R',rna=file.path(inp,'sample_rna128.csv.gz'),metadata=file.path(inp,'sample_metadata.csv'))
s<-merge(fread(paths['rna']),fread(paths['metadata'])[,.(SampleID,DonorKey,DonorKnown,ConnGroup,ObservedPercent,MaterialBalancedWeight)],by='SampleID');setorder(s,ReferenceKey,SampleID)
genes<-sort(grep('^ENSG',names(s),value=TRUE));x<-as.data.frame(s[,..genes]);y<-s$ObservedPercent
stopifnot(nrow(s)==1898,length(genes)==128,!anyDuplicated(s$SampleID))
# One fixed donor split; all known records of a donor stay together, even across materials.
donors<-sort(unique(s$DonorKey));testdonors<-sample(donors,floor(length(donors)*.2))
# Never leave a material without a training record. Singleton materials are training-only.
for(k in sort(unique(s$ReferenceKey))){v<-unique(s[ReferenceKey==k,DonorKey]);if(all(v%in%testdonors))testdonors<-setdiff(testdonors,v[1])}
te<-which(s$DonorKey%in%testdonors);tr<-setdiff(seq_len(nrow(s)),te)
stopifnot(!length(intersect(s$DonorKey[tr],s$DonorKey[te])),all(s$ReferenceKey[te]%in%s$ReferenceKey[tr]))
s[,Split:=ifelse(seq_len(.N)%in%te,'Within-material test','Within-material training')]
fwrite(s[,.(SampleID,ReferenceKey,DonorKey,DonorKnown,ConnGroup,Split)],file.path(out,'tables/split.csv'))
train<-function(ii,mtry){
 z<-s[ii,.(ReferenceKey,DonorKey)];z[,W:=1/(.N),by=.(ReferenceKey,DonorKey)];z[,W:=W/uniqueN(DonorKey),by=ReferenceKey]
 ranger(x=x[ii,,drop=FALSE],y=y[ii],num.trees=500,mtry=mtry,min.node.size=1,
  case.weights=z$W,seed=25,num.threads=4,write.forest=TRUE,importance='none',oob.error=FALSE)
}
rows<-list();add<-function(model,ev,ii,pr){rows[[length(rows)+1]]<<-data.table(Model=model,Evaluation=ev,SampleID=s$SampleID[ii],ReferenceKey=s$ReferenceKey[ii],ConnGroup=s$ConnGroup[ii],Label=y[ii],Predicted=pr)}
for(mt in c(11L,64L)){
 model<-paste0('RF_mtry',mt);f<-train(seq_len(nrow(s)),mt)
 saveRDS(list(genes=genes,model=f,target='Shared material-level detected Kla percentage'),file.path(out,paste0(model,'.rds')))
 add(model,'Training fit',seq_len(nrow(s)),predict(f,x)$predictions)
 f<-train(tr,mt);add(model,'Within-material test',te,predict(f,x[te,,drop=FALSE])$predictions)
 for(b in sort(unique(s$ConnGroup))){ii<-which(s$ConnGroup!=b);jj<-which(s$ConnGroup==b);f<-train(ii,mt);add(model,'Source held out',jj,predict(f,x[jj,,drop=FALSE])$predictions)}
 cat('Completed',model,'\n')
}
p<-rbindlist(rows);p<-merge(p,s[,.(SampleID,DonorKey)],by='SampleID')
# Re-normalize donor/material weights separately within each evaluation subset.
p[,W:=1/.N,by=.(Model,Evaluation,ReferenceKey,DonorKey)];p[,W:=W/uniqueN(DonorKey),by=.(Model,Evaluation,ReferenceKey)]
metrics<-rbindlist(lapply(c('Record','Material'),function(level)p[,{
 w<-if(level=='Record')rep(1,.N)else W;w<-w/sum(w);e<-Predicted-Label
 list(Weighting=level,N=.N,Materials=uniqueN(ReferenceKey),MAE=sum(w*abs(e)),RMSE=sqrt(sum(w*e^2)),R2=1-sum(w*e^2)/sum(w*(Label-sum(w*Label))^2),Within2=100*sum(w*(abs(e)<=2)),Within5=100*sum(w*(abs(e)<=5)))
},by=.(Model,Evaluation)]))
fwrite(p,file.path(out,'tables/predictions.csv.gz'));fwrite(metrics,file.path(out,'tables/metrics.csv'))
# No-RNA reference for each corresponding evaluation; never uses held-out labels.
wm<-function(ii){z<-s[ii,.(Y=ObservedPercent),by=ReferenceKey][,.(Y=Y[1]),by=ReferenceKey];median(z$Y)}
bas<-copy(p[Model=='RF_mtry11']);bas[,Model:='NoRNA'];bas[Evaluation=='Training fit',Predicted:=wm(seq_len(nrow(s)))];bas[Evaluation=='Within-material test',Predicted:=wm(tr)]
for(b in unique(s$ConnGroup))bas[Evaluation=='Source held out' & ConnGroup==b,Predicted:=wm(which(s$ConnGroup!=b))]
fwrite(bas,file.path(out,'tables/baseline_predictions.csv.gz'))
fwrite(bas[,.(MAE=sum(W*abs(Predicted-Label))/sum(W),Within5=100*sum(W*(abs(Predicted-Label)<=5))/sum(W)),by=Evaluation],file.path(out,'tables/baseline_metrics.csv'))
for(model in unique(p$Model)){
 d<-merge(p[Model==model],metrics[Model==model & Weighting=='Material',.(Evaluation,MAE,Within5,N)],by='Evaluation')
 evs<-c('Training fit','Within-material test','Source held out');mm<-metrics[Model==model & Weighting=='Material'][match(evs,Evaluation)]
 panel<-function(e,n,ma,acc)sprintf('%s (n=%d)\nMAE %.2f pp | within +/-5 pp: %.1f%%',e,n,ma,acc)
 d[,Panel:=factor(panel(Evaluation,N,MAE,Within5),levels=panel(mm$Evaluation,mm$N,mm$MAE,mm$Within5))]
 g<-ggplot(d,aes(Label,Predicted))+geom_abline(slope=1,intercept=0,linetype=2,colour='grey55')+geom_point(size=.9,alpha=.3,colour='#2166AC')+facet_wrap(~Panel,nrow=1)+coord_equal(xlim=c(0,40),ylim=c(0,40))+labs(title=paste('Nonlinear prediction for individual RNA samples:',model),subtitle='128 RNA genes only | material-balanced random forest | 500 trees',x='Assigned material Kla detection percentage (%)',y='Predicted percentage (%)',caption='Each dot is an RNA record. Labels are shared material measurements, not individual protein measurements.\nMiddle: held-out RNA records from represented materials; right: entire unseen source blocks. Metrics weight materials equally.')+theme_minimal(base_size=10,base_family='Arial')+theme(panel.grid.minor=element_blank(),plot.caption=element_text(hjust=0))
 for(ext in c('pdf','png'))ggsave(file.path(out,'figures',paste0(model,'.',ext)),g,width=13,height=5.5,device=if(ext=='pdf')cairo_pdf else agg_png,dpi=180)
}
stopifnot(all(is.finite(p$Predicted)),!anyDuplicated(p,by=c('Model','Evaluation','SampleID')),all(p$Predicted>=0 & p$Predicted<=100))
fwrite(data.table(Role=names(paths),Path=unname(paths),SHA256=vapply(paths,digest::digest,character(1),algo='sha256',file=TRUE)),file.path(out,'input_sha256.csv'))
writeLines(trimws(capture.output(sessionInfo()),which='right'),file.path(out,'sessionInfo.txt'))
fs<-sort(list.files(out,recursive=TRUE,full.names=TRUE));fwrite(data.table(Path=substring(fs,nchar(out)+2),SHA256=vapply(fs,digest::digest,character(1),algo='sha256',file=TRUE)),file.path(out,'release_sha256.csv'))
print(metrics[Weighting=='Material'])
