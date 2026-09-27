#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table);library(ranger);library(ggplot2);library(ragg)})
set.seed(25)
a<-commandArgs(TRUE);stopifnot(length(a)==3);inp<-a[1];out<-a[2];target<-a[3];stopifnot(target %in% c('DDR','DNA_repair'));if(dir.exists(out))stop('Fresh output required')
dir.create(file.path(out,'tables'),recursive=TRUE);dir.create(file.path(out,'figures'))
paths<-c(script='workflow/ddr_fraction/train.R',labels='outputs/20260927_ddr_fraction_labels/material_labels_28.csv',rna=file.path(inp,'full_rna.rds'),split='outputs/20260926_nonlinear_samples_final/tables/split.csv')
bundle<-readRDS(paths['rna']);s<-bundle$metadata;x<-as.data.frame(bundle$x);genes<-colnames(x);lab<-fread(paths['labels'])[Target==target];s[,ObservedPercent:=lab$ObservedPercent[match(ReferenceKey,lab$ReferenceKey)]];y<-s$ObservedPercent;stopifnot(all(is.finite(y)));rm(bundle);gc(FALSE)
stopifnot(nrow(s)==1898,length(genes)>128,!anyDuplicated(s$SampleID),identical(rownames(x),s$SampleID))
frozen<-fread(paths['split']);frozen<-frozen[match(s$SampleID,SampleID)]
stopifnot(identical(frozen$SampleID,s$SampleID),identical(frozen$ReferenceKey,s$ReferenceKey),all(frozen$ConnGroup==s$ConnGroup))
te<-which(frozen$Split=='Within-material test');tr<-setdiff(seq_len(nrow(s)),te)
stopifnot(length(te)==374,length(tr)==1524,!length(intersect(s$DonorKey[tr],s$DonorKey[te])))
s[,Split:=ifelse(seq_len(.N)%in%te,'Within-material test','Within-material training')]
fwrite(s[,.(SampleID,ReferenceKey,DonorKey,DonorKnown,ConnGroup,Split)],file.path(out,'tables/split.csv'))
train<-function(ii,mtry){
 z<-s[ii,.(ReferenceKey,DonorKey)];z[,W:=1/(.N),by=.(ReferenceKey,DonorKey)];z[,W:=W/uniqueN(DonorKey),by=ReferenceKey]
 ranger(x=x[ii,,drop=FALSE],y=y[ii],num.trees=500,mtry=mtry,min.node.size=1,
  case.weights=z$W,seed=25,num.threads=4,write.forest=TRUE,importance='impurity',oob.error=FALSE)
}
rows<-list();add<-function(model,ev,ii,pr){rows[[length(rows)+1]]<<-data.table(Model=model,Evaluation=ev,SampleID=s$SampleID[ii],ReferenceKey=s$ReferenceKey[ii],ConnGroup=s$ConnGroup[ii],Label=y[ii],Predicted=pr)}
for(mt in 135L){
 model<-paste0('RF_mtry',mt);f<-train(seq_len(nrow(s)),mt)
 saveRDS(list(genes=genes,model=f,target=paste('Shared material-level Kla',target,'percentage')),file.path(out,paste0(model,'.rds')))
 add(model,'Training fit',seq_len(nrow(s)),predict(f,x)$predictions)
 f<-train(tr,mt);add(model,'Training split',tr,predict(f,x[tr,,drop=FALSE])$predictions);
 imp<-data.table(Ensembl=names(f$variable.importance),Importance=as.numeric(f$variable.importance));setorder(imp,-Importance,Ensembl);imp[,Rank:=.I];fwrite(imp,file.path(out,'tables/training_importance.csv'));
 add(model,'Within-material test',te,predict(f,x[te,,drop=FALSE])$predictions)
 for(b in sort(unique(s$ConnGroup))){ii<-which(s$ConnGroup!=b);jj<-which(s$ConnGroup==b);f<-train(ii,mt);add(model,'Source held out',jj,predict(f,x[jj,,drop=FALSE])$predictions);cat(model,'source block',b,'completed\n')}
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
bas<-copy(p[Model=='RF_mtry135']);bas[,Model:='NoRNA'];bas[Evaluation=='Training fit',Predicted:=wm(seq_len(nrow(s)))];bas[Evaluation %in% c('Within-material test','Training split'),Predicted:=wm(tr)]
for(b in unique(s$ConnGroup))bas[Evaluation=='Source held out' & ConnGroup==b,Predicted:=wm(which(s$ConnGroup!=b))]
fwrite(bas,file.path(out,'tables/baseline_predictions.csv.gz'))
fwrite(bas[,.(MAE=sum(W*abs(Predicted-Label))/sum(W),Within5=100*sum(W*(abs(Predicted-Label)<=5))/sum(W)),by=Evaluation],file.path(out,'tables/baseline_metrics.csv'))
for(model in unique(p$Model)){
 d<-merge(p[Model==model],metrics[Model==model & Weighting=='Material',.(Evaluation,MAE,Within5,N)],by='Evaluation')
 evs<-c('Training split','Within-material test','Training fit','Source held out');mm<-metrics[Model==model & Weighting=='Material'][match(evs,Evaluation)]
 panel<-function(e,n,ma,acc)sprintf('%s (n=%d)\nMAE %.2f pp | within +/-5 pp: %.1f%%',e,n,ma,acc)
 d[,Panel:=factor(panel(Evaluation,N,MAE,Within5),levels=panel(mm$Evaluation,mm$N,mm$MAE,mm$Within5))]
 g<-ggplot(d,aes(Label,Predicted))+geom_abline(slope=1,intercept=0,linetype=2,colour='grey55')+geom_point(size=.9,alpha=.3,colour='#2166AC')+facet_wrap(~Panel,nrow=1)+coord_equal(xlim=c(0,100),ylim=c(0,100))+labs(title=paste('RNA prediction of Kla-positive',target,'/ detected',target),subtitle=paste(length(genes),'common RNA genes | material-balanced random forest | 500 trees'),x=paste('Assigned material Kla-',target,' percentage (%)',sep=''),y='Predicted percentage (%)',caption='Each dot is an RNA record. Labels are shared material measurements, not individual protein measurements.\nInternal test: held-out RNA records from represented materials; source test: entire unseen source blocks. Metrics weight materials equally.')+theme_minimal(base_size=10,base_family='Arial')+theme(panel.grid.minor=element_blank(),plot.caption=element_text(hjust=0))
 for(ext in c('pdf','png'))ggsave(file.path(out,'figures',paste0(model,'.',ext)),g,width=16,height=5.8,device=if(ext=='pdf')cairo_pdf else agg_png,dpi=180)
}
stopifnot(all(is.finite(p$Predicted)),!anyDuplicated(p,by=c('Model','Evaluation','SampleID')),all(p$Predicted>=0 & p$Predicted<=100))
fwrite(lab,file.path(out,'tables/material_target_labels.csv'))
pp<-copy(p);pp[,AbsError:=abs(Predicted-Label)]
fwrite(pp[Evaluation=='Within-material test',.(N=.N,Label=Label[1],MeanPrediction=mean(Predicted),MAE=mean(AbsError),Within5=100*mean(AbsError<=5)),by=ReferenceKey],file.path(out,'tables/internal_test_by_material.csv'))
imp<-fread(file.path(out,'tables/training_importance.csv'))[1:25]
g<-ggplot(imp,aes(Importance,reorder(Ensembl,Importance)))+geom_col(fill='#2166AC')+labs(title=paste(target,'model: top 25 training-split features'),subtitle='Ranked across all 18,332 input genes; regression variance-reduction importance',x='Impurity importance (not causal effect)',y='Ensembl ID')+theme_minimal(base_size=11,base_family='Arial')
ggsave(file.path(out,'figures/top25.png'),g,width=9,height=8,device=agg_png,dpi=180)
ggsave(file.path(out,'figures/top25.pdf'),g,width=9,height=8,device=cairo_pdf)
fwrite(data.table(Role=names(paths),Path=unname(paths),SHA256=vapply(paths,digest::digest,character(1),algo='sha256',file=TRUE)),file.path(out,'input_sha256.csv'))
writeLines(trimws(capture.output(sessionInfo()),which='right'),file.path(out,'sessionInfo.txt'))
fs<-sort(list.files(out,recursive=TRUE,full.names=TRUE));fwrite(data.table(Path=substring(fs,nchar(out)+2),SHA256=vapply(fs,digest::digest,character(1),algo='sha256',file=TRUE)),file.path(out,'release_sha256.csv'))
print(metrics[Weighting=='Material'])
