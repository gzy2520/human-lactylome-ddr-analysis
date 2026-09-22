# Post-main-result robustness check; no candidate grid or fitting changes.
source('workflow/lactylation_ml_optimize/engine.R')
setDTthreads(1)
d <- fread('outputs/20260921_br_v4/BRef_features_v4.csv.gz')
base <- fread('outputs/20260922_br_v5/oos_predictions.csv.gz')
gs <- sort(unique(d$ConnGroup));records <- choices <- list()
for(drop in gs){
 cat('Delete group',drop,'\n');flush.console()
 for(test in setdiff(gs,drop)){
  ev<-d[ConnGroup==test];f<-v5_fit(d[!ConnGroup%in%c(drop,test)],ev)
  original<-merge(ev[,.(GroupID,BaseAccession)],base[ConnGroup==test],by=c('GroupID','BaseAccession'),sort=FALSE)
  j<-match(paste(ev$GroupID,ev$BaseAccession),paste(original$GroupID,original$BaseAccession))
  for(m in c('SHRINK','SHRINK_CAL')){
   col<-paste0('p_',m)
   records[[paste(drop,test,m)]]<-data.table(Dropped=drop,Test=test,Model=m,LogLoss=logloss(ev$y,f$pred[[col]]),MPLogLoss=logloss(ev$y,f$pred$p_MP_check),OriginalLogLoss=logloss(ev$y,original[[col]][j]))
  }
  f$selected[,`:=`(Dropped=drop,Test=test)];choices[[paste(drop,test)]]<-f$selected
 }
}
r<-rbindlist(records);r[,`:=`(DeltaVsMP=LogLoss-MPLogLoss,ChangeVsOriginal=LogLoss-OriginalLogLoss)]
fwrite(r,'outputs/20260922_br_v5/s2_per_group.csv')
fwrite(r[,.(MeanLogLoss=mean(LogLoss),MeanDeltaVsMP=mean(DeltaVsMP),MeanChangeVsOriginal=mean(ChangeVsOriginal)),by=.(Dropped,Model)],'outputs/20260922_br_v5/s2_summary.csv')
fwrite(rbindlist(choices),'outputs/20260922_br_v5/s2_selected_parameters.csv')
