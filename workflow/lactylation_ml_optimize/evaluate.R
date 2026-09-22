source('workflow/lactylation_ml_optimize/engine.R')
setDTthreads(1)
out <- 'outputs/20260922_br_v5'
p <- fread(file.path(out,'oos_predictions.csv.gz'))
models <- c('M0','MP','MP_cal','ME_P','ME2P','SHRINK','SHRINK_CAL')
# Explicit loop avoids duplicate group columns from the metric helper.
fold <- rbindlist(lapply(models,function(m)rbindlist(lapply(sort(unique(p$ConnGroup)),function(g){z<-p[ConnGroup==g];metrics_row_v3(z$y,z[[paste0('p_',m)]],m,g)}))))
summary <- fold[,lapply(.SD,mean),by=Model,.SDcols=c('LogLoss','AP','ROC','Brier')]
cal <- rbindlist(lapply(models,function(m)p[,.(N=.N,Observed=mean(y),Predicted=mean(get(paste0('p_',m)))),by=ReferenceKey][,Model:=m]))
cal[,AbsoluteError:=abs(Predicted-Observed)]
summary <- merge(summary,cal[,.(MaterialMAE=mean(AbsoluteError)),by=Model],by='Model',sort=FALSE)
delta <- merge(fold[Model%in%c('SHRINK','SHRINK_CAL'),.(Model,ConnGroup,LogLoss)],fold[Model=='MP_cal',.(ConnGroup,Baseline=LogLoss)],by='ConnGroup')
delta[,Delta:=LogLoss-Baseline]
set.seed(25)
uncertainty <- delta[,{
  draws<-replicate(10000,mean(sample(Delta,length(Delta),replace=TRUE)))
  list(MeanDelta=mean(Delta),Better=sum(Delta<0),Worse=sum(Delta>0),BootstrapLow=unname(quantile(draws,.025)),BootstrapHigh=unname(quantile(draws,.975)),DeleteOneLow=min((sum(Delta)-Delta)/(length(Delta)-1)),DeleteOneHigh=max((sum(Delta)-Delta)/(length(Delta)-1)))
},by=Model]
fwrite(fold,file.path(out,'per_group_metrics.csv'));fwrite(summary,file.path(out,'performance_summary.csv'))
fwrite(cal,file.path(out,'material_calibration.csv'));fwrite(delta,file.path(out,'paired_logloss.csv'));fwrite(uncertainty,file.path(out,'conditional_uncertainty.csv'))
draw <- function(){
 par(mfrow=c(1,2),mar=c(5,5,3,1))
 for(m in c('SHRINK','SHRINK_CAL')){
  z<-delta[Model==m][order(ConnGroup)]
  plot(seq_len(nrow(z)),z$Delta,pch=19,col=ifelse(z$Delta<0,'#238B45','#CB181D'),xlab='Held-out connected group',ylab='LogLoss difference vs MP_cal',main=m,xaxt='n');axis(1,seq_len(nrow(z)),z$ConnGroup,cex.axis=.7);abline(h=0,lty=2,col='gray40')
 }
}
pdf(file.path(out,'fig_paired_logloss.pdf'),width=11,height=4.5);draw();dev.off()
png(file.path(out,'fig_paired_logloss.png'),width=1650,height=675,res=150);draw();dev.off()
drawcal <- function(){par(mfrow=c(1,3),mar=c(4,4,3,1));for(m in c('MP_cal','SHRINK','SHRINK_CAL')){z<-cal[Model==m];plot(z$Observed,z$Predicted,xlim=c(0,1),ylim=c(0,1),pch=19,col='#2C7FB8',xlab='Observed material detection rate',ylab='Mean predicted probability',main=m);abline(0,1,lty=2,col='gray40')}}
pdf(file.path(out,'fig_material_calibration.pdf'),width=11,height=4);drawcal();dev.off()
png(file.path(out,'fig_material_calibration.png'),width=1650,height=600,res=150);drawcal();dev.off()
print(summary);print(uncertainty)
