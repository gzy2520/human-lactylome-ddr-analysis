library(data.table)
source('workflow/lactylation_ml_review/05c_metrics_and_diagnostics_v3.R')
v5_grid<-function(){rbind(data.table(ID='original',Weight='row',Center='half',Strength=2),CJ(Weight=c('row','group'),Strength=c(.5,2,8,32))[,.(ID=paste(Weight,Strength,sep='_'),Weight,Center='train_group',Strength)])}
v5_folds<-function(groups,k=5L){set.seed(25);g<-sample(sort(unique(groups)));split(g,rep(seq_len(min(k,length(g))),length.out=length(g)))}
v5_weights<-function(d){n<-table(d$ConnGroup);as.numeric(1/n[as.character(d$ConnGroup)])/length(n)}
v5_freq<-function(tr,ev,cfg){
 stopifnot(nrow(cfg)==1,nrow(tr)>0)
 group_rate<-tr[,.(rate=mean(y)),by=ConnGroup][,mean(rate)]
 center<-if(cfg$Center=='half') .5 else group_rate
 if(cfg$Weight=='row') st<-tr[,.(n=.N,k=sum(y)),by=BaseAccession] else {
  votes<-tr[,.(v=mean(y)),by=.(ConnGroup,BaseAccession)]
  st<-votes[,.(n=.N,k=sum(v)),by=BaseAccession]
 }
 j<-match(ev$BaseAccession,st$BaseAccession)
 p<-(st$k[j]+cfg$Strength*center)/(st$n[j]+cfg$Strength)
 fallback<-if(cfg$ID=='original') (sum(tr$y)+1)/(nrow(tr)+2) else group_rate
 unseen<-is.na(j);p[unseen]<-fallback
 list(p=pmin(pmax(p,1e-9),1-1e-9),unseen=unseen)
}
v5_cf<-function(d,cfg,k=4L){
 p<-rep(NA_real_,nrow(d));folds<-v5_folds(d$ConnGroup,k)
 for(gs in folds){idx<-which(d$ConnGroup%in%gs);p[idx]<-v5_freq(d[-idx],d[idx,.(BaseAccession)],cfg)$p}
 stopifnot(all(is.finite(p)));p
}
v5_cal<-function(d,p){
 x<-qlogis(p);w<-v5_weights(d);y<-d$y
 obj<-function(b){eta<-b[1]+b[2]*x;sum(w*(pmax(eta,0)+log1p(exp(-abs(eta)))-y*eta))+.01*(b[2]-1)^2}
 grad<-function(b){r<-w*(plogis(b[1]+b[2]*x)-y);c(sum(r),sum(r*x)+.02*(b[2]-1))}
 f<-optim(c(0,1),obj,gr=grad,method='BFGS',control=list(reltol=1e-10,maxit=500))
 stopifnot(f$convergence==0,all(is.finite(f$par)))
 list(par=f$par,objective=f$value,gradient=max(abs(grad(f$par))))
}
v5_apply<-function(f,p)pmin(pmax(plogis(f$par[1]+f$par[2]*qlogis(p)),1e-9),1-1e-9)
v5_group_loss<-function(d,p){data.table(ConnGroup=d$ConnGroup,l= -d$y*log(p)-(1-d$y)*log1p(-p))[,.(LogLoss=mean(l)),by=ConnGroup]}
v5_fit<-function(tr,ev,grid=v5_grid()){
 folds<-v5_folds(tr$ConnGroup);logs<-list()
 for(i in seq_len(nrow(grid))){cfg<-grid[i]
  for(f in seq_along(folds)){
   a<-tr[!ConnGroup%in%folds[[f]]];v<-tr[ConnGroup%in%folds[[f]]]
   pv<-v5_freq(a,v[,.(BaseAccession)],cfg)$p
   fit<-v5_cal(a,v5_cf(a,cfg))
   for(m in c('SHRINK','SHRINK_CAL')){
    z<-v5_group_loss(v,if(m=='SHRINK') pv else v5_apply(fit,pv));z[,`:=`(Model=m,Candidate=cfg$ID,InnerFold=f)]
    logs[[paste(i,f,m)]]<-z
   }
  }
 }
 logs<-rbindlist(logs)
 risks<-logs[,.(Risk=mean(LogLoss),Groups=.N),by=.(Model,Candidate)]
 stopifnot(all(risks$Groups==uniqueN(tr$ConnGroup)))
 pred<-ev[,.(GroupID,ReferenceKey,ConnGroup,BaseAccession,EnsemblGeneID)]
 selected<-list()
 for(m in c('SHRINK','SHRINK_CAL')){
  id<-risks[Model==m][order(Risk,Candidate)][1,Candidate];cfg<-grid[ID==id]
  raw<-v5_freq(tr,ev[,.(BaseAccession)],cfg)
  cal<-if(m=='SHRINK_CAL') v5_cal(tr,v5_cf(tr,cfg,k=5L)) else list(par=c(0,1),gradient=NA_real_)
  pred[,(paste0('p_',m)):=if(m=='SHRINK')raw$p else v5_apply(cal,raw$p)]
  selected[[m]]<-data.table(Model=m,Candidate=id,Intercept=cal$par[1],Slope=cal$par[2],Gradient=cal$gradient,UnseenRate=mean(raw$unseen))
 }
 pred[,p_MP_check:=v5_freq(tr,ev[,.(BaseAccession)],grid[ID=='original'])$p]
 list(pred=pred,selected=rbindlist(selected),inner=logs,risks=risks)
}
