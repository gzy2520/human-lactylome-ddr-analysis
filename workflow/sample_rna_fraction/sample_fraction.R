#!/usr/bin/env Rscript
# Sample-level weak supervision and distribution regression for frozen group labels.
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(ragg)})
set.seed(25)
args<-commandArgs(TRUE);stopifnot(length(args)==2); inp<-args[1];out<-args[2]
if(dir.exists(out)&&length(list.files(out,all.files=TRUE,no..=TRUE)))stop('Output exists')
dir.create(file.path(out,'tables'),recursive=TRUE);dir.create(file.path(out,'figures'))
put<-function(x,n)fwrite(x,file.path(out,'tables',paste0(n,'.csv')))
paths<-c(script='workflow/sample_rna_fraction/sample_fraction.R',
  rna=file.path(inp,'sample_rna128.csv.gz'),metadata=file.path(inp,'sample_metadata.csv'),
  labels=file.path(inp,'material_labels_28.csv'),
  original_labels='outputs/20260923_rna_only_kla_fraction/protein_group_labels_31.csv',
  modules='config/lactate_metabolism/gene_module_membership.csv',
  catalog='outputs/20260923_rna_only_kla_fraction_optimization/module_catalog.csv')
d<-fread(paths['labels']);setorder(d,ReferenceKey)
s<-fread(paths['rna']);meta<-fread(paths['metadata'])
s<-merge(s,meta[,.(SampleID,ConnGroup,DonorKey,MaterialBalancedWeight,SourceBalancedWeight)],by='SampleID',sort=FALSE)
setorder(s,ReferenceKey,SampleID)
genes<-sort(grep('^ENSG',names(s),value=TRUE));X<-as.matrix(s[,..genes])
bi<-match(s$ReferenceKey,d$ReferenceKey);blocks<-d$ConnGroup;y<-d$ObservedPercent
stopifnot(length(genes)==128,nrow(s)==1898,nrow(d)==28,uniqueN(blocks)==14,
 all(s$ConnGroup==blocks[bi]),all(is.finite(X)))
counts<-tabulate(bi,nbins=nrow(d));group_weight<-1/as.numeric(table(blocks)[as.character(blocks)])
sample_weight<-s$SourceBalancedWeight; within_weight<-s$MaterialBalancedWeight
stopifnot(max(abs(as.vector(rowsum(sample_weight,bi))-group_weight))<1e-12)
mem<-unique(fread(paths['modules'])[Ensembl %chin% genes,.(Origin,Module,Ensembl)])
catlog<-fread(paths['catalog']);mem<-merge(mem,catlog[,.(Origin,Module,FeatureID)],by=c('Origin','Module'))
modules<-lapply(split(mem$Ensembl,mem$FeatureID),function(z)match(z,genes))
stopifnot(length(modules)==16)
di<-as.integer(factor(paste(s$ReferenceKey,s$DonorKey,sep='::')))
donorX<-rowsum(X,di,reorder=TRUE)/as.numeric(table(di))
donorbi<-as.integer(tapply(bi,di,function(v)v[1]))
medX<-do.call(rbind,lapply(seq_len(nrow(d)),function(i)apply(donorX[donorbi==i,,drop=FALSE],2,median)))
# Fixed random Fourier feature map; no outcomes influence the nonlinear basis.
W<-matrix(rnorm(16*32),16,32)/sqrt(16);phase<-runif(32,0,2*pi)
wmean<-function(x,w)as.vector(crossprod(w/sum(w),x))
standard<-function(a,w){mu<-wmean(a,w);z<-sweep(a,2,mu,'-');sd<-sqrt(wmean(z*z,w));sd[sd<1e-8]<-1;list(mu=mu,sd=sd)}
transform<-function(x,z)sweep(sweep(x,2,z$mu,'-'),2,z$sd,'/')
bagmean<-function(x)rowsum(x*within_weight,bi,reorder=TRUE)
features<-function(train){
  ii<-which(bi %in% train);z<-standard(X[ii,,drop=FALSE],sample_weight[ii]);
  sx<-transform(X,z);raw128<-sx;gx<-transform(medX,z)
  score<-function(a)do.call(cbind,lapply(modules,function(ix)rowMeans(a[,ix,drop=FALSE])))
  sx<-score(sx);gx<-score(gx);z2<-standard(sx[ii,,drop=FALSE],sample_weight[ii])
  sx<-transform(sx,z2);gx<-transform(gx,z2)
  nonlinear<-sqrt(2/32)*cos(sweep(sx%*%W,2,phase,'+'))
  list(sample=sx,raw128=raw128,median=gx,mean=bagmean(sx),
       distribution=bagmean(cbind(sx/sqrt(16),nonlinear)/sqrt(2)))
}
fit_ridge<-function(a,b,w,lambda){
  w<-w/sum(w);cx<-wmean(a,w);cy<-sum(w*b);az<-sweep(a,2,cx,'-')
  e<-eigen(crossprod(az*sqrt(w)),symmetric=TRUE)
  rhs<-crossprod(az,w*(b-cy));beta<-e$vectors%*%((crossprod(e$vectors,rhs))/(pmax(e$values,0)+lambda))
  list(cx=cx,cy=cy,beta=beta)
}
predict_ridge<-function(f,a)as.vector(f$cy+sweep(a,2,f$cx,'-')%*%f$beta)
models<-c('MaterialMedian','SampleUnweighted','SampleMaterialBalanced','SampleBalanced',
          'SampleBalanced128','BagMean','BagDistribution')
lambdas<-c(100,10,1,.1,.01)
make_fit<-function(f,train,model,lambda){
  if(grepl('^Sample',model)){
    ii<-which(bi %in% train)
    ww<-if(model=='SampleUnweighted')rep(1,length(ii)) else if(model=='SampleMaterialBalanced')within_weight[ii] else sample_weight[ii]
    a<-if(model=='SampleBalanced128')f$raw128 else f$sample
    return(fit_ridge(a[ii,,drop=FALSE],y[bi[ii]],ww,lambda))
  }
  key<-switch(model,MaterialMedian='median',BagMean='mean',BagDistribution='distribution')
  fit_ridge(f[[key]][train,,drop=FALSE],y[train],group_weight[train],lambda)
}
pred_all<-function(f,fit,model){
  if(grepl('^Sample',model)){
    a<-if(model=='SampleBalanced128')f$raw128 else f$sample
    p<-predict_ridge(fit,a);return(list(bag=as.vector(rowsum(p*within_weight,bi,reorder=TRUE)),sample=p))
  }
  key<-switch(model,MaterialMedian='median',BagMean='mean',BagDistribution='distribution')
  list(bag=predict_ridge(fit,f[[key]]),sample=NULL)
}
clip<-function(p)pmin(100,pmax(0,p))
weighted_median<-function(v,w){o<-order(v);v[o][which(cumsum(w[o])>=sum(w)/2)[1]]}
preds<-list();tuning<-list();sample_preds<-list();fullmodels<-list()
for(outer in sort(unique(blocks))){
  train<-which(blocks!=outer);test<-which(blocks==outer)
  loss<-matrix(0,length(models),length(lambdas),dimnames=list(models,as.character(lambdas)))
  for(inner in sort(unique(blocks[train]))){
    tr<-train[blocks[train]!=inner];te<-train[blocks[train]==inner];f<-features(tr)
    for(m in models)for(j in seq_along(lambdas)){
      fit<-make_fit(f,tr,m,lambdas[j]);p<-clip(pred_all(f,fit,m)$bag[te])
      loss[m,j]<-loss[m,j]+mean(abs(p-y[te]))
    }
  }
  f<-features(train)
  for(m in models){
    j<-which.min(loss[m,]);fit<-make_fit(f,train,m,lambdas[j]);p<-pred_all(f,fit,m)
    preds[[length(preds)+1]]<-data.table(ReferenceKey=d$ReferenceKey[test],ConnGroup=outer,
      Model=m,Observed=y[test],Predicted=clip(p$bag[test]))
    tuning[[length(tuning)+1]]<-data.table(OuterBlock=outer,Model=m,Lambda=lambdas,
      InnerSourceMAE=loss[m,]/length(unique(blocks[train])),Selected=seq_along(lambdas)==j)
    if(!is.null(p$sample)){
      ii<-which(bi %in% test)
      sample_preds[[length(sample_preds)+1]]<-data.table(Model=m,ReferenceKey=s$ReferenceKey[ii],
        SampleID=s$SampleID[ii],ConnGroup=outer,SharedGroupLabel=y[bi[ii]],RawPrediction=p$sample[ii])
    }
  }
  preds[[length(preds)+1]]<-data.table(ReferenceKey=d$ReferenceKey[test],ConnGroup=outer,
    Model='NoRNA',Observed=y[test],Predicted=weighted_median(y[train],group_weight[train]))
  cat('Completed source',outer,'\n')
}
p<-rbindlist(preds);tun<-rbindlist(tuning)
stopifnot(nrow(p)==28*(length(models)+1),!anyDuplicated(p,by=c('Model','ReferenceKey')),all(is.finite(p$Predicted)))
put(p,'heldout_material_predictions');put(tun,'inner_tuning');put(rbindlist(sample_preds),'heldout_sample_weak_label_predictions')
by_source<-p[,.(MAE=mean(abs(Predicted-Observed))),by=.(Model,ConnGroup)];put(by_source,'source_errors')
metrics<-p[,.(MaterialMAE=mean(abs(Predicted-Observed)),RMSE=sqrt(mean((Predicted-Observed)^2))),by=Model]
metrics<-merge(metrics,by_source[,.(SourceMAE=mean(MAE)),by=Model],by='Model')
put(metrics,'performance')
put(data.table(ReferenceKey=d$ReferenceKey,ConnGroup=blocks,Samples=counts,GroupLabel=y,
              GroupWeight=group_weight),'weight_contract')
original31<-fread(paths['original_labels'])
original31<-merge(original31[,.(GroupID,ReferenceKey,OriginalPercent=DetectedKlaPercent)],
                  p[,.(ReferenceKey,Model,Predicted)],by='ReferenceKey',allow.cartesian=TRUE)
put(original31,'original31_label_predictions')
f<-features(seq_len(nrow(d)));fitrows<-list()
for(m in models){
  # Median selected penalty is only for a descriptive full-data fit; no test metric is recomputed.
  lam<-median(tun[Model==m & Selected==TRUE,Lambda]);fit<-make_fit(f,seq_len(nrow(d)),m,lam)
  fitrows[[m]]<-data.table(Model=m,ReferenceKey=d$ReferenceKey,Observed=y,
    Fitted=clip(pred_all(f,fit,m)$bag),Lambda=lam)
  fullmodels[[m]]<-fit
}
put(rbindlist(fitrows),'descriptive_full_fit')
labels<-c(NoRNA='No RNA',MaterialMedian='Material median',SampleUnweighted='Samples: unweighted',
          SampleMaterialBalanced='Samples: material balanced',SampleBalanced='Samples: source balanced',
          SampleBalanced128='Samples: 128 genes balanced',BagMean='Bags: linear mean',BagDistribution='Bags: nonlinear distribution')
p<-merge(p,metrics[,.(Model,SourceMAE)],by='Model')
p[,Panel:=factor(sprintf('%s\nSource MAE %.2f pp',labels[Model],SourceMAE),
  levels=sprintf('%s\nSource MAE %.2f pp',labels[c('NoRNA',models)],metrics$SourceMAE[match(c('NoRNA',models),metrics$Model)]))]
g<-ggplot(p,aes(Observed,Predicted))+geom_abline(slope=1,intercept=0,colour='grey65')+
  geom_point(colour='#2166AC',size=1.7)+facet_wrap(~Panel,ncol=3)+coord_equal(xlim=c(0,40),ylim=c(0,40))+
  labs(title='Does sample-level RNA improve group-percentage prediction?',
       subtitle='1,898 RNA records | 28 group labels | 14 donor-aware source-held-out blocks',
       x='Observed group detection percentage',y='Held-out group prediction',
       caption='All models use the same source splits and target labels. Sample scores are averaged within material before evaluation.\nNo individual RNA sample has a matched protein percentage; bag regression uses sample distributions.')+
  theme_minimal(base_size=10,base_family='Arial')+theme(panel.grid.minor=element_blank(),plot.caption=element_text(hjust=0))
ggsave(file.path(out,'figures/sample_vs_material.pdf'),g,width=11,height=10,device=cairo_pdf)
ggsave(file.path(out,'figures/sample_vs_material.png'),g,width=11,height=10,device=agg_png,dpi=160)
put(data.table(Role=names(paths),Path=unname(paths),SHA256=vapply(paths,digest::digest,character(1),algo='sha256',file=TRUE,USE.NAMES=FALSE)),'input_sha256')
cat('PASS: all 1,898 records evaluated by held-out material/source; no patient-level accuracy claimed.\n')

writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
