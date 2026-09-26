#!/usr/bin/env Rscript
# Reconstruct the existing direct-sample full fits; reuse frozen source-held-out predictions.
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(ragg)})
set.seed(25)
a<-commandArgs(TRUE);stopifnot(length(a)==3);inp<-a[1];old<-a[2];out<-a[3]
if(dir.exists(out))stop('Use a fresh output directory')
dir.create(file.path(out,'tables'),recursive=TRUE);dir.create(file.path(out,'figures'))
paths<-c(script='workflow/sample_rna_fraction/render_individual_samples.R',
 rna=file.path(inp,'sample_rna128.csv.gz'),metadata=file.path(inp,'sample_metadata.csv'),
 heldout=file.path(old,'tables/heldout_sample_weak_label_predictions.csv'),
 tuning=file.path(old,'tables/inner_tuning.csv'),fullfit=file.path(old,'tables/descriptive_full_fit.csv'),
 modules='config/lactate_metabolism/gene_module_membership.csv',catalog='outputs/20260923_rna_only_kla_fraction_optimization/module_catalog.csv')
m<-fread(paths['metadata']);s<-merge(fread(paths['rna']),m[,.(SampleID,ConnGroup,ObservedPercent,DonorKey,Category,MaterialBalancedWeight,SourceBalancedWeight)],by='SampleID')
setorder(s,ReferenceKey,SampleID);genes<-sort(grep('^ENSG',names(s),value=TRUE));x<-as.matrix(s[,..genes])
stopifnot(nrow(s)==1898,length(genes)==128,!anyDuplicated(s$SampleID))
wmean<-function(x,w)as.vector(crossprod(w/sum(w),x))
standard<-function(x,w){mu<-wmean(x,w);z<-sweep(x,2,mu,'-');sd<-sqrt(wmean(z*z,w));sd[sd<1e-8]<-1;list(mu=mu,sd=sd)}
trans<-function(x,z)sweep(sweep(x,2,z$mu,'-'),2,z$sd,'/')
z<-standard(x,s$SourceBalancedWeight);gx<-trans(x,z)
mem<-unique(fread(paths['modules'])[Ensembl %chin% genes,.(Origin,Module,Ensembl)])
mem<-merge(mem,fread(paths['catalog'])[,.(Origin,Module,FeatureID)],by=c('Origin','Module'))
modules<-lapply(split(mem$Ensembl,mem$FeatureID),function(v)match(v,genes))
mx<-do.call(cbind,lapply(modules,function(i)rowMeans(gx[,i,drop=FALSE])))
z2<-standard(mx,s$SourceBalancedWeight);mx<-trans(mx,z2)
models<-c('SampleBalanced128','SampleMaterialBalanced','SampleBalanced','SampleUnweighted')
names0<-c(SampleBalanced128='128 genes | source-balanced training',SampleMaterialBalanced='16 modules | material-balanced training',SampleBalanced='16 modules | source-balanced training',SampleUnweighted='16 modules | unweighted training')
tun<-fread(paths['tuning']);rows<-list();fits<-list()
for(model in models){
 xx<-if(model=='SampleBalanced128')gx else mx
 w<-switch(model,SampleUnweighted=rep(1,nrow(s)),SampleMaterialBalanced=s$MaterialBalancedWeight,s$SourceBalancedWeight);w<-w/sum(w)
 lam<-median(tun[Model==model & Selected==TRUE,Lambda]);cx<-wmean(xx,w);cy<-sum(w*s$ObservedPercent);az<-sweep(xx,2,cx,'-')
 e<-eigen(crossprod(az*sqrt(w)),symmetric=TRUE);rhs<-crossprod(az,w*(s$ObservedPercent-cy))
 beta<-e$vectors%*%(crossprod(e$vectors,rhs)/(pmax(e$values,0)+lam));pred<-as.vector(cy+az%*%beta)
 rows[[model]]<-data.table(Model=model,ReferenceKey=s$ReferenceKey,SampleID=s$SampleID,ConnGroup=s$ConnGroup,SharedGroupLabel=s$ObservedPercent,RawPrediction=pred,Evaluation='Training fit')
 fits[[model]]<-list(cx=cx,cy=cy,beta=beta,lambda=lam)
}
p<-rbindlist(rows);prior<-fread(paths['fullfit'])
check<-merge(p,s[,.(SampleID,MaterialBalancedWeight)],by='SampleID')[,.(Rebuilt=pmin(100,pmax(0,sum(RawPrediction*MaterialBalancedWeight)))),by=.(Model,ReferenceKey)]
check<-merge(check,prior,by=c('Model','ReferenceKey'));stopifnot(nrow(check)==112,max(abs(check$Rebuilt-check$Fitted))<1e-10)
h<-fread(paths['heldout']);h[,Evaluation:='Source held out'];p<-rbind(p,h,use.names=TRUE)
p<-merge(p,m[,.(SampleID,Category,MaterialBalancedWeight,SourceBalancedWeight)],by='SampleID')
p[,Prediction:=pmin(100,pmax(0,RawPrediction))]
stopifnot(all(p[,.(N=.N,U=uniqueN(SampleID)),by=.(Model,Evaluation)]$N==1898),!anyDuplicated(p,by=c('Model','Evaluation','SampleID')))
p[,Evaluation:=factor(Evaluation,levels=c('Training fit','Source held out'))]
metrics<-rbindlist(lapply(c('Record','Material','Source'),function(level){
 p[,{
  w<-switch(level,Record=rep(1,.N),Material=MaterialBalancedWeight,Source=SourceBalancedWeight);w<-w/sum(w)
  list(Weighting=level,MAE=sum(w*abs(Prediction-SharedGroupLabel)),R2=1-sum(w*(Prediction-SharedGroupLabel)^2)/sum(w*(SharedGroupLabel-sum(w*SharedGroupLabel))^2))
 },by=.(Model,Evaluation)]
}))
fwrite(p,file.path(out,'tables/individual_predictions.csv.gz'));fwrite(metrics,file.path(out,'tables/shared_label_agreement.csv'))
saveRDS(list(genes=genes,gene_scaling=z,module_indices=modules,module_scaling=z2,models=fits),file.path(out,'direct_sample_models.rds'))
base<-theme_minimal(base_size=11,base_family='Arial')+theme(panel.grid.minor=element_blank(),plot.caption=element_text(hjust=0))
cols<-c(normal_tissue='#2166AC',cancer_tissue='#B2182B',cancer_cells='#D89000',normal_cells='#1B9E77')
stopifnot(all(unique(p$Category)%in%names(cols)))
export<-function(g,name,w,h){ggsave(file.path(out,'figures',paste0(name,'.pdf')),g,width=w,height=h,device=cairo_pdf);ggsave(file.path(out,'figures',paste0(name,'.png')),g,width=w,height=h,device=agg_png,dpi=160)}
for(model in models){
 dd<-p[Model==model];mm<-metrics[Model==model & Weighting=='Material']
 dd<-merge(dd,mm[,.(Evaluation,MAE)],by='Evaluation');dd[,Panel:=factor(sprintf('%s | shared-label MAE %.2f pp',Evaluation,MAE),levels=sprintf('%s | shared-label MAE %.2f pp',mm$Evaluation,mm$MAE))]
 g<-ggplot(dd,aes(SharedGroupLabel,Prediction,colour=Category))+geom_abline(slope=1,intercept=0,linetype=2,colour='grey50')+geom_point(alpha=.3,size=1.1)+facet_wrap(~Panel,nrow=1)+coord_equal(xlim=c(0,40),ylim=c(0,40))+scale_colour_manual(values=cols)+labs(title='One prediction per RNA sample',subtitle=paste(names0[model], '| 1,898 sample records in EACH panel'),x='Assigned material-level Kla detection percentage (%)',y='Predicted percentage for each RNA sample (%)',colour=NULL,caption='Each point is one RNA record; exact overlaps can hide points. No horizontal jitter.\nSamples share their material label, not individually measured protein percentages. MAE gives each material equal total weight.')+base+theme(legend.position='bottom')
 export(g,paste0(model,'_sample_scatter'),11,6.5)
 ord<-unique(dd[,.(ReferenceKey,SharedGroupLabel)])[order(SharedGroupLabel),ReferenceKey]
 dd[,Material:=factor(ReferenceKey,levels=rev(ord))]
 lab<-dd[,.(N=.N),by=.(Evaluation,Material)][,.(N=N[1]),by=Material]
 ticks<-setNames(paste0(as.character(lab$Material),' (n=',lab$N,')'),as.character(lab$Material))
 targets<-unique(dd[,.(Evaluation,Material,SharedGroupLabel)])
 g<-ggplot(dd,aes(Prediction,Material))+geom_boxplot(outlier.shape=NA,width=.5,fill='#DCE9F3',colour='#52758E')+geom_point(position=position_jitter(width=0,height=.15,seed=25),size=.65,alpha=.25,colour='#2166AC')+geom_point(data=targets,aes(x=SharedGroupLabel),shape=18,size=2.5,colour='#B2182B')+facet_wrap(~Evaluation,nrow=1)+scale_y_discrete(labels=ticks)+labs(title='RNA sample predictions within each material',subtitle=names0[model],x='Predicted Kla detection percentage (%)',y=NULL,caption='Blue dots: individual RNA predictions; boxes: median and interquartile range; red diamonds: shared material labels.\nVertical jitter separates records visually; predicted values are not shifted. n counts RNA records.')+base
 export(g,paste0(model,'_material_distributions'),13,10)
}
fwrite(data.table(Role=names(paths),Path=unname(paths),SHA256=vapply(paths,digest::digest,character(1),algo='sha256',file=TRUE)),file.path(out,'input_sha256.csv'))
writeLines(c('PASS: four direct-sample models, 1,898 rows per model per evaluation','PASS: reconstructed training predictions reproduce all 112 frozen material fits within 1e-10','PASS: source-held-out sample predictions reused unchanged','PASS: no material aggregation before prediction; clipping only for percentage display and sample agreement metrics'),file.path(out,'validation.txt'))
writeLines(trimws(capture.output(sessionInfo()),which='right'),file.path(out,'sessionInfo.txt'))
fs<-sort(list.files(out,recursive=TRUE,full.names=TRUE));fwrite(data.table(Path=substring(fs,nchar(out)+2),SHA256=vapply(fs,digest::digest,character(1),algo='sha256',file=TRUE)),file.path(out,'release_sha256.csv'))
print(metrics[Weighting=='Material']);cat('PASS: individual sample figures exported\n')
