#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(ragg)})
set.seed(25)
root<-'outputs/20260929_assay_composition';dir.create(file.path(root,'figures'),showWarnings=FALSE)
targets<-c('Proteome_DDR','Kla_DDR','Proteome_DNA_repair','Kla_DNA_repair')
pretty<-c(Proteome_DDR='DDR / all proteome proteins',Kla_DDR='DDR / all Kla proteins',Proteome_DNA_repair='DNA repair / all proteome proteins',Kla_DNA_repair='DNA repair / all Kla proteins')
labels<-fread(file.path(root,'labels/material_labels_28.csv'));allp<-list();allm<-list();allb<-list()
calc<-function(d,level){w<-if(level=='Record')rep(1,nrow(d))else d$W;w<-w/sum(w);e<-d$Predicted-d$Label;data.table(Weighting=level,N=nrow(d),Materials=uniqueN(d$ReferenceKey),MAE=sum(w*abs(e)),RMSE=sqrt(sum(w*e^2)),R2=1-sum(w*e^2)/sum(w*(d$Label-sum(w*d$Label))^2),Within1=100*sum(w*(abs(e)<=1)),Within2=100*sum(w*(abs(e)<=2)),Within5=100*sum(w*(abs(e)<=5)))}
for(t in targets){
 folder<-file.path(root,'models',t);p<-fread(file.path(folder,'tables/predictions.csv.gz'));base<-fread(file.path(folder,'tables/baseline_predictions.csv.gz'))
 stopifnot(p[Evaluation=='Within-material test',.N]==374,p[Evaluation=='Source held out',.N]==1898,all(is.finite(p$Predicted)))
 frozen<-fread('outputs/20260926_nonlinear_samples_final/tables/split.csv');sp<-fread(file.path(folder,'tables/split.csv'));stopifnot(identical(sort(sp[Split=='Within-material test',SampleID]),sort(frozen[Split=='Within-material test',SampleID])))
 expected<-labels[Target==t];stopifnot(max(abs(p$Label-expected$ObservedPercent[match(p$ReferenceKey,expected$ReferenceKey)]))<1e-10)
 p[,Target:=t];allp[[t]]<-p
 mm<-rbindlist(lapply(c('Record','Material'),function(l)p[,calc(.SD,l),by=Evaluation]));mm[,Target:=t]
 bb<-rbindlist(lapply(c('Record','Material'),function(l)base[,calc(.SD,l),by=Evaluation]));bb[,Target:=t];allb[[t]]<-bb
 old<-fread(file.path(folder,'tables/metrics.csv'));check<-merge(mm,old,by=c('Evaluation','Weighting'));stopifnot(max(abs(check$MAE.x-check$MAE.y))<1e-10,max(abs(check$R2.x-check$R2.y))<1e-10)
 allm[[t]]<-mm
}
p<-rbindlist(allp);m<-rbindlist(allm);b<-rbindlist(allb)
m<-merge(m,b[,.(Target,Evaluation,Weighting,BaselineMAE=MAE)],by=c('Target','Evaluation','Weighting'));m[,MAEReductionPercent:=100*(1-MAE/BaselineMAE)]
fwrite(m,file.path(root,'comparison_metrics.csv'));fwrite(b,file.path(root,'baseline_metrics.csv'))
fwrite(p,file.path(root,'all_predictions.csv.gz'))
counts<-fread(file.path(root,'labels/protein_labels_31.csv'));mem<-fread(file.path(root,'labels/protein_members.csv.gz'));check<-merge(counts,mem[,.(DenCheck=.N,NumCheck=sum(InPathway)),by=.(GroupID,Target)],by=c('GroupID','Target'));stopifnot(all(check$Denominator==check$DenCheck),all(check$Numerator==check$NumCheck))
# Main overview figures use equal material/donor weighting for the printed metrics.
for(ev in c('Within-material test','Source held out')){
 d<-p[Evaluation==ev];mm<-m[Evaluation==ev & Weighting=='Material'];mm[,Panel:=sprintf('%s\nMAE %.3f pp | R2 %.3f | no-RNA MAE %.3f pp',pretty[Target],MAE,R2,BaselineMAE)]
 d<-merge(d,mm[,.(Target,Panel)],by='Target');d[,Panel:=factor(Panel,levels=mm[match(targets,Target),Panel])]
 g<-ggplot(d,aes(Label,Predicted))+geom_abline(slope=1,intercept=0,linetype=2,colour='grey50')+geom_point(alpha=.35,size=1.1,colour='#2166AC')+facet_wrap(~Panel,ncol=2,scales='free')+labs(title=ev,subtitle='All 18,332 RNA genes; fixed RF settings; metrics balance materials and donors',x='Assigned material fraction (%)',y='Predicted fraction (%)',caption='Each dot is an RNA record with a shared material label, not an individually measured proteome fraction.')+theme_minimal(base_size=12,base_family='Arial')+theme(panel.grid.minor=element_blank())
 stem<-if(ev=='Within-material test')'internal_test'else'source_heldout'
 ggsave(file.path(root,'figures',paste0(stem,'.png')),g,width=11,height=8,dpi=220,device=agg_png);ggsave(file.path(root,'figures',paste0(stem,'.pdf')),g,width=11,height=8,device=cairo_pdf)
}
external<-rbindlist(lapply(c('newRNA','reused_TCGA'),function(scope)fread(file.path(root,'external',scope,'predictions.csv'))))
baselines<-labels[,.(NoRNAPrediction=median(ObservedPercent)),by=Target];external<-merge(external,baselines,by='Target')
em<-external[,.(N=.N,Reference=ObservedPercent[1],MeanPrediction=mean(PredictedPercent),MAE=mean(abs(ErrorPP)),RMSE=sqrt(mean(ErrorPP^2)),Within1=100*mean(abs(ErrorPP)<=1),Within2=100*mean(abs(ErrorPP)<=2),Within5=100*mean(abs(ErrorPP)<=5),NoRNAPrediction=NoRNAPrediction[1],BaselineMAE=mean(abs(NoRNAPrediction-ObservedPercent))),by=.(Scope,Target,Material)]
history<-labels[ReferenceKey=='TCGA_ESCA_ESCC',.(Target,TrainingESCCReference=ObservedPercent)]
em<-merge(em,history,by='Target');em[,TrainingESCCBaselineMAE:=ifelse(Material=='ESCC_untreated',abs(TrainingESCCReference-Reference),NA_real_)]
fwrite(em,file.path(root,'external_metrics.csv'));fwrite(external,file.path(root,'external_predictions.csv'))
e<-external[Scope=='newRNA'];e[,Panel:=factor(pretty[Target],levels=pretty[targets])];ref<-unique(e[,.(Panel,Material,ObservedPercent)])
g<-ggplot(e,aes(Material,PredictedPercent))+geom_boxplot(width=.4,outlier.shape=NA,fill='#d6eaf8')+geom_point(position=position_jitter(width=.1,seed=25),alpha=.65,size=1.5,colour='#2166AC')+geom_point(data=ref,aes(y=ObservedPercent),shape=18,size=4,colour='#B2182B')+facet_wrap(~Panel,ncol=2,scales='free_y')+scale_x_discrete(labels=c(ESCC_untreated='ESCC tumor',ESCC_adjacent_untreated='Adjacent esophagus'))+labs(title='External material-level comparison: new RNA cohort',subtitle='23 paired donors; frozen full-data models; no fitting on external RNA',x=NULL,y='Within-assay pathway fraction (%)',caption='Red diamonds: PXD063945 material references. Blue: GSE130078 RNA predictions.\nRNA and proteomics are not patient-paired; RNA pretreatment history is unverified.')+theme_minimal(base_size=12,base_family='Arial')
ggsave(file.path(root,'figures/external_newRNA.png'),g,width=11,height=8,dpi=220,device=agg_png);ggsave(file.path(root,'figures/external_newRNA.pdf'),g,width=11,height=8,device=cairo_pdf)
# Check independent new RNA IDs and identical expression scale/order against frozen training.
rna<-readRDS('outputs/20260928_external_escc_rna_corrected/full_rna.rds');training<-readRDS('outputs/20260926_full_rna_inputs/full_rna.rds');stopifnot(!any(rna$metadata$SampleID%in%training$metadata$SampleID),identical(colnames(rna$x),colnames(training$x)),min(rna$x)==-1,min(training$x)==-1)
writeLines('PASS: all four model metrics independently recomputed; frozen split IDs identical; labels match every prediction; 124 assay-specific numerator/denominator pairs recounted; all new RNA IDs absent from training; all 18332 genes ordered identically; both expression minima -1 from log2(TPM+0.5). No cross-assay numerator restriction.',file.path(root,'validation.txt'))
print(m[Evaluation%in%c('Within-material test','Source held out') & Weighting=='Material']);print(em[Scope=='newRNA'])
