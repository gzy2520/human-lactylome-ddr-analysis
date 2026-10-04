#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(ragg);library(uwot)})
set.seed(25)
out<-'outputs/20261004_all_figures/ML';dir.create(out,recursive=TRUE,showWarnings=FALSE)
meta<-fread('outputs/20260925_sample_fraction_inputs_final/sample_metadata.csv')
b<-readRDS('outputs/20260926_full_rna_inputs/full_rna.rds');v<-apply(b$x,2,var);pc<-prcomp(b$x[,order(v,decreasing=TRUE)[1:2000]],scale.=TRUE,rank.=10)
set.seed(25);u<-uwot::umap(pc$x[,1:10],n_neighbors=20,min_dist=.3,n_threads=1);u<-sweep(u,2,apply(u,2,min),'-')
um<-data.table(SampleID=rownames(b$x),UMAP1=u[,1],UMAP2=u[,2]);rm(b,pc);gc(FALSE)
configs<-list()
for(t in c('Proteome_DDR','Kla_DDR','Proteome_DNA_repair','Kla_DNA_repair'))configs[[paste0('01_current_four_targets/',t)]]<-list(folder=paste0('outputs/20260929_assay_composition/models/',t),pred='tables/predictions.csv.gz')
for(t in c('DDR','DNA_repair'))configs[[paste0('02_previous_cross_assay/',t)]]<-list(folder=paste0('outputs/20260927_ddr_fraction_models/',t),pred='tables/predictions.csv.gz')
configs[['03_previous_total_Kla_fraction']]<-list(folder='outputs/20260926_full_rna_models',pred='tables/extended_predictions.csv.gz')
axes<-list();save<-function(g,folder,stem,w=11,h=7){for(ext in c('png','pdf'))ggsave(file.path(folder,paste0(stem,'.',ext)),g,width=w,height=h,dpi=240,device=if(ext=='png')agg_png else cairo_pdf,bg='white')}
zscale<-function()scale_x_continuous(limits=c(0,NA),expand=expansion(mult=c(0,.05)))
yscale<-function()scale_y_continuous(limits=c(0,NA),expand=expansion(mult=c(0,.05)))
base<-theme_minimal(base_size=11,base_family='Arial')+theme(panel.grid.minor=element_blank())
for(n in names(configs)){
 c<-configs[[n]];folder<-file.path(out,n);dir.create(folder,recursive=TRUE,showWarnings=FALSE)
 p<-fread(file.path(c$folder,c$pred));if('Model'%in%names(p))p<-p[Model=='RF_mtry135'];if(!'Label'%in%names(p))stop('Inspect prediction schema')
 if('Category'%in%names(p))p[,Category:=NULL]
 p[Evaluation=='Full dataset fit',Evaluation:='Training fit']
 p<-merge(p,meta[,.(SampleID,Category)],by='SampleID',all.x=TRUE)
 p[,Error:=Predicted-Label];p[,AbsError:=abs(Error)]
 g<-ggplot(p,aes(Label,Predicted))+geom_abline(slope=1,intercept=0,linetype=2,colour='grey50')+geom_point(alpha=.35,size=1,colour='#2166AC')+facet_wrap(~Evaluation,ncol=2)+zscale()+yscale()+labs(title=n,x='Shared material reference (%)',y='Prediction (%)',caption='Each point is an RNA record; reference labels are shared by material.')+base
 save(g,folder,'01_prediction_all_evaluations',12,9)
 for(e in intersect(c('Within-material test','Source held out'),unique(p$Evaluation))){d<-p[Evaluation==e];g<-ggplot(d,aes(Label,Predicted))+geom_abline(slope=1,intercept=0,linetype=2,colour='grey50')+geom_point(alpha=.4,size=1,colour='#2166AC')+facet_wrap(~Category)+zscale()+yscale()+labs(title=paste(n,e),x='Shared material reference (%)',y='Prediction (%)')+base;save(g,folder,paste0('02_',gsub(' ','_',e)))}
 d<-p[Evaluation%in%c('Training split','Within-material test')];if(nrow(d)==0)d<-p[Evaluation!='Source held out']
 g<-ggplot(d,aes(ReferenceKey,Predicted,colour=Evaluation))+geom_boxplot(outlier.shape=NA)+geom_point(position=position_jitter(width=.12,seed=25),alpha=.35,size=.7)+yscale()+labs(title=n,x='Material',y='Prediction (%)')+base+theme(axis.text.x=element_text(angle=55,hjust=1,size=7));save(g,folder,'03_material_distributions',14,7)
 impfile<-file.path(c$folder,'tables/training_importance.csv');if(!file.exists(impfile))impfile<-file.path(c$folder,'tables/feature_importance_top25.csv');if(file.exists(impfile)){imp<-fread(impfile)[1:25];g<-ggplot(imp,aes(Importance,reorder(Ensembl,Importance)))+geom_col(fill='#2166AC')+zscale()+labs(title=paste(n,'Top 25 features'),x='Regression importance (variance reduction)',y='Ensembl ID')+base;save(g,folder,'04_top25',10,9)}
 a<-merge(um,p[Evaluation=='Training fit',.(SampleID,Label,Predicted)],by='SampleID');if(nrow(a)){a<-melt(a,id.vars=c('SampleID','UMAP1','UMAP2'),measure.vars=c('Label','Predicted'));g<-ggplot(a,aes(UMAP1,UMAP2,colour=value))+geom_point(size=.9,alpha=.7)+facet_wrap(~variable)+scale_colour_gradient(low='#2166AC',high='#B2182B')+zscale()+yscale()+labs(title=paste(n,'UMAP: shared labels and full-data fit'),x='UMAP 1 (translated minimum = 0)',y='UMAP 2 (translated minimum = 0)',colour='Percent',caption='Coordinates translated only; distances unchanged. Full-data fitted values are descriptive, not validation.')+base;save(g,folder,'05_UMAP_zero_origin',12,6)}
 g<-ggplot(d,aes(Label,AbsError,colour=Evaluation))+geom_point(alpha=.4,size=1)+zscale()+yscale()+labs(title=paste(n,'Absolute prediction error'),x='Shared material reference (%)',y='Absolute error (percentage points)')+base;save(g,folder,'06_absolute_error_zero_origin')
 g<-ggplot(d,aes(AbsError,fill=Evaluation))+geom_histogram(bins=40,boundary=0,alpha=.55,position='identity')+zscale()+yscale()+labs(title=paste(n,'Absolute error distribution'),x='Absolute error (percentage points)',y='RNA records')+base;save(g,folder,'07_absolute_error_distribution_zero_origin')
 # Signed residuals must preserve negative errors. They are placed in a separately labelled diagnostic appendix.
 diag<-file.path(folder,'signed_diagnostic_appendix');dir.create(diag)
 g<-ggplot(d,aes(Label,Error,colour=Evaluation))+geom_hline(yintercept=0)+geom_point(alpha=.4,size=1)+zscale()+labs(title=paste(n,'Signed residuals'),x='Shared material reference (%)',y='Prediction - reference (percentage points)',caption='Negative residuals retained to show underprediction; main error figures start at zero.')+base;save(g,diag,'signed_residual')
 fwrite(p[,.(SampleID,Evaluation,Label,Predicted,Error,AbsError)],file.path(folder,'plot_data.csv.gz'))
 cat('Completed',n,'\n')
}
# External comparison uses only corrected +0.5 new RNA and explicitly labelled reused RNA.
for(scope in c('newRNA','reused_TCGA')){
 p<-fread(file.path('outputs/20260929_assay_composition/external',scope,'predictions.csv'));ref<-unique(p[,.(Target,Material,ObservedPercent)])
 folder<-file.path(out,'01_current_four_targets','external',scope);dir.create(folder,recursive=TRUE)
 g<-ggplot(p,aes(Material,PredictedPercent))+geom_boxplot(outlier.shape=NA,fill='#d6eaf8')+geom_point(position=position_jitter(width=.1,seed=25),alpha=.6,size=1)+geom_point(data=ref,aes(y=ObservedPercent),colour='#B2182B',shape=18,size=3)+facet_wrap(~Target,scales='free_y')+yscale()+labs(title=paste('External material comparison:',scope),x='Material',y='Within-assay fraction (%)',caption='Red: proteome material reference; points: RNA predictions. Cohorts are not patient-paired.')+base;save(g,folder,'external_comparison_zero_origin',12,8)
}
writeLines('Nonnegative numeric axes start at zero with no lower padding. Categorical axes retain labels. UMAP is translated to zero minima with distances preserved. Signed residual appendix retains negative values; absolute-error main figures start at zero. Models and predictions were not refitted.',file.path(out,'axis_policy.txt'))

# Corrected prior cross-assay external predictions; historical target, retained for comparison.
f<-'outputs/20260928_report_revision/external/predictions.csv'
if(file.exists(f)){
 p<-fread(f);ref<-unique(p[,.(Target,Material,ObservedPercent)])
 folder<-file.path(out,'02_previous_cross_assay/external_corrected');dir.create(folder,recursive=TRUE,showWarnings=FALSE)
 g<-ggplot(p,aes(Material,PredictedPercent))+geom_boxplot(outlier.shape=NA,fill='#d6eaf8')+geom_point(position=position_jitter(width=.1,seed=25),alpha=.6,size=1)+geom_point(data=ref,aes(y=ObservedPercent),colour='#B2182B',shape=18,size=3)+facet_wrap(~Target,scales='free_y')+yscale()+labs(title='Previous cross-assay targets: corrected external RNA scale',x='Material',y='Kla / detected pathway proteins (%)')+base;save(g,folder,'external_comparison_zero_origin',12,7)
}
