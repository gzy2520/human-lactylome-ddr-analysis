#!/usr/bin/env Rscript
# Preserve the latest figure content while changing nonnegative axis origins.
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(ragg);library(patchwork);library(uwot)})
set.seed(25)
out <- Sys.getenv('KLA_ML_FIGURE_OUTPUT', 'outputs/20261006_figures_latest/ML')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
meta <- fread('outputs/20260925_sample_fraction_inputs_final/sample_metadata.csv')
setorder(meta,ReferenceKey,SampleID)
cache <- Sys.getenv('KLA_UMAP_CACHE', '')
if (nzchar(cache)) {
 um <- fread(cache)
 stopifnot(identical(um$SampleID,meta$SampleID))
} else {
 b <- readRDS('outputs/20260926_full_rna_inputs/full_rna.rds')
 stopifnot(identical(rownames(b$x),meta$SampleID))
 v <- apply(b$x,2,var)
 pc <- prcomp(b$x[,order(v,decreasing=TRUE)[1:2000]],scale.=TRUE,rank.=10)
 set.seed(25)
 u <- uwot::umap(pc$x[,1:10],n_neighbors=20,min_dist=.3,seed=25,n_threads=1)
 um <- data.table(SampleID=meta$SampleID,Raw1=u[,1],Raw2=u[,2],
                  UMAP1=u[,1]-min(u[,1]),UMAP2=u[,2]-min(u[,2]))
 rm(b,pc);gc(FALSE)
}
stopifnot(max(abs(diff(um$Raw1)-diff(um$UMAP1)))<1e-10,
          max(abs(diff(um$Raw2)-diff(um$UMAP2)))<1e-10)
fwrite(um,file.path(out,'umap_coordinates.csv.gz'))
titles <- c(Proteome_DDR='DDR / all proteome proteins',Kla_DDR='DDR / all Kla proteins',
 Proteome_DNA_repair='DNA repair / all proteome proteins',Kla_DNA_repair='DNA repair / all Kla proteins',
 DDR='Historical: Kla-DDR / proteome DDR',DNA_repair='Historical: Kla-DNA repair / proteome DNA repair',
 Total_Kla='Historical: Kla overlap / all proteome proteins')
targets <- c('Proteome_DDR','Kla_DDR','Proteome_DNA_repair','Kla_DNA_repair')
cat_cols <- c(normal_tissue='#377EB8',cancer_tissue='#E41A1C',normal_cells='#4DAF4A',cancer_cells='#984EA3')
cat_names <- c(normal_tissue='Normal tissue',cancer_tissue='Cancer tissue',normal_cells='Normal cells',cancer_cells='Cancer cells')
split_cols <- c('Training split'='#2B5C8F','Within-material test'='#E6550D')
zx <- function(pad=.05)scale_x_continuous(limits=c(0,NA),expand=expansion(mult=c(0,pad)))
zy <- function()scale_y_continuous(limits=c(0,NA),expand=expansion(mult=c(0,.05)))
base <- theme_minimal(base_size=11,base_family='Arial')+theme(panel.grid.minor=element_blank())
if (Sys.getenv('KLA_FIGURE_NO_GRID','0') == '1')
 base <- base + theme(panel.grid=element_blank(),panel.grid.major=element_blank(),
                     axis.line=element_line(colour='grey55',linewidth=.3))
axis_audit <- list()
save <- function(g,folder,stem,w=11,h=7,zero=c('x','y')){
 dir.create(folder,recursive=TRUE,showWarnings=FALSE)
 if(inherits(g,'ggplot') && !inherits(g,'patchwork')){
  panels <- ggplot_build(g)$layout$panel_params
  for(a in zero)for(i in seq_along(panels)){
   lower <- panels[[i]][[paste0(a,'.range')]][1]
   stopifnot(is.finite(lower),abs(lower)<1e-10)
   axis_audit[[length(axis_audit)+1]] <<- data.table(Figure=file.path(folder,stem),Panel=i,Axis=a,Lower=lower)
  }
 }
 for(ext in c('png','pdf'))ggsave(file.path(folder,paste0(stem,'.',ext)),g,width=w,height=h,dpi=240,
   device=if(ext=='png')agg_png else cairo_pdf,bg='white')
 if (Sys.getenv('KLA_SAVE_PLOTS','0') == '1') saveRDS(g,file.path(folder,paste0(stem,'.rds')))
}
configs <- list()
for(t in targets)configs[[t]] <- list(folder=paste0('outputs/20260929_assay_composition/models/',t),
 dest=paste0('current_four_targets/',t),pred='tables/predictions.csv.gz',metrics='tables/metrics.csv')
for(t in c('DDR','DNA_repair'))configs[[t]] <- list(folder=paste0('outputs/20260927_ddr_fraction_models/',t),
 dest=paste0('history_cross_assay/',t),pred='tables/predictions.csv.gz',metrics='tables/metrics.csv')
configs[['Total_Kla']] <- list(folder='outputs/20260926_full_rna_models',dest='history_total_Kla',
 pred='tables/extended_predictions.csv.gz',metrics='tables/extended_metrics.csv')
for(t in names(configs)){
 c <- configs[[t]];folder <- file.path(out,c$dest)
 p <- fread(file.path(c$folder,c$pred));if('Model'%in%names(p))p<-p[Model=='RF_mtry135']
 if('Category'%in%names(p))p[,Category:=NULL]
 p[Evaluation=='Full dataset fit',Evaluation:='Training fit']
 p <- merge(p,meta[,.(SampleID,Category)],by='SampleID',all.x=TRUE)
 p[,`:=`(Error=Predicted-Label,AbsError=abs(Predicted-Label))]
 stopifnot(nrow(p)==5694,all(is.finite(p$Predicted)),all(p$Predicted>=0),all(p$Label>=0))
 m <- fread(file.path(c$folder,c$metrics));if('Model'%in%names(m))m<-m[Model=='RF_mtry135']
 m[Evaluation=='Full dataset fit',Evaluation:='Training fit']
 m <- m[Weighting=='Material'];evorder<-c('Training fit','Training split','Within-material test','Source held out')
 m[,Panel:=sprintf('%s (n = %d)\nMaterial-balanced MAE %.3f pp | R2 %.3f',Evaluation,N,MAE,R2)]
 a <- merge(p,m[,.(Evaluation,Panel)],by='Evaluation')
 a[,Panel:=factor(Panel,levels=m[match(evorder,Evaluation),Panel])]
 g <- ggplot(a,aes(Label,Predicted))+geom_abline(slope=1,intercept=0,linetype=2,colour='grey50')+
 geom_point(alpha=.35,size=1,colour='#2166AC')+facet_wrap(~Panel,ncol=2,scales='free')+zx()+zy()+
 labs(title=titles[t],x='Shared material reference (%)',y='Prediction (%)',
 caption='Each point is an RNA record with a shared material label; full-data and training fits are not validation.')+base
 save(g,folder,'01_prediction_all_evaluations',12,9)
 for(e in c('Within-material test','Source held out')){
  a<-p[Evaluation==e];a[,Category:=factor(Category,levels=names(cat_names),labels=cat_names)]
  g<-ggplot(a,aes(Label,Predicted))+geom_abline(slope=1,intercept=0,linetype=2,colour='grey50')+
  geom_point(alpha=.4,size=1,colour='#2166AC')+facet_wrap(~Category,drop=FALSE)+zx()+zy()+
  labs(title=paste(titles[t],e),x='Shared material reference (%)',y='Prediction (%)')+base
  save(g,folder,paste0('02_',gsub(' ','_',e)),12,8)
 }
 d <- copy(p[Evaluation%in%names(split_cols)])
 cm<-d[Evaluation=='Within-material test',.(N=.N,MAE=mean(AbsError)),by=Category]
 cm[,Panel:=sprintf('%s\nInternal test: n = %d | record MAE %.3f pp',cat_names[Category],N,MAE)]
 a<-merge(d,cm[,.(Category,Panel)],by='Category',all.x=TRUE)
 a[is.na(Panel),Panel:=paste0(cat_names[Category],'\nNo internal-test records')]
 panel_order<-unique(a[,.(Category,Panel)])[match(names(cat_names),Category),Panel]
 a[,Panel:=factor(Panel,levels=panel_order)]
 g<-ggplot(a,aes(Label,Predicted,colour=Evaluation))+geom_abline(slope=1,intercept=0,linetype=2,colour='grey50')+
 geom_point(alpha=.45,size=1)+scale_colour_manual(values=split_cols)+facet_wrap(~Panel,ncol=2,scales='free')+zx()+zy()+
 labs(title=paste(titles[t],'Category-stratified predictions'),x='Shared material reference (%)',y='Prediction (%)',
 caption='Blue: training split. Orange: internal test. Category metrics use test records only; labels remain material-level references.')+base
 save(g,folder,'02_category_training_and_test',12,9)
 refs <- unique(d[,.(ReferenceKey,Category,Label)]);setorder(refs,Label,ReferenceKey)
 d[,ReferenceKey:=factor(ReferenceKey,levels=refs$ReferenceKey)]
 g<-ggplot(d,aes(ReferenceKey,Predicted))+geom_boxplot(aes(fill=Category),alpha=.25,outlier.shape=NA,width=.55)+
 geom_point(aes(colour=Evaluation),position=position_jitter(width=.18,seed=25),alpha=.55,size=.9)+
 geom_point(data=refs,aes(y=Label),shape=23,size=3.2,fill='gold',colour='black')+
 scale_fill_manual(values=cat_cols,labels=cat_names,name='Material category')+
 scale_colour_manual(values=split_cols,name='Sample split')+zy()+
 labs(title=titles[t],subtitle='Diamonds: shared material reference; dots: predictions from the split-trained model',
 x='Material (ordered by reference percentage)',y='Prediction (%)')+base+
 theme(axis.text.x=element_text(angle=50,hjust=1,size=8),legend.position='top')
 save(g,folder,'03_material_distributions',14.5,7,zero='y')
 impfile<-file.path(c$folder,'tables/training_importance_annotated.csv')
 if(!file.exists(impfile))impfile<-file.path(c$folder,'tables/feature_importance_top25.csv')
 imp<-fread(impfile)[1:25];stopifnot(all(c('Ensembl','Symbol','Importance')%in%names(imp)))
 if(t!='Total_Kla'){
  raw<-fread(file.path(c$folder,'tables/training_importance.csv'))[1:25]
  stopifnot(identical(imp$Ensembl,raw$Ensembl),max(abs(imp$Importance-raw$Importance))<1e-10)
 }
 imp[,GeneLabel:=factor(paste0(Symbol,' (',Ensembl,')'),levels=rev(paste0(Symbol,' (',Ensembl,')')))]
 g<-ggplot(imp,aes(Importance,GeneLabel))+geom_col(fill='#2166AC',alpha=.88,width=.72)+
 geom_text(aes(label=sprintf('%.1f',Importance)),hjust=-.15,size=3.1)+zx(.16)+
 labs(title=paste(titles[t],'Top 25 features'),
 subtitle=if(t=='Total_Kla')'Full-data fitted model; ranked across all 18,332 input genes' else 'Training-split model; ranked across all 18,332 input genes',
 x='Regression importance (variance reduction)',y='Gene: Symbol (Ensembl ID)',
 caption='Symbols are display annotations; stable Ensembl IDs remain the analysis keys.')+base+
 theme(panel.grid.major.y=element_blank(),axis.text.y=element_text(size=9))
 save(g,folder,'04_top25',10,9,zero='x')
 a<-merge(um,p[Evaluation=='Training fit',.(SampleID,Category,Label,Predicted)],by='SampleID')
 g1<-ggplot(a,aes(UMAP1,UMAP2,colour=Category))+geom_point(size=1,alpha=.65)+
 scale_colour_manual(values=cat_cols,labels=cat_names,name='Category')+zx()+zy()+
 labs(title='A. Material category',x='UMAP 1 (translated)',y='UMAP 2 (translated)')+base+theme(legend.position='bottom')
 limits<-range(c(a$Label,a$Predicted));mid<-median(a$Label)
 g2<-ggplot(a,aes(UMAP1,UMAP2,colour=Label))+geom_point(size=1,alpha=.7)+
 scale_colour_gradient2(low='#2166AC',mid='#FFFFBF',high='#B2182B',midpoint=mid,limits=limits,name='Reference %')+
 zx()+zy()+labs(title='B. Shared material reference',x='UMAP 1 (translated)',y='UMAP 2 (translated)')+base+theme(legend.position='bottom')
 g3<-ggplot(a,aes(UMAP1,UMAP2,colour=Predicted))+geom_point(size=1,alpha=.7)+
 scale_colour_gradient2(low='#2166AC',mid='#FFFFBF',high='#B2182B',midpoint=mid,limits=limits,name='Fitted %')+
 zx()+zy()+labs(title='C. Full-data fitted values',x='UMAP 1 (translated)',y='UMAP 2 (translated)')+base+theme(legend.position='bottom')
 for(gp in list(g1,g2,g3)){
  pp<-ggplot_build(gp)$layout$panel_params[[1]]
  stopifnot(pp$x.range[1]==0,pp$y.range[1]==0)
 }
 g<-(g1|g2|g3)+plot_annotation(title=titles[t],
 subtitle='UMAP from top 2,000 variable genes / 10 PCs; axes translated to zero minima. Full-data fit is descriptive, not validation.')
 save(g,folder,'05_UMAP_zero_origin',16,6,zero=character())
 g<-ggplot(d,aes(Label,AbsError,colour=Evaluation))+geom_point(alpha=.45,size=1)+
 scale_colour_manual(values=split_cols)+zx()+zy()+labs(title=titles[t],x='Shared material reference (%)',y='Absolute error (percentage points)')+base
 save(g,folder,'06_absolute_error_zero_origin')
 g<-ggplot(d,aes(AbsError,fill=Evaluation))+geom_histogram(bins=40,boundary=0,alpha=.55,position='identity')+
 scale_fill_manual(values=split_cols)+zx()+zy()+labs(title=titles[t],x='Absolute error (percentage points)',y='RNA records')+base
 stopifnot(sum(ggplot_build(g)$data[[1]]$count)==nrow(d))
 save(g,folder,'07_absolute_error_distribution_zero_origin')
 g<-ggplot(d,aes(Label,Error,colour=Evaluation))+geom_hline(yintercept=0)+geom_point(alpha=.45,size=1)+
 scale_colour_manual(values=split_cols)+zx()+labs(title=paste(titles[t],'Signed residuals'),
 x='Shared material reference (%)',y='Prediction - reference (percentage points)',caption='Signed errors retain negative values to show underprediction.')+base
 save(g,file.path(folder,'signed_diagnostic_appendix'),'signed_residual',zero='x')
 fwrite(p[,.(SampleID,Evaluation,Label,Predicted,Error,AbsError)],file.path(folder,'plot_data.csv.gz'))
 cat('Completed',t,'\n')
}
# Latest four-target overview figures were missing from the previous package.
p<-fread('outputs/20260929_assay_composition/all_predictions.csv.gz')
m<-fread('outputs/20260929_assay_composition/comparison_metrics.csv')
for(e in c('Within-material test','Source held out')){
 mm<-copy(m[Evaluation==e & Weighting=='Material']);mm<-mm[match(targets,Target)]
 mm[,Panel:=sprintf('%s\nMAE %.3f pp | R2 %.3f | no-RNA MAE %.3f pp',titles[Target],MAE,R2,BaselineMAE)]
 d<-merge(p[Evaluation==e],mm[,.(Target,Panel)],by='Target');d[,Panel:=factor(Panel,levels=mm$Panel)]
 g<-ggplot(d,aes(Label,Predicted))+geom_abline(slope=1,intercept=0,linetype=2,colour='grey50')+
 geom_point(alpha=.35,size=1.1,colour='#2166AC')+facet_wrap(~Panel,ncol=2,scales='free')+zx()+zy()+
 labs(title=e,subtitle='All 18,332 RNA genes; fixed RF settings; metrics balance materials and donors',
 x='Shared material reference (%)',y='Prediction (%)',caption='RNA records share material labels; these are not individual proteomic measurements.')+base
 save(g,file.path(out,'current_four_targets/overview'),if(e=='Within-material test')'internal_test'else'source_heldout',12,9)
}
for(scope in c('newRNA','reused_TCGA')){
 p<-fread(file.path('outputs/20260929_assay_composition/external',scope,'predictions.csv'))
 p[,Panel:=factor(titles[Target],levels=titles[targets])];ref<-unique(p[,.(Panel,Material,ObservedPercent)])
 g<-ggplot(p,aes(Material,PredictedPercent))+geom_boxplot(width=.4,outlier.shape=NA,fill='#d6eaf8')+
 geom_point(position=position_jitter(width=.1,seed=25),alpha=.65,size=1.3,colour='#2166AC')+
 geom_point(data=ref,aes(y=ObservedPercent),shape=18,size=4,colour='#B2182B')+facet_wrap(~Panel,ncol=2,scales='free_y')+
 scale_x_discrete(labels=c(ESCC_untreated='ESCC tumor',ESCC_adjacent_untreated='Adjacent esophagus'))+zy()+
 labs(title=paste('External material comparison:',scope),x=NULL,y='Within-assay pathway fraction (%)',
 caption=if(scope=='newRNA')'Red: PXD063945 material references. Blue: GSE130078 RNA predictions.\n23 paired RNA donors; no patient pairing to proteomics; RNA pretreatment history is unverified.' else
 'Previously used TCGA RNA compared with new proteomic references; this is not an independent external RNA test.')+base
 save(g,file.path(out,'current_four_targets/external',scope),'external_comparison_zero_origin',12,9,zero='y')
}
p<-fread('outputs/20260928_report_revision/external/predictions.csv')
ref<-unique(p[,.(Target,Material,ObservedPercent)])
g<-ggplot(p,aes(Material,PredictedPercent))+geom_boxplot(outlier.shape=NA,fill='#d6eaf8')+
 geom_point(position=position_jitter(width=.1,seed=25),alpha=.6,size=1)+
 geom_point(data=ref,aes(y=ObservedPercent),shape=18,size=3,colour='#B2182B')+facet_wrap(~Target,scales='free_y')+zy()+
 labs(title='Historical cross-assay targets: corrected external RNA scale',x='Material',y='Historical cross-assay fraction (%)',
 caption='Corrected log2(TPM+0.5) RNA; not the current Figure 1a target; cohorts are not patient-paired.')+base
save(g,file.path(out,'history_cross_assay/external_corrected'),'external_comparison_zero_origin',12,8,zero='y')
fwrite(rbindlist(axis_audit),file.path(out,'axis_lower_bound_audit.csv'))
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
