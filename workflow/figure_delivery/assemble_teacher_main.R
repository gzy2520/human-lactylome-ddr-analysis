#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(patchwork);library(ragg)})
set.seed(25)
out<-Sys.getenv('KLA_TEACHER_OUTPUT','outputs/20261006_teacher_figures')
main<-file.path(out,'main_figures');dir.create(main,recursive=TRUE,showWarnings=FALSE)
read_plot<-function(folder,stem)readRDS(file.path(out,'plot_objects',paste0(folder,'__',stem,'.rds')))
save<-function(g,stem,w,h){
 for(ext in c('png','pdf'))ggsave(file.path(main,paste0(stem,'.',ext)),g,width=w,height=h,
  dpi=300,device=if(ext=='png')agg_png else cairo_pdf,bg='white',limitsize=FALSE)
 saveRDS(g,file.path(main,paste0(stem,'.rds')))
}
f1a<-read_plot('fraction','Figure_1_DDR_fraction_candidate_category_boxplot_refined')
f1b<-read_plot('core','Figure_1b_RNA_proliferation_hallmark_28materials_boxplot')
f1<-(f1a|f1b)+plot_annotation(tag_levels='a')&theme(plot.tag=element_text(size=24,face='bold'))
save(f1,'Figure_1_proteomic_DDR_and_RNA_proliferation',17,8.5)
up1<-read_plot('publication','Figure_2b_Kla_DDR_UpSet')
up2<-read_plot('publication','Figure_2a_whole_proteome_DDR_UpSet')
paths<-lapply(c('BER','NER','MMR','FA','HR','AEJ','NHEJ'),function(pa){
 read_plot('pathways',paste0('Figure_2_DDR_pathway_summary_',pa,'_barplot'))+
  labs(title=paste0('Kla-',ifelse(pa=='AEJ','A-EJ',pa)),caption=NULL)+
  scale_x_discrete(labels=c('non-tumor tissues'='non-tumor\ntissues','tumor tissues'='tumor\ntissues',
                           'normal cell lines'='normal\ncell lines','cancer cell lines'='cancer\ncell lines'))+
  theme(axis.text.x=element_text(size=10),axis.text.y=element_text(size=10),
                           axis.title.y=element_text(size=11),plot.title=element_text(size=14))
})
# The quantitative panels keep their reviewed order. The separate schematic
# in the WeChat montage has no native source in this repository.
f2<-wrap_plots(c(list(wrap_elements(full=patchworkGrob(up1)),wrap_elements(full=patchworkGrob(up2))),paths),
 design='AAAABBBB\nCCDDEEFF\nGGHHII##',heights=c(1.25,1,1))+plot_annotation(tag_levels='a',
 caption='Prespecified source-balanced contrasts; study-clustered CR2, BH correction across 56 pathway tests.\nns: q >= 0.05. Pro/Inh are tested separately; zero-variance contrasts have no bracket.')&
 theme(plot.tag=element_text(size=23,face='bold'),plot.caption=element_text(size=10,hjust=0))
save(f2,'Figure_2_DDR_sets_and_Kla_pathways',18,14)

targets<-c('Proteome_DDR','Kla_DDR','Proteome_DNA_repair','Kla_DNA_repair')
names<-c(Proteome_DDR='Proteome DDR',Kla_DDR='Lactylome DDR',
         Proteome_DNA_repair='Proteome DNA repair',Kla_DNA_repair='Lactylome DNA repair')
p<-fread('outputs/20260929_assay_composition/all_predictions.csv.gz')
m<-fread('outputs/20260929_assay_composition/comparison_metrics.csv')[Weighting=='Material']
ex<-fread('outputs/20260929_assay_composition/external/newRNA/predictions.csv')
em<-fread('outputs/20260929_assay_composition/external_metrics.csv')[Scope=='newRNA']
base<-theme_classic(base_family='Arial',base_size=11)+theme(panel.grid=element_blank(),
 plot.title=element_text(size=12,face='bold'),plot.subtitle=element_text(size=10),
 plot.margin=margin(8,10,8,8))
zx<-function()scale_x_continuous(limits=c(0,NA),expand=expansion(mult=c(0,.06)))
zy<-function()scale_y_continuous(limits=c(0,NA),expand=expansion(mult=c(0,.06)))
figs<-list();table<-list()
for(t in targets){
 for(e in c('Within-material test','Source held out')){
  d<-p[Target==t&Evaluation==e];mm<-m[Target==t&Evaluation==e];stopifnot(nrow(mm)==1)
  g<-ggplot(d,aes(Label,Predicted))+geom_abline(slope=1,intercept=0,linetype=2,colour='grey50')+
   geom_point(alpha=.35,size=.9,colour='#2166AC')+zx()+zy()+labs(title=names[t],
    subtitle=sprintf('%s | MAE %.3f pp | R2 %.3f',ifelse(e=='Source held out','Source held out','Internal test'),mm$MAE,mm$R2),
    x='Shared material reference (%)',y='RNA prediction (%)')+base
  figs[[length(figs)+1]]<-g
 }
 d<-ex[Target==t];ref<-unique(d[,.(Material,ObservedPercent)])
 ma<-em[Target==t,mean(MAE)]
 g<-ggplot(d,aes(Material,PredictedPercent))+geom_boxplot(width=.42,outlier.shape=NA,fill='#d6eaf8')+
  geom_point(position=position_jitter(width=.09,seed=25),alpha=.6,size=1,colour='#2166AC')+
  geom_point(data=ref,aes(y=ObservedPercent),shape=18,size=3.5,colour='#B2182B')+zy()+
  scale_x_discrete(labels=c(ESCC_untreated='ESCC tumor',ESCC_adjacent_untreated='Adjacent esophagus'))+
  labs(title=names[t],subtitle=sprintf('External material comparison | MAE %.3f pp',ma),x=NULL,
       y='Within-assay fraction (%)')+base+theme(axis.text.x=element_text(size=10))
 figs[[length(figs)+1]]<-g
 mi<-m[Target==t&Evaluation=='Within-material test'];ms<-m[Target==t&Evaluation=='Source held out']
 table[[length(table)+1]]<-data.table(Target=t,TargetLabel=names[t],InternalMAE=mi$MAE,InternalR2=mi$R2,
 SourceMAE=ms$MAE,SourceR2=ms$R2,SourceBaselineMAE=ms$BaselineMAE,ExternalMaterialMAE=ma,
 ExternalTumorMAE=em[Target==t&Material=='ESCC_untreated',MAE],
 ExternalAdjacentMAE=em[Target==t&Material=='ESCC_adjacent_untreated',MAE])
}
f3<-wrap_plots(figs,ncol=3)+plot_annotation(tag_levels='a',
 title='RNA prediction of within-assay DDR and DNA repair fractions',
 caption='All 18,332 input RNA genes; fixed random-forest settings. Metrics balance materials and donors.\nInternal test: 374 RNA records. Source holdout: 14 connected source blocks, 1,898 records.\nExternal comparison: 46 RNA records / 23 paired RNA donors; two material references, not patient-matched proteomics. Red diamonds: references.')&
 theme(plot.tag=element_text(size=18,face='bold'),plot.caption=element_text(size=10,hjust=0),
       plot.title=element_text(size=12,face='bold'))
save(f3,'Figure_3_RNA_machine_learning_four_targets',16,14)
fwrite(rbindlist(table),file.path(main,'Table_1_ML_performance.csv'))
