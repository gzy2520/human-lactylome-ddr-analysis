#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(clubSandwich);library(ragg)})
set.seed(25)
out<-'outputs/20261008_all_category_statistics'
cats<-c('normal_tissue','cancer_tissue','normal_cells','cancer_cells')
frozen<-'corrected_final_result_20260905_sample_only/Corrected_Data/sample_only_inputs'
objects<-'outputs/20261006_teacher_figures/plot_objects'
units<-list()
fit<-function(d,family,panel,direction=''){
 d<-copy(d);d[,Category:=factor(Category,levels=cats)]
 stopifnot(all(is.finite(d$Value)),!anyNA(d$Source),all(cats%in%d$Category))
 d[,Weight:=1/.N,by=.(Source,Category)]
 units[[length(units)+1]]<<-cbind(data.table(Family=family,Panel=panel,Direction=direction),d)
 m<-lm(Value~0+Category,data=d,weights=Weight)
 V<-vcovCR(m,cluster=d$Source,type='CR2',inverse_var=FALSE)
 rbindlist(lapply(cats[-2],function(cmp){
  C<-matrix(0,1,length(coef(m)));C[1,match('Categorycancer_tissue',names(coef(m)))]<- -1
  C[1,match(paste0('Category',cmp),names(coef(m)))]<-1
  z<-as.data.table(linear_contrast(m,V,C,p_values=TRUE))
  data.table(Family=family,Panel=panel,Direction=direction,Control='cancer_tissue',Comparison=cmp,
   Estimate=z$Est,SE=z$SE,DF=z$df,LowerCI=z$CI_L,UpperCI=z$CI_U,P=z$p_val,
   ControlRows=sum(d$Category=='cancer_tissue'),ComparisonRows=sum(d$Category==cmp),
   ControlSources=uniqueN(d[Category=='cancer_tissue',Source]),ComparisonSources=uniqueN(d[Category==cmp,Source]))
 }))
}
get_source<-function(x)sub('^.*?(PXD[0-9]+).*','\\1',x,perl=TRUE)
f<-fread(file.path(frozen,'figure1_sample_boxplot_values.csv'))
f[,Source:=get_source(SourceFile)]
keys<-c('Dataset','SourceFile','SampleID','Category')
stopifnot(f[,all(uniqueN(DdrFractionPercentage)==1),by=keys]$V1)
fu<-unique(f,by=keys)
r<-rbindlist(lapply(unique(fu$Dataset),function(ds)
 fit(fu[Dataset==ds,.(Source,Category,Unit=paste(SourceFile,SampleID,Category),Value=DdrFractionPercentage)],'ProteinDDR',ds)))
p<-fread(file.path(frozen,'figure1_pathway_summary_sample_boxplot_values.csv'))
p[,Source:=get_source(SourceFile)]
stopifnot(nrow(p)==504,!anyDuplicated(p,by=c('SourceFile','SampleID','Category','Pathway')))
for(pa in c('BER','NER','MMR','FA','HR','AEJ','NHEJ'))for(di in c('Pro','Inh'))
 r<-rbind(r,fit(p[Pathway==pa,.(Source,Category,Unit=paste(SourceFile,SampleID,Category),
             Value=100*get(if(di=='Pro')'PositiveFraction'else'NegativeFraction'))],'KlaPathways',pa,di))
rat<-fread(file.path(frozen,'figure1_mki67_ratio_sample_values.csv'))[ObservationType=='sample']
rat[,Source:=get_source(SourceFile)]
rk<-c('Denominator','SourceFile','SampleID','Category')
stopifnot(rat[,all(uniqueN(Ratio)==1),by=rk]$V1)
rat<-unique(rat,by=rk)
for(den in c('H3C1','ACTB','TUBB'))r<-rbind(r,fit(rat[Denominator==den,
 .(Source,Category,Unit=paste(SourceFile,SampleID,Category),Value=log10(Ratio))],'MKI67Ratios',den))
rna<-fread('outputs/20260919_repair_release/core/figure1a_fixed_ddr_28materials.csv')
src<-unique(fread('outputs/20260925_sample_fraction_inputs_final/sample_metadata.csv')[,.(ReferenceKey,ConnGroup)])
rna<-merge(rna,src,by='ReferenceKey');stopifnot(nrow(rna)==28,!anyNA(rna$ConnGroup))
r<-rbind(r,fit(rna[,.(Source=as.character(ConnGroup),Category,Unit=ReferenceKey,Value=DdrFractionMedian)],'RNADDRFraction','RNA DDR fraction'))
r[,Status:=ifelse(is.finite(P)&is.finite(SE)&SE>1e-10,'estimable','zero variance / not estimable')]
r[Status!='estimable',P:=NA_real_]
r[,Q:=p.adjust(P,'BH',n=.N),by=Family]
r[,Label:=ifelse(is.na(Q),'NE',ifelse(Q<.0001,'****',ifelse(Q<.001,'***',ifelse(Q<.01,'**',ifelse(Q<.05,'*','ns')))))]
fwrite(r,file.path(out,'statistics/tumor_reference_comparisons.csv'))
fwrite(rbindlist(units),file.path(out,'statistics/inference_units.csv'))
# Reuse the exact already-reviewed RNA expression comparison family (six tests).
re<-fread('audit/20261008_rna_figure2_tumor_reference/tables/RNA_expression_pairwise_source_clustered.csv')
fwrite(re,file.path(out,'statistics/RNA_expression_existing_comparisons.csv'))

checks<-list();annotations<-list()
grid_lines<-function(g){
 n<-as.integer(grepl('panel.grid',g$name,fixed=TRUE)&&any(c('polyline','segments')%in%class(g)))
 if(length(g$grobs))n<-n+sum(vapply(g$grobs,grid_lines,integer(1)))
 if(length(g$children))n<-n+sum(vapply(g$children,grid_lines,integer(1)))
 n
}
agg_png(tempfile(fileext='.png'),width=1000,height=800)
render<-function(rds,st,stem,group_field=NULL,offsets=NULL){
 g<-readRDS(file.path(objects,rds));st<-copy(st)
 # Remove prior comparison annotations only, preserving all data and summary layers.
 old<-which(vapply(g$layers,function(l)is.data.frame(l$data)&&all(c('Control','Comparison','Q')%in%names(l$data)),logical(1)))
 if(length(old))g$layers<-g$layers[seq_len(min(old)-1)]
 base<-ggplot_build(g);k<-length(g$layers)
 if(is.null(group_field)){st[,Offset:=0];st[,Text:=Label]}else{
  st[,Offset:=unname(offsets[get(group_field)])]
  st[,Text:=paste(ifelse(get(group_field)=='Whole proteome','Proteome',
                ifelse(get(group_field)=='Lactylome (Kla)','Kla',get(group_field))),Label)]
 }
 setorder(st,Comparison,Offset)
 ys<-g$scales$get_scales('y');trans<-ys$trans
 ymax<-max(unlist(lapply(base$data,function(d)unlist(d[,intersect(c('y','ymax','yend'),names(d)),drop=FALSE]))),na.rm=TRUE)
 ymin<-ys$limits[1];step<-(ymax-ymin)*.075
 st[,`:=`(Left=pmin(match(Control,cats),match(Comparison,cats))+Offset,
           Right=pmax(match(Control,cats),match(Comparison,cats))+Offset,
           Y=ymax+seq_len(.N)*step)]
 tips<-rbind(st[,.(X=Left,Y,Yend=Y-step*.23)],st[,.(X=Right,Y,Yend=Y-step*.23)])
 st[,`:=`(PlotY=trans$inverse(Y),TextY=trans$inverse(Y+step*.08))]
 tips[,`:=`(PlotY=trans$inverse(Y),PlotYend=trans$inverse(Yend))]
 g<-g+geom_segment(data=st,aes(x=Left,xend=Right,y=PlotY,yend=PlotY),inherit.aes=FALSE,linewidth=.35)+
  geom_segment(data=tips,aes(x=X,xend=X,y=PlotY,yend=PlotYend),inherit.aes=FALSE,linewidth=.35)+
  geom_text(data=st,aes(x=(Left+Right)/2,y=TextY,label=Text),inherit.aes=FALSE,size=3.5,
            family='Arial Unicode MS',vjust=0)+
  labs(subtitle=NULL,caption='Tumor-tissue reference; source-balanced CR2 contrasts, BH q.\nns: q >= 0.05; NE: contrast not estimable.')+
  theme(panel.grid=element_blank(),plot.caption=element_text(size=9,hjust=0),plot.margin=margin(8,12,8,8))
 ys<-g$scales$get_scales('y');ys$limits<-c(ymin,max(st$Y)+step*.85)
 built<-ggplot_build(g)
 same<-all(vapply(seq_len(k),function(i){
  a<-base$data[[i]];b<-built$data[[i]]
  cols<-intersect(c('x','y','xend','yend','ymin','ymax','lower','upper','middle','count'),names(a))
  isTRUE(all.equal(a[,cols,drop=FALSE],b[,cols,drop=FALSE],tolerance=1e-12))
 },logical(1)))
 checks[[length(checks)+1]]<<-data.table(Figure=stem,OriginalLayers=k,ValuesUnchanged=same,
  Brackets=nrow(st),LabelsMatch=identical(built$data[[k+3]]$label,st$Text),
  TumorReference=all(st$Control=='cancer_tissue'),GridLines=grid_lines(ggplotGrob(g)))
 st[,Figure:=stem];annotations[[length(annotations)+1]]<<-st
 for(ext in c('png','pdf'))ggsave(file.path(out,'figures',paste0(stem,'.',ext)),g,
  width=8.8,height=7.4,dpi=300,bg='white',device=if(ext=='png')agg_png else cairo_pdf)
 cat('RENDERED',stem,'\n')
}
render('fraction__Figure_1_DDR_fraction_candidate_category_boxplot_refined.rds',r[Family=='ProteinDDR'],
 'Figure_1a_DDR_fraction_boxplot','Panel',c('Whole proteome'=-.165,'Lactylome (Kla)'=.165))
for(pa in c('BER','NER','MMR','FA','HR','AEJ','NHEJ'))render(
 paste0('pathways__Figure_2_DDR_pathway_summary_',pa,'_barplot.rds'),r[Family=='KlaPathways'&Panel==pa],
 paste0('Kla_',pa,'_barplot'),'Direction',c(Pro=-.18,Inh=.18))
for(den in c('H3C1','ACTB','TUBB'))render(paste0('ratios__Figure_1_MKI67_over_',den,'_boxplot.rds'),
 r[Family=='MKI67Ratios'&Panel==den],paste0('MKI67_over_',den,'_boxplot'))
render('core__Figure_1a_RNA_DDR_expression_28materials_boxplot.rds',re[GeneSet=='DDR genes'],
 'RNA_DDR_expression_28materials_boxplot')
render('core__Figure_1a_RNA_DDR_fraction_28materials_boxplot.rds',r[Family=='RNADDRFraction'],
 'RNA_DDR_fraction_28materials_boxplot')
dev.off()
v<-rbindlist(checks);fwrite(v,file.path(out,'qa/render_validation.csv'))
stopifnot(nrow(v)==13,all(v$ValuesUnchanged),all(v$LabelsMatch),all(v$TumorReference),all(v$GridLines==0))
fwrite(rbindlist(annotations,fill=TRUE),file.path(out,'statistics/displayed_annotations.csv'))
writeLines(trimws(capture.output(sessionInfo()),which='right'),file.path(out,'sessionInfo.txt'))
print(r[,.(Tests=.N,Estimable=sum(Status=='estimable'),Significant=sum(Q<.05,na.rm=TRUE)),by=Family])
