#!/usr/bin/env Rscript
# Presentation-only export for figure selection; no statistical fitting.
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(ragg)})
set.seed(25)
out<-'outputs/20261008_figure_selection_no_statistics'
plan<-fread('audit/20261008_figure_selection_no_statistics/render_plan.csv')
checks<-list()
agg_png(tempfile(fileext='.png'),width=1000,height=800)
for(i in seq_len(nrow(plan))){
 g<-readRDS(plan$RDS[i])
 old<-which(vapply(g$layers,function(l)is.data.frame(l$data)&&all(c('Control','Comparison','Q')%in%names(l$data)),logical(1)))
 if(length(old))g$layers<-g$layers[seq_len(min(old)-1)]
 before<-ggplot_build(g)
 g<-g+labs(subtitle=NULL,caption=NULL)+theme(panel.grid=element_blank())
 if(length(old)){
  ys<-g$scales$get_scales('y');lower<-ys$limits[1]
  ymax<-max(unlist(lapply(before$data,function(d)unlist(d[,intersect(c('y','ymax','yend'),names(d)),drop=FALSE]))),na.rm=TRUE)
  ys$limits<-c(lower,ymax+(ymax-lower)*.06)
 }
 after<-ggplot_build(g)
 same<-all(vapply(seq_along(before$data),function(j){
  cols<-intersect(c('x','y','xend','yend','ymin','ymax','lower','upper','middle','count'),names(before$data[[j]]))
  isTRUE(all.equal(before$data[[j]][,cols,drop=FALSE],after$data[[j]][,cols,drop=FALSE],tolerance=1e-12))
 },logical(1)))
 checks[[i]]<-data.table(Stem=plan$Stem[i],DataLayersUnchanged=same,
  ComparisonLayers=length(which(vapply(g$layers,function(l)is.data.frame(l$data)&&'Control'%in%names(l$data),logical(1)))),
  SubtitleRemoved=is.null(g$labels$subtitle),CaptionRemoved=is.null(g$labels$caption))
 for(ext in c('png','pdf'))ggsave(file.path(out,'figures',paste0(plan$Stem[i],'.',ext)),g,
  width=plan$Width[i],height=plan$Height[i],dpi=300,bg='white',device=if(ext=='png')agg_png else cairo_pdf)
 cat('EXPORTED',plan$Stem[i],'\n')
}
dev.off()
v<-rbindlist(checks);stopifnot(all(v$DataLayersUnchanged),all(v$ComparisonLayers==0),all(v$SubtitleRemoved),all(v$CaptionRemoved))
fwrite(v,file.path(out,'qa/no_statistics_validation.csv'))
