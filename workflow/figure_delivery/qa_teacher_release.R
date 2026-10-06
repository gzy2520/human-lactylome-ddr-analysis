#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(patchwork)})
out<-'outputs/20261006_teacher_figures';dir.create(file.path(out,'qa'),showWarnings=FALSE)
checks<-list();put<-function(test,pass,detail='')checks[[length(checks)+1]]<<-data.table(Test=test,Pass=pass,Detail=detail)
same_table<-function(a,b){
 x<-fread(a);y<-fread(b)
 isTRUE(all.equal(as.data.frame(x),as.data.frame(y),tolerance=1e-12,check.attributes=FALSE))
}
for(group in c('core','reference')){
 old<-file.path('outputs/20260919_repair_release',group)
 for(f in list.files(old,pattern='\\.(csv|csv.gz)$',full.names=TRUE)){
  new<-file.path(out,'RNA',group,basename(f))
  if(file.exists(new)&&!basename(f)%in%c('render_input_fingerprints.csv','figure_manifest.csv'))
   put(paste('RNA unchanged',group,basename(f)),same_table(f,new))
 }
}
for(f in list.files('outputs/20260920_lactate_metabolism/tables',pattern='\\.(csv|csv.gz)$',full.names=TRUE)){
 new<-file.path(out,'lactate/tables',basename(f))
 if(file.exists(new))put(paste('Lactate unchanged',basename(f)),same_table(f,new))
}
for(f in list.files(file.path(out,'ML'),pattern='plot_data.csv.gz',recursive=TRUE,full.names=TRUE)){
 old<-sub('20261006_teacher_figures','20261006_figures_latest',f,fixed=TRUE)
 put(paste('ML predictions unchanged',f),same_table(old,f))
}
put('UMAP coordinates unchanged',same_table('outputs/20261006_figures_latest/ML/umap_coordinates.csv.gz',
                                            file.path(out,'ML/umap_coordinates.csv.gz')))
axis<-fread(file.path(out,'ML/axis_lower_bound_audit.csv'))
put('ML nonnegative-axis lower bounds',all(abs(axis$Lower)<1e-10),paste(nrow(axis),'panel axes'))
st<-fread(file.path(out,'statistics/pairwise_source_clustered.csv'))
put('Pairwise tests count',nrow(st)==68,'8 Figure1a + 4 Figure1b + 56 Figure2')
put('BH matches prescribed families',all(sapply(split(st,st$Family),function(d){
 isTRUE(all.equal(d$Q,p.adjust(d$P,'BH',n=nrow(d)),check.attributes=FALSE))})))
put('No stars assigned to unavailable contrasts',all(is.na(st[Status!='estimable',Q])))
# Inspect actual built grid grobs, keeping scientific lines and heatmap cells.
grid_lines<-function(g){
 n<-as.integer(grepl('panel.grid',g$name,fixed=TRUE)&&any(c('polyline','segments')%in%class(g)))
 if(length(g$grobs))n<-n+sum(vapply(g$grobs,grid_lines,integer(1)))
 if(length(g$children))n<-n+sum(vapply(g$children,grid_lines,integer(1)))
 n
}
figs<-c(list.files(file.path(out,'plot_objects'),pattern='\\.rds$',full.names=TRUE),
 list.files(file.path(out,'ML'),pattern='\\.rds$',full.names=TRUE,recursive=TRUE),
 list.files(file.path(out,'main_figures'),pattern='\\.rds$',full.names=TRUE))
ga<-rbindlist(lapply(figs,function(f){
 g<-readRDS(f);gr<-if(inherits(g,'patchwork'))patchworkGrob(g)else if(inherits(g,'ggplot'))ggplotGrob(g)else g
 count<-grid_lines(gr);cat('GRID_QA',basename(f),count,'\n')
 data.table(File=f,GridPolylines=count)
}))
fwrite(ga,file.path(out,'qa/grid_audit.csv'))
put('All built plots have no background grid lines',all(ga$GridPolylines==0),paste(nrow(ga),'plots'))
v<-rbindlist(checks);fwrite(v,file.path(out,'qa/numerical_validation.csv'))
print(v[,.(Checks=.N,Passed=sum(Pass))]);stopifnot(all(v$Pass))

# Final export checks for the panels corrected during visual review.
final<-list();add_final<-function(test,pass,detail='')final[[length(final)+1]]<<-data.table(Test=test,Pass=pass,Detail=detail)
no_dividers<-function(p){
 ok<-!any(vapply(p$layers,function(l)inherits(l$geom,'GeomVline')&&identical(l$aes_params$linetype,'dashed'),logical(1)))
 if(inherits(p,'patchwork'))ok<-ok&&all(vapply(p$patches$plots,no_dividers,logical(1)))
 ok
}
for(stem in c('Figure_S1a_DDR_fraction_by_PXD','Figure_S1b_MKI67_over_H3C1_by_PXD','Figure_S1b_MKI67_over_H3C1_11_detected_only')){
 g<-readRDS(file.path(out,'plot_objects',paste0('s1__',stem,'.rds')))
 add_final(paste('No decorative dividers',stem),no_dividers(g))
}
png_size<-function(f){
 con<-file(f,'rb');on.exit(close(con));seek(con,16);readBin(con,integer(),n=2,size=4,endian='big')
}
add_final('Lactate route original canvas',identical(png_size(file.path(out,'lactate/figures/Lactate_routes.png')),c(2080L,960L)),'13 x 6 inches at 160 dpi; white background visually checked')
mat<-c(Figure_S3a_DDR_pathway_matrix_tumor_tissues=192,Figure_S3b_DDR_pathway_matrix_non_tumor_tissues=183,
       Figure_S4a_DDR_pathway_matrix_cancer_cell_lines=381,Figure_S4b_DDR_pathway_matrix_normal_cell_lines=292)
for(stem in names(mat)){
 folder<-if(startsWith(stem,'Figure_S3'))'Supplementary_Figure_S3' else 'Supplementary_Figure_S4'
 base<-file.path(out,'protein/matrix_root/results/final_figures_and_tables',folder,stem)
 txt<-paste(system2('pdftotext',c(shQuote(paste0(base,'.pdf')),'-'),stdout=TRUE),collapse=' ')
 add_final(paste('Matrix title and original canvas',stem),grepl(paste0('(',mat[[stem]],' proteins)'),txt,fixed=TRUE)&&identical(png_size(paste0(base,'.png')),c(3300L,1440L)),'PNG rendered from checked Cairo PDF; title also visually checked')
}
final<-rbindlist(final);fwrite(final,file.path(out,'qa/final_export_validation.csv'));stopifnot(all(final$Pass))
