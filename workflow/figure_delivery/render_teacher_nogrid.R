#!/usr/bin/env Rscript
# Re-evaluate the reviewed source renderers with only grid/separator styling
# removed, then annotate prespecified source-aware comparisons on Figure 1/2.
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(patchwork)})
set.seed(25)
root<-normalizePath('.')
out<-normalizePath(Sys.getenv('KLA_TEACHER_OUTPUT','outputs/20261006_teacher_figures'))
rna_root<-Sys.getenv('KLA_RNA_MATRIX_ROOT','/Users/gzy2520/Desktop/Research/kla/.worktrees/kla31-rnaseq-20260912')
s<-fread(file.path(out,'statistics/pairwise_source_clustered.csv'),na.strings='')
plot_records<-list()
cat_order<-c('normal_tissue','cancer_tissue','normal_cells','cancer_cells')

style_expression<-function(x){
 if(missing(x))return(quote(expr=))
 if(is.expression(x))return(as.expression(lapply(x,style_expression)))
 if(!is.call(x))return(x)
 head<-paste(deparse(x[[1]]),collapse='')
 if(head=='geom_vline' && 'xintercept'%in%names(x)){
  v<-tryCatch(eval(x[['xintercept']]),error=function(e)NULL)
  dashed<-'linetype'%in%names(x) && identical(tryCatch(eval(x[['linetype']]),error=function(e)NULL),'dashed')
  if(identical(v,c(1.5,2.5,3.5)) || dashed)return(quote(geom_blank()))
 }
 if(head=='theme')for(i in seq_along(x)[-1])
  if(!is.null(names(x)) && grepl('^panel\\.grid',names(x)[i]))x[[i]]<-quote(element_blank())
 if(head=='ggplot2::ggsave')x[[1]]<-as.name('ggsave')
 for(i in seq_along(x)[-1])x[i]<-list(style_expression(x[[i]]))
 if(head=='theme_minimal')x<-call('+',x,quote(theme(panel.grid=element_blank())))
 x
}

annotate_current<-function(g,stem){
 st<-NULL
 if(stem=='Figure_1_DDR_fraction_candidate_category_boxplot_refined'){
  st<-s[Family=='Figure1a' & (MainBracket | (!is.na(Q)&Q<.05))]
  st[,Offset:=ifelse(Panel=='Whole proteome',-.165,.165)]
  g<-g+labs(title='DDR fraction within each proteomic assay',subtitle=NULL,
   y='DDR proteins / all detected proteins (%)')
 }else if(stem=='Figure_1b_RNA_proliferation_hallmark_28materials_boxplot'){
  st<-s[Family=='Figure1b'];st[,Offset:=0]
  g<-g+labs(title='Cell proliferation score\n28 independent RNA materials',subtitle=NULL,
   y='Proliferation score (log2 scale)')
 }else if(grepl('^Figure_2_DDR_pathway_summary_.*_(barplot|boxplot)$',stem)){
  pa<-sub('^Figure_2_DDR_pathway_summary_(.*)_(barplot|boxplot)$','\\1',stem)
  st<-s[Family=='Figure2' & Panel==pa & MainBracket & Status=='estimable']
  st[,Offset:=ifelse(Direction=='Pro',-.18,.18)]
  g<-g+labs(subtitle=NULL)
 }
 if(!is.null(st)){
  st<-copy(st);setorder(st,Control,Comparison,Offset)
  ylim<-g$scales$get_scales('y')$limits
  stopifnot(is.numeric(ylim),length(ylim)==2)
  step<-(ylim[2]-ylim[1])*.075
  st[,`:=`(Left=match(Control,cat_order)+Offset,Right=match(Comparison,cat_order)+Offset,
            Y=ylim[2]+seq_len(.N)*step)]
  st[,Text:=paste0(ifelse(nzchar(Direction),paste0(Direction,' '),
                     ifelse(Family=='Figure1a',ifelse(Panel=='Whole proteome','Proteome ','Kla '),'')),Label)]
  if(nrow(st)){
   tips<-rbind(st[,.(X=Left,Y,Yend=Y-step*.25)],st[,.(X=Right,Y,Yend=Y-step*.25)])
   g<-g+geom_segment(data=st,aes(x=Left,xend=Right,y=Y,yend=Y),inherit.aes=FALSE,linewidth=.35)+
    geom_segment(data=tips,aes(x=X,xend=X,y=Y,yend=Yend),inherit.aes=FALSE,linewidth=.35)+
    geom_text(data=st,aes(x=(Left+Right)/2,y=Y+step*.12,label=Text),inherit.aes=FALSE,
              size=3.3,family='Arial',vjust=0)
   ys<-g$scales$get_scales('y')
   ys$limits<-c(ylim[1],max(st$Y)+step*.8)
   ys$breaks<-scales::pretty_breaks(6)
  }
  cap<-'Source-balanced contrasts; study-clustered CR2 / BH q. ns: q >= 0.05.'
  if(any(st$Family=='Figure1a'))cap<-paste(cap,'Dots and n retain frozen display rows; reused references count once in tests.',sep='\n')
  if(grepl('pathway_summary',stem))cap<-paste(cap,'Pro and Inh tested separately; zero-variance contrasts are not annotated.',sep='\n')
  if(stem=='Figure_1b_RNA_proliferation_hallmark_28materials_boxplot')cap<-paste(cap,
    'Material median of the 196-gene Hallmark G2M mean log2(TPM + 0.5).',sep='\n')
  g<-g+labs(caption=cap)+theme(plot.caption=element_text(size=9,hjust=0,colour='#404040'),
                            plot.title=element_text(size=15,face='bold'))
 }
 g
}

run_renderer<-function(script,args,variables=list()){
 marker<-file.path(out,'render_completed',basename(script))
 if(file.exists(marker)){cat('ALREADY_RENDERED',script,'\n');return(invisible(NULL))}
 e<-new.env(parent=globalenv());e$commandArgs<-function(...)args
 previous<-setNames(Sys.getenv(names(variables),unset=NA_character_),names(variables))
 if(length(variables))do.call(Sys.setenv,variables)
 on.exit({for(n in names(variables))if(is.na(previous[[n]]))Sys.unsetenv(n)else do.call(Sys.setenv,setNames(list(previous[[n]]),n))},add=TRUE)
 e$ggsave<-function(filename,plot=last_plot(),...){
  stem<-sub('\\.(png|pdf)$','',basename(filename))
  if(inherits(plot,'ggplot')&&!inherits(plot,'patchwork')){
   # Save the effective global theme too; otherwise rebuilding an RDS in a
   # fresh R process can restore the factory grid for globally styled plots.
   plot$theme<-theme_get()+plot$theme
   plot<-annotate_current(plot,stem)
  }
  dir.create(dirname(filename),recursive=TRUE,showWarnings=FALSE)
  if(endsWith(filename,'.png')){
   rds<-file.path(out,'plot_objects',paste0(basename(dirname(filename)),'__',stem,'.rds'))
   dir.create(dirname(rds),recursive=TRUE,showWarnings=FALSE);saveRDS(plot,rds)
   plot_records[[length(plot_records)+1L]]<<-data.table(Renderer=script,File=filename,PlotRDS=rds)
  }
  options<-list(...)
  if(!'bg'%in%names(options))options$bg<-'white'
  do.call(ggplot2::ggsave,c(list(filename=filename,plot=plot),options))
 }
 expressions<-style_expression(parse(script))
 for(x in expressions)eval(x,envir=e)
 w<-warnings();if(length(w))print(head(w,6))
 dir.create(dirname(marker),recursive=TRUE,showWarnings=FALSE);writeLines(script,marker)
 if(length(plot_records))fwrite(rbindlist(plot_records),file.path(out,'renderer_outputs.csv'))
 cat('RENDERED',script,'\n')
 invisible(e)
}

frozen<-file.path(root,'corrected_final_result_20260905_sample_only/Corrected_Data/sample_only_inputs')
run_renderer('R/publication/build_publication_outputs.R',root,list(
 KLA_PUBLICATION_INPUT=file.path(root,'data/publication_input'),KLA_PUBLICATION_OUTPUT=file.path(out,'protein/publication')))
run_renderer('workflow/corrected_renderers/build_figure1_category_boxplot.R',root,list(
 KLA_CANDIDATE_INPUT=frozen,KLA_CANDIDATE_OUTPUT=file.path(out,'protein/fraction')))
run_renderer('workflow/corrected_renderers/build_figure1_mki67_ratio_boxplots.R',root,list(
 KLA_MKI67_INPUT=frozen,KLA_MKI67_OUTPUT=file.path(out,'protein/ratios')))
run_renderer('workflow/corrected_renderers/build_ddr_pathway_summary_boxplots.R',root,list(
 KLA_CANDIDATE_INPUT=frozen,KLA_CANDIDATE_OUTPUT=file.path(out,'protein/pathways')))
run_renderer('workflow/corrected_renderers/build_figure_s1_dataset_plots.R',root,list(
 KLA_CANDIDATE_INPUT=frozen,KLA_S1_OUTPUT=file.path(out,'protein/s1'),KLA_SAMPLE_ONLY='TRUE',
 KLA_FULL_FIGURE1_INPUT=file.path(root,'data/candidate/figure1_sample_boxplot_values.csv'),
 KLA_S1A_EXPECTED_ROWS='272',KLA_S1B_EXPECTED_ROWS='79'))
matrix_root<-file.path(out,'protein/matrix_root');dir.create(matrix_root,recursive=TRUE,showWarnings=FALSE)
run_renderer('R/candidate/build_split_pathway_matrices.R',matrix_root,list(
 KLA_PUBLICATION_INPUT=file.path(root,'data/publication_input')))
run_renderer('workflow/render_teacher_questions_rna_exact.R',c(root,file.path(out,'RNA/core'),
 file.path(rna_root,'outputs/20260919_repair_qsmooth'),file.path(rna_root,'outputs/20260919_expression_corrected'),
 file.path(rna_root,'outputs/20260919_repair_panel')))
run_renderer('workflow/plot_rna_reference_31group_20260917.R',c(root,file.path(out,'RNA/reference'),
 file.path(rna_root,'outputs/20260919_repair_qsmooth'),file.path(rna_root,'outputs/20260919_repair_panel'),
 file.path(rna_root,'outputs/20260919_repair_assisted')))
# The lactate renderer resolves its release config relative to the working
# directory. Read the frozen local matrix checkout; all writes use absolute out.
oldwd<-getwd();setwd(rna_root)
run_renderer(file.path(root,'workflow/lactate_metabolism/analyze.R'),file.path(out,'lactate'))
setwd(oldwd)
if(length(plot_records))fwrite(rbindlist(plot_records),file.path(out,'renderer_outputs.csv'))
# Match the earlier release's title fix: the Cairo PDF contains all digits,
# whereas the macOS raster font backend can omit title digits in these panels.
for(folder in c('Supplementary_Figure_S3','Supplementary_Figure_S4')){
 d<-file.path(matrix_root,'results/final_figures_and_tables',folder)
 for(p in list.files(d,pattern='\\.pdf$',full.names=TRUE)){
  status<-system2('pdftoppm',c('-png','-r','300','-scale-to-x','3300','-scale-to-y','1440',
                              '-singlefile',shQuote(p),shQuote(sub('\\.pdf$','',p))))
  stopifnot(status==0)
 }
}
