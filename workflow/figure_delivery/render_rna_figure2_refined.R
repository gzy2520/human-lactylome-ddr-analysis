#!/usr/bin/env Rscript
# Separate RNA panels for the supervisor's Figure 2; frozen expression scores.
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(ragg);library(clubSandwich)})
set.seed(25)
out<-Sys.getenv('KLA_RNA_FIGURE2_OUTPUT','outputs/20261007_rna_figure2_refined')
dir.create(file.path(out,'tables'),recursive=TRUE,showWarnings=FALSE)
dir.create(file.path(out,'figures'),recursive=TRUE,showWarnings=FALSE)
frozen<-'outputs/20260919_repair_release/core'
a_file<-file.path(frozen,'figure1a_rna_ddr_and_lactylated_28materials.csv')
b_file<-file.path(frozen,'figure1b_hallmark_proliferation_28materials.csv')
a<-fread(a_file);b<-fread(b_file)
cats<-c('normal_tissue','cancer_tissue','normal_cells','cancer_cells')
gs<-c('DDR genes','Lactylated protein genes')
stopifnot(nrow(a)==56,nrow(b)==28,uniqueN(b$ReferenceKey)==28,
          !anyDuplicated(a,by=c('GroupID','GeneSet')),all(is.finite(a$MedianExpression)),
          identical(as.integer(b[, .N,by=Category][match(cats,Category),N]),c(9L,3L,6L,10L)))
src<-unique(fread('outputs/20260925_sample_fraction_inputs_final/sample_metadata.csv')[,.(ReferenceKey,ConnGroup)])
map<-b[,.(GroupID,ReferenceKey)]
ai<-merge(merge(a,map,by='GroupID'),src,by='ReferenceKey')
stopifnot(nrow(ai)==56,!anyNA(ai$ConnGroup))
co<-data.table(Control=c(cats[1],cats[1],cats[1],cats[3]),
               Comparison=c(cats[2],cats[3],cats[4],cats[4]),MainBracket=c(TRUE,FALSE,FALSE,TRUE))
tests<-rbindlist(lapply(gs,function(gene_set){
 d<-copy(ai[GeneSet==gene_set]);d[,Category:=factor(Category,levels=cats)]
 d[,Weight:=1/.N,by=.(ConnGroup,Category)]
 m<-lm(MedianExpression~0+Category,data=d,weights=Weight)
 V<-vcovCR(m,cluster=d$ConnGroup,type='CR2',inverse_var=FALSE)
 rbindlist(lapply(seq_len(nrow(co)),function(i){
  C<-matrix(0,1,length(coef(m)))
  C[1,match(paste0('Category',co$Control[i]),names(coef(m)))]<- -1
  C[1,match(paste0('Category',co$Comparison[i]),names(coef(m)))]<-1
  z<-as.data.table(linear_contrast(m,V,C,p_values=TRUE))
  data.table(GeneSet=gene_set,Control=co$Control[i],Comparison=co$Comparison[i],
   MainBracket=co$MainBracket[i],Estimate=z$Est,SE=z$SE,DF=z$df,LowerCI=z$CI_L,
   UpperCI=z$CI_U,P=z$p_val,ControlMaterials=sum(d$Category==co$Control[i]),
   ComparisonMaterials=sum(d$Category==co$Comparison[i]))
 }))
}))
tests[,Q:=p.adjust(P,'BH',n=8)]
label_q<-function(q)ifelse(q<1e-4,'****',ifelse(q<.001,'***',ifelse(q<.01,'**',ifelse(q<.05,'*','ns'))))
tests[,Label:=label_q(Q)]
fwrite(tests,file.path(out,'tables/RNA_expression_pairwise_source_clustered.csv'))
bt<-fread('outputs/20261006_teacher_figures/statistics/pairwise_source_clustered.csv')[Family=='Figure1b']
fwrite(bt,file.path(out,'tables/RNA_proliferation_pairwise_frozen.csv'))
stopifnot(file.copy(a_file,file.path(out,'tables',basename(a_file)),overwrite=TRUE),
          file.copy(b_file,file.path(out,'tables',basename(b_file)),overwrite=TRUE))

font<-'Arial Unicode MS';ink<-'#2F3437';mean_red<-'#C0392B'
pal<-c(normal_tissue='#0072B2',cancer_tissue='#D55E00',normal_cells='#009E73',cancer_cells='#CC79A7')
counts<-b[, .N,by=Category][match(cats,Category)]
xlabs<-paste0(c('Non-tumor\ntissues','Tumor\ntissues','Normal\ncell lines','Cancer\ncell lines'),
              '\n(n = ',counts$N,')')
clean<-theme_classic(base_family=font,base_size=9)+theme(
 panel.grid=element_blank(),axis.line=element_line(linewidth=.35,colour=ink),
 axis.ticks=element_line(linewidth=.3,colour=ink),axis.ticks.length=unit(1.1,'mm'),
 axis.text=element_text(size=8,colour=ink),axis.text.x=element_text(lineheight=1.08,margin=margin(t=5)),
 axis.title.y=element_text(size=9,margin=margin(r=7)),
 plot.title=element_text(size=10,face='bold',hjust=0,margin=margin(b=4)),
 plot.subtitle=element_text(size=8,colour='#65717D',hjust=0,margin=margin(b=5)),
 legend.position='top',legend.justification='left',legend.text=element_text(size=8),
 legend.title=element_blank(),legend.key.width=unit(3.7,'mm'),legend.key.height=unit(3,'mm'),
 legend.spacing.x=unit(1.5,'mm'),legend.margin=margin(0,0,3,0),
 plot.margin=margin(7,8,7,7),plot.background=element_rect(fill='white',colour=NA))
xscale<-function()scale_x_continuous(breaks=1:4,labels=xlabs,limits=c(.5,4.5),expand=c(0,0))

a[,`:=`(GeneSet=factor(GeneSet,levels=gs),X=match(Category,cats))]
sa<-a[,.(Mean=mean(MedianExpression),Median=median(MedianExpression),N=.N),by=.(X,GeneSet)]
sa[,Pos:=X+ifelse(GeneSet==gs[1],-.16,.16)]
pa<-ggplot(a,aes(X,MedianExpression,fill=GeneSet))+
 geom_boxplot(aes(group=interaction(X,GeneSet)),position=position_dodge(.64),width=.48,
              outlier.shape=NA,linewidth=.32,median.linewidth=.45,alpha=.8,colour=ink)+
 geom_segment(data=sa,aes(x=Pos-.105,xend=Pos+.105,y=Mean,yend=Mean),inherit.aes=FALSE,
              colour=mean_red,linewidth=.42)+
 geom_point(position=position_jitterdodge(jitter.width=.10,dodge.width=.64,seed=25),
            shape=21,size=1.55,stroke=.22,alpha=.9,colour='white')+
 scale_fill_manual(values=setNames(c('#4E79A7','#F28E2B'),gs),
                   breaks=gs,labels=c('DDR genes','Kla protein genes'))+
 xscale()+scale_y_continuous(limits=c(2,6),breaks=seq(2,6,1),expand=c(0,0))+
 labs(title='RNA gene-set expression',subtitle='28 RNA reference materials',
                      x=NULL,y='Median RNA expression\nlog2(qsmooth TPM + 0.5)')+clean+
 guides(fill=guide_legend(nrow=1,override.aes=list(size=1.5)))

b[,X:=match(Category,cats)]
sb<-b[,.(Mean=mean(ProliferationMedian),Median=median(ProliferationMedian),N=.N),by=X]
pb<-ggplot(b,aes(X,ProliferationMedian,fill=Category))+
 geom_boxplot(aes(group=Category),width=.46,outlier.shape=NA,linewidth=.32,
              median.linewidth=.45,alpha=.8,colour=ink)+
 geom_segment(data=sb,aes(x=X-.18,xend=X+.18,y=Mean,yend=Mean),inherit.aes=FALSE,
              colour=mean_red,linewidth=.42)+
 geom_point(position=position_jitter(width=.075,height=0,seed=25),shape=21,size=1.85,
            stroke=.25,alpha=.94,colour='white')+
 scale_fill_manual(values=pal,guide='none')+xscale()+
 scale_y_continuous(limits=c(1,7),breaks=seq(1,7,1),expand=c(0,0))+
 labs(title='RNA proliferation score',subtitle='Hallmark G2M Checkpoint · 28 RNA materials',
                      x=NULL,y='Proliferation score\n(log2 expression)')+clean

grid_lines<-function(g){
 n<-as.integer(grepl('panel.grid',g$name,fixed=TRUE)&&any(c('polyline','segments')%in%class(g)))
 if(length(g$grobs))n<-n+sum(vapply(g$grobs,grid_lines,integer(1)))
 if(length(g$children))n<-n+sum(vapply(g$children,grid_lines,integer(1)))
 n
}
checks<-list()
save_panel<-function(p,stem,values){
 # Use the same font backend for the grob audit as for the delivered raster.
 agg_png(tempfile(fileext='.png'),width=3.8,height=3.45,units='in',res=100)
 built<-ggplot_build(p)
 checks[[length(checks)+1]]<<-data.table(Panel=stem,InputPoints=length(values),
  DrawnPoints=nrow(built$data[[3]]),Boxplots=nrow(built$data[[1]]),
  GridLines=grid_lines(ggplotGrob(p)),FinitePoints=all(is.finite(built$data[[3]]$y)),
  PointValuesUnchanged=isTRUE(all.equal(sort(built$data[[3]]$y),sort(values),tolerance=1e-12)),
  NoComparisonAnnotations=length(p$layers)==3)
 dev.off()
 for(ext in c('png','pdf','svg'))ggsave(file.path(out,'figures',paste0(stem,'.',ext)),p,
  width=3.8,height=3.45,dpi=600,bg='white',device=if(ext=='png')agg_png else if(ext=='pdf')cairo_pdf else grDevices::svg)
}
save_panel(pa,'Figure_2a_RNA_DDR_and_Kla_gene_expression',a$MedianExpression)
save_panel(pb,'Figure_2b_RNA_proliferation',b$ProliferationMedian)
v<-rbindlist(checks);stopifnot(all(v$InputPoints==v$DrawnPoints),all(v$GridLines==0),all(v$FinitePoints),
                            all(v$PointValuesUnchanged),all(v$NoComparisonAnnotations),identical(v$Boxplots,c(8L,4L)))
fwrite(v,file.path(out,'tables/render_validation.csv'))
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
print(v)
