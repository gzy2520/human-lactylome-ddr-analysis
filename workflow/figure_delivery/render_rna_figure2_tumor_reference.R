#!/usr/bin/env Rscript
# Separate RNA panels with tumor tissue as the reference category.
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(ragg);library(clubSandwich)})
set.seed(25)
out<-Sys.getenv('KLA_RNA_FIGURE2_OUTPUT','outputs/20261008_rna_figure2_tumor_reference')
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
stopifnot(!anyDuplicated(src$ReferenceKey))
bi<-merge(b,src,by='ReferenceKey')
co<-data.table(Control=cats[2],Comparison=cats[-2])
fit_contrasts<-function(d,value,panel){
 d<-copy(d);d[,Category:=factor(Category,levels=cats)]
 d[,Weight:=1/.N,by=.(ConnGroup,Category)]
 m<-lm(reformulate('Category',response=value,intercept=FALSE),data=d,weights=Weight)
 V<-vcovCR(m,cluster=d$ConnGroup,type='CR2',inverse_var=FALSE)
 rbindlist(lapply(seq_len(nrow(co)),function(i){
  C<-matrix(0,1,length(coef(m)))
  C[1,match(paste0('Category',co$Control[i]),names(coef(m)))]<- -1
  C[1,match(paste0('Category',co$Comparison[i]),names(coef(m)))]<-1
  z<-as.data.table(linear_contrast(m,V,C,p_values=TRUE))
  data.table(Panel=panel,Control=co$Control[i],Comparison=co$Comparison[i],
   Estimate=z$Est,SE=z$SE,DF=z$df,LowerCI=z$CI_L,
   UpperCI=z$CI_U,P=z$p_val,ControlMaterials=sum(d$Category==co$Control[i]),
   ComparisonMaterials=sum(d$Category==co$Comparison[i]),
   ControlSources=uniqueN(d[Category==co$Control[i],ConnGroup]),
   ComparisonSources=uniqueN(d[Category==co$Comparison[i],ConnGroup]))
 }))
}
tests<-rbindlist(lapply(gs,function(gene_set){
 cbind(data.table(GeneSet=gene_set),fit_contrasts(ai[GeneSet==gene_set],'MedianExpression','RNA expression'))
}))
tests[,Q:=p.adjust(P,'BH',n=6)]
label_q<-function(q)ifelse(q<1e-4,'****',ifelse(q<.001,'***',ifelse(q<.01,'**',ifelse(q<.05,'*','ns'))))
tests[,Label:=label_q(Q)]
fwrite(tests,file.path(out,'tables/RNA_expression_pairwise_source_clustered.csv'))
bt<-fit_contrasts(bi,'ProliferationMedian','RNA proliferation')
bt[,`:=`(Q=p.adjust(P,'BH',n=3))];bt[,Label:=label_q(Q)]
stopifnot(nrow(tests)==6,nrow(bt)==3,all(is.finite(c(tests$Q,bt$Q))))
fwrite(bt,file.path(out,'tables/RNA_proliferation_pairwise_source_clustered.csv'))
fwrite(ai[,.(GroupID,ReferenceKey,GeneSet,Category,ConnGroup,Value=MedianExpression)],
       file.path(out,'tables/RNA_expression_inference_units.csv'))
fwrite(bi[,.(GroupID,ReferenceKey,Category,ConnGroup,Value=ProliferationMedian)],
       file.path(out,'tables/RNA_proliferation_inference_units.csv'))
stopifnot(file.copy(a_file,file.path(out,'tables',basename(a_file)),overwrite=TRUE),
          file.copy(b_file,file.path(out,'tables',basename(b_file)),overwrite=TRUE))

font<-'Arial Unicode MS';ink<-'#2F3437';mean_red<-'#C0392B'
gene_pal<-setNames(c('#4E79A7','#F28E2B'),gs)
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

# Bracket endpoints align with the corresponding boxes; colors identify gene sets.
ba<-copy(tests)
ba[,`:=`(Offset=ifelse(GeneSet==gs[1],-.16,.16),
          Y=5.86+(match(Comparison,cats[-2])-1)*.60+ifelse(GeneSet==gs[1],0,.27),
          Colour=gene_pal[GeneSet])]
ba[,`:=`(X1=pmin(match(Control,cats),match(Comparison,cats))+Offset,
          X2=pmax(match(Control,cats),match(Comparison,cats))+Offset)]
bb<-copy(bt);bb[,`:=`(X1=pmin(match(Control,cats),match(Comparison,cats)),
                      X2=pmax(match(Control,cats),match(Comparison,cats)),
                      Y=7.02+(match(Comparison,cats[-2])-1)*.45,Colour=ink)]
fwrite(ba,file.path(out,'tables/Figure_2a_displayed_brackets.csv'))
fwrite(bb,file.path(out,'tables/Figure_2b_displayed_brackets.csv'))
add_brackets<-function(p,d){
 lines<-rbind(d[,.(X=X1,Xend=X2,Y,Yend=Y,Colour)],
              d[,.(X=X1,Xend=X1,Y=Y-.09,Yend=Y,Colour)],
              d[,.(X=X2,Xend=X2,Y=Y-.09,Yend=Y,Colour)])
 p+geom_segment(data=lines,aes(x=X,xend=Xend,y=Y,yend=Yend,colour=Colour),
                 inherit.aes=FALSE,linewidth=.30,show.legend=FALSE)+
   geom_text(data=d,aes(x=(X1+X2)/2,y=Y+.06,label=Label,colour=Colour),
             inherit.aes=FALSE,family=font,size=2.6,vjust=0,show.legend=FALSE)+
   scale_colour_identity()
}

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
 scale_fill_manual(values=gene_pal,
                   breaks=gs,labels=c('DDR genes','Kla protein genes'))+
 xscale()+scale_y_continuous(limits=c(2,7.65),breaks=seq(2,7,1),expand=c(0,0))+
 labs(title='RNA gene-set expression',subtitle='28 RNA reference materials',
                      x=NULL,y='Median RNA expression\nlog2(qsmooth TPM + 0.5)')+clean+
 guides(fill=guide_legend(nrow=1,override.aes=list(size=1.5)))
pa<-add_brackets(pa,ba)

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
 scale_y_continuous(limits=c(1,8.35),breaks=seq(1,8,1),expand=c(0,0))+
 labs(title='RNA proliferation score',subtitle='Hallmark G2M Checkpoint · 28 RNA materials',
                      x=NULL,y='Proliferation score\n(log2 expression)')+clean
pb<-add_brackets(pb,bb)

grid_lines<-function(g){
 n<-as.integer(grepl('panel.grid',g$name,fixed=TRUE)&&any(c('polyline','segments')%in%class(g)))
 if(length(g$grobs))n<-n+sum(vapply(g$grobs,grid_lines,integer(1)))
 if(length(g$children))n<-n+sum(vapply(g$children,grid_lines,integer(1)))
 n
}
checks<-list()
save_panel<-function(p,stem,values,brackets){
 # Use the same font backend for the grob audit as for the delivered raster.
 agg_png(tempfile(fileext='.png'),width=4,height=4.1,units='in',res=100)
 built<-ggplot_build(p)
 expected_offset<-if('Offset'%in%names(brackets))brackets$Offset else 0
 endpoints_match<-all(brackets$X1==pmin(match(brackets$Control,cats),match(brackets$Comparison,cats))+expected_offset)&
                  all(brackets$X2==pmax(match(brackets$Control,cats),match(brackets$Comparison,cats))+expected_offset)
 checks[[length(checks)+1]]<<-data.table(Panel=stem,InputPoints=length(values),
  DrawnPoints=nrow(built$data[[3]]),Boxplots=nrow(built$data[[1]]),
  GridLines=grid_lines(ggplotGrob(p)),FinitePoints=all(is.finite(built$data[[3]]$y)),
  PointValuesUnchanged=isTRUE(all.equal(sort(built$data[[3]]$y),sort(values),tolerance=1e-12)),
  ComparisonBrackets=nrow(built$data[[4]])/3,
  ComparisonLabels=nrow(built$data[[5]]),
  LabelsMatchStatistics=identical(built$data[[5]]$label,brackets$Label),
  BracketEndpointsMatch=endpoints_match,
  AnnotationColoursMatch=identical(built$data[[5]]$colour,brackets$Colour),
  AllControlsTumor=all(brackets$Control=='cancer_tissue'))
 dev.off()
 for(ext in c('png','pdf','svg'))ggsave(file.path(out,'figures',paste0(stem,'.',ext)),p,
  width=4,height=4.1,dpi=600,bg='white',device=if(ext=='png')agg_png else if(ext=='pdf')cairo_pdf else grDevices::svg)
}
save_panel(pa,'Figure_2a_RNA_DDR_and_Kla_gene_expression',a$MedianExpression,ba)
save_panel(pb,'Figure_2b_RNA_proliferation',b$ProliferationMedian,bb)
v<-rbindlist(checks);stopifnot(all(v$InputPoints==v$DrawnPoints),all(v$GridLines==0),all(v$FinitePoints),
  all(v$PointValuesUnchanged),identical(v$Boxplots,c(8L,4L)),
  identical(v$ComparisonBrackets,c(6,3)),identical(v$ComparisonLabels,c(6L,3L)),
  all(v$LabelsMatchStatistics),all(v$BracketEndpointsMatch),all(v$AnnotationColoursMatch),all(v$AllControlsTumor))
fwrite(v,file.path(out,'tables/render_validation.csv'))
writeLines(trimws(capture.output(sessionInfo()),which='right'),file.path(out,'sessionInfo.txt'))
print(v)
