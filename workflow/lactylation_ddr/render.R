#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(patchwork);library(ragg)})
set.seed(25)
args<-commandArgs(TRUE)
release<-if(length(args))args[1] else 'deliverables/20261011_lactylation_ddr'
figdir<-file.path(release,'figures');tabdir<-file.path(release,'tables')
dir.create(figdir,recursive=TRUE,showWarnings=FALSE)
m<-fread(file.path(tabdir,'display_mapping_31.csv'))
genes<-fread(file.path(tabdir,'gene_panel.csv'))
rna<-fread(cmd=paste('gzip -dc',shQuote(file.path(tabdir,'rna_values_28_references.csv.gz'))))
stopifnot(nrow(rna)==12*1898,uniqueN(rna$ReferenceKey)==28,!anyDuplicated(rna[,.(Ensembl,SampleID)]),nrow(m)==31)
# Reference expansion changes display positions only; source sample IDs remain unchanged.
d<-merge(rna,m[,.(ReferenceKey,GroupID,X)],by='ReferenceKey',allow.cartesian=TRUE)
prot<-fread(file.path(release,'inputs','proteome_samples.csv'))
prot<-merge(prot,m[,.(PXD,SampleGroup,X,GroupID)],by=c('PXD','SampleGroup'),all.x=TRUE)
stopifnot(nrow(prot)==272,!anyNA(prot$X),all(prot$ObservationType=='sample'))
ps<-prot[,.(N=.N,Mean=mean(DdrFractionPercentage),SEM=if(.N>1)sd(DdrFractionPercentage)/sqrt(.N) else NA_real_),by=.(GroupID,X,Dataset)]
ps[,Position:=X+ifelse(Dataset=='Whole proteome',-.19,.19)]
fwrite(ps,file.path(tabdir,'ddr_bars_mean_sem.csv'))
means<-rna[,.(Mean=mean(Expression),N=.N),by=.(Ensembl,Symbol,ReferenceKey)]
fwrite(means,file.path(tabdir,'rna_line_means.csv'))
ink<-'#262B33';font<-'Arial Unicode MS'
gene_col<-c(LDHA='#DE6B00',LDHB='#785AF2',EP300='#252525',CREBBP='#007DB3',KAT5='#B84D76',AARS1='#159878',AARS2='#A27719',HDAC1='#1B5D9D',HDAC2='#AD342B',HDAC3='#658B28',SIRT2='#886093',SIRT3='#319FBB')
catcol<-c(normal_tissue='#0072B2',cancer_tissue='#D55E00',normal_cells='#009E73',cancer_cells='#CC79A7')
catlabels<-c('Non-tumor tissues','Tumor tissues','Normal cell lines','Cancer cell lines')
ranges<-m[,.(xmin=min(X)-.47,xmax=max(X)+.47,xmid=mean(range(X)),Label=catlabels[first(CategoryOrder)]),by=Category]
bg<-m[X%%2==1,.(xmin=X-.5,xmax=X+.5)]
base<-theme_minimal(base_family=font,base_size=13)+theme(panel.grid.minor=element_blank(),panel.grid.major.x=element_blank(),panel.grid.major.y=element_line(colour='#E5E9EF',linewidth=.35),axis.text=element_text(colour=ink),axis.title.y=element_text(face='bold',size=13,margin=margin(r=10)),plot.margin=margin(4,10,4,10),legend.title=element_blank(),legend.position='top',legend.text=element_text(size=12),plot.background=element_rect(fill='white',colour=NA))
xscale<-function(labels=FALSE)scale_x_continuous(limits=c(.5,31.5),breaks=1:31,labels=if(labels)m$PlotLabel else NULL,expand=c(0,0))
backdrop<-function()geom_rect(data=bg,aes(xmin=xmin,xmax=xmax,ymin=-Inf,ymax=Inf),inherit.aes=FALSE,fill='#F4F6F9',colour=NA)
track<-ggplot(ranges)+geom_rect(aes(xmin=xmin,xmax=xmax,ymin=.08,ymax=.92,fill=Category),colour=NA)+geom_text(aes(xmid,.5,label=Label),family=font,colour='white',size=4.5)+scale_fill_manual(values=catcol,guide='none')+xscale()+scale_y_continuous(limits=c(0,1),expand=c(0,0))+theme_void(base_family=font)+theme(plot.margin=margin(0,10,0,10))
# Positive percentages increase downwards from the central material axis.
pbar<-ggplot(ps,aes(Position,Mean,fill=Dataset))+backdrop()+geom_col(width=.33,colour=ink,linewidth=.25)+geom_errorbar(data=ps[N>1],aes(ymin=pmax(0,Mean-SEM),ymax=Mean+SEM),width=.14,linewidth=.4)+scale_fill_manual(values=c('Whole proteome'='#4E79A7','Lactylome (Kla)'='#F28E2B'),breaks=c('Whole proteome','Lactylome (Kla)'))+xscale()+scale_y_reverse(limits=c(22,0),breaks=seq(0,20,5),labels=function(z)paste0(z,'%'),expand=c(0,0))+labs(x=NULL,y='GO-DDR protein fraction (%)')+base+theme(legend.position='bottom',axis.text.x=element_blank(),axis.line.x=element_blank(),legend.margin=margin(0,0,0,0))
make_plot<-function(selected,title){
 dd<-d[Symbol %in% selected]
 dd[,Offset:=(match(Symbol,selected)-(length(selected)+1)/2)*.8/length(selected)]
 dd[,Position:=X+Offset]
 sm<-dd[,.(Mean=mean(Expression)),by=.(Ensembl,Symbol,X,Position)]
 ng<-length(selected)
 prna<-ggplot(dd,aes(Position,Expression,colour=Symbol,fill=Symbol))+backdrop()+geom_boxplot(aes(group=interaction(Ensembl,GroupID)),width=.62/ng,outlier.shape=NA,linewidth=.38,alpha=.24)+geom_point(position=position_jitter(width=.07/ng,height=0,seed=25),size=if(ng>7).25 else .4,alpha=.28,show.legend=FALSE)+geom_line(data=sm,aes(Position,Mean,colour=Symbol,group=Ensembl),inherit.aes=FALSE,linetype='dashed',linewidth=.55,alpha=.95,show.legend=FALSE)+scale_colour_manual(values=gene_col,breaks=selected)+scale_fill_manual(values=gene_col,breaks=selected,guide='none')+xscale(TRUE)+scale_y_continuous(limits=c(-1.2,13),breaks=c(0,4,8,12),expand=c(0,0))+labs(x=NULL,y='RNA expression\nqsmooth B, log2(TPM + 0.5)')+base+guides(colour=guide_legend(nrow=if(ng>7)2 else 1,byrow=TRUE,override.aes=list(alpha=1,linewidth=.7)))+theme(axis.text.x=element_text(angle=55,hjust=1,vjust=1,size=11,face='bold'),axis.ticks.x=element_line(colour='#B3BAC3',linewidth=.3),axis.ticks.length=grid::unit(3,'pt'),legend.margin=margin(0,0,4,0),plot.margin=margin(2,10,12,10))
 caption<-paste('RNA: boxes show median and IQR; whiskers extend to 1.5 x IQR. Points are sample records; dashed lines connect material-level means.',
 'DDR: downward bars show mean percentages; error bars show SEM (n > 1). Blank positions indicate no sample-level observation for that modality.',
 '* Shared RNA reference: HCT116 (3 columns), HK-2 (2 columns). 31 display columns correspond to 28 RNA references; reference cohorts and conditions are listed in the mapping table.',sep='\n')
 wrap_plots(track,prna,pbar,ncol=1,heights=c(.25,6.3,2.6))+plot_annotation(title=title,caption=caption,theme=theme(plot.title=element_text(family=font,size=21,face='bold',colour=ink,margin=margin(b=9)),plot.caption=element_text(family=font,size=10,hjust=0,colour='#65717D',lineheight=1.25),plot.margin=margin(12,14,10,14)))
}
sets<-list(
 Overview_12genes=list(genes$Symbol,'Lactylation regulator expression and DDR protein fractions'),
 Lactate_and_KATs=list(c('LDHA','LDHB','EP300','CREBBP','KAT5'),'Lactate metabolism / KAT-family lactyltransferases and DDR protein fractions'),
 Lactate_and_AARS=list(c('LDHA','LDHB','AARS1','AARS2'),'Lactate metabolism / AARS lactyltransferases and DDR protein fractions'),
 Delactylation_regulators=list(c('HDAC1','HDAC2','HDAC3','SIRT2','SIRT3'),'Delactylation-associated regulators and DDR protein fractions'))
plots<-lapply(sets,function(s)make_plot(s[[1]],s[[2]]))
for(n in names(plots)){
 ggsave(file.path(figdir,paste0(n,'.png')),plots[[n]],width=24,height=12,device=agg_png,dpi=300,limitsize=FALSE)
 ggsave(file.path(figdir,paste0(n,'.pdf')),plots[[n]],width=24,height=12,device=cairo_pdf,limitsize=FALSE)
}
cairo_pdf(file.path(release,'Combined_RNA_DDR.pdf'),width=24,height=12,onefile=TRUE,family=font)
for(p in plots)print(p)
dev.off()
display_means<-d[,.(DisplayMean=mean(Expression)),by=.(Ensembl,Symbol,ReferenceKey,GroupID)]
mean_check<-merge(display_means,means,by=c('Ensembl','Symbol','ReferenceKey'))
checks<-data.table(Check=c('All RNA points within axes','All DDR means and SEM within axes','Mean lines equal original RNA means','31 material positions retained','272 protein observations retained'),Pass=c(all(rna$Expression>=-1.2 & rna$Expression<=13),all(ps$Mean+fifelse(is.na(ps$SEM),0,ps$SEM)<=22),nrow(mean_check)==12*31 && max(abs(mean_check$DisplayMean-mean_check$Mean))<1e-12,nrow(m)==31,nrow(prot)==272))
stopifnot(all(checks$Pass));fwrite(checks,file.path(tabdir,'revision_validation.csv'))
writeLines(capture.output(sessionInfo()),file.path(release,'sessionInfo.txt'))
cat('Four paired-panel figures exported; middle material axis, upper RNA boxes, lower DDR bars.\n')
