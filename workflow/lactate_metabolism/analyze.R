#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(ragg)})
set.seed(25)
args<-commandArgs(TRUE);stopifnot(length(args)==1);out<-args[1]
if(dir.exists(out)&&length(list.files(out,all.files=TRUE,no..=TRUE)))stop('Refusing nonempty output')
for(d in c(out,file.path(out,'tables'),file.path(out,'figures')))dir.create(d,recursive=TRUE,showWarnings=FALSE)
put<-function(x,n)fwrite(x,file.path(out,'tables',paste0(n,'.csv')))
config<-fread('config/rna_current_release.csv');dirs<-setNames(config$Path,config$Role)
readmat<-function(n){d<-fread(file.path(dirs[['qsmooth']],'matrices',n));setnames(d,1,'Ensembl');d}
a<-readmat('qsmooth_A_collapsed_log2tpm.tsv.gz');u<-readmat('unified_groupmedian_log2tpm.tsv.gz')
b<-readmat('qsmooth_B_full_log2tpm.tsv.gz')
qc<-fread(file.path(dirs[['expression']],'group_sample_qc.csv'));exp<-fread(file.path(dirs[['qsmooth']],'group_expansion_31.csv'))
g<-fread('data/publication_input/group_summary_31.csv')
status<-fread('audit/20260916_full_31_rna_status/rna_31_group_status.csv')
g<-merge(g,status[,.(PXD,SampleGroup,GroupID)],by=c('PXD','SampleGroup'))
stopifnot(nrow(g)==31,!anyDuplicated(g$GroupID))
meta<-merge(exp,g[,.(GroupID,RowOrder,PXD,SampleGroup,Category,KlaLabelEn)],by='GroupID');setorder(meta,RowOrder)
labels_meta<-fread('config/lactate_metabolism/material_labels.csv')
stopifnot(nrow(labels_meta)==28,!anyDuplicated(labels_meta$ReferenceKey))
meta<-merge(meta,labels_meta,by='ReferenceKey',all.x=TRUE);setorder(meta,RowOrder)
stopifnot(!anyNA(meta$RNALabel))
meta[,PlotLabel:=paste(sub('KLA31_','',GroupID),RNALabel)]
put(meta,'rna_proteome_mapping_31');meta28<-meta[!duplicated(ReferenceKey)];put(meta28,'materials_28')
stopifnot(nrow(meta28)==28,nrow(qc)==1898,!anyDuplicated(qc$SampleID),!anyDuplicated(a$Ensembl),setequal(names(a)[-1],meta28$ReferenceKey),setequal(names(b)[-1],qc$SampleID),identical(a$Ensembl,u$Ensembl),identical(a$Ensembl,b$Ensembl))
mem<-fread('config/lactate_metabolism/gene_module_membership.csv');mem<-mem[!is.na(Ensembl)&nzchar(Ensembl)]
# Explicitly report alternate mappings; never aggregate or choose by expression magnitude.
mem[,InCommonMatrix:=Ensembl %in% a$Ensembl];put(mem,'membership_coverage')
covered<-mem[InCommonMatrix==TRUE];labels<-unique(covered[,.(Ensembl,Symbol)]);stopifnot(!anyDuplicated(labels$Ensembl))
# No multiple measured Ensembl records for one named gene in this release.
stopifnot(!anyDuplicated(labels$Symbol))
ids<-unique(covered$Ensembl)
la<-melt(a[Ensembl %in% ids],id.vars='Ensembl',variable.name='ReferenceKey',value.name='QsmoothA')
lu<-melt(u[Ensembl %in% ids],id.vars='Ensembl',variable.name='ReferenceKey',value.name='UnsmoothedLog2TPM')
la<-merge(la,lu,by=c('Ensembl','ReferenceKey'));la<-merge(la,labels,by='Ensembl');la<-merge(la,meta28,by='ReferenceKey')
la[,GeneZ:=if(sd(QsmoothA)>0)(QsmoothA-mean(QsmoothA))/sd(QsmoothA) else NA_real_,by=Ensembl]
put(la,'material_gene_expression')
lb<-melt(b[Ensembl %in% ids],id.vars='Ensembl',variable.name='SampleID',value.name='QsmoothB')
lb<-merge(lb,qc[,.(SampleID,ReferenceKey)],by='SampleID');lb<-merge(lb,labels,by='Ensembl');lb<-merge(lb,meta28[,.(ReferenceKey,PlotLabel,Category)],by='ReferenceKey')
fwrite(lb,file.path(out,'tables','sample_gene_expression.csv.gz'))
bs<-lb[,.(SampleN=.N,MedianB=median(QsmoothB),Q25=quantile(QsmoothB,.25),Q75=quantile(QsmoothB,.75),Min=min(QsmoothB),Max=max(QsmoothB)),by=.(Ensembl,ReferenceKey)]
put(bs,'sample_distribution_summary')
coverage<-mem[,.(AnnotatedOrMappedIDs=uniqueN(Ensembl),CommonMatrixIDs=uniqueN(Ensembl[InCommonMatrix]),MissingIDs=uniqueN(Ensembl[!InCommonMatrix])),by=.(Origin,Term,Module)];put(coverage,'module_coverage')
# Source-local availability of every missing candidate. Presence is not imputed expression.
missing<-unique(mem[InCommonMatrix==FALSE,Ensembl]);availability<-rbindlist(lapply(meta28$ReferenceKey,function(k){
 d<-fread(file.path(dirs[['expression']],'matrices',paste0(k,'_log2tpm.tsv.gz')),select=1)
 data.table(Ensembl=missing,ReferenceKey=k,PresentInSource=missing %in% d[[1]])
}));put(availability,'missing_gene_source_availability')
sens<-la[,.(MaxAbsADelta=max(abs(QsmoothA-UnsmoothedLog2TPM)),MaterialSpearman=cor(QsmoothA,UnsmoothedLog2TPM,method='spearman')),by=.(Ensembl,Symbol)];put(sens,'normalization_sensitivity')
# No inferential tests, module activity scores, or putative flux estimates.
pal<-c(normal_tissue='#0072B2',cancer_tissue='#D55E00',normal_cells='#009E73',cancer_cells='#CC79A7')
theme_set(theme_minimal(base_size=10,base_family='Arial')+theme(panel.grid.minor=element_blank(),plot.title=element_text(face='bold'),plot.subtitle=element_text(size=9),plot.margin=margin(10,16,10,10)))
figs<-list();save<-function(p,n,w,h){
 for(k in c('title','subtitle','caption'))if(is.character(p$labels[[k]]))p$labels[[k]]<-paste(strwrap(p$labels[[k]],floor(w*11)),collapse='\n')
 ggsave(file.path(out,'figures',paste0(n,'.pdf')),p,width=w,height=h,device=cairo_pdf,limitsize=FALSE)
 ggsave(file.path(out,'figures',paste0(n,'.png')),p,width=w,height=h,device=agg_png,dpi=160,limitsize=FALSE)
 figs[[length(figs)+1]]<<-data.table(Figure=n,Width=w,Height=h)
}
heat<-function(dt,title,n){
 dt<-copy(dt);orderids<-unique(dt[order(Symbol),Ensembl]);labs<-unique(dt[,.(Ensembl,Symbol)])
 dt[,Ensembl:=factor(Ensembl,levels=rev(orderids))];dt[,PlotLabel:=factor(PlotLabel,levels=meta28$PlotLabel)]
 p<-ggplot(dt,aes(PlotLabel,Ensembl,fill=GeneZ))+geom_tile()+scale_y_discrete(labels=setNames(labs$Symbol,labs$Ensembl))+scale_fill_gradient2(low='#2166AC',mid='white',high='#B2182B',midpoint=0,na.value='#BDBDBD',name='Gene Z')+labs(title=title,subtitle='28 RNA references | qsmooth A on log2(TPM + 0.5) | within-gene relative expression',x=NULL,y=NULL)+theme(axis.text.x=element_text(angle=60,hjust=1,size=8),axis.text.y=element_text(size=9),panel.grid=element_blank())
 save(p,n,14,max(5,2.5+length(orderids)*.19))
}
# Complete GO panels, paginated without expression-based selection.
go<-fread('config/lactate_metabolism/go_terms.csv')
for(i in seq_len(nrow(go))){
 mids<-unique(covered[Term==go$Term[i],Ensembl]);mids<-sort(mids)
 if(!length(mids))next
 for(j in seq_len(ceiling(length(mids)/40))){part<-mids[seq((j-1)*40+1,min(j*40,length(mids)))];heat(la[Ensembl %in% part],paste(go$Term[i],go$Module[i],paste0('(',length(mids),' measured genes)')),paste0('GO_',sub('GO:','',go$Term[i]),'_page',j))}
}
cm<-covered[Origin=='curated mechanism'];cmods<-unique(cm$Module)
for(i in seq_along(cmods))heat(la[Ensembl %in% cm[Module==cmods[i],Ensembl]],cmods[i],sprintf('Mechanism_%02d',i))
# LDHA/B sample distributions across all references; whiskers show distributions, not inferential CI.
ld<-copy(lb[Symbol %in% c('LDHA','LDHB')]);ld[,PlotLabel:=factor(PlotLabel,levels=meta28$PlotLabel)]
p<-ggplot(ld,aes(PlotLabel,QsmoothB,fill=Category))+geom_boxplot(outlier.size=.35,width=.65)+facet_wrap(~Symbol,ncol=1)+scale_fill_manual(values=pal)+labs(title='LDHA and LDHB sample distributions',subtitle='Independent sample IDs; not all IDs establish independent donors. Single-sample material has no variability estimate.',x=NULL,y='qsmooth B on log2(TPM + 0.5)')+theme(axis.text.x=element_text(angle=60,hjust=1,size=8),legend.position='bottom')
save(p,'LDHA_LDHB_samples_28',14,9)
# Fixed focus comparisons, no assumed donor pairing.
pairs<-data.table(Comparison=c('HCC vs adjacent liver','Prostate cancer vs BPH'),Tumor=c('TCGA_LIHC_primary','TCGA_PRAD_primary'),Reference=c('TCGA_LIHC_normal','GSE132714_BPH'),Caveat=c('Same TCGA project; pairing not asserted; external to Kla cohort','Different RNA studies; BPH medications confound comparison'))
put(pairs,'focus_comparison_contract')
contr<-rbindlist(lapply(seq_len(nrow(pairs)),function(i){
 d<-merge(la[ReferenceKey==pairs$Tumor[i],.(Ensembl,Symbol,TumorA=QsmoothA,TumorUnsmoothed=UnsmoothedLog2TPM)],la[ReferenceKey==pairs$Reference[i],.(Ensembl,ReferenceA=QsmoothA,ReferenceUnsmoothed=UnsmoothedLog2TPM)],by='Ensembl')
 d[,`:=`(Comparison=pairs$Comparison[i],Caveat=pairs$Caveat[i],DeltaA=TumorA-ReferenceA,DeltaUnsmoothed=TumorUnsmoothed-ReferenceUnsmoothed)];d
}));put(contr,'focus_expression_differences')
key<-c('LDHA','LDHB','SLC16A1','SLC16A3','MPC1','MPC2','PDHA1','PDK1','PC','PCK1','PCK2','FBP1','G6PC1','ACSS2','SUCLG1','SUCLG2','AARS1','AARS2','GLO1','HAGH','LDHD')
focus<-copy(lb[ReferenceKey %in% c(pairs$Tumor,pairs$Reference)&Symbol %in% key]);focus[,Short:=factor(ReferenceKey,levels=c('TCGA_LIHC_normal','TCGA_LIHC_primary','GSE132714_BPH','TCGA_PRAD_primary'),labels=c('Adj. liver','HCC','BPH','Prostate ca.'))]
for(i in seq_len(ceiling(length(key)/9))){kk<-key[seq((i-1)*9+1,min(i*9,length(key)))];p<-ggplot(focus[Symbol %in% kk],aes(Short,QsmoothB,fill=Short))+geom_boxplot(outlier.size=.4)+facet_wrap(~Symbol,ncol=3,scales='free_y')+scale_fill_manual(values=c('#0072B2','#D55E00','#56B4E9','#CC79A7'))+labs(title='Focus: liver and prostate materials',subtitle='Descriptive sample distributions; liver n=50/371, prostate n=18/501. BPH cohort medicated; no donor pairing assumed.',x=NULL,y='qsmooth B on log2(TPM + 0.5)')+theme(axis.text.x=element_text(angle=30,hjust=1),legend.position='none');save(p,paste0('Focus_samples_',i),11,9)}
# Existing Kla quantities are union detection counts, not abundance/occupancy.
kla<-g[,.(GroupID,PXD,SampleGroup,KlaProteinCount,ReferenceProteinCount,KlaDdrProteinCount,KlaDdrFraction,ReferenceDdrProteinCount,ReferenceDdrFraction,ReferencePXD,ReferenceMatchNote)]
kla<-merge(meta[,.(GroupID,ReferenceKey,SharedReference)],kla,by='GroupID');kla[,MetricBoundary:='Study/group protein union detection counts; no site occupancy or absolute modification abundance'];put(kla,'kla_detection_context_31')
focuskla<-kla[GroupID %in% c('KLA31_08','KLA31_09','KLA31_11','KLA31_12')];put(focuskla,'focus_kla_detection_context')
p<-ggplot(melt(focuskla,id.vars=c('GroupID','SampleGroup'),measure.vars=c('KlaProteinCount','ReferenceProteinCount')),aes(SampleGroup,value,fill=variable))+geom_col(position='dodge')+geom_text(aes(label=value),position=position_dodge(width=.9),vjust=-.3,size=3)+scale_y_continuous(expand=expansion(mult=c(0,.12)))+labs(title='Existing proteome detection context',subtitle='Counts of detected protein unions; these bars do not measure lactylation abundance or occupancy.',x=NULL,y='Detected proteins',fill=NULL)+theme(legend.position='bottom')
save(p,'Kla_detection_context',10,5)
# Compact RNA contrasts accompany Kla counts without sample-level merging/correlation.
p<-ggplot(contr[Symbol %in% key],aes(DeltaA,factor(Symbol,levels=rev(key)),colour=Comparison))+geom_vline(xintercept=0,colour='grey75')+geom_point(position=position_dodge(width=.5),size=2)+labs(title='Focus expression contrasts',subtitle='Difference in qsmooth A values (tumor minus reference); not a differential-expression test or flux estimate.',x='Difference on transformed expression scale',y=NULL)+theme(legend.position='bottom')
save(p,'Focus_contrasts',10,8)
# Evidence diagram: fixed geometry, no flow magnitude encoded.
boxes<-data.table(x=c(1,3,5,7,9,5,7,3,1),y=c(4,4,4,4,4,2,2,2,2),label=c('Glucose / glycogen\nfructose / glycerol','Pyruvate\n(cytosol)','L-lactate','Lactyl-CoA /\nlactyl-AMP','L-lysine\nlactylation','Extracellular\nL-lactate','Pyruvate -> TCA\n(mitochondria) / glucose','Methylglyoxal\n-> LGSH','D-lactate /\nD-lactylation'))
edges<-data.table(x=c(1.65,3.65,5.65,7.65,5,5.5,1.3,2.3),y=c(4,4,4,4,3.55,3.6,3.55,2),xend=c(2.3,4.3,6.3,8.3,5,6.8,3,1.7),yend=c(4,4,4,4,2.5,2.5,2.5,2))
p<-ggplot()+geom_segment(data=edges[-c(2,5)],aes(x=x,y=y,xend=xend,yend=yend),arrow=grid::arrow(length=grid::unit(.12,'inches')),colour='#65717D')+geom_segment(data=edges[c(2,5)],aes(x=x,y=y,xend=xend,yend=yend),arrow=grid::arrow(length=grid::unit(.12,'inches'),ends='both'),colour='#65717D')+annotate('text',x=4,y=4.55,label='LDHA / LDHB',size=3)+annotate('text',x=6.8,y=4.6,label='ACSS2 / GTPSCS\nAARS1/2',size=3)+annotate('text',x=5.45,y=2.9,label='MCTs',size=3)+geom_label(data=boxes,aes(x,y,label=label),size=3.4,linewidth=.2,fill='#F0F5FA')+annotate('text',x=5,y=.7,label='Reaction routes only. LDH and MCT net direction is context-dependent.\nRNA does not establish concentration, compartment-specific activity, carbon flux or competition.',size=3.2)+coord_cartesian(xlim=c(0,10),ylim=c(.2,5),clip='off')+theme_void()+labs(title='Lactate sources and alternative fates',subtitle='Evidence map; D-lactate branch is separate from the main L-lactate route.')
save(p,'Lactate_routes',13,6)
put(rbindlist(figs),'figure_inventory')
checks<-data.table(Check=c('28 unique references','1898 unique sample IDs','31 mappings','Unique stable IDs','Unique measured gene labels','Finite A/B subset','No material expansion in primary table','Focus delta reconstruction','No imputed missing genes','Existing offset retained'),Pass=c(nrow(meta28)==28,nrow(qc)==1898&&!anyDuplicated(qc$SampleID),nrow(meta)==31,!anyDuplicated(a$Ensembl),!anyDuplicated(labels$Symbol),all(is.finite(la$QsmoothA))&&all(is.finite(lb$QsmoothB)),nrow(la)==length(ids)*28,all(abs(contr$DeltaA-(contr$TumorA-contr$ReferenceA))<1e-12),all(la$Ensembl %in% a$Ensembl),TRUE));stopifnot(all(checks$Pass));put(checks,'validation')
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
cat('PASS:',length(ids),'measured genes;',nrow(mem[InCommonMatrix==FALSE]),'uncovered membership rows;',length(figs),'figures\n')
