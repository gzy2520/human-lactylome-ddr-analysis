# Prespecified exploratory comparisons for the current frozen Figure 1/2 inputs.
# Replicate rows remain on figures; inference balances sources and clusters by study.
suppressPackageStartupMessages({library(data.table);library(clubSandwich)})
set.seed(25)
out <- Sys.getenv('KLA_TEACHER_OUTPUT','outputs/20261006_teacher_figures')
dir.create(file.path(out,'statistics'),recursive=TRUE,showWarnings=FALSE)
cats <- c('normal_tissue','cancer_tissue','normal_cells','cancer_cells')
contrasts <- data.table(Control=c(cats[1],cats[1],cats[1],cats[3]),
                       Comparison=c(cats[2],cats[3],cats[4],cats[4]),
                       MainBracket=c(TRUE,FALSE,FALSE,TRUE))
frozen <- 'corrected_final_result_20260905_sample_only/Corrected_Data/sample_only_inputs'
stats_inputs <- list()
fit_contrasts <- function(d,family,panel,direction='') {
 d <- copy(d);stopifnot(all(is.finite(d$Value)),all(d$Category %in% cats),!anyNA(d$Source))
 d[,Category:=factor(Category,levels=cats)]
 # Equal total weight per source/category. Study clusters retain dependencies
 # across categories, such as adjacent tissue and tumor from the same study.
 d[,Weight:=1/.N,by=.(Source,Category)]
 stats_inputs[[length(stats_inputs)+1L]] <<- cbind(data.table(Family=family,Panel=panel,Direction=direction),d)
 m <- lm(Value~0+Category,data=d,weights=Weight)
 V <- vcovCR(m,cluster=d$Source,type='CR2',inverse_var=FALSE)
 rbindlist(lapply(seq_len(nrow(contrasts)),function(i){
  co<-contrasts$Control[i];cmp<-contrasts$Comparison[i]
  C<-matrix(0,nrow=1,ncol=length(coef(m)));C[1,match(paste0('Category',co),names(coef(m)))]<- -1
  C[1,match(paste0('Category',cmp),names(coef(m)))]<-1
  z<-as.data.table(linear_contrast(m,V,C,p_values=TRUE))
  data.table(Family=family,Panel=panel,Direction=direction,Control=co,Comparison=cmp,
   MainBracket=contrasts$MainBracket[i],Estimate=z$Est,SE=z$SE,DF=z$df,
   LowerCI=z$CI_L,UpperCI=z$CI_U,P=z$p_val,
   ControlRows=sum(d$Category==co),ComparisonRows=sum(d$Category==cmp),
   ControlSources=uniqueN(d[Category==co,Source]),ComparisonSources=uniqueN(d[Category==cmp,Source]))
 }))
}
f <- fread(file.path(frozen,'figure1_sample_boxplot_values.csv'))
f[,Source:=sub('^.*?(PXD[0-9]+).*','\\1',SourceFile,perl=TRUE)]
stopifnot(all(grepl('^PXD[0-9]+$',f$Source)))
# Same source/sample/category may be reused by multiple registry rows.
keys<-c('Dataset','SourceFile','SampleID','Category')
dups<-f[duplicated(f,by=keys)]
fwrite(dups,file.path(out,'statistics/reused_reference_rows.csv'))
stopifnot(f[,all(uniqueN(DdrFractionPercentage)==1),by=keys]$V1)
fu<-unique(f,by=keys)
r<-rbindlist(lapply(unique(fu$Dataset),function(ds){
 d<-fu[Dataset==ds,.(Source,Category,Unit=paste(SourceFile,SampleID,Category),Value=DdrFractionPercentage)]
 fit_contrasts(d,'Figure1a',ds)
}))
rna<-fread('outputs/20260919_repair_release/core/figure1b_hallmark_proliferation_28materials.csv')
meta<-fread('outputs/20260925_sample_fraction_inputs_final/sample_metadata.csv')
src<-unique(meta[,.(ReferenceKey,ConnGroup)])
stopifnot(!anyDuplicated(src$ReferenceKey))
rna<-merge(rna,src,by='ReferenceKey',all.x=TRUE)
stopifnot(nrow(rna)==28,!anyNA(rna$ConnGroup))
r<-rbind(r,fit_contrasts(rna[,.(Source=paste0('RNA_block_',ConnGroup),Category,Unit=ReferenceKey,
                              Value=ProliferationMedian)],'Figure1b','RNA proliferation'))
p<-fread(file.path(frozen,'figure1_pathway_summary_sample_boxplot_values.csv'))
p[,Source:=sub('^.*?(PXD[0-9]+).*','\\1',SourceFile,perl=TRUE)]
stopifnot(nrow(p)==504,!anyDuplicated(p,by=c('SourceFile','SampleID','Category','Pathway')))
for(pa in c('BER','NER','MMR','FA','HR','AEJ','NHEJ'))for(di in c('Pro','Inh')){
 d<-p[Pathway==pa,.(Source,Category,Unit=paste(SourceFile,SampleID,Category),
                   Value=100*get(if(di=='Pro')'PositiveFraction'else'NegativeFraction'))]
 r<-rbind(r,fit_contrasts(d,'Figure2',pa,di))
}
r[,Status:=fifelse(is.finite(P)&is.finite(SE)&SE>1e-10,'estimable','zero variance / not estimable')]
r[Status!='estimable',P:=NA_real_]
r[,Q:=p.adjust(P,method='BH',n=.N),by=Family]
r[,Label:=fifelse(is.na(Q),'NA',fifelse(Q<1e-4,'****',fifelse(Q<.001,'***',fifelse(Q<.01,'**',fifelse(Q<.05,'*','ns')))))]
fwrite(r,file.path(out,'statistics/pairwise_source_clustered.csv'))
fwrite(rbindlist(stats_inputs),file.path(out,'statistics/inference_units.csv'))
writeLines(capture.output(sessionInfo()),file.path(out,'statistics/sessionInfo.txt'))
print(r[,.(Tests=.N,Estimable=sum(Status=='estimable'),Qbelow05=sum(Q<.05,na.rm=TRUE)),by=Family])
