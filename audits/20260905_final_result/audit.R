#!/usr/bin/env Rscript
# Read-only audit of publication inputs and final artifacts. Writes only here.
suppressPackageStartupMessages({library(data.table); library(readxl); library(digest)})
set.seed(25)
args <- commandArgs(TRUE)
out <- if(length(args)) normalizePath(args[1],mustWork=TRUE) else normalizePath("audits/20260905_final_result")
root <- "/Users/gzy2520/Desktop/Research/kla"
scripts <- "/Users/gzy2520/orca/workspaces/kla/figure_refine"
res <- file.path(root,"results")
final <- file.path(res,"final_figures_and_tables")
expanded <- file.path(root,"data/candidate/escc_inclusion_20260903_pxd065830_tumor_reference")
inputs <- c(baseline30=file.path(root,"data/candidate"),expanded31=file.path(expanded,"candidate_input"))
sha <- function(f) vapply(f,function(x) digest(file=x,algo="sha256"),character(1))
rel <- function(f,r) substring(f,nchar(r)+2L)
save_csv <- function(x,name) fwrite(x,file.path(out,paste0(name,".csv")),na="")
final_files <- list.files(final,recursive=TRUE,full.names=TRUE,all.files=TRUE)
final_files <- final_files[basename(final_files)!=".DS_Store"]
before <- data.table(File=rel(final_files,final),SHA256=sha(final_files),Bytes=file.info(final_files)$size)
save_csv(before,"final_files_sha256")
sources <- list.files(res,pattern="\\.(png|pdf)$",recursive=TRUE,full.names=TRUE)
sources <- sources[!grepl("/(final_figures_and_tables|final_result)/",sources)]
hashes <- data.table(SourceFile=rel(sources,res),SHA256=sha(sources))
matches <- merge(before[grepl("\\.(png|pdf)$",File)],hashes,by="SHA256",all.x=TRUE,allow.cartesian=TRUE)
save_csv(matches,"figure_hash_provenance")
index<-before[grepl("[.]png$",File),.(File,SHA256)]
index[,`:=`(Scope="expanded31",Unit="protein/group summary",Assessment="Review caveats",Evidence="Hash provenance + source tables + PDF text + visual inspection")]
index[grepl("^Figure_1/",File),`:=`(Scope="baseline30 subset",Unit="mixed sample/aggregate",Assessment="Stale inputs; non-sample observations")]
index[grepl("Figure_1a_",File),`:=`(Scope="baseline30",Unit="210 observations: 172 sample + 38 other",Evidence="Matches baseline candidate PNG; F=142.81; tumor n=13/6 instead of 107/12")]
index[grepl("Figure_1b_",File),`:=`(Unit="35 observations: 34 sample + 1 aggregate",Evidence="Matches baseline H3C1 PNG; q=3.14e-6; missing 45 ESCC ratios")]
index[grepl("companion_controls/Figure_1_companion_",File),`:=`(Scope="expanded31 quantifiable subset",Evidence="Matches expanded candidate; ACTB n=49, TUBB n=145")]
index[grepl("^Figure_2/.*_barplot",File),`:=`(Scope="baseline30",Unit="92 observations/direction: 66 sample + 26 other",Assessment="Stale inputs; repeated-measures ANOVA error",Evidence="Baseline candidate hash; tumor n=6 not 12; final statistics use NPoint=98")]
index[grepl("UpSet",File),`:=`(Unit="unique UniProt BaseAccession membership",Assessment="31-group scope verified; not sample plots")]
index[grepl("^Figure_3/",File),`:=`(Unit="group-level median within-sample percentile",Assessment="31 groups verified; missing/zero encoding caveat",Evidence="ESCC row present; 1488 reference/1519 Kla input rows; 58/31 zero values have detections")]
index[grepl("S1a_",File),`:=`(Unit="310 observations: 272 sample + 38 other",Assessment="31 groups verified; mixed observation units",Evidence="Expanded default source; 31 labels; 94 WP + 6 Kla ESCC points")]
index[grepl("S1b_",File),`:=`(Scope="expanded31 quantifiable subset",Unit="80 observations: 79 sample + 1 aggregate across 11 groups",Assessment="Current subset; non-sample/reused observations; ND label too broad",Evidence="80 rows; ESCC n=45; remaining 20 groups lack valid ratio")]
index[grepl("pathway_matrix",File),`:=`(Unit="unique Kla-DDR protein columns",Assessment="31-group source verified; PNG title digits missing",Evidence="Workbook membership/state checks; tumor192/non-tumor183/cancer381/normal292; PDF title intact")]
index[grepl("(legacy_tissue_summaries|companion_summaries)",File),`:=`(Unit="category union protein counts/fractions; no sample points",Assessment="31-group source verified; clarify union denominator and direction",Evidence="Hash matches expanded formal figures; signed counts independently recomputed")]
save_csv(index,"per_figure_audit")
for(f in final_files[grepl("[.]pdf$",final_files)]) {
  txt <- system2("pdftotext",c("-layout",shQuote(f),"-"),stdout=TRUE)
  writeLines(txt,file.path(out,paste0(gsub("/","__",rel(f,final)),".txt")))
}
sourcefiles <- c(list.files(inputs[1],pattern="[.]csv$",full.names=TRUE),list.files(expanded,pattern="[.](csv|xlsx)$",recursive=TRUE,full.names=TRUE),list.files(file.path(scripts,"R"),pattern="[.]R$",recursive=TRUE,full.names=TRUE))
save_csv(data.table(File=sourcefiles,SHA256=sha(sourcefiles)),"inputs_and_scripts_sha256")
source(file.path(scripts,"R/candidate/boxplot_significance.R"))
categories <- c("normal_tissue","cancer_tissue","normal_cells","cancer_cells")
pathways <- c("BER","NER","MMR","FA","HR","AEJ","NHEJ")
profiles <- list(); exceptions <- list(); summaries <- list(); tests <- list(); order_checks <- list()
for(scope in names(inputs)) {
  for(nm in c("figure1_sample_boxplot_values","figure1_pathway_summary_sample_boxplot_values","figure1_mki67_ratio_sample_values")) {
    d<-fread(file.path(inputs[[scope]],paste0(nm,".csv")))
    d[,Scope:=scope]; d[,Input:=nm]
    bycols<-intersect(c("Scope","Input","Dataset","Category","Denominator","Pathway","ObservationType"),names(d))
    profiles[[length(profiles)+1L]]<-d[,.(Rows=.N,Groups=uniqueN(paste(PXD,SampleGroup))),by=bycols]
    exceptions[[length(exceptions)+1L]]<-d[ObservationType!="sample"]
    if(nm=="figure1_sample_boxplot_values") {
      summaries[[scope]]<-d[,.(N=.N,Mean=mean(DdrFractionPercentage),Median=median(DdrFractionPercentage)),by=.(Scope,Category,Dataset)]
      a<-compute_figure1_category_omnibus_anova(d,categories); a[,Scope:=scope];tests[[length(tests)+1L]]<-a
      d[,CategoryFactor:=factor(Category,levels=categories)]; d[,DatasetFactor:=factor(Dataset)]
      a1<-summary(aov(DdrFractionPercentage~CategoryFactor*DatasetFactor,d))[[1]]
      a2<-summary(aov(DdrFractionPercentage~DatasetFactor*CategoryFactor,d))[[1]]
      order_checks[[scope]]<-data.table(Scope=scope,CategoryFirstF=a1[1,"F value"],CategoryAfterModalityF=a2[2,"F value"])
      stopifnot(max(abs(d$DdrFractionPercentage-100*d$DdrProteinCount/d$ProteinCount))<1e-8)
      save_csv(d[Dataset=="Whole proteome",.(Rows=.N,Groups=paste(unique(paste(PXD,SampleGroup)),collapse=";"),Values=paste(unique(DdrFractionPercentage),collapse=";")),by=.(ReferencePXD,SourceFile,SampleID)][Rows>1],paste0(scope,"_reused_reference_identifiers"))
    }
    if(nm=="figure1_pathway_summary_sample_boxplot_values") {
      a<-compute_pathway_sample_two_way_anova(d,categories,pathways);a[,Scope:=scope];tests[[length(tests)+1L]]<-a
      stopifnot(max(abs(d$PositiveFraction-d$PositiveProteinCount/d$KlaDdrProteinCount))<1e-10,max(abs(d$NegativeFraction-d$NegativeProteinCount/d$KlaDdrProteinCount))<1e-10)
      save_csv(d[,.(N=.N,PositiveMean=mean(PositiveFraction),NegativeMean=mean(NegativeFraction),PositiveSEM=sd(PositiveFraction)/sqrt(.N),NegativeSEM=sd(NegativeFraction)/sqrt(.N)),by=.(Category,Pathway)],paste0(scope,"_pathway_stats"))
    }
    if(nm=="figure1_mki67_ratio_sample_values") {
      a<-compute_ratio_global_significance(d,c("ACTB","TUBB","H3C1"),categories);a[,Scope:=scope];tests[[length(tests)+1L]]<-a
    }
  }
}
save_csv(rbindlist(profiles,fill=TRUE),"observation_type_counts")
save_csv(rbindlist(exceptions,fill=TRUE),"non_sample_observations")
save_csv(rbindlist(summaries),"figure1_recomputed_summary")
save_csv(rbindlist(tests,fill=TRUE),"recomputed_existing_anova_not_endorsed")
save_csv(rbindlist(order_checks),"anova_term_order_sensitivity")
pub<-file.path(expanded,"publication_input")
groups<-fread(file.path(pub,"group_summary_30.csv"))
stopifnot(nrow(groups)==31L,!anyDuplicated(groups[,.(PXD,SampleGroup)]))
save_csv(groups,"expanded_31_group_registry")
upsets<-list()
for(f in c("venn_reference_ddr.csv","venn_kla_ddr.csv","venn_reference.csv","venn_all_kla.csv")) {
 d<-fread(file.path(pub,f)); cols<-paste0("In_",c("normal_tissue","cancer_tissue","cancer_cells","normal_cells"))
 m<-sapply(d[,..cols],function(x) tolower(as.character(x))%in%c("true","1","yes"))
 mask<-as.integer(m%*%c(1,2,4,8));stopifnot(!anyDuplicated(d$BaseAccession),all(mask>0))
 u<-data.table(Input=f,Mask=1:15,IntersectionCount=tabulate(mask,15),UnionCount=nrow(d))
 for(j in seq_along(cols)) u[,(cols[j]):=sum(m[,j])]
 upsets[[f]]<-u
}
save_csv(rbindlist(upsets),"upset_recomputed_counts")
ranking<-file.path(pub,"Supplementary_Table_S4_Pathway_Protein_Ranking.xlsx")
panels<-c(TumorTissues="cancer_tissue",NonTumorTissues="normal_tissue",CancerCellLines="cancer_cells",NormalCellLines="normal_cells")
membership<-fread(file.path(pub,"venn_kla_ddr.csv"));matrix_counts<-list()
for(s in names(panels)) {
 d<-as.data.table(read_excel(ranking,sheet=s)); v<-as.matrix(d[,..pathways]); flag<-paste0("In_",panels[[s]])
 expected<-membership[tolower(as.character(get(flag)))%in%c("true","1"),BaseAccession]
 stopifnot(!anyDuplicated(d$BaseAccession),setequal(expected,d$BaseAccession),all(v%in%c(-1,0,1)))
 matrix_counts[[s]]<-data.table(Sheet=s,Proteins=nrow(d),Pathway=pathways,Positive=colSums(v==1),Negative=colSums(v== -1))
}
save_csv(rbindlist(matrix_counts),"matrix_and_union_summary_counts")
hm<-list()
for(f in c("regulator_reference_percentiles_30.csv","regulator_kla_percentiles_30.csv")) {
 d<-fread(file.path(pub,f)); col<-if(grepl("reference",f)) "WholeProteomeRelativePercentile" else "RelativeKlaPercentile"
 stopifnot(uniqueN(paste(d$PXD,d$SampleGroup))==31L,all(is.na(d[[col]])|d[[col]]>=0&d[[col]]<=100))
 hm[[f]]<-d[,.(Input=f,PXD,SampleGroup,RegulatorBaseAccession,GeneSymbol,Value=get(col),DetectedSampleCount,QuantSampleCount)]
}
hm<-rbindlist(hm); save_csv(hm[Value==0&DetectedSampleCount>0],"heatmap_detected_but_zero_percentile")
save_csv(hm[,.(Rows=.N,Groups=uniqueN(paste(PXD,SampleGroup)),MissingValue=sum(is.na(Value)),ZeroNoDetection=sum(Value==0&DetectedSampleCount==0,na.rm=TRUE),ZeroWithDetection=sum(Value==0&DetectedSampleCount>0,na.rm=TRUE)),by=Input],"heatmap_summary")
# Independently verify the new ESCC ratios in the deposited workbook by accession.
rawpath<-file.path(root,"data/PXD065830/supplementary/Dataset1.xlsx")
raw<-as.data.table(suppressMessages(read_excel(rawpath,sheet="2.a protein raw information",skip=1)))
tcols<-grep("^ESCC-.*T$",names(raw),value=TRUE);stopifnot(length(tcols)==94L)
acc<-sub("-[0-9]+$","",as.character(raw[[1]]));targets<-c(MKI67="P46013",ACTB="P60709",TUBB="P07437",H3C1="P68431")
signal<-sapply(targets,function(a){ix<-which(acc==a);stopifnot(length(ix)==1L);as.numeric(raw[ix,..tcols])})
ratios<-fread(file.path(inputs[["expanded31"]],"figure1_mki67_ratio_sample_values.csv"))[PXD=="PXD064038"]
rawchecks<-rbindlist(lapply(c("ACTB","TUBB","H3C1"),function(a){
 r<-signal[,"MKI67"]/signal[,a];ok<-is.finite(r)&r>0&signal[,a]>0&signal[,"MKI67"]>0
 stored<-ratios[Denominator==a];idx<-match(stored$SampleID,tcols)
 data.table(Denominator=a,PositiveDenominator=sum(signal[,a]>0,na.rm=TRUE),ValidRatios=sum(ok),StoredRatios=nrow(stored),MaxAbsoluteRatioError=max(abs(stored$Ratio-r[idx])))
}))
stopifnot(all(rawchecks$ValidRatios==rawchecks$StoredRatios),all(rawchecks$MaxAbsoluteRatioError<1e-9))
save_csv(rawchecks,"escc_raw_ratio_verification")
# Contact sheets are a visual index only, not new scientific figures.
pngs<-before[grepl("[.]png$",File)][!duplicated(SHA256),File]
for(start in seq(1,length(pngs),by=6)) {
 batch<-pngs[start:min(start+5,length(pngs))]
 png(file.path(out,sprintf("contact_%02d.png",ceiling(start/6))),width=2400,height=2100,res=130)
 grid::grid.newpage();grid::pushViewport(grid::viewport(layout=grid::grid.layout(3,2)))
 for(i in seq_along(batch)) {
  im<-png::readPNG(file.path(final,batch[i]));aspect<-dim(im)[2]/dim(im)[1]
  grid::pushViewport(grid::viewport(layout.pos.row=ceiling(i/2),layout.pos.col=(i-1)%%2+1))
  grid::grid.text(paste(strwrap(batch[i],65),collapse="\n"),y=.97,gp=grid::gpar(fontsize=8))
  grid::grid.raster(im,y=.45,width=grid::unit(min(.96,.82*aspect/(12/7)),"npc"),height=grid::unit(min(.82,.96*(12/7)/aspect),"npc"))
  grid::popViewport()
 }
 dev.off()
}
after<-sha(final_files);stopifnot(identical(unname(before$SHA256),unname(after)))
save_csv(data.table(Check="Final artifact hashes unchanged during audit",Pass=TRUE,Files=length(final_files)),"preservation_check")
capture.output(sessionInfo(),file=file.path(out,"sessionInfo.txt"))
cat("AUDIT_COMPLETE",nrow(before),"files;",length(pngs),"distinct PNGs; final artifacts unchanged\n")
