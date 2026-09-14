#!/usr/bin/env Rscript
# Read-only inspection of registered source files. No analysis output is regenerated.
suppressPackageStartupMessages({library(data.table); library(readxl); library(digest)})
set.seed(25)
args <- commandArgs(trailingOnly=TRUE)
source_root <- if(length(args)) args[1] else '/Users/gzy2520/Desktop/Research/kla'
out <- 'outputs/20260914_reference_material_audit'
dir.create(out, recursive=TRUE, showWarnings=FALSE)
g <- fread('audit/20260912_workspace_rnaseq/groups_31.csv')
stopifnot(nrow(g)==31, !anyDuplicated(g$GroupID))
paths <- unique(g$ReferenceEvidenceFile)
inventory <- rbindlist(lapply(paths, function(p) {
  f <- file.path(source_root,p)
  data.table(Path=p, Exists=file.exists(f), Bytes=if(file.exists(f)) file.info(f)$size else NA_real_,
             SHA256=if(file.exists(f)) digest(file=f,algo='sha256') else NA_character_)
}))
fwrite(inventory,file.path(out,'reference_source_files.csv'))
headers <- list()
for(p in paths) {
  f <- file.path(source_root,p)
  if(!file.exists(f)) next
  if(grepl('xlsx$',p)) {
    for(s in excel_sheets(f)) {
      x <- suppressMessages(read_excel(f,sheet=s,n_max=0,.name_repair='minimal'))
      headers[[length(headers)+1]] <- data.table(Path=p,Part=s,Field=names(x))
    }
  } else if(grepl('zip$',p)) {
    z <- unzip(f,list=TRUE)
    headers[[length(headers)+1]] <- data.table(Path=p,Part='archive members',Field=z$Name)
  } else {
    x <- fread(f,nrows=0,showProgress=FALSE,skip=if(grepl('PXD072220',p)) 2L else 0L)
    headers[[length(headers)+1]] <- data.table(Path=p,Part='table header',Field=names(x))
  }
}
fwrite(rbindlist(headers),file.path(out,'source_headers.csv'))
checks <- list()
check <- function(label,ok,evidence) {
  checks[[length(checks)+1]] <<- data.table(Check=label,Pass=isTRUE(ok),Evidence=evidence)
}
base <- function(p) file.path(source_root,'data',p)
f <- base('PXD065775/search_results/20170330_01-24_patients_iTRAQ.xlsx')
check('Liver workbook distinguishes ANTs / DNTs / CISs',setequal(excel_sheets(f),c('ANTs','DNTs','CISs')),paste(excel_sheets(f),collapse=';'))
for(s in excel_sheets(f)) {
  x <- as.data.table(read_excel(f,sheet=s))
  checks[[length(checks)+1]] <- data.table(Check=paste('Liver sheet profile',s),Pass=NA,
    Evidence=paste('rows',nrow(x),'p<0.05',sum(x$p<0.05,na.rm=TRUE),'Change',paste(unique(x$Change),collapse=';')))
}
x <- fread(base('PXD046800/search_results/HFX2_LFQ_QB002_Proteins.txt'),nrows=0)
h <- grep('^Abundance:',names(x),value=TRUE)
check('Scar and adjacent-skin ordinary-proteome columns separated',sum(grepl('HSP',h))==4 && sum(grepl('NSP',h))==4,paste(h,collapse=';'))
x <- fread(base('PXD050147/search_results/SIRT_proteinGroups.txt'),nrows=0)
h <- grep('^LFQ intensity',names(x),value=TRUE)
check('HepG2 protein file has nine pro replicates',length(h)==9 && all(grepl('_pro_rep[123]$',h)),paste(h,collapse=';'))
x <- fread(base('PXD030304/search_results/ProCan-DepMapSanger_protein_matrix_6692_averaged.txt'),select=1)
targets <- c('MCF7','HCT-116','A549','MDA-MB-468','T47D','RKO')
for(t in targets) {
  v <- x[[1]][grepl(paste0(';',t,'$'),x[[1]])]
  check(paste('Atlas exact cell row',t),length(v)==1,paste(v,collapse=';'))
}
x <- fread(base('PXD030304/search_results/ProCan-DepMapSanger_mapping_file_replicates.txt'))
check('HEK293T control lysate identity',sum(x$Cell_line=='Control_HEK293T_lys')>0,
      paste('Control_HEK293T_lys run entries',sum(x$Cell_line=='Control_HEK293T_lys')))
x <- fread(base('PXD072220/search_results/HK-2_Spectronaut-report_PG_Quantity.txt'),skip=2,nrows=0)
h <- grep('PG.Log2Quantity$',names(x),value=TRUE)
check('HK2 selector takes first three log2 columns',length(h)==9 && all(grepl('amostra(1|3|4)\\.raw',h[1:3])),paste(h[1:3],collapse=';'))
checks[[length(checks)+1]] <- data.table(Check='HK2 raw filename alias evidence',Pass=NA,
  Evidence='Processed headers use amostra1/3/4, public raw names use Control_1/3/4; correspondence inferred from numbering and historical registry, explicit rename manifest not located.')
fwrite(rbindlist(checks),file.path(out,'targeted_checks.csv'))
for(p in c('PXD002400/search_results/msms.zip','PXD022005/search_results/txt_proteomics.zip')) {
  con <- unz(base(p),'parameters.txt')
  pars <- readLines(con,warn=FALSE); close(con)
  writeLines(pars[grepl('FDR|Fasta file|Match between runs',pars,ignore.case=TRUE)],
             file.path(out,paste0(sub('/.*','',p),'_parameter_evidence.txt')))
}
a <- g[,.(GroupID,PXD,SampleGroup,ReferencePXD,ReferenceEvidenceFile,NamingBasis)]
a[,MaterialAssessment := 'consistent_with_registered_material']
a[,ReadinessCaveat := 'Material review only; not a new statistical or protein-inference validation.']
a[ReferencePXD=='PXD030304',ReadinessCaveat := 'Atlas technical replicates/averaged model profiles; not independent biological replicates.']
a[SampleGroup=='HEK293T',ReadinessCaveat := '401 process-control lysate run entries, not 401 biological samples; detection background only.']
a[SampleGroup=='MCF10A',ReadinessCaveat := 'Source-specific parameter: parameters.txt reports Protein FDR=1 and PSM FDR=0.01. Per project instruction, retain the source protocol without imposing a stricter protein-FDR cutoff; disclose the parameter and do not claim protein-level 1% FDR.']
a[ReferencePXD=='PXD072220',MaterialAssessment := 'cell_identity_consistent_sample_alias_pending']
a[ReferencePXD=='PXD072220',ReadinessCaveat := 'amostra1/3/4 assumed to be Control_1/3/4; explicit rename manifest not found. Both Kla groups reuse the same reference.']
a[SampleGroup=='human sperm',ReadinessCaveat := 'Source is sperm, but current table includes L and N samples; donor phenotype of L/N remains unverified. Do not label all as healthy.']
a[SampleGroup=='PC-3M',ReadinessCaveat := 'Heavy SILAC PC-3M channel; not PC-3 light channel. Aggregate rather than independent sample replication.']
a[SampleGroup=='HCT116' & PXD=='PXD053474',ReadinessCaveat := 'Same HCT116 model; whole-cell background is not matched subcellular fractionation.']
a[SampleGroup=='HepG2 WT and SIRT1 or SIRT3 KO',ReadinessCaveat := 'Nine pro-labelled protein replicates and paper protein normalization support whole-proteome arm; source raw-to-pro mapping not independently reconstructed.']
a[SampleGroup %in% c('normal human lung','normal pregnancy placenta'),ReadinessCaveat := 'Organ-matched healthy-tissue atlas; anatomical subregion/placental layer not locked.']
a[,EvidenceURL := paste0('https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=',ReferencePXD)]
fwrite(a,file.path(out,'proteome_material_review_31.csv'))
r <- fread('audit/20260912_workspace_rnaseq/rnaseq_reference_registry_31.csv')
r[,`:=`(PriorSource=Source,ReviewDate='2026-09-14',InclusionBasis='same anatomical tissue and sampling class; condition differences recorded separately')]
r[GroupID %in% c('KLA31_03','KLA31_04'),`:=`(Source='GSE181540',
  ReferenceKey=ifelse(GroupID=='KLA31_03','GSE181540_HS_input','GSE181540_adjacent_NS_input'),
  SampleSelection=ifelse(GroupID=='KLA31_03','GSM5505079;GSM5505081;GSM5505083','GSM5505085;GSM5505087;GSM5505089'),
  EvidenceStatus='material_verified',AssayAndUnits='rRNA-depleted total RNA input; no IP; PE150; reported FPKM; hg19/Ensembl annotation',
  Limitations='Paper confirms paired HS and adjacent full-thickness skin. Include only input libraries (antibody none), exclude all IP. GEO library_strategy RIP-Seq conflicts with input protocol; document exception. Paper Gallus-gallus alignment sentence conflicts with GEO human hg19; verify actual reads/alignment and complete matrix before analysis.',
  EvidenceURL='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE181540;https://doi.org/10.3389/fcell.2021.748703')]
r[GroupID=='KLA31_22',`:=`(Source='GSE235595',ReferenceKey='GSE235595_PC3M_DMSO',
  SampleSelection='GSM7506018;GSM7506019;GSM7506020',EvidenceStatus='material_verified',
  AssayAndUnits='bulk RNA-seq; hg19; processed transcript counts indexed by RefSeq NM',
  Limitations='Exact PC-3M cell identity, DMSO 4-day controls. Exclude KMI169 and KMI169Ctrl. RefSeq NM must map through versioned stable identifiers to Ensembl/Entrez; not a symbol join. Matrix not yet validated.',
  EvidenceURL='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE235595')]
r[GroupID=='KLA31_23',`:=`(Source='GSE266884',ReferenceKey='GSE266884_MES28_shMCT1Control',
  SampleSelection='GSM8255530;GSM8255531',EvidenceStatus='material_verified',
  AssayAndUnits='bulk RNA-seq; SMARTER mRNA; hg38; reported FPKM (verify quantifier units)',
  Limitations='Same paper as PXD069969/PXD070007; exact MES28 model only, not all six GSC models. shMCT1 control vector; separate shKU70 controls GSM8255534/35 are alternatives, not automatically four independent replicates. Do not use MES28 RNA for NSC. No matrix acceptance yet.',
  EvidenceURL='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE266884;https://doi.org/10.1038/s41556-025-01839-y')]
stopifnot(identical(r$GroupID,g$GroupID),all(!r$AnalysisReady))
fwrite(r,file.path(out,'rnaseq_reference_candidate_31.csv'))
cat('Registered groups:',nrow(g),'; unique reference files:',nrow(inventory),'; present:',sum(inventory$Exists),'\n')
print(rbindlist(checks)[,.(Check,Pass)])
cat('Scope: material identity, source headers and selector audit; no new protein inference or statistical validation.\n')
