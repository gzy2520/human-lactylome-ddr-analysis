suppressPackageStartupMessages({library(data.table);library(digest)})
root <- getwd()
cat('CURRENT HEAD\n'); system('git rev-parse HEAD')
frozen <- 'corrected_final_result_20260905_sample_only'
h <- fread(file.path(frozen,'Corrected_Data/corrected_package_sha256.csv'))
p <- file.path(frozen,h$File)
cat('FROZEN missing',sum(!file.exists(p)),'of',length(p),'\n'); print(head(h$File[!file.exists(p)],8))
h2 <- fread('config/kla31_local_inputs_sha256.csv'); pp <- h2$Path
cat('LOCAL INPUT missing',sum(!file.exists(pp)),'of',length(pp),'\n')
ex <- file.exists(pp); cat('LOCAL INPUT hash mismatch',sum(vapply(pp[ex],function(f) digest(file=f,algo='sha256'), '') != h2$SHA256[ex]),'\n')
q <- fread('outputs/20260916_expression_extraction/group_sample_qc.csv')
cat('QC samples',nrow(q),'unique',uniqueN(q$SampleID),'refs',uniqueN(q$ReferenceKey),'\n')
baseline_results <- if (dir.exists('results/archive/rna_teacher_questions_before_20260919')) 'results/archive/rna_teacher_questions_before_20260919' else 'results/rna_teacher_questions'
for (d in c(baseline_results,'outputs/20260917_delivery_31group/rna')) {
 cat('DIRECTORY',d,'\n')
 for (f in c('figure1a_fixed_ddr_31groups.csv','figure1a_fixed_ddr_28materials.csv','figure1b_hallmark_proliferation_31groups.csv','figure1b_hallmark_proliferation_28materials.csv','figure1a_fixed_ddr_sample_values.csv')) {
  p <- file.path(d,f); if(!file.exists(p))next
  x<-fread(p); cat(f,'rows',nrow(x),'unique sample',if('SampleID'%in%names(x))uniqueN(x$SampleID) else NA,'unique refs',if('ReferenceKey'%in%names(x))uniqueN(x$ReferenceKey) else NA,'mtime',as.character(file.info(p)$mtime),'\n')
 }
}
d <- 'outputs/20260918_qsmooth_31group_hgnc'
rd<-function(f){x<-fread(file.path(d,'matrices',f)); m<-as.matrix(x[,-1]);rownames(m)<-x[[1]];m}
a<-rd('qsmooth_A_collapsed_log2tpm.tsv.gz'); u<-rd('unified_groupmedian_log2tpm.tsv.gz'); b<-rd('qsmooth_B_full_log2tpm.tsv.gz')
cat('MATRICES A',dim(a),'B',dim(b),'stable',all(grepl('^ENSG[0-9]+',rownames(b))),'duplicate ids',anyDuplicated(rownames(b)),'duplicate samples',anyDuplicated(colnames(b)),'finite',all(is.finite(b)),'metadata match',setequal(colnames(b),q$SampleID),'\n')
cat('A deviation from unsmoothed: max',max(abs(a-u)),'changed entries',sum(abs(a-u)>1e-10),'of',length(a),'\n')
for (w in c('A','B')) {v<-fread(file.path(d,paste0('qsmooth_weights_',w,'.csv.gz')));cat('WEIGHTS',w,'\n');print(summary(v))}
cat('AB comparison\n'); print(fread(file.path(d,'AB_comparison.csv')))
cat('GROUP CATEGORY\n'); g<-fread(file.path(d,'group_expansion_31.csv'));print(g[duplicated(ReferenceKey)|duplicated(ReferenceKey,fromLast=TRUE)])
cat('PROTEOME VALIDATION\n');print(fread('outputs/20260917_corrected_final_31group/Corrected_Data/correction_validation.csv'))
# Test local Stage2 matrices against the matrices actually fed into qsmooth.
v<-rd('unified_full_log2tpm.tsv.gz')
for (key in unique(q$ReferenceKey)) {
 x<-fread(file.path('outputs/20260916_expression_extraction/matrices',paste0(key,'_log2tpm.tsv.gz')))
 m<-as.matrix(x[match(rownames(v),x[[1]]),-1]); delta<-max(abs(m-v[,colnames(m),drop=FALSE]));
 cat('SOURCE_CHAIN',key,'max_delta',delta,'\n')
}
