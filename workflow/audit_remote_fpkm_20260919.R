root <- '/home/user/gzy/kla31-rnaseq-20260914'
args <- commandArgs(TRUE)
expression_dir <- if (length(args)) args[[1L]] else file.path(root, 'outputs/20260916_expression_extraction')
source(file.path(root,'workflow/lib_kla31_expression_20260916.R'))
expr <- parse(file.path(root,'workflow/build_31_expression_matrices_20260916.R'))
for(e in expr) if(is.call(e)&&identical(e[[1]],as.name('<-'))&&as.character(e[[2]])%in%c('h_single_ensembl_fpkm','h_per_sample_rpkm'))eval(e)
ct<-read.csv(file.path(root,'config/rnaseq_expression_extraction_contract_20260916.csv'))
l<-read.delim(file.path(root,'metadata/annotation/human_gene_lengths_ensembl111.tsv'),header=FALSE); lens<-setNames(l[[2]],l[[1]])
for(key in c('GSE132714_BPH','GSE163787_TALL104')) {
 s<-as.list(ct[ct$ReferenceKey==key,]);s$selector<-strsplit(s$SelectorValue,';',fixed=TRUE)[[1]]
 f<-get(paste0('h_',s$SourceKind))(s)$fpkm
 o<-readRDS(file.path(expression_dir,'objects',paste0(key,'.rds')))
 common<-intersect(rownames(f),names(lens)); r<-f[common,,drop=FALSE]
 correct<-sweep(r,2,colSums(r),'/')*1e6
 old_rate<-sweep(r,1,lens[common]/1000,'/')
 wrong<-sweep(old_rate,2,colSums(old_rate),'/')*1e6
 if (length(args)) stopifnot(max(abs(o$tpm[common,]-correct)) < 1e-8)
 cat(key,'stored_old_formula_max_error',max(abs(o$tpm[common,]-wrong)),'correct_log2_max_change',max(abs(log2(correct+0.5)-log2(o$tpm[common,]+0.5))),'correct_log2_median_change',median(abs(log2(correct+0.5)-log2(o$tpm[common,]+0.5))),'\n')
 cat('matrix_md5',unname(tools::md5sum(file.path(expression_dir,'matrices',paste0(key,'_log2tpm.tsv.gz')))),'\n')
}
