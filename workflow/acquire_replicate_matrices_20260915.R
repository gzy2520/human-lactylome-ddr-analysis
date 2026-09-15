#!/usr/bin/env Rscript
args<-commandArgs(TRUE); stopifnot(length(args)==1); root<-normalizePath(args[1])
meta<-file.path(root,'metadata/replicate_sources_20260915')
gse_files <- list(
 GSE157383='GSE157383_MCF7-2_M453_M468_DMSOvsAbema_7d_RNA-seq_readcounts.txt.gz',
 GSE171750='GSE171750_rsem_count.txt.gz',GSE203551='GSE203551_rawCounts.txt.gz',
 GSE253699='GSE253699_raw_counts.txt.gz',GSE318640='GSE318640_Expression_matrix.txt.gz',
 GSE283812='GSE283812_C48_RNAseq_raw_counts.txt.gz')
urls<-character(); studies<-character()
for(g in names(gse_files)) {
 urls<-c(urls,sprintf('https://ftp.ncbi.nlm.nih.gov/geo/series/%snnn/%s/suppl/%s',sub('[0-9]{3}$','',g),g,gse_files[[g]]));studies<-c(studies,g)
}
d<-read.csv(file.path(meta,'GSE158552_samples.csv'))
urls<-c(urls,sub('^ftp:','https:',d$Supplement[d$GSM %in% paste0('GSM',4802741:4802743)]));studies<-c(studies,rep('GSE158552',3))
records<-list()
for(i in seq_along(urls)) {
 dir<-file.path(root,'raw/replicate_replacements',studies[i]);dir.create(dir,recursive=TRUE,showWarnings=FALSE)
 dest<-file.path(dir,basename(urls[i]))
 code<-0L
 if(!file.exists(dest)) {
  code<-system2('curl',c('-fLsS','--retry','2','--max-time','180','-o',shQuote(paste0(dest,'.part')),shQuote(urls[i])))
  if(code==0)stopifnot(file.rename(paste0(dest,'.part'),dest))
 }
 ok<-code==0 && file.exists(dest) && system2('gzip',c('-t',shQuote(dest)))==0
 records[[i]]<-data.frame(GSE=studies[i],URL=urls[i],File=dest,DownloadValid=ok,
   MD5=if(ok)unname(tools::md5sum(dest)) else '')
 write.csv(do.call(rbind,records),file.path(meta,'matrix_download_status.csv'),row.names=FALSE)
 if(ok) {
  con<-gzfile(dest); lines<-readLines(con,n=3,warn=FALSE);close(con)
  cat(studies[i],basename(dest),'\n',paste(substr(lines,1,2500),collapse='\n'),'\n');flush.console()
 }
}
