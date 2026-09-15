#!/usr/bin/env Rscript
args<-commandArgs(TRUE);stopifnot(length(args)>1)
out<-args[1];dir.create(out,recursive=TRUE,showWarnings=FALSE);records<-list()
for(gse in args[-1]) {
 page<-file.path(out,paste0(gse,'_download.html'))
 url<-paste0('https://www.ncbi.nlm.nih.gov/geo/download/?acc=',gse,'&type=rnaseq_counts')
 code<-system2('curl',c('-fLsS','--retry','1','--max-time','35','-o',shQuote(page),shQuote(url)))
 html<-if(code==0)paste(readLines(page,warn=FALSE),collapse=' ') else ''
 pattern<-paste0(gse,'_raw_counts_[A-Za-z0-9_.]+_NCBI.tsv.gz')
 hits<-unique(regmatches(html,gregexpr(pattern,html))[[1]])
 state<-'not_listed_or_fetch_failed';dest<-'';nr<-nc<-NA_integer_;ids<-''
 if(length(hits)==1 && nzchar(hits[1])) {
  dest<-file.path(out,hits[1]);download<-paste0('https://www.ncbi.nlm.nih.gov/geo/download/?type=rnaseq_counts&acc=',gse,'&format=file&file=',hits[1])
  code<-if(file.exists(dest))0 else system2('curl',c('-fLsS','--retry','2','--max-time','120','-o',shQuote(paste0(dest,'.part')),shQuote(download)))
  if(code==0 && !file.exists(dest))stopifnot(file.rename(paste0(dest,'.part'),dest))
  if(code==0 && system2('gzip',c('-t',shQuote(dest)))==0) {
    d<-read.delim(gzfile(dest),check.names=FALSE)
    nr<-nrow(d);nc<-ncol(d)-1L;ids<-paste(head(as.character(d[[1]]),3),collapse=';')
    state<-'downloaded_table';writeLines(names(d),file.path(out,paste0(gse,'_columns.txt')))
  } else state<-'download_failed'
 }
 records[[gse]]<-data.frame(GSE=gse,State=state,File=dest,Genes=nr,Samples=nc,ExampleID=ids,
  MD5=if(state=='downloaded_table')unname(tools::md5sum(dest)) else '')
 write.csv(do.call(rbind,records),file.path(out,'ncbi_count_inventory.csv'),row.names=FALSE)
 cat(gse,state,nr,nc,ids,'\n');flush.console()
}
