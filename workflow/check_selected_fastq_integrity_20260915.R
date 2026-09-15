#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
stopifnot(length(args)==2)
root <- normalizePath(args[1]); out <- args[2]
dir.create(out,recursive=TRUE,showWarnings=FALSE)
manifest <- read.delim(file.path(root,'metadata/ena_fastq_active_manifest.tsv'))
results <- list()
for(run in unique(manifest$SRR)) {
  rows <- manifest[manifest$SRR==run,]
  rows <- rows[order(rows$ENAFastqURL),]
  stopifnot(nrow(rows)==2)
  paths <- file.path(root,'raw/ena_fastq',rows$GSE,paste0(run,'_',basename(rows$ENAFastqURL)))
  # gzip -t traverses the entire file; structural pairing checks below sample 10k reads.
  integrity <- vapply(paths,function(p)system2('gzip',c('-t',shQuote(p))),integer(1))
  reads <- lapply(paths,function(p){con<-gzfile(p,'rt');on.exit(close(con));readLines(con,n=40000,warn=FALSE)})
  valid <- vapply(reads,function(x)length(x)>0 && length(x)%%4==0 &&
    all(startsWith(x[seq(1,length(x),4)],'@')) && all(startsWith(x[seq(3,length(x),4)],'+')) &&
    all(nchar(x[seq(2,length(x),4)])==nchar(x[seq(4,length(x),4)])),logical(1))
  names <- lapply(reads,function(x)sub('/[12]$','',sub(' .*','',x[seq(1,length(x),4)])))
  results[[run]] <- data.frame(SRR=run,GroupID=rows$GroupIDs[1],GzipR1=integrity[1],GzipR2=integrity[2],
    SampledReads=length(names[[1]]),StructurePass=all(valid),SampledPairNamesMatch=identical(names[[1]],names[[2]]),
    Scope='gzip full traversal; first 10000 read pairs structural QC only; not full sequencing QC')
  write.table(do.call(rbind,results),file.path(out,'fastq_integrity.tsv'),sep='\t',row.names=FALSE,quote=TRUE)
  cat(run,'checked\n'); flush.console()
}
stopifnot(all(vapply(results,function(x)x$GzipR1==0 && x$GzipR2==0 && x$StructurePass && x$SampledPairNamesMatch,logical(1))))
writeLines('PASS',file.path(out,'COMPLETE.txt'))
