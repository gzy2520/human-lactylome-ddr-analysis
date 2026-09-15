#!/usr/bin/env Rscript
out <- 'audit/20260915_teacher_review'
dir.create(file.path(out,'hk2_source'),showWarnings=FALSE)
tables <- list(); evidence <- list()
for(i in 1:3) {
  gsm <- paste0('GSM',7708661+i)
  meta <- readLines(file.path(out,'evidence',paste0(gsm,'.soft.txt')),warn=FALSE)
  stopifnot(any(grepl('genotype: WT',meta)),any(grepl('treatment: Control',meta)),
            any(grepl('HK2 cells were untreated or treated with TGF',meta)))
  url <- sub('^!Sample_supplementary_file_1 = ftp:','https:',meta[startsWith(meta,'!Sample_supplementary_file_1 = ')])
  dest <- file.path(out,'hk2_source',basename(url))
  if(!file.exists(dest)) stopifnot(system2('curl',c('-fLsS','--retry','2','--max-time','60','-o',shQuote(dest),shQuote(url)))==0)
  stopifnot(system2('gzip',c('-t',shQuote(dest)))==0)
  d <- read.delim(gzfile(dest),check.names=FALSE)
  stopifnot(all(c('Counts','TPM','FPKM') %in% names(d)))
  values <- as.matrix(d[,c('Counts','TPM','FPKM')])
  stopifnot(is.numeric(values))
  tables[[gsm]] <- d
  evidence[[gsm]] <- data.frame(SampleID=gsm,Rows=nrow(d),Columns=paste(names(d),collapse=';'),
     FirstID=as.character(d[[1]][1]),EnsemblRows=sum(grepl('^ENSG[0-9]+',d[[1]])),
     DuplicatedFirstID=sum(duplicated(d[[1]])),MissingValues=sum(is.na(values)),
     NonfiniteValues=sum(!is.finite(values)),NegativeValues=sum(values<0,na.rm=TRUE),
     StableIDMappingReady=FALSE,URL=url,MD5=unname(tools::md5sum(dest)))
}
stopifnot(identical(tables[[1]]$Gene_id,tables[[2]]$Gene_id),identical(tables[[1]]$Gene_id,tables[[3]]$Gene_id))
write.csv(do.call(rbind,evidence),file.path(out,'hk2_replacement_file_qc.csv'),row.names=FALSE)
saveRDS(tables,file.path(out,'hk2_replacement_source_tables.rds'))
cat('HK2_SOURCE_FILES_INSPECTED\n');print(do.call(rbind,evidence)[,1:6])
