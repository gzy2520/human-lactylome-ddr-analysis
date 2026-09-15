#!/usr/bin/env Rscript
# Retrieve source metadata only. Candidate selection is a subsequent explicit step.
args <- commandArgs(TRUE)
stopifnot(length(args)>=2)
out <- args[1]; series <- args[-1]
dir.create(out,recursive=TRUE,showWarnings=FALSE)
for(gse in series) {
  url <- sprintf('https://ftp.ncbi.nlm.nih.gov/geo/series/%snnn/%s/soft/%s_family.soft.gz',sub('[0-9]{3}$','',gse),gse,gse)
  file <- file.path(out,paste0(gse,'_family.soft.gz'))
  if(!file.exists(file)) {
    code <- system2('curl',c('-fLsS','--retry','2','--max-time','180','-o',shQuote(paste0(file,'.part')),shQuote(url)))
    if(code!=0) {cat('FETCH_FAILED',gse,'\n');next}
    stopifnot(file.rename(paste0(file,'.part'),file))
  }
  con <- gzfile(file,'rt'); records <- list(); current <- NULL
  repeat {
    chunk <- readLines(con,n=10000,warn=FALSE)
    if(!length(chunk))break
    for(line in chunk) {
      if(startsWith(line,'^SAMPLE = ')) {
        if(!is.null(current)) records[[length(records)+1L]] <- current
        current <- list(GSE=gse,GSM=sub('^\\^SAMPLE = ','',line),Title='',Characteristics='',Treatment='',Growth='',Extraction='',Strategy='',Processing='',Supplement='',Relation='')
      }
      if(is.null(current))next
      fields <- c('!Sample_title = '='Title','!Sample_characteristics_ch1 = '='Characteristics',
        '!Sample_treatment_protocol_ch1 = '='Treatment','!Sample_growth_protocol_ch1 = '='Growth',
        '!Sample_extract_protocol_ch1 = '='Extraction','!Sample_library_strategy = '='Strategy',
        '!Sample_data_processing = '='Processing','!Sample_supplementary_file_1 = '='Supplement',
        '!Sample_relation = '='Relation')
      for(prefix in names(fields)) if(startsWith(line,prefix)) {
        key <- fields[[prefix]]; value <- substring(line,nchar(prefix)+1L)
        current[[key]] <- paste(c(current[[key]][nzchar(current[[key]])],value),collapse=' | ')
      }
    }
  }
  close(con)
  if(!is.null(current))records[[length(records)+1L]]<-current
  stopifnot(length(records)>0)
  d <- do.call(rbind,lapply(records,as.data.frame,stringsAsFactors=FALSE))
  write.csv(d,file.path(out,paste0(gse,'_samples.csv')),row.names=FALSE)
  cat('METADATA_COMPLETE',gse,nrow(d),'\n');flush.console()
}
