#!/usr/bin/env Rscript
# Record actual input files for directory-backed sources; never MD5 a directory itself.
args<-commandArgs(TRUE);stopifnot(length(args)==3L)
root<-normalizePath(args[[1]],mustWork=TRUE);contract<-read.csv(args[[2]],stringsAsFactors=FALSE)
out<-normalizePath(args[[3]],mustWork=TRUE)
contract<-contract[!duplicated(contract$ReferenceKey),]
gdc<-read.delim(file.path(root,'config/gdc_star_counts_manifest_20260914.tsv'),check.names=FALSE)
rows<-list()
for(i in seq_len(nrow(contract))) {
 s<-contract[i,];p<-file.path(root,s$SourcePath)
 if(s$SourceKind=='gtex_v8') {
  p<-file.path(root,'raw/gtex_v8',c('GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz',
                                  'GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt'))
 } else if(s$SourceKind=='gdc_star') {
  m<-gdc[gdc$CohortLabel==s$SelectorValue,]
  stopifnot(nrow(m)>0L)
  p<-file.path(root,'raw/gdc_star_counts',s$SelectorValue,paste0(m$FileUUID,'_',m$FileName))
 }
 stopifnot(all(file.exists(p)),!any(dir.exists(p)))
 rows[[s$ReferenceKey]]<-data.frame(ReferenceKey=s$ReferenceKey,Path=p,stringsAsFactors=FALSE)
}
x<-do.call(rbind,rows);paths<-unique(x$Path)
md5<-tools::md5sum(paths);stopifnot(!anyNA(md5))
x$MD5<-unname(md5[x$Path]);rownames(x)<-NULL
write.csv(x,file.path(out,'source_file_fingerprints.csv'),row.names=FALSE)
cat('SOURCE_FILE_FINGERPRINTS_PASS references=',length(unique(x$ReferenceKey)),
    ' unique_files=',length(paths),'\n',sep='')
