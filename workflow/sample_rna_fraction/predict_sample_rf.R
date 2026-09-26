#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table);library(ranger)})
a<-commandArgs(TRUE);stopifnot(length(a)==3)
if(file.exists(a[3]))stop('Output already exists')
m<-readRDS(a[1]);d<-fread(a[2]);stopifnot('SampleID'%in%names(d),!anyDuplicated(d$SampleID),all(m$genes%in%names(d)))
x<-as.data.frame(d[,m$genes,with=FALSE]);stopifnot(all(is.finite(as.matrix(x))))
p<-predict(m$model,x)$predictions
fwrite(data.table(SampleID=d$SampleID,PredictedMaterialKlaPercent=p),a[3])
cat('Predicted',nrow(d),'RNA records. Input must use the same corrected log2(TPM + 0.5) scale.\n')
