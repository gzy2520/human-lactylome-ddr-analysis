#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table);library(ranger)})
a<-commandArgs(TRUE);stopifnot(length(a)==3,!file.exists(a[3]));m<-readRDS(a[1]);b<-readRDS(a[2])
stopifnot(all(m$genes%in%colnames(b$x)),identical(rownames(b$x),b$metadata$SampleID),all(is.finite(b$x)))
p<-predict(m$model,as.data.frame(b$x[,m$genes,drop=FALSE]))$predictions
fwrite(data.table(SampleID=b$metadata$SampleID,PredictedMaterialKlaPercent=p),a[3])
