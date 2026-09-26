#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
a<-commandArgs(TRUE);stopifnot(length(a)==2);old<-a[1];new<-a[2]
m0<-fread(file.path(old,'tables/metrics.csv'));m0<-m0[Model=='RF_mtry64'];m0[,Input:='128 genes']
m1<-fread(file.path(new,'tables/metrics.csv'));m1[,Input:='18332 genes']
fwrite(rbind(m0,m1),file.path(new,'tables/comparison_to_128.csv'))
s0<-fread(file.path(old,'tables/split.csv'));s1<-fread(file.path(new,'tables/split.csv'));stopifnot(identical(s0,s1))
p<-fread(file.path(new,'tables/predictions.csv.gz'))
stopifnot(all(p[Evaluation=='Within-material test',uniqueN(SampleID),by=Model]$V1==374),all(p[Evaluation!='Within-material test',uniqueN(SampleID),by=.(Model,Evaluation)]$V1==1898))
for(model in unique(p$Model)){
 z<-p[Model==model & Evaluation=='Within-material test'];mt<-m1[Model==model & Evaluation=='Within-material test' & Weighting=='Record']
 stopifnot(abs(mean(abs(z$Predicted-z$Label))-mt$MAE)<1e-10,abs(100*mean(abs(z$Predicted-z$Label)<=5)-mt$Within5)<1e-10)
}
h<-fread(file.path(new,'input_sha256.csv'));stopifnot(all(h$SHA256==vapply(h$Path,digest::digest,character(1),algo='sha256',file=TRUE)))
writeLines(c('PASS: exact same sample, donor, material, source and split assignment as 128-gene experiment','PASS: all models cover all 374 internal test records and 1,898 source-held-out records','PASS: internal test errors and within-five-point fractions independently recomputed','PASS: model input hashes'),file.path(new,'validation.txt'))
fwrite(data.table(Path='workflow/full_rna_fraction/compare.R',SHA256=digest::digest(file='workflow/full_rna_fraction/compare.R',algo='sha256')),file.path(new,'comparison_script_sha256.csv'))
fs<-sort(list.files(new,recursive=TRUE,full.names=TRUE));fs<-fs[basename(fs)!='release_sha256.csv'];fwrite(data.table(Path=substring(fs,nchar(new)+2),SHA256=vapply(fs,digest::digest,character(1),algo='sha256',file=TRUE)),file.path(new,'release_sha256.csv'))
print(rbind(m0,m1)[Evaluation!='Training fit' & Weighting=='Record',.(Input,Model,Evaluation,MAE,R2,Within5)])
