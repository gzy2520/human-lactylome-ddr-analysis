#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table);library(ranger)})
lab<-fread('outputs/20260928_external_escc_labels/external_labels.csv');mem<-fread('outputs/20260928_external_escc_labels/protein_members.csv')
z<-merge(lab,mem[,.(D=.N,K=sum(InKla)),by=.(Target,Material)],by=c('Target','Material'));stopifnot(all(z$Denominator==z$D),all(z$Numerator==z$K),max(abs(z$ObservedPercent-100*z$K/z$D))<1e-10)
sel<-fread('outputs/20260928_external_escc_labels/selected_columns.csv');stopifnot(!any(grepl('NAC',sel$Column)))
b<-readRDS('outputs/20260928_external_escc_rna/full_rna.rds');old<-readRDS('outputs/20260926_full_rna_inputs/full_rna.rds');stopifnot(nrow(b$x)==46,ncol(b$x)==18332,!any(b$metadata$SampleID%in%old$metadata$SampleID),identical(rownames(b$x),b$metadata$SampleID))
qc<-fread('outputs/20260928_external_escc_rna/qc.csv');stopifnot(all(qc$MissingModelGenes==0),all(qc$FractionReadsWithoutLength==0),max(abs(qc$TPMSum-1e6))<1e-6)
for(scope in c('newRNA','reused_TCGA')){
 root<-file.path('outputs/20260928_external_escc_evaluation',scope);p<-fread(file.path(root,'predictions.csv'));m<-fread(file.path(root,'summary.csv'))
 for(i in seq_len(nrow(m))){d<-p[Target==m$Target[i] & Material==m$Material[i]];stopifnot(abs(mean(abs(d$PredictedPercent-d$ObservedPercent))-m$MAE[i])<1e-9,abs(mean(d$PredictedPercent)-m$MeanPrediction[i])<1e-9)}
 for(t in c('DDR','DNA_repair')){
  rfpath<-file.path('outputs/20260927_ddr_fraction_models',t,'RF_mtry135.rds');h<-fread(file.path('outputs/20260927_ddr_fraction_models',t,'release_sha256.csv'));stopifnot(digest::digest(file=rfpath,algo='sha256')==h[Path=='RF_mtry135.rds',SHA256]);rf<-readRDS(rfpath)
  inp<-if(scope=='newRNA')b else old;d<-p[Target==t];ii<-match(d$SampleID,inp$metadata$SampleID);pr<-predict(rf$model,as.data.frame(inp$x[ii,rf$genes,drop=FALSE]))$predictions
  stopifnot(max(abs(pr-d$PredictedPercent))<1e-10)
 }
}
writeLines(c('PASS: eight material/target label rows recounted from protein membership.','PASS: all NAC columns excluded.','PASS: 46 new RNA IDs, zero overlap with frozen 1898 training records.','PASS: 18332 model genes all present; every counted gene has GENCODE19 length; TPM sums to 1e6.','PASS: both frozen model hashes unchanged; every new and reused RNA prediction reproduced.','PASS: summary MAE and prediction means independently recomputed.','SCOPE: independent RNA cohort + new proteome cohort, tissue-matched only; two shared material labels, not 46 paired individual labels.','LIMIT: clinical pretreatment history of RNA cohort not verified; proteome and Kla column counts differ (25 vs 22 per untreated tissue group).'),'outputs/20260928_external_escc_evaluation/validation.txt')
cat('External material validation PASS\n')
