#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
release <- 'deliverables/20261011_lactylation_ddr'
out <- file.path(release, 'inputs')
dir.create(out, recursive=TRUE, showWarnings=FALSE)
paths <- c(
 RNA='outputs/20260919_repair_qsmooth/matrices/qsmooth_B_full_log2tpm.tsv.gz',
 Samples='outputs/20260919_expression_corrected/group_sample_qc.csv',
 Mapping='outputs/20260920_lactate_metabolism/tables/rna_proteome_mapping_31.csv',
 Proteome='corrected_final_result_20260905_sample_only/Corrected_Data/sample_only_inputs/figure1_sample_boxplot_values.csv',
 OriginalLDH='outputs/20260920_lactate_metabolism/tables/sample_gene_expression.csv.gz')
h <- fread(file.path(out,'hgnc_selected.tsv'))
stopifnot(nrow(h)==12,all(h$status=='Approved'),!anyDuplicated(h$ensembl_gene_id))
b <- fread(cmd=paste('gzip -dc',shQuote(paths[['RNA']])));setnames(b,1,'Ensembl')
stopifnot(nrow(b)==18332,ncol(b)==1899,!anyDuplicated(b$Ensembl),all(h$ensembl_gene_id %in% b$Ensembl))
sub <- b[Ensembl %in% h$ensembl_gene_id]
fwrite(sub,file.path(out,'selected_qsmooth_B.tsv.gz'),sep='\t')
qc <- fread(paths[['Samples']],select=c('SampleID','ReferenceKey'))
m <- fread(paths[['Mapping']]);p<-fread(paths[['Proteome']])
stopifnot(nrow(qc)==1898,!anyDuplicated(qc$SampleID),setequal(qc$SampleID,names(sub)[-1]),nrow(m)==31,uniqueN(m$ReferenceKey)==28,nrow(p)==272,all(p$ObservationType=='sample'))
fwrite(qc,file.path(out,'rna_samples.csv'));fwrite(m,file.path(out,'material_mapping_31.csv'));fwrite(p,file.path(out,'proteome_samples.csv'))
# Verify the retained LDHA/B data against the exact original figure input.
old <- fread(cmd=paste('gzip -dc',shQuote(paths[['OriginalLDH']])))[Symbol %in% c('LDHA','LDHB')]
new <- melt(sub[Ensembl %in% h[symbol %in% c('LDHA','LDHB'),ensembl_gene_id]],id.vars='Ensembl',variable.name='SampleID',value.name='NewValue')
v <- merge(old[,.(Ensembl,SampleID,QsmoothB)],new,by=c('Ensembl','SampleID'))
stopifnot(nrow(v)==2*1898,max(abs(v$QsmoothB-v$NewValue))<1e-12)
fwrite(data.table(Check='LDHA_LDHB_vs_original_input',Rows=nrow(v),MaxAbsoluteDifference=max(abs(v$QsmoothB-v$NewValue))),file.path(release,'tables','ldh_preservation.csv'))
hashes <- vapply(paths,function(p)strsplit(system2('shasum',c('-a','256',shQuote(p)),stdout=TRUE),' ')[[1]][1],character(1))
fwrite(data.table(Role=names(paths),SourcePath=unname(paths),SHA256=hashes),file.path(out,'original_source_manifest.csv'))
cat('12 genes; 1898 RNA samples; 28 references; 31 columns; 272 proteome observations. Original LDHA/B values preserved.\n')
