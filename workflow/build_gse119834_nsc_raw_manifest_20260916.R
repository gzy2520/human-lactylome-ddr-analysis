#!/usr/bin/env Rscript
# Freeze selected, untreated NSC RNA-seq libraries before raw-read re-quantification.
args <- commandArgs(TRUE)
stopifnot(length(args) == 2L)
suppressPackageStartupMessages(library(data.table))

identity_file <- normalizePath(args[[1L]], mustWork=TRUE)
out_file <- args[[2L]]
dir.create(dirname(out_file), recursive=TRUE, showWarnings=FALSE)

selected <- fread(identity_file)
stopifnot(nrow(selected) == 9L, all(selected$NonPerturbationPass),
  all(selected$RNA_GSE == 'GSE119834'))
annotations <- fread('audit/20260915_replicate_resolution/GSE119834_samples.csv')
annotations <- annotations[GSM %in% selected$RNA_GSM, .(GSM, Relation, Strategy)]
out <- merge(selected, annotations, by.x='RNA_GSM', by.y='GSM', all.x=TRUE, sort=FALSE)
out[, SRA_Experiment := sub('.*SRA: https://www\\.ncbi\\.nlm\\.nih\\.gov/sra\\?term=([^ |]+).*', '\\1', Relation)]
out[, SRA_Experiment := fifelse(grepl('^SRX[0-9]+$', SRA_Experiment), SRA_Experiment, NA_character_)]
out[, RequantificationReference := 'Homo sapiens Ensembl GRCh37-compatible release, to be frozen before quantification']
out[, RawRequantificationStatus := 'PENDING_RUN_ACCESSION_AND_FASTQ']
setorder(out, RNA_GSM)
stopifnot(all(out$Strategy == 'RNA-Seq'), !anyNA(out$SRA_Experiment), !anyDuplicated(out$RNA_GSM))
setcolorder(out, c('GroupID','RNA_GSE','RNA_GSM','SRA_Experiment','RNA_CellModel','RNA_CellType',
  'NonPerturbationPass','RNA_TreatmentStatement','ExactProteomeModelMatch','IdentityClass',
  'RequantificationReference','RawRequantificationStatus','ProteomePXD','ProteomeReferencePXD'))
fwrite(out, out_file, sep='\t')
cat('GSE119834_NSC_RAW_MANIFEST_PASS: 9 selected untreated libraries with SRX accessions.\n')
