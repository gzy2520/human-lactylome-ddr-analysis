#!/usr/bin/env Rscript
# Audit the neural-stem-cell RNA candidate against the registered proteome models.
args <- commandArgs(TRUE)
stopifnot(length(args) == 1L)
suppressPackageStartupMessages(library(data.table))

out_dir <- args[[1L]]
if (dir.exists(out_dir) && length(list.files(out_dir, all.files=TRUE, no..=TRUE))) {
  stop('Output directory already exists and is nonempty: ', out_dir, call.=FALSE)
}
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

rna <- fread('audit/20260915_replicate_resolution/GSE119834_samples.csv')
rna <- rna[grepl('cell type: neural stem cells', Characteristics, fixed=TRUE)]
rna[, RNA_CellModel := sub('.*cell line: ([^|]+) .*', '\\1', Characteristics)]
rna[, RNA_CellModel := trimws(RNA_CellModel)]
setorder(rna, GSM)

prot <- fread('corrected_final_result_20260905_sample_only/Corrected_Data/source_provenance/biological_sample_count_record.csv')
prot <- prot[SampleGroup == 'neural stem cells']
stopifnot(nrow(prot) == 1L)
proteome_models <- trimws(unlist(strsplit(prot$KlaSampleIDs[[1L]], ';', fixed=TRUE)))
stopifnot(identical(sort(proteome_models), c('ENSA', 'HMP1')))

out <- rna[, .(
  GroupID='KLA31_30',
  RNA_GSE=GSE,
  RNA_GSM=GSM,
  RNA_Title=Title,
  RNA_CellModel,
  RNA_CellType='neural stem cells',
  RNA_TreatmentStatement=Treatment,
  RNA_Assay='bulk RNA-seq',
  RNA_AuthorQuantification='FPKM; UCSC RefSeq annotation; hg19/GRCh37',
  ProteomePXD='PXD070007',
  ProteomeReferencePXD='PXD069969',
  ProteomeCellModels=paste(proteome_models, collapse=';'),
  ExactProteomeModelMatch=RNA_CellModel %in% proteome_models,
  IdentityClass=fifelse(RNA_CellModel %in% proteome_models,
    'exact_model_match', 'category_matched_external_NSC_model'),
  NonPerturbationPass=grepl('not treated', Treatment, ignore.case=TRUE)
)]

stopifnot(nrow(out) == 9L, sum(out$ExactProteomeModelMatch) == 1L,
  out[ExactProteomeModelMatch == TRUE, RNA_CellModel][[1L]] == 'ENSA',
  all(out$NonPerturbationPass))
fwrite(out, file.path(out_dir, 'nsc_rna_vs_proteome_model_map.csv'))
fwrite(out, file.path(out_dir, 'nsc_rna_vs_proteome_model_map.tsv'), sep='\t')

writeLines(c(
  '# NSC identity-resolution decision (2026-09-16)', '',
  '## Established evidence', '',
  '- KLA31_30 is registered as **neural stem cells**. Its proteome comparison uses the two cell models `ENSA` and `HMP1` (PXD070007; reference PXD069969).',
  '- GSE119834 has nine bulk-RNA-seq neural-stem-cell libraries. The source annotation states that these RNA-seq cells were not treated.',
  '- One RNA library is the same cell model as the proteome: `GSM3384849`, model `ENSA`.',
  '- The other proteome model, `HMP1`, is absent from GSE119834. The remaining eight libraries are distinct NSC models; their names and the exact match flag are in `nsc_rna_vs_proteome_model_map.csv`.', '',
  '## What this means for the current inclusion rules', '',
  'GSE119834 passes the stated material category (neural stem cells) and no-perturbation requirement. It does **not** establish a two-model ENSA/HMP1 RNA match. Its nine libraries should be described as independent external NSC models, not technical replicates of ENSA or HMP1.', '',
  '## Inclusion decision (recorded 2026-09-16)', '',
  'The user approved category-matched NSC libraries when they have no extra treatment. Therefore all nine GSE119834 NSC libraries are selected as the RNA reference for KLA31_30.',
  'The next gate is source-matrix QC: validate the downloaded author FPKM matrix, retain its RefSeq stable identifiers, verify all nine selected columns and finite numeric values, and then run source-local analysis. The ENSA/HMP1 distinction remains a provenance caveat only.', '',
  'No gene symbols are used as analytical join keys in either path. The model labels above are sample-identity annotations only.'
), file.path(out_dir, 'MANUAL_DECISION.md'))

cat('NSC_IDENTITY_RESOLUTION_PASS: 9 untreated NSC RNA libraries; 1 ENSA exact model match; HMP1 absent.\n')
