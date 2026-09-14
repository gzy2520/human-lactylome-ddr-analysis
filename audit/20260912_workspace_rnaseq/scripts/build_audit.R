#!/usr/bin/env Rscript
# Run from repository root; read-only toward historical analysis inputs/results.
suppressPackageStartupMessages(library(data.table))
set.seed(25)
out <- 'audit/20260912_workspace_rnaseq'
src <- 'data/candidate/escc_inclusion_20260903_pxd065830_tumor_reference/publication_input/group_summary_30.csv'
g <- fread(src)
setorder(g, RowOrder)
stopifnot(nrow(g)==31L, !anyDuplicated(g[, .(PXD, SampleGroup)]),
          !anyDuplicated(g$RowOrder), sum(g$PXD=='PXD064038')==1)
# Historical RowOrder has gaps (13, 26, 27); preserve it, use separate audit ID.
g[, GroupID := sprintf('KLA31_%02d', seq_len(.N))]
fwrite(g, file.path(out, 'groups_31.csv'))
counts <- g[, .N, by=Category]
stopifnot(counts[Category=='normal_tissue',N]==9,
          counts[Category=='cancer_tissue',N]==3,
          counts[Category=='cancer_cells',N]==12,
          counts[Category=='normal_cells',N]==7)
fwrite(counts, file.path(out,'group_counts.csv'))
fwrite(g[, .(Groups=.N, GroupIDs=paste(GroupID,collapse=';')), by=ReferencePXD],
       file.path(out,'proteome_reference_reuse.csv'))
base <- dirname(dirname(src))
obs <- rbindlist(lapply(c('figure1_sample_boxplot_values.csv',
                        'figure1_pathway_summary_sample_boxplot_values.csv'), function(f) {
  x <- fread(file.path(base,'candidate_input',f))
  x[, .N, by=ObservationType][, InputFile:=f]
}))
fwrite(obs, file.path(out,'observation_unit_counts.csv'))
rna <- fread(file.path(out,'rnaseq_reference_registry_31.csv'))
stopifnot(nrow(rna)==31L, identical(rna$GroupID,g$GroupID),
          identical(rna$PXD,g$PXD), identical(rna$SampleGroup,g$SampleGroup),
          all(rna$AnalysisReady==FALSE), all(nzchar(rna$EvidenceURL)))
fwrite(rna[, .(Groups=.N, GroupIDs=paste(GroupID,collapse=';')), by=ReferenceKey],
       file.path(out,'rnaseq_reference_reuse.csv'))
writeLines(c('AUDIT_REGISTRY_PASS',
 sprintf('31 groups; %d unique Kla PXD; %d unique reference PXD',uniqueN(g$PXD),uniqueN(g$ReferencePXD)),
 'This validates registry structure, not RNA matrices or biological eligibility.'),file.path(out,'validation.txt'))
