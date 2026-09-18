# 31-group RNA expression matrices (server-side extraction, 2026-09-16)

> **Frozen snapshot of 2026-09-16.** The live matrices are
> [`outputs/20260916_expression_extraction/`](../../outputs/20260916_expression_extraction/), which
> differ in two places after the 2026-09-18 HGNC-rename fix: `GSE171750_A549_untreated` is
> 23,771 genes (here 22,360) and `GSE283812_T47D_vehicle` is 19,154 (here 18,815), and their
> `IDRule` now records the HGNC fallback. Cite the `outputs/` copy for anything current; this
> directory is kept only as the pre-fix record.

Groups: 28 unique reference matrices covering all 31 proteome/Kla group rows (three rows share the HCT116 reference, two share the HK-2 reference).
Samples: 1898 in total. Genes per matrix: 18815 to 61365.

## What this is
Per-group expression profiles, not a differential-expression comparison. Each of the 31 group
rows was reduced to the expression profile of its own selected, unperturbed samples. No group
was compared against another group here; the cross-tissue comparison (qsmooth) is a separate
downstream step that runs on these matrices.

## Conventions
* Primary metric: log2(TPM + 0.5). Cleaned integer counts are archived alongside wherever the
  source actually provides counts.
* TPM is computed from counts against one shared merged-exon length table
  (Ensembl release 111, metadata/annotation/human_gene_lengths_ensembl111.tsv), so every
  count-based group is normalised the same way. Matrices built from a source that ships only
  FPKM/RPKM are converted with the same table and are marked as converted in the manifest.
* Row keys are stable Ensembl gene identifiers with the version suffix removed (ENSG...,
  ENSG..._PAR_Y for pseudoautosomal duplicates). Gene symbols are never a join key: the two
  author matrices whose rows are symbols (GSE171750, GSE283812) go through the official NCBI
  Symbol -> GeneID -> Ensembl route, and symbols with more than one official GeneID are dropped.
* Columns are named by the source's own sample identifier (GSM, GTEx SAMPID, TCGA barcode).
* Gene sets differ between matrices because the sources were quantified against different
  annotations; intersect on the stable gene ID before combining them.

## Files
* `matrices/<ReferenceKey>_log2tpm.tsv.gz` - genes x samples, the primary metric.
* `matrices/<ReferenceKey>_counts.tsv.gz` - archived integer counts, where available.
* `group_index.csv` - group row -> reference key -> sample and gene counts.
* `group_manifest.csv` - full provenance per matrix: source file, selector, value route, ID rule,
  mapping yield and source MD5.
* `group_sample_qc.csv` - per sample: library size, TPM total, detected genes, zero fraction and
  within-group pairwise Spearman range.
* `group_gene_summary.csv` - per gene per group: mean log2(TPM+0.5) and detection.
* `objects/<ReferenceKey>.rds` - per-group object with counts, TPM, log2 TPM, QC and provenance.

## Selection rule applied
Selected samples carry no experimental drug, knockdown, overexpression or infection. DMSO/vehicle
arms are marked as vehicle and are never merged with untreated samples. One group row accepted an
n = 1 descriptive reference (KLA31_17 TALL-104) where no second public unperturbed library exists;
no within-group variance is estimable for that row.
Known limitations that were not resolvable from public metadata are recorded per group in
`group_manifest.csv` (`SourceNote` and `IDNote`) and in the project ledger.

## Session

