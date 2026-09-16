# Cross-material reference profile, 31 groups (2026-09-16)

Smooth quantile normalization (`qsmooth` 1.22.0) of the per-group RNA profiles, plus the DDR
annotation overlaid on the result. This is the step that turns 31 separately quantified
datasets into one profile that can be read across material classes.

## What was smoothed

| | input | shape |
|---|---|---|
| A | one median profile per reference matrix | 17,340 genes x 28 |
| B | every sample, reference matrix as group factor | 17,340 genes x 1,898 |

Both start from the 17,340 Ensembl genes present in all 28 matrices (the limit is T-47D at
18,815). 31 group rows map onto 28 matrices, because the three HCT116 rows share one matrix and
the two HK-2 rows share another.

## Result

`summary.csv` carries the numbers; the two that matter:

* **Median smoothing weight is 0.** qsmooth weights each gene between global quantile
  normalization (1) and within-group normalization (0). A median of 0 means that for most genes
  it kept the groups on their own quantiles, which is the correct behaviour when the between-group
  differences are real biology rather than technical spread - exactly the multi-tissue case
  qsmooth exists for, and the reason it is the method YARN wraps as `normalizeTissueAware`.
* **A and B agree.** Median per-gene correlation between the collapsed and the full run is
  **0.9898**, the median absolute difference is 0.21 log2 units, and no gene differs by more
  than 0.5. Collapsing to one profile per material class therefore does not materially change the
  outcome, so the simpler A form is safe to use for cross-tissue reading.

## Files

* `matrices/unified_full_log2tpm.tsv.gz`, `unified_groupmedian_log2tpm.tsv.gz` - the two inputs.
* `matrices/qsmooth_A_collapsed_log2tpm.tsv.gz` - **the reference profile to read across groups**
  (17,340 x 28, one column per reference matrix).
* `matrices/qsmooth_B_full_log2tpm.tsv.gz` - the same smoothing with within-class variation kept.
* `qsmooth_weights_A.csv.gz`, `qsmooth_weights_B.csv.gz` - one weight per gene.
* `AB_comparison.csv`, `AB_delta_by_gene_and_group.csv.gz` - the A-versus-B cross-check.
* `group_expansion_31.csv` - which group rows share a reference matrix.

The DDR overlay lives in `outputs/20260916_ddr_panel_31group/`:

* `kla_ddr_expression_31groups.tsv.gz` - the **401 Kla n DDR proteins**, 377 of which resolve to a
  gene present in the profile, across all 31 group rows.
* `reference_ddr_expression_31groups.tsv.gz` - the 836 reference DDR proteins (675 resolved).
* `ddr_go_universe_expression_31groups.tsv.gz` - the wider DDR GO universe, 2,719 accessions
  (1,950 resolved).
* `kla_ddr_annotation.csv` - the panel with its seven-pathway states and GO terms.
* `ddr_uniprot_to_ensembl.tsv` - the UniProt -> Ensembl bridge, built twice over: directly from
  UniProt's idmapping file and again through Entrez GeneID via NCBI gene2ensembl. Where both
  routes apply they agree for 879 of 882 accessions; the three exceptions are recorded, not
  silently resolved.

Gene symbols appear nowhere as a key: the DDR side is UniProt BaseAccession, the RNA side is
Ensembl gene ID, and the two are joined only through the official mappings above.

## Limitations, stated rather than buried

* **This is a normalization, not a test.** No group was tested against another. A low weight in
  some region means the groups genuinely differ there, not that something failed.
* **Platform and batch are not handled.** qsmooth corrects between-group versus within-group
  quantile behaviour; it does not correct for the fact that these 28 matrices come from roughly
  20 studies on different platforms and annotations. Its `batch` argument is unusable here
  because each group is its own batch, which is perfectly collinear with the group factor.
* **Two group rows are n = 1 by construction**: KLA31_14/15/16 are one HCT116 matrix and
  KLA31_27/28 are one HK-2 matrix. Their columns are identical in the output and are flagged in
  `group_expansion_31.csv` rather than silently duplicated.
* **24 Kla n DDR proteins are not in the profile** (unmapped accession, or the gene is absent
  from the 17,340 common set). They are listed by omission from the expression table; the count
  is in `overlay_summary.csv`.
