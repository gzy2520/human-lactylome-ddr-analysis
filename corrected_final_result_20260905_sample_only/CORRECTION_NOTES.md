# final_result data correction — 2026-09-05

This is an isolated corrected package. The frozen package at
`/Users/gzy2520/Desktop/Research/kla/results/final_figures_and_tables` was
copied first and was not overwritten.

## Data contract applied

The latest expanded 31-group candidate tables under
`data/candidate/escc_inclusion_20260903_pxd065830_tumor_reference/candidate_input`
were used. For every boxplot/barplot point and every statistic supporting those
plots, only `ObservationType == "sample"` was retained:

| Input table | Original rows | Sample rows used | Removed non-sample rows |
|---|---:|---:|---:|
| Figure 1 DDR fraction | 310 | 272 | 38 |
| Figure 1 MKI67 ratios | 274 | 264 | 10 |
| Figure 2 seven-pathway values | 686 | 504 | 182 |

The 31-group registry is retained separately. Four registry entries have no
source-level fraction sample rows and therefore remain as empty aligned slots
in S1a; PC-3M's H3C1 aggregate is excluded from sample-only S1b and the
sample-level detected-only panel consequently contains 10 datasets.
The ratio protein key is recorded by stable UniProt accessions in
`Corrected_Data/source_provenance/figure1_mki67_ratio_protein_key.csv`.

## Re-rendered files

- Figure 1a DDR-fraction category boxplot;
- Figure 1b MKI67/H3C1 and both ACTB/TUBB companion aliases;
- Figure 2 BER, NER, MMR, FA, HR, AEJ and NHEJ barplot aliases;
- Supplementary Figure S1a, S1b and the legacy-named
  `Figure_S1b_MKI67_over_H3C1_11_detected_only` file (the title reports the
  actual 10 sample-level detected datasets).

The reviewed plotting methods, colour mapping, ordering, seed (25), and
significance functions were retained. Corrected statistical sidecars are in
`Tables/statistical_tests/`; sample-only inputs, the 31-group registry,
modality coverage, source registries and renderer manifest are in
`Corrected_Data/`; the package SHA-256 inventory is
`Corrected_Data/corrected_package_sha256.csv`.

For the four unchanged S3/S4 matrix panels, the PNGs were re-rasterized from
their intact PDFs so the protein-count digits that were missing in the old PNG
titles are visible in both delivery formats.

The whole-proteome reference reuse identified by the audit was not silently
deduplicated because it is part of the existing per-dataset comparison. The
22 affected rows/11 source keys are listed in
`Corrected_Data/whole_proteome_reference_reuse_audit.csv` and should not be
interpreted as 22 independent biological samples in downstream inference.

Reproduce to a new output directory from the workspace root with:

```sh
KLA_CORRECTED_OUTPUT="$(pwd)/corrected_final_result_20260905_rebuild" \
  Rscript --vanilla workflow/build_corrected_final_result_20260905.R "$(pwd)"
```
