# Matched-tissue RNA reference expansion for fixed Kla protein labels

Date: 2026-09-24. This is a sensitivity analysis for the RNA-only predictor of
the detected Kla-protein fraction. The protein catalog, denominator, 28
independent-material labels, 128 stable-Ensembl-ID gene panel, and 15
source-connected held-out folds are unchanged. No RNA sample is treated as a
new protein-labeled training observation.

## Source decisions

| GEO source | Target protein material | Included RNA | Decision |
| --- | --- | ---: | --- |
| [GSE178411](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE178411) | KLA31_03 hypertrophic scar | 28 scar sections from 26 patients | Include as an independent bulk RNA reference. The GEO series summary says 30 HTS, but sample-level fields identify only 28 as `HTS`; 2 other scar samples are `Normal scar` and are excluded. Uninjured skin and acute wounds are also excluded. |
| [GSE180836](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE180836) | KLA31_01 torn rotator-cuff tendon | 31 patients, 16 traumatic and 15 degenerative tears | Include as a sensitivity reference for the same pathological tissue. The GEO metadata show human tendon RNA-seq and no experimental gene perturbation. Diabetes status is not supplied, unlike the original non-diabetic RNA reference, so the source is not an exact clinical-condition match. |
| [GSE188952](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE188952) | KLA31_03 hypertrophic scar | 0 | Exclude: FFPE and RNA Access target capture do not meet the preference for broad whole-transcriptome RNA coverage. |
| [GSE178562](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE178562) | KLA31_03 hypertrophic scar | 0 | Exclude: source includes experimentally perturbed/cultured fibroblasts; it is not a clean untreated bulk-scar substitute. |

The included sources are human pathological tissues, not adjacent normal tissue
or healthy-tendon substitutes. “No experimental perturbation” is based on
public study/sample descriptions; unrecorded treatments or clinical covariates
cannot be ruled out.

## Processing and checks

- Author-provided gene counts and GEO per-sample metadata were downloaded and
  SHA-256 checked by `workflow/rna_expansion/download_sources.sh`. GSE180836's
  116 MiB original tar archive remains in ignored candidate data because it
  exceeds GitHub's ordinary file limit. Its SHA-256 is
  `a7cb6e33131c7edfaab4877c3592d6dae0c021148da889b947617d9388ae7368`.
- GSE178411 Entrez counts were mapped unambiguously to Ensembl genes using the
  saved human-only NCBI gene2ensembl snapshot; ambiguous Entrez IDs were
  discarded. GSE180836 author counts already use Ensembl gene IDs. Gene symbols
  were not analysis keys.
- Both sources use the same Ensembl 111 merged-exon gene lengths and existing
  `tpm_from_counts`/`log2_tpm` helper as the project RNA pipeline. Each library
  was normalized independently to TPM and transformed to log2(TPM + 0.5).
  GSE178411's repeated sections were summarized first within donor, then
  across 26 donors, so one donor did not gain extra weight. GSE180836 has one
  sample per subject and was summarized by material median. Traumatic and
  degenerative subset profiles were also retained.
- GSE178411: 25,375 genes after TPM calculation, all 128 model genes covered;
  minimum raw library 23,645,053 counts and 19,067 detected genes. GSE180836:
  60,422 genes after TPM calculation, all 128 model genes covered; minimum raw
  library 8,113,716 counts and 35,945 detected genes. Comparison of the 128
  genes with the existing RNA references gave Spearman correlations 0.888 and
  0.928 respectively. This is a compatibility check, not a paired proteomic
  validation.
- The source processors, independent output verifiers, and model-output
  verifier pass. `workflow/rna_expansion/compare_expanded_rna_models.R` also
  confirms that every scenario preserves the same 28 labels, 15 held-out source
  components, and no-RNA baseline predictions.

## Fixed-design held-out performance

Errors are percentage points. The RNA model is the prespecified 16-module
ridge candidate; the no-RNA comparison is the source-balanced median. Updating
the reference means replacing or averaging **one material-level RNA vector**,
not creating additional labeled rows. These exploratory scenarios were chosen
after inspecting the available datasets, so the best number should not be
presented as prospective generalization performance.

| RNA reference scenario | All16 material MAE | All16 source-macro MAE | All16 RMSE |
| --- | ---: | ---: | ---: |
| Original | 8.597 | 9.150 | 10.548 |
| New scar replaces old scar | 8.504 | 9.184 | 9.900 |
| Old/new scar equal-source mean | 8.460 | 9.146 | 9.893 |
| New tendon replaces old tendon | 8.563 | 9.119 | 10.507 |
| Old/new tendon equal-source mean | 8.611 | 9.164 | 10.575 |
| Degenerative tendon only replaces old tendon | 8.730 | 9.271 | 10.640 |
| Equal-source means for both scar and tendon | 8.470 | 9.157 | 9.908 |
| No-RNA balanced median, same in every run | 8.636 | 8.744 | 10.137 |

The combined scenario reduces material MAE by only 0.165 pp versus the no-RNA
baseline while increasing source-macro MAE by 0.413 pp. The actual scar
held-out absolute error increases from 0.549 pp with original RNA to 1.772 pp
with both references averaged; the tendon error remains about 6.3 pp. The
nested, fold-selected model also remains worse than the no-RNA baseline.
Therefore these data improve source coverage and auditability, but do not
establish a reliable RNA-only predictor. Preserve the original model as the
historical reference and treat these replacements/averages as sensitivity
analyses. Do not promote a scenario based on the smallest observed MAE.

The target itself is the fraction of *detected* Kla-bearing proteins among the
reference proteome's detected proteins. It is not biochemical Kla occupancy,
lactate concentration, or a matched donor-level measurement. To truly grow
training N, one needs new independent paired RNA and Kla/proteome labels; more
RNA samples for the same fixed protein material cannot supply them.

The source-held-out prediction plots, material-level comparison, and reusable
captions are in [FIGURES.md](FIGURES.md). They display the modest material-MAE
change alongside source-level errors and the narrow prediction range.

## Reproduction from the repository root

```bash
bash workflow/rna_expansion/download_sources.sh
bash workflow/rna_expansion/prepare_annotation.sh \
  data/candidate/rna_expansion/annotation/Homo_sapiens.GRCh38.111.gtf.gz \
  data/candidate/rna_expansion/annotation/gene2ensembl.gz \
  data/candidate/rna_expansion/annotation
Rscript workflow/rna_expansion/process_GSE178411_HTS.R outputs/20260924_GSE178411_HTS
Rscript workflow/rna_expansion/verify_GSE178411_HTS.R outputs/20260924_GSE178411_HTS
Rscript workflow/rna_expansion/process_GSE180836_tendon.R outputs/20260924_GSE180836_tendon
Rscript workflow/rna_expansion/verify_GSE180836_tendon.R outputs/20260924_GSE180836_tendon
Rscript workflow/lactylation_ml_review/15_optimize_rna_only_kla_fraction.R \
  outputs/20260924_rna_ref_combined_equal \
  outputs/20260924_GSE180836_tendon/tendon_material_profile.csv.gz \
  equal_source tendon_nondiabetic \
  outputs/20260924_GSE178411_HTS/HTS_material_profile.csv.gz \
  equal_source GSE181540_HS_input
Rscript workflow/lactylation_ml_review/16_verify_rna_only_optimization.R \
  outputs/20260924_rna_ref_combined_equal
Rscript workflow/rna_expansion/compare_expanded_rna_models.R \
  outputs/20260924_expanded_rna_model_comparison
```

The scripts refuse to overwrite nonempty output directories. Use fresh output
directories for a new run. Saved derived matrices and model outputs are under
`outputs/20260924_*`; raw source archives and large annotation downloads are
in ignored `data/candidate/rna_expansion/` and can be recovered by SHA-verified
download. The NCBI mapping source is a live snapshot, so the frozen derived
human mapping is also preserved in the repository.
