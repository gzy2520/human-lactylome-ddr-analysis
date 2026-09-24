# RNA–proteome interpretation bridge, 2026-09-24

The RNA-only model was not promoted. This analysis is a descriptive, article-oriented bridge between the frozen lactate-route RNA panel and existing Kla-protein membership, with a separate within-study ordinary-proteome check for PXD066054. No accepted RNA, protein, or manuscript DOCX artifact was overwritten.

## Reproducible results

- The prior 21-gene focus panel yields 42 gene–comparison rows, each joined by Ensembl ID. Protein joins use exact UniProt accession; 18/21 prostate genes have an unambiguous exact ordinary-protein group.
- HCC versus adjacent liver: 1,775 versus 1,798 Kla-protein species, 1,668 shared (Jaccard 0.876). The fraction of DDR proteins **among detected Kla proteins** is 6.37% versus 6.73%. The RNA references are TCGA-LIHC; the Kla and ordinary-proteome cohorts are independent of them and of each other.
- Prostate cancer versus BPH: 2,118 versus 2,125 Kla-protein species, 2,116 shared (Jaccard 0.995). Raw PXD066054 site-report extraction reproduced both frozen accession sets exactly. The DDR fraction among detected Kla proteins is 6.70% versus 6.73%. `NAT1–NAT5` in this dataset refers to BPH, not healthy adjacent prostate.
- In external prostate RNA, LDHA is +1.23 and LDHB −4.69 on the log2(TPM + 0.5) expression-difference scale. Within PXD066054, their ordinary-protein median quantity ratios are +0.05 and −1.17 on the log2 scale, respectively. These are different cohorts and cannot be treated as donor-matched changes.
- The fixed 21 genes have zero binary Kla-detection changes in either tissue pair. This says nothing about site intensity or occupancy. Five LDHA site quantities have mixed median directions; six LDHB site quantities are lower in cancer. They are separately measured, non-normalized modified-site signals and are exploratory.

## Reproduction

From this worktree root, with the original repository's downloaded PXD066054 raw exports available:

```sh
R_PROFILE_USER=/dev/null Rscript workflow/rna_protein_bridge/analyze.R \
  outputs/20260924_rna_protein_bridge \
  /Users/gzy2520/Desktop/Research/kla
```

The script refuses a nonempty output directory. Use a fresh output path for any new run. The release includes input and output SHA-256 manifests. The local run passed its row, key, cohort, raw accession-set, and protein-group assertions; all 10 released files matched the output manifest. Both PDFs have one page, and both PNGs were visually inspected.

## Inference boundary

The RNA cohorts are external to the mass-spectrometry patients. The prostate RNA comparison spans TCGA-PRAD and a BPH GEO cohort with medication/background differences. The protein-quantity comparison is within PXD066054 and limited to exact single-accession groups; no cross-study protein-intensity comparison is made. Detection sets do not measure total Kla percentage, site occupancy, writer/eraser activity, lactate level, or metabolic flux. The mechanistic explanation is therefore a set of plausible control points to test, not a demonstrated balance or diversion pathway.

Candidate article text and figure legends: `manuscript/RNA_PROTEOME_BRIDGE_CANDIDATE_20260924.md`.
