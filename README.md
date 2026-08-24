# Human lactylome DDR publication reproduction package

This repository reproduces only the analyses reported in the final manuscript
on human lactylation (Kla) and the DNA-damage response (DDR). It is a frozen
publication package: it does not discover new data, download public
repositories, or rerun exploratory analyses.

All biological joins, deduplication and set membership use isoform-stripped
human UniProt `BaseAccession`. Gene symbols and protein names are display-only
annotations.

## Fixed study scope

The input boundary starts with exactly 30 pan-Kla-enriched groups and their
whole-proteome references.

| Biological category | Groups |
| --- | ---: |
| Non-tumor tissues | 9 |
| Tumor tissues | 2 |
| Cancer cell lines | 12 |
| Normal cell lines | 7 |

The Kla-DDR union contains 399 `BaseAccession` values. The curated seven-pathway
score table uses `-1`, `0` and `+1` states for BER, NER, MMR, FA, HR, A-EJ and
NHEJ. Its four reported set sizes are 183, 178, 381 and 292 proteins.

## Exact output contract

Each publication run first removes only `results/figures/` and
`results/supplementary/`, then creates the following files and no others. Every
panel is exported as both PNG and PDF. Tables S1–S3 and S6 are rebuilt as XLSX;
the author-approved S4 ranking workbook and S5 regulator workbook are copied
byte-for-byte from the frozen input boundary.

| Manuscript item | Rebuilt output |
| --- | --- |
| Figure 1 | Kla and whole-proteome DDR fractions |
| Figure 2a–b | Whole-proteome DDR and Kla-DDR four-category Venn diagrams |
| Figure 2c–e | Tissue pathway matrices and tumor/non-tumor pathway summaries |
| Figure 3a–b | Reference-proteome and Kla regulator percentile heatmaps |
| Figure S1a–b | Whole-proteome and Kla four-category Venn diagrams |
| Figure S2a–c | Cell-line pathway matrices and cancer/normal cell-line summaries |
| Tables S1–S3, S6 | Kla data, reference data, human DDR GO annotations and Venn membership |
| Tables S4–S5 | Frozen pathway-ranking and regulator workbooks, preserved unchanged |

The final validation checks the input manifest, 30-group scope, 9/2/12/7
category counts, 399-protein union, 183/178/381/292 pathway panels, all 15
Venn regions for each analysis, all figure filenames and the exact S1–S6
workbook layouts.

## Reproduce

Use **R 4.4.3**. The committed `renv.lock` fixes the exact package versions
used for this release, including the figure and spreadsheet writers. From the
repository root:

```bash
Rscript workflow/install_r_dependencies.R
Rscript workflow/run_pipeline.R publication
```

Outputs are written to `results/figures/` and `results/supplementary/`. A
successful run ends with `PASS: final 30-group publication workflow completed.`
The installation step restores the locked project library; it does not upgrade
packages to current CRAN releases.

### Verified rendering environment

The release was verified on **macOS 15.6 (arm64)** with R 4.4.3, the committed
lockfile and the system font **Arial Unicode MS**. This font is mandatory: the
workflow stops before creating figures if it is absent. Other operating systems
or font substitutions have not been validated and can change figure typography
or pixels; they are outside this release's exact-reproduction claim.

## Frozen input boundary and provenance

`data/publication_input/` is the complete versioned input boundary. Its
`INPUT_MANIFEST.csv` records an MD5 checksum and byte count for every frozen
file; the workflow refuses changed, missing or extra-input substitutions.

The package retains standardized final evidence, sample-level memberships,
regulator percentiles, exact Venn memberships, the author-approved S4/S5
workbooks and frozen human DDR GO annotations. Original public PXD files are
not committed because of their size. Their ProteomeXchange accessions and
source-file locators remain in `group_summary_30.csv` and Table S1.

## Repository layout

```text
data/publication_input/     frozen final inputs and integrity manifest
R/publication/              manuscript-figure generation
workflow/                   one command entry point; builds S1–S3/S6 and copies S4/S5
tests/                      executable publication contract
results/                    regenerated outputs (ignored by Git)
```

Historical scopes, raw-download parsers, exploratory embeddings, alternative
pathway analyses and superseded score sources are intentionally absent.
