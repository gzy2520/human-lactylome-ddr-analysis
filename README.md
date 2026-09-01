# Reproducing the human lactylome–DDR analysis

This repository contains the code and release inputs used for the analyses in
the final manuscript on human lysine lactylation (Kla) and the DNA-damage
response (DDR).

This repository includes the step that turns downloaded, processed source
files into the analysis tables. The release-table workflow remains available,
so a reader can choose either of the following:

1. Start with the release tables already included in the repository and
   reproduce the figures and supplementary tables.
2. Start with the processed files from the public datasets, build the same
   analysis tables, and then use the same downstream workflow.

Both choices use the same 30-group study definition and the same final figure
and supplementary-table code. Files from exploratory or superseded analyses
are not part of this release.

All joins, deduplication and set membership use the isoform-stripped human
UniProt `BaseAccession`. Gene symbols and protein names are used only as
display labels.

## Study scope

The study contains 30 Kla-enriched sample groups and their whole-proteome
references.

| Biological category | Groups |
| --- | ---: |
| Non-tumor tissues | 9 |
| Tumor tissues | 2 |
| Cancer cell lines | 12 |
| Normal cell lines | 7 |

The Kla-DDR union contains 399 `BaseAccession` values. The curated
seven-pathway score table uses `-1`, `0` and `+1` states for BER, NER, MMR,
FA, HR, A-EJ and NHEJ. The four reported pathway panels contain 183, 178,
381 and 292 proteins.

## Figures and supplementary tables

A run writes to `results/figures/` and `results/supplementary/`. Before it
starts, the workflow clears only those two output folders so that files from a
previous run cannot be mistaken for current results. Every figure is saved as
both PNG and PDF.

| Manuscript item | Output |
| --- | --- |
| Figure 1 | Kla and whole-proteome DDR fractions |
| Figure 2a–b | Whole-proteome DDR and Kla-DDR four-category Venn diagrams |
| Figure 2c–e | Tissue pathway matrices and tumor/non-tumor pathway summaries |
| Figure 3a–b | Reference-proteome and Kla regulator percentile heatmaps |
| Figure S1a–b | Whole-proteome and Kla four-category Venn diagrams |
| Figure S2a–c | Cell-line pathway matrices and cancer/normal cell-line summaries |
| Tables S1–S3, S6 | Kla data, reference data, human DDR GO annotations and Venn membership |
| Tables S4–S5 | Author-provided pathway-ranking and regulator workbooks, copied unchanged |

The accompanying checks confirm the 30-group scope, the 9/2/12/7 category
counts, the 399-protein union, the four pathway panels, the Venn membership
fields, the figure filenames and the S1–S6 workbook layout.

## Reproduce the publication outputs

Use **R 4.4.3**. The committed `renv.lock` records the package versions used
for this release. From the repository root, run:

```bash
Rscript workflow/install_r_dependencies.R
Rscript workflow/run_pipeline.R publication
```

This uses the release tables in `data/publication_input/` and produces the
figures and Tables S1–S6. A successful run ends with:

```text
PASS: final 30-group publication workflow completed.
```

The installation command restores the recorded project library; it does not
upgrade packages to the latest CRAN versions.

## Start from processed source files

The source-to-table workflow expects the processed search-result and
supplementary files for the final 30 groups. Raw spectra are not needed. Put
the files in PXD folders under a separate source directory, then run:

```bash
Rscript workflow/run_pipeline.R source /absolute/path/to/source_cache
```

The preparation step reads those files and writes the 30-group summary,
Kla-protein membership, reference-protein membership and four analytical Venn
tables. It checks the source-derived fields against the release tables before
the common figure and supplementary-table steps are run. The release tables
remain unchanged.

S4 pathway ranking, S5 regulator annotations, regulator percentiles, pathway
display settings and Venn display annotations are author-provided inputs.
They are used as supplied; the source workflow does not invent replacements
for them.

The large PXD files are not stored in Git. Most are available through the
ProteomeXchange or iProX accessions recorded in the release tables and Table
S1. A small number of article supplements are exposed only through publisher
download pages. Download those files from the cited source and place them in
their recorded PXD paths before running the source workflow.

## Rendering environment

The release was checked on **macOS 15.6 (arm64)** with R 4.4.3 and the system
font **Arial Unicode MS**. This font is required for the figure labels. Other
operating systems or font substitutions may change typography and have not
been checked for pixel-level agreement.

## Release inputs and provenance

`data/publication_input/` is the versioned input set for the publication
workflow. `INPUT_MANIFEST.csv` records the MD5 checksum and byte count of each
file, and the workflow stops if a file is missing, changed or unexpectedly
added.

The release includes the standardized evidence tables, sample-level
memberships, regulator percentiles, Venn memberships, the author-provided S4
and S5 workbooks and the frozen human DDR GO annotations. The original public
PXD files are not committed because of their size. Their accessions and
source-file locations are retained in `group_summary_30.csv` and Table S1.

## Sample-level Figure 1 candidate

This branch contains an optional Figure 1 candidate for reviewing within-PXD
sample variation. It is separate from the publication workflow: it does not
change the frozen figures or tables, and it writes only to
`results/candidate/`.

The plot uses 88 source-defined observations across the same 30 publication
groups. Each point represents one sample, condition or model when that
identity can be recovered from the deposited processed file. Technical
channels, structural runs and pooled measurements are kept as transparent
single observations when independent biological sample identities cannot be
established; they are not counted as biological replicates. The `n=` label on
each row shows the number of plotted observations.

The upper block contains non-tumor and tumor tissues side by side. Cancer cell
lines and normal cell lines remain in separate rows below, with the original
group order retained. Orange boxes and points show sample-level Kla-DDR
fractions. Blue diamonds show the frozen whole-proteome group reference used
in the manuscript; they are not a second sample-level distribution.

The committed candidate tables are enough to inspect the figure directly:

```bash
Rscript R/candidate/build_figure1_sample_boxplot.R .
```

To rebuild those tables from a local cache of the deposited processed files,
run the preparation step first and then render the same figure:

```bash
Rscript R/candidate/prepare_sample_boxplot_inputs.R . /absolute/path/to/source_cache
Rscript R/candidate/build_figure1_sample_boxplot.R .
```

The preparation step also writes a reconciliation table. Every group must
match the frozen Kla and Kla-DDR protein counts before the figure is rendered.
The sample-count and provenance checks can be run with:

```bash
Rscript tests/validate_sample_boxplot_contract.R .
```

This is a review candidate for the manuscript revision, not a replacement for
the approved publication output.

## Repository layout

```text
data/publication_input/     release tables, workbooks and input manifest
R/publication/              manuscript-figure generation
R/candidate/                optional sample-level Figure 1 review code
workflow/                   publication workflow and source-to-table preparation
tests/                      checks for the study scope and outputs
results/                    regenerated figures and supplementary tables
```
