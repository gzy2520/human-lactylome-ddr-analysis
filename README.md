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
| Figure 2a–b | Whole-proteome DDR and Kla-DDR four-set UpSet plots |
| Figure 2c–e | Tissue pathway matrices and tumor/non-tumor pathway summaries |
| Figure 3a–b | Reference-proteome and Kla regulator percentile heatmaps |
| Figure S1a–b | Whole-proteome and Kla four-set UpSet plots |
| Figure S2a–c | Cell-line pathway matrices and cancer/normal cell-line summaries |
| Tables S1–S3, S6 | Kla data, reference data, human DDR GO annotations and four-set membership |
| Tables S4–S5 | Author-provided pathway-ranking and regulator workbooks, copied unchanged |

The accompanying checks confirm the 30-group scope, the 9/2/12/7 category
counts, the 399-protein union, the four pathway panels, the four-set membership
fields, the UpSet figure filenames and the S1–S6 workbook layout.

The four UpSet plots read the same four membership files used by Table S6.
They rank the exact intersections by size and keep all fifteen possible
combinations, including combinations whose count is zero. This changes only
the display of the frozen membership fields; it does not recalculate the
underlying tables.

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
Kla-protein membership, reference-protein membership and four analytical
set-membership tables. It checks the source-derived fields against the release tables before
the common figure and supplementary-table steps are run. The release tables
remain unchanged.

S4 pathway ranking, S5 regulator annotations, regulator percentiles, pathway
display settings and set-membership display annotations are author-provided inputs.
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
memberships, regulator percentiles, four-set memberships, the author-provided S4
and S5 workbooks and the frozen human DDR GO annotations. The original public
PXD files are not committed because of their size. Their accessions and
source-file locations are retained in `group_summary_30.csv` and Table S1.

## Figure 1 boxplot review plots

This branch contains a separate Figure 1 review plot under
`results/candidate/`. The approved publication figures and tables are left
untouched.

There are two isolated Figure 1 views. The category-level view uses the four
biological categories on the x axis and two adjacent boxes per category—whole
proteome and Kla. Each point is one PXD/sample-group dataset union, so it is
not inflated by the number of source samples. The original-layout view keeps
one row per publication group and modality, but retains source-resolved sample
points because that layout is intended to show the underlying observations.
Both views show the median as the dark box center line and the mean as a red
line.

The current category totals are:

| Category | Whole proteome | Kla |
| --- | ---: | ---: |
| non-tumor tissues | 9 | 9 |
| tumor tissues | 2 | 2 |
| cancer cell lines | 12 | 12 |
| normal cell lines | 7 | 7 |

Technical channels, fractions, acquisition modes and pooled measurements are
not counted as separate biological samples. If a source table contains only
one recoverable material or profile, that observation is retained rather than
inventing replicates. The whole-proteome and Kla sides use the same rule.

The candidate inputs are:

```text
data/candidate/biological_sample_count_record.csv
data/candidate/figure1_sample_boxplot_values.csv
data/candidate/figure1_sample_boxplot_source_registry.csv
data/candidate/figure1_dataset_boxplot_values.csv
```

Render the candidate directly with:

```bash
Rscript R/candidate/build_figure1_category_boxplot.R .
```

The renderer writes the refined review image to
`results/candidate/Figure_1_DDR_fraction_candidate_category_boxplot_refined.png`
and the matching PDF.

The renderer also writes `figure1_category_one_way_anova.csv`. The ANOVA is
run separately for each modality with category as the factor, and the two
omnibus p values are BH-adjusted. The original-layout renderer additionally
writes `figure1_original_boxplot_mean_median.csv`.

To rebuild the sample-level inputs from a local cache of the deposited
processed files, run:

```bash
Rscript R/candidate/prepare_sample_boxplot_inputs.R . /absolute/path/to/source_cache
Rscript R/candidate/build_figure1_category_boxplot.R .
```

The Kla-only reconciliation records two known source-scope differences
(PXD033146 and TALL-104); all other groups are checked against the frozen
group-level counts. Run the input checks with:

```bash
Rscript tests/validate_sample_boxplot_contract.R .
Rscript tests/validate_figure1_sample_boxplot_contract.R .
Rscript tests/validate_figure1_category_boxplot_contract.R .
```

This is a review candidate for the manuscript revision, not a replacement for
the approved publication output.

The category-level Figure 1 plot reports one-way ANOVA by modality across the
four categories. The original-layout view is descriptive and has no pooled
protein-level test. The pathway-summary plots use three-term two-way ANOVA
(category, direction and their interaction) on the dataset-level points.
Displayed q values are BH-adjusted within the stated review-test family.

## Optional MKI67-normalized Figure 1 views

For an additional view of the whole-proteome data, the candidate workflow can
plot Ki-67 (MKI67) relative to ACTB, TUBB or H3C1. Each point is one
source-resolved whole-proteome observation, and the three ratios are rendered
as three separate figures in the same four-category layout used for the
Figure 1 review plot.

The calculation uses the exact human UniProt accessions MKI67 P46013, ACTB
P60709, TUBB P07437 and H3C1 P68431. Values are taken from the deposited
whole-proteome tables. A source profile is left out of a particular ratio when
one of the required proteins is absent or non-positive; there is no
imputation. Protein-group values are kept as reported by the source, and the
ProCan matrix's log2 abundances are back-transformed before the ratio is
calculated. The input ratios themselves are not log-transformed; the y-axis is
shown on a log10 scale so that the observed range remains readable.

With a local source cache, prepare and render the three plots with:

~~~bash
Rscript R/candidate/prepare_figure1_mki67_ratio_inputs.R . /absolute/path/to/source_cache
Rscript R/candidate/build_figure1_mki67_ratio_boxplots.R .
~~~

The images are written separately to
results/candidate/mki67_ratio_boxplot/:

~~~
Figure_1_MKI67_over_ACTB_boxplot.png/.pdf
Figure_1_MKI67_over_TUBB_boxplot.png/.pdf
Figure_1_MKI67_over_H3C1_boxplot.png/.pdf
~~~

The current source cache produces 44 ACTB-normalized, 51 TUBB-normalized and
35 H3C1-normalized observations. The exact source values, missingness reasons,
coverage summary and statistics are recorded in
data/candidate/figure1_mki67_ratio_*.csv. For each denominator, the four broad
categories are first compared with a one-way ANOVA. Its p-value is adjusted
across the three denominators. Pairwise brackets are shown only when that
omnibus q-value is below 0.05; the brackets then use two-sided Wilcoxon
rank-sum tests with BH adjustment across the six category pairs for that
denominator. The statistics file retains both the three omnibus tests and all
18 pairwise tests. Each box also has a red mean line in addition to the dark
median line.

The source-audit contract for these three plots is checked with:

~~~bash
Rscript tests/validate_figure1_mki67_ratio_contract.R .
~~~

## DDR pathway-summary dataset-level review plots

The pathway summary is rendered as 14 isolated boxplots: seven pathways for
the lactylome (Kla) and the same seven pathways for whole proteome. Each plot
uses the four biological categories on the x axis. Within a category,
Up/positive and Down/negative are adjacent boxes; the Down/negative values are
not reflected to the opposite side of the axis. Each point is one
PXD/sample-group dataset union. The dark line is the median and the red line
is the mean.

The source-preparation script writes the dataset-level input (one point per
PXD/sample-group, two modalities and seven pathways) to:

```text
data/candidate/pathway_summary_dataset_boxplot_values.csv
```

Render all 14 candidate plots with:

```bash
Rscript R/candidate/build_ddr_pathway_summary_boxplots.R .
```

The PNG and PDF files are written under
`results/candidate/pathway_summary_dataset_boxplot/`. The directory also
contains `pathway_summary_dataset_boxplot_manifest.csv` and
`pathway_summary_two_way_anova.csv`. The original publication renderer and
its formal figure files are unchanged.

The pathway-summary significance labels use two-way ANOVA for each
dataset/pathway plot. The statistics retain the Category, Direction and
Category × Direction terms, with a single BH adjustment across the 42 terms
from 14 plots. Exact p values, adjusted q values and dataset point counts can
be regenerated from the same inputs by sourcing
`R/candidate/boxplot_significance.R`.

## Isolated PXD064038 ESCC inclusion scope

For the teacher-directed scope update, only the representative PXD064038
group `MEC and NEC ESCC groups` is added. MEC and NEC are the two ESCC tumor
subgroups defined by lymph-node metastasis status; PXD048995 and PXD063945
remain excluded. The six PXD064038 observations contribute 1,239 Kla
BaseAccessions and 92 Kla-DDR BaseAccessions. Its non-lactylated control is
the ordinary whole-proteome tumor arm of iProX/PXD065830: only the 94 `T`
ESCC tumor columns from Dataset1 sheet `2.a protein raw information` are
used, giving 8,083 reference BaseAccessions and 420 DDR BaseAccessions. The
24 `N` non-tumor columns are excluded. This is an independent ESCC tumor
cohort, so the comparison is by BaseAccession and is not a same-specimen
paired measurement; the healthy-esophagus PXD010154 file remains
background-only. PXD053809 is not selected for this reproducible run because
its deposited processed XLS is encrypted and is not currently auditable.

The update is isolated under:

```text
data/candidate/escc_inclusion_20260903_pxd065830_tumor_reference/
results/candidate/escc_inclusion_20260903_pxd065830_tumor_reference/
results/escc_inclusion_20260903_pxd065830_tumor_reference/formal_figures/
```

Rebuild the scope and figures with:

```bash
Rscript R/candidate/prepare_escc_inclusion_inputs.R . /absolute/path/to/source_cache
KLA_PUBLICATION_INPUT=data/candidate/escc_inclusion_20260903_pxd065830_tumor_reference/publication_input \
KLA_PUBLICATION_OUTPUT=results/escc_inclusion_20260903_pxd065830_tumor_reference/formal_figures \
KLA_PUBLICATION_EXPECTED_GROUPS=31 \
KLA_PUBLICATION_CATEGORY_COUNTS='normal_tissue=9;cancer_tissue=3;cancer_cells=12;normal_cells=7' \
Rscript R/publication/build_publication_outputs.R .
```

The candidate Figure 1 and pathway-summary renders use the corresponding
`candidate_input/` directory. Validate the isolated scope with:

```bash
Rscript tests/validate_escc_inclusion_scope.R .
```

For the expanded boxplot candidates, use:

```bash
KLA_CANDIDATE_INPUT=data/candidate/escc_inclusion_20260903_pxd065830_tumor_reference/candidate_input \
KLA_CANDIDATE_OUTPUT=results/candidate/escc_inclusion_20260903_pxd065830_tumor_reference \
KLA_CANDIDATE_EXPECTED_GROUPS=31 \
Rscript R/candidate/build_figure1_category_boxplot.R .
KLA_CANDIDATE_INPUT=data/candidate/escc_inclusion_20260903_pxd065830_tumor_reference/candidate_input \
KLA_PUBLICATION_INPUT=data/candidate/escc_inclusion_20260903_pxd065830_tumor_reference/publication_input \
KLA_CANDIDATE_OUTPUT=results/candidate/escc_inclusion_20260903_pxd065830_tumor_reference/pathway_summary_dataset_boxplot \
Rscript R/candidate/build_ddr_pathway_summary_boxplots.R .
KLA_CANDIDATE_INPUT=data/candidate/escc_inclusion_20260903_pxd065830_tumor_reference/candidate_input \
KLA_CANDIDATE_OUTPUT=results/candidate/escc_inclusion_20260903_pxd065830_tumor_reference \
KLA_CANDIDATE_EXPECTED_GROUPS=31 \
Rscript R/candidate/build_figure1_original_boxplot.R .
```

## Repository layout

```text
data/publication_input/     release tables, workbooks and input manifest
R/publication/              manuscript-figure generation
R/candidate/                optional sample-level Figure 1 review code
workflow/                   publication workflow and source-to-table preparation
tests/                      checks for the study scope and outputs
results/                    regenerated figures and supplementary tables
```
