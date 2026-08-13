# Methods

## Study design

We integrated publicly available human lysine lactylation (Kla) mass
spectrometry datasets and paired each eligible tissue or cell sample group with
a material- and condition-matched, non-Kla-enriched whole-proteome dataset.
The Methods are organized following the reporting logic of a recent
bioinformatics study in *Genome Biology*: data and preprocessing, analytical
rules, quantitative metrics, visualization, computing environment, and
reproducibility. All procedures and wording below describe the present
proteomics study.

The analysis unit was a unique `PXD+SampleGroup`. Thirty-seven sample groups
had traceable Kla quantification. Of these, 30 were retained for paired
heatmaps, DDR fraction comparisons, and four-set Venn analyses. The final
scope comprised nine normal/non-tumor tissues, two cancer tissues, seven
normal/non-tumor cell groups, and 12 cancer cell groups.

## Public data acquisition and organization

Lactylome datasets and article metadata were obtained from
ProteomeXchange/ProteomeCentral, PRIDE, iProX, and open article supplements.
Files for each PXD were organized into `raw`, `search_results`,
`supplementary`, and `metadata` directories. Raw spectra were retained for
provenance but were not treated as Kla evidence. The primary analysis used
author-generated database-search site, modified-peptide, and protein tables or
author supplementary tables. File-level provenance, fields, sample columns,
DOIs, and download sources were recorded in the pairing configuration and
teacher-review table.

Candidate data generated before the recognition of Kla, or enriched with a
generic anti-acyl antibody, were not classified as lactylomes based on
enrichment alone. Inclusion required explicit lactyl-lysine annotations,
modified peptides, or localized Kla sites in the search results or author
tables.

## Kla evidence extraction

Dataset-specific parsers converted heterogeneous outputs into a common
sample-level long table. Accepted evidence included MaxQuant Kla site and
modification-specific peptide tables; PEAKS Lactylation, `K(+72.02)`, or
`K(Lactyl)` annotations; Spectronaut Kla PTM or precursor quantities;
Proteome Discoverer lactylated peptide tables; and explicitly annotated Kla
sites or peptides in author supplements. The output retained the dataset,
sample, experimental group, modified peptide, site, localization probability,
PEP, score, source file, source row or site identifier, evidence mode, and
source confidence.

Samples, conditions, and replicates were not merged during parsing. Decoy,
reverse, potential-contaminant, non-human, invalid-protein, non-Kla, and
unlocalized records were excluded. Diagnostic ions were retained as supporting
evidence but could not independently classify an unmodified peptide as Kla.

No universal localization-probability threshold was imposed across
heterogeneous author outputs. Valid author-reported Kla evidence was retained
after the quality controls above. PXD014870 was excluded from the final main
analysis because the experiment was considered technically unsuccessful and
yielded only 193 Kla proteins, far fewer than the other MCF7 Kla dataset.

## Protein identifier harmonization

All analytical matching, deduplication, merging, GO intersection, and Venn
membership used human UniProt `BaseAccession`. Database prefixes and UniProt
isoform suffixes were removed. Ensembl protein identifiers in whole-proteome
sources were converted through the cached project Ensembl-to-UniProt mapping.
Reviewed status, gene symbols, and protein names were retained only for display
and manual audit and were never used for matching or missing-value fallback.
The final GeneSymbol fallback count was zero.

Repeated observations of the same protein arising from isoforms, technical
records, DDA/DIA acquisitions, or subcellular fractions were consolidated at
the traceable feature level before BaseAccession protein sets were constructed.

## Selection of whole-proteome references

Whole-proteome references were required to be non-Kla-enriched quantitative
protein matrices. References were prioritized in the following order:
non-PTM measurements of the same biological samples; the same study, material,
and treatment; or an independent dataset with exact material identity,
anatomical granularity, disease state, and treatment status. Similar tissues or
related cell types were not considered exact matches.

PXD062720, PXD063047/severe preeclampsia placenta, PXD064038, and PXD075014
were excluded from paired analyses because no exact, auditable quantitative
whole-proteome reference was available. Their source Kla data and exclusion
audits were retained. Three PXD037371 clinical liver groups were excluded from
Kla quantification because TMT channels could not be reliably assigned to the
clinical groups.

In addition, PXD055230 and PXD057709 were excluded together because they
belonged to the same herpesvirus-infection publication family, infection
strongly perturbed lactate biology, and these deep datasets dominated the
normal-cell Kla union. Their raw and processed files were retained for audit.

The 30 Kla groups mapped to 28 unique whole-proteome display rows. Pairs of
HK-2 and HCT116 Kla studies shared identical reference matrices.
Sample-level pairing remained at 30 rows, whereas the whole-proteome heatmap
displayed each identical reference once. The 30-row Kla heatmap was ordered by
the corresponding whole-proteome reference axis.

## GO-DDR annotation and fractions

The human `GO-repair+damage(human).tsv` table defined the DNA damage response
(DDR) set. Annotations with the `NOT` qualifier were excluded. GO evidence
codes were retained but were not used to restrict the primary analysis.
Protein accessions in the GO table were normalized to BaseAccession.

For every sample group, Kla and whole-proteome protein sets were intersected
with the DDR set using BaseAccession only. The DDR fraction was calculated as:

```text
DDR fraction (%) = 100 x number of unique DDR-intersecting BaseAccessions
                         / number of all unique BaseAccessions in that sample group.
```

Kla and whole-proteome denominators were calculated independently. Updating a
reference dataset therefore changed only the whole-proteome counts and
fractions.

## Lactylation regulator annotations

Writer, Eraser, and Reader proteins were obtained from
`乳酸化调控因子_Writer-Eraser-Reader.xlsx` and mapped to human UniProt
BaseAccession. A regulator could participate in more than one role. Such
multi-role annotations were retained in separate display facets, while protein
set statistics remained deduplicated by BaseAccession.

## Within-sample Kla percentile heatmap

Raw Kla intensities are not directly comparable across studies because of
differences in enrichment, instrumentation, acquisition, search software, and
measurement level. Within each quantitative sample, signals for the same
normalized feature were summed, transformed as `log2(signal+1)`, and ranked
among all quantitative Kla features. Percentiles were calculated as:

```text
percentile = 100 x (average rank - 1) / (number of features - 1).
```

A single-feature sample was assigned 100. Regulator percentiles were extracted,
and the median was used across conditions or replicates within a sample group.
An undetected regulator in an otherwise quantitative sample was assigned 0.
Values represent relative positions within each sample and are neither
cross-study abundance estimates nor Kla fold changes.

## Whole-proteome percentile heatmap

The whole-proteome heatmap used protein-level quantities from non-Kla-enriched
matrices. Decoy, contaminant, and invalid records were excluded, and duplicate
BaseAccessions were consolidated. Raw intensities were transformed as
`log2(signal+1)`; quantities already reported on a log2 scale were ranked
directly. Percentiles were calculated among all valid proteins in each sample
before regulators were extracted. The median percentile was used across
controls or technical replicates. An undetected regulator in a usable
whole-proteome sample was assigned 0. Kla-enriched intensities were never used
as substitutes.

Both heatmaps were ordered as normal tissues, cancer tissues, normal cells, and
cancer cells, with Writer, Eraser, Writer-Eraser, and Reader column facets.
The color scale ranged from white through yellow and orange to dark red, with
warmer colors indicating higher percentiles. Chinese and English figures used
identical data and axis orders.

## Four-set Venn/Euler analyses

Four category sets were generated separately for all Kla proteins, Kla-DDR
proteins, whole-proteome proteins, and whole-proteome DDR proteins. Each set was
deduplicated by BaseAccession, and every exact membership region was computed.
Area-proportional Euler diagrams were fitted with `eulerr`. Displayed counts
came from exact membership tables; ellipse geometry was an approximate
area-proportional representation. Each analysis exported membership, region
count, and set count tables for independent reconstruction.

## Linear pathway-state matrices

The union of Kla-DDR proteins retained in the 30-group analysis contained 399
unique BaseAccessions. Four category-specific sets and the 399-protein union
were visualized separately. Seven manually curated DDR pathway states (BER,
NER, MMR, FA, HR, AEJ, and NHEJ) were retained from the score workbook.
State values of +1, -1, and 0 were displayed as solid pathway color, dark
charcoal, and light gray, respectively. A weighted score was used only to sort
proteins within each panel and did not alter pathway state or color.

## Software and reproducibility

The main analyses used R 4.4.3 with dplyr 1.2.0, readr 2.2.0, tidyr 1.3.2,
ggplot2 4.0.2, eulerr 7.0.4, readxl 1.4.5, and digest 0.6.39. Rule-based
extraction of heterogeneous search outputs used Python 3.12.10. Project
configuration, source provenance, intermediate tables, commands, and automated
tests were retained.

Tests verified the 30-row Kla axis, 28-row unique reference axis, category
counts of 9/2/7/12, exclusion of seven groups, zero GeneSymbol
fallbacks, non-collapsed within-sample percentiles, deduplicated shared
references, and exact reconstruction of all Venn counts from membership and
region tables. UMAP, t-SNE, PCA, and Cytoscape outputs were not regenerated
for this analysis revision. Final deliverables were recorded in a SHA256
manifest.

## Data and code availability

PXD accessions, article DOIs, repository links, and local evidence paths are
listed in `kla_and_reference_teacher_review_zh.csv`. The active code, execution
order, and required intermediate files are documented in
`NEW_CHAT_PROJECT_PROMPT.md`. Final bilingual figures and tables are available
under `results/figures/` and `results/tables/`.

## Methods writing reference

Hu Y, Xie M, Li Y, et al. Benchmarking clustering, alignment, and integration
methods for spatial transcriptomics. *Genome Biology*. 2024;25:212.
doi:10.1186/s13059-024-03361-0.
