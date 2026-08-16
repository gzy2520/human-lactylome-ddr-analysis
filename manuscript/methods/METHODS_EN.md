# Methods

## Study design

We integrated publicly available human lysine lactylation (Kla) mass
spectrometry datasets and paired each eligible tissue or cell sample group with
a non-Kla-enriched whole-proteome dataset of the same biological material.
Reference relationships ranged from the same biospecimens or matched
same-study samples to material-matched baselines and independent cohorts;
experimental-state matching was therefore recorded separately.
The Methods are organized following the reporting logic of a recent
bioinformatics study in *Genome Biology*: data and preprocessing, analytical
rules, quantitative metrics, visualization, computing environment, and
reproducibility. All procedures and wording below describe the present
proteomics study.

The analysis unit was a unique `PXD+SampleGroup`. Thirty-seven sample groups
had traceable Kla quantification. Of these, 30 were retained for paired
heatmaps, DDR fraction comparisons, and four-set Venn analyses. In display
order, the final scope comprised nine non-tumor tissues, two tumor tissues,
12 cancer cell lines, and seven normal cell lines.

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
non-PTM measurements of the same biological samples; the same study and
material; a material-matched baseline; or an independent dataset with the same
material identity and anatomical granularity. Similar tissues or related cell
types were not considered material matches. Material identity, donor/cohort
relationship, disease state, and treatment state were audited as separate
fields; a material match was not described as an experimental-state-matched
control unless the source data supported that claim.

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

Both heatmaps were ordered as non-tumor tissues, tumor tissues, cancer cell
lines, and normal cell lines, with Writer, Eraser, Writer-Eraser, and Reader
column facets. Kla rows were labelled with the biological material and
experimental condition of the Kla dataset, whereas whole-proteome rows were
labelled independently using the material and condition of the reference
dataset, followed by the corresponding PXD accession. The Kla heatmap used a
white-to-blue scale and the whole-proteome heatmap used a white-to-orange-red
scale, with darker colors indicating higher percentiles. Chinese and English
figures used identical data and axis orders.

## Four-set Venn analyses

Four category sets were generated separately for all Kla proteins, Kla-DDR
proteins, whole-proteome proteins, and whole-proteome DDR proteins. Each set was
deduplicated by BaseAccession, and every exact membership region was computed.
Four-set diagrams were drawn with `ggVennDiagram` from the exact membership
tables. Fixed schematic geometry retained all 15 logical four-set regions,
including zero-count regions; region area was not used as a quantitative
encoding. Displayed counts therefore came from exact membership rather than an
area fit. `eulerr` was retained only to export historical area-fit audit tables
and was not used for the final displayed diagrams. Each analysis exported
membership, region-count, and set-count tables for independent reconstruction.

## Direct-GO-term pathway matrices

The union of Kla-DDR proteins retained in the 30-group analysis contained 399
unique BaseAccessions. Four category-specific sets and the 399-protein union
were visualized separately. All direct UniProt BP, CC, and MF annotations for
the 399 proteins were retained, yielding 10,605 unique protein-term pairs and
2,785 unique direct GO terms. Curated GO seeds and their pathway-specific
descendants assigned each term to BER, NER, MMR, FA, HR, NHEJ, AEJ, or
`Others`. A term could map to more than one repair pathway. Broad DNA repair,
damage-response, or binding terms without pathway-specific evidence remained
`Others`.

For each protein and pathway, the score was the number of distinct direct GO
terms assigned to that pathway. The sum across the seven pathways was used
only to order proteins in ascending order within each panel; BaseAccession was
the deterministic tie-breaker. A pathway was shown in its solid pathway color
when at least one direct term was assigned and in light gray otherwise.
`Others` was retained for completeness auditing but was not plotted as an
eighth pathway. The former manual ±1 score workbook and pathway coefficients
were not used.

## Software and reproducibility

The main analyses used R 4.4.3 with dplyr 1.2.0, readr 2.2.0, tidyr 1.3.2,
ggplot2 4.0.2, ggVennDiagram 1.5.7, eulerr 7.0.4, GO.db 3.20.0,
readxl 1.4.5, and digest 0.6.39. Rule-based extraction of heterogeneous search
outputs used Python 3.12.10. Project configuration, source provenance,
intermediate tables, commands, and automated tests were retained.

Tests verified the 30-row Kla axis, 28-row unique reference axis, category
counts of 9/2/12/7 in display order, exclusion of seven groups, zero GeneSymbol
fallbacks, non-collapsed within-sample percentiles, deduplicated shared
references, and exact reconstruction of all Venn counts from membership and
region tables. The pathway contract additionally verified 399 proteins,
10,605 direct protein-GO pairs, 2,785 unique terms, 103 terms assigned to at
least one of the seven pathways, eight multi-pathway terms, exhaustive
seven-pathway/`Others` decisions, and 35 five-set pathway summary rows. UMAP,
t-SNE, PCA, and Cytoscape outputs were not regenerated for this analysis
revision. Final deliverables were recorded in a SHA256 manifest.

## Data and code availability

PXD accessions, article DOIs, repository links, and local evidence paths are
listed in `kla_and_reference_teacher_review_zh.csv`. The active code, execution
order, and required intermediate files are documented in
`PROJECT_INDEX.md` and `NEW_CHAT_PROJECT_PROMPT.md`. Final bilingual figures
and tables are available under `results/figures/` and `results/tables/`.
Archived 33-group, 507-protein, embedding, Cytoscape, and manual-score
materials are historical and are not part of the current result contract.

## Methods writing reference

Hu Y, Xie M, Li Y, et al. Benchmarking clustering, alignment, and integration
methods for spatial transcriptomics. *Genome Biology*. 2024;25:212.
doi:10.1186/s13059-024-03361-0.
