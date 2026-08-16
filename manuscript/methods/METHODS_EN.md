# Methods

<!--
Provisional manuscript labels used in this draft:
Supplementary Table S1, dataset catalogue and inclusion decisions;
Supplementary Table S2, whole-proteome reference pairing and audit;
Supplementary Table S3, lactylation regulator annotations;
Supplementary Table S4, exact Venn membership and region counts;
Supplementary Table S5, GO-term-to-pathway decisions;
Supplementary Table S6, manually curated signed pathway annotations;
Supplementary Data S1, sample-level Kla evidence;
Supplementary Data S2, human GO-DDR annotations;
Supplementary Data S3, direct GO annotations and pathway matrices.
Table numbers, figure numbers and repository links should be updated when the
manuscript layout is finalized.
-->

## Study design and analysis scope

We performed an integrative secondary analysis of publicly available human
lysine lactylation (Kla) mass-spectrometry datasets. The analytical unit was a
unique ProteomeXchange accession and biological sample-group combination.
Each eligible Kla group was linked, where possible, to a quantitative
whole-proteome dataset representing the same biological material. Material
identity and experimental-state matching were evaluated separately because
some references represented the same biospecimens or matched samples, whereas
others represented material-matched baselines or independent cohorts.

The source catalogue comprised 40 candidate sample groups, of which 37
provided traceable protein-level Kla quantification. Thirty-three of these
groups had an acceptable non-Kla-enriched whole-proteome reference. Before the
final comparative analysis, three additional groups were excluded according to
the predefined study-scope and quality-control decisions described below,
leaving 30 Kla groups with an associated whole-proteome reference. In the
order used throughout the figures, these comprised nine non-tumor tissues, two
tumor tissues, 12 cancer cell lines, and seven normal cell lines. The complete
dataset catalogue, eligibility decisions, sample annotations and source
publications are summarized in Supplementary Table S1.

## Public proteomics data acquisition

Lactylome datasets, study metadata and associated publications were retrieved
from ProteomeXchange/ProteomeCentral, PRIDE, iProX and openly available article
supplements. Deposited database-search results and author-curated
supplementary tables were used as the primary evidence sources. Raw
mass-spectrometry files were retained for provenance but were not uniformly
re-searched and were not treated as evidence of Kla in the absence of an
explicit lactyl-lysine assignment.

Studies performed before Kla was recognized, datasets enriched with a generic
anti-acyl antibody, and non-PTM proteomes were not classified as lactylomes on
the basis of enrichment or protein detection alone. Inclusion required an
explicit Kla site, a lactylated peptide, a Kla-specific modification
annotation, or a Kla-specific quantitative feature in the deposited search
results or author-curated tables. File provenance, source fields, sample
columns, digital object identifiers and repository accessions were recorded
for every included group (Supplementary Table S1).

## Extraction and harmonization of Kla evidence

Dataset-specific parsers converted heterogeneous search outputs into a common
sample-level evidence table (Supplementary Data S1). Supported inputs included
MaxQuant modification-site and evidence tables, PEAKS peptide and
protein-peptide tables, Spectronaut PTM and precursor reports, Proteome
Discoverer modification-specific peptide tables, and author-curated Kla site
or peptide tables. Lactylation annotations included software-specific Kla or
Lac labels and lysine modifications corresponding to the lactyl mass shift.

For each retained record, the harmonized evidence table preserved the study,
sample, biological material, experimental group, replicate, modified peptide,
protein site, localization probability where available, posterior error
probability, identification score, evidence level and source record. Samples,
experimental conditions and replicates were not merged during parsing. When
quantitative columns could contain imputed values, sample-level detection was
defined using explicit identification flags rather than intensity alone.

Decoy, reverse, potential-contaminant, non-human, invalid-protein, non-Kla and
unlocalized records were excluded. Diagnostic ions were retained as supporting
information but were not sufficient to classify an otherwise unmodified
peptide as Kla. Because the contributing studies used different search
software, acquisition strategies and reporting conventions, a universal
localization-probability threshold was not imposed across all datasets.
Instead, valid author-reported Kla evidence was retained after the
source-specific quality controls documented in Supplementary Table S1.

Where both author-curated site tables and independently processed search
results were available, the integration rule was specified before constructing
protein sets. For the HCT116 subcellular dataset, author-reported sites were
retained, and additional search-only sites were admitted only when supported
independently by both data-dependent and data-independent acquisitions.
Analogous source-specific integration rules are reported in Supplementary
Table S1.

## Protein identifier harmonization

All analytical matching, deduplication, merging, Gene Ontology (GO)
intersection and set membership used human UniProt accessions. Database
prefixes were removed and UniProt isoform suffixes were collapsed to the
canonical accession stem, hereafter referred to as the BaseAccession.
Whole-proteome features reported as Ensembl protein identifiers were converted
to UniProt through an explicit cached mapping. Gene symbols, protein names and
UniProt reviewed status were retained for display and audit only and were not
used as matching fallbacks. No final protein match depended on a gene-symbol
fallback.

Repeated observations of the same BaseAccession arising from isoforms,
technical records, replicate acquisitions, DDA/DIA measurements or
subcellular fractions were consolidated only after the underlying
sample-level evidence had been retained. Protein-set counts therefore
represent unique BaseAccessions within each biological sample group.

## Selection of whole-proteome references

Whole-proteome references were required to be quantitative protein matrices
generated without Kla or other PTM enrichment. References were prioritized in
the following order: non-PTM measurements of the same biospecimens; matched
samples from the same study; a baseline dataset for the same biological
material; and an independent cohort with the same material identity and
anatomical granularity. A related tissue or cell type was not accepted as a
material match solely on the basis of biological similarity.

Material identity, donor or cohort relationship, disease state, treatment
state, source measurements and sample subsets were reviewed separately
(Supplementary Table S2). Consequently, a material-matched reference was not
described as a condition-matched control unless the source data supported that
interpretation. Kla and whole-proteome measurements were treated as parallel
protein-set references rather than as paired treatment-effect measurements
when they originated from different cohorts or experimental states.

Three clinical liver groups from PXD037371 were not included among the 37
quantifiable Kla groups because the deposited TMT channels could not be
reliably assigned to the reported clinical groups. Four otherwise
quantifiable Kla groups were excluded from paired analysis because an
acceptable quantitative whole-proteome reference was unavailable:
PXD062720, the severe-preeclampsia group from PXD063047, PXD064038 and
PXD075014. Their source data and exclusion records were retained.

Three further datasets were excluded before the final 30-group analysis.
PXD014870 was excluded because the Kla experiment yielded 193 proteins, far
fewer than the independent MCF7 lactylome included in the analysis. The
related PXD055230 and PXD057709 herpesvirus-infection datasets were excluded
together because infection substantially perturbs lactate biology and the
depth of this single publication family contributed disproportionately to the
normal-cell Kla union. These decisions and their quantitative effects are
documented in Supplementary Tables S1 and S2.

The final 30 Kla groups mapped to 28 unique whole-proteome display rows because
two pairs of Kla studies used identical reference matrices and sample subsets.
Sample-level Kla–reference linkage was retained as 30 analytical rows, whereas
a duplicated reference was displayed once in the whole-proteome heatmap. The
Kla heatmap was ordered relative to the corresponding whole-proteome reference
axis.

## Definition of the DNA damage-response protein set

The DNA damage-response (DDR) set was defined using a precompiled human
GO-derived annotation resource (Supplementary Data S2). Only human records
(taxon 9606) linked to UniProtKB accessions were retained. Annotations
containing the GO qualifier `NOT` were removed. GO evidence codes were
preserved for audit but were not used to restrict the primary analysis.
Protein identifiers in this resource were normalized to BaseAccession before
intersection with the Kla and whole-proteome sets.

## Calculation of Kla and whole-proteome DDR fractions

For each of the 30 sample groups, Kla and whole-proteome protein sets were
intersected independently with the DDR set using BaseAccession. The DDR
fraction was calculated as

$$
\mathrm{DDR\ fraction} =
\frac{\text{number of unique DDR-intersecting BaseAccessions}}
{\text{number of all unique BaseAccessions in the corresponding protein set}}.
$$

Kla and whole-proteome denominators were calculated independently. A
whole-proteome protein count was never inferred from Kla-enriched signal, and
Kla signal was never substituted for missing whole-proteome abundance. The
cross-study DDR fractions were summarized descriptively; they were not
interpreted as direct measurements of DNA damage, pathway activation or
treatment response.

## Lactylation-regulator annotations

A curated list of Kla Writers, Erasers and Readers was mapped to human UniProt
BaseAccessions (Supplementary Table S3). Proteins with more than one reported
regulatory role were retained in each corresponding display facet. Protein-set
statistics remained deduplicated by BaseAccession even when the same regulator
appeared under multiple functional labels.

## Within-sample Kla percentile heatmap

Raw Kla intensities were not compared directly across studies because the
datasets differed in enrichment strategy, instrumentation, acquisition mode,
search software and quantitative level. Within each quantitative sample,
signals belonging to the same normalized Kla feature were summed and
transformed as $\log_2(x+1)$. Features were ranked within that sample using
average ranks for ties, and the percentile was calculated as

$$
\mathrm{percentile} =
100 \times
\frac{\mathrm{average\ rank}-1}{n_{\mathrm{features}}-1}.
$$

A sample containing a single quantitative feature was assigned a percentile
of 100. Regulator percentiles were then extracted, and the median percentile
was used across conditions or replicates belonging to the same sample group.
A regulator not detected in an otherwise usable quantitative sample was
assigned 0 for visualization. This value denotes an undetected feature under
the plotting convention and not a measured biological zero. The resulting
heatmap represents relative positions within individual lactylomes rather than
cross-study abundance or Kla fold change.

## Whole-proteome percentile heatmap

Whole-proteome percentiles were calculated from protein-level quantities in
the selected non-PTM reference matrices. Decoys, contaminants and invalid
protein records were excluded, and duplicate BaseAccessions were consolidated.
Raw intensities were transformed as $\log_2(x+1)$; quantities already reported
on a log2 scale were ranked without an additional transformation. Percentiles
were calculated among all valid proteins within each quantitative sample
before regulator values were extracted. The median percentile was used across
replicates or linked control samples. An undetected regulator in an otherwise
usable whole-proteome sample was assigned 0 under the same visualization
convention.

Both heatmaps were ordered as non-tumor tissues, tumor tissues, cancer cell
lines and normal cell lines, with Writer, Eraser, Writer-Eraser and Reader
column facets. Kla rows described the material and treatment state of the Kla
dataset. Whole-proteome rows independently described the material and state of
the reference dataset and included the corresponding PXD accession. The Kla
heatmap used a white-to-blue scale, and the whole-proteome heatmap used a
white-to-orange-red scale. Chinese and English figures used identical data and
axis orders.

## Exact four-set Venn analyses

Four-category protein sets were constructed separately for all Kla proteins,
Kla-DDR proteins, whole-proteome proteins and whole-proteome DDR proteins.
Within each analysis, proteins were deduplicated by BaseAccession and assigned
to one of the 15 mutually exclusive logical regions defined by membership in
the four biological categories. Exact set memberships and region counts are
provided in Supplementary Table S4.

The four-set diagrams were drawn using fixed schematic geometry. Region area
was not used to encode protein abundance or intersection size; the printed
region counts were the only quantitative encoding. Zero-count regions were
retained, preventing small or empty intersections from being hidden by an
area-fitting procedure. The displayed counts were reconstructed directly from
the underlying membership table.

## Direct-GO-term assignment to DNA-repair pathways

The union of Kla-DDR proteins in the final 30-group analysis contained 399
unique BaseAccessions. Direct GO annotations for these proteins were obtained
from UniProt release 2026_02 (taxon 9606). Biological-process, cellular-
component and molecular-function annotations were all retained, producing
10,605 unique protein-GO pairs and 2,785 unique direct GO terms
(Supplementary Data S3).

Seven DNA-repair pathways were considered: base-excision repair (BER),
nucleotide-excision repair (NER), mismatch repair (MMR), the Fanconi-anemia
pathway (FA), homologous recombination (HR), classical non-homologous end
joining (NHEJ), and alternative end joining (AEJ). For each pathway, a
manually reviewed set of pathway-specific GO seeds was defined
(Supplementary Table S5). A rule either matched the seed term exactly or,
when biologically justified, included all descendants of the seed according
to GO.db version 3.20.0. Descendant expansion was not applied to broad or
mechanistically ambiguous seeds.

A direct GO term could be assigned to more than one pathway when its
definition explicitly spanned multiple repair mechanisms. For example,
interstrand cross-link repair terms could map to FA together with HR or NER,
and generic non-homologous end-joining terms could map to both NHEJ and AEJ.
Terms that did not match any of the seven curated rule sets were assigned to
`Others`. Broad terms such as generic DNA repair, DNA-damage response, DNA
binding or DNA-ligase activity were not assigned to a specific pathway without
pathway-specific evidence. Positive and negative regulatory terms both
indicated involvement in a pathway and were not interpreted as pathway
activation or inhibition.

Recent UniProt terms absent from GO.db were retained in the annotation audit.
Such terms could enter a pathway only through an exact manually reviewed rule;
otherwise, they remained `Others`. Every one of the 2,785 direct terms
received either at least one seven-pathway assignment or the mutually
exclusive `Others` designation. In total, 103 direct terms mapped to at least
one of the seven pathways, and eight terms mapped to two pathways.

## Linear pathway matrices and pathway summaries

The 399-protein union and the four category-specific subsets were plotted
separately. The subset sizes were 183 proteins for non-tumor tissues, 178 for
tumor tissues, 381 for cancer cell lines, 292 for normal cell lines and 399 for
the complete union.

For each protein and pathway, the pathway annotation count was the number of
distinct direct GO terms assigned to that pathway. Counts were summed across
the seven pathways solely to order proteins in ascending order within each
panel; BaseAccession was used as the deterministic tie-breaker. The displayed
matrix was binary: a pathway was shown in its fixed pathway color when at
least one corresponding direct GO term was present and in light gray
otherwise. `Others` was retained for completeness auditing but was not plotted
as an eighth pathway and did not contribute to the ordering score. Separate
summary panels reported the number and proportion of proteins with at least
one direct term assigned to each pathway. The former signed manual score and
pathway coefficients were not used in this primary analysis.

## Software, quality control and reproducibility

Data processing, statistical summaries and visualization were performed in R
4.4.3 using data.table 1.18.2.1, dplyr 1.2.0, readr 2.2.0, tidyr 1.3.2,
ggplot2 4.0.2, ggVennDiagram 1.5.7, GO.db 3.20.0, readxl 1.4.5, ragg 1.5.1
and digest 0.6.39. Dataset-specific evidence extraction also used Python
3.12.10. Figure labels used Arial Unicode MS to maintain consistent rendering
between the English and Chinese versions.

Automated checks verified the 40-source/37-quantifiable/30-paired-group scope,
28 unique whole-proteome display rows, category counts of 9/2/12/7, 399 unique
Kla-DDR proteins, five subset sizes of 183/178/381/292/399, and exact
reconstruction of all Venn regions. The pathway checks verified 10,605 direct
protein-GO pairs, 2,785 unique terms, exhaustive seven-pathway/`Others`
decisions, eight multi-pathway terms and 35 five-set pathway-summary rows.
Final deliverables and software versions were recorded using SHA256 checksums.
UMAP, t-SNE, PCA and Cytoscape analyses were not part of the present analysis.

## Data and code availability

ProteomeXchange accessions, source publications, inclusion decisions and
whole-proteome reference relationships are listed in Supplementary Tables S1
and S2. Harmonized Kla evidence, GO-DDR annotations, exact set memberships and
direct-GO-term pathway decisions are provided in the corresponding
Supplementary Data and Tables described above. Raw and deposited processed
files remain available from their original public repositories. Analysis code
and the reproducible workflow will be made available at
`[repository URL to be added]`.

---

## Alternative pathway-matrix method based on the revised curated score table

<!--
This section is an alternative manuscript version for the revised signed-score
workbook. Retain either this section or the direct-GO-term pathway sections as
the primary method after the final manuscript decision.
-->

As an alternative pathway annotation strategy, we used a manually curated
signed pathway table covering 507 unique UniProt BaseAccessions
(Supplementary Table S6). The table assigned each protein a state of `+1`,
`0` or `-1` for BER, NER, MMR, FA, HR, NHEJ and AEJ. A value of `+1`
represented a curated promoting association, `-1` represented a suppressing
association, and `0` indicated that no directional assignment was made.
Chromatin-interaction and other transcriptional, RNA-processing or
proteostasis annotations were retained in the source table for audit but were
not included in the seven-row pathway matrix.

The curated table was intersected with the current 30-group Kla-DDR membership
using BaseAccession. All 399 proteins in the current union were covered; the
remaining 108 scored proteins were excluded from this analysis because they
were not members of the current Kla-DDR union. The same four category-specific
sets and complete union were used as in the direct-GO-term analysis, yielding
set sizes of 183, 178, 381, 292 and 399 proteins.

For visualization, a weighted signed score was calculated for each protein as

$$
\mathrm{weighted\ score} =
1s_{\mathrm{BER}} + 2s_{\mathrm{NER}} + 3s_{\mathrm{MMR}} +
4s_{\mathrm{FA}} + 5s_{\mathrm{HR}} + 6s_{\mathrm{AEJ}} +
7s_{\mathrm{NHEJ}},
$$

where \(s_{\mathrm{pathway}}\) denotes the corresponding signed state
(`-1`, `0` or `+1`). The weighted score was used only to order proteins in
ascending order within each panel and did not alter pathway states or colors.
BaseAccession was used as the deterministic tie-breaker. In the linear matrix,
`+1` states were displayed using the fixed pathway color, `-1` states were
displayed in dark charcoal, and `0` states were displayed in light gray.
Separate summary panels reported promoting, suppressing and unassigned counts
and proportions for each pathway. English and Chinese versions used identical
protein orders, pathway orders and color assignments.
