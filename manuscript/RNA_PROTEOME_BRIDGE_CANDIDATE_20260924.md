# Candidate manuscript addition: metabolic RNA context for the Kla proteome

This is a proposed addition to the manuscript, not an edit to the frozen DOCX. The safest placement is a short Results subsection after the writer/eraser/reader landscape, with the two new graphics as exploratory supplementary figures. Symbols below are display names; all computational joins used Ensembl gene IDs and exact UniProt accessions.

## Results — proposed new subsection

### Tissue-level transcriptional context of lactate metabolism

To place the lactylated-protein atlas in a metabolic context, we examined a prespecified 21-gene panel spanning lactate production, transport, pyruvate utilization, gluconeogenesis, and candidate lactyl-donor formation. RNA-seq references were matched to the proteomics groups by tissue category, not by donor. In TCGA-LIHC, the HCC-to-adjacent-liver difference was small for LDHA (−0.07 on the log2[TPM + 0.5] scale), while LDHB was modestly higher (+0.26). PCK1 (−2.31), FBP1 (−2.00), and MPC1 (−1.07) were lower in the HCC RNA reference. The corresponding Kla datasets identified 1,775 and 1,798 lactylated-protein species in adjacent liver and HCC, respectively, including 113 and 121 DDR proteins (6.37% and 6.73% of the detected Kla-protein sets). These are set-composition measures; they do not quantify total cellular lactylation or modification occupancy.

In the prostate comparison, the external RNA references showed higher LDHA (+1.23) and lower LDHB (−4.69) in TCGA-PRAD than in the GSE132714 BPH cohort. We then inspected the separately measured ordinary proteome from PXD066054, which contains five BPH and five prostate-cancer specimens. Among the 21 genes, 18 had an unambiguous exact UniProt protein group. The LDHA protein-quantity median was similar between groups (log2 cancer/BPH ratio +0.05), whereas LDHB was lower (−1.17). Although the group-level Kla-protein sets overlapped extensively (2,116 proteins shared among 2,118 BPH and 2,125 cancer detections), this binary overlap does not establish equality of site-specific modified-peptide signals. As an exploratory illustration, the exported LDHA site quantities had mixed directions across five sites detected in at least three samples per group, whereas LDHB sites tended to be lower in cancer. Thus, the tissue-level RNA contrasts and the lactylome describe related but distinct biological and measurement layers; no direct RNA-to-Kla relationship can be estimated from these unpaired cohorts.

## Discussion — proposed paragraph after the tumor-versus-non-tumor Kla comparison

The metabolic RNA results offer a context for interpreting the broad overlap of lactylated-protein identities without invoking a demonstrated homeostatic balance between “writers” and “erasers.” Transcripts describe potential enzyme supply; ordinary proteomics measures protein abundance; the detected Kla-protein list records whether any eligible modified peptide was observed in a group; and modified-site quantities may vary even when protein identity remains shared. In addition, LDH catalysis is reversible, lactate transport depends on gradients, and lactyl-donor availability and enzyme activity are not determined by RNA abundance alone. Mechanistic studies have established routes involving lactate-dependent AARS1/2 activity (Li et al., 2024), ACSS2-linked lactyl-CoA production (Zhu et al., 2025), and delactylation by HDAC1–3 (Moreno-Yruela et al., 2022) in specific experimental systems; our data do not show which route dominates in HCC or prostate cancer. A hypothesis of increased synthesis counterbalanced by increased removal, or of lactate diversion into the TCA cycle, would require matched lactate and metabolite measurements, quantitative site-resolved Kla normalized to its parent protein, and ideally isotope-tracing or perturbation experiments. We therefore present the RNA data as pathway context rather than a predictor or mechanistic explanation of Kla abundance.

## Methods — proposed addition

The fixed 21-gene panel and the TCGA-LIHC primary/adjacent and TCGA-PRAD/GSE132714-BPH RNA reference contrasts were taken from the independently validated lactate-metabolism analysis. The group-level RNA difference is the difference of qsmooth-A expression values on the log2(TPM + 0.5) scale; it is descriptive and is not a differential-expression test or an exact log2 fold change. The prostate ordinary-proteome quantities were read from the PXD066054 Spectronaut `Protein_Quant.tsv` export. Genes were linked by measured Ensembl ID to a curated exact UniProt accession; ambiguous semicolon-delimited protein groups were excluded. For each exact group, the median positive exported quantity within the five BPH and five cancer samples was calculated only when at least three samples per group were detected. The log2 ratio of these medians was used for display, without cross-study protein-quantity comparison. Existing Kla-protein membership was linked by exact accession and checked against the original PXD066054 site-report selection, which reproduced both frozen group-level protein sets exactly. Exploratory LDHA/LDHB site summaries used only L-Lac(K) at lysine, localization probability ≥0.75, positive quantity, and at least three detected samples per group. They are modified-site signal summaries, not estimates of site occupancy. No RNA specimen was matched to a mass-spectrometry specimen. The prostate RNA comparison spans different studies, and the BPH cohort has medication/background differences; it cannot be interpreted as a controlled normal-versus-tumor molecular contrast.

## Proposed supplementary figure legends

**Figure Sx. Lactate-route RNA and ordinary-protein contrasts.** The prespecified 21-gene panel is shown for HCC versus adjacent liver RNA, prostate cancer versus BPH RNA, and prostate cancer versus BPH ordinary protein quantities. RNA points are differences in qsmooth-A log2(TPM + 0.5); protein points are log2 ratios of group medians in PXD066054 (five samples per group). Blank protein rows lacked an unambiguous exact UniProt protein group. Comparisons are descriptive and across independent RNA and protein cohorts; neither matched-patient RNA–protein correlation nor lactate flux is implied.

**Figure Sy. LDHA/LDHB lactylated-site quantities in PXD066054.** Each point represents the log2 cancer/BPH ratio of the median positive modified-site quantity among at least three detected samples per group (localization probability ≥0.75). Opposing LDHA site directions illustrate why group-level presence of a lactylated protein is not a quantitative site-level result. These peptide-derived signals are not protein-normalized modification occupancies, and the panel is exploratory.

## One existing sentence to soften

In Results section 4, the current final sentence says that most regulator components “are themselves regulated by lactylation.” Detection alone establishes that they **were detected as lactylated**, not that their activity was regulated. Suggested replacement:

> In summary, many candidate lactylation-related proteins were present across tissue and cell-line categories, and several were also detected in the lactylated-protein catalog. Whether their lactylation changes enzyme or reader activity remains to be tested.

## Interpretation notes for coauthors

- The prostate `NAT1–NAT5` filenames denote the five BPH samples in this dataset; BPH is not healthy adjacent prostate tissue.
- The 2,116-protein overlap was independently reproduced from the original site report using the frozen selection rule. It is not evidence of unchanged total Kla or site occupancy.
- The 21-gene panel has zero **binary Kla-detection switches** in either tissue comparison. That is a limitation of this measurement, not proof that modification levels are identical.
- The RNA-only prediction model did not show reliable predictive utility. Its negative result can be reported briefly as a feasibility limitation, without making it the central story of the article.
- The existing manuscript's DDR percentage refers to **DDR among detected Kla proteins**. It must not be called the percentage of lactylated protein in the total proteome.

## Primary references to support the mechanism boundaries

- Li et al. 2024, AARS1/2 lactate-sensing and lactyltransferase activity: https://www.nature.com/articles/s41586-024-07992-y
- Zhu et al. 2025, ACSS2/KAT2A lactyl-CoA route in an experimental tumor model: https://pubmed.ncbi.nlm.nih.gov/39561764/
- Moreno-Yruela et al. 2022, HDAC1–3 delactylase activity: https://pubmed.ncbi.nlm.nih.gov/35044827/
- PXD066054 original deposited study and cohort description: https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD066054

## Local, reproducible evidence

- Analysis: `workflow/rna_protein_bridge/analyze.R`
- Figure and table release: `outputs/20260924_rna_protein_bridge/`
- Input and output SHA-256 manifests: `tables/input_manifest.csv`, `release_sha256.csv`

The original manuscript DOCX and the prior accepted RNA/proteomics releases were not changed.
