# Candidate manuscript text: a two-step RNA–proteome association

This is a proposed addition for coauthor review, not an edit to the frozen manuscript DOCX. The model explains **detection across the assembled experiments**; it does not predict the percentage of lactylated proteins in the total proteome. The earlier RNA/proteome bridge text remains a separate, descriptive tissue-contrast candidate.

## Results — candidate subsection

### RNA expression and sequential protein detection across lactate-related genes

We asked whether the tissue-matched RNA references and the ordinary proteome could help explain which lactate-related proteins were identified in the Kla atlas. The previously curated panel contained 51 measured genes with one-to-one Ensembl-to-UniProt mappings. Across 31 proteomics material groups, this yielded 1,581 gene–material observations; 928 had an ordinary-proteome detection and 175 had a Kla-protein detection. RNA references were matched by material category rather than by patient.

In a first-stage logistic model, RNA expression and the total number of ordinary-proteome detections distinguished whether the corresponding ordinary protein was observed (source-held-out ROC AUC 0.854, average precision 0.875), compared with AUC 0.613 and average precision 0.664 for detection depth alone. In a second-stage model, RNA expression, ordinary-protein detection, and detection depth distinguished Kla-protein detections (source-held-out AUC 0.733, average precision 0.249). The corresponding model without RNA had AUC 0.574 and average precision 0.126; the Kla-positive fraction was 11.1%. Source-held-out predictions were assembled by leaving out each of 15 connected study-source blocks in turn. An identity-adjusted sensitivity analysis showed a much smaller incremental contribution from RNA: adding RNA to gene identity, ordinary-protein detection, and depth changed the Kla AUC from 0.760 to 0.772. Thus RNA expression provides a visible association with protein detection in these datasets, while stable gene-specific detectability explains some of the same signal.

This analysis has a different outcome from the originally proposed material-level percentage. For that percentage, the ordinary-proteome detection profile alone gave a source-mean absolute error of 7.99 percentage points, whereas adding RNA raised it to 11.83 points. We therefore do not present the two-step detection AUC as evidence that RNA accurately predicts the fraction of Kla proteins among all detected proteins.

## Methods — candidate addition

The fixed 51-gene panel was taken from the curated lactate-metabolism evidence table. Joins used Ensembl gene IDs and exact UniProt accessions; symbols were used only as display labels. Each proteomics material group was paired with its existing RNA `ReferenceKey`. The expression input was unsmoothed log2(TPM + 0.5), summarized in the earlier RNA pipeline; these RNA samples were not collected from the proteomics patients. Ordinary-protein and Kla-protein outcomes were binary accession-level detections in the group's published membership lists. An undetected protein was coded zero for the respective assay; this is a detection label, not biological absence. Ordinary-proteome depth was log1p of the number of unique accessions detected in the group.

Binomial logistic regression was fitted for (i) ordinary-protein detection using RNA and depth and (ii) Kla-protein detection using RNA, ordinary-protein detection, and depth. Baselines omitted RNA. Sensitivity models also included stable gene identity as a categorical covariate. Source-blocked validation left out one of 15 connected groups at a time; groups sharing a Kla study, ordinary-proteome study, or RNA reference belonged to the same block. Predictions from the held-out blocks were pooled for ROC AUC, average precision, and Brier score. The relationship chart also displays full-data fitted values and explicitly labels them as descriptive. The random seed was 25. The 1,581 gene–material rows were not treated as 1,581 independent patients, and the 31 groups were not represented as independent RNA references; there were 28 distinct RNA references.

The separate percentage benchmark compared a source-balanced median, the previously fixed 16 RNA modules, ordinary-proteome detection depth, a 24,397-accession binary ordinary-proteome profile, and a joint RNA-plus-profile model. It used a nested source-blocked penalty choice and outer leave-one-source-block-out evaluation. The percentage is the fraction of ordinary-proteome detected proteins that were also detected in the Kla-protein list; it is not site occupancy or a quantitative Kla-to-total-protein ratio. The ordinary-proteome profile contains information about the denominator and detection coverage, so its benchmark performance should not be interpreted as a biological mechanism.

## Discussion — candidate paragraph

The sequential detection model offers a limited but direct bridge between the RNA and proteomics layers. More highly expressed genes were more likely to yield an ordinary-protein identification, and RNA added discrimination for which gene–material pairs yielded a Kla-protein identification within the assembled datasets. However, the RNA and mass-spectrometry specimens were not patient matched, protein abundance was not harmonized quantitatively across studies, and the Kla outcome records detection of at least one modified peptide rather than modification intensity or occupancy. The association could reflect protein availability, assay coverage, and gene-specific detectability; it cannot establish that lactate metabolism caused a change in lactylation. Paired quantitative RNA, protein, Kla-site and lactate measurements would be required to test that mechanism or construct a patient-level predictor.

## Candidate figure legends

**Figure Sx. RNA expression and sequential protein detection.** Blue points show observed detection fractions across RNA deciles; red points show full-data fitted fractions from the two logistic models. The left panel concerns ordinary-protein detection across all 1,581 gene–material records. The right panel concerns Kla detection among the 928 records with ordinary-protein detection. Deciles were formed separately for each panel. This chart describes the assembled data and is not a cross-validated calibration plot. RNA was linked by material category, not patient.

**Figure Sy. Source-held-out discrimination for two detection tasks.** ROC curves use predictions obtained after withholding each of 15 connected source blocks. The left model uses RNA and ordinary-proteome depth to classify ordinary-protein detection (AUC 0.854); the right uses RNA, ordinary-protein detection, and depth to classify Kla-protein detection (AUC 0.733). The training-fit AUC values shown in the panel headers are labeled separately. A Kla-positive result denotes mass-spectrometry detection, not protein-normalized Kla occupancy.

## Suggested message to the supervisor (Chinese)

我们把现有数据做成了一个“两步检出关联”模型：先用 RNA 及蛋白组检测深度解释普通蛋白是否检出，再结合普通蛋白检出状态解释 Kla 蛋白是否检出。按研究来源整组留出，两个任务的 AUC 分别为 0.854 和 0.733；有两张可以展示的关联图和验证图。但这不是“输入 RNA 预测乳酸化蛋白占全蛋白百分比”的成功模型，那个原目标加入 RNA 后仍未改善。RNA 与质谱不是同一患者，所以文章里应称为跨组学关联，不能称为患者级预测或乳酸化机制证明。

## Reproducible evidence

- Scripts: `workflow/rna_proteome_ml_benchmark/explain_association.R` and `workflow/rna_proteome_ml_benchmark/run.R`.
- New results: `outputs/20260924_rna_protein_explanatory/` and `outputs/20260924_rna_proteome_ml_benchmark/`.
- The result directories include source-held-out predictions, model specifications, input SHA-256, and release SHA-256 manifests.
