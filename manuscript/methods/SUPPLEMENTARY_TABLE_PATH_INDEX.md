# Supplementary table and analysis-table path index

本文件用于标识 `METHODS_EN.md` 中暂定的 Supplementary Table/Data 编号与
项目真实文件路径。编号仍为投稿前占位编号；最终投稿时可将同一编号下的多个
CSV/TSV合并为一个含多个工作表的Excel文件，并统一改为英文列名。

## 已生成的当前30组工作簿

- `results/supplementary/Supplementary_Tables_S1-S6_Draft_20260817.xlsx`
- `results/supplementary/Supplementary_Data_S1_Kla_Evidence_Draft_20260817.xlsx`
- `results/supplementary/Supplementary_Data_S2_Human_GO_DDR_Draft_20260817.xlsx`
- `results/supplementary/Supplementary_Data_S3_GO_Pathway_Draft_20260817.xlsx`
- `results/supplementary/supplementary_workbook_manifest.csv`

## Methods中暂定的补充表

### Supplementary Table S1 — dataset catalogue and quantitative summary

主来源：

- `results/tables/cell_type_kla_vs_reference_ddr_statistics_accession_only_paired_30.csv`：
  当前30组的类别、材料、Kla与DDR蛋白数、普通全蛋白参照及对应比例。
- `config/sample_group_catalog.csv`：仅用于补充当前30组的论文年份、DOI、
  数据集标题和ProteomeXchange链接。

生成表只保留当前30组，并使用英文列名。

### Supplementary Table S2 — whole-proteome reference pairing and audit

主来源：

- `config/lactylome_reference_pairing.csv`：Kla样本组与普通全蛋白参照的配置
  关系。
- `config/strict_reference_material_identity_review.csv`：材料身份和实验状态
  匹配的人工复核。
- `results/tables/strict_reference_material_identity_audit.csv`：程序生成的完整
  参照配对审计结果。

生成表以程序审计结果为主表，仅保留当前30组；保留
`MaterialIdentityMatch`、`ExperimentalStateMatch`、参照PXD、实际文件和样本
子集，避免把独立队列或基线参照写成严格条件匹配对照。

### Supplementary Table S3 — lactylation-regulator annotations

主来源：

- `data/identifier/乳酸化调控因子_Writer-Eraser-Reader.xlsx`：Writer、
  Eraser和Reader的原始整理表。
- `config/lactylation_regulator_uniprot_mapping.csv`：调控因子到UniProt
  BaseAccession的固定映射。

相关分析审计：

- `results/tables/kla_regulator_dataset_audit.csv`
- `results/tables/kla_regulator_intensity_id_mapping_audit.csv`

### Supplementary Table S4 — exact Venn membership and region counts

四类Venn概览：

- `results/tables/four_class_venn/four_venn_set_counts_4x4.csv`
- `results/tables/four_class_venn/venn_sample_group_scope.csv`

四套精确成员表及区域计数：

- `results/tables/four_class_venn/all_kla_four_class_venn/membership.csv`
- `results/tables/four_class_venn/all_kla_four_class_venn/region_counts.csv`
- `results/tables/four_class_venn/all_kla_four_class_venn/set_counts.csv`
- `results/tables/four_class_venn/kla_ddr_four_class_venn/membership.csv`
- `results/tables/four_class_venn/kla_ddr_four_class_venn/region_counts.csv`
- `results/tables/four_class_venn/kla_ddr_four_class_venn/set_counts.csv`
- `results/tables/four_class_venn/reference_proteome_four_class_venn/membership.csv`
- `results/tables/four_class_venn/reference_proteome_four_class_venn/region_counts.csv`
- `results/tables/four_class_venn/reference_proteome_four_class_venn/set_counts.csv`
- `results/tables/four_class_venn/reference_proteome_ddr_four_class_venn/membership.csv`
- `results/tables/four_class_venn/reference_proteome_ddr_four_class_venn/region_counts.csv`
- `results/tables/four_class_venn/reference_proteome_ddr_four_class_venn/set_counts.csv`

投稿处理建议：一个工作簿设置四个分析分组，每组至少包含`membership`和
`region_counts`两个工作表。图中数字以这些精确表为准，不以图形面积反推。

### Supplementary Table S5 — GO-term-to-pathway decisions

投稿主表：

- `results/tables/go_term_pathway_scoring_30groups/go_term_to_pathway_long.csv`

规则及完整性审计：

- `config/go_term_pathway_seed_rules.csv`
- `results/tables/go_term_pathway_scoring_30groups/go_term_decision_audit_2785.csv`
- `results/tables/go_term_pathway_scoring_30groups/term_pathway_rule_hit_evidence.csv`
- `results/tables/go_term_pathway_scoring_30groups/multi_pathway_go_terms.csv`
- `results/tables/go_term_pathway_scoring_30groups/go_db_missing_terms.csv`

`go_term_to_pathway_long.csv`允许同一GO term出现多行，以保留一个term映射至
多个DNA修复通路的情况；`go_term_decision_audit_2785.csv`用于证明2,785个
直接GO term均已得到七通路或`Others`判定。

### Supplementary Table S6 — revised manually curated signed pathway annotations

当前评分表：

- `data/identifier/乳酸化DDR基因评分表_Revised_20260816.xlsx`

当前399蛋白及4+1作图结果：

- `results/tables/five_set_pathway_matrix_revised_excel_20260816/protein_order_and_seven_pathway_matrix_5sets.csv`
- `results/tables/five_set_pathway_matrix_revised_excel_20260816/pathway_state_summary_5sets_35rows.csv`

S6只输出当前30组Kla-DDR并集中的399个蛋白。

## Methods中暂定的补充数据

### Supplementary Data S1 — sample-level Kla evidence

- `results/tables/cell_type_kla_vs_reference_ddr_statistics_accession_only_paired_30.csv`：
  当前30组概览。
- `work/intermediate/expanded_ddr_by_accession/kla_proteins_by_sample_group.csv`：
  按当前30个`PXD+SampleGroup`过滤后生成Kla蛋白成员表和Kla-DDR成员表。

### Supplementary Data S2 — human GO-DDR annotations

- `data/annotations/GO-repair+damage(human).tsv`

### Supplementary Data S3 — direct GO annotations and pathway matrices

- `results/tables/go_term_pathway_scoring_30groups/direct_go_annotations_399proteins.csv`
- `results/tables/go_term_pathway_scoring_30groups/protein_go_term_pathway_long.csv`
- `results/tables/go_term_pathway_scoring_30groups/protein_pathway_direct_term_counts_long.csv`
- `results/tables/go_term_pathway_scoring_30groups/protein_pathway_direct_term_count_matrix.csv`
- `results/tables/go_term_pathway_scoring_30groups/protein_seven_pathway_binary_matrix.csv`

投稿处理建议：将以上文件作为一个多工作表工作簿或一个压缩数据包，分别保留
直接蛋白–GO配对、蛋白–GO–通路长表、直接term计数矩阵和七通路二值矩阵。

## 当前四类正式图对应的结果表

### DDR fraction柱状图

- `results/tables/cell_type_kla_vs_reference_ddr_statistics_accession_only_paired_30.csv`
- `results/tables/cell_type_kla_vs_reference_ddr_plot_rows.csv`

### 两张调控因子强度热图

- `results/tables/kla_regulator_normalized_intensity_long.csv`
- `results/tables/kla_regulator_whole_proteome_normalized_long.csv`
- `results/tables/kla_vs_whole_proteome_heatmap_axis_alignment.csv`

### 四类Venn图

- `results/tables/four_class_venn/four_venn_set_counts_4x4.csv`
- 详细成员和区域计数见Supplementary Table S4对应路径。

### GO-term主版本4+1线性图及summary

- `results/tables/five_set_pathway_matrix_go_term_30groups/protein_order_and_direct_term_counts_5sets.csv`
- `results/tables/five_set_pathway_matrix_go_term_30groups/linear_matrix_plot_data.csv`
- `results/tables/five_set_pathway_matrix_go_term_30groups/pathway_presence_summary_5sets_35rows.csv`
- `results/tables/five_set_pathway_matrix_go_term_30groups/protein_set_counts.csv`

### 老师修订评分表4+1线性图及summary

- `results/tables/five_set_pathway_matrix_revised_excel_20260816/protein_order_and_seven_pathway_matrix_5sets.csv`
- `results/tables/five_set_pathway_matrix_revised_excel_20260816/linear_matrix_plot_data.csv`
- `results/tables/five_set_pathway_matrix_revised_excel_20260816/pathway_state_summary_5sets_35rows.csv`
- `results/tables/five_set_pathway_matrix_revised_excel_20260816/protein_set_counts.csv`

## 可复现性说明

`results/`和`work/intermediate/`中的生成文件不由Git跟踪。缺失时应在项目根
目录运行：

```bash
Rscript workflow/run_pipeline.R selected_figures
Rscript workflow/run_pipeline.R revised_score_preview
Rscript workflow/build_manifest.R /Users/gzy2520/Desktop/Research/kla
```

生成后使用`results/reports/publication_manifest_sha256.csv`核对文件完整性。
