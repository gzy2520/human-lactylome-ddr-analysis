# Kla 项目交接 Prompt

请继续处理以下项目：

`/Users/gzy2520/Desktop/Research/kla`

这是一个人源蛋白乳酸化质谱数据整合、严格普通全蛋白参照配对、GO-DDR
交集和乳酸化调控因子热图项目。开始工作前先读取：

1. `PROJECT_INDEX.md`
2. `reanalysis/README.md`
3. `reanalysis/config/lactylome_reference_pairing.csv`
4. `reanalysis/results/tables/kla_and_reference_teacher_review_zh.csv`
5. `reanalysis/reports/KLA_REGULATOR_INTENSITY_HEATMAP_METHODS.md`
6. `reanalysis/reports/EXPANDED_DDR_FRACTION_ACCESSION_ONLY.md`
7. `reanalysis/reports/METHODS_CURRENT_33GROUP_ZH.md`

## 不可改变的分析规则

- 生物信息学分析优先使用 R。
- 蛋白分析、去重、合并、DDR 匹配和 Venn 集合只能使用去除 isoform
  后缀的 UniProt `BaseAccession`。
- GeneSymbol 和 Protein name 只能用于显示和人工审计，禁止作为命中回退。
- Ensembl protein ID 必须先通过显式映射转换为 UniProt。
- GO 输入为 `data/annotations/GO-repair+damage(human).tsv`；排除 `NOT`
  注释，主分析暂不按 evidence code 缩减。
- Kla 证据必须来自作者检索结果中的 Kla 位点/修饰肽段/蛋白表或作者补充表，
  不能直接把原始谱图当作 Kla 蛋白证据。
- 普通全蛋白强度只能来自未做 Kla/PTM 富集的全蛋白矩阵，不能用 Kla 信号替代。
- 不删除原始数据；废弃结果移动到 `archive/` 并记录迁移清单。
- 图必须同时输出中文和英文版本。热图色阶为白色到黄色、橙色、深红色。

## 当前最终范围

来源审计保留 37 个具有 Kla 定量的样本组。最终配对热图、DDR 比较和四分类
Venn 只使用其中 33 组，分类数量固定为：

- 正常/非肿瘤组织：9
- 癌症组织：2
- 正常/非肿瘤细胞：9
- 癌症细胞：13

四个无完全匹配普通全蛋白强度参照的 Kla 组不进入最终比较：

- `PXD062720 / bladder cancer cells treated with EPI`
- `PXD063047 / severe preeclampsia placenta`
- `PXD064038 / MEC and NEC ESCC groups`
- `PXD075014 / AC16 control and hypoxia`

`PXD063047 / normal pregnancy placenta` 仍保留。`PXD037371` 的 normal
liver、nonmetastatic HCC、lung-metastatic HCC 因 TMT 通道与临床组无法可靠
映射，仅保留来源审计，不进入 Kla 定量范围。

33 个 Kla 组对应 30 个唯一普通全蛋白参照行。HK-2、MCF7、HCT116 各有两项
Kla 研究共享完全相同的普通全蛋白参照，统计表保留 33 条配对，普通全蛋白热图
只显示 30 条唯一参照。Kla 热图保留 33 行，并严格按这 30 条普通全蛋白参照的
顺序展开共享参照。

## 当前有效输入和配置

- Kla 与普通全蛋白配对：
  `reanalysis/config/lactylome_reference_pairing.csv`
- 四分类和显示顺序：
  `reanalysis/config/four_class_sample_grouping.csv`
- 严格参照排除：
  `reanalysis/config/strict_reference_exclusions.csv`
- 材料身份审计：
  `reanalysis/config/strict_reference_material_identity_review.csv`
- 调控因子：
  `data/identifier/乳酸化调控因子_Writer-Eraser-Reader.xlsx`
- 调控因子稳定 ID 映射：
  `reanalysis/config/lactylation_regulator_uniprot_mapping.csv`
- GO-DDR：
  `data/annotations/GO-repair+damage(human).tsv`
- Ensembl 到 UniProt 映射：
  `reanalysis/config/ensembl_protein_to_uniprot_biomart.tsv`

每个 PXD 的实际原始文件、检索结果、补充表和元数据分别放在
`data/PXD*/raw`、`search_results`、`supplementary`、`metadata`。不要移动或
重复复制大型原始数据，除非用户明确要求。

## 当前有效代码

数据盘点、下载和解压：

- `reanalysis/scripts/build_human_lactylome_inventory.R`
- `reanalysis/scripts/build_lactylome_acquisition_manifests.R`
- `reanalysis/scripts/download_lactylome_pair_files.R`
- `reanalysis/scripts/extract_lactylome_pair_archives.R`
- `reanalysis/scripts/summarize_acquired_lactylome_data.R`

Kla 主提取和 GO-DDR：

- `reanalysis/scripts/run_pipeline.py`
- `reanalysis/scripts/common.py`
- `reanalysis/scripts/extractors.py`
- `reanalysis/scripts/analyze_ddr.R`

严格普通全蛋白参照和老师审阅表：

- `reanalysis/scripts/apply_exact_reference_policy.R`
- `reanalysis/scripts/build_lactylome_reference_pairing.R`
- `reanalysis/scripts/build_strict_reference_material_identity_audit.R`
- `reanalysis/scripts/build_teacher_review_table.R`

最终比较和作图：

- `reanalysis/scripts/analyze_kla_regulator_whole_proteome_intensity.R`
- `reanalysis/scripts/analyze_kla_regulator_intensity.R`
- `reanalysis/scripts/analyze_expanded_ddr_fraction_by_accession.R`
- `reanalysis/scripts/plot_four_class_area_proportional_venn.R`
- `reanalysis/scripts/build_final_manifest.R`

## 必要中间表

- Kla sample-level 长表：
  `reanalysis/results/tables/kla_regulator_intensity_sample_level_long.csv`
- 普通全蛋白调控因子 sample-level 长表：
  `reanalysis/results/tables/kla_regulator_whole_proteome_regulator_sample_level_long.csv`
- 33 行 Kla 轴：
  `reanalysis/results/tables/kla_regulator_heatmap_axis_order.csv`
- Kla/普通全蛋白轴对齐：
  `reanalysis/results/tables/kla_vs_whole_proteome_heatmap_axis_alignment.csv`
- 30 行普通全蛋白轴：
  `reanalysis/results/tables/kla_regulator_whole_proteome_heatmap_rows.csv`
- Kla 蛋白集合：
  `reanalysis/intermediate/expanded_ddr_by_accession/kla_proteins_by_sample_group.csv`
- 普通全蛋白集合：
  `reanalysis/intermediate/expanded_ddr_by_accession/reference_proteins_by_sample_group.csv`
- Venn 范围审计：
  `reanalysis/results/tables/four_class_venn/venn_sample_group_scope.csv`

## 推荐重跑顺序

```bash
cd /Users/gzy2520/Desktop/Research/kla

Rscript reanalysis/scripts/analyze_kla_regulator_whole_proteome_intensity.R .
Rscript reanalysis/scripts/analyze_kla_regulator_intensity.R .
Rscript reanalysis/scripts/analyze_expanded_ddr_fraction_by_accession.R .
Rscript reanalysis/scripts/plot_four_class_area_proportional_venn.R .

for test_file in $(rg --files reanalysis/tests | rg 'test_.*\.R$' | sort); do
  Rscript "$test_file" . || exit 1
done

PYTHONPATH=reanalysis/scripts python3 \
  -m unittest discover -s reanalysis/tests -p 'test_*.py' -v

Rscript reanalysis/scripts/build_final_manifest.R .
```

普通全蛋白热图必须先运行，因为 Kla 热图从其 30 行参照轴读取显示顺序。

## 当前正式成品

热图：

- `reanalysis/results/figures/kla_regulator_cross_study_relative_intensity_heatmap_zh.png`
- `reanalysis/results/figures/kla_regulator_cross_study_relative_intensity_heatmap_en.png`
- `reanalysis/results/figures/kla_regulator_whole_proteome_relative_intensity_heatmap_zh.png`
- `reanalysis/results/figures/kla_regulator_whole_proteome_relative_intensity_heatmap_en.png`

DDR 柱状图：

- `reanalysis/results/figures/cell_type_kla_vs_reference_ddr_fraction_zh.png`
- `reanalysis/results/figures/cell_type_kla_vs_reference_ddr_fraction_en.png`

四分类 Venn：

- `reanalysis/results/figures/four_class_venn/*_33groups_zh.png`
- `reanalysis/results/figures/four_class_venn/*_33groups_en.png`

33 组统计表：

- `reanalysis/results/tables/cell_type_kla_vs_reference_ddr_statistics_accession_only_paired_33.csv`
- `reanalysis/results/tables/cell_type_kla_vs_reference_ddr_statistics_accession_only_paired_33_zh.csv`

每套 Venn 的 `membership.csv`、`region_counts.csv`、`set_counts.csv` 和稳定
ID 映射审计位于 `reanalysis/results/tables/four_class_venn/<analysis>/`。

最终交付前必须确认：

1. Kla 轴 33 行、普通全蛋白轴 30 行。
2. 33 个 Kla `PXD+SampleGroup` 与严格配对表完全一致。
3. 四类数量为 9/2/9/13。
4. 四个无参照组不进入热图、DDR 比较或 Venn。
5. `GeneSymbolFallbackCount=0`。
6. 全部 R 和 Python 测试通过。
7. 重新生成 `reanalysis/reports/final_file_manifest_sha256.csv`。

可直接用于论文/汇报的方法稿：

- 中文：`reanalysis/reports/METHODS_CURRENT_33GROUP_ZH.md`
- 英文：`reanalysis/reports/METHODS_CURRENT_33GROUP_EN.md`
- Genome Biology 范文结构说明：
  `reanalysis/reports/METHODS_STYLE_REFERENCE.md`
