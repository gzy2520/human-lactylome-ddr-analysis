# 37组普通全蛋白参照更新报告

更新日期：2026-08-07

## 本次替换

原 HUVEC 参照 `PXD009687` 的 `prot.xml` 只有蛋白鉴定，没有逐蛋白定量强度，
不适合进入跨样本强度热图。因此重新采用同研究的普通全蛋白定量矩阵：

`data/PXD073311/search_results/extracted_pairing/IPX0015307001_Database_search_result/Database_search_result/report.pg_matrix.tsv`

只使用 `A0h_1`、`A0h_2`、`A0h_3` 三个基线/对照列，排除 `A6h`。三列共有
7,709 个阳性蛋白组行，拆分并去除 isoform 后为 7,794 个唯一 UniProt
BaseAccession。

## 最终范围

普通全蛋白热图、DDR占比图和四分类 Venn 都严格跟随：

`reanalysis/results/tables/kla_regulator_intensity_availability_audit.csv`

中 `定量可用=TRUE` 的 37 个 PXD+样本组。PXD037371 的
`normal liver`、`nonmetastatic HCC` 和 `lung-metastatic HCC` 不在最终范围，
但原始数据和旧结果均保留在 archive。

最终四分类样本组数为：

| 分类 | 样本组数 |
|---|---:|
| 正常/非肿瘤组织 | 10 |
| 癌症组织 | 3 |
| 正常/非肿瘤细胞 | 10 |
| 癌症细胞 | 14 |

## 结果变化

- Kla DDR统计的 37 行蛋白数和 DDR 数与归档基线完全一致。
- Kla 四分类 Venn：全部 Kla 为 8,719 个 BaseAccession，Kla-DDR 为 510 个；
  新旧 membership 和区域逐 ID 一致。
- 普通全蛋白四分类 Venn：普通全蛋白新增 64、移除 50 个 BaseAccession；
  普通全蛋白 DDR 新增 6 个 BaseAccession。
- HUVEC 普通全蛋白 DDR 为 512/7,794，比例为 6.569%。
- 全部37组普通全蛋白调控因子热图都有定量来源，未使用 Kla 富集强度替代，
  未出现整行 `?`。
- 普通全蛋白热图原先10行整行纯白是 Ensembl protein ID 未转换为 UniProt
  BaseAccession 的代码问题，现已修正；10个 PXD010154 健康组织参照组均恢复调控因子命中。
- Kla 热图的 5 个纯白样本组已单独核查：PXD028488 的 3 组和 PXD064912
  在 Kla 来源中确实没有目标调控蛋白 accession 命中；PXD014870/MCF7 有 1 个条件检出
  AARS2，但 4 条件百分位取中位数后为 0。详见
  `reanalysis/results/tables/kla_regulator_intensity_pure_white_audit.csv`。

## 方法约束

所有分析匹配都使用去 isoform 后缀的 UniProt BaseAccession。GeneSymbol 和
Protein name 只用于显示及人工审计，GeneSymbol 回退数为0。普通全蛋白热图在
每个样本的全部普通蛋白特征中计算 `log2(signal+1)` 后的平均秩百分位：

`100 × (average_rank - 1) / (n - 1)`

重复或对照列在样本组内取百分位中位数。白色表示0或未检出，颜色从白色依次过渡到
黄色、橙色和深红色。

## 主要输出

- `reanalysis/results/figures/kla_regulator_whole_proteome_relative_intensity_heatmap.png`
- `reanalysis/results/figures/cell_type_kla_vs_reference_ddr_fraction.png`
- `reanalysis/results/figures/cell_type_kla_vs_reference_ddr_fraction_accession_only.png`
- `reanalysis/results/figures/four_class_venn/reference_proteome_four_class_venn_zh.png`
- `reanalysis/results/figures/four_class_venn/reference_proteome_ddr_four_class_venn_zh.png`
- `reanalysis/results/tables/cell_type_kla_vs_reference_ddr_statistics_accession_only.csv`
- `reanalysis/results/tables/cell_type_kla_vs_reference_ddr_statistics_accession_only_zh.csv`
- `reanalysis/results/tables/four_class_venn/`
- `reanalysis/results/tables/cell_type_kla_vs_reference_ddr_statistics.xlsx`
- `reanalysis/results/tables/lactylome_and_reference_proteome_pairing_zh.xlsx`
- `reanalysis/results/tables/four_class_venn_tables.xlsx`
- 完整绝对路径清单：`reanalysis/reports/FINAL_OUTPUT_FILE_LIST_2026-08-07.txt`
