# 本地数据目录

`data/`保存公开质谱项目的原始文件、作者检索结果、补充表和本项目的小型注释输入。
大型PXD文件不进入Git。某个PXD存在于本目录，不代表其进入当前30组主分析。

当前样本范围只能由以下配置确定：

- `../config/sample_group_catalog.csv`
- `../config/four_class_sample_grouping.csv`
- `../config/lactylome_reference_pairing.csv`
- `../config/main_analysis_scope_exclusions.csv`

## 跟踪的小型输入

- `annotations/GO-repair+damage(human).tsv`：人源GO-DDR集合，分析时排除`NOT`。
- `identifier/乳酸化调控因子_Writer-Eraser-Reader.xlsx`：调控蛋白显示列表。

`identifier/260810乳酸化DDR基因评分表.xlsx`属于旧人工评分版本，保留仅供历史
审计，不参与当前GO-term七通路评分或图形。

`identifier/乳酸化DDR基因评分表_Revised_20260816.xlsx`是老师修订的507蛋白
人工评分表。它只用于独立的399蛋白4+1试绘，不替代当前直接GO-term版本。

## 数据使用原则

- 原始谱图用于来源追溯，不直接等同于Kla证据。
- 主分析读取作者检索后的Kla位点、修饰肽段、蛋白定量或作者补充表。
- 普通全蛋白参照必须来自非Kla/PTM富集矩阵。
- 所有蛋白集合分析使用去isoform的UniProt`BaseAccession`。
- 实际使用的文件、字段、样本列和配对关系保存在结果审计表及
  `../docs/DATA_PROVENANCE.md`。
