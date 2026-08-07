# Kla 项目索引

更新时间：2026-08-06  
项目根目录：`/Users/gzy2520/Desktop/Research/kla`

## 项目目标

本项目用于整理人源乳酸化质谱数据，提取可信的 Kla 蛋白/位点，与人源
GO repair/damage（DDR）注释取交集，并为各细胞系或组织寻找普通全蛋白组参照。
当前还包括调控因子在 Kla 数据和普通全蛋白数据中的相对信号热图。

## 最高优先级 ID 规则

- 蛋白分析、去重、GO 匹配和统计必须使用稳定 ID。
- 当前主键为去除 isoform 后缀的 `UniProt BaseAccession`。
- 有 Ensembl protein ID 的参照数据，先通过明确的 Ensembl→UniProt 映射再分析。
- `GeneSymbol` 不得作为分析主键、GO 回退键或蛋白去重键；它只能保留为来源元数据或图中显示标签。
- PXD043880 的原始全蛋白表只有 symbol，是输入格式限制；先转换为人源 reviewed UniProt
  BaseAccession，转换完成后只按 BaseAccession 分析。
- 不得通过 GeneSymbol 回退来扩大交集，也不得用 Kla 信号替代普通全蛋白强度。

## 目录结构

| 目录 | 用途 |
|---|---|
| `data/` | PXD 原始文件、检索结果、补充表、文章和元数据 |
| `previous_umap/` | 最新纯 GO repair/damage UMAP 及其输入、代码、结果 |
| `reanalysis/` | 当前重分析脚本、配置、中间表、结果、报告和测试 |
| `archive/` | 旧脚本、旧结果、迁移清单和不再使用的内容；不得永久删除 |

## Kla 主分析数据

当前主纳入范围：

- `PXD014870`
- `PXD028488`
- `PXD050470`
- `PXD053474`
- `PXD060185`
- `PXD078013`
- `PXD078736`

暂不分析：

- `PXD038880`
- `PXD050906`

这两个 PXD 只能保留原始下载文件、文章/元数据和 `excluded` 说明，不能进入主
Kla 蛋白集合、GO-DDR 交集或 Venn 图。

## 当前分析状态

### Kla 与 GO-DDR

- 扩展比较表包含 40 个细胞系/组织样本组。
- 37 个样本组完成 Kla 与普通全蛋白组的 ID 配对和 DDR 占比计算。
- 3 个 `PXD037371` 样本组因 TMT 通道无法可靠对应临床分组而排除。
- GO 输入唯一使用：
  `data/annotations/GO-repair+damage(human).tsv`
- GO `NOT` 注释已排除；主分析暂不按 evidence code 缩减。
- GO 交集按 UniProt BaseAccession 判定。

主结果：

- [cell_type_kla_vs_reference_ddr_statistics_accession_only_zh.csv](/Users/gzy2520/Desktop/Research/kla/reanalysis/results/tables/cell_type_kla_vs_reference_ddr_statistics_accession_only_zh.csv)
- [cell_type_kla_vs_reference_ddr_accession_only_audit.csv](/Users/gzy2520/Desktop/Research/kla/reanalysis/results/tables/cell_type_kla_vs_reference_ddr_accession_only_audit.csv)
- [cell_type_kla_vs_reference_ddr_fraction.png](/Users/gzy2520/Desktop/Research/kla/reanalysis/results/figures/cell_type_kla_vs_reference_ddr_fraction.png)

### 海马体

Kla 数据：`PXD050470`，作者补充表 S3/S12，853 个 Kla 蛋白。  
普通参照：`PXD043880`，正常 CA1 海马全蛋白 LFQ 矩阵。

- 保留 74 名正常 CA1 供体。
- 原矩阵包含 2,092 个蛋白特征。
- 复合 symbol 特征拆分并转换后得到 2,105 个 reviewed UniProt BaseAccession。
- Kla DDR：29/853，3.40%。
- 普通参照 DDR：83/2105，3.94%。
- 海马 symbol 转换审计：
  [PXD043880_hippocampus_symbol_to_reviewed_uniprot_mapping_audit.csv](/Users/gzy2520/Desktop/Research/kla/reanalysis/results/tables/PXD043880_hippocampus_symbol_to_reviewed_uniprot_mapping_audit.csv)
- reviewed 人源 UniProt 快照：
  [uniprot_human_reviewed_2026-08-05.tsv](/Users/gzy2520/Desktop/Research/kla/reanalysis/config/uniprot_human_reviewed_2026-08-05.tsv)

### 调控因子热图

调控因子来源：

[乳酸化调控因子_Writer-Eraser-Reader.xlsx](/Users/gzy2520/Desktop/Research/kla/data/identifier/乳酸化调控因子_Writer-Eraser-Reader.xlsx)

1. Kla 信号热图
   - 使用 Kla 位点、修饰肽段或作者定量表中的信号。
   - 每个样本内计算相对 Kla 信号百分位。
   - 37/40 个样本组有可审计 Kla 定量值。
   - 不是 Log FC，也不是普通蛋白表达量。

   [kla_regulator_cross_study_relative_intensity_heatmap.png](/Users/gzy2520/Desktop/Research/kla/reanalysis/results/figures/kla_regulator_cross_study_relative_intensity_heatmap.png)

2. 普通全蛋白组热图
   - 使用普通全蛋白组的强度或 LFQ/PG.Quantity/Area 等定量字段。
   - 不使用 Kla 富集信号替代。
   - 每个样本内部计算普通全蛋白信号百分位。
   - 16/40 个样本组有可追溯普通全蛋白强度；没有强度的样本显示 `?`。
   - PXD050470 使用 PXD043880 的 74 个 CA1 供体强度矩阵。

   [kla_regulator_whole_proteome_relative_intensity_heatmap.png](/Users/gzy2520/Desktop/Research/kla/reanalysis/results/figures/kla_regulator_whole_proteome_relative_intensity_heatmap.png)

两张主热图均使用白色到暖色：低值接近白色，高值依次为浅黄、橙色和深红色。

## 细胞/组织分类

当前 Venn 分组配置：

- `hippocampus_tissue`
- `normal_immortalized_cell_lines`
  - `HEK293T`
  - `HK-2`
  - `MCF10A`
- `tumor_cell_lines`
  - 其余按老师当前方案归入肿瘤组的细胞系

MCF10A 按老师当前要求作为第三个永生化模型。所有分组配置见：

- `reanalysis/config/grouping_schemes.csv`
- `reanalysis/config/sample_map.csv`
- `reanalysis/results/tables/venn_regions/`

## 主要脚本

- 主 Kla/GO/Venn 流程：
  `reanalysis/scripts/run_pipeline.py`
- 参考蛋白组 DDR 占比：
  `reanalysis/scripts/analyze_expanded_ddr_fraction_by_accession.R`
- Kla 调控因子信号热图：
  `reanalysis/scripts/analyze_kla_regulator_intensity.R`
- 普通全蛋白组调控因子热图：
  `reanalysis/scripts/analyze_kla_regulator_whole_proteome_intensity.R`
- 最终文件清单：
  `reanalysis/scripts/build_final_manifest.py`

## 推荐运行顺序

```bash
cd /Users/gzy2520/Desktop/Research/kla

Rscript reanalysis/scripts/analyze_expanded_ddr_fraction_by_accession.R .
Rscript reanalysis/scripts/analyze_kla_regulator_intensity.R .
Rscript reanalysis/scripts/analyze_kla_regulator_whole_proteome_intensity.R .

Rscript reanalysis/tests/test_expanded_ddr_fraction_by_accession.R .
Rscript reanalysis/tests/test_kla_regulator_intensity.R .
Rscript reanalysis/tests/test_kla_regulator_landscape.R .
Rscript reanalysis/tests/test_kla_regulator_whole_proteome_intensity.R .

python3 reanalysis/scripts/build_final_manifest.py
```

## 重要报告

- [FINAL_REPORT.md](/Users/gzy2520/Desktop/Research/kla/reanalysis/reports/FINAL_REPORT.md)
- [EXPANDED_DDR_FRACTION_ACCESSION_ONLY.md](/Users/gzy2520/Desktop/Research/kla/reanalysis/reports/EXPANDED_DDR_FRACTION_ACCESSION_ONLY.md)
- [KLA_REGULATOR_INTENSITY_HEATMAP_METHODS.md](/Users/gzy2520/Desktop/Research/kla/reanalysis/reports/KLA_REGULATOR_INTENSITY_HEATMAP_METHODS.md)
- [REFERENCE_PROTEOME_SELECTION.md](/Users/gzy2520/Desktop/Research/kla/reanalysis/reports/REFERENCE_PROTEOME_SELECTION.md)
- [DATA_SOURCE_AND_ANALYSIS_ALGORITHM.md](/Users/gzy2520/Desktop/Research/kla/reanalysis/reports/DATA_SOURCE_AND_ANALYSIS_ALGORITHM.md)

## 当前限制

- 普通全蛋白强度热图目前只有 16/40 个样本组可用；不能为缺失强度的组伪造连续信号。
- PXD037371 仍未进入 DDR 配对统计。
- PXD038880/PXD050906 仍按老师要求排除。
- PXD043880 的 symbol 是来源文件本身的格式限制；后续不得再直接以 symbol 作为分析键。
- `reanalysis/reports/final_file_manifest_sha256.csv` 是当前文件状态校验清单。
