# 当前30组Kla与普通全蛋白DDR占比审计

## 计算范围

DDR占比仅对当前30个配对`PXD+SampleGroup`计算。四分类顺序为非肿瘤组织9组、
肿瘤组织2组、癌细胞系12组和正常细胞系7组。

Kla与普通全蛋白分别计算：

```text
DDR占比 = DDR交集中的唯一BaseAccession数 / 本组全部唯一BaseAccession数
```

两侧分母独立。普通全蛋白数据不能用Kla/PTM富集信号替代。

## ID和GO规则

- 唯一分析键为去isoform的UniProt`BaseAccession`。
- DDR来自`data/annotations/GO-repair+damage(human).tsv`。
- 排除`NOT`限定符；主分析保留全部evidence code。
- GeneSymbol不用于匹配或回退。

## 当前输出

- 30组英文统计：
  `results/tables/cell_type_kla_vs_reference_ddr_statistics_accession_only_paired_30.csv`
- 30组中文统计：
  `results/tables/cell_type_kla_vs_reference_ddr_statistics_accession_only_paired_30_zh.csv`
- 绘图数据：
  `results/tables/cell_type_kla_vs_reference_ddr_plot_rows.csv`
- 英文图：
  `results/figures/cell_type_kla_vs_reference_ddr_fraction_en.png`
- 中文图：
  `results/figures/cell_type_kla_vs_reference_ddr_fraction_zh.png`

该图比较蛋白集合中的DDR比例，不是差异丰度分析，也不能直接解释为DNA损伤程度、
修复活性或通路激活。
