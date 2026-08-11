# 五集合UMAP、t-SNE和PCA九宫格三版本对照

## 共同分析范围

- 三版均使用固定33组产生的507个Kla∩DDR蛋白。
- 五个展示集合均为：全部507、正常/非肿瘤组织183、正常/非肿瘤细胞471、癌症组织178、癌症细胞383。
- 分析键为去isoform的UniProt `BaseAccession`；GeneSymbol仅用于显示。
- 坐标拟合均基于BP语义特征，通路评分不参与降维。
- 通路色实心为促进`+1`，通路色空心为抑制`-1`，中灰为未分配至当前通路。

## 三版的核心区别

| 版本 | 坐标拟合范围 | 参数策略 | 集合外蛋白 | 主要用途 |
|---|---|---|---|---|
| V1 | 每种算法仅对全部507拟合一次，五个集合复用同一坐标 | 固定初始配置 | 四个类别图中以极浅灰保留 | 在严格相同坐标系中比较类别覆盖位置 |
| V2 | 五个集合分别重新拟合UMAP、t-SNE和PCA | 每种算法在五个集合中使用同一套整体最优参数 | 不绘制 | 比较独立拟合效果，控制参数一致 |
| V3 | 五个集合分别重新拟合UMAP、t-SNE和PCA | 每个集合、每种算法分别调参 | 不绘制 | 优先改善各集合内部的结构保持和可读性 |

## 解释边界

- V1同一算法的五张图复用坐标，因此可以比较不同集合在同一全局点云中的覆盖位置。
- V2和V3的五个集合均独立拟合，不能直接比较集合之间的绝对轴、方向、簇位置或单个蛋白绝对坐标。
- 三版不是三次独立生物学分析，而是同一数据和通路分配在不同坐标策略下的可视化比较。
- V3参数选择兼顾局部邻域保持、全局距离、点间距、主体面积和轴向均衡；没有人为移动蛋白坐标。

## 图形目录

- V1：`reanalysis/results/figures/kla_ddr_five_set_three_embedding_pathway_grids_33groups_v1/`
- V2：`reanalysis/results/figures/kla_ddr_five_set_three_embedding_pathway_grids_33groups_v2/`
- V3：`reanalysis/results/figures/kla_ddr_five_set_three_embedding_pathway_grids_33groups_v3/`

每个目录均包含5个集合 × 3种算法，共15张九宫格，并同时输出PNG、PDF和SVG。
