# Kla-DDR通路分面UMAP（33组，V1）

## 固定范围

- 使用固定33组范围产生的507个Kla∩DDR蛋白。
- 分析与合并键为去isoform的UniProt `BaseAccession`；GeneSymbol仅用于显示。
- 所有面板复用V4 BP semantic UMAP的原始坐标，未执行防碰撞位移。
- UMAP坐标由直接UniProt BP term及GO.db中可获得的BP祖先构建；通路评分不参与UMAP。

## 图形编码

- 每个面板只突出一个通路/功能，并保留全部507个蛋白作为共同背景。
- `+1`：通路颜色实心点。
- `-1`：白色空心点，边框使用该通路颜色。
- `0`：浅灰色背景点，表示未分配至当前面板的通路/功能。
- 九宫格和九张单独通路图使用相同原始UMAP坐标范围，因此可以直接比较空间分布。

## 产物

- 九宫格主图：`results/figures/pathway_specific_umap/kla_ddr_pathway_specific_umap_3x3_v1.png`
- 九张单独通路图目录：`results/figures/pathway_specific_umap/individual_pathways`
- 完整4,563行绘图表：`results/tables/pathway_specific_umap/pathway_specific_umap_plot_data.csv`
- 通路计数表：`results/tables/pathway_specific_umap/pathway_specific_umap_summary.csv`
- 输入文件审计：`results/tables/pathway_specific_umap/input_file_audit.csv`
