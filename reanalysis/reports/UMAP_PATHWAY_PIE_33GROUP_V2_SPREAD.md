# 33组Kla∩DDR通路饼图UMAP V2：扩展间距与负向空心扇形

## 本版修改

- 不再使用V1固定坐标；从同一507 × 66原始GO-term二值矩阵重新拟合UMAP。
- 通路评分仍不参与UMAP，仅在坐标生成后用于绘制饼图。
- 新参数：cosine距离、`n_neighbors = 15`、`min_dist = 3`、`spread = 10`、`repulsion_strength = 1.8`、`negative_sample_rate = 10`、750轮、随机种子25、单线程。
- 饼图半径为0.95坐标单位。以该直径计算，V1中有5330对蛋白距离小于一个饼图直径；V2降为9对。
- 102个仅命中`GO:0006974`的蛋白拥有完全相同的输入特征；它们在V2中的相对摊开仅用于显示，不能解释为生物学距离。

## 扇形编码

- `+1`促进性评分：使用对应通路颜色的实心扇形，扇形之间以细白线分隔。
- `-1`抑制性评分：使用白色空心扇形，并以对应通路颜色描边。
- `0`：不绘制该通路扇形；9列全为0的22个蛋白显示为灰色圆点。
- 每个非零通路扇形等权；扇形面积不表示表达量、蛋白丰度或证据强弱。

## 配色与版式

- 使用NPG离散色板中9个颜色，该色板源自Nature Reviews Cancer图形配色。
- 调色板来源：`https://github.com/nanxstats/ggsci/blob/master/R/palettes.R`。
- 图中文字使用统一无衬线字体，并同时导出600 dpi PNG及PDF/SVG矢量文件，便于按Nature双栏宽度缩放。
- Nature作者指南：`https://www.nature.com/nature/for-authors/final-submission`。

## 输出

- 推荐饼图：`reanalysis/results/figures/kla_ddr_pathway_pie_umap_33groups_v2_spread/kla_ddr_pathway_hollow_negative_pie_umap_v2_spread.{png,pdf,svg}`
- 组合图：`reanalysis/results/figures/kla_ddr_pathway_pie_umap_33groups_v2_spread/kla_ddr_pathway_hollow_negative_umap_and_summary_v2_spread.{png,pdf,svg}`
- V2坐标：`reanalysis/results/tables/kla_ddr_pathway_pie_umap_33groups_v2_spread/umap_coordinates_v2_spread.csv`
- 参数：`reanalysis/results/tables/kla_ddr_pathway_pie_umap_33groups_v2_spread/umap_v2_parameters.csv`
- 重叠比较：`reanalysis/results/tables/kla_ddr_pathway_pie_umap_33groups_v2_spread/pie_overlap_comparison_v1_v2.csv`
