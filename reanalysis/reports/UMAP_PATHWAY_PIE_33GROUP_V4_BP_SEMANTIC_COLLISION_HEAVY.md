# 33组Kla∩DDR通路饼图UMAP V4：BP层级语义与防碰撞展示

## V3问题与V4修正

- V3直接混合BP/CC/MF并使用IDF；3,461个直接GO term中1,864个只命中一个蛋白，稀有term权重使个体区分压过共享功能，因此图上几乎没有清楚的功能簇。
- V4仍完整保留UniProt原始全GO表、13,738条蛋白–GO长表以及507 × 3,461直接GO二值矩阵，但最终通路布局只使用BP，因为当前饼图表达的是DNA修复通路/生物过程。
- CC和MF未被删除，仍保存在全GO归档表中；它们不再参与最终坐标，以免宽泛定位或分子功能跨过程连接不同BP簇。

## BP语义特征

- 注释来源：UniProt 2026_02（10-June-2026）；GO层级来源：GO.db 3.20.0。
- 每个直接BP term扩展到当前GO.db版本可提供的BP祖先节点，删除BP根`GO:0008150`。
- 有4个UniProt直接BP ID不在当前GO.db祖先映射中；这些直接term仍保留，但无法补充祖先。明细保存在`bp_direct_terms_missing_from_go_db.csv`。
- 只有命中至少2个、至多405个蛋白的共享特征参与UMAP；单蛋白term仍保存在原始表中，但不能贡献蛋白间相似性；覆盖>80%蛋白的近根term被排除。
- 最终UMAP矩阵为507 × 3008，使用binary cosine、`n_neighbors = 8`、`min_dist = 0.2`、`spread = 1.5`、1,000轮、种子25、单线程。
- 通路评分和33组检出模式均不参与UMAP。

## 饼图防碰撞

- 原始UMAP坐标单独保存，是功能结构解释使用的坐标。
- 主饼图先对原始坐标统一缩放，再执行确定性防碰撞；饼半径为0.75、直径为1.5、最小中心距为2.55，饼边缘最小净间距为1.05。间距越大，展示坐标对原始UMAP局部结构的改变越强。
- 防碰撞后小于设定中心距的蛋白对为0；原始与展示坐标全部蛋白对距离的Spearman相关为0.9533。
- 防碰撞会改变稠密区域中的局部邻居关系；主图只用于读取通路组合和宽尺度几何。局部功能拓扑必须使用另存的原始UMAP坐标/点图解释。

## 编码与配色

- `+1`为通路颜色实心扇形，并用完全相同的通路颜色描边以覆盖polygon闭合缝；扇形之间不加白色分隔线。`-1`为白色空心扇形并用通路色描边。
- 22个全0蛋白使用浅灰圈叉；Chromatin interaction为`#7E6148`，Other support为`#6F6F6F`。

## 保存的GO表和矩阵

- UniProt原始全GO返回表：`reanalysis/config/uniprot_kla_ddr_507_all_go_2026-08-10.tsv`
- UniProt下载元数据：`reanalysis/config/uniprot_kla_ddr_507_all_go_2026-08-10_metadata.tsv`
- 全部直接蛋白–GO长表：`reanalysis/results/tables/kla_ddr_pathway_pie_umap_33groups_v3_all_go/uniprot_direct_go_annotation_long.csv`
- 507 × 3,461全部直接GO矩阵：`reanalysis/results/tables/kla_ddr_pathway_pie_umap_33groups_v3_all_go/protein_all_go_direct_binary_matrix.csv`
- BP直接term到祖先节点展开表：`reanalysis/results/tables/kla_ddr_pathway_pie_umap_33groups_v4_bp_semantic/bp_direct_to_ancestor_expansion_long.csv`
- V4实际使用的BP语义矩阵：`reanalysis/results/tables/kla_ddr_pathway_pie_umap_33groups_v4_bp_semantic/protein_bp_semantic_shared_binary_matrix.csv`

## 图形输出

- 推荐饼图：`reanalysis/results/figures/kla_ddr_pathway_pie_umap_33groups_v4_bp_semantic/kla_ddr_pathway_pie_umap_v4_bp_semantic.{png,pdf,svg}`
- 未位移原始点图：`reanalysis/results/figures/kla_ddr_pathway_pie_umap_33groups_v4_bp_semantic/kla_ddr_raw_bp_semantic_umap_v4.{png,pdf,svg}`
- 原始与展示对照：`reanalysis/results/figures/kla_ddr_pathway_pie_umap_33groups_v4_bp_semantic/kla_ddr_raw_and_pie_bp_semantic_umap_v4.{png,pdf,svg}`
