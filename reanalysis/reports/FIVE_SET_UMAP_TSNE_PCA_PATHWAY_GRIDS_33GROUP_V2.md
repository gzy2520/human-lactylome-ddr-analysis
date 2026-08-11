# 五个蛋白集合 × 三种独立降维 × 九通路图（固定33组，V2）

## 分析范围

- 第一个集合为固定33组产生的全部507个Kla∩DDR蛋白；其余四个集合分别为正常/非肿瘤组织183、正常/非肿瘤细胞471、癌症组织178、癌症细胞383。
- 每个集合仅使用自身包含的蛋白重新学习UMAP、t-SNE和PCA坐标，共5 × 3 = 15张互不共用坐标的九宫格。
- 分析键为去isoform的UniProt `BaseAccession`；GeneSymbol仅用于显示。
- 所有集合使用相同的3,008个BP语义特征定义；通路评分不参与任何降维。

## 参数筛选与正式配置

- 使用五个蛋白集合上的10近邻保持率、高低维全局距离Spearman相关及归一化最近邻距离第5百分位联合评价候选参数；综合评分权重依次为0.70、0.20、0.10。
- UMAP正式配置：cosine、`n_neighbors = 12`、`min_dist = 0.35`、`spread = 1.8`、`repulsion_strength = 1.5`、`negative_sample_rate = 10`、1,000轮、seed 25、单线程。
- t-SNE正式配置：预计算cosine距离、perplexity 25、theta 0.5、1,000轮、seed 25、单线程。
- PCA：对每个蛋白的BP二元向量先作L2归一化，再按特征中心化但不按方差缩放；使用中心化蛋白×蛋白Gram矩阵的特征分解求取与标准PCA等价的样本得分，PC符号确定性定向。该处理使PCA与UMAP/t-SNE采用的cosine语义更一致，并减弱GO注释条目数量的影响。
- 三种方法均记录seed 25；PCA本身为确定性计算。

## 图形解释

- 通路色实心：评分`+1`；通路色空心：评分`-1`；中灰：未分配至当前面板通路。
- 每张图只绘制该集合自身蛋白，不再保留集合外蛋白作为浅灰背景。
- 同一张九宫格内九个面板严格复用同一套坐标；不同蛋白集合之间为独立拟合，不能比较绝对坐标轴、方向或单个蛋白的绝对位置。

## 产物

- 15张九宫格图清单：`reanalysis/results/tables/kla_ddr_five_set_three_embedding_pathway_grids_33groups_v2/figure_manifest_15_grids.csv`
- 五套三种降维坐标：`reanalysis/results/tables/kla_ddr_five_set_three_embedding_pathway_grids_33groups_v2/embedding_coordinates_5sets_long.csv`
- 各集合完整参数：`reanalysis/results/tables/kla_ddr_five_set_three_embedding_pathway_grids_33groups_v2/embedding_parameters_5sets.csv`
- 五集合蛋白成员表：`reanalysis/results/tables/kla_ddr_five_set_three_embedding_pathway_grids_33groups_v2/embedding_set_membership.csv`
- 集合×通路绘图表：`reanalysis/results/tables/kla_ddr_five_set_three_embedding_pathway_grids_33groups_v2/pathway_plot_data_5sets.csv`
- 集合×通路计数：`reanalysis/results/tables/kla_ddr_five_set_three_embedding_pathway_grids_33groups_v2/pathway_summary_5sets.csv`
- 参数候选评价：`reanalysis/results/tables/kla_ddr_five_set_three_embedding_pathway_grids_33groups_v2/parameter_tuning_candidate_summary.csv`
- 参数筛选原始指标：`reanalysis/results/tables/kla_ddr_five_set_three_embedding_pathway_grids_33groups_v2/parameter_tuning_metrics.csv`
- 输入文件审计：`reanalysis/results/tables/kla_ddr_five_set_three_embedding_pathway_grids_33groups_v2/input_file_audit.csv`
