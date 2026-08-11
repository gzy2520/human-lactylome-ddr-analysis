# 五集合 × 三种降维 × 九通路图（固定33组，V1）

## 固定分析范围

- 使用固定33组产生的507个Kla∩DDR蛋白。
- 分析键为去isoform的UniProt `BaseAccession`；GeneSymbol仅用于显示。
- 五个展示集合为全部507、正常/非肿瘤组织183、正常/非肿瘤细胞471、癌症组织178、癌症细胞383；四类成员来自现有`kla_ddr_four_class_venn/membership.csv`，同一蛋白可在多个类别中被检出。
- 三种降维均使用完全相同的507 × 3,008二元BP语义特征矩阵。通路评分和类别检出均不参与坐标计算。

## 降维配置

- UMAP：复用已审计的V4原始坐标；cosine、`n_neighbors = 8`、`min_dist = 0.2`、`spread = 1.5`、1,000轮、seed 25。
- t-SNE：同一矩阵的预计算cosine距离；perplexity 30、theta 0.5、1,000轮、单线程、seed 25。
- PCA：同一二元矩阵中心化但不按方差缩放；PC符号按最大绝对蛋白得分确定性定向。PCA本身无随机性，流程仍记录seed 25。
- 每种算法只计算一套507蛋白坐标，全部507和四个类别的九通路面板均复用；因此同一算法内可直接比较五个集合。

## 图形编码

- 通路色实心：该类别检出且评分`+1`。
- 通路色空心：该类别检出且评分`-1`。
- 中灰：该类别检出，但未分配至当前面板通路。
- 极浅灰：未在当前生物类别中检出，仅作为全局坐标背景。

## 产物

- 15张九宫格图清单：`reanalysis/results/tables/kla_ddr_five_set_three_embedding_pathway_grids_33groups_v1/figure_manifest_15_grids.csv`
- 三种降维坐标：`reanalysis/results/tables/kla_ddr_five_set_three_embedding_pathway_grids_33groups_v1/embedding_coordinates_507.csv`
- 降维配置：`reanalysis/results/tables/kla_ddr_five_set_three_embedding_pathway_grids_33groups_v1/embedding_parameters.csv`
- 五集合成员长表：`reanalysis/results/tables/kla_ddr_five_set_three_embedding_pathway_grids_33groups_v1/embedding_set_membership_507_long.csv`
- 22,815行集合×通路绘图表：`reanalysis/results/tables/kla_ddr_five_set_three_embedding_pathway_grids_33groups_v1/embedding_set_pathway_plot_data.csv`
- 集合×通路计数：`reanalysis/results/tables/kla_ddr_five_set_three_embedding_pathway_grids_33groups_v1/embedding_set_pathway_summary.csv`
- 输入文件审计：`reanalysis/results/tables/kla_ddr_five_set_three_embedding_pathway_grids_33groups_v1/input_file_audit.csv`
