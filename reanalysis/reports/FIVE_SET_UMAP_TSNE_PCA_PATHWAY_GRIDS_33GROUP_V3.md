# 五集合逐算法独立调参的九通路图（固定33组，V3）

## 与V2的区别

- V2为五个集合独立计算坐标，但同一算法仍使用统一参数；V3进一步为每个集合分别选择UMAP、t-SNE和PCA配置。
- 共5个集合 × 3种算法 = 15张九宫格；每张图只包含本集合蛋白，九个通路面板复用该图坐标。
- V1和V2结果均保留，V3不覆盖旧产物。

## 参数选择原则

- UMAP和t-SNE的逐集合综合评分：10近邻保持率0.50、全局距离Spearman相关0.15、最近邻距离第5百分位0.20、核心点云面积占比0.10、稳健轴向均衡度0.05。
- PCA的逐集合综合评分：10近邻保持率0.30、全局距离Spearman相关0.40、最近邻距离第5百分位0.15、核心点云面积占比0.10、稳健轴向均衡度0.05。
- 核心点云面积指标用于避免少数远端蛋白把主体压缩在画布很小区域；它只用于选择算法参数，不移动任何蛋白坐标。
- 所有随机算法固定seed 25并使用单线程；通路评分不参与坐标计算。

## 最终逐集合配置

- `all_507` UMAP：n_neighbors=18，min_dist=0.25，repulsion=1.5。
- `all_507` t-SNE：perplexity=20，theta=0.3。
- `all_507` PCA：l2_all_features，保留3008个特征。
- `normal_tissue` UMAP：n_neighbors=8，min_dist=0.55，repulsion=2.0。
- `normal_tissue` t-SNE：perplexity=40，theta=0.5。
- `normal_tissue` PCA：l2_all_features，保留3008个特征。
- `normal_cells` UMAP：n_neighbors=18，min_dist=0.25，repulsion=1.5。
- `normal_cells` t-SNE：perplexity=25，theta=0.3。
- `normal_cells` PCA：l2_subset_filter_02_90，保留2968个特征。
- `cancer_tissue` UMAP：n_neighbors=12，min_dist=0.55，repulsion=1.5。
- `cancer_tissue` t-SNE：perplexity=30，theta=0.5。
- `cancer_tissue` PCA：l2_subset_filter_02_90，保留1688个特征。
- `cancer_cells` UMAP：n_neighbors=12，min_dist=0.70，repulsion=2.0。
- `cancer_cells` t-SNE：perplexity=25，theta=0.3。
- `cancer_cells` PCA：l2_subset_filter_02_90，保留2418个特征。

## 图形编码与解释边界

- 通路色实心：评分`+1`；通路色空心：评分`-1`；中灰：未分配至当前面板通路。
- 分析键为去isoform的UniProt `BaseAccession`；GeneSymbol仅用于显示。
- 不同集合使用不同参数且独立拟合，不能直接比较绝对坐标轴、方向或单个蛋白的绝对位置。

## 产物

- 15张九宫格清单：`reanalysis/results/tables/kla_ddr_five_set_three_embedding_pathway_grids_33groups_v3/figure_manifest_15_grids.csv`
- 五套调参坐标：`reanalysis/results/tables/kla_ddr_five_set_three_embedding_pathway_grids_33groups_v3/embedding_coordinates_5sets_tuned_long.csv`
- 最终参数长表：`reanalysis/results/tables/kla_ddr_five_set_three_embedding_pathway_grids_33groups_v3/embedding_parameters_5sets_tuned.csv`
- 逐集合推荐：`reanalysis/results/tables/kla_ddr_five_set_three_embedding_pathway_grids_33groups_v3/parameter_tuning_recommendation_by_set.csv`
- 候选参数排名：`reanalysis/results/tables/kla_ddr_five_set_three_embedding_pathway_grids_33groups_v3/parameter_tuning_ranked_candidates_by_set.csv`
- 参数筛选原始指标：`reanalysis/results/tables/kla_ddr_five_set_three_embedding_pathway_grids_33groups_v3/parameter_tuning_metrics_by_set.csv`
- 输入审计：`reanalysis/results/tables/kla_ddr_five_set_three_embedding_pathway_grids_33groups_v3/input_file_audit.csv`
