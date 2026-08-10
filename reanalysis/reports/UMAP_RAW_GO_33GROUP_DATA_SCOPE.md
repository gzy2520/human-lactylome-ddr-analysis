# 33组项目范围与Kla∩DDR原始GO-term UMAP数据审计

## 当前项目正式范围

- 原始Kla成员表含37个可量化样本组。
- 本次及后续正式分析按`PXD + SampleGroup`精确排除4个不可用组，得到33组。
- 正式33组与`four_class_sample_grouping.csv`逐组完全一致：正常细胞9组、癌细胞13组、正常组织9组、癌组织2组。
- `PXD037371`的3个临床组此前已因TMT通道无法可靠映射而排除；它们不在上述37组成员源表中，也不计入本次“37减4等于33”。

本次从37组中排除的4组为：

1. `PXD062720 / bladder cancer cells treated with EPI`
2. `PXD063047 / severe preeclampsia placenta`
3. `PXD064038 / MEC and NEC ESCC groups`
4. `PXD075014 / AC16 control and hypoxia`

## 本次UMAP选择的数据

- 分析对象：33组中所有`IsDdr == TRUE`的Kla蛋白并集。
- 唯一分析键：去除isoform后缀的UniProt `BaseAccession`；不以Gene Symbol匹配、去重或建模。
- 一个点代表一个唯一`BaseAccession`，共507个点。
- GO来源：`data/annotations/GO-repair+damage(human).tsv`。
- GO过滤：仅保留`GENE PRODUCT DB == UniProtKB`、人源`TAXON ID == 9606`、不含`NOT`限定词且GO term非空的记录。
- 同一蛋白与同一GO term的重复记录折叠为一个二值命中。
- 507个蛋白全部至少有1个合格GO term；最终为66个原始GO term、1,029个唯一蛋白–GO term配对。
- 507个蛋白形成205种不同的原始GO成员模式；最大的一组为102个仅命中`GO:0006974`的蛋白。
- UMAP的唯一输入为507 × 66的“蛋白 × 原始GO term”0/1矩阵。
- 不对GO term做HR、NHEJ、BER等类别归并；样本类别、蛋白出现组数、表达/强度和Gene Symbol均不进入UMAP特征矩阵。

## UMAP参数与后续使用

- R包：`uwot`。
- 距离：cosine，适用于当前稀疏二值成员矩阵。
- 参数：`n_neighbors = 30`、`min_dist = 0.8`、`spread = 3.0`、`n_epochs = 500`。
- 使用随机初始化并固定随机种子25；最近邻与SGD均单线程，以保证当前环境内可重复。
- 当前图中所有点为同一种蓝色，不编码通路。
- 后续获得通路归属后，应直接按`BaseAccession`合并到`umap_coordinates_fixed.csv`并在固定坐标上将点替换为饼图；不得重新拟合UMAP。
- 具有完全相同GO成员模式的蛋白在输入空间中不可区分；它们在二维图内的小幅分散不代表额外生物学差异。

## 输出边界

- 本流程只新建`kla_ddr_raw_go_umap_33groups`目录及本审计文档。
- 未重跑、未覆盖既有Venn图、DDR柱状图或`previous_umap`。
- UMAP是基于GO注释相似性的二维可视化；坐标轴及点间全局距离不应解释为定量生物学效应。

## 关键文件

- 正式33组：`reanalysis/results/tables/kla_ddr_raw_go_umap_33groups/project_sample_group_scope_33.csv`
- 4个排除组：`reanalysis/results/tables/kla_ddr_raw_go_umap_33groups/removed_sample_groups_4.csv`
- 507个蛋白：`reanalysis/results/tables/kla_ddr_raw_go_umap_33groups/kla_ddr_protein_scope.csv`
- 原始GO长表：`reanalysis/results/tables/kla_ddr_raw_go_umap_33groups/protein_raw_go_term_long.csv`
- 原始GO模式汇总：`reanalysis/results/tables/kla_ddr_raw_go_umap_33groups/raw_go_pattern_summary.csv`
- 507 × 66二值矩阵：`reanalysis/results/tables/kla_ddr_raw_go_umap_33groups/protein_raw_go_term_binary_matrix.csv`
- 固定坐标：`reanalysis/results/tables/kla_ddr_raw_go_umap_33groups/umap_coordinates_fixed.csv`
- UMAP参数：`reanalysis/results/tables/kla_ddr_raw_go_umap_33groups/umap_parameters.csv`
- 单色图：`reanalysis/results/figures/kla_ddr_raw_go_umap_33groups/kla_ddr_raw_go_umap_33groups_solid_blue.{png,pdf,svg}`
