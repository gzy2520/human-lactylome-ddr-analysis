# 33组Kla∩DDR通路评分与固定UMAP展示审计

## 数据合并

- 评分来源：`data/identifier/260810乳酸化DDR基因评分表.xlsx`的“评分表”工作表。
- 工作表含507个有`BaseAccession`的正式评分行，另有7个仅保留ID的空白模板行；空白行不参与分析。
- 507个评分蛋白与`umap_coordinates_fixed.csv`中的507个蛋白按去isoform的UniProt `BaseAccession`逐一完全匹配。
- `GeneSymbol`和蛋白名称仅保留用于显示与审计，不作为合并、去重或分析键。
- UMAP坐标直接读取既有固定坐标，未重新拟合。

## 通路评分

- 9个通路/功能类别：HR、NHEJ、AEJ、BER、NER、MMR、FA、Chromatin Interaction和Other support。
- 分值仅为`-1/0/+1`：`+1`表示促进性证据，`-1`表示抑制性证据，`0`表示未分配该通路。
- 共485个蛋白至少有1个非零分配；22个蛋白9列均为0。
- 共1,175个非零蛋白–通路分配，其中1,108个为`+1`，67个为`-1`；53个蛋白至少含1个`-1`分配。
- “参考文献表”恰好有1,175行，并可按`ID + Pathway + Score`与每个非零分配一对一匹配。

## 展示方案

- 主图保留固定UMAP坐标，并把每个有分配的蛋白画成等权饼图；扇区面积表示通路成员关系，不表示表达量、强度或评分大小。
- 实色扇区表示`+1`，同色浅色扇区表示`-1`；全零蛋白显示为灰点。
- 由于饼图中的浅色负向扇区在密集区域较小，推荐使用“UMAP饼图 + 正负分配计数条形图”的组合图作为主展示；同时保留纯UMAP饼图。
- Chromatin Interaction和Other support是功能类别，不应在文字中误称为经典DNA修复通路。

## 输出

- 推荐组合图：`reanalysis/results/figures/kla_ddr_pathway_pie_umap_33groups/kla_ddr_pathway_umap_and_signed_summary_33groups.{png,pdf,svg}`
- 纯UMAP饼图：`reanalysis/results/figures/kla_ddr_pathway_pie_umap_33groups/kla_ddr_pathway_signed_pie_umap_33groups.{png,pdf,svg}`
- 固定坐标绘图数据：`reanalysis/results/tables/kla_ddr_pathway_pie_umap_33groups/pathway_umap_plot_data_fixed_coordinates.csv`
- 带文献证据的1,175条分配：`reanalysis/results/tables/kla_ddr_pathway_pie_umap_33groups/pathway_assignments_with_evidence_1175.csv`
