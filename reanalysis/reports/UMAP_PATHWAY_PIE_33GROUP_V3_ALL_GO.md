# 33组Kla∩DDR通路饼图UMAP V3：全GO功能空间

## 为什么改为全GO

- 正式分析对象仍是33组数据得到的507个Kla∩DDR蛋白；一个饼代表一个去isoform的UniProt `BaseAccession`。
- V1/V2仅使用66个DDR相关GO term，205种输入模式中最大重复组包含102个蛋白，难以支撑可读的功能布局。
- V3从UniProt 2026_02（10-June-2026）取得这507个蛋白的全部直接BP/CC/MF注释：13,738条蛋白–GO关系、3,461个不同GO term。
- 全GO直接注释形成506种模式；仅RFC2（P35250）与RFC5（P40937）完全相同。二者属于同一Replication factor C复合体近缘小亚基，保留接近具有生物学合理性。

## UMAP输入和边界

- 使用507 × 3,461蛋白×GO二值矩阵；不对GO term预先归入HR/NHEJ等通路。
- 纳入BP、CC和MF直接注释，不扩展GO祖先节点，避免人为重复加入宽泛上位term。
- 对每个GO term使用逆文档频率权重：`IDF = log((507 + 1) / (注释蛋白数 + 1)) + 1`，再以cosine距离拟合UMAP；常见的nucleus、cytoplasm等term因此影响较低。
- 9列通路评分和33组检出模式均不参与UMAP，只在坐标生成后用于饼图编码或范围核验。
- UMAP邻近表示全GO注释相似性，不代表物理互作、共同表达、乳酸化丰度或通路激活强度。

## 扇形、全0蛋白与配色

- `+1`：对应通路颜色的实心扇形；`-1`：白色空心扇形并以对应通路颜色描边。
- 9列全为0的22个蛋白改为浅灰色空心圆内带叉（circled cross），明确表示未获当前评分表的通路/功能分配，不再画成类似第10类通路的实心灰点。
- Chromatin interaction保留棕色`#7E6148`；Other support改为中性炭灰`#6F6F6F`，消除原两个棕色过近的问题。
- 配色以NPG/Nature Reviews Cancer风格为基础，并使用中性炭灰扩展；导出600 dpi PNG和PDF/SVG矢量图。

## 可读性

- 饼半径为0.95坐标单位；V3有0对蛋白中心距离小于一个饼直径，涉及0个蛋白。

## 主要输出

- 推荐图：`reanalysis/results/figures/kla_ddr_pathway_pie_umap_33groups_v3_all_go/kla_ddr_pathway_hollow_negative_pie_umap_v3_all_go.{png,pdf,svg}`
- 组合图：`reanalysis/results/figures/kla_ddr_pathway_pie_umap_33groups_v3_all_go/kla_ddr_pathway_hollow_negative_umap_and_summary_v3_all_go.{png,pdf,svg}`
- 坐标：`reanalysis/results/tables/kla_ddr_pathway_pie_umap_33groups_v3_all_go/umap_coordinates_v3_all_go.csv`
- 全GO长表：`reanalysis/results/tables/kla_ddr_pathway_pie_umap_33groups_v3_all_go/uniprot_direct_go_annotation_long.csv`
- 全GO二值矩阵：`reanalysis/results/tables/kla_ddr_pathway_pie_umap_33groups_v3_all_go/protein_all_go_direct_binary_matrix.csv`
- 注释与IDF：`reanalysis/results/tables/kla_ddr_pathway_pie_umap_33groups_v3_all_go/all_go_term_dictionary_and_idf.csv`
- 参数：`reanalysis/results/tables/kla_ddr_pathway_pie_umap_33groups_v3_all_go/umap_v3_parameters.csv`
