# 乳酸化调控基因表达与 DDR 蛋白占比组合图

## 图件

`Combined_RNA_DDR.pdf` 包含四页，均采用中间材料横轴、上方 RNA 箱线图、下方 DDR 占比柱状图：

1. `Overview_12genes`：12 基因总览。
2. `Lactate_and_KATs`：LDHA、LDHB、EP300、CREBBP、KAT5。
3. `Lactate_and_AARS`：LDHA、LDHB、AARS1、AARS2。
4. `Delactylation_regulators`：HDAC1、HDAC2、HDAC3、SIRT2、SIRT3。

各页的独立 PNG（300 dpi）与矢量 PDF 位于 `figures/`。全部图使用相同的 31 列材料顺序，便于对照。分组图使每种材料中的基因箱线图更易分辨。

## 图注与计算

上方面板显示 qsmooth B 校正后的 log2(TPM + 0.5) 表达值。箱体为第一至第三四分位数，中线为中位数，须延伸至 1.5 倍 IQR 范围内最远观测；点为 RNA 样本记录。不同颜色表示不同基因，虚线连接每个材料的表达均值。

下方面板的蓝色、橙色柱分别表示全蛋白组、乳酸化蛋白组的 GO-DDR 注释蛋白占比均值，误差线为 SEM（样本数大于 1 时）。单个样本的占比为 `100 × DDR 注释蛋白数 / 该样本检出的蛋白数`。柱从中间向下延伸，刻度保持正百分比。原图缺少样本级观测的位置继续留空；每列每种组学的样本数见 `tables/ddr_bars_mean_sem.csv`。

RNA 与蛋白面板通过材料对应表对齐。RNA 共 1,898 个唯一样本记录、28 个 ReferenceKey；HCT116 的三列、HK-2 的两列分别复用同一套 RNA 参考数据，并以星号标识。RNA 来源包括 GTEx、TCGA 和其他参考队列，具体来源与处理条件见 `tables/display_mapping_31.csv` 的 ReferenceKey 和 RNALabel。例如两列 HK-2 使用 untreated RNA 参考；BPH RNA 来自用药队列。图用于材料间表达与 DDR 组成的并列比较，未计算跨组学样本相关系数。

## 基因选择依据

| 基因 | 纳入理由 | 原始研究 |
|---|---|---|
| LDHA、LDHB | 乳酸代谢，保留原图两基因 | [Cell 2024](https://doi.org/10.1016/j.cell.2023.11.022) |
| EP300 | 经典组蛋白乳酸化转移酶 | [Nature 2019](https://doi.org/10.1038/s41586-019-1678-1) |
| CREBBP（CBP） | MRE11 乳酸化，与同源重组修复直接相关 | [Cell 2024](https://doi.org/10.1016/j.cell.2023.11.022) |
| KAT5（TIP60） | NBS1 乳酸化，与 MRN 复合物和 DNA 修复相关 | [Nature 2024](https://doi.org/10.1038/s41586-024-07620-9) |
| AARS1、AARS2 | ATP 依赖的乳酸化转移酶 | [Cell 2024](https://doi.org/10.1016/j.cell.2024.04.002)、[Nature 2024](https://doi.org/10.1038/s41586-024-07992-y) |
| HDAC1、HDAC2、HDAC3 | 经典去乳酸化酶；HDAC3 有 NBS1 底物证据 | [Science Advances 2022](https://doi.org/10.1126/sciadv.abi6696)、[Nature 2024](https://doi.org/10.1038/s41586-024-07620-9) |
| SIRT2 | 组蛋白去乳酸化 | [Cell Discovery 2022](https://doi.org/10.1038/s41421-022-00398-y) |
| SIRT3 | CCNE2 去乳酸化，具有结构与生化证据 | [EMBO Reports 2023](https://doi.org/10.15252/embr.202256052) |

HDAC1–3 在特定乳酸条件下还被报道具有反向乳酸化催化能力，因此在图中使用调控酶名称，功能解释需结合具体底物与条件。[JBC 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC11888385/)

## 数据核对与复现

- 12 个基因均通过 HGNC 映射到唯一 Ensembl ID；矩阵提取和样本汇总按 Ensembl ID 完成，符号用于显示。
- 原 LDHA、LDHB 的 3,796 个表达值逐项保留，最大绝对差为 0。
- 原修正版蛋白数据 272 条全部保留，均为 `ObservationType=sample`；占比重算误差小于 1e-10。
- RNA 原始长表保留 12 × 1,898 = 22,776 条记录；31 列展开仅用于绘图。
- 原始输入路径与 SHA-256 位于 `inputs/original_source_manifest.csv`；本次提取的基因矩阵、对应表和蛋白输入均保存在 `inputs/`。
- 配色按基因区分；材料分类条保留历史四分类配色；全蛋白组与乳酸化组保持蓝、橙色。随机种子 25，字体 Arial Unicode MS。

在项目根目录运行 `Rscript workflow/lactylation_ddr/render.R` 可从已保存的输入及表格重建图件。`prepare_inputs.R` 用于从原始本地发布目录重新提取所需数据。R 依赖：data.table、ggplot2、patchwork、ragg；压缩表通过系统 gzip 读取。软件版本见 `sessionInfo.txt`。
