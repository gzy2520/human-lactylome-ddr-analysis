# 乳酸化调控基因与 DDR 蛋白占比：逐样本重建版

后续核验发现本版TALL的12.71%漏拆PEAKS多accession，完整解析均值为11.421069%。当前更正图见`deliverables/20261011_lactylation_ddr_tall_verified/`，算法对照见`audit/20261011_tall_algorithm_history/README.md`。本文件保留历史记录。

2026-10-11。交付主文件：`Combined_RNA_DDR.pdf`，共4页；各页独立PDF及300 dpi PNG位于`figures/`。

## 图形内容

沿用确认的布局：上方为RNA分组箱线图与均值虚线，中间为31个材料名称，下方为向下延伸的DDR占比柱及误差线。4页分别为12基因总览、乳酸代谢与KAT、乳酸代谢与AARS、去乳酸化调控基因。

RNA包含1898条唯一样本记录、28个参考来源、12个Ensembl基因，共22776个表达值。HCT116三列和HK-2两列复用RNA参考，以星号标注；具体来源及实验条件见`tables/display_mapping_31.csv`。LDHA、LDHB的3796个原表达值全部保持不变。RNA与蛋白按材料对齐，用于比较材料层面的表达和蛋白组成。

## DDR计算

1. 从当前项目的原始处理表提取每个样本或来源单位检出的蛋白；分级、重复采集先按样本或技术制备合并。
2. 按去除isoform后缀的稳定UniProt accession去重，逐accession与冻结GO-DDR集合匹配。冻结集合含2660个人源、非NOT的UniProt accession。
3. 每单位计算`100 × DDR蛋白数 / 检出蛋白数`；乳酸化组的分母仅包括有乳酸化鉴定证据的蛋白。
4. 每个显示位置的柱高为单位百分比的算术均值，SD由单位百分比计算，`SEM = SD / sqrt(n)`。该均值由逐样本比例汇总获得。

最终325条来源单位显示记录覆盖62个材料×检测类型位置；HCT116的两个全蛋白组位置共用3个技术制备，因此显示记录数包含参考复用。40个位置绘制SEM，22个位置保留均值柱而不绘制SEM。

- `t`：7个全蛋白组位置，误差线来自每细胞系3个技术制备；每制备的跨仪器重复采集先合并。
- `m`：NSC与GBM两材料、两检测类型共4个位置，柱高为不同细胞模型的均值；模型内重复未提供，图中不画SEM。
- `q`：HEK293T全蛋白组QC lys来源，401次采集合并为1个单位。
- 单一来源、单一样本，以及每条件仅1条记录的材料不画SEM。其他误差线描述纳入来源记录的分布；TALL的3条记录对应3个原始分析文件，独立供体信息未提供。

保留31组原有材料和实验条件范围。包含多个有重复条件的材料按所有纳入样本汇总，例如HepG2的WT/KO、HUVEC的0/6小时以及RKO的WT/KO；条件分别列于`ddr_source_unit_values.csv`，条件内汇总另存`ddr_condition_summary.csv`。

## 与最初图的核对

最初图来源表有272条记录、39个材料×检测类型位置。269条记录的蛋白分母和DDR分子逐项完全一致；只有TALL乳酸化的3条记录需要修正。新增53条显示记录补齐其余23个位置。

TALL旧解析从富集蛋白表计入全部有Area的蛋白，包含未乳酸化蛋白。当前按Lac PSM、-10lgP≥20、排除decoy，并校验同一accession的肽段位置映射，三个样本分别为21/173、24/178、22/176，均值12.707291%。旧均值约1.79%。按照原论文AScore≥20定位阈值另作敏感性核对，均值12.328074%；该结果保存在审计表，当前图保留历史解析的定位筛选规则以便逐项比较。

原PC-3M乳酸化均值13.765432%和ESCC全蛋白组均值5.533467%均复算一致。逐样本比较见`tables/original_sample_comparison.csv`；逐柱比较见`tables/original_figure_comparison.csv`。新增来源的标识符修正和采集合并依据见审计目录：肺/胎盘使用逐UniProt直接GO判定；MCF10A保留81个旧正则遗漏的合法10字符accession。

旧交接文档中的部分分类已由原表和元数据修正：adjacent liver/HCC的Kla是同一患者的两类单样本组织；HEK的HCD24–40为碰撞能量；HUVEC与RKO原表存在可拆分样本。当前版与历史组并集柱图分别保存。

## 数据与复现

- `tables/ddr_source_unit_values.csv`：每单位的样本ID、材料、检测类型、稳定ID计数、百分比、重复类型、条件及来源文件。
- `inputs/ddr_protein_membership.csv.gz`：逐单位UniProt和DDR判断，可直接复算全部柱。
- `tables/ddr_mean_sem_62.csv`：62根柱的n、均值、SD、SEM及实际绘图误差`PlotSEM`。
- `tables/ddr_condition_summary.csv`：按条件汇总，保留样本分组。
- `tables/gene_panel.csv`：12基因的Ensembl、HGNC、Entrez及显示名称。
- `SHA256SUMS.tsv`：本交付文件校验值。

从Git工作区根目录运行：

```sh
Rscript workflow/lactylation_ddr/rebuild_sample_inputs.R
Rscript workflow/lactylation_ddr/render.R
```

原表恢复代码和证据位于`audit/20261011_ddr_sample_recovery/`，详见其中文README。合并脚本从这些逐蛋白恢复结果重新计数；不会读取旧组并集比例作为柱高。绘图依赖data.table、ggplot2、patchwork、ragg；使用Arial Unicode MS、随机种子25、Cairo PDF和300 dpi PNG。软件版本见`sessionInfo.txt`。

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
