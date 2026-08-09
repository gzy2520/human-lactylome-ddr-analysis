# 方法

## 研究设计与方法组织

本研究整合公开的人源赖氨酸乳酸化（lysine lactylation, Kla）质谱数据，并为
每个可用的组织或细胞样本组配对材料身份和实验状态相符的非 Kla 富集普通全蛋白
质谱数据。方法部分的组织方式参考 Hu 等发表于 *Genome Biology* 的生信方法学
文章，即依次说明数据来源与预处理、分析判定规则、定量指标、可视化、计算环境
和可复现性；具体算法和文字均根据本项目重新撰写。

分析以 `PXD+SampleGroup` 为最小研究单元。来源审计中共有 37 个样本组具有可
追溯的 Kla 定量，其中 33 个具有完全匹配且可定量的普通全蛋白参照，进入最终
Kla/普通全蛋白热图、DDR 占比比较和四集合 Venn 分析。最终 33 组依次分为
正常/非肿瘤组织 9 组、癌症组织 2 组、正常/非肿瘤细胞 9 组和癌症细胞 13 组。

## 公共数据获取与文件管理

乳酸化数据及对应论文信息来自 ProteomeXchange/ProteomeCentral、PRIDE、iProX
和论文开放补充材料。每个 PXD 的文件分别保存于 `raw`、`search_results`、
`supplementary` 和 `metadata` 子目录。原始谱图仅作为来源数据保存，不直接
作为 Kla 蛋白证据；主分析优先读取作者已经完成数据库检索后导出的位点表、修饰
肽段表、蛋白表或作者补充表。所有实际使用文件、读取字段、样本列、文章 DOI 和
下载来源均记录在 `lactylome_reference_pairing.csv` 和老师审阅总表中。

对于 2019 年之前生成或使用泛乙酰化抗体富集的候选数据，不因富集方式本身将
蛋白判定为 Kla。只有检索结果或作者补充表明确报告 lactyl-lysine 修饰、修饰肽段
或可定位 Kla 位点时才纳入。

## Kla 证据提取与质量控制

不同搜索软件的 Kla 表示被转换为统一的 sample-level 长表。可接受证据包括
MaxQuant 的 `La/Lactyl (K)Sites` 和修饰特异肽段、PEAKS 中的 Lactylation、
`K(+72.02)` 或 `K(Lactyl)`、Spectronaut 的 Kla PTM/precursor quantity、
Proteome Discoverer 的乳酸化修饰肽，以及作者补充表中的明确 Kla 位点或肽段。
数据集特异的文件和字段规则由提取器固定，并保留来源文件、来源行或 site ID、
样本、实验组、修饰肽段、位点、定位概率、PEP、score、证据模式和来源置信度。

在读取阶段不合并不同样本、条件或重复。反库/decoy、reverse、potential
contaminant、非人源记录、无有效蛋白 ID、无明确 Kla 修饰或无法定位到蛋白的
记录被排除。PEAKS 数据另保留肽段鉴定分数和 diagnostic-ion 审计；diagnostic
ion 只能增强来源可信度，不能单独把未修饰肽判为 Kla。

主分析不对所有异质数据强制施加统一的位点定位概率阈值，而是保留作者报告的
有效 Kla 证据并执行上述基础质量控制。PXD014870 按老师要求使用额外定位概率
阈值 0 作为主分析，并单独生成 0.75 阈值的敏感性结果；敏感性结果不混入主集合。

## 蛋白标识符标准化

全部分析、去重、合并、GO 交集和 Venn 集合均使用人源 UniProt
`BaseAccession`。首先去除数据库前缀和 UniProt isoform 后缀，例如
`P12345-2` 转换为 `P12345`。来源为 Ensembl protein ID 的普通全蛋白数据先
通过项目缓存的 Ensembl-to-UniProt 映射表转换。reviewed/non-reviewed 状态、
GeneSymbol 和 Protein name 仅保留用于显示和人工审计，不参与命中、合并、
去重或缺失回退。最终分析的 GeneSymbol 回退数为 0。

同一蛋白在同一样本中因 isoform、重复记录、DDA/DIA 或亚细胞组分重复出现时，
先在可追溯的特征层级汇总，再按 BaseAccession 生成蛋白集合，避免重复计数。

## 普通全蛋白参照选择

普通全蛋白参照必须来自未进行 Kla/PTM 富集的蛋白定量矩阵。选择顺序为：

1. 同一研究、同一生物样本的非 PTM 普通全蛋白数据；
2. 同一研究、同一细胞或组织及相同处理状态的数据；
3. 独立研究中材料身份、解剖粒度、疾病状态和处理状态完全匹配的数据。

“组织相近”或“细胞来源相似”不视为严格匹配。每个 Kla 样本组的具体材料、
处理、参照 PXD、实际读取文件和样本列均经过材料身份审计。PXD062720、
PXD063047 的 severe preeclampsia placenta、PXD064038 和 PXD075014 因缺少
完全匹配且可审计的普通全蛋白强度参照，不进入最终配对分析；原始 Kla 数据和
排除记录仍保留。PXD037371 的三个肝组织临床组因 TMT 通道与临床分组无法可靠
对应，不进入 Kla 定量范围。

33 个 Kla 样本组对应 30 个唯一普通全蛋白参照行。HK-2、MCF7 和 HCT116
各有两个 Kla 研究共享完全相同的普通全蛋白数据；样本级配对表保留 33 条关系，
普通全蛋白热图只显示 30 条唯一参照，Kla 热图保留 33 行并按相应参照行连续排列。

## GO-DDR 蛋白集合与占比

DDR 注释来源为人源 `GO-repair+damage(human).tsv`。排除带 `NOT` 限定符的
注释，主分析不再按 GO evidence code 缩减，但保留 evidence code 供后续敏感性
分析。GO 中的 protein accession 同样转换为 BaseAccession。

对每个样本组分别生成 Kla 蛋白集合和普通全蛋白集合，并仅按 BaseAccession
与 GO-DDR 集合取交集。DDR 占比定义为：

```text
DDR蛋白占比（%） = 100 × DDR交集中的唯一BaseAccession数
                         / 该样本组全部唯一BaseAccession数
```

Kla 和普通全蛋白的分母分别独立计算。普通全蛋白参照改变时只重新计算参照侧
总蛋白数、DDR 数和占比，不改变 Kla 侧蛋白集合。

## 调控因子列表与多角色注释

Writer、Eraser 和 Reader 蛋白来自
`乳酸化调控因子_Writer-Eraser-Reader.xlsx`，并通过项目映射表转换为人源
UniProt BaseAccession。一个蛋白可同时参与多个调控环节；此类多角色注释按
identifier 表原意保留，在不同角色分面中分别显示，但在蛋白集合统计中仍按
BaseAccession 去重。

## Kla 相对信号热图

不同研究的原始 Kla 强度受富集方式、仪器、采集模式、搜索软件和定量层级影响，
不能直接比较。因此对每个定量样本分别进行组内排名：

1. 在同一样本内汇总同一规范化特征的信号；
2. 对原始信号计算 `log2(signal+1)`；
3. 在该样本全部可定量 Kla 特征中，以平均秩计算 0--100 百分位；
4. 提取调控蛋白的百分位；
5. 同一样本组有多个条件或重复时取百分位中位数；
6. 样本可定量但调控蛋白未检出时记为 0。

百分位定义为：

```text
百分位 = 100 × (平均秩 - 1) / (样本内特征数 - 1)
```

若样本仅有一个有效特征，则记为 100。该数值只表示调控蛋白在本样本 Kla
特征中的相对位置，不是跨研究表达量、Kla fold change 或差异分析结果。

## 普通全蛋白相对信号热图

普通全蛋白热图使用非 Kla 富集矩阵中的逐蛋白定量。反库、污染物和无效蛋白
记录被排除，同一 BaseAccession 的重复记录先合并。原始强度计算
`log2(signal+1)`；来源已经提供 log2 quantity 时直接排名，不再次 log2。
随后在每个普通全蛋白样本的全部有效蛋白中按上述平均秩公式计算百分位，再提取
调控蛋白。多个对照或技术重复取百分位中位数；普通全蛋白样本可用但目标蛋白
未检出时记为 0。Kla 富集强度从不用于替代普通全蛋白强度。

两套热图均按正常组织、癌症组织、正常细胞、癌症细胞排列，列按
Writer、Eraser、Writer-Eraser 和 Reader 分面。颜色由白色经黄色、橙色到
深红色，百分位越高颜色越暖。中文和英文版本使用相同数据和轴顺序。

## 四集合 Venn/Euler 分析

分别对全部 Kla 蛋白、Kla-DDR 蛋白、普通全蛋白和普通全蛋白-DDR 蛋白构建
四个类别集合。每个集合先按 BaseAccession 去重，再计算所有精确交集区域。
图形使用 `eulerr` 生成面积近似与集合大小成比例的四集合 Euler 表示；图中数字
来自精确 membership，椭圆几何仅为比例拟合。每套分析同时输出
`membership.csv`、`region_counts.csv` 和 `set_counts.csv`，因此所有图中数字
均可由表格独立复算。

## 软件、复现与质量控制

主要分析在 R 4.4.3 中完成，使用 dplyr 1.2.0、readr 2.2.0、tidyr 1.3.2、
ggplot2 4.0.2、eulerr 7.0.4、readxl 1.4.5 和 digest 0.6.39。异质检索结果
的规则提取和基础 pipeline 使用 Python 3.13.12。所有代码、配置、输入来源、
中间表、命令记录和测试均保存在项目目录。

自动测试验证最终 Kla 轴为 33 行、普通全蛋白轴为 30 行、四分类数量为
9/2/9/13、四个无严格参照组未进入最终分析、GeneSymbol 回退数为 0、样本内
百分位具有真实变化、共享参照未被重复显示，以及所有 Venn 数字均可由 membership
和 region 表精确复算。最终文件使用 SHA256 manifest 记录。

## 数据和代码可得性

公开质谱数据的 PXD、文章 DOI、下载链接和本地文件路径见
`kla_and_reference_teacher_review_zh.csv`。当前有效代码、重跑顺序和关键中间
文件见项目根目录 `NEW_CHAT_PROJECT_PROMPT.md`。最终双语图位于
`reanalysis/results/figures/`，表格位于 `reanalysis/results/tables/`。

## 方法写作参考

Hu Y, Xie M, Li Y, et al. Benchmarking clustering, alignment, and integration
methods for spatial transcriptomics. *Genome Biology*. 2024;25:212.
doi:10.1186/s13059-024-03361-0.
