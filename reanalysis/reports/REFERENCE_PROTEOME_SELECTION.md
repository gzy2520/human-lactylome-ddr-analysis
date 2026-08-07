# 常规全蛋白组对照数据选型说明

核验日期：2026-07-30

## 1. 纳入范围

本次不再重新推断细胞类型，严格使用现有统计表
`reanalysis/results/tables/cell_type_kla_ddr_statistics.csv` 中的 10 类：

1. Human hippocampus
2. HEK293T
3. HK-2
4. MCF10A
5. MCF7
6. HCT116
7. T-ALL
8. MDA-MB-468
9. T-47D
10. RKO

选型表中的 Kla 蛋白数、Kla∩DDR 蛋白数和比例均与该统计表逐行核对一致。

## 2. 选型目标

每类优先选择：

- 相同细胞系或相同组织；
- 未处理、野生型或明确的 control/baseline 子集；
- 不进行 Kla 或其他 PTM 富集、免疫沉淀、BioID、膜蛋白富集；
- 蛋白鉴定深度较高；
- 原始文件、处理矩阵、样本映射和方法信息尽量完整。

普通肽段分级可用于增加蛋白组深度，不视为 PTM 偏好富集。

癌细胞系在本报告中的“正常/基线”仅表示该细胞系的 untreated baseline，
不表示其为生物学意义上的正常细胞。

## 3. 最终推荐

| 当前模型 | 推荐对照 | 使用子集 | 对照蛋白数 | 匹配等级 |
|---|---|---|---:|---|
| Human hippocampus | PXD043880 | 74 名神经学正常的 CA1 海马供者 | 2,092 | A，精确组织匹配 |
| HEK293T | PXD030304 | Control_HEK293T_lys；排除 standard-QC | 6,441 | B，精确细胞系匹配 |
| HK-2 | PXD072220 | Control_1、Control_3、Control_4 | 4,897 | A，精确细胞系匹配 |
| MCF10A | PXD002400 | 未处理 MCF-10A；20 个分级/技术进样 raw | 4,839 | B，精确细胞系匹配 |
| MCF7 | PXD030304 | MCF7 / SIDM00148 | 3,099 | A，精确细胞系匹配 |
| HCT116 | PXD030304 | HCT-116 / SIDM00783 | 4,010 | A，精确细胞系匹配 |
| T-ALL | PXD030304 | TALL-1 / SIDM00370；Jurkat / SIDM01016 作敏感性替代 | 3,383（TALL-1 主分析） | C，疾病类别替代 |
| MDA-MB-468 | PXD030304 | MDA-MB-468 / SIDM00628 | 3,760 | A，精确细胞系匹配 |
| T-47D | PXD030304 | T47D / SIDM00097 | 3,516 | A，精确细胞系匹配 |
| RKO | PXD030304 | RKO / SIDM01090 | 3,615 | A，精确细胞系匹配 |

4 个主要参考 PXD 覆盖全部 10 类：

- PXD043880：Human hippocampus
- PXD072220：HK-2
- PXD002400：MCF10A
- PXD030304：HEK293T、MCF7、HCT116、T-ALL 替代、MDA-MB-468、T-47D、RKO

其中 9 类为精确细胞系/组织匹配，T-ALL 为替代匹配。

## 4. 关键限制

### T-ALL

当前 Kla 项目标签 `T-ALL` 对应的原始模型实际是 `TALL-104`。
未找到独立、未处理、精确匹配 TALL-104 的外部常规全蛋白组。

主分析应将 TALL-1 和 Jurkat 分开输出。只有二者均检出的蛋白，
才可标记为“由两个 T-ALL 替代背景共同支持”，不能写成 TALL-104 精确背景。

PXD028488 中的 TALL-104 non-enrichment 数据仍经历乳酸暴露，
且来自同一 Kla 研究，只能作为同细胞系技术敏感性检查，不能作为正常主对照。

### 技术采集不等于生物学重复

- PXD030304 的每个癌细胞系 6 次 DIA 采集不能自动视为 6 个独立生物学重复。
- HEK293T 的 401 个 lysate run 是长期过程控制采集，不能作为 401 个生物学重复。
- 663 个 HEK293T standard-QC run 不进入主要蛋白集合。
- PXD002400 的 IEF 组分和重复进样用于增加检出深度，不能用于生物重复统计。

### 比较解释

常规全蛋白组只作为“在相应模型中可检出的蛋白背景”。

- reference 中检出、Kla 中未检出：只能表述为“未在当前 Kla 实验中鉴定”。
- Kla 中检出、reference 中未检出：可能来自低丰度、批次、仪器、搜索库差异或富集增敏，不能直接删除。
- 常规全蛋白组没有 Kla 位点证据，不能用来证明某个蛋白或位点“不乳酸化”。
- 不同平台之间不直接比较绝对强度。

## 5. 后续分析规则

1. 统一用去除 isoform 后缀后的 `BaseAccession` 匹配，`GeneSymbol` 仅作辅助。
2. PXD030304 主分析使用 6,692 蛋白高置信矩阵。
3. PXD030304 同时使用 8,498 蛋白矩阵生成敏感性结果。
4. 每种模型分别输出：
   - `reference_detected`
   - `Kla_intersection_reference`
   - `Kla_not_in_reference`
   - DDR 相关归一化比例
   - accession/gene 映射失败记录
5. TALL-1 与 Jurkat 必须分别输出，再生成共同检出集合。

## 6. 数据来源

- PXD030304：https://proteomecentral.proteomexchange.org/?pxid=PXD030304
- PXD030304 处理矩阵：https://doi.org/10.6084/m9.figshare.19345397
- PXD030304 论文：https://doi.org/10.1016/j.ccell.2022.06.010
- PXD043880：https://proteomecentral.proteomexchange.org/?pxid=PXD043880
- PXD043880 论文：https://doi.org/10.1186/s13024-023-00650-3
- PXD072220：https://proteomecentral.proteomexchange.org/?pxid=PXD072220
- PXD072220 论文：https://doi.org/10.1152/ajpcell.00311.2026
- PXD002400：https://proteomecentral.proteomexchange.org/?pxid=PXD002400

## 7. 当前完成状态

已经完成 10 类模型的参考数据选型、必要处理结果下载、样本子集核验、
蛋白集合构建、GO-DDR 交集和 Kla-reference DDR 占比比较。

本次没有下载全部大型 raw 文件；分析优先使用作者发布的处理矩阵、检索结果和
补充表。正式分析结果见 `reanalysis/reports/REFERENCE_PROTEOME_DDR_ANALYSIS.md`。

正式输出：

- `reanalysis/results/tables/cell_type_reference_proteome_selection.xlsx`
- `reanalysis/results/tables/cell_type_reference_proteome_selection.csv`
- `reanalysis/results/tables/reference_proteome_pxd030304_sample_audit.csv`
- `reanalysis/reports/REFERENCE_PROTEOME_SELECTION.md`
