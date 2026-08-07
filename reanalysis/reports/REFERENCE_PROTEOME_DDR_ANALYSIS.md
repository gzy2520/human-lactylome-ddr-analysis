# 常规蛋白组与 Kla 蛋白组的 DDR 占比分析

分析日期：2026-08-05

## 1. 分析问题

对现有 10 类组织或细胞模型，分别使用已经选定的常规全蛋白组作为检出背景，
计算 DNA repair/damage response（DDR）相关蛋白占全部检出蛋白的比例，并与
同一模型的 Kla 蛋白集合中 DDR 蛋白比例并排比较。

计算公式：

`DDR占比 = GO repair/damage 相交蛋白数 / 对应蛋白集合总数`

GO 注释使用 `data/annotations/GO-repair+damage(human).tsv`。仅保留人源注释，
排除 qualifier 中含 `NOT` 的记录。匹配时优先使用去除 isoform 后缀的
UniProt `BaseAccession`，来源数据具有 GeneSymbol 时才使用基因符号辅助匹配。

## 2. 每个常规蛋白组如何选择

| 模型 | 对照数据 | 年份 | 选用内容 | 选择理由 | 主要限制 |
|---|---|---:|---|---|---|
| Human hippocampus | PXD043880 | 2023 | 74 名神经学正常 CA1 供者 | 组织精确匹配、年龄接近 Kla 海马样本、无 PTM 富集且有发布矩阵 | 为老年 CA1，不代表年轻或完整海马 |
| HEK293T | PXD030304 | 2022 | 401 个 Control_HEK293T_lys；排除 standard-QC | 细胞系精确匹配、未处理裂解物、标准化 DIA 平台且数据量大 | 属于过程控制运行，不能作为 401 个生物重复 |
| HK-2 | PXD072220 | 2026 | Control_1、Control_3、Control_4 | 细胞系精确匹配、明确未处理、3 个对照且完整度高 | 只能使用 control，不能混入 Cu 或 TTM |
| MCF10A | PXD002400 | 2015 | 未处理 MCF-10A，10 个 IEF 组分×2 次进样 | 细胞系精确匹配、无 PTM 富集、分级深度较高且有 MaxQuant evidence | 平台较旧；分级和进样是技术深度，不是生物重复 |
| MCF7 | PXD030304 | 2022 | MCF7 / SIDM00148 | 精确未处理细胞系，统一 DIA 图谱和高置信矩阵 | 只作检出背景，不与 Kla 富集数据直接比较绝对强度 |
| HCT116 | PXD030304 | 2022 | HCT-116 / SIDM00783 | 精确未处理细胞系，可作为独立全细胞背景 | Kla 数据含亚细胞组分，不能据此推断亚细胞缺失 |
| T-ALL | PXD030304 | 2022 | TALL-1 主替代；Jurkat 敏感性替代 | 未找到独立未处理 TALL-104 全蛋白组，只能使用疾病类别相近模型 | 不是 TALL-104 精确细胞系，结果等级为 C |
| MDA-MB-468 | PXD030304 | 2022 | MDA-MB-468 / SIDM00628 | 精确未处理细胞系，统一 DIA 图谱 | 只作蛋白检出背景 |
| T-47D | PXD030304 | 2022 | T47D / SIDM00097 | 精确未处理细胞系，名称已映射到项目标签 T-47D | 只作蛋白检出背景 |
| RKO | PXD030304 | 2022 | RKO / SIDM01090 | 精确未处理细胞系，可作为 WT/KO Kla 实验的共同背景 | 不能替代 WT 与 KO 的组内比较 |

完整字段还包括样本数、蛋白数、仪器、检索软件、PTM 富集状态、实际分析文件、
仓库完整度、数据集链接和论文链接，保存在
`cell_type_reference_control_information.csv` 及分析工作簿的“对照选择信息”页。

## 3. 结果

| 组织或细胞系 | 常规蛋白数 | 常规DDR蛋白数 | 常规DDR占比 | Kla蛋白数 | Kla∩DDR | Kla DDR占比 | Kla/常规倍数 |
|---|---:|---:|---:|---:|---:|---:|---:|
| Human hippocampus | 2,092 | 77 | 3.68% | 853 | 30 | 3.52% | 0.96 |
| HEK293T | 6,441 | 381 | 5.92% | 438 | 44 | 10.05% | 1.70 |
| HK-2 | 4,897 | 308 | 6.29% | 1,250 | 134 | 10.72% | 1.70 |
| MCF10A | 4,839 | 312 | 6.45% | 947 | 116 | 12.25% | 1.90 |
| MCF7 | 3,099 | 209 | 6.74% | 1,316 | 139 | 10.56% | 1.57 |
| HCT116 | 4,010 | 273 | 6.81% | 643 | 74 | 11.51% | 1.69 |
| T-ALL | 3,383 | 267 | 7.89% | 368 | 37 | 10.05% | 1.27 |
| MDA-MB-468 | 3,760 | 230 | 6.12% | 1,612 | 199 | 12.34% | 2.02 |
| T-47D | 3,516 | 228 | 6.48% | 1,301 | 161 | 12.38% | 1.91 |
| RKO | 3,615 | 252 | 6.97% | 1,217 | 116 | 9.53% | 1.37 |

## 4. 主要观察

1. 10 类模型中有 9 类的 Kla DDR 占比高于常规蛋白组背景。
2. Human hippocampus 基本相当：常规背景为 3.68%，Kla 为 3.52%，相差
   -0.16 个百分点。
3. 描述性倍数最高的是 MDA-MB-468：Kla 为 12.34%，常规背景为 6.12%，
   约为 2.02 倍。
4. T-47D 和 MCF10A 也接近 1.9 倍；HEK293T、HK-2 和 HCT116 约为 1.7 倍。
5. T-ALL 的常规背景 DDR 占比在 10 类中最高，为 7.89%，因此其 Kla/背景
   倍数相对较低，为 1.27。

这些结果支持“多数细胞系的 Kla 蛋白集合中 DDR 相关蛋白所占比例高于常规
可检出蛋白背景”这一描述性结论，但不能单凭该比例证明 Kla 对 DDR 蛋白存在
因果性或统计学显著富集。

## 5. T-ALL 替代模型敏感性

原 Kla 模型是 TALL-104，未找到精确匹配的独立未处理外部全蛋白组。
主表使用 TALL-1；Jurkat 仅作替代敏感性分析：

| 替代集合 | 蛋白数 | DDR蛋白数 | DDR占比 |
|---|---:|---:|---:|
| TALL-1 primary | 3,383 | 267 | 7.89% |
| Jurkat sensitivity | 3,363 | 245 | 7.29% |
| TALL-1/Jurkat 并集 | 3,835 | 285 | 7.43% |
| TALL-1/Jurkat 交集 | 2,911 | 227 | 7.80% |

不同替代集合的 DDR 背景占比为 7.29%-7.89%，方向稳定，均低于 TALL-104
Kla 集合的 10.05%。但这些结果仍不能替代精确 TALL-104 外部对照。

## 6. 分母和匹配注意事项

- Human hippocampus 的 2,092 个分母是论文发布的 gene/protein feature，
  通过 GeneSymbol 与 GO-DDR 匹配；77 个 DDR feature 均为基因符号辅助匹配。
- HEK293T 使用 401 个 `Control_HEK293T_lys` 运行中检出的 6,692 高置信
  蛋白范围内的 BaseAccession 并集；standard-QC 运行排除。
- HK-2 使用三个未处理 control 的蛋白并集。
- MCF10A 因仓库未提供 `proteinGroups.txt`，分母为 20 个 MCF-10A
  分级/进样 evidence 中的唯一 leading-razor BaseAccession。
- PXD030304 的癌细胞系使用 6,692 高置信 averaged matrix 中对应细胞系行的
  非缺失 BaseAccession。

## 7. 统计解释限制

这些数据来自不同研究、仪器、数据库检索流程和采集设计。多数参考数据中的
运行数是技术采集、过程控制或肽段分级，不是独立生物学重复。因此本次只进行
蛋白集合层面的描述性比例比较，不进行跨项目丰度显著性检验。

常规蛋白组未检出某蛋白，不能证明该蛋白在样本中不存在；Kla 数据未检出某
蛋白，也不能证明该蛋白不发生乳酸化。

## 8. 输出文件

- `reanalysis/results/tables/cell_type_reference_control_information.csv`
- `reanalysis/results/tables/cell_type_reference_ddr_statistics.csv`
- `reanalysis/results/tables/cell_type_kla_vs_reference_ddr_statistics.csv`
- `reanalysis/results/tables/cell_type_kla_vs_reference_ddr_statistics.xlsx`
- `reanalysis/results/tables/reference_proteome_all_proteins.csv`
- `reanalysis/results/tables/reference_proteome_ddr_proteins.csv`
- `reanalysis/results/tables/reference_proteome_go_unmatched.csv`
- `reanalysis/results/tables/tall104_surrogate_ddr_sensitivity.csv`
- `reanalysis/results/figures/cell_type_kla_vs_reference_ddr_fraction.png`
- `reanalysis/results/figures/cell_type_kla_vs_reference_ddr_fraction.pdf`
- `reanalysis/results/figures/cell_type_kla_vs_reference_ddr_fraction_v2.png`
- `reanalysis/results/figures/cell_type_kla_vs_reference_ddr_fraction_v2.pdf`

新版图使用中文标题和图例，并在每根柱末端同时标注
`DDR蛋白数/总蛋白数（DDR占比）`。图中仅纳入同时具有完整 Kla 蛋白集合、
DDR 交集和可计数常规蛋白组的 10 类模型；其余新增样本不会在分母口径尚未
统一时提前加入。
