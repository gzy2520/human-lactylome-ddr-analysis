# 各数据集的数据来源、Kla 提取算法与分析结果说明

## 1. 本文件解决什么问题

本文件专门说明以下四件事：

1. 每个 PXD 数据从哪里获得，本地保存了什么。
2. 程序实际读取了哪个文件，没有读取哪些文件。
3. 程序根据什么字段判断一个位点是赖氨酸乳酸化（Kla）。
4. Kla 位点怎样转换成蛋白集合、怎样与 GO repair/damage 取交集，以及怎样生成 Venn 图。

核心原则是：**原始质谱文件只作为原始数据保存，不直接当作 Kla 证据。** 本次分析优先使用作者已经完成数据库检索后导出的位点表、肽段表、蛋白表或论文补充表。这样可以避免在没有统一搜索参数、谱图库和 FDR 设置的情况下，自行重新搜索 raw 数据而引入新的不可比性。

## 2. 数据纳入范围

### 2.1 正式纳入

- PXD014870
- PXD028488
- PXD050470
- PXD053474
- PXD060185
- PXD078013
- PXD078736

### 2.2 暂不分析

- PXD038880
- PXD050906

PXD038880/PXD050906 的文章、下载文件和元数据仍保存在 `data` 中，但按照老师要求，不进入 Kla 蛋白集合、GO 交集或 Venn 图。程序还包含自动检查：如果这两个 PXD 出现在主分析长表中，pipeline 会直接报错停止。

## 3. 一句话概括整个算法

```text
作者检索结果/补充表
    -> 识别明确的 Kla 位点
    -> 排除反库、污染物、非人源和无法定位的记录
    -> 保留样本、实验组、重复和来源行号
    -> 统一 UniProt accession 并去除 isoform 后缀
    -> 在最后一步按 BaseAccession 生成蛋白集合
    -> 与 human GO repair/damage 取交集
    -> 按海马体、永生化模型、肿瘤细胞系生成精确 Venn 区域
```

这不是机器学习分类算法，也不是根据蛋白名称猜测 Kla。它是一个**基于作者检索证据、规则明确且可追溯的集合分析流程**。

## 4. “阈值改成 0”到底是什么意思

老师所说的阈值 0，在当前流程中指的是：

> 对作者已经报告为 Kla 的位点，主分析不再额外施加 `LocalizationProb >= 0.75` 之类的定位概率截断，即额外定位概率阈值为 0。

这不等于取消所有质量控制。以下规则仍然保留：

- reverse/decoy 记录排除。
- contaminant 记录排除。
- 非人源目录排除。
- 没有明确 Kla 修饰或无法映射到蛋白位置的记录排除。
- PEAKS PSM/DDA 仍要求 `-10lgP >= 20`，因为这是肽段鉴定可信度控制，不是定位概率阈值。
- PXD078013 必须同时满足 evidence 的 Kla 阳性和 proteinGroups 的位点映射。

PXD014870 另外保存了 `LocalizationProb >= 0.75` 的敏感性分析，但它不进入主结果：

| PXD014870 结果 | Kla 位点 | Kla 蛋白 |
|---|---:|---:|
| 主分析，额外定位概率阈值 0 | 442 | 193 |
| 敏感性分析，定位概率阈值 0.75 | 407 | 185 |

因此，阈值 0 的作用是避免在数据量已经较少时再次删除作者报告位点，同时仍保留基本的证据有效性检查。

## 5. 每个 PXD 是怎样处理的

### 5.1 PXD014870

**文章和数据来源**

- 相关方法学文章：[Cyclic immonium ion of lactyllysine reveals widespread lactylation in the human proteome](https://www.nature.com/articles/s41592-022-01523-1)
- DOI：`10.1038/s41592-022-01523-1`
- PXD014870 是较早的 MCF7 数据来源，原始实验发表于 2019 年，2022 年文章又使用该数据进行分析。
- 本地目录：`data/PXD014870/search_results/`

**细胞和处理**

- 细胞：MCF7。
- 条件：DCA、Hypoxia、Rotenone、U-13C6 glucose。
- 数据库检索软件：MaxQuant。

**程序读取的文件**

- 每个条件的 `txt/Lactyl (K)Sites.txt`
- 每个条件的 `txt/modificationSpecificPeptides.txt`

**Kla 判断规则**

1. 从 `Lactyl (K)Sites.txt` 读取作者定位的赖氨酸乳酸化位点。
2. 要求该位点的 `id` 同时存在于 `modificationSpecificPeptides.txt` 的有效 `Lactyl (K) site IDs` 中。
3. `Amino acid` 必须为 K。
4. 必须存在有效蛋白 accession 和蛋白位置。
5. 排除 reverse 和 potential contaminant。
6. 主分析不要求定位概率达到 0.75；定位概率只保留为证据字段，并另做 0.75 敏感性分析。

**得到的结果**

- sample-level 证据行：675。
- 去重 accession-site：442。
- 去重 Kla 蛋白：193。
- Kla 与 GO repair/damage 交集蛋白：16。

**为什么一个文章会对应两个 PXD**

PXD014870 和 PXD028488 都与 2022 Nature Methods 文章有关，但不是重复下载。PXD014870 是早期 MCF7 数据的再利用；PXD028488 是该方法学文章新沉积的 PEAKS、碰撞能和 diagnostic ion 实验数据。两者的实验设计和检索结果文件不同，所以应分别读取，最后再按蛋白集合合并。

### 5.2 PXD028488

**文章和数据来源**

- 文章：[Cyclic immonium ion of lactyllysine reveals widespread lactylation in the human proteome](https://www.nature.com/articles/s41592-022-01523-1)
- DOI：`10.1038/s41592-022-01523-1`
- 本地目录：`data/PXD028488/search_results/`
- 检索结果包含 enrichment、non-enrichment 和不同 HCD 碰撞能目录。

**细胞和目录盘点**

- 人源：HCT116、T-ALL、HEK293T。
- 非人源：BV2、RAW264.7，均排除。
- 共盘点 19 个搜索目录：13 个进入主分析，5 个排除，1 个 all-HCD 聚合目录只用于验证。
- HCT116 non-enrichment 因缺少 `DB search psm.csv`，无法按同一规则提取，因此排除并记录原因。

**程序读取的文件**

- `DB search psm.csv`
- `protein-peptides.csv`
- `proteins.csv`
- `PSM ions.csv`，只用于检查 marker 156 辅助证据。

**Kla 判断规则**

1. PSM 的 `-10lgP` 必须大于等于 20。
2. 乳酸化必须由下列信息之一明确给出：
   - `AScore` 中出现 K 位点的 Lac/Lactylation；
   - 修饰肽段中出现 `K(+72.02)`；
   - 修饰肽段中出现 `K(Lactyl)` 或等价 Lac 标记；
   - `PTM` 字段标注 Lactylation，同时还必须能确定具体 K 位置。
3. 从 `protein-peptides.csv` 读取肽段在蛋白中的 `Start`，计算蛋白位点：

```text
蛋白 Kla 位置 = 肽段起始位置 + 肽段内 Kla 的 K 位置 - 1
```

4. 肽段映射的蛋白 accession 必须与 PSM accession 一致；共享肽段不能无条件映射到无关蛋白。
5. reverse、decoy、contaminant 排除。
6. `PSM ions.csv` 中理论 m/z 约为 156.1025 的 marker(156) 只提高来源置信度，不能单独把普通肽段判定为 Kla。原因是 marker 156 也可能出现在非乳酸化肽段记录中。
7. `HEK293T-Enrichment-all HCD` 是多个 HCD 条件的聚合结果，只作为验证，避免与 HCD24/28/32/36/40 重复计数。

**得到的结果**

- 主分析 sample-level 证据行：4,888。
- 去重 accession-site：1,278。
- 去重 Kla 蛋白：665。
- Kla 与 GO repair/damage 交集蛋白：65。

**旧流程遗漏**

旧流程只使用三个目录，而且包含 all-HCD 聚合目录，未系统读取其余可分析的人源 HCD 和 non-enrichment 目录。新流程盘点全部 19 个目录后，PXD028488 比旧流程新增 638 个蛋白。这个变化来自数据覆盖范围扩大，不是降低定位概率阈值得到的。

### 5.3 PXD050470

**文章和数据来源**

- 文章：[Global Profiling of Protein Lactylation in Human Hippocampi](https://onlinelibrary.wiley.com/doi/10.1002/prca.202400061)
- DOI：`10.1002/prca.202400061`
- 组织：人海马体生理组织。
- 本地补充表目录：`data/PXD050470/supplementary/`

**本地数据情况**

仓库记录提到 27 个 raw 文件，但本地没有这些 raw。根据“先使用作者已经完成的检索结果、不自行重新搜索 raw”的原则，本数据集使用论文补充表：

- `prca2331-sup-0005-tables3.xlsx`，Supplementary Table S3。
- `prca2331-sup-0014-tables12.xlsx`，Supplementary Table S12。

**Kla 判断规则**

1. S3/S12 本身就是作者整理的乳酸化位点表，因此 `Protein accession + Positions within proteins` 作为 Kla 位点证据。
2. 保留作者的 localization probability、PEP、Score、modified sequence 和 H072/H081/H187 样本信息。
3. S3 的 Class I 信息保留为高置信作者证据。
4. S3 和 S12 中重复出现的相同 accession-site 在蛋白集合阶段只计一次，但来源文件和样本证据都保留。

**曾经出现的格式问题**

旧读取代码将 Excel 表头位置偏移了一行，导致有效数据被读成 0 行。根据工作簿实际结构修正表头后，能够正常得到 2,600 个 accession-site。因此，PXD050470 原始补充表是完整的，问题出在旧代码对 Excel 格式的假设，不是文章没有提供可处理的数据。

**得到的结果**

- sample-level 证据行：11,537。
- 去重 accession-site：2,600。
- 去重 Kla 蛋白：853。
- Kla 与 GO repair/damage 交集蛋白：30。

### 5.4 PXD053474

**文章和数据来源**

- 文章：[Subcellular Proteomic Mapping of Lysine Lactylation](https://pubs.acs.org/doi/10.1021/jasms.4c00366)
- DOI：`10.1021/jasms.4c00366`
- 细胞：HCT116。
- 本地目录：`data/PXD053474/search_results/` 和 `data/PXD053474/supplementary/`

**程序读取的四套检索结果**

- Enriched DDA。
- Unenriched DDA。
- Enriched DIA。
- Unenriched DIA。
- 另外读取作者 Supplementary Table S3：`js4c00366_si_003.xlsx`。

**DDA Kla 判断规则**

1. 读取 PEAKS `peptide.csv`、`protein-peptides.csv` 和 `proteins.csv`。
2. `-10lgP >= 20`。
3. AScore、PTM 或修饰肽段必须明确指向 Lac/Lactylation，并能定位到 K。
4. 通过 `protein-peptides.csv` 的肽段起始位置计算蛋白 Kla 位点。
5. 要求肽段映射 accession 与鉴定 accession 一致。

**DIA Kla 判断规则**

1. 读取 `ptm-site_Report.tsv`。
2. `PTM.ModificationTitle` 必须等于 `Lac`。
3. `PTM.SiteAA` 必须等于 K。
4. 物种必须为 Homo sapiens。
5. 读取 `PTM.ProteinId` 和 `PTM.SiteLocation` 作为蛋白和位点。
6. 保留每个 run 的 `PG.IsIdentified`、site probability、PEP 和 Cscore。

**检索结果和 S3 怎样整合**

以 `BaseAccession + KlaSite` 作为比较键：

| 比较类别 | 位点数 |
|---|---:|
| search results 与 S3 都存在 | 1,275 |
| 仅 search results 存在 | 826 |
| 仅 S3 存在 | 3 |

826 个 search-only 位点中，只有 20 个同时获得 DDA 和 DIA 两种独立采集模式支持；其余 806 个只得到单一模式支持。

主结果采用以下可追溯规则：

```text
PXD053474 主位点 = 作者 S3 的全部位点
                  + 不在 S3 中但同时被 DDA 和 DIA 支持的位点
```

所以主结果为 `1,278 + 20 = 1,298` 个 accession-site。806 个单模式 search-only 候选保留在审计表中，不进入主蛋白集合，也没有被删除。

**得到的结果**

- 主分析 sample-level 证据行：9,655。
- 主分析 accession-site：1,298。
- 去重 Kla 蛋白：590。
- Kla 与 GO repair/damage 交集蛋白：67。

### 5.5 PXD060185

**文章和数据来源**

- DOI：`10.1186/s13046-025-03512-6`
- 本地检索结果：`data/PXD060185/search_results/RESULT/combined/txt/`
- 读取主文件：`La (K)Sites.txt`

**细胞和样本映射**

- A -> MCF7。
- B -> MDA-MB-468。
- C -> MCF10A。
- D -> T-47D。

**Kla 判断规则**

1. `La (K)Sites.txt` 是 MaxQuant 已定位的乳酸化位点表。
2. `Amino acid` 必须为 K。
3. 排除 reverse 和 potential contaminant。
4. 必须有蛋白 accession 和蛋白位置。
5. 一个 Kla 位点是否属于某个样本，依据该样本的 `Identification type A/B/C/D`，而不是仅凭 intensity。原因是强度列可能经过缺失值填补，不能把“有强度”简单解释成“该样本鉴定到该位点”。
6. 每个样本分别读取自己的 localization probability、PEP 和 Score。
7. 样本在读取阶段不合并；最后求细胞类别蛋白集合时才取并集。

**得到的结果**

- sample-level 证据行：12,113。
- 去重 accession-site：5,583。
- 去重 Kla 蛋白：1,991。
- Kla 与 GO repair/damage 交集蛋白：223。

### 5.6 PXD078013

**数据来源**

- 细胞：RKO 结直肠癌细胞。
- 处理：GSK3B WT 与 KO，各两个重复。
- 本地目录：`data/PXD078013/search_results/`
- 仓库没有独立的 `La(K)Sites.txt`。

**程序读取的文件**

- `evidence.txt`
- `proteinGroups.txt`
- `peptides.txt` 保留为辅助检索结果，但主位点判断由前两个文件联合完成。

**Kla 判断规则**

1. evidence 行的 `La (K)` 必须大于 0。
2. evidence 行必须存在 `La (K) site IDs`。
3. 排除 reverse 和 potential contaminant。
4. 从 `proteinGroups.txt` 获取相同 site ID 对应的 `La (K) site positions`。
5. site ID 数量必须与 position 数量相同，然后按顺序一一配对；数量不一致时不做笛卡尔组合，避免制造不存在的位点。
6. proteinGroups 中还要存在有效 majority protein ID 或 protein ID。
7. 保留 sg_1、sg_2、sg737_1、sg737_2，映射为 WT/KO 和具体重复。

**为什么这种联合判断是合理的**

`evidence.txt` 能证明某条 PSM/肽段包含 La(K)，但单独使用 evidence 不一定能稳定得到蛋白绝对位置；`proteinGroups.txt` 提供蛋白和 site position。只有两者通过 site ID 对上，才生成最终 Kla 位点，因此不是从 proteinGroups 的普通蛋白存在信息直接推断 Kla。

**得到的结果**

- sample-level 证据行：13,836。
- 去重 accession-site：3,655。
- 去重 Kla 蛋白：1,217。
- Kla 与 GO repair/damage 交集蛋白：115。

### 5.7 PXD078736

**数据来源**

- 细胞：HK-2 永生化肾小管上皮细胞。
- 处理：ctrl 与 man，各三个重复。
- 本地目录：`data/PXD078736/search_results/txt/`
- 读取主文件：`La(K)Sites.txt`

**Kla 判断规则**

1. 读取 MaxQuant `La(K)Sites.txt` 中已经定位的 La(K) 位点。
2. `Amino acid` 必须为 K。
3. 排除 reverse 和 potential contaminant。
4. 必须有蛋白 accession 和蛋白位置。
5. 依据 `Identification type ctr_1` 至 `Identification type man_3` 判断位点在哪个样本被鉴定。
6. 不使用所有样本都有值的 intensity 列来强行复制位点。
7. ctrl/man 和六个重复在 sample-level 表中分别保留，最后才按 HK-2 所属类别求蛋白并集。

**得到的结果**

- sample-level 证据行：11,142。
- 去重 accession-site：3,661。
- 去重 Kla 蛋白：1,250。
- Kla 与 GO repair/damage 交集蛋白：133。

### 5.8 PXD038880 和 PXD050906

这两个 PXD 与 AGS 亲本/顺铂耐药和 NBS1 lactylation 研究有关。当前只保留：

- raw 下载文件。
- 检索结果。
- 文章和元数据。
- excluded 说明。

由于老师认为数据可能有问题，当前 pipeline 不调用这两个数据集的提取器。它们不会影响 3,112 个 Kla 蛋白、275 个 Kla-GO 蛋白或任何 Venn 区域。

## 6. sample-level 长表保存了什么

每条证据在读取阶段都保留以下信息：

- PXD 和 DOI。
- SampleName、CellOrTissueType、ExperimentalGroup、Replicate。
- enrichment 状态和 DDA/DIA 模式。
- 原始 Accession 和去除 isoform 后缀的 BaseAccession。
- GeneSymbol、ProteinName。
- KlaSite、ModifiedPeptide。
- LocalizationProb、PEP、Score。
- diagnostic ion、Class I 等辅助信息。
- SourceFile、SourceRow、SiteID。
- EvidenceMode、SourceConfidence。
- 是否进入主结果及未进入的原因。

所以最终蛋白表中的每个蛋白都能追溯到具体 PXD、样本、来源文件和来源行。样本、实验组和重复不会在读取阶段被提前合并。

## 7. UniProt 统一和去重算法

不同文件可能把同一个蛋白写成：

```text
sp|P49959|MRE11_HUMAN
P49959
P49959-2
```

程序先去除 `sp|`、`tr|`、reverse/contaminant 前缀等数据库格式，再去除 isoform 后缀 `-数字`，得到统一的：

```text
BaseAccession = P49959
```

去重分两个层次：

1. 位点层：`BaseAccession + KlaSite`，例如 `P49959 + K510`。
2. 蛋白层：只按 `BaseAccession` 计数，同一个蛋白在多个样本、多个重复或 DDA/DIA 中出现仍只算一个蛋白。

样本证据没有丢失，而是以分号合并保存在蛋白汇总表的 Sample、PXD、KlaSites、SourceFile 和 EvidenceMode 字段中。

## 8. GO repair/damage 交集算法

GO 输入文件为：

`data/annotations/GO-repair+damage(human).tsv`

处理步骤：

1. 只保留人源 taxon 9606。
2. 排除 qualifier 中含 `NOT` 的否定注释。
3. 主匹配键为 UniProt BaseAccession。
4. BaseAccession 未匹配时，才使用 GeneSymbol 做辅助匹配。
5. GeneSymbol fallback 也只使用人源且非 NOT 的 GO 注释。
6. 主分析不按 GO evidence code 再次缩减，但保留 IDA、IMP、IEA 等 evidence code，后续可生成高置信 GO 版本。

这里的 GO 交集不会决定一个蛋白是不是 Kla。正确顺序是：

```text
先由质谱/补充表确定 Kla 蛋白
再判断该 Kla 蛋白是否属于 GO repair/damage
```

因此，未进入 GO 交集不等于 Kla 证据被删除；它只表示该蛋白没有匹配到当前 GO repair/damage 表。

## 9. 细胞分类和 Venn 算法

最终采用的分类为：

| 类别 | 成员 |
|---|---|
| hippocampus_tissue | Human hippocampus |
| normal_immortalized_cell_lines | HEK293T、HK-2、MCF10A |
| tumor_cell_lines | MCF7、HCT116、T-ALL、MDA-MB-468、T-47D、RKO |

注意：HEK293T 和 HK-2 是转化/永生化模型，不是正常组织；MCF10A 是本项目确认的第三个非肿瘤/永生化模型；MCF7 归入肿瘤细胞系。

每个类别先对其所有样本的 BaseAccession 取并集。设海马体集合为 H、永生化模型为 N、肿瘤细胞为 T，则七个精确区域为：

```text
hippocampus_only              = H - N - T
normal_only                   = N - H - T
tumor_only                    = T - H - N
hippocampus_and_normal_only   = (H ∩ N) - T
hippocampus_and_tumor_only    = (H ∩ T) - N
normal_and_tumor_only         = (N ∩ T) - H
all_three                     = H ∩ N ∩ T
```

分别对两类数据做 Venn：

1. 全部 Kla 蛋白。
2. Kla 与 GO repair/damage 的交集蛋白。

当前精确区域计数为：

| Venn 区域 | 全部 Kla | Kla 与 GO repair/damage |
|---|---:|---:|
| hippocampus only | 397 | 4 |
| immortalized only | 231 | 15 |
| tumor only | 923 | 99 |
| hippocampus + immortalized only | 32 | 0 |
| hippocampus + tumor only | 157 | 6 |
| immortalized + tumor only | 1,105 | 131 |
| all three | 267 | 20 |
| 七区域合计 | 3,112 | 275 |

`tumor_specific_kla_proteins.csv` 对应严格的 `tumor_only`，不是“所有在肿瘤中出现过的蛋白”；`tumor_specific_kla_ddr_proteins.csv` 同理。因此肿瘤特异结果分别为 923 和 99 个蛋白。

## 10. MRE11、XLF 和 NBS1 为什么在旧数据中没有

目标蛋白统一标识为：

- MRE11：GeneSymbol `MRE11`，UniProt `P49959`。
- XLF：正式基因名 `NHEJ1`，UniProt `Q9H9Q4`。
- NBS1：正式基因名 `NBN`，UniProt `O60934`。

源表级审计结果为：

| 数据集 | MRE11/XLF/NBS1 的情况 | 结论 |
|---|---|---|
| PXD014870 | 可用作者检索/Kla 表中未找到目标记录 | 不是阈值或 GO 删除 |
| PXD028488 | 普通蛋白或未乳酸化肽段中能找到目标蛋白，但目标 Kla 记录为 0 | 蛋白被鉴定不等于检测到 Kla |
| PXD050470 | S3/S12 中未找到目标 Kla 记录 | 作者补充表本身没有这些 Kla 蛋白 |
| PXD053474 | 普通蛋白、普通肽段或非 Lac PTM 中能找到目标蛋白，但 PTM=Lac 的目标记录为 0 | 不是 DDA/DIA 整合规则删除 |
| PXD060185 | 检出 MRE11、NHEJ1、NBN 的 Kla | 新数据中确实存在 Kla 证据 |
| PXD078736 | 检出 MRE11 Kla | 新数据中确实存在 Kla 证据 |

PXD060185 中的目标证据包括：

- MRE11：K510、K609、K625、K673。
- XLF/NHEJ1：K288。
- NBS1/NBN：K435、K441。

PXD078736 中检出 MRE11 K510。

所以旧结果缺少这些蛋白的主要原因是：**旧文章对应的作者检索结果没有把这些蛋白报告为 Kla，而不是定位概率阈值 0、GO 交集或去重算法误删。**

## 11. 当前算法是否合理

### 11.1 合理之处

- 不把 raw 文件名称或普通蛋白存在当作 Kla 证据。
- 针对 MaxQuant、PEAKS、DIA 和作者补充表分别使用与格式匹配的解析规则。
- 保留样本级证据，最后才合并，避免不同条件和重复被过早混合。
- 统一 BaseAccession，避免 isoform 和重复实验造成蛋白重复计数。
- PXD053474 不盲目接收所有 search-only 结果，而采用 S3 或 DDA+DIA 双模式支持的规则。
- marker 156 只作辅助，不单独定义 Kla。
- GO 在 Kla 判定之后进行，不会反向影响 Kla 提取。
- 每条记录保留来源文件和行号，可回到原表核查。
- 14 项自动测试已通过，包括 GO NOT 排除、人源 fallback、PEAKS accession 映射、site ID-position 配对、Venn 区域求和和肿瘤特异表一致性。

### 11.2 仍需说明的限制

- 不同 PXD 使用的检索软件、FDR 设置和实验富集策略不同，蛋白数量不能简单解释为生物学丰度差异。
- PXD050470 使用作者补充表，证据粒度与直接搜索结果不同。
- PXD050470 的 27 个 raw 未在本地重新搜索；本项目不能声称独立重现其谱图鉴定。
- 缺失定位概率时，只有在来源文件本身已经明确报告 Kla 位点的情况下才保留。
- 当前 Venn 是“是否检出”的集合比较，不考虑定量强度、效应方向或统计显著性。
- 对未来新增数据不能只凭文件扩展名自动判断；必须先确认软件格式和 Kla 字段，再选择现有解析器或增加新的适配器。

总体判断：**算法适合作为跨数据集的可追溯 Kla 蛋白集合整理和 GO/Venn 比较，但不等同于对所有 raw 数据进行统一参数的独立质谱重搜索。**

## 12. 主要结果文件位置

- 数据集总表：`reanalysis/results/tables/dataset_analysis_summary.csv`
- 每个 PXD 的 sample-level Kla：`reanalysis/intermediate/kla_by_dataset/`
- 每个 PXD 的去重 Kla 蛋白：`reanalysis/results/tables/per_pxd/`
- 每个 PXD 的 Kla 与 GO 交集：`reanalysis/intermediate/go_intersection/`
- 全部 Kla 蛋白：`reanalysis/results/tables/all_unique_kla_proteins.csv`
- GO 未匹配蛋白：`reanalysis/results/tables/go_unmatched_kla_proteins.csv`
- 精确 Venn 区域：`reanalysis/results/tables/venn_regions/all_kla/` 和 `reanalysis/results/tables/venn_regions/kla_go_ddr/`
- 肿瘤特异 Kla：`reanalysis/results/tables/tumor_specific_kla_proteins.csv`
- 肿瘤特异 Kla-GO：`reanalysis/results/tables/tumor_specific_kla_ddr_proteins.csv`
- PXD028488 目录审计：`reanalysis/results/tables/pxd028488/`
- PXD053474 检索与 S3 对照：`reanalysis/results/tables/pxd053474/`
- MRE11/XLF/NBS1 源表追踪：`reanalysis/results/tables/target_protein_source_level_audit_MRE11_XLF_NBS1.csv`
- 排除日志：`reanalysis/results/tables/exclusion_log.csv`
- 自动测试日志：`reanalysis/logs/test_run_2026-07-22.log`

## 13. 可复现命令

```bash
cd /Users/gzy2520/Desktop/Research/kla
PYTHONPATH=reanalysis/scripts python3 reanalysis/scripts/run_pipeline.py --project-root /Users/gzy2520/Desktop/Research/kla
PYTHONPATH=reanalysis/scripts python3 -m unittest discover -s reanalysis/tests -p 'test_*.py' -v
```

正式分析入口为 `reanalysis/scripts/run_pipeline.py`。不同 PXD 的具体 Kla 解析规则位于 `reanalysis/scripts/extractors.py`，通用 accession、Kla 编码和 GO 处理函数位于 `reanalysis/scripts/common.py`。
