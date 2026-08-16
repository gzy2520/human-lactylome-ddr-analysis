# 当前30组分析的数据来源与证据链

## 范围来源

当前项目以`PXD+SampleGroup`为分析单元。范围由
`config/sample_group_catalog.csv`、`four_class_sample_grouping.csv`、
`lactylome_reference_pairing.csv`和`main_analysis_scope_exclusions.csv`
共同确定，而不是由`data/`中存在的文件自动决定。

范围审计为：

- 40个候选样本组；
- 37组具有可追溯的逐蛋白Kla定量；
- 33组原本具有可用普通全蛋白参照；
- 老师排除3组后，当前主分析为30组；
- 30组对应28条唯一普通全蛋白参照展示行。

## Kla证据

Python核心证据构建器对7个历史核心PXD执行数据集特异解析，并输出统一的
sample-level Kla位点长表。R端的`build_kla_regulator_landscape.R`继续读取
配置中其他项目的MaxQuant、PEAKS、Spectronaut、Proteome Discoverer或作者补充表，
形成40/37组目录。

只有明确报告Kla位点、Kla修饰肽段或Kla定量的记录可作为Kla证据。原始谱图、
泛酰化富集、诊断离子或未标记修饰的普通肽段不能单独证明Kla。反库、污染物、
非人源、无有效蛋白ID和无法定位到蛋白的记录被排除。

异质研究不强制使用单一定位概率阈值；保留作者报告的有效Kla证据，并在来源表中
记录定位概率、PEP、score、数据层级、源文件和源行。PXD014870的`>=0.75`
结果仅为敏感性输出，且该MCF7组已按老师要求排除出当前30组。

## 蛋白标识符

集合分析统一使用去数据库前缀和isoform后缀的人源UniProt
`BaseAccession`。Ensembl protein ID必须通过缓存映射表显式转换。
GeneSymbol、ProteinName和reviewed状态仅供显示与审计，不参与命中回退。

## 普通全蛋白参照

参照必须为非Kla富集的普通蛋白定量矩阵。优先顺序为同一批样本、同研究匹配材料、
同材料基线或独立队列的相同材料/表型。材料身份匹配并不自动代表实验处理状态
完全相同，因此审计表分别记录：

- `MaterialIdentityMatch`
- `ExperimentalStateMatch`
- `ReferenceMatchQuality`
- `PairingCaveat`

当前30条配对关系和28条展示行分别见：

- `results/tables/strict_reference_material_identity_audit.csv`
- `results/tables/kla_regulator_whole_proteome_heatmap_rows.csv`
- `results/tables/kla_and_reference_teacher_review_zh.csv`

## GO和DDR

DDR集合来自`data/annotations/GO-repair+damage(human).tsv`。排除`NOT`限定符，
不按GO evidence code进一步缩减。Kla集合和普通全蛋白集合分别与DDR集合取交集，
并分别使用自己的全部蛋白数作为分母。

线性通路图使用另一条独立证据链：UniProt release 2026_02中399个当前蛋白的
全部直接BP、CC和MF term。七通路种子及边界见
`config/go_term_pathway_seed_rules.csv`；完整判定见
`docs/GO_TERM_PATHWAY_SCORING_30GROUPS.md`。

## 当前结果与历史材料

当前可发表结果由`workflow/run_pipeline.R selected_figures`生成，并由
`tests/validate_publication_contract.R`验证。`results/`可重新生成且不进入Git。
`archive/`和旧33组/507蛋白/人工评分/降维报告只用于历史追溯，不是当前来源说明。
