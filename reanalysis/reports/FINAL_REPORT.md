# Kla 项目重构与重新分析报告

## 结论

项目已重构为 `data`、`previous_umap`、`reanalysis` 三个主要工作区，旧脚本和旧结果均
保留在 `archive`，没有永久删除数据。最终分析获得 3,112 个去重 Kla 蛋白，其中 275 个
与人 GO repair/damage 相交。相较第一阶段旧结果 2,983/263，分别增加 129/12；变化来自
更完整的数据目录读取和修正的证据规则，没有为追求旧数量而修改阈值。

## 数据和证据

| PXD | 年份 | 细胞或组织 | 主要处理 | 最终证据 | Kla 蛋白 | Kla-GO DDR |
|---|---:|---|---|---|---:|---:|
| PXD014870 | 2019/2022 | MCF7 | DCA、缺氧、rotenone、U-13C6 glucose | MaxQuant site + modified peptide | 193 | 16 |
| PXD028488 | 2022 | HCT116、T-ALL、HEK293T | enrichment/non-enrichment、HCD 24-40 | PEAKS PSM + protein-peptide | 665 | 65 |
| PXD050470 | 2025 | 人海马体 | 三个生物样本 | 作者 S3/S12 | 853 | 30 |
| PXD053474 | 2024 | HCT116 | 亚细胞/全细胞、DDA/DIA、富集/非富集 | 四套检索 + S3 | 590 | 67 |
| PXD060185 | 2025 | MCF7、MDA-MB-468、MCF10A、T-47D | A/B/C/D Kla 富集 | MaxQuant site table | 1,991 | 223 |
| PXD078013 | 2026 | RKO | GSK3B WT/KO | evidence + proteinGroups | 1,217 | 115 |
| PXD078736 | 2026 | HK-2 | ctrl/man | MaxQuant site table | 1,250 | 133 |

PXD038880/PXD050906 保留原始下载、文章和元数据，但按老师要求完全排除。PXD050470 本地
没有仓库所列 27 个 raw，本次未自行重新搜索，只使用作者补充表。

PXD014870 和 PXD028488 对应同一篇 2022 方法学文章，但不是重复下载：PXD014870 是 2019
MCF7 数据在该文章中的重新分析来源，PXD028488 是 2022 文章新沉积的 PEAKS/诊断离子和
碰撞能实验数据。

## 关键修正

PXD050470 旧解析器把 Excel 表头偏移了一行，导致读取为 0 行；按工作簿实际表头修复后，
S3/S12 并集为 2,600 个 accession-site。原文件完整，问题在旧读取代码而非下载格式。

PXD028488 共有 19 个目录：13 个进入主分析，5 个排除，1 个 all-HCD 聚合目录只用于验证。
旧流程只覆盖三个目录，且把聚合目录当主来源；新流程增加 638 个蛋白，最终为 665 个。
156.1025 diagnostic ion 只标记支持，不单独判定 Kla。

PXD053474 的 accession-site 对照为：search+S3 一致 1,275、search-only 826、S3-only 3。
search-only 中仅 20 个同时获得 DDA 和 DIA 支持，因此主结果为 S3 的 1,278 个位点加这
20 个位点，共 1,298 个；另外 806 个单模式候选保留在审计表。

## 目标蛋白

最初三篇文章对应的可用 Kla 数据中，MRE11、XLF/NHEJ1、NBS1/NBN 缺失不是阈值 0、
GO 交集或 isoform 去重造成的：

- PXD014870 和 PXD050470 的可用作者 Kla 表中没有这些目标记录。
- PXD028488 能找到 MRE11/NHEJ1/NBN 的普通蛋白或未乳酸化肽段，但 0 条带 Lactylation、
  `K(+72.02)` 或 `K(Lactyl)` 的目标记录。
- PXD053474 能找到这些蛋白/普通 PTM 肽段，但 0 条目标记录的 PTM 为 Lac。
- 新纳入的 PXD060185 明确检出 MRE11、XLF/NHEJ1、NBS1/NBN 的 Kla；PXD078736 检出
  MRE11 Kla。

因此，最初缺失主要是作者检索结果没有给出目标蛋白的 Kla 证据，而不是算法把有效 Kla
记录过滤掉。逐文件证据见 `target_protein_source_level_audit_MRE11_XLF_NBS1.csv`。

## Venn 结果

主分组使用 HEK293T、HK-2、MCF10A 作为永生化/非肿瘤模型，MCF7 归肿瘤组。

| 区域 | 全部 Kla | Kla ∩ GO-DDR |
|---|---:|---:|
| hippocampus only | 397 | 4 |
| normal/immortalized only | 231 | 15 |
| tumor only | 923 | 99 |
| hippocampus + normal only | 32 | 0 |
| hippocampus + tumor only | 157 | 6 |
| normal + tumor only | 1,105 | 131 |
| all three | 267 | 20 |

老师分组和常规生物学分组在最终确认后使用相同成员，因此两套 Venn 数量差异为 0。仍保留
两套输出与 classification warning，明确 HEK293T/HK-2 不是正常组织。

## 质量控制

14 项自动测试全部通过：排除数据隔离、GO `NOT`/人源过滤、BaseAccession 规范化、PEAKS
accession 映射、site ID/position 配对、Kla 编码识别、PXD053474
对照规则、PXD014870 0.75 敏感性、七区域精确求和、区域内唯一性和肿瘤特异表一致性。
PXD014870 主分析为 193 蛋白/442 位点；0.75 敏感性为 185 蛋白/407 位点，敏感性结果未
混入主结果。严格化共享肽段和 site-position 边界规则后，主计数未发生变化。

## 主要交付路径

- 数据说明：`reanalysis/results/tables/dataset_analysis_summary.csv`
- 样本级 Kla：`reanalysis/intermediate/kla_by_dataset/`
- 每 PXD Kla-GO：`reanalysis/intermediate/go_intersection/`
- Venn 区域表：`reanalysis/results/tables/venn_regions/`
- Venn 图片：`reanalysis/results/figures/`
- 肿瘤特异表：`reanalysis/results/tables/tumor_specific_kla_proteins.csv`
- 肿瘤特异 DDR 表：`reanalysis/results/tables/tumor_specific_kla_ddr_proteins.csv`
- 目标蛋白追踪：`reanalysis/results/tables/target_protein_source_level_audit_MRE11_XLF_NBS1.csv`
- 迁移清单：`archive/migration_manifest_reconstructed_2026-07-22.csv`
