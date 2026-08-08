# 严格普通全蛋白参照更新报告

更新日期：2026-08-08

## 最终范围

- Kla 定量母集仍为 40 个候选样本组，其中 37 组可定量。
- PXD037371 的 3 个临床组因 TMT 通道不能可靠映射，仍不进入最终 37 组。
- 37 个 Kla 组中，33 组具有生物材料匹配且含逐蛋白强度的普通全蛋白参照。
- 普通全蛋白热图、DDR 配对柱状图和普通全蛋白 Venn 只使用这 33 组。
- Kla 热图和两类 Kla Venn 不因普通参照排除而删除 Kla 证据。
- 参照材料身份和实验状态现已分列审计；逐行结果见
  `strict_reference_material_identity_audit.csv` 及中文版。

33 个严格配对组的四分类数量为：

- 正常/非肿瘤组织：9
- 癌症组织：2
- 正常/非肿瘤细胞：9
- 癌症细胞：13

## 排除的普通参照

以下 4 个 Kla 组保留在 37 组 Kla 母集中，但不进入普通全蛋白配对分析：

| Kla PXD | 样本组 | 排除原因 |
|---|---|---|
| PXD062720 | bladder cancer cells treated with EPI | 下载的 proteinGroups 是 Kla 富集检索结果；2223 条有效蛋白中 2060 条有 La(K) site IDs，只有 487 条有 A_1/A_3 强度 |
| PXD063047 | severe preeclampsia placenta | 未找到重度子痫前期胎盘逐蛋白强度矩阵；早发型子痫前期不能静默等同 |
| PXD064038 | MEC and NEC ESCC groups | 可解析候选只有蛋白鉴定或肿瘤/正常比值，没有可审计逐蛋白绝对强度 |
| PXD075014 | AC16 control and hypoxia | PXD059985 是增殖/分化状态，不能替代 control/hypoxia；PXD075014 补充表是 Kla 富集肽段 |

ESCC 候选大包的未完成下载保留在
`data/PXD037295/metadata/txt.rar.partial_download_2026-08-08`，不进入分析。

## 新接入的精确参照

| Kla 样本组 | 普通全蛋白参照 | 强度字段或样本 |
|---|---|---|
| TALL-104 | PXD028488 | non-enrichment PEAKS Area，3 个样本 |
| HMC3 | PXD028737 | MaxQuant LFQ intensity H0/H24 |
| PC-3M | PXD022005 | SILAC `Intensity H` |
| glioblastoma stem cells | PXD069969 | 相同 6 个 GSC 模型 LFQ |
| neural stem cells | PXD069969 | 相同 2 个 NSC 模型 LFQ |
| HCC | PXD065775 | CISs 肿瘤组织 iTRAQ |
| adjacent liver | PXD065775 | ANTs 邻近组织 iTRAQ |

HK-2 两组继续按既定要求使用 PXD072220 的 amostra1、amostra3、amostra4
未处理对照。HUVEC 使用 PXD073311 同研究 A0h_1、A0h_2、A0h_3 普通
PG 矩阵。正常妊娠胎盘使用 PXD010154 组织图谱中的独立胎盘矩阵，不是肺矩阵。

人海马已从独立研究的 CA1 亚区 PXD043880 改为 Kla 论文同研究 Table S4。
Table S4 含同一批 H072、H081、H0187 三份海马样本的普通全蛋白强度，
共 6,082 个唯一有效 UniProt BaseAccession，不再进行 GeneSymbol 转换。

HCC 与癌旁肝虽然都引用 PXD065775，但实际读取子集不同：HCC 只读 CISs，
癌旁肝只读 ANTs，DNTs 不进入二者。其他共享工作簿的组织和细胞子集也已逐项
写入材料粒度审计表。

## 算法与质量控制

- 所有匹配仅使用去除 isoform 后缀的 UniProt BaseAccession。
- GeneSymbol 只用于显示和人工审计，分析回退数为 0。
- 普通全蛋白强度按每个样本的全蛋白排名换算为 0–100 百分位，再对重复取中位数。
- 未检出的调控蛋白记为 0；颜色为白、黄、橙、深红，数值越高越暖。
- 最终 33 行普通全蛋白热图没有全白行，也没有整行百分位塌缩；每行至少有 4 个调控蛋白为正值。
- 修复了单竖线蛋白 ID（例如 `P49327|FAS_HUMAN`）被清空的问题。
- 修复了全缺失强度分组使用 `max(..., na.rm=TRUE)` 产生 `-Inf` 和大量 warning 的问题。
- 六个 R 测试和 19 个 Python 测试全部通过。

## 与归档基线比较

- 全部 Kla Venn：8719 个蛋白，集合与区域均保持不变。
- Kla-DDR Venn：510 个蛋白，集合与区域均保持不变。
- 本轮海马材料粒度修正前后，普通全蛋白 Venn 并集由 24566 变为 24735；
  新增 175，移除 6。
- 普通全蛋白 DDR Venn 并集由 764 变为 765；新增 1，移除 0。
- 当前普通全蛋白四类集合数依次为 18468、8756、13285、14989。
- 当前普通全蛋白 DDR 四类集合数依次为 649、426、631、616。

旧图表保留在
`archive/reanalysis_2026-08-07_pre_37group_reference_fix`。
本轮迁移记录位于
`archive/2026-08-08_strict_reference_policy/MIGRATION_MANIFEST.csv`。
旧的“37 组均具有普通全蛋白参照”报告和约 191 MB 的工作簿检查缓存也已
移入本轮归档；当前三份交付工作簿已重建并通过公式错误及视觉检查。

## 主要输出

- `reanalysis/results/figures/kla_regulator_whole_proteome_relative_intensity_heatmap.png`
- `reanalysis/results/figures/cell_type_kla_vs_reference_ddr_fraction.png`
- `reanalysis/results/figures/four_class_venn/reference_proteome_four_class_venn_zh.png`
- `reanalysis/results/figures/four_class_venn/reference_proteome_ddr_four_class_venn_zh.png`
- `reanalysis/results/tables/kla_regulator_whole_proteome_intensity_availability_audit.csv`
- `reanalysis/results/tables/kla_regulator_whole_proteome_normalized_long.csv`
- `reanalysis/results/tables/cell_type_kla_vs_reference_ddr_statistics_accession_only.csv`
- `reanalysis/results/tables/cell_type_kla_vs_reference_ddr_statistics_accession_only_zh.csv`
- `reanalysis/results/tables/strict_reference_exclusion_audit.csv`
- `reanalysis/results/tables/exact_reference_selection_audit.csv`
