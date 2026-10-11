# 现有样本 DDR 计算核查

2026-10-11。核查对象为 `deliverables/20261011_lactylation_ddr/inputs/proteome_samples.csv`，独立复算 TALL-104 乳酸化、PC-3M 乳酸化、ESCC 全蛋白组。DDR 注释取当前冻结 `human_ddr_go_annotations.tsv` 中人类 UniProtKB、非 NOT 注释，以去除 isoform 后的 UniProt accession 运算。逐样本去重，不以基因符号计数。

## 结论

- 272 条现有记录共引用 29 条不同 SourceFile；全部可在当前项目数据目录或工作区解析。详见 `source_path_audit.csv`。
- ESCC 全蛋白组 94 个 `ESCC-数字T` 样本从 `Dataset1.xlsx` 的 `2.a protein raw information` 按各列正值检出，分母、DDR 分子均与原记录逐项完全一致。均值 5.533467%，SD 0.164882 个百分点。来源是 PXD065830 外部 ESCC 肿瘤参考队列；没有逐样本配对到 PXD064038 的 6 个 Kla 样本。检测列为实际 T 样本列，不把 N 列或 pooled 列并入。
- PC-3M 原样本值复算通过：1=41/348=11.781609%；2=24/145=16.551724%；3=28/216=12.962963%。均值 13.765432%，SD 2.484244 个百分点。与组汇总 10.68% 的差异来自集合及汇总方式，未发现这三个样本 DDR 标注遗漏。`mqpar.xml` 将三个 raw 分别映射 experiment 1、2、3，fraction 均为 32767，不应因文件编号直接合并。单样本定位概率大于零的条目有有效样本 Score，支持作为已鉴定位点计数；若另要求 Intensity>0，均值会变成14.113727%，这是附加量化检出条件，不应在其他组仍保留定性鉴定时单独改变。
- TALL-104 原值错误：旧脚本 `read_tall_by_sample` 从富集 `proteins.csv` 对所有 Area>0 蛋白计数，包含大量没有 Lac 证据的蛋白；又以该组冻结 Kla membership 给这些蛋白赋 DDR，其余全部置 FALSE。3 个旧分母为1662/1686/1678，旧分子29/30/31。仅将全部富集蛋白正确赋 DDR 后分子为114/112/117，但这仍是富集表全部蛋白而非乳酸化蛋白。

## TALL-104 推荐修正

按 `DB search psm.csv` 中 Lac PSM、-10lgP≥20、排除 decoy/contaminant，并要求与 `protein-peptides.csv` 同一 accession 的肽段位置匹配，按 `Source File` 取 accession 去重集合。此判据复现 `workflow/source_legacy/python/utils/extractors.py` 的原始解析规则，获得：

| source raw | Kla 蛋白 | DDR-Kla 蛋白 | 占比 |
|---|---:|---:|---:|
| TALL-NaLAC-enrichment-1.raw | 173 | 21 | 12.138728% |
| TALL-NaLAC-enrichment-2.raw | 178 | 24 | 13.483146% |
| TALL-NaLAC-enrichment-3.raw | 176 | 22 | 12.500000% |

均值 12.707291%，SD 0.695767 个百分点。选择 `Case=TALL_Kla_legacy_mapping`。仅依据PSM accession而不校验肽段映射会得到175/22、184/25、182/23，保存在 `Case=TALL_Kla_corrected` 供比对，不作为最终推荐值。Source File 与 Sample 1–3 的末尾编号对应，但本次核查的原始表和文章未找到逐 raw 的生物学独立性表；图注可将其称为 3 个源样本/分析文件观测，避免沿用旧交接文档仅由列名推定的“biological samples”。

作为检测方法敏感性记录，`protein-peptides.csv` 先筛 Lac 肽段再以每 Sample Area>0 判定，结果为184/20、190/22、195/22，均值11.243521%。此表含跨样本量化与蛋白映射，其集合与 PSM 直接鉴定集合不同。当前建议采用直接 PSM 判定，与其他 Kla 位点鉴定表的观察对象更一致。

## 可审计产物

- `audit_existing.R`：独立复算脚本。
- `recomputed_counts.csv`：每样本分子分母、比例、计算方法、实际来源。
- `recomputed_membership.csv.gz`：每样本 UniProt membership 和 DDR 标记；用 Case 字段选定计算方法。
- `recomputed_summary.csv`：均值和 SD，单位为百分比/百分点。
- `escc_old_vs_recomputed.csv`：94 例 ESCC 与原结果逐项核对。
- `source_path_audit.csv`：所有 SourceFile 的当前存在性与解析路径。
- `tall_source_files.csv`：PSM 原始文件及记录分布。

所有复算文件仅写入本审计目录；未改动原数据、图或共享脚本。

## PSM 质量字段复核

TALL PSM 表全部1033条均标记PEAKS DB；-10lgP最小21.20，没有decoy/contaminant accession。原论文Methods（本地文本3383–3387行）说明PSM使用1% FDR，Lac要求AScore≥20并进行手工光谱检查。导出表仍包含AScore<20记录；旧解析器提取AScore位点却没有执行≥20的阈值。追加该阈值的敏感性结果为139/15、147/20、143/18，均值12.328074%、SD1.424850，保存在 `Case=TALL_Kla_AScore20`。若后续选择论文严格定位阈值，应同步核查全部PEAKS来源；当前复现冻存分析原则的推荐列为legacy_mapping。已导出PSM的FDR遵循作者流程，本表没有独立q-value列。
