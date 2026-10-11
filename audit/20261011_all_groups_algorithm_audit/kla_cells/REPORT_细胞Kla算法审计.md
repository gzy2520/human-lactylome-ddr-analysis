# 19 个细胞材料 Kla 算法审计

审计对象：`deliverables/20261011_lactylation_ddr_tall_verified` 的 19 个细胞材料位置、54 个样本／模型／技术合并单位。来源为本地原始处理表、作者补充表、历史 Python 解析器与旧样本 R 脚本。

## 结论

15 个材料组的当前样本蛋白集合与原表独立重建完全一致。4 个材料组存在未标记反向序列记录，影响 7 个样本；每个受影响样本的分母多计 1 个蛋白，DDR 分子不变。

| GroupID | 材料 | 原始记录证据 | 受影响样本 |
|---|---|---|---|
| KLA31_18 | HepG2 | `Lactyl_K_Sites.txt` 文件行 8238，Protein=`REV__O14732`，Reverse 与 Potential contaminant 空 | SIRT3KO_rep1、SIRT3KO_rep2 |
| KLA31_26 | HMC3 | `La(K)Sites.txt` 文件行 1887，Protein=`REV__Q14093`，排除标志空 | H0、H24 |
| KLA31_27 | pretreated HK-2 | `La (K)Sites.txt` 文件行 2589，Protein=`REV__Q14093`，排除标志空 | HK2 |
| KLA31_28 | HK-2 ctr/man | `La(K)Sites.txt` 文件行 1969、id=1967，Protein=`REV__Q14093`，排除标志空 | ctr_1、ctr_2 |

修正从原始行开始：先排除 Protein 或首个 Leading proteins 明确带 `REV__` / `CON__` 前缀的记录，再按原有检出规则构建蛋白集合。随后重新去重、标注冻结 DDR 并计算比例。7 个受影响样本均无被删除 accession 的其他合格正向记录支持，见 `affected_id_forward_support.csv`。

问题有两条进入路径：`first_accession_current` 和历史 `normalize_accession` 会先剥离 REV/CON 前缀；另一些旧样本函数读取 Proteins 列，而该列已经将反向条目的 accession 写成普通 Q14093/O14732。仅检查 Reverse 字段或仅改最终 ID 分隔规则均不能完整解决问题。

| 组别 | 当前均值 ± SEM（%） | 原行修正后均值 ± SEM（%） |
|---|---:|---:|
| HepG2（9 个源样本） | 11.696698 ± 0.033664 | 11.698228 ± 0.033907 |
| HMC3（H0/H24 两条件） | 11.341986 ± 0.065557 | 11.351657 ± 0.066391 |
| pretreated HK2（n=1） | 10.166359 | 10.171059 |
| HK2 ctr/man（6 个源样本） | 10.596034 ± 0.152656 | 10.599823 ± 0.154014 |

这些整体均值沿用当前图在各原材料组内汇总全部源单位的计算方法；Condition 字段继续保留各处理条件。

## 原表和算法核查

独立解析覆盖 MaxQuant、Spectronaut 和作者工作簿，以原始检出字段筛选行；输出 156,787 条合格原始记录的 ID 解析证据。PXD028488、PXD053474 另重新运行完整原始 Python 位点解析器；MCF 四线、HK2 ctr/man 和 RKO 同时与本轮重新运行的 Python 原始解析器逐样本比对精确集合。最终独立 membership 有 66,271 条样本—蛋白记录。

- UniProt ID：检查 6/10 字符 accession、分号/冒号/逗号、`sp|` / `tr|`、isoform 后缀；稳定 accession 为计算键。
- 蛋白组：MCF 四线与 HK2 ctr/man 使用源 Protein 代表；RKO 使用原顺序第一 Majority protein ID。没有按字母排序另选代表，也没有将该类蛋白组任意扩展成所有成员。
- RKO：独立展开阳性 Kla evidence 的 site IDs，与 proteinGroups 的 site ID—位置配对关联。结果与原始 Python parser 完全一致。
- 分子与分母：先构建每个源单位的 Kla 蛋白集合，DDR 分子为该集合中冻结 GO 注释命中的子集。所有 54 个单位均使用相同集合计算分子/分母。
- 全字段污染扫描：另外读取 72 份底层原表／映射文件，保存 REV/CON 及排除标志。存在合法代表带有混合污染物次级成员的蛋白组；其合法正向代表保留，污染物 token 不计入集合。除上述 4 条反向代表记录外，没有发现导致当前 19 组数值进一步改变的前缀污染记录。
- 技术与材料单位：HEK293T HCD 是碰撞能量，合为一个技术单位；PXD053474 亚细胞组分、分级及 DDA/DIA 合并为一个单位；NSC/GBM stem 分别保留 2/6 个模型；HUVEC 6、RKO 4、HepG2 9 样本保留各自条件。Roseburia 和 HMC3 两条件各一个源样本。没有把 fractions 当生物学重复。

## TALL 复核

重新运行完整原始 PEAKS 算法：PSM 分数、乳酸化位置、protein-peptides 起点与 accession 候选映射共同确定蛋白；不是直接以 proteins.csv 的正 Area 作为 Kla 成员。当前三个富集样本分别为 196/201/198 蛋白、21/25/22 DDR，与本轮原始解析器精确一致。

旧 R `read_tall_by_sample()` 直接读取富集 proteins.csv 正 Area，产生 1662/1686/1678 分母；旧冻结 membership 标注提供 29/30/31 的 DDR 分子，两者选集不同。该问题已在 tall_verified 版本修正，本轮确认修正后的完整原始 parser 路径。非富集臂的 Kla 位点另有来源证据，未混入三个富集样本的分母。

## 条件释义

PXD078736 元数据仅明确 `ctrl versus man with three replicates each`，并写明 `The biological meaning of man remains pending publication.` 现 SampleGroup 的 mannitol 与冻结英文标签的 hypoxia-reoxygenation 释义不一致。数值审计保留原始 ctr/man 样本名，条件释义需根据正式来源核定。

## 输出与复现

- `group_algorithm_audit.csv`：每个 GroupID + Dataset 的结论、变化幅度和原始来源。
- `old_R_current_independent_threeway.csv`：旧 R、当前交付与本轮独立重建的逐样本对照。
- `exact_set_comparison.csv`、`fresh_legacy_set_comparison.csv`：精确集合差异。
- `independent_membership.csv.gz`：全 19 组独立修正版 membership。
- `corrected_affected_membership.csv.gz`：仅 4 个受影响组，便于定向替换。
- `corrected_affected_sample_counts.csv`、`affected_group_mean_sem.csv`：7 个受影响样本及 4 个组前后数值。
- `confirmed_all_reverse_rows.csv`、`affected_id_forward_support.csv`：原始反向记录与其他正向支持核查。
- `raw_row_identifier_evidence.csv.gz`、`raw_headers.json`、`raw_source_rules.csv`：原表行、字段、检出规则。
- `raw_prefix_scan_coverage.csv`、`all_raw_prefixed_identifiers.csv.gz`、`unflagged_raw_prefixed_identifiers.csv`：额外污染前缀扫描。

在工作区根目录依次运行 `reparse_legacy.py`、`independent_audit.py`、`audit_raw_prefixes.py`、`finalize_audit.py`。Python 环境需要 pandas、numpy、openpyxl。本轮可用运行时：`/Users/gzy2520/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3`。新输出仅写入本目录，没有修改共享图件、输入或脚本。
