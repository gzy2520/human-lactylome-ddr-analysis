# 31 个全蛋白组位置的原表算法审计

本报告比较旧柱状图的组内去重集合、旧样本脚本和本次按实际源表重建的结果。基准图包为 `deliverables/20261011_lactylation_ddr_tall_verified`，原表位于项目主目录的 `data/`。逐位置结论见 `all_31_whole_proteome_review.csv`；其 Historical 字段指冻结的 `group_summary_31.csv`，具体最早提交历史另见 `audit/20261011_tall_algorithm_history/earliest_barplot/REPORT.md`。

## 原表重建结果

19 个位置直接重读作者原表得到 200 个单位；其余 12 个位置重读 MaxQuant、ProCan 和归档 evidence 表得到 27 个主分析单位，另保存 2 个 HEK QC 替代来源。共 31 个位置、227 个原图单位的逐单位 UniProt 集合和 DDR 集合与核查时图包逐项相同（`raw_rebuild_vs_current_unit_comparison.csv`）。这项结果确认数字可由原表复现；下述完整矩阵与条件混合仍需改变展示方式。

分母为选定检测单位内去重的 canonical UniProt，分子为其中直接匹配冻结 human、非 NOT GO 注释的蛋白。有效 6 字符和 10 字符 ID 均保留，isoform 后缀归一后去重。多数 MaxQuant 源使用 Majority protein IDs 并展开其有效成员；MCF10A evidence 使用 Leading Razor Protein，避免把每个共享肽的全部候选 Proteins 展开成额外检测蛋白。注释基表为 `data/publication_input/human_ddr_go_annotations.tsv`，直接 DDR ID 共 2660 个。

## 需要改为作者完整处理矩阵的三个全蛋白组位置

- Hippocampus：S4 原表 6082 行、H072/H081/H187 三列全部为正；表头说明为三个生物重复的相对定量强度。精确最小值 0.00262465560808778 分别重复 84、69、594 次。由该处理矩阵的正值无法恢复原始检出 mask。建议保留该 panel 蛋白集合及比例，按一个处理矩阵单位显示，省略样本 SEM；来源列数仍记录 3。
- Adjacent liver：ANTs sheet 的 4254 行在 8 列全部为正，Non-rec1 整列恰为 1。
- HCC：CISs sheet 的 5198 行在 8 列全部为正，Non-rec1 整列也恰为 1。这两个是相对比值矩阵，建议各按一个 panel 集合展示，省略由重复集合产生的零 SEM，来源列数各记录 8。同名 Non-rec/Rec 列来自不同 sheet、不同材料，来源键应包含 Sheet。

以上证据在 `complete_channel_minima_audit.csv`、`hippocampus_workbook_preface.txt`。本审计未推定具体缺失值填补算法。

Tendon 则有实际零值：有效行 3709，三个 RCT reporter 列各 3336 个正值、373 个零值，正值 mask 完全相同。六个 reporter 通道并集也相同。其零 SD 是该源表定量集合一致，保留这一事实即可。

## 旧组集合与当前样本正值集合的差异

|来源|旧全部有效鉴定集合（分母/DDR）|当前选定列正值并集|原表原因|
|---|---:|---:|---|
|Tendon|3984/154|3574/140|373 个有效 protein-group 行在全部六个 reporter 通道均无正值，展开后减少 410 ID、14 DDR。并非只筛选 RCT 而漏掉正常通道。|
|HepG2|6772/450|6713/448|59 ID 在九个 WT/SIRT1KO/SIRT3KO 强度列均无正值。|
|Sperm|10680/461|9492/385|旧算法取全部 PG.ProteinAccessions；当前要求 33 个 L/N 样本至少一列正值，减少 1188 ID、76 DDR。|

这些数值均由实际原表重算，见 `raw_selection_variant_counts.csv` 和逐 ID 的 `identified_but_not_selected_ids.csv`。旧 `workflow/source_legacy/R/analyze_ddr_fraction.R` 对前两种 MaxQuant 来源调用未指定 abundance_pattern 的鉴定集合提取，对 sperm 直接提取全部 PG.ProteinAccessions；旧样本脚本则逐列使用正值筛选。

Lung 与 placenta 的变化还涉及计数对象：旧冻结分母分别为 14022、12395 个 SourceProteinID，主要为 ENSP；当前映射后去重 UniProt 分母为 12358、10847，DDR 为 535、434。一个 ENSP 映射多 UniProt 时，逐 UniProt 查 GO；若把旧源 ID 的 any-DDR 布尔值传播给所有映射对象，DDR 会变为 606、502，分别多 71、68。映射表保留 1426/1301 个一对多源 ID 和 595/546 个未映射源 ID，见 `identifier_mapping_issues.csv`、`legacy_method_comparison.csv`。各一个组织 profile 的 36 个采集在作者 proteinGroups 内已整合。

MCF10A 原 6 字符规则为 4839/242；保留有效 10 字符 Leading Razor ID 后为 4920/242，多 81 个真实有效 ID。一个制备的 10 个 fraction、每个 fraction 两次采集先合并，总共 20 个采集，分析单位为 1。PC3M 以 heavy Intensity H 1st/2nd 分别汇总，两个生物重复的 10/26 个采集合并；集合 8049/262。

ProCan 新取实际 run peptide-count 检出集合并在同一制备内合并两个仪器采集，旧柱状图则用平均表达 feature matrix 的集合。此数据选择解释多个细胞系的新旧集合差异，不应为追求旧数限制到平均矩阵的 feature universe。独立读取 replicate mapping 确认六个细胞系各三个 technical preparations，每个两个 unique runs；metadata 原文为技术三重复并在不同质谱仪重复测量。HCT116 两个显示位置共享同三制备，因此 registry 跨位置有重复行，单位置合并仍为两个 unique runs。HEK 的 lys、H002、H003 保持分开；主展示 lys 的 401 个 QC 采集合并为一个 source pool，7817/446。

## 标识符、污染物与代表字段

19 个直接来源采用原解析与严格解析双跑，严格解析检查冒号、分号、逗号、pipe、isoform、10 字符、REV/CON；全部 200 个单位集合一致。Whole TALL 源为 nonenrichment proteins.csv 的 Accession 与 Area Sample1–3，未发现 Kla PSM 源的冒号解析问题。

另外独立扫描全部原表的 Protein、Proteins、Majority、Leading Razor、Accession、Reverse/contaminant 字段。直接来源发现 75 条 flag 空的前缀记录，全部位于附属 Protein IDs 的非首位混合成员；未发现代表字段为隐藏 reverse 的情况。恢复来源 134 条类似记录均为非首位混合成员，其中 placenta 的一条 Majority 为 `ENSP00000340761.7_RNA;REV__ENSP00000263674.3_RNA`，只有 forward ENSP 进入映射。污染物/反向候选不进入稳定 ID 集合，不能因共享组中含一个此类候选就删除其有效 forward 代表。字段级证据见 `all_direct_identifier_fields.csv`、`all_direct_prefixed_rows.csv`、`recovered_all_identifier_fields.csv`、`recovered_prefix_without_flags.csv`。

恢复来源额外逐字段核查原始输入，独立检查 ProCan header ID 集合覆盖输出 ID、映射文件到 run/preparation 的对应，以及源 metadata；检查结果保存在 `independent_procan_registry_checks.csv`、`procan_metadata_evidence.txt`。原始文件校验和在 `source_sha256.csv`。

## 单位与分散度解释

HMC3 的 H0/H24 和 HCT116–Roseburia 的 NC116/R116 各是两个条件、各一列，建议保留条件点而省略跨条件的重复 SEM。HepG2 含 WT、SIRT1KO、SIRT3KO 各三列；sperm 有 L/N 材料组；GSC/NSC 的单位为不同细胞模型。相应分散度需标明来自条件或模型。HK2 两位置和 HCT116 两位置是真实参考复用。

HMC3 现算法要求 LFQ intensity 正值：H0 为 5422/387，H24 为 5946/412。若改用未做 LFQ 的 Intensity 正值，则分别为 6072/423、6079/419（`hmc_intensity_vs_lfq.csv`）。因此现比例具体描述 LFQ 可定量蛋白集合；是否改为原始强度检出属于选择指标的变化，本轮保留既有选列并明确记录。

## 文件与复现

`all_31_whole_proteome_review.csv` 是全 31 位置的结论主表，含 GroupID、Dataset、Status、ID/分子分母/单位结论，以及旧集合、旧样本、重建数值。`raw_rebuilt_whole_membership.csv.gz` 保存从原表重建的全部主分析成员。建议更改展示单位的三个完整矩阵尚在此文件保留各源列，用于追溯；后续图包应在 panel 层合并。

执行顺序：`audit_direct_sources.R`、`recover_uncovered_sources.py`、`audit_variants.R`、`scan_source_identifier_fields.R`、`independent_recovered_fields.py`、`hmc_quantity_comparison.R`、`assemble_review.py`。全部输出限定在本审计目录，原始数据和共享图脚本未修改。
