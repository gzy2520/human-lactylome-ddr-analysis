# 前12个组织位置的Kla样本算法审计

## 结果

从9个实际来源文件独立重读和解析了12个材料位置、44个源样本列；没有直接使用现有membership作为输入。逐样本完成检测判据、accession清洗、去重、固定2660个GO-DDR注释、分子分母和集合对照。44个样本的独立membership与当前图数据逐键完全一致，未发现TALL式冒号漏拆导致的数值错误。**其中10个位置为原表复算一致，tendon存在统计对象变化，hippocampus来源量化表不能恢复原始样本检出。** 后两个问题需要在结果对象与误差线处理中体现。

| GroupID | 材料 | 源样本列 | 样本占比均值% | 审计状态 |
|---|---|---:|---:|---|
| 01 | Tendon | 3 | 4.960465 | 统计对象变化：全鉴定union→reporter正值检出 |
| 02 | Lung | 6 | 6.163456 | 原表复算一致 |
| 03 | Hypertrophic scar | 4 | 4.220028 | 原表复算一致 |
| 04 | Adjacent skin | 4 | 5.196381 | 原表复算一致 |
| 05 | Hippocampus | 3 | 3.399766 | 源表不能恢复原始检出；推荐量化panel单值 |
| 06 | Placenta | 3 | 3.914570 | 原表复算一致 |
| 07 | Sperm | 3 | 3.109480 | 原表复算一致 |
| 08 | BPH | 5 | 6.563144 | 原表复算一致 |
| 09 | Adjacent liver | 1 | 6.366197 | 原表复算一致；单患者组织 |
| 10 | ESCC | 6 | 6.998222 | 原表复算一致；MEC/NEC各3 |
| 11 | Prostate cancer | 5 | 6.785623 | 原表复算一致 |
| 12 | HCC | 1 | 6.729700 | 原表复算一致；单患者组织 |

## Tendon：571/30变成样本union464/23的原因

最初 `workflow/source_legacy/R/analyze_ddr_fraction.R` 的 `kla_set(PXD033146)` 调用 `extract_maxquant_sites`，使用所有通过reverse/contaminant筛选且全局localization>0的位点。该全部鉴定集合为571个蛋白、30个DDR。

当前样本流程增加各个reporter通道强度>0的要求：通道1为463/23，通道2与3为464/23。**通道4–6也各为464/23，全部6通道正值并集仍只有464/23**。原571集合中107个蛋白（其中7个DDR）在全部6通道均无正reporter。差异来自全鉴定集合与可定量reporter检测集合，而非简单的RCT/NC分组排除。数值及移除集合详见 `tendon_channel_scope.csv`、`tendon_legacy_union_removed.csv`。

RCT2与RCT3检出集合相同，但源位点正值行分别998与997；量化值不同，去重后protein存在集合相同。两者DDR占比相同由源表决定。当前仅在已出现正reporter的同一集合中数DDR，分子分母一致。

旧脚本注释称每channel的`___1/2/3`是3个technical acquisitions，应修改为位点表中不同修饰multiplicity量化列；实际fraction acquisition见源`summary.txt`的Fr1等。当前对三列取OR合并没有制造3倍样本。通道1–3→RCT1–3的标签沿用既定样本映射，当前源表本身只标通道号，本次文件复核没有新增外部临床通道对应证据。

## Hippocampus：SD=0来自作者处理后的量化panel

当前来源为 `prca2331-sup-0005-tables3.xlsx`。表头原文为“Table S3 Lactylated proteins with quantifiable Class I lactylated sites in hippocampi”，Intensity列解释为“Relative quantitative values of intensity for three biological replicates”。表中2579位点、853蛋白，三列2579行均为正数，无NA或0。

| 列 | 正值行 | NA/0 | 最小值 | 最小值出现次数 |
|---|---:|---:|---:|---:|
| H072 | 2579 | 0/0 | 0.09011405 | 73 |
| H081 | 2579 | 0/0 | 0.09011405 | 70 |
| H187 | 2579 | 0/0 | 0.09011405 | 1279 |

原论文说明：全部鉴定3650个lactylated sites中1167个位点在三个样本均鉴定；S3选的是“detected in at least two biological replicates”的2579个Class I量化位点。证据在本目录 `hippocampus_article.txt` 第168–182行。**S3各列全正值的模式没有保留原始每样本鉴定是否缺失的信息**。同一极小值重复出现是来源数值特征；本次材料不足以确定作者采用何种填补或下限处理，不能指定一种未核实的插补方法。

因此“按Intensity>0”得到三个完全相同853/29集合的计算虽可重复，但不是三个独立原始检出蛋白集合。建议将其作为一个 `quantification_panel` 值29/853=3.3997655334%，不画样本SEM；在来源元数据保留3个生物样本。推荐表见 `recommended_hippocampus_panel.csv`。若要恢复每样本检出，需原始site/PSM sample-identification表或作者Table S2；当前本地该项目只有S3、S4、S12和部分mqpar，没有S2或原始Kla检出表。

S12含656蛋白，全部是S3的子集，没有新增蛋白，无法补回原始样本检出。原legacy Python提取器同时读取S3/S12，当前只取S3，因此蛋白集合未因漏S12改变。`mqpar.xml`将H072的Fr1–Fr6归同一Experiment H072；当前把fraction合并成来源样本，没有把6个fraction当6个人。

## 其他来源的实际判据与识别层级

- Lung、Placenta、ESCC：MaxQuant Kla位点表，排除带reverse/contaminant/site-only标记的行，按对应样本localization>0识别位点。每行`Proteins`多ID以accession展开并集，与原legacy R定义一致。保留由同一肽支持的全部候选accession，未改为排序后选first。ESCC是MEC_1–3与NEC_1–3两个临床亚组，现有材料均值包括两组。
- Scar、Skin：Proteome Discoverer中Modifications匹配Lacty/Lactyl/La(K)，并按样本`Found in Sample`分类筛选。实际分类只有High、Peak Found、Not Found；当前同时保留High与Peak Found。此处检出包括峰级匹配记录，保留原作者/旧R量化检测定义，不改称每样本都由独立MS/MS确认。Master Protein Accessions按分隔符展开。
- Sperm：作者PTM表`PTM.ModificationTitle=lactylation`、`PTM.SiteAA=K`，分别使用Rep 1/2/3的SiteProbability>0；三个真实列独立解析。PTM.ProteinId在所选记录无多分隔符，6/10字符ID按标准保留。
- BPH、Prostate cancer：Spectronaut `R.Condition`对应NAT1–5、PCa1–5；每个condition对应一个不同`R.FileName`，全表10个不同raw；R.Replicate均0。按L-Lac(K)、K、SiteProbability>0建立该condition的Kla蛋白集合。当地metadata明确5例BPH与5例PCa组织。
- Adjacent liver、HCC：同一患者的Control和HCC强度列，按Intensity>0及Localization probability>0，各算一个组织样本；不将两条件作为同组重复，不画样本误差。

## ID与污染物检查

稳健解析与旧规则分别独立实现。扩展分隔符`;,:`，处理标准sp/tr pipe及NX_前缀，移除isoform后缀，保留合法6/10字符UniProt；带REV__/CON__/DECOY等前缀的token排除。

选中原记录中出现的冒号来自`CON__ENSEMBL:`和`CON__H-INV:`污染物标识，而非TALL那样的多个合法人UniProt列表。部分未整行标contaminant的混合蛋白行含CON__token和正常人ID；只丢污染物token、保留有效人ID，与旧正则实际结果一致。所有实际异常字符串及源行号保存在 `special_id_source_rows.csv`。扩展冒号解析没有增加或删除任何样本蛋白。当前每样本GO-DDR分子来自其自身最终蛋白集合，与分母使用同一检测规则，未再次套用该组冻存membership限制注释。

## 审计文件

`group_review.csv`列出12位置审计状态及建议；`sample_review.csv`提供44个样本的计数和判据；`independent_membership.csv.gz`是从原表重新提取的集合；`membership_differences.csv`为空差异表（含表头）。`source_evidence.csv`及`source_field_values.csv`记录来源字段、实际取值、ID结构；tendon与hippocampus另有专项CSV。`audit.R`、`source_followup.R`、`finalize.R`可复算全部审计输出。

本次只写审计目录，未改动现有共享脚本、输入及图。

## 所有原始ID字段的隐藏decoy补查

另行读取9个来源文件并扫描实际存在的Protein、Proteins、Leading proteins、Leading razor protein、Master Protein Accessions及PTM ID字段。Tendon、Lung、Placenta、ESCC中有REV__/CON__的`Protein`代表字段分别68、27、44、145行，但它们全部已由原始Reverse或Potential contaminant标记排除；Leading proteins中的带前缀记录也全部有对应排除标记。未发现代表ID带REV__/CON__而标志为空的隐藏decoy进入这12个位置。

`Proteins`多ID字段中确有带CON__token但未整行标污染的混合行（上述四来源分别6、3、1、8行）；其正常人ID仍有原始位点支持，而污染token被ID验证丢弃。按token过滤与旧结果一致。完整证据为 `all_id_fields_flagged_rows.csv`，复核脚本为 `hidden_flags.R`。
