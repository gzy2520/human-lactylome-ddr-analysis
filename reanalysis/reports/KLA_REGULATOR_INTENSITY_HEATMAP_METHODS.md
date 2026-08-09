# 乳酸化调控蛋白定量热图说明

## 结论

现有质谱结果中确实存在可用的定量信息，包括 MaxQuant 位点 `Intensity`、
PEAKS 修饰肽段 `Area`、Spectronaut `PTM.Quantity`/`Precursor.Quantity`、
Proteome Discoverer 乳酸化修饰肽 normalized abundance，以及作者补充表中的位点强度。

40 个候选样本组中，37 个可以把定量值可靠关联到 Kla 蛋白或位点；其中 4 个没有完全匹配且
可审计的普通全蛋白强度参照，因此不进入配对 Kla 热图、普通全蛋白热图或后续对照比较。
最终配对热图范围为 33 个 PXD+样本组，并按普通全蛋白热图的行、列顺序排列。
PXD037371 的三个组织组因无法可靠拆分 TMT 通道而从配对图中排除；四个无参照组和这三个
无定量组均保留原始文件及审计记录，不从磁盘删除。

调控蛋白身份统一使用
`reanalysis/config/lactylation_regulator_uniprot_mapping.csv` 中的人源
UniProt BaseAccession。reviewed 状态只作注释，不参与纳入或排除；
GeneSymbol 只作为图上的显示标签，不参与蛋白命中、强度汇总或缺失判断。
当前身份匹配流程的 GeneSymbol 回退数为 0。

## 跨研究热图

不同 PXD 的原始强度不能直接拼接比较，因为仪器、富集方法、采集模式、搜索软件和定量层级不同。
主热图采用以下流程：

1. 在每个定量样本内，将同一蛋白的 Kla 位点、修饰肽段或蛋白层级信号求和。
2. 对信号计算 `log2(signal + 1)`。
3. 在该样本的全部可定量 Kla 特征中计算百分位。
4. 同一样本组有多个条件或重复时，对百分位取中位数。
5. 调控蛋白在定量样本中未检出时记为 0；没有可审计定量值时保持缺失并显示 `?`。

主热图采用白色到暖色的连续色阶：低百分位接近白色，高百分位为橙红色；
颜色方向只表示组内相对 Kla 信号强弱，不表示跨研究的原始强度或 Log FC。

不同软件的蛋白字段统一处理如下：

- MaxQuant：`Proteins` 或明确的 protein accession 字段。
- PEAKS：`Protein Accession`。
- Spectronaut：`Protein.Group`、`PTM.ProteinId` 或 `PG.ProteinGroups`。
- Proteome Discoverer：`Accession` 或 `Master Protein Accessions`。
- 作者补充表：优先使用表中的 protein accession，不使用 GeneSymbol 回退。

所有 accession 去除 isoform 后缀后与调控蛋白 BaseAccession 比较。

两套数据在复核后收紧了证据范围：

- PXD046800 只读取 `Lacty_PeptideGroups.txt` 中 `Modifications` 明确含
  `Lacty` 且在对应 HSP/NSP 样本检出的修饰肽，并使用
  `Master Protein Accessions`；不再把乳酸化富集蛋白表中的全部蛋白视为 Kla。
- PXD066351 只读取 `DPLa-MSstats_Input.csv` 中
  `PTM.ModificationTitle = Lac (K)`、`PTM.SiteAA = K`、
  `PTM.SiteProbability > 0` 且 `PTM.Quantity > 0` 的记录；不再读取混有未修饰肽段的
  整张 DIA result 蛋白定量表。

该热图只能解释为“蛋白在各自样本乳酸化蛋白组中的相对 Kla 信号位置”，不是蛋白表达量，
也不是 Log FC，不能据此判断上调或下调。

最终 Kla 热图中有 5 个定量可用样本组整行显示为白色，逐组原因已记录在
`results/tables/kla_regulator_intensity_pure_white_audit.csv`：

- PXD028488 的 HEK293T、HCT116、TALL-104，以及 PXD064912 的 human sperm：
  Kla 来源文件中没有任何 identifier 调控蛋白的 UniProt BaseAccession 命中，因此白色表示
  “有可用 Kla 定量数据，但目标调控蛋白未检出”，不是解析失败。
- PXD014870 的 MCF7：AARS2 在 4 个条件中仅 1 个条件有信号，组内百分位为
  `0, 0, 0, 54.90`，取中位数后为 0；因此该行白色是稀疏检出被当前组内中位数规则压到 0，
  不能解释为所有条件都未检出。

这 5 行与普通全蛋白热图此前的 Ensembl ID 漏映射不是同一个问题。Kla 版使用的来源字段
已经能解析到 UniProt BaseAccession，且 GeneSymbol 回退数仍为 0。

## 普通全蛋白组热图

第二张热图不使用 Kla 富集信号，而读取同研究普通全蛋白定量文件或独立匹配的正常组织
全蛋白参照。每个样本内部先对全部普通蛋白信号计算 `log2(signal + 1)` 和百分位，再提取
identifier 表中调控蛋白的百分位。蛋白命中仅按人源 UniProt BaseAccession；
普通蛋白源文件中的 `ENSP...`、`ENSP..._RNA` 和带版本号的 Ensembl protein ID
先通过 `config/ensembl_protein_to_uniprot_biomart.tsv` 转换为 UniProt
BaseAccession；GeneSymbol 只用于显示，不能作为命中回退。

本次普通全蛋白组版本使用的算法审计文件为
`results/tables/kla_regulator_whole_proteome_algorithm_audit.csv`。它明确固定了：

1. 排名对象是每个普通全蛋白定量样本中全部保留的蛋白特征，而不是调控蛋白子集。
2. 信号先计算 `log2(signal + 1)`，再用平均秩换算为 0--100 百分位。
3. 同一乳酸化样本组的重复或条件取百分位中位数。
4. 在有可用普通全蛋白组但未检出某调控蛋白时记为 0（白色）；整个样本组没有可用普通全蛋白组时保留缺失并显示 `?`。
5. 只用 UniProt BaseAccession 匹配；Ensembl protein ID 先按项目映射表转换；
   GeneSymbol 仅作为显示和人工审计字段。
6. 不以 Kla 富集强度替代普通全蛋白强度。

与旧版乳酸化位点图相比，保留了相同的“样本内排名、组内中位数、白色到暖色”的数值逻辑，
但输入表已完全切换为普通全蛋白定量表，并加入了重复特征检查和输入来源审计。
因此这张图表示“调控蛋白在普通全蛋白组中的相对信号位置”，不能解释为乳酸化强度或
跨研究可比较的表达量。

另外，旧版普通全蛋白脚本曾使用分组 `ifelse()` 计算百分位；由于条件是单个分组值，
会把整组百分位错误地压成同一个数。当前版本已改为逐组 `if (...) ... else ...`，
并在自动测试中检查同一样本内的百分位不能全部塌缩为一个值。

37个乳酸化定量样本组中，33组具有材料身份严格匹配且可解析的普通全蛋白强度；
另外4组保留在排除审计中，不进入配对热图、普通全蛋白热图或DDR对照比较。
此前有10行整行纯白，原因是PXD010154健康组织参照使用Ensembl protein ID，
旧脚本没有执行Ensembl到UniProt转换，并非这些组织没有普通蛋白。当前版本已补入该映射，
并输出 `results/tables/kla_regulator_whole_proteome_ensembl_mapping_audit.csv`。
参照关系以
`reanalysis/results/tables/cell_type_kla_vs_reference_ddr_statistics_accession_only_zh.csv`
和
`reanalysis/config/lactylome_reference_pairing.csv`
为准：乳酸化 PXD 只定义行对应的细胞或组织，普通全蛋白参照 PXD 才提供热图颜色。
因此最终普通全蛋白热图不再有整行 `?`；如果后续新增没有强度的候选数据，
必须先进入待处理审计表，不能用 Kla 信号替代。

例如 PXD014870 的 MCF7 乳酸化数据使用 PXD030304 的 MCF7 普通全蛋白矩阵；
PXD060185 的 MCF7 使用同一参照。显示层将两行合并为
`MCF7 · Ref:PXD030304 · Kla:PXD014870/PXD060185`，但样本级配对表仍保留
两条记录。HK-2 和 HCT116 的共享参照采用相同规则；因此33条严格配对在热图中
显示为30个唯一普通全蛋白参照行。

identifier 表中的 `KAT2A` 即常用名 `GCN5`，稳定分析 ID 为 UniProt
`Q92830`。为避免老师按常用名审图时误判为遗漏，热图显示名改为
`GCN5 (KAT2A)`；没有重复添加第二个蛋白，也没有改用 GeneSymbol 分析。

PXD073311 的 HUVEC 参照已改为同研究的普通全蛋白 PG 矩阵，仅使用
`A0h_1`、`A0h_2`、`A0h_3`，A6h不混入正常参照；该矩阵有 7,709 个阳性蛋白组行，
去重后得到 7,794 个 UniProt BaseAccession。
PXD050470 的海马普通全蛋白参照为同研究 Supplementary Table S4：

- 使用 H072、H081、H0187 三个样本。
- 使用 6,082 个普通全蛋白特征的强度。
- 表内直接提供 UniProt accession，不使用 symbol 回退。
- 旧 PXD043880 CA1 参照不再进入活动分析。

该图采用白色到暖色的连续色阶，只表示调控蛋白在各自普通全蛋白样本中的相对信号百分位，
不是 Kla 水平、跨研究表达量或 Log FC。

## PXD 内热图

对于同一 PXD 内具有多个可比较定量样本的数据，另行生成蛋白 Z 分数热图：

1. `log2(signal + 1)`。
2. 每个样本按全部 Kla 特征的中位数中心化。
3. 每个调控蛋白仅在至少两个样本中具有真实定量信号时，跨该 PXD 的样本计算 Z 分数。

Z 分数只允许在同一 PXD 内解释，不跨 PXD 比较。最终共有 21 个 PXD 满足绘图条件。

## PXD050470 海马体

作者 Supplementary Table S3 已经包含 2579 个 Class I Kla 位点、853 个蛋白以及三份生物学样本的
`Intensity_H072`、`Intensity_H081` 和 `Intensity_H187`。因此海马体不需要重新搜索原始谱图，
即可纳入计数和组内相对定量热图。

仓库另外提供 27 个 `.raw` 和一个 `mqpar.xml`。这些文件只在需要从原始谱图完整重现作者
MaxQuant 搜索时使用，不是当前计数和热图的必要输入。

## 重新核对后恢复的数据

- PXD028737 HMC3：旧 `txt.zip` 是不含 Kla 的普通搜索结果；补下载 `combined.zip` 后获得
  `La(K)Sites.txt`，包含 H0/H24 强度。
- PXD073311 HUVEC：第二套 Spectronaut 矩阵使用 `K(UniMod:2114)` 表示 L-lactyllysine，
  可恢复 1881 个蛋白和 A0h/A6h 各三个重复。
- PXD075014 AC16：iProX 的 `.pdResult` 下载地址失效，但论文开放补充表 Table2 含 Kla 肽段和
  6 个 TMT 通道，可直接用于计数和定量。

## 从图中排除的数据

- PXD037371 三个肝组织组：位点表包含 TMT reporter intensity，但缺少临床组与通道的可靠映射，
  因此不能拆分到 normal liver、nonmetastatic HCC 和 lung-metastatic HCC。
- PXD062720 / bladder cancer cells treated with EPI：当前 Kla 结果不能提供与该组严格匹配的
  普通全蛋白强度参照。
- PXD063047 / severe preeclampsia placenta：没有同时满足该具体胎盘状态和逐蛋白强度的严格参照；
  同一 PXD 的 normal pregnancy placenta 仍保留。
- PXD064038 / MEC and NEC ESCC groups：当前候选文件未提供可审计的逐蛋白绝对强度。
- PXD075014 / AC16 control and hypoxia：现有候选参照与该研究的处理状态不严格匹配。

四个有 Kla 定量但无严格普通全蛋白参照的组，在
`results/tables/kla_regulator_intensity_availability_audit.csv` 中标记为
`严格配对分析纳入=FALSE`，并在 `kla_regulator_intensity_plot_exclusions.csv` 集中列出。

## 输出

- `results/figures/kla_regulator_cross_study_relative_intensity_heatmap.png`
- `results/figures/kla_regulator_cross_study_relative_intensity_heatmap.pdf`
- `results/figures/kla_regulator_whole_proteome_relative_intensity_heatmap.png`
- `results/figures/kla_regulator_whole_proteome_relative_intensity_heatmap.pdf`
- `results/figures/kla_regulator_within_pxd_zscore_heatmaps.pdf`
- `results/tables/kla_regulator_intensity_availability_audit.csv`
- `results/tables/kla_regulator_intensity_plot_exclusions.csv`
- `results/tables/kla_regulator_intensity_id_mapping_audit.csv`
- `results/tables/kla_regulator_intensity_pure_white_audit.csv`
- `results/tables/kla_regulator_intensity_sample_level_long.csv`
- `results/tables/kla_regulator_normalized_intensity_long.csv`
- `results/tables/cell_type_kla_vs_reference_ddr_statistics_accession_only_paired_33.csv`
- `results/tables/cell_type_kla_vs_reference_ddr_statistics_accession_only_paired_33_zh.csv`
- `results/tables/four_class_venn/venn_sample_group_scope.csv`
- `results/tables/kla_regulator_within_pxd_zscore_long.csv`
- `results/tables/kla_regulator_whole_proteome_intensity_availability_audit.csv`
- `results/tables/kla_regulator_whole_proteome_normalized_long.csv`
- `results/tables/kla_regulator_whole_proteome_heatmap_rows.csv`
- `results/tables/kla_regulator_whole_proteome_heatmap_display_long.csv`
- `results/tables/kla_regulator_whole_proteome_hippocampus_id_mapping_audit.csv`
