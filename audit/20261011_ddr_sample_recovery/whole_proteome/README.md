# 缺失全蛋白组观测的原始数据恢复

## 输出与计算

本次从本地原始处理表恢复12个显示位置，共27条主分析单位。`recovered_unit_values.csv`以`UseRole=primary`筛选这27条；另2条HEK293T标准品QC合并结果标为`qc_alternative`，仅供来源审计。

每个单位先按下述设计合并分级或重复采集，再以去除亚型后缀的UniProt accession去重。DDR直接取冻结`data/publication_input/human_ddr_go_annotations.tsv`中人源、非NOT、UniProtKB条目的2660个唯一accession。百分比为`100 × DDR蛋白数 / 检测蛋白数`。源ENSP仅使用冻结映射转为UniProt；一对多映射全部保留并逐accession判定DDR，未映射条目列在`identifier_mapping_issues.csv`。

|来源|单位定义|主分析重复数|
|---|---|---|
|PXD030304：MCF7、HCT116、A549、MDA-MB-468、T47D、RKO|每个技术制备的两次跨仪器采集先取检测蛋白并集|每细胞系3个技术制备；HCT116同一数据对应两个显示位置|
|PXD030304：HEK293T|`Control_HEK293T_lys`的401次QC采集取并集；H002/H003标准品另存|1个lys参考来源，无SD/SEM|
|PXD010154：肺、胎盘|P013163与P013680分别为单一命名组织profile，各含36次分级采集，使用已按实验合并的proteinGroups|各1个来源profile，无SD/SEM|
|PXD022005：PC-3M|`Intensity H 1st`、`Intensity H 2nd`分别对应两个生物学重复；每重复的分级已在proteinGroups中合并|2个生物学重复|
|PXD002400：MCF10A|10个IEF fraction各2次LC-MS采集，20个raw file合并；使用`Leading Razor Protein`代表蛋白|1份来源材料，无SD/SEM|

## 来源依据与技术重复识别

PXD030304的本地`metadata/project.json`明确说明细胞系采用technical triplicate，并在不同质谱仪重复采集。`ProCan-DepMapSanger_mapping_file_replicates.txt`给出每个run的细胞系；同一制备标识在两次采集中保持相同。旧式文件名（含`_e0022_`）使用细胞系编号后的1/2/3，较新文件名使用第三个下划线分隔字段（如MCF7的00556、005a4、005cs；HCT116的00diz、00dkr、00mxq）。每个技术制备均核对为两次采集。具体对应保存在`run_to_unit_registry.csv`，不把采集重复当作新的制备。

ProCan检测来自`ProCan-DepMapSanger_peptide_counts_per_protein_per_sample.txt`中正肽段计数。使用实际run检测集，未限制在旧averaged表达矩阵的6692特征内，因此检测并集可能大于旧均值矩阵的非缺失蛋白集。

PC-3M通道与重复类型经[原论文](https://pmc.ncbi.nlm.nih.gov/articles/PMC10372916/)的“Succinylated proteins are identified by a global succinylome analysis”结果段核对：PC-3M为heavy标记，实验包含两个生物学重复。zip内`summary.txt`的`Experiment`字段将所有分级采集连接至`1st`或`2nd`。

PXD002400本地元数据说明IEF分级且每份分析重复采集；`summary.txt`确认MCF10A为20111109的F1–F10，每fraction有1/2两个采集。`evidence.txt`的`Proteins`包含共享肽的众多候选蛋白，因此采用`Leading Razor Protein`代表赋值，排除Reverse/Contaminant。原搜索的`Protein FDR`参数为1，当前恢复保留作者已有肽证据，未重新进行蛋白水平FDR估计。

## 与旧解析的具体差异

- 肺与胎盘恢复的UniProt并集分别为12358、10847，与冻结映射并集一致。直接GO判定DDR分别为535、434。旧解析将ENSP的IsDdr布尔传播给其所有映射UniProt，得到606、502；旧算法数字保存在`legacy_method_comparison.csv`，推荐使用直接GO结果。
- MCF10A恢复4920个有效代表UniProt，其中242个DDR。旧正则只接受6字符accession，得到4839个蛋白；新增81个合法10字符accession，DDR数保持242。差异来自标识符规则。
- PC-3M两个重复的并集为8049个蛋白、262个DDR，与冻结原始集合完全一致；本次分解为两个实际生物学重复。
- ProCan从实际肽检测run计算每制备的集合；旧均值强度矩阵经过特征保留和汇总，其分母与检测集合有所不同。`recovery_summary.csv`记录各组新旧集合差异。
- HEK293T使用已明确标为lys的QC参考来源；H002/H003为不同标准品来源，分别列为替代审计记录，未与lys一起计算SD/SEM。

`recovery_summary.csv`中的Mean、SD和SEM由单位百分比计算。技术制备的误差线反映技术重复，PC-3M误差线反映两个生物学重复；n=1不计算SD/SEM。

## 文件与复现

- `recover.py`：从项目根运行`python3 audit/20261011_ddr_sample_recovery/whole_proteome/recover.py`。
- `recovered_unit_values.csv`：每单位的分母、DDR数、比例、单位类型、来源与筛选字段。
- `recovered_protein_membership.csv.gz`：159848行逐单位UniProt及直接GO判定。
- `run_to_unit_registry.csv`：raw file、制备、显示组的完整对应；同一HCT116参考在两个显示位置重复列出。
- `source_inventory.csv`：每来源的实际列与恢复单位数。
- `source_sha256.csv`：读取的主要源文件和冻结注释的校验值。
- `identifier_mapping_issues.csv`：未映射和一对多ENSP，保留原始标识符。
- `legacy_method_comparison.csv`：旧标识符与DDR传播规则的数值影响。

独立回算已核对29条主/替代单位的ProteinCount与DdrProteinCount均与membership逐项相符，run登记数与AcquisitionCount逐项相符，membership的`PXD + SampleGroup + SampleID + BaseAccession`无重复。
