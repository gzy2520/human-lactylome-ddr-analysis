# 乳酸化数据与常规蛋白组配对状态

更新日期：2026-08-07

## 口径

- “乳酸化候选”不等于“已下载且可进入分析”。镜像、PRM、机制验证、方法特异数据和只有 raw 的项目均单独标记。
- “常规蛋白组”指未做 Kla 富集的普通全蛋白组，不等同于健康正常组织。匹配时优先同研究同样本，其次精确细胞系或组织。
- 邻癌、BPH、疾病组织、原代细胞和永生化细胞系均保留原始生物学身份，不互相冒充。

## 当前结果

- 共检索到 92 个相关人源 PXD，其中 46 个全局乳酸化候选已全部给出数据类别和分析资格。
- 这不等于全球所有可用乳酸化文件均已下载；超大处理包、仓库缺失结果和专有格式均保留真实状态。
- 配对表按样本组拆为 56 行。
- 已取得全局 Kla 蛋白明细的 40 个样本组，普通非富集蛋白组和健康组织基线两列均已补齐并可计数。
- 6 个纳入样本组仍缺可审计 Kla 蛋白明细，涉及 5 个 PXD：PXD037530、PXD046344、PXD053029、PXD063945、PXD065831。
- 46 个要求纳入的样本组均已配置可计数的健康组织基线。

## 主要缺口

- 当前缺口集中在远程超大处理包、只提供 Spectronaut SNE 的项目，以及仓库没有可用处理结果的项目。
- 10 个可成对样本组使用健康器官代理普通蛋白组；这些记录均标成“健康器官代理”，不是同细胞系或同患者精确匹配。
- 健康组织列已经补齐；T-ALL 使用健康脾脏淋巴组织、HUVEC 使用健康动脉组织、宫颈使用阴道相近组织代理，均明确标注不是精确样本匹配。
- PXD038880/PXD050906 继续 hold；PXD077426 是 PXD078736 镜像；PXD058173 和 PXD065104 不属于全局 Kla。
- 超大处理包已登记远程大小和来源，但未伪装成已下载。

## 健康组织参考来源

- PXD010154：12种健康器官的MaxQuant proteinGroups，包括肺、胎盘、肝、胃、脑、膀胱、食管、心、子宫内膜、结肠、肾和前列腺。
- PXD016999：GTEx 32种正常组织定量图谱；本项目使用乳腺、未暴露皮肤、脾脏、主动脉和阴道组织列。
- PXD018212：40个健康人跟腱/胫骨前肌腱mzTab文件，唯一BaseAccession并集为648。
- PXD037660：4名健康口腔黏膜对照的MaxQuant蛋白组，唯一leading BaseAccession为1050。
- PXD043880：正常人CA1海马组织；PXD066517：正常人精子DIA蛋白组。
- PXD073311：同研究非PTM普通全蛋白PG矩阵；仅使用A0h_1、A0h_2、A0h_3基线重复，A6h不进入参照，共7794个唯一UniProt BaseAccession。

## 计数规则补充

- 病毒感染成纤维细胞Spectronaut乳酸化表按K(UniMod:378)、precursor/蛋白组q值不高于0.01、位点置信度大于0提取；这里没有使用0.75定位阈值。
- 正常组织蛋白数的口径随来源保留在表中；BaseAccession、leading protein和蛋白编码基因数不能静默混称。

## 输出

- reanalysis/results/tables/lactylome_and_reference_proteome_pairing_zh.csv
- reanalysis/results/tables/lactylome_and_reference_proteome_gaps_zh.csv
- reanalysis/results/tables/lactylome_reference_all_gaps_zh.csv
- reanalysis/results/tables/lactylome_group_two_reference_columns_complete_zh.csv
- reanalysis/results/tables/lactylome_dataset_decisions_zh.csv
- reanalysis/results/tables/lactylome_reference_pairing_summary_zh.csv
- reanalysis/results/tables/lactylome_and_reference_proteome_pairing_zh.xlsx
