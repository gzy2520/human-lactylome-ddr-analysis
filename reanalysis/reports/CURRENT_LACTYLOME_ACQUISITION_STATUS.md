# 人源乳酸化质谱数据获取进展

更新日期：2026-08-06

## 当前范围

- ProteomeXchange 多关键词去重后共 92 个相关 PXD，其中 46 个列为全局乳酸化蛋白组或候选。
- 当前阶段优先获取作者检索结果、位点表和论文补充表，不批量下载大型 raw，也暂不把新增数据并入旧 Venn/DDR。
- PXD038880/PXD050906 继续保留 hold，不进入分析。

## 已获取并完成基础校验的数据

- PXD036307：正常人肺组织；MaxQuant La (K)Sites + evidence + proteinGroups；两个检索压缩包完整且内容相同；已解压一份
- PXD054919：A549 肺腺癌细胞；论文补充表 MOESM2；仓库 Results.zip 无效；补充表可用；Results.zip 校验值匹配但内容全零
- PXD063047：重度子痫前期胎盘与正常妊娠胎盘；MaxQuant La (K)Sites + evidence + proteinGroups；检索压缩包可列出并解压；第二仓库副本已确认重复
- PXD064912：人精子；Spectronaut .sne；仓库 SHA-1 匹配
- PXD066054：前列腺癌组织与良性前列腺增生组织；Spectronaut PTMSiteReport + Identification；两套结果 ZIP 完整且仓库 SHA-1 匹配
- PXD075377：肝细胞癌组织与邻癌组织；MS_identified_information 逐位点表；结果 ZIP 完整

## 正常材料解释

- PXD036307 是健康正常生理人肺组织，也是当前最明确的新增正常组织乳酸化数据。
- PXD063047 含 3 例正常妊娠胎盘和 3 例重度子痫前期胎盘，可拆分健康对照和疾病组。
- PXD064912 是三个正常人精子样本，属于正常人源生物样本，但不是组织。
- PXD066054 的 BPH 是良性病变对照，不能写成健康正常组织。
- PXD075377 的邻癌组织是疾病研究对照，不能写成健康正常组织。
- PXD054919 是 A549 肺腺癌细胞乳酸化质谱，不是正常材料。

## PXD054919 原稿与数据

- 老师提供的 PDF 与项目归档论文 SHA-256 完全一致。
- 仓库 Results.zip 与提交 SHA-1 一致，但文件全部为 0x00，无法解压，属于仓库提交质量问题。
- 论文补充表 MOESM2 可用，包含 A549 三个重复、3,110 个 Kla 位点和 1,220 个唯一蛋白。

## 主要交付

- 人源乳酸化质谱总表：reanalysis/results/tables/human_lactylome_mass_spectrometry_inventory.xlsx
- 92 个 PXD 的机器可读总表：reanalysis/results/tables/human_lactylome_mass_spectrometry_inventory.csv
- 仓库文件与来源 URL：reanalysis/results/tables/human_lactylome_repository_file_manifest.csv
- 当前已获取数据 QC：reanalysis/results/tables/priority_dataset_acquisition_qc_summary_zh.csv
- 当前下载文件登记：reanalysis/results/tables/priority_dataset_acquisition_manifest.csv

