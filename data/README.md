# 数据目录

## datasets

每个 PXD 使用独立子目录，保留原始压缩包、已解压检索结果和论文补充表。

- PXD014870、PXD028488、PXD050470、PXD053474：此前文章数据。
- PXD060185、PXD078013、PXD078736：本次纳入重分析的数据。
- PXD038880、PXD050906：数据保留，但因记录关系尚需确认，本次明确排除。

## articles

按 PXD 保存对应论文 PDF。PXD078013 和 PXD078736 尚无正式论文 PDF。

## annotations

`GO-repair+damage(human).tsv` 是本次 Kla 蛋白取 GO repair/damage 交集的
唯一 GO 输入文件，不按 GO evidence code 额外过滤。

## metadata

保存数据集年份、样本类型、处理方法、数据库记录和下载清单。实际分析分组
以 `../config/sample_group_catalog.csv`、`../config/lactylome_reference_pairing.csv`
和 `../config/four_class_sample_grouping.csv` 为准。
