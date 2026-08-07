# 2019 年前乳酸化数据专项复核

## 判定原则

2019 年前采集的普通蛋白组不能仅因数据库搜索出现 `+72.021 Da` 就作为乳酸化数据。
本项目要求此类旧数据至少具有实验时明确的乳酸化富集或正交验证。乳酸化通用抗体应写作
`pan anti-Kla`；`pan anti-Kac` 是乙酰化抗体，不能作为乳酸化富集证据。

## PXD014870

当前纳入的 Kla 数据中，只有 PXD014870 需要按 2019 年前数据进行专项审核。

- DCA 和 rotenone 检索结果中的原始文件名包含 `20181127`，说明相关质谱在
  2018 年 11 月 27 日采集。
- 2019 年原始论文 Methods 明确使用 `pan anti-Kla (PTM-1401)`：
  胰酶消化后，将 pan anti-Kla 抗体与 Protein A Sepharose beads 偶联，
  对乳酸化肽段进行免疫沉淀和洗脱，再进行 LC-MS/MS。
- 本地结果目录名称为 `MCF7_*_SILAC_Kla_IP`，与 Kla 免疫沉淀实验一致。
- 实验还包含 SILAC、`13C3` 乳酸或 `U-13C6` 葡萄糖示踪。
- 2022 年研究使用乳酸赖氨酸循环亚胺诊断离子重新检查 PXD014870，
  支持这些肽段属于真实乳酸化，而不是普通未富集蛋白的质量偏移误判。

结论：PXD014870 可以保留，但必须标记为
`pre-2019, pan anti-Kla enriched and later CycIm-revalidated`。

## 普通参照蛋白组

DDR 占比图中的蓝色柱是普通全蛋白组参照。部分参照质谱可能早于 2019 年，
但这些蛋白只用于普通蛋白组分母，不会被标记为 Kla，因此不受上述 Kla 纳入规则影响。

## 输出

- `results/tables/pre_2019_lactylation_dataset_review.csv`
- `results/tables/cell_type_kla_ddr_lactylation_evidence_audit.csv`
