# 扩展版 Kla 与常规蛋白组 DDR 占比图

## 分析范围

新版图从现有 40 个候选样本组中纳入 37 个。纳入条件为：

1. Kla 证据文件能提取蛋白 ID。
2. 对应常规蛋白组能提取蛋白 ID。
3. 两侧都能通过 ID 映射与 GO repair/damage 注释比较。

GO-DDR 的最终交集判定没有使用 GeneSymbol 回退。PXD050470 对应的正常
CA1 海马参照 PXD043880 原表只有 symbol，因此先使用人源 reviewed UniProt
映射快照将 symbol 转换为 BaseAccession，再按 ID 与 GO-DDR 表匹配。

## ID 处理

- UniProt isoform 后缀先去除，得到 BaseAccession。
- PXD010154 健康组织蛋白组使用 Ensembl protein ID。
- 这些 `ENSP` ID 通过 Ensembl BioMart 显式映射到 UniProt。
- 映射缓存保存在：
  `reanalysis/config/ensembl_protein_to_uniprot_biomart.tsv`
- 25,726 个 Ensembl protein ID 中有 24,462 个获得 UniProt 映射。
- 未映射 ID 保存在：
  `reanalysis/config/ensembl_protein_unmapped_biomart.tsv`
- 分母保留原始蛋白 ID 数，避免多个 Ensembl isoform 汇总到一个 UniProt
  accession 后人为缩小总蛋白数。
- DDR 分子由蛋白 ID 显式映射到 GO 表中的 UniProt BaseAccession 后判定。
- PXD043880 的 2,092 个蛋白特征拆分为 2,153 个 symbol token；唯一映射到
  2,105 个 reviewed UniProt BaseAccession。歧义和未映射 symbol 不进入 ID 分母，
  但全部保存在映射审计表中。

## 排除数据

- PXD037371 的 normal liver、nonmetastatic HCC 和 lung-metastatic HCC：
  TMT 通道无法可靠对应到三个临床组。

## 海马结果

- Kla PXD050470：853 个 UniProt BaseAccession，其中 29 个属于 GO-DDR，
  占 3.40%。
- 正常 CA1 参照 PXD043880：2,105 个映射后的 reviewed UniProt
  BaseAccession，其中 83 个属于 GO-DDR，占 3.94%。
- 两侧最终均按 BaseAccession 与去除 `NOT` 注释后的人源 GO-DDR 表匹配。

## 输出

- 主图：
  `reanalysis/results/figures/cell_type_kla_vs_reference_ddr_fraction.png`
- 独立命名图：
  `reanalysis/results/figures/cell_type_kla_vs_reference_ddr_fraction_accession_only.png`
- 英文统计表：
  `reanalysis/results/tables/cell_type_kla_vs_reference_ddr_statistics_accession_only.csv`
- 中文统计表：
  `reanalysis/results/tables/cell_type_kla_vs_reference_ddr_statistics_accession_only_zh.csv`
- 纳入和排除审计：
  `reanalysis/results/tables/cell_type_kla_vs_reference_ddr_accession_only_audit.csv`
- Kla 蛋白 ID 明细：
  `reanalysis/intermediate/expanded_ddr_by_accession/kla_proteins_by_sample_group.csv`
- 常规蛋白组 ID 与映射明细：
  `reanalysis/intermediate/expanded_ddr_by_accession/reference_proteins_by_sample_group.csv`
- 海马 symbol 到 reviewed UniProt 的映射审计：
  `reanalysis/results/tables/PXD043880_hippocampus_symbol_to_reviewed_uniprot_mapping_audit.csv`

## 运行命令

```bash
Rscript reanalysis/scripts/build_ensembl_uniprot_mapping.R .
Rscript reanalysis/scripts/analyze_expanded_ddr_fraction_by_accession.R .
Rscript reanalysis/tests/test_expanded_ddr_fraction_by_accession.R .
```
