# 扩展版 Kla 与常规蛋白组 DDR 占比图

## 分析范围

新版图的来源审计从现有 40 个候选样本组中保留 37 个；当前严格配对比较和四分类
Venn 使用其中 33 个。纳入严格配对的条件为：

1. Kla 证据文件能提取蛋白 ID。
2. 对应常规蛋白组能提取蛋白 ID。
3. 两侧都能通过 ID 映射与 GO repair/damage 注释比较，并且普通全蛋白参照与 Kla
   材料和处理状态严格匹配。

GO-DDR 的最终交集判定没有使用 GeneSymbol 回退。PXD050470 对应的正常
PXD050470 海马使用同研究 Table S4 的 H072/H081/H0187 普通全蛋白强度，按
UniProt BaseAccession 与 GO-DDR 表匹配；不再使用旧的 CA1 替代参照。

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
- PXD062720、PXD063047/severe preeclampsia placenta、PXD064038 和 PXD075014：
  没有完全匹配且可审计的普通全蛋白强度参照，因此不进入当前 DDR 对照和 Venn；
  原始数据与排除审计仍保留。

## 海马结果

- Kla PXD050470：853 个 UniProt BaseAccession，其中 29 个属于 GO-DDR，
  占 3.40%。
- 同研究普通全蛋白参照 PXD050470：6,082 个 UniProt BaseAccession，其中 219 个
  属于 GO-DDR，占 3.60%。
- 两侧最终均按 BaseAccession 与去除 `NOT` 注释后的人源 GO-DDR 表匹配。

## 输出

- 主图：
  `reanalysis/results/figures/cell_type_kla_vs_reference_ddr_fraction.png`
- 中文主图：
  `reanalysis/results/figures/cell_type_kla_vs_reference_ddr_fraction_zh.png`
- 英文主图：
  `reanalysis/results/figures/cell_type_kla_vs_reference_ddr_fraction_en.png`
- 独立命名图：
  `reanalysis/results/figures/cell_type_kla_vs_reference_ddr_fraction_accession_only.png`
- 37组来源审计统计表：
  `reanalysis/results/tables/cell_type_kla_vs_reference_ddr_statistics_accession_only.csv`
- 37组来源审计中文表：
  `reanalysis/results/tables/cell_type_kla_vs_reference_ddr_statistics_accession_only_zh.csv`
- 33组严格配对英文表：
  `reanalysis/results/tables/cell_type_kla_vs_reference_ddr_statistics_accession_only_paired_33.csv`
- 33组严格配对中文表：
  `reanalysis/results/tables/cell_type_kla_vs_reference_ddr_statistics_accession_only_paired_33_zh.csv`
- 纳入和排除审计：
  `reanalysis/results/tables/cell_type_kla_vs_reference_ddr_accession_only_audit.csv`
- Kla 蛋白 ID 明细：
  `reanalysis/intermediate/expanded_ddr_by_accession/kla_proteins_by_sample_group.csv`
- 常规蛋白组 ID 与映射明细：
  `reanalysis/intermediate/expanded_ddr_by_accession/reference_proteins_by_sample_group.csv`
- 四分类 Venn 样本范围审计：
  `reanalysis/results/tables/four_class_venn/venn_sample_group_scope.csv`
- 四分类 Venn 集合、区域和 membership 表：
  `reanalysis/results/tables/four_class_venn/`
- 四套 33 组中文 Venn：
  `reanalysis/results/figures/four_class_venn/*_33groups_zh.png`
- 四套 33 组英文 Venn：
  `reanalysis/results/figures/four_class_venn/*_33groups_en.png`

## 运行命令

```bash
Rscript reanalysis/scripts/build_ensembl_uniprot_mapping.R .
Rscript reanalysis/scripts/analyze_expanded_ddr_fraction_by_accession.R .
Rscript reanalysis/scripts/plot_four_class_area_proportional_venn.R .
Rscript reanalysis/tests/test_expanded_ddr_fraction_by_accession.R .
Rscript reanalysis/tests/test_four_class_area_proportional_venn.R .
```

Venn 图标题直接标明严格配对 33 组、四分类数量 9/2/9/13 和排除 4 个无严格
参照组。四套分析分别为全部 Kla、Kla DDR、普通全蛋白、普通全蛋白 DDR。
