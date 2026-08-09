# Kla 项目索引

更新时间：2026-08-09
项目根目录：`/Users/gzy2520/Desktop/Research/kla`

## 目录

| 目录 | 当前用途 |
|---|---|
| `data/` | PXD 原始文件、检索结果、补充表、文章、元数据和 GO/identifier 输入 |
| `previous_umap/` | 保留的最新纯 GO repair/damage 版本 |
| `reanalysis/` | 当前有效代码、配置、中间表、结果、报告、日志和测试 |
| `archive/` | 已弃用脚本、旧结果和迁移清单；不得永久删除 |

## 当前最终分析范围

- 37 个样本组具有可审计 Kla 定量。
- 33 个样本组具有完全匹配且可定量的普通全蛋白参照，进入最终热图、DDR
  对照和四分类 Venn。
- 四分类数量：正常组织 9、癌症组织 2、正常细胞 9、癌症细胞 13。
- 普通全蛋白热图为 30 条唯一参照行；Kla 热图为 33 行，顺序按普通全蛋白轴展开。
- 排除的四个无严格参照组为 PXD062720、PXD063047/severe preeclampsia、
  PXD064038、PXD075014。
- PXD037371 三个临床组因 TMT 通道映射不可靠，不进入 Kla 定量范围。

范围依据：

- `reanalysis/config/lactylome_reference_pairing.csv`
- `reanalysis/results/tables/kla_regulator_intensity_availability_audit.csv`
- `reanalysis/results/tables/four_class_venn/venn_sample_group_scope.csv`

## 分析主键

蛋白分析统一使用去除 isoform 后缀的 UniProt `BaseAccession`。GeneSymbol 和
Protein name 只用于显示/审计，禁止参与命中、去重、合并或 GO 回退。

GO-DDR 输入：

`data/annotations/GO-repair+damage(human).tsv`

排除 `NOT` 注释；主分析暂不按 evidence code 缩减。

## 当前主流程

1. `run_pipeline.py` 与 `extractors.py` 从检索结果/作者补充表提取 Kla。
2. `build_lactylome_reference_pairing.R` 和
   `build_strict_reference_material_identity_audit.R` 固定一一对应的普通全蛋白参照。
3. `analyze_kla_regulator_whole_proteome_intensity.R` 生成 30 行全蛋白热图轴。
4. `analyze_kla_regulator_intensity.R` 按该轴生成 33 行 Kla 热图。
5. `analyze_expanded_ddr_fraction_by_accession.R` 计算 Kla/普通全蛋白 DDR 占比。
6. `plot_four_class_area_proportional_venn.R` 生成四套 33 组面积比例 Venn。
7. `build_final_manifest.R` 生成 SHA256 清单。

完整交接说明和重跑命令：

`NEW_CHAT_PROJECT_PROMPT.md`

## 当前正式结果

- 双语 Kla 热图：
  `reanalysis/results/figures/kla_regulator_cross_study_relative_intensity_heatmap_{zh,en}.png`
- 双语普通全蛋白热图：
  `reanalysis/results/figures/kla_regulator_whole_proteome_relative_intensity_heatmap_{zh,en}.png`
- 双语 DDR 柱状图：
  `reanalysis/results/figures/cell_type_kla_vs_reference_ddr_fraction_{zh,en}.png`
- 双语 33 组 Venn：
  `reanalysis/results/figures/four_class_venn/*_33groups_{zh,en}.png`
- 33 组统计表：
  `reanalysis/results/tables/cell_type_kla_vs_reference_ddr_statistics_accession_only_paired_33{,_zh}.csv`
- 老师审阅表：
  `reanalysis/results/tables/kla_and_reference_teacher_review_zh.csv`
- 最终 SHA256：
  `reanalysis/reports/final_file_manifest_sha256.csv`

方法细节：

- `reanalysis/reports/METHODS_CURRENT_33GROUP_ZH.md`
- `reanalysis/reports/METHODS_CURRENT_33GROUP_EN.md`
- `reanalysis/reports/METHODS_STYLE_REFERENCE.md`
- `reanalysis/reports/KLA_REGULATOR_INTENSITY_HEATMAP_METHODS.md`
- `reanalysis/reports/EXPANDED_DDR_FRACTION_ACCESSION_ONLY.md`
- `reanalysis/reports/DATA_SOURCE_AND_ANALYSIS_ALGORITHM.md`
- `reanalysis/reports/COMMANDS.md`
