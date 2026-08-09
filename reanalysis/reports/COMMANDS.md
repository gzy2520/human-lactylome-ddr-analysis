# Reanalysis command log

Run on 2026-07-22 from `/Users/gzy2520/Desktop/Research/kla`:

```bash
python3 reanalysis/scripts/build_project_metadata.py --project-root /Users/gzy2520/Desktop/Research/kla
PYTHONPATH=reanalysis/scripts python3 reanalysis/scripts/run_pipeline.py --project-root /Users/gzy2520/Desktop/Research/kla
PYTHONPATH=reanalysis/scripts python3 -m unittest discover -s reanalysis/tests -p 'test_*.py' -v
```

Reference-proteome DDR comparison added on 2026-07-30:

```bash
# Historical command only. The script is archived because it used the retired
# PXD043880 CA1 hippocampus reference:
# PYTHONPATH=reanalysis/scripts python3 reanalysis/scripts/analyze_reference_proteome_ddr.py

# Spreadsheet builders use the bundled @oai/artifact-tool runtime.
node reanalysis/scripts/build_reference_proteome_selection_workbook.mjs
node reanalysis/scripts/build_reference_ddr_comparison_workbook.mjs

# 扩展版 DDR 占比图，仅使用蛋白 ID 映射，不使用 GeneSymbol 回退
Rscript reanalysis/scripts/build_ensembl_uniprot_mapping.R .
Rscript reanalysis/scripts/analyze_expanded_ddr_fraction_by_accession.R .
Rscript reanalysis/tests/test_expanded_ddr_fraction_by_accession.R .

PYTHONPATH=reanalysis/scripts /Users/gzy2520/miniconda3/bin/python3 \
  -m unittest discover -s reanalysis/tests -p 'test_*.py' -v

# 交付工作簿使用 Codex bundled @oai/artifact-tool 运行时
CODEX_DEPS=/Users/gzy2520/.cache/codex-runtimes/codex-primary-runtime/dependencies
ln -s "$CODEX_DEPS/node/node_modules" node_modules
"$CODEX_DEPS/node/bin/node" \
  reanalysis/scripts/build_37group_deliverable_workbooks.mjs \
  /Users/gzy2520/Desktop/Research/kla
unlink node_modules

Rscript reanalysis/scripts/build_final_manifest.R \
  /Users/gzy2520/Desktop/Research/kla
```

Final pipeline summary:

```text
Primary sample-level rows: 63,846
Unique Kla proteins: 3,112
Kla GO-DDR proteins: 275
Automated tests: 14 passed
```

Kla regulator quantitative heatmaps added on 2026-08-05:

```bash
Rscript reanalysis/scripts/analyze_kla_regulator_intensity.R \
  /Users/gzy2520/Desktop/Research/kla

Rscript reanalysis/tests/test_kla_regulator_intensity.R \
  /Users/gzy2520/Desktop/Research/kla
```

Quantitative heatmap summary:

```text
Sample groups with auditable Kla quantitative values: 37/40
PXD datasets with valid within-dataset Z-score heatmaps: 21
```

37-group ordinary-proteome reference update on 2026-08-07:

```bash
Rscript reanalysis/scripts/build_lactylome_reference_pairing.R \
  /Users/gzy2520/Desktop/Research/kla

Rscript reanalysis/scripts/analyze_expanded_ddr_fraction_by_accession.R \
  /Users/gzy2520/Desktop/Research/kla

Rscript reanalysis/scripts/analyze_kla_regulator_whole_proteome_intensity.R \
  /Users/gzy2520/Desktop/Research/kla

Rscript reanalysis/scripts/plot_four_class_area_proportional_venn.R \
  /Users/gzy2520/Desktop/Research/kla

Rscript reanalysis/tests/test_lactylome_acquisition.R \
  /Users/gzy2520/Desktop/Research/kla
Rscript reanalysis/tests/test_expanded_ddr_fraction_by_accession.R \
  /Users/gzy2520/Desktop/Research/kla
Rscript reanalysis/tests/test_kla_regulator_whole_proteome_intensity.R \
  /Users/gzy2520/Desktop/Research/kla
Rscript reanalysis/tests/test_four_class_area_proportional_venn.R \
  /Users/gzy2520/Desktop/Research/kla

node reanalysis/scripts/build_37group_deliverable_workbooks.mjs \
  /Users/gzy2520/Desktop/Research/kla
Rscript reanalysis/scripts/build_final_manifest.R \
  /Users/gzy2520/Desktop/Research/kla
```

本次普通全蛋白参照修正使用 PXD073311 `report.pg_matrix.tsv` 的
`A0h_1/A0h_2/A0h_3`，不使用旧 PXD009687 `prot.xml`。

Kla 纯白行审计与普通全蛋白热图复核于 2026-08-07 运行：

```bash
Rscript reanalysis/scripts/analyze_kla_regulator_intensity.R \
  /Users/gzy2520/Desktop/Research/kla
Rscript reanalysis/tests/test_kla_regulator_intensity.R \
  /Users/gzy2520/Desktop/Research/kla
Rscript reanalysis/tests/test_kla_regulator_whole_proteome_intensity.R \
  /Users/gzy2520/Desktop/Research/kla
Rscript reanalysis/scripts/build_final_manifest.R \
  /Users/gzy2520/Desktop/Research/kla
```

严格普通全蛋白参照策略于 2026-08-08 运行：

```bash
Rscript reanalysis/scripts/apply_exact_reference_policy.R \
  /Users/gzy2520/Desktop/Research/kla
Rscript reanalysis/scripts/build_lactylome_reference_pairing.R \
  /Users/gzy2520/Desktop/Research/kla
Rscript reanalysis/scripts/plot_kla_regulator_landscape.R \
  /Users/gzy2520/Desktop/Research/kla
Rscript reanalysis/scripts/analyze_kla_regulator_intensity.R \
  /Users/gzy2520/Desktop/Research/kla
Rscript reanalysis/scripts/analyze_kla_regulator_whole_proteome_intensity.R \
  /Users/gzy2520/Desktop/Research/kla
Rscript reanalysis/scripts/analyze_expanded_ddr_fraction_by_accession.R \
  /Users/gzy2520/Desktop/Research/kla
Rscript reanalysis/scripts/plot_four_class_area_proportional_venn.R \
  /Users/gzy2520/Desktop/Research/kla

Rscript reanalysis/tests/test_lactylome_acquisition.R .
Rscript reanalysis/tests/test_expanded_ddr_fraction_by_accession.R .
Rscript reanalysis/tests/test_four_class_area_proportional_venn.R .
Rscript reanalysis/tests/test_kla_regulator_intensity.R .
Rscript reanalysis/tests/test_kla_regulator_landscape.R .
Rscript reanalysis/tests/test_kla_regulator_whole_proteome_intensity.R .

PYTHONPATH=reanalysis/scripts /Users/gzy2520/miniconda3/bin/python3 \
  -m unittest discover -s reanalysis/tests -p 'test_*.py' -v
```

结果：Kla 定量 37/40 组；严格普通全蛋白参照 33/37 组；
6 个 R 测试和 19 个 Python 测试全部通过。

严格材料身份和取材粒度审计于 2026-08-08 运行：

```bash
Rscript reanalysis/scripts/apply_exact_reference_policy.R \
  /Users/gzy2520/Desktop/Research/kla
Rscript reanalysis/scripts/build_lactylome_reference_pairing.R \
  /Users/gzy2520/Desktop/Research/kla
Rscript reanalysis/scripts/build_strict_reference_material_identity_audit.R \
  /Users/gzy2520/Desktop/Research/kla
Rscript reanalysis/scripts/analyze_kla_regulator_whole_proteome_intensity.R \
  /Users/gzy2520/Desktop/Research/kla
Rscript reanalysis/scripts/analyze_expanded_ddr_fraction_by_accession.R \
  /Users/gzy2520/Desktop/Research/kla
Rscript reanalysis/scripts/plot_four_class_area_proportional_venn.R \
  /Users/gzy2520/Desktop/Research/kla

Rscript reanalysis/tests/test_strict_reference_material_identity.R .
```

PXD050470 使用同研究 Table S4 的 H072、H081、H0187 普通全蛋白强度；
旧 PXD043880 CA1 参照不再进入活动分析。HCC 与癌旁肝分别读取
PXD065775 的 CISs 和 ANTs 工作表。

老师审阅总表于 2026-08-08 生成并复核：

```bash
Rscript reanalysis/scripts/build_lactylome_reference_pairing.R \
  /Users/gzy2520/Desktop/Research/kla
Rscript reanalysis/scripts/build_strict_reference_material_identity_audit.R \
  /Users/gzy2520/Desktop/Research/kla
Rscript reanalysis/scripts/build_teacher_review_table.R \
  /Users/gzy2520/Desktop/Research/kla

/Users/gzy2520/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/bin/node \
  reanalysis/scripts/build_teacher_review_workbook.mjs \
  /Users/gzy2520/Desktop/Research/kla

Rscript reanalysis/tests/test_teacher_review_table.R \
  /Users/gzy2520/Desktop/Research/kla
Rscript reanalysis/scripts/build_final_manifest.R \
  /Users/gzy2520/Desktop/Research/kla
```

老师审阅总表固定为37个Kla组，其中33组有严格材料匹配且可定量的
普通全蛋白参照，4组不进入普通全蛋白配对分析。两个HK-2组实际读取
PXD072220 的 `amostra1/amostra3/amostra4 PG.Log2Quantity`，并跳过
文件前两行说明。

DDR占比柱状图按老师要求进行细胞系内聚类及参照显示去重：

```bash
Rscript reanalysis/scripts/analyze_expanded_ddr_fraction_by_accession.R \
  /Users/gzy2520/Desktop/Research/kla
Rscript reanalysis/tests/test_expanded_ddr_fraction_by_accession.R \
  /Users/gzy2520/Desktop/Research/kla

/Users/gzy2520/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/bin/node \
  reanalysis/scripts/build_37group_deliverable_workbooks.mjs \
  /Users/gzy2520/Desktop/Research/kla
```

33个Kla研究组全部保留。普通全蛋白参照按参照PXD、来源文件、总蛋白数
和DDR数生成唯一显示键，共显示30根蓝柱。HK-2、MCF7和HCT116各有一组
完全重复的共享参照，分别只显示一次；统计表中的样本级配对关系不删除。

普通全蛋白调控因子热图按老师意见去重并补充GCN5显示名：

```bash
Rscript reanalysis/scripts/analyze_kla_regulator_whole_proteome_intensity.R \
  /Users/gzy2520/Desktop/Research/kla
Rscript reanalysis/tests/test_kla_regulator_whole_proteome_intensity.R \
  /Users/gzy2520/Desktop/Research/kla
```

33条严格Kla到普通全蛋白配对仍保留；最终热图显示30个唯一普通全蛋白
参照行。HK-2、MCF7和HCT116各合并一条重复显示行。identifier中的KAT2A
即GCN5，分析仍使用UniProt Q92830，图上显示为`GCN5 (KAT2A)`。

最终回归测试：

```bash
for test_file in $(rg --files reanalysis/tests | rg 'test_.*\.R$' | sort); do
  Rscript "$test_file" /Users/gzy2520/Desktop/Research/kla || exit 1
done

PYTHONPATH=reanalysis/scripts /Users/gzy2520/miniconda3/bin/python3 \
  -m unittest discover -s reanalysis/tests -p 'test_*.py' -v
```

结果：8个R测试和19个Python测试全部通过。

四个无严格普通全蛋白参照的Kla组于 2026-08-09 从配对Venn和对应比较集合中排除：

```bash
Rscript reanalysis/scripts/analyze_kla_regulator_intensity.R \
  /Users/gzy2520/Desktop/Research/kla
Rscript reanalysis/scripts/plot_four_class_area_proportional_venn.R \
  /Users/gzy2520/Desktop/Research/kla
Rscript reanalysis/tests/test_kla_regulator_intensity.R \
  /Users/gzy2520/Desktop/Research/kla
Rscript reanalysis/tests/test_four_class_area_proportional_venn.R \
  /Users/gzy2520/Desktop/Research/kla
```

本次保留 37 组来源审计，但严格配对及 Venn 使用 33 组；排除的具体组为
PXD062720、PXD063047/severe preeclampsia placenta、PXD064038 和 PXD075014。
PXD063047/normal pregnancy placenta 不在排除列表中。Venn 四类数量改为
9、2、9、13；`venn_sample_group_scope.csv` 保存了 37 组审计范围与 33 组实际纳入范围。
Venn 集合、区域表、membership 表和 `four_class_venn_tables.xlsx` 已重新生成。

2026-08-09 双语热图、33组显式Venn和交接文档更新：

```bash
Rscript reanalysis/scripts/analyze_kla_regulator_whole_proteome_intensity.R \
  /Users/gzy2520/Desktop/Research/kla
Rscript reanalysis/scripts/analyze_kla_regulator_intensity.R \
  /Users/gzy2520/Desktop/Research/kla
Rscript reanalysis/scripts/analyze_expanded_ddr_fraction_by_accession.R \
  /Users/gzy2520/Desktop/Research/kla
Rscript reanalysis/scripts/plot_four_class_area_proportional_venn.R \
  /Users/gzy2520/Desktop/Research/kla
```

Kla 与普通全蛋白热图均输出 `_zh` 和 `_en` PNG/PDF。Kla 热图33行按普通
全蛋白30行唯一参照轴排列，并显示9/2/9/13四分类分区。四套Venn同时输出
`_33groups_zh` 和 `_33groups_en`，标题直接注明33组严格配对范围。

当前项目交接Prompt：

`/Users/gzy2520/Desktop/Research/kla/NEW_CHAT_PROJECT_PROMPT.md`

当前正文方法稿：

- 中文：`reanalysis/reports/METHODS_CURRENT_33GROUP_ZH.md`
- 英文：`reanalysis/reports/METHODS_CURRENT_33GROUP_EN.md`
- 写作结构参考：`reanalysis/reports/METHODS_STYLE_REFERENCE.md`

写作结构参考文献：

Hu Y, Xie M, Li Y, et al. Benchmarking clustering, alignment, and integration
methods for spatial transcriptomics. Genome Biology. 2024;25:212.
doi:10.1186/s13059-024-03361-0.
