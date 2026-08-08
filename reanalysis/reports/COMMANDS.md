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
