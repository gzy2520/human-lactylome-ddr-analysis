### Task 8: 迁移参照 DDR 分析 → `analyze_ddr.R`（第三部分）

**Files:**
- Modify: `reanalysis/scripts/analyze_ddr.R`（追加 `--stage reference` 段）
- 参照（只读）：`reanalysis/scripts/analyze_reference_proteome_ddr.py`

**Interfaces:**
- Produces: `--stage reference` 输出与 analyze_reference_proteome_ddr.py 一致：`reference_proteome_all_proteins.csv`、`reference_proteome_ddr_proteins.csv`、`reference_proteome_go_unmatched.csv`、`reference_proteome_source_manifest.csv`、`cell_type_reference_ddr_statistics.csv`、`cell_type_reference_control_information.csv`、`tall104_surrogate_ddr_sensitivity.csv`、`cell_type_kla_vs_reference_ddr_statistics.csv`、`cell_type_kla_vs_reference_ddr_fraction_v2.png/pdf`、`pxd030304_sample_audit.csv`（路径均在 tables/figures 原位）。

- [ ] **Step 1: 翻译 8 个 Python 函数**

`load_pxd030304`、`load_pxd072220`、`load_pxd002400`、`load_pxd043880`（保留 Python 版 symbol→reviewed UniProt 转换逻辑）、`build_go_maps`、`annotate_go`、`build_reference_statistics`、`build_control_information`、`build_tall_sensitivity`、`plot_comparison`（ggplot 重绘 v2 图，允许渲染差异但数字一致）、`build_source_manifest` 逐函数翻译并入 `run_stage("reference", {...})`。注意 `load_pxd043880` 与 Task 7 的 `extract_pxd043880_reference` 是不同的函数（一个在 R 版 expanded 里、一个在 Python reference 里），两者输出不同文件，都保留。

- [ ] **Step 2: 运行并 diff + 测试点 18/19 验证**

Run:
```bash
cd /Users/gzy2520/Desktop/Research/kla
Rscript reanalysis/scripts/analyze_ddr.R . --stage reference
Rscript reanalysis/scripts/lib/verify_outputs.R . archive/refactor_baseline_2026-08-06/output_sha256_baseline.csv
Rscript reanalysis/tests/test_pipeline.R .
```
Expected: verify_outputs "OK"；test_pipeline.R 断言 18/19（reference 分母、tall surrogate）通过。

- [ ] **Step 3: 检查点**

---

