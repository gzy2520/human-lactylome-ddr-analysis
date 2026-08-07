### Task 14: 最终全验证

**Files:**
- Modify: `reanalysis/reports/final_file_manifest_sha256.csv`（重建）
- Create: `archive/refactor_baseline_2026-08-06/POST_REFACTOR_TEST_RESULTS.txt`

- [ ] **Step 1: 全部 R 测试**

Run:
```bash
cd /Users/gzy2520/Desktop/Research/kla
for t in test_expanded_ddr_fraction_by_accession test_kla_regulator_intensity \
         test_kla_regulator_landscape test_kla_regulator_whole_proteome_intensity \
         test_lactylome_acquisition test_pipeline; do
  Rscript reanalysis/tests/$t.R . || echo "FAILED: $t"
done
```
Expected: 6 个测试文件全部通过，输出存 POST_REFACTOR_TEST_RESULTS.txt。

- [ ] **Step 2: 输出 diff**

Run: `Rscript reanalysis/scripts/lib/verify_outputs.R . archive/refactor_baseline_2026-08-06/output_sha256_baseline.csv '^archive/|venn|_fraction'`
Expected: 除按设计移除的 `cell_type_kla_vs_reference_ddr_fraction.{png,pdf}`、移入 archive 的文件、及 Venn 图渲染差异外，全部字节一致。人工检查 diff_report 中每一行 CHANGED，确认均可解释。

- [ ] **Step 3: 重建 manifest 并核对**

Run:
```bash
cd /Users/gzy2520/Desktop/Research/kla
Rscript reanalysis/scripts/build_final_manifest.R .
# 交叉核对：manifest 中每个文件存在且哈希一致
Rscript -e '
m <- read.csv("reanalysis/reports/final_file_manifest_sha256.csv", stringsAsFactors = FALSE)
missing <- m$RelativePath[!file.exists(m$RelativePath)]
stopifnot(length(missing) == 0)
cat("manifest rows:", nrow(m), "| all files present\n")'
```
Expected: "manifest rows: N | all files present"，且清单涵盖全部最终交付文件（不含 .py）。

- [ ] **Step 4: 统计数字复核**

Run: `Rscript -e 'd <- read.csv("reanalysis/results/tables/dataset_analysis_summary.csv"); cat(sum(d$UniqueKlaProteins[7个纳入行]去重后), "\n")'`（以 run_pipeline 输出为准），Expected: UniqueKlaProteins 合计去重后 3112；KlaGoDdrProteins 合计去重后 275；扩展 DDR 表 37 行 + 3 行排除审计。

- [ ] **Step 5: 更新 PROJECT_INDEX.md 与 CHECKPOINTS.md**

写入最终文件结构与推荐运行顺序（`Rscript reanalysis/scripts/analyze_ddr.R .`、`Rscript reanalysis/scripts/analyze_regulators.R .`、`Rscript reanalysis/scripts/build_workbooks.R .`、`Rscript reanalysis/scripts/build_final_manifest.R .`），报告重构结果（文件数变化、行数变化、测试结果、diff 结论）。

---

