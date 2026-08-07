### Task 13: 归档被合并/迁移的原文件

**Files:**
- Move: `reanalysis/scripts/*.py`、`*.mjs`、`*.js`（build_final_manifest.py 等全部）与已合并的 `.R` 原脚本 → `archive/2026-08-06_full_r_refactor/scripts/`
- Move: `reanalysis/tests/test_pipeline.py` → 同 archive
- Clean: `reanalysis/scripts/__pycache__/`、`reanalysis/scripts/reanalysis/`（仅含 .DS_Store 的空目录）

**Interfaces:**
- Produces: `reanalysis/scripts/` 最终只含：`lib/`（accession_utils.R、io_utils.R、extractors.R、verify_outputs.R）+ `acquire_data.R` + `analyze_ddr.R` + `analyze_regulators.R` + `build_workbooks.R` + `build_final_manifest.R`。

- [ ] **Step 1: 核对脚本引用**

`grep -rn "run_pipeline\|extractors\|common\|\.mjs\|\.js\b\|\.py" --include="*.md" --include="*.R" --include="*.csv" /Users/gzy2520/Desktop/Research/kla/reanalysis /Users/gzy2520/Desktop/Research/kla/PROJECT_INDEX.md /Users/gzy2520/Desktop/Research/kla/NEW_CHAT_PROJECT_PROMPT.md | grep -v archive | grep -v "\.pyc"` —— 找到所有残留引用，逐处更新为 R 版路径/命令（Task 12 已处理文档，此处处理漏网引用）。

- [ ] **Step 2: 移动文件**

```bash
cd /Users/gzy2520/Desktop/Research/kla
mkdir -p archive/2026-08-06_full_r_refactor/scripts
mv reanalysis/scripts/*.py reanalysis/scripts/*.mjs reanalysis/scripts/*.js archive/2026-08-06_full_r_refactor/scripts/ 2>/dev/null
mv reanalysis/scripts/analyze_expanded_ddr_fraction_by_accession.R \
   reanalysis/scripts/analyze_kla_regulator_intensity.R \
   reanalysis/scripts/analyze_kla_regulator_whole_proteome_intensity.R \
   reanalysis/scripts/plot_kla_regulator_landscape.R \
   reanalysis/scripts/build_lactylome_reference_pairing.R \
   reanalysis/scripts/build_human_lactylome_inventory.R \
   reanalysis/scripts/build_lactylome_acquisition_manifests.R \
   reanalysis/scripts/build_ensembl_uniprot_mapping.R \
   reanalysis/scripts/build_healthy_special_reference_manifest.R \
   reanalysis/scripts/download_healthy_tissue_references.R \
   reanalysis/scripts/download_lactylome_pair_files.R \
   reanalysis/scripts/extract_lactylome_pair_archives.R \
   reanalysis/scripts/probe_lactylome_pair_files.R \
   reanalysis/scripts/register_additional_lactylome_pair_files.R \
   reanalysis/scripts/summarize_acquired_lactylome_data.R \
   reanalysis/scripts/build_venn_combined_workbook.mjs archive/2026-08-06_full_r_refactor/scripts/ 2>/dev/null
mv reanalysis/tests/test_pipeline.py archive/2026-08-06_full_r_refactor/scripts/
rm -rf reanalysis/scripts/__pycache__ reanalysis/scripts/reanalysis
```
注意：`.mjs/.js` 与部分 `.R` 已在上一行 glob 覆盖，具体以实际文件为准；移动前用 `ls reanalysis/scripts/` 核对清单，**确认 build_workbooks.R 需要的 xlsx 写入逻辑已全部在 R 版内**。

- [ ] **Step 3: 验证 scripts/ 终态**

Run: `ls /Users/gzy2520/Desktop/Research/kla/reanalysis/scripts/`
Expected: 只有 `lib/`、`acquire_data.R`、`analyze_ddr.R`、`analyze_regulators.R`、`build_workbooks.R`、`build_final_manifest.R`（+ .DS_Store 可忽略）。

- [ ] **Step 4: 检查点**

---

