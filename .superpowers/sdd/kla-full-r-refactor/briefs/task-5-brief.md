### Task 5: 合并数据获取脚本 → `acquire_data.R`

**Files:**
- Create: `reanalysis/scripts/acquire_data.R`
- Test: `reanalysis/tests/test_lactylome_acquisition.R`（不改）

**Interfaces:**
- Produces: CLI `Rscript acquire_data.R <project_root> [--stage pairing|inventory|manifests|ensembl_map|download|summarize]`；**download 段默认跳过**（避免重复联网下载），仅当显式 `--stage download` 时执行。

- [ ] **Step 1: 合并纯计算段**

把 `build_ensembl_uniprot_mapping.R`（161 行，函数 `connect_mart` + main）、`build_lactylome_reference_pairing.R`（1239 行）、`build_human_lactylome_inventory.R`（496 行）、`build_lactylome_acquisition_manifests.R`（218 行）、`build_healthy_special_reference_manifest.R`（66 行）、`summarize_acquired_lactylome_data.R`（214 行）按 `run_stage()` 模式（复用 Task 4 定义的骨架模式，在本脚本内复制该模式，不依赖 analyze_regulators.R）合并入 `acquire_data.R`。各段内部函数保持原逻辑，删除与 lib 重复的工具函数。

- [ ] **Step 2: 合并 IO 段（默认跳过）**

`download_healthy_tissue_references.R`（157）、`download_lactylome_pair_files.R`（175）、`extract_lactylome_pair_archives.R`（70）、`probe_lactylome_pair_files.R`（104）、`register_additional_lactylome_pair_files.R`（73）并入 `--stage download` 段（默认不执行）。

- [ ] **Step 3: 只重跑计算段并验证**

Run:
```bash
cd /Users/gzy2520/Desktop/Research/kla
Rscript reanalysis/scripts/acquire_data.R .
Rscript reanalysis/scripts/lib/verify_outputs.R . archive/refactor_baseline_2026-08-06/output_sha256_baseline.csv
Rscript reanalysis/tests/test_lactylome_acquisition.R .
```
Expected: verify_outputs "OK"；test_lactylome_acquisition.R 通过（注意该测试可能引用 `reanalysis/config/` 下清单文件，合并后文件路径未变应仍通过）。

- [ ] **Step 4: 检查点**

---

