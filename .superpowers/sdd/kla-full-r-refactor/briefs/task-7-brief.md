### Task 7: 并入 DDR 占比分析 → `analyze_ddr.R`（第二部分）

**Files:**
- Modify: `reanalysis/scripts/analyze_ddr.R`（追加 `--stage expanded` 段）
- Test: `reanalysis/tests/test_expanded_ddr_fraction_by_accession.R`（不改）

**Interfaces:**
- Produces: `--stage expanded` 输出与 `analyze_expanded_ddr_fraction_by_accession.R` 一致；**同时移除重复输出文件**：`cell_type_kla_vs_reference_ddr_fraction.png/.pdf` 与 `_accession_only` 版本字节相同，只写 `_accession_only` 名（在 Task 12 把同名旧副本移 archive 前，本任务先只写主名；顺序：Task 7 改脚本只写一份 → Task 12 处理旧文件）。

- [ ] **Step 1: 迁入 expanded 段**

把 `analyze_expanded_ddr_fraction_by_accession.R` 的 main 流程（pairing 读取、37 组循环、kla_set/reference_set、统计/审计组装、evidence_audit、绘图）包进 `run_stage("expanded", {...})`。删除与 lib 重复的 470 行工具函数（原 41-87、344-572 行等），改用 lib；保留 `biomart_mapping`、`map_ensembl_proteins`（私有，依赖 `ensembl_mapping_records` 状态，保留原样）、`extract_pxd043880_reference`、`reviewed_uniprot` 处理段、`kla_set`/`reference_set` 全部分支。
绘图段：4 个 `ggsave`（原 1429-1466 行，两组文件名 × png/pdf）改为 1 次 `save_figure(file.path(figure_dir, "cell_type_kla_vs_reference_ddr_fraction_accession_only"), p, ...)`——即只保留 `_accession_only` 命名。

- [ ] **Step 2: 运行并 diff**

Run:
```bash
cd /Users/gzy2520/Desktop/Research/kla
Rscript reanalysis/scripts/analyze_ddr.R . --stage expanded
Rscript reanalysis/scripts/lib/verify_outputs.R . archive/refactor_baseline_2026-08-06/output_sha256_baseline.csv
Rscript reanalysis/tests/test_expanded_ddr_fraction_by_accession.R .
```
Expected: verify_outputs "OK"（注意：本任务起 `cell_type_kla_vs_reference_ddr_fraction.png/.pdf` 不再被重写——它们仍是旧内容，verify 会通过；Task 12 才移除它们。`_accession_only` 两组文件必须字节一致）。

- [ ] **Step 3: 检查点**

---

