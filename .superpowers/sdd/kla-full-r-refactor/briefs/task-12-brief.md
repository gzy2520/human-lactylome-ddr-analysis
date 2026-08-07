### Task 12: 文档精简与重复输出清理

**Files:**
- Modify: `PROJECT_INDEX.md`、`reanalysis/README.md`、`NEW_CHAT_PROJECT_PROMPT.md`
- Move to archive: 被替代的 reports（先核对，候选 `reports/METHODS.md`、`reports/CURRENT_LACTYLOME_ACQUISITION_STATUS.md`、`reports/HUMAN_HIPPOCAMPUS_25_DATASET_REVIEW_TABLE.md`、`reports/HUMAN_HIPPOCAMPUS_25_DATASET_SCREENING.md`、`reports/HUMAN_LACTYLOME_DATA_ACQUISITION.md`；核对方法：读被替代者与其继任者，确认继任者包含其独有事实，否则保留）
- Move to archive: `reanalysis/results/figures/cell_type_kla_vs_reference_ddr_fraction.png`、`cell_type_kla_vs_reference_ddr_fraction.pdf`（与 `_accession_only` 字节相同的重复件）

**Global Constraints 提醒:** 不删除数据——全部移入 `archive/2026-08-06_full_r_refactor/docs/` 与 `archive/2026-08-06_full_r_refactor/figures/`，并记录迁移清单。

- [ ] **Step 1: 核对报告新旧关系**

逐个读候选报告与其继任者（如 METHODS.md → DATA_SOURCE_AND_ANALYSIS_ALGORITHM.md；CURRENT_LACTYLOME_ACQUISITION_STATUS.md → LACTYLOME_REFERENCE_PAIRING_STATUS.md），确认继任者覆盖其独有内容后才归档；无法确认的保留原位。归档清单写入 `archive/2026-08-06_full_r_refactor/archive_manifest.csv`。

- [ ] **Step 2: 精简入口三文档**

- PROJECT_INDEX.md：更新"主要脚本"节为合并后的 6 个脚本名；更新推荐运行顺序；删除与 README/启动 prompt 重复的段落（保留：ID 规则、数据范围、分析状态、结果链接、限制）。
- reanalysis/README.md：更新脚本引用与运行命令；保留目录说明与报告指引。
- NEW_CHAT_PROJECT_PROMPT.md：更新脚本名与命令（`python3 run_pipeline.py` → `Rscript analyze_ddr.R .`；`python3 build_final_manifest.py` → `Rscript build_final_manifest.R .` 等）；保留分析原则与结果状态。
- 验证：三个文档互不重复同一命令块。

- [ ] **Step 3: 移动重复图到 archive**

`mkdir -p archive/2026-08-06_full_r_refactor/figures && mv reanalysis/results/figures/cell_type_kla_vs_reference_ddr_fraction.{png,pdf} archive/2026-08-06_full_r_refactor/figures/`

- [ ] **Step 4: 检查点**

---

