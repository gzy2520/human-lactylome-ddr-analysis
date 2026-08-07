### Task 10: 合并交付工具 → `build_workbooks.R`

**Files:**
- Create: `reanalysis/scripts/build_workbooks.R`
- 参照（只读）：`build_cell_type_statistics_workbook.mjs`、`build_hippocampus_review_md.mjs`、`build_reference_ddr_comparison_workbook.mjs`、`build_reference_proteome_selection_workbook.mjs`、`build_venn_combined_workbook.mjs`、`create_bilingual_figure_legends_docx.js`、`build_project_metadata.py`

**Interfaces:**
- Produces: CLI `Rscript build_workbooks.R <project_root> [--stage workbooks|metadata]`；输出文件与原工具一致：`cell_type_kla_ddr_statistics.xlsx`、`lactylome_and_reference_proteome_pairing_zh.xlsx`、`cell_type_kla_vs_reference_ddr_statistics.xlsx`、`cell_type_reference_proteome_selection.xlsx`、`venn 相关 workbook`（以原 .mjs 输出清单为准）、`Kla_Venn_figure_legends_bilingual.docx`、`reanalysis/reports/project_file_inventory.csv` 等 metadata 输出。

- [ ] **Step 1: 核对原 .mjs 的输出文件清单与 sheet 结构**

读 4 个 .mjs + 1 个 .js 的 main 段，列出各自输出路径与 sheet 名/列结构，写入本任务记录（`docs/superpowers/plans/task10_output_spec.md`），作为 R 实现与验证依据。

- [ ] **Step 2: 实现 build_workbooks.R**

xlsx 用 `writexl::write_xlsx(list(sheet1 = df1, ...), path)`（sheet 名与列结构照抄）；docx 用 `officer`（标题/段落/表格与 .js 生成内容一致）；`build_project_metadata.py` 的 inventory csv 逻辑照译（用 data.table 快速扫描文件树）。

- [ ] **Step 3: 运行并验证内容一致**

Run:
```bash
cd /Users/gzy2520/Desktop/Research/kla
Rscript reanalysis/scripts/build_workbooks.R .
Rscript -e '
suppressPackageStartupMessages(library(readxl))
# 对比新旧 xlsx：读回列名与单元格值
for (f in c("reanalysis/results/tables/cell_type_kla_ddr_statistics.xlsx",
            "reanalysis/results/tables/cell_type_kla_vs_reference_ddr_statistics.xlsx")) {
  sheets <- excel_sheets(f)
  for (s in sheets) {
    d <- read_excel(f, sheet = s)
    cat(f, "|", s, "| rows:", nrow(d), "cols:", ncol(d), "\n")
  }
}
'
Rscript reanalysis/scripts/lib/verify_outputs.R . archive/refactor_baseline_2026-08-06/output_sha256_baseline.csv
```
Expected: xlsx 行/列数与原版一致（xlsx 为二进制，字节不必一致，但**内容必须一致**：sheet 名、列名、行数、单元格值抽样对比）；docx 用 `unzip -p` 提取 document.xml 文本对比图注文字一致；metadata csv 字节一致。

- [ ] **Step 4: 检查点**

---

