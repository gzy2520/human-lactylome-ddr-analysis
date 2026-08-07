### Task 4: 合并调控因子三图 → `analyze_regulators.R`

**Files:**
- Create: `reanalysis/scripts/analyze_regulators.R`
- Modify: 无现有脚本（原三个脚本在 Task 13 归档前保持可用）
- Test: `reanalysis/tests/test_kla_regulator_intensity.R`、`test_kla_regulator_landscape.R`、`test_kla_regulator_whole_proteome_intensity.R`（不改，仍指向原输出路径）

**Interfaces:**
- Produces: CLI `Rscript analyze_regulators.R <project_root> [--stage heatmap|whole_proteome|landscape]`（缺省三段全跑）；输出文件与三个原脚本完全一致。

- [ ] **Step 1: 建立合并骨架**

Create: `reanalysis/scripts/analyze_regulators.R`，头部：
```r
#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(ggplot2)
  library(readxl); library(stringr); library(tidyr)
})
args <- commandArgs(trailingOnly = TRUE)
project_root <- normalizePath(if (length(args)) args[[1]] else ".")
stage <- if (length(args) >= 2) args[[2]] else "all"
lib <- file.path(project_root, "reanalysis", "scripts", "lib")
for (f in c("accession_utils.R", "io_utils.R", "extractors.R"))
  source(file.path(lib, f))
run_stage <- function(name, code) {
  if (stage == "all" || stage == name) { message("[stage] ", name); code }
}
```
（stage 参数放在 project_root 之后，与现有 `Rscript x.R <root>` 调用兼容。）

- [ ] **Step 2: 迁入 landscape 段（`plot_kla_regulator_landscape.R`）**

把原脚本 128 行之后的函数（`split_target_genes`、`empty_evidence`、`extract_maxquant_sites`（本脚本私有变体，签名不同，保留私有名）、`extract_spectronaut`（同样私有）、`add_evidence`、`set_audit`、`make_role_panel`）与 main 执行段包进 `run_stage("landscape", {...})`。脚本私有函数与 lib 同名者（`base_accession` 等）**从原脚本删除**改用 lib 版本（先用自检/测试确认行为一致）。若私有 `extract_maxquant_sites` 与 lib 版签名/行为不同，保留私有定义且不改名（R 中后 source 的遮蔽先 source 的，属预期，注释说明）。

- [ ] **Step 3: 迁入 heatmap 段（`analyze_kla_regulator_intensity.R`）**

同样包进 `run_stage("heatmap", {...})`，删除与 lib 重复的定义，保留私有函数（`add_quant`、`extract_spectronaut_precursor`、`unavailable_reason` 等）。注意原脚本 182-194 行的 `read_delimited` 与 lib 版行为须一致后才删除（Step 3/4 自检已覆盖）。

- [ ] **Step 4: 迁入 whole_proteome 段（`analyze_kla_regulator_whole_proteome_intensity.R`）**

包进 `run_stage("whole_proteome", {...})`，同上处理私有函数（`map_reviewed_symbol_features`、`update_audit`、`add_total_quant`、`add_maxquant_proteome`、`add_pd_proteome`、`add_spectronaut_report`、`add_spectronaut_standard_report`、`add_spectronaut_matrix`、`add_peaks_proteins`）。该脚本私有 `add_spectronaut_matrix` 与 lib `extract_spectronaut_matrix` 不同，保留私有。

- [ ] **Step 5: 全跑合并脚本并做输出 diff**

Run:
```bash
cd /Users/gzy2520/Desktop/Research/kla
Rscript reanalysis/scripts/analyze_regulators.R . 
Rscript reanalysis/scripts/lib/verify_outputs.R . archive/refactor_baseline_2026-08-06/output_sha256_baseline.csv
Rscript reanalysis/tests/test_kla_regulator_intensity.R .
Rscript reanalysis/tests/test_kla_regulator_landscape.R .
Rscript reanalysis/tests/test_kla_regulator_whole_proteome_intensity.R .
```
Expected: verify_outputs 输出 "OK"（所有图与表字节不变）；3 个 R 测试全部通过。若某图字节变化，检查是否是 ggVennDiagram/字体差异导致——图只允许"数字与标签一致"的渲染差异，须人工目检确认。

- [ ] **Step 6: 检查点**

---

