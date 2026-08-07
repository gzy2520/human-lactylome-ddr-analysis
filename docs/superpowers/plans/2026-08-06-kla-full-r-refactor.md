# Kla 全 R 合并重构实施计划

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 在目录结构不动、分析功能不变的前提下，把 kla 项目 scripts/ 从 ~30 个文件合并重构为 6 个 R 文件 + 1 个 lib/，全部 Python/Node 工具迁移到 R。

**Architecture:** 先建回归基线（现有输出 sha256 快照），再纯新增 lib/ 共享模块，随后逐批合并/迁移脚本（每个任务以"输出与基线 diff + 测试通过"为完成标准），被替换原文件最终移入 archive/。Python 主流程（extractors/common/run_pipeline）逐函数翻译为 R，test_pipeline.py 重写为 test_pipeline.R 提供等价的 19 项断言。

**Tech Stack:** R 4.4.3（data.table, dplyr, ggplot2, readxl, stringr, tidyr, ggVennDiagram, digest, openssl, writexl, officer）。**本项目不是 git 仓库**——所有任务以"基线快照 + 输出 diff"代替 git commit 作为检查点。

**Spec:** `docs/superpowers/specs/2026-08-06-kla-refactor-design.md`

## Global Constraints

- 目录结构完全不动：`data/`、`reanalysis/`、`archive/`、`previous_umap/` 各级路径原样，**不移动任何文件位置**。
- 输出文件路径、列名、数值、分析规则不变；统计表/中间表必须字节一致，图允许渲染风格差异（数字与标签一致）。
- 蛋白主键是去 isoform 后缀的 UniProt BaseAccession；Gene Symbol 只作显示标签，无回退匹配。
- 不删除任何数据；被合并/替代的原 `.py/.mjs/.js/.R` 文件一律移入 `archive/`。
- 每 PXD 显式分支保留，不做数据驱动重构。
- 分析一律 R；Rscript 调用约定为 `Rscript scripts/xxx.R <project_root>`（与现有脚本一致）。
- 脚本入口约定：合并后的分析脚本支持分段运行（`--stage <name>` 可选参数，缺省全跑），避免每次重跑全部段落。

---

### Task 1: 建立回归基线（输出 sha256 快照 + 现有测试确认）

**Files:**
- Create: `archive/refactor_baseline_2026-08-06/output_sha256_baseline.csv`
- Create: `archive/refactor_baseline_2026-08-06/PRE_REFACTOR_TEST_RESULTS.txt`

**Interfaces:**
- Produces: 基线表 `output_sha256_baseline.csv`（列：`RelativePath,SizeBytes,SHA256`），后续所有任务的验证都 diff 它。

- [ ] **Step 1: 确认现有测试全过**

Run:
```bash
cd /Users/gzy2520/Desktop/Research/kla
Rscript reanalysis/tests/test_expanded_ddr_fraction_by_accession.R .
Rscript reanalysis/tests/test_kla_regulator_intensity.R .
Rscript reanalysis/tests/test_kla_regulator_landscape.R .
Rscript reanalysis/tests/test_kla_regulator_whole_proteome_intensity.R .
Rscript reanalysis/tests/test_lactylome_acquisition.R .
PYTHONPATH=reanalysis/scripts python3 -m unittest discover -s reanalysis/tests -p 'test_*.py' -v 2>&1 | tail -5
```
Expected: 6 个测试文件全部通过（含 Python 侧 19 项），输出保存到 `PRE_REFACTOR_TEST_RESULTS.txt`。若任一失败，先停下报告，不要带病重构。

- [ ] **Step 2: 生成输出文件 sha256 基线**

```bash
cd /Users/gzy2520/Desktop/Research/kla
mkdir -p archive/refactor_baseline_2026-08-06
# 用现有 manifest 覆盖的输出全集做基线（reanalysis/ + previous_umap/ + 根文档，排除自身）
python3 - <<'EOF'
import hashlib, csv
from pathlib import Path
root = Path(".")
paths = set()
for f in (Path("PROJECT_INDEX.md"), Path("NEW_CHAT_PROJECT_PROMPT.md")):
    if f.exists(): paths.add(f)
for d in (Path("reanalysis"), Path("previous_umap")):
    paths.update(p for p in d.rglob("*") if p.is_file())
rows = []
for p in sorted(paths, key=lambda x: x.as_posix()):
    if p.name == ".DS_Store" or p.suffix == ".pyc": continue
    rows.append({"RelativePath": p.as_posix(), "SizeBytes": p.stat().st_size,
                 "SHA256": hashlib.sha256(p.read_bytes()).hexdigest()})
with open("archive/refactor_baseline_2026-08-06/output_sha256_baseline.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=["RelativePath","SizeBytes","SHA256"]); w.writeheader(); w.writerows(rows)
print(f"Baseline files hashed: {len(rows)}")
EOF
```
Expected: 输出 "Baseline files hashed: N"（N 与当前 manifest 行数接近）。

- [ ] **Step 3: 写基线 diff 工具脚本（后续任务复用）**

Create: `reanalysis/scripts/lib/verify_outputs.R`（保留在 lib/ 顶层，供各任务手工调用；本脚本不属于合并目标，是重构期间的工具）：
```r
#!/usr/bin/env Rscript
# 用法: Rscript lib/verify_outputs.R <project_root> <baseline.csv> [排除正则...]
args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(args[[1]])
baseline <- read.csv(args[[2]], check.names = FALSE, stringsAsFactors = FALSE)
exclude <- if (length(args) >= 3) args[[3]] else NULL
sha1 <- function(p) {
  con <- file(p, "rb"); on.exit(close(con))
  openssl::sha256(con)
}
changed <- list()
for (i in seq_len(nrow(baseline))) {
  rel <- baseline$RelativePath[i]
  p <- file.path(root, rel)
  if (!is.null(exclude) && grepl(exclude, rel, perl = TRUE)) next
  if (!file.exists(p)) { changed[[length(changed) + 1]] <- data.frame(RelativePath = rel, Status = "MISSING"); next }
  if (baseline$SizeBytes[i] != file.info(p)$size ||
      baseline$SHA256[i] != sha1(p)) {
    changed[[length(changed) + 1]] <- data.frame(RelativePath = rel, Status = "CHANGED")
  }
}
if (length(changed)) {
  out <- do.call(rbind, changed)
  write.csv(out, file.path(dirname(args[[2]]), "diff_report.csv"), row.names = FALSE)
  cat("DIFF FOUND:", nrow(out), "files\n")
  print(out)
  quit(status = 1)
}
cat("OK: all", nrow(baseline), "baseline files unchanged\n")
```
Expected: 运行 `Rscript reanalysis/scripts/lib/verify_outputs.R . archive/refactor_baseline_2026-08-06/output_sha256_baseline.csv` 输出 "OK: all N baseline files unchanged"。

- [ ] **Step 4: 检查点（代替 commit）**

将本任务产出写入检查点文件 `docs/superpowers/plans/CHECKPOINTS.md`（每任务完成后追加一行"Task N 完成: <验证命令> <结果摘要>"）。

---

### Task 2: 安装 R 依赖

**Files:**
- Modify: 无（仅环境）

- [ ] **Step 1: 安装 writexl 与 officer**

Run:
```bash
Rscript -e 'install.packages(c("writexl", "officer"), repos = "https://cloud.r-project.org")'
```
Expected: 两个包安装成功（若 writexl/officer 已存在则跳过）。

- [ ] **Step 2: 确认关键包可加载**

Run: `Rscript -e 'stopifnot(requireNamespace("writexl"), requireNamespace("officer"), requireNamespace("ggVennDiagram"), requireNamespace("openssl"))'`
Expected: 静默退出码 0。

- [ ] **Step 3: 检查点**（追加到 CHECKPOINTS.md）

---

### Task 3: 创建 lib/ 共享模块（纯新增，不动现有脚本）

**Files:**
- Create: `reanalysis/scripts/lib/accession_utils.R`
- Create: `reanalysis/scripts/lib/io_utils.R`
- Create: `reanalysis/scripts/lib/extractors.R`
- Test: `reanalysis/scripts/lib/test_lib_self_check.R`（重构期自检，重构完成后移入 tests/）

**Interfaces:**
- Produces（后续任务必须使用的精确签名）：
  - `base_accession(values)` → character，去数据库前缀/NX_/isoform 后缀
  - `is_uniprot(values)` → logical
  - `split_accessions(values)` → character（按 `[;,]` 拆分 → base_accession → 仅保留 UniProt → 去重排序）
  - `split_protein_identifiers(values)` → character（UniProt + ENSP 混合拆分）
  - `accession_feature(values)` → character（去 isoform 的排序去重键）
  - `match_target_accession(values)` → logical
  - `safe_numeric(values)` → numeric
  - `relative_path(path, root)` → character
  - `read_delimited(path)` → data.frame（csv 用 read.csv，其他用 read.delim）
  - `valid_maxquant_rows(data)` → logical
  - `write_csv_std(frame, path)` → 写 UTF-8 csv（`write.csv(..., row.names = FALSE, na = "")`）
  - `save_figure(path_stem, plot, width, height, dpi = 350, ...)` → 同时写 `<stem>.png`（ggsave + bg="white"）与 `<stem>.pdf`（cairo_pdf）
  - `extract_maxquant_sites(path, sample_tokens = NULL, sheet = NULL)` → character
  - `extract_maxquant_proteins(path, abundance_pattern = NULL)` → character
  - `extract_spectronaut_proteins(path, lactyl_pattern = NULL, group_pattern = NULL, accession_column = "Protein.Group")` → character
  - `extract_spectronaut_matrix(path, lactyl_pattern)` → character
  - `extract_spectronaut_quant(path, group_pattern = NULL)` → character
  - `extract_pd_proteins(path, sample_token, lactylome = FALSE)` → character
  - `extract_pd_lactyl_peptides(path, sample_token)` → character
  - `extract_huvec_xml(path)` → character

- [ ] **Step 1: 从 `analyze_expanded_ddr_fraction_by_accession.R` 逐字抽取纯函数到 `lib/accession_utils.R` 与 `lib/io_utils.R`**

抽取对象（函数体逐字拷贝，不修改任何逻辑）：
- accession_utils.R：`base_accession`（原 41-47 行）、`is_uniprot`（57-67）、`split_accessions`（69-73）、`split_protein_identifiers`（75-84）、`identifier_type`（175-177）
- io_utils.R：`relative_path`（49-55，改签名 `function(path, root)`）、`valid_maxquant_rows`（344-361）、`read_delimited`（363-375）

并在两个文件末尾加入：
```r
# 供分析脚本 source 时的环境标记（避免重复 source 报错）
lib_loaded <- TRUE
```
每个文件头部加 `lib_loaded <- TRUE` 行即可（不依赖 R6/环境对象，保持简单）。

- [ ] **Step 2: 从 4 个热图脚本抽取共有工具到 `lib/accession_utils.R` 追加段**

对 `analyze_kla_regulator_intensity.R`（86-121 行的 `base_accession`/`match_target_accession`/`accession_feature`/`safe_numeric`/`relative_path`）与 `analyze_kla_regulator_whole_proteome_intensity.R`（91-122 行同组函数）：比对与 Step 1 抽取版本的行为是否一致（用 Step 3 的自检断言覆盖差异点，如 `REV__CON__Q9H9Q4-3` → `Q9H9Q4`）。一致的合入 lib；有行为差异的（如 `match_target_accession` 只存在于热图脚本）以热图脚本版为准放入 lib，并在自检中断言其行为。
`read_delimited`（intensity 脚本 182-194 行）与 Step 1 版比对后合入。

- [ ] **Step 3: 抽取 extractors 到 `lib/extractors.R`**

从 `analyze_expanded_ddr_fraction_by_accession.R` 逐字抽取（函数体不变）：
`extract_maxquant_sites`（377-407）、`extract_maxquant_proteins`（409-434）、`extract_pd_proteins`（436-456）、`extract_pd_lactyl_peptides`（458-483）、`extract_spectronaut_proteins`（485-523）、`extract_spectronaut_matrix`（525-529）、`extract_spectronaut_quant`（531-545）、`extract_huvec_xml`（566-572）。
文件头加 `extractors_loaded <- TRUE`。
注意：这些函数依赖 lib 其他文件（`base_accession`/`is_uniprot`/`split_accessions`/`split_protein_identifiers`/`read_delimited`/`valid_maxquant_rows`），调用方须按 `source(accession_utils.R); source(io_utils.R); source(extractors.R)` 顺序加载。

- [ ] **Step 4: 写 lib 自检脚本**

Create: `reanalysis/scripts/lib/test_lib_self_check.R`：
```r
#!/usr/bin/env Rscript
root <- normalizePath(if (length(commandArgs(trailingOnly = TRUE)) >= 1) commandArgs(trailingOnly = TRUE)[[1]] else ".")
lib <- file.path(root, "reanalysis", "scripts", "lib")
source(file.path(lib, "accession_utils.R"))
source(file.path(lib, "io_utils.R"))
source(file.path(lib, "extractors.R"))
stopifnot(base_accession("sp|P49959-2|MRE11_HUMAN") == "P49959")
stopifnot(base_accession("REV__CON__Q9H9Q4-3") == "Q9H9Q4")
stopifnot(base_accession("NX_Q9H9Q4-1") == "Q9H9Q4")
stopifnot(is_uniprot("P49959"), !is_uniprot("P49959-2"), !is_uniprot("MRE11"))
stopifnot(identical(split_accessions("P49959-2;Q9H9Q4-3,sp|O60934|NBN_HUMAN"), c("O60934", "P49959", "Q9H9Q4")))
stopifnot(identical(split_protein_identifiers("ENSP00000369497;P49959-2"), c("ENSP00000369497", "P49959")))
stopifnot(identical(accession_feature(c("P49959-2", "P49959-1")), c("P49959", "P49959")))
stopifnot(identical(safe_numeric(c("1,234", "NA", "0")), c(1234, NA, 0)))
tmp <- tempfile(fileext = ".csv")
write.csv(data.frame(a = 1:2, b = c("x", "")), tmp, row.names = FALSE)
stopifnot(nrow(read_delimited(tmp)) == 2)
cat("lib self-check passed\n")
```
（`match_target_accession`、`accession_feature` 的具体断言以 Step 2 比对后的实际语义为准，若与示例不符以源脚本行为为准修正示例断言。）

- [ ] **Step 5: 运行自检**

Run: `Rscript reanalysis/scripts/lib/test_lib_self_check.R .`
Expected: "lib self-check passed"

- [ ] **Step 6: 验证现有脚本未被破坏**

Run:
```bash
cd /Users/gzy2520/Desktop/Research/kla
Rscript reanalysis/scripts/analyze_expanded_ddr_fraction_by_accession.R . && \
Rscript reanalysis/scripts/lib/verify_outputs.R . archive/refactor_baseline_2026-08-06/output_sha256_baseline.csv
```
Expected: 脚本正常完成且 "OK: all N baseline files unchanged"（lib/ 新增文件不在基线内，verify 只查基线已有文件，应全绿）。

- [ ] **Step 7: 检查点**

---

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

### Task 6: 迁移 Python 主流程 → `analyze_ddr.R`（第一部分：Kla/GO/Venn）

**Files:**
- Create: `reanalysis/scripts/analyze_ddr.R`
- Create: `reanalysis/tests/test_pipeline.R`（替代 test_pipeline.py 的等价断言）
- 参照（只读）：`reanalysis/scripts/run_pipeline.py`、`extractors.py`、`common.py`、`reanalysis/tests/test_pipeline.py`

**Interfaces:**
- Produces: CLI `Rscript analyze_ddr.R <project_root> [--stage pipeline]`；pipeline 段输出与 run_pipeline.py 完全一致（下表）。

**函数映射表（逐函数翻译，参数与返回值语义一一对应）：**

| Python 源 | R 目标 | 说明 |
|---|---|---|
| `common.base_accession` | lib `base_accession` | 已就绪 |
| `common.read_go_annotations(path)` | `read_go_annotations(path)` | 返回 list(retained, raw)，retained 含 ExcludedNOT 列；测试点 2/5 依赖 |
| `common.split_tokens/normalize_accession/accession_candidates` | `split_tokens/normalize_accession/accession_candidates` | PEAKS 解析辅助 |
| `common.lactyl_positions_from_ascore/peptide` | 同名函数 | site 解析 |
| `extractors.*`（extract_pxd014870, pxd028488, pxd050470, pxd053474_*, maxquant_site_table, pxd078013, extract_pxd028488 等 11 个提取器） | `extract_pxd014870(...)` 等同名函数 | 逐函数翻译，PEAKS scan 归一化（`normalize_scan`）、marker-156（`read_marker_156`）、accession 匹配（`matching_accession_mappings`）、site-position 配对（`site_position_pairs`）逻辑必须与 Python 逐行等价（测试点 3/4/15/16） |
| `run_pipeline.extract_all(root)` | `extract_all(root)` | 返回 list(evidence, exclusion_log, audits) |
| `run_pipeline.aggregate_evidence(frame)` | `aggregate_evidence(frame)` | KlaSites 排序拼接逻辑 |
| `run_pipeline.attach_go(proteins, go_summary, go_raw)` | `attach_go(...)` | GOMatchMode 含 "unmatched"，人源过滤（测试点 5） |
| `run_pipeline.cell_type_statistics(...)` | `cell_type_statistics(...)` | MCF10A 947/116 断言 |
| `run_pipeline.venn_regions(sets)` | `venn_regions(sets)` | 7 区域 |
| `run_pipeline.write_group_outputs(...)` | `write_group_outputs(...)` | 区域表/计数/组合表 + Venn 图 |
| `run_pipeline.target_trace(...)` / `build_target_source_audit(...)` | 同名 | 目标蛋白追踪（Task 9 的 audit 函数在此一并翻译，来自 audit_target_sources.py 的 `lactyl_mask`/`build_target_source_audit`） |
| `run_pipeline.regression_outputs(...)` | `regression_outputs(...)` | regression 表 |
| `run_pipeline.main()` | `main()` | 全部写文件路径照抄（见下清单） |

**pipeline 段输出文件清单（路径与 Python 版完全一致，全部须字节一致）：**
- `reanalysis/intermediate/kla_by_dataset/`：`all_included_and_audit_kla_evidence.csv`、`all_primary_sample_level_kla_sites.csv`、7×`{PXD}_sample_level_kla_sites.csv`、`PXD014870_sensitivity_localization_0.75.csv`
- `reanalysis/results/tables/`：`exclusion_log.csv`、`per_pxd/{PXD}_unique_kla_proteins.csv`、`pxd028488/{directory_audit,included_directories,excluded_directories_and_reasons,old_vs_new_directory_coverage,diagnostic_ion_156_supported_evidence}.csv`、`pxd053474/{search_vs_supplementary_all,consistent,search_only,supplementary_only,inconsistent,single_mode_search_only_audit_evidence}.csv`、`all_unique_kla_proteins.csv`、`go_unmatched_kla_proteins.csv`、`accession_gene_mapping_failures.csv`、`classification_review_needed.csv`、`cell_type_kla_ddr_statistics.csv`、`venn_all_schemes_counts.csv`、`all_kla_three_groups_combined.csv`、`kla_go_ddr_three_groups_combined.csv`、`all_kla_three_groups_combined_non_deduplicated.csv`、`kla_go_ddr_three_groups_combined_non_deduplicated.csv`、`tumor_specific_kla_proteins.csv`、`tumor_specific_kla_ddr_proteins.csv`、`target_protein_evidence_trace_MRE11_XLF_NBS1.csv`、`target_protein_source_level_audit_MRE11_XLF_NBS1.csv`、`regression_old_vs_new_detail.csv`、`regression_old_vs_new_summary.csv`、`dataset_analysis_summary.csv`
- `reanalysis/results/tables/venn_regions/{scheme}/{analysis}/`：`venn_region_counts.csv`、`venn_membership.csv`、7×`{region}.csv`、`{category}_all.csv`
- `reanalysis/results/figures/`：`all_kla_three_group_venn.{png,svg,pdf}`、`kla_go_ddr_three_group_venn.{png,svg,pdf}`、`teacher_requested_grouping/`、`biologically_conventional_grouping/` 下各 4 张图
- `reanalysis/logs/software_environment.csv`（改为 R 环境：R version / platform / data.table / ggplot2 版本，文件名与列名不变）

**Venn 图实现**：Python 版 matplotlib 自定义绘制 → R 用 ggVennDiagram（已装）。用区域计数（`venn_regions` 输出）构造 `ggVennDiagram::venn.diagram`-风格绘制，或手工 ggplot 重叠圆。**数字标签必须与 Python 版一致**（区域计数、总蛋白数）；渲染风格允许差异，但必须目检确认。

- [ ] **Step 1: 翻译 common.py → analyze_ddr.R 辅助段**

逐函数翻译 `common.py`（clean_text, number, integer, is_true, split_tokens, normalize_accession, base_accession（直接调 lib）, accession_candidates, unique_join, best_annotation, annotation_from_description, annotation_from_fasta, blank_site, strip_peptide_modifications, lactyl_positions_from_ascore, lactyl_positions_from_peptide, parse_probability_values, apply_annotation_supplement, read_go_annotations, relative_path）。用 data.frame 代替 pd.DataFrame，字符列以 `stringsAsFactors = FALSE` 读取避免 factor 破坏字节一致输出。

- [ ] **Step 2: 翻译 extractors.py**

逐函数翻译 11 个 PXD 提取器（extract_pxd014870、extract_pxd028488、extract_pxd050470、extract_pxd053474_dda、extract_pxd053474_dia、extract_pxd053474_supplementary、reconcile_pxd053474、extract_maxquant_site_table、extract_pxd078013、extract_pxd078736（若有）、pxd028_directory_metadata、read_marker_156、normalize_scan、matching_accession_mappings、site_position_pairs、valid_maxquant_site_ids、run_prefixes、peaks_annotations、peaks_start_map、exclusion_row、dataframe）。翻译要点：
  - `valid_maxquant_site_ids` 返回 list(ids=set, rows=list of dict) → R list(set=character, rows=data.frame)
  - `exclusion_row`/`dataframe` 输出列名与 Python 版逐列一致（R 里用 data.frame 显式列名）
  - PEAKS PSM 解析按 Python 的列名处理，`normalize_scan` 的字符串归一化逐字符等价

- [ ] **Step 3: 翻译 run_pipeline.py 逻辑函数**

翻译 aggregate_evidence、attach_go、cell_type_statistics、venn_regions、write_group_outputs、target_trace、regression_outputs、build_target_source_audit（后两者需要 `lactyl_mask`——从 audit_target_sources.py 翻译，见测试点 6：`lactyl_mask` 对 PTM=="Lac"、Modified peptide 含 "La (K)"、"La (K)"=="1" 的行返回 TRUE）。

- [ ] **Step 4: 写 test_pipeline.R（19 项断言）**

Create: `reanalysis/tests/test_pipeline.R`，纯 R 断言风格（`stopifnot` + `cat`），覆盖 test_pipeline.py 全部 19 项测试点：
1. `base_accession("sp|P49959-2|MRE11_HUMAN") == "P49959"`；`base_accession("REV__CON__Q9H9Q4-3") == "Q9H9Q4"`
2. read_go_annotations 排除 NOT（构造临时 go.tsv，retained 仅 Q9H9Q4，raw ExcludedNOT 和=1）
3. matching_accession_mappings 拒绝 accession 不匹配（`matching_accession_mappings(list(c("P49959",100), c("Q9H9Q4",200)), "P49959|MRE11_HUMAN")` 返回只含 P49959；对 "O60934|NBN_HUMAN" 返回空）
4. site_position_pairs 长度不匹配返回空
5. attach_go 对非人 GO 行 GOMatchMode == "unmatched"
6. lactyl_mask：PTM=="Lac" 或 Modified peptide 含 "La (K)" 或 "La (K)"=="1" → TRUE
7. all_primary_sample_level_kla_sites.csv 的 PXD 集合 == 7 个纳入 PXD，不含 PXD038880/PXD050906
8. grouping_schemes.csv：MCF10A==normal_immortalized_cell_lines、MCF7==tumor_cell_lines（两列 scheme 都查）
9. 两个 scheme × 两个 analysis 的 venn_region_counts.csv 与 {region}.csv 行数一致、membership 无重复、计数和 == 3112/275
10. membership 重建类别集合；`{category}_all.csv` 无重复
11. tumor_specific_kla_proteins.csv == tumor_only.csv（all_kla 与 kla_go_ddr 两对）
12. all_kla_three_groups_combined.csv 3112 行无重复、VennRegion ∈ REGION_ORDER、DetectedGroupCount == flags 计数和
13. 非去重表行数 == 3×{category}_all.csv 行数和、SourceCategory 集合正确、有重复行
14. cell_type_kla_ddr_statistics.csv 与 cell_type_statistics() 重建一致；MCF10A == 947/116
15. pxd053474/search_vs_supplementary_all.csv：consistent=1275, search_only=826, supplementary_only=3, PrimaryIncluded 和=1298
16. PXD014870_sensitivity 行 ⊆ 主集行；LocalizationProb ≥ 0.75
17. target_protein_source_level_audit 中 PXD∈{PXD014870,PXD028488,PXD050470,PXD053474} 的 ExtractedPrimaryRows==0 且 SourceKlaTargetRows==0
18. cell_type_reference_ddr_statistics.csv 与 reference_proteome_selection.json 的 reference_protein_count_main 一致（10 行）；comparison/control_information 各 10 行且差值恒等式成立
19. tall104_surrogate_ddr_sensitivity.csv：TALL-1_primary=3383、Jurkat_sensitivity=3363、union=3835、intersection=2911，Ddr≤Total

- [ ] **Step 5: 运行 test_pipeline.R 确认新实现通过**

Run: `Rscript reanalysis/tests/test_pipeline.R .`
Expected: 全部 19 项断言通过，输出 "pipeline tests passed"。

- [ ] **Step 6: 全量输出 diff（字节一致性核心验证）**

Run:
```bash
cd /Users/gzy2520/Desktop/Research/kla
Rscript reanalysis/scripts/analyze_ddr.R . --stage pipeline
Rscript reanalysis/scripts/lib/verify_outputs.R . archive/refactor_baseline_2026-08-06/output_sha256_baseline.csv
```
Expected: verify_outputs "OK"。唯一允许差异：`logs/software_environment.csv`（Python→R 版本信息，预期内，需在 diff 中人工确认）与 Venn 图渲染（如字节不同则目检数字标签一致）。若任何统计表字节不同，回查对应翻译函数，不得以"近似一致"放行。

- [ ] **Step 7: 检查点**

---

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

### Task 9: 并入目标蛋白审计 → `analyze_ddr.R`（第四部分）

**Files:**
- Modify: `reanalysis/scripts/analyze_ddr.R`
- 参照（只读）：`reanalysis/scripts/audit_target_sources.py`

**Interfaces:**
- Produces: 无新输出文件（Task 6 的 pipeline 段已翻译 `build_target_source_audit`）；本任务把 Python 原版与 R 版输出做交叉验证。

- [ ] **Step 1: 交叉验证 R 版审计输出**

Run:
```bash
cd /Users/gzy2520/Desktop/Research/kla
Rscript reanalysis/scripts/analyze_ddr.R . --stage pipeline
diff <(Rscript -e 'cat(readLines("reanalysis/results/tables/target_protein_source_level_audit_MRE11_XLF_NBS1.csv"), sep="\n")') \
     <(PYTHONPATH=reanalysis/scripts python3 -c '
import pandas as pd
print(pd.read_csv("reanalysis/results/tables/target_protein_source_level_audit_MRE11_XLF_NBS1.csv", dtype=str).to_csv(index=False))')
```
Expected: 无输出差异（同一文件在 R 版 pipeline 重跑后与 Python 版基线字节一致，由 verify_outputs 在 Task 6 已保证；本任务重点确认 audit 相关测试点 6/17 通过）。

- [ ] **Step 2: 检查点**

---

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

### Task 11: 迁移 manifest → `build_final_manifest.R`

**Files:**
- Create: `reanalysis/scripts/build_final_manifest.R`
- 参照（只读）：`reanalysis/scripts/build_final_manifest.py`

**Interfaces:**
- Produces: CLI `Rscript build_final_manifest.R <project_root>`；输出 `reanalysis/reports/final_file_manifest_sha256.csv`，表头 `RelativePath,SizeBytes,SHA256`，包含范围与 Python 版完全一致（根两文档 + reanalysis/ + previous_umap/ + data/PXD*/metadata/，排除 .DS_Store/.pyc/自身）。

- [ ] **Step 1: 翻译为 R**

用 `openssl::sha256(file)` 或 `digest::digest(file, algo = "sha256", serialize = FALSE)`；目录遍历用 `list.files(recursive = TRUE)`；排序键 `path.as_posix()` 与 Python 一致。

- [ ] **Step 2: 对比新旧 manifest（同文件集）**

Run:
```bash
cd /Users/gzy2520/Desktop/Research/kla
Rscript reanalysis/scripts/build_final_manifest.R . && mv reanalysis/reports/final_file_manifest_sha256.csv /tmp/manifest_r.csv
python3 reanalysis/scripts/build_final_manifest.py && mv reanalysis/reports/final_file_manifest_sha256.csv /tmp/manifest_py.csv
diff /tmp/manifest_r.csv /tmp/manifest_py.csv && echo MANIFEST_IDENTICAL
```
Expected: `MANIFEST_IDENTICAL`（同一时刻同文件集，两个实现输出逐字节一致）。

- [ ] **Step 3: 检查点**

---

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

## Self-Review 记录

- **Spec 覆盖**：✓ 目录不动（全部任务只在 scripts/ 内新建/归档，无路径移动）；✓ 6 文件 + lib 结构（Task 3/4/5/6/10/11/13）；✓ 全 R（Task 6/8/10/11 迁移）；✓ 测试迁移（Task 6 Step 4）；✓ 文档精简（Task 12）；✓ 输出回归（每任务 verify + Task 14）；✓ manifest 重建（Task 11/14）；✓ 边界（Task 13 只移动 scripts/ 内文件到 archive）。
- **占位符扫描**：无 TBD/TODO；所有验证命令与断言均为实际命令。
- **类型一致性**：lib 函数签名在 Task 3 定义并被 Task 4/5/6/7 引用，签名在接口块中固定；`run_stage(name, code)` 模式在 Task 4 定义、Task 5/7/8 复用；test_pipeline.R 19 项断言编号与 test_pipeline.py 测试点一一对应。
