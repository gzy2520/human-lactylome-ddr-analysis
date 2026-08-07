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

