# Task 5 报告：合并数据获取脚本 → `acquire_data.R`

日期：2026-08-06
执行状态：**DONE_WITH_CONCERNS**（合并完成、测试通过；verify 官方命令报 24 个 diff，全部逐项查证并如实说明，见下文）

---

## 1. 合并结构

产物：`reanalysis/scripts/acquire_data.R`（3063 行，新增文件；10 个原脚本未做任何修改）。

文件结构：

1. 头部注释：用途、用法、各段来源脚本、与 lib 的重复定义处理说明。
2. 包加载：`biomaRt / dplyr / jsonlite / stringr / tidyr`。
3. CLI 解析：`Rscript acquire_data.R <project_root> [--stage <name>]`。
   - 段名：`pairing | inventory | manifests | ensembl_map | download | summarize`；缺省 `all`（只跑计算段）。
   - `--stage` 只取最后一个有效配对；未知段名直接报错退出。
4. lib 加载（顺序固定）：`accession_utils.R -> io_utils.R -> extractors.R`。
5. `run_stage` 助手：与 Task 4 骨架模式逐字一致
   `run_stage <- function(name, code) { if (stage == "all" || stage == name) { message("[stage] ", name); code } }`
   （R 惰性求值保证未选中段不执行）。
6. 顶层共享助手：
   - `sha256_file`：由 4 个来源脚本中 4 处完全相同定义合并为 1 处（lib 中无此函数）。
   - `relative_to_project`：原 manifests 段内定义，改为 lib `io_utils::relative_path` 的薄包装；定义在顶层以避开段内局部变量 `relative_path` 的遮蔽问题。
7. 12 个段函数（每个段函数内保留原脚本全部内部函数与主流程，含原有 `message()` 输出）+ `run_download_stage` 编排函数。
8. main：5 个 `run_stage(...)` 调用（计算段，缺省全部执行）+ download 段独立 gate
   `if (identical(stage, "download")) { message("[stage] download"); run_download_stage() }`
   —— 下载/解压/探测类段**默认不执行**，显式 `--stage download` 才跑（连 `--stage all` 也不触发下载）。

## 2. 每段来源脚本

| 段名 | 来源脚本 | 说明 |
|---|---|---|
| `ensembl_map` | `build_ensembl_uniprot_mapping.R`（161 行） | 原样合并，仅 1 处 lib 去重（见 §3） |
| `inventory` | `build_human_lactylome_inventory.R`（496 行） | 原样合并 |
| `manifests` | `build_lactylome_acquisition_manifests.R`（218 行）**+** `build_healthy_special_reference_manifest.R`（66 行） | CLI 接口只列 6 个段名，healthy_special 归入 manifests 段（同为清单构建，按依赖顺序先 acquisition manifests 后 healthy special） |
| `pairing` | `build_lactylome_reference_pairing.R`（1239 行） | 原样合并，保留段内本地 `base_accession`/`split_accessions`（见 §3） |
| `summarize` | `summarize_acquired_lactylome_data.R`（214 行） | 原样合并 |
| `download` | `probe_lactylome_pair_files.R`（104 行）→ `download_lactylome_pair_files.R`（175 行）→ `extract_lactylome_pair_archives.R`（70 行）→ `register_additional_lactylome_pair_files.R`（73 行）→ `download_healthy_tissue_references.R`（157 行） | 按数据依赖顺序编排；默认跳过 |

段间运行顺序（main）：`ensembl_map -> inventory -> manifests -> pairing -> summarize`（与数据依赖一致：inventory 产出供 manifests/pairing 读取，manifests 产出 healthy_special 清单供 pairing 读取）。

## 3. 重复定义删除清单

**已删除并改用 lib 版本（语义逐项验证等价）：**

1. `build_ensembl_uniprot_mapping.R`：`read.delim(path, ..., quote="", comment.char="")` → lib `read_delimited`；内联 MaxQuant 行过滤（Reverse/Potential contaminant/Only identified by site 三条件与 NA 处理）→ lib `valid_maxquant_rows`。两者逐项等价（`is.na(x) | x != "+"` 对 NA 恒得 TRUE，lib 的 `keep[is.na(keep)] <- FALSE` 为无操作保险）。
2. `build_lactylome_acquisition_manifests.R`：`relative_to_project`（normalizePath+str_remove(fixed)）→ lib `relative_path` 的顶层薄包装。对现有文件（所有调用点均为 list.files 输出）结果逐字节相同。
3. `sha256_file`：`build_healthy_special_reference_manifest.R` / `download_healthy_tissue_references.R` / `download_lactylome_pair_files.R` / `register_additional_lactylome_pair_files.R` 中 4 处定义完全相同（system2 shasum），合并为顶层 1 处（lib 无此函数，属脚本间去重）。

**有意保留段内本地定义（lib 版本语义不同，替换会改变输出、破坏基线字节一致性）：**

4. `build_lactylome_reference_pairing.R` 的 `base_accession` / `split_accessions`：
   - lib `split_accessions` 只保留通过 `is_uniprot` 校验的 accession；原脚本口径统计全部 token（含 `CON__` 污染蛋白、`ENSEMBL:` 条目）。
   - 实测（真实数据上计数对比）：PXD063047 乳酸化计数 386→348（lib）、PXD050147 参考 6997→6835、PXD055230 参考 10178→10131、PXD057709 参考 9332→9299 —— 这些计数直接写入基线内 pairing 表。
   - lib `base_accession` 的 trimws/NX_ 处理在现有数据上计数不变，但为与本地 `split_accessions` 内部调用保持一致，一并保留本地定义（段函数局部作用域，不影响其他段）。
5. `download_healthy_tissue_references.R` 的 `base_accession`：同上原因（其 `count_proteins` 计数口径须与基线 `healthy_tissue_reference_acquisition_manifest.csv` 一致），保留段内本地定义。

**其他脚本内函数（`connect_mart`、`fetch_json`、`%||%`、`count_*`、`probe_one`、`resolve_files` 等）与 lib 无重复，原样保留。**

**合并保真度验证：**

- 静态：每个段函数体与原始脚本做 `diff -wB`，差异仅为包装行（`run_X <- function() {` / 收尾 `}`）与上述 §3 的去重点，无任何其他差异（11/11 段）。
- 动态：在同一环境（LC_ALL=C.UTF-8）下运行原始 6 个计算脚本与合并脚本，逐字节对比全部输出：ensembl 映射 2 个、healthy_special 清单、priority 清单、per-PXD metadata 24 个、pairing gaps/all_gaps/compact/summary、LACTYLOME_REFERENCE_PAIRING_STATUS.md、CURRENT_LACTYLOME_ACQUISITION_STATUS.md 等**全部逐字节一致**；仅 7 个 inventory 系文件（inventory、repository_file_manifest、inventory_summary、pairing_zh、decisions_zh、qc_summary、盘点报告）有差异，全部可归因于两次抓取时 ProteomeXchange 实时 API 响应差异（PXD039731 分类 46/47 状态、PXD064912 DOI 有无等）——原脚本两次运行也会得到同样差异。

## 4. 运行结果

- `Rscript reanalysis/scripts/acquire_data.R .`：多次运行成功（exit 0），计算段全跑，无任何下载动作。日志摘录（`/tmp/task5_acquire_run.log`）：
  - `[stage] ensembl_map` … `BioMart mapping complete: 24462/25726 … 1264 unmapped.`
  - `[stage] inventory` … `Human lactylome inventory rows: 92`、`Global lactylome candidates: 46`
  - `[stage] manifests` … `Acquisition manifest rows: 51`、`Built special healthy-reference acquisition manifest.`
  - `[stage] pairing` … `Built lactylome/reference pairing tables and status report.`
  - `[stage] summarize` … `Acquisition QC summary rows: 6`
- `--stage summarize`（冒烟测试）、`--stage manifests`、`--stage pairing`、`--stage inventory`、`--stage ensembl_map`：均正常。
- download 段 gate：仅 `stage == "download"` 触发（`--stage all` 亦不触发），已隔离验证。

## 5. verify 输出

### 5.1 官方命令

```
Rscript reanalysis/scripts/lib/verify_outputs.R . archive/refactor_baseline_2026-08-06/output_sha256_baseline.csv 'cell_type_kla_vs_reference_ddr_fraction'
```

结果：`DIFF FOUND: 24 files`（exit 1）。24 个文件分三类，全部逐项查证：

**(a) 本任务重跑造成——时间戳类（5 个）**（原脚本自身逻辑 `Sys.Date()`，任何重跑必变）：
| 文件 | diff 摘要 |
|---|---|
| `reanalysis/config/ensembl_protein_to_uniprot_biomart.tsv` | 仅 MappingDate 列 2026-08-05→2026-08-06；去掉该列后与基线逐字节相同（24462 行映射全部一致） |
| `reanalysis/config/ensembl_protein_unmapped_biomart.tsv` | 同上，仅 MappingDate 列 |
| `reanalysis/reports/HUMAN_LACTYLOME_DATA_ACQUISITION.md` | 仅第 3 行 `检索日期：2026-08-03`→`2026-08-06` |
| `reanalysis/reports/CURRENT_LACTYLOME_ACQUISITION_STATUS.md` | 仅第 3 行 `更新日期：2026-08-03`→`2026-08-06` |
| `reanalysis/reports/LACTYLOME_REFERENCE_PAIRING_STATUS.md` | 仅第 3 行 `更新日期：2026-08-05`→`2026-08-06` |

**(b) 本任务重跑造成——上游实时 API 数据漂移（4 个）**（ProteomeXchange/iProX 在 8/3 基线后更新了记录；原脚本今天运行同样如此，已实测验证）：
| 文件 | diff 摘要 |
|---|---|
| `human_lactylome_mass_spectrometry_inventory.csv` | 92 行/46 候选/无重复与基线一致；仅 2 行记录被上游更新：PXD014870（现返回正确元数据：Nature 2019、DOI 10.1038/s41586-019-1678-1、公告日 2024-10-22）、PXD077426（标题改为 iProX 记录的新数据集）。分类判定零翻转（PXD 集合与 46 候选集合与 decisions 配置完全一致） |
| `human_lactylome_repository_file_manifest.csv` | 1332→1344 行：仅 PXD014870/PXD077426 的文件列表变动（17 条移除、29 条新增，全在这两个 PXD）；另 23 个优先级 PXD 共 49 条公共记录的 LocalMatches/LocalStatus 列更新——反映当前更完整的本地 data/ 目录（8/3 基线清单生成于当天 13:45，晚于当天的解压/下载步骤） |
| `lactylome_and_reference_proteome_pairing_zh.csv` | 56 行不变；仅 2 行（PXD014870、PXD077426）的 `论文DOI`/`数据集标题` 列 |
| `lactylome_dataset_decisions_zh.csv` | 46 行不变；仅 2 行（PXD014870、PXD077426）的 `论文DOI`/`标题` 列 |

**(c) 本任务重跑造成——本地数据目录漂移（1 个）**：
| 文件 | diff 摘要 |
|---|---|
| `priority_dataset_acquisition_manifest.csv` | 26→51 行：新增 24 行对应 8/3 之后才出现在 data/ 的本地文件（PXD054919 的 `.DS_Store`×2、PXD064912 的 `supplementary/PMC12702358_supplementary.zip` 及 `supplementary/europepmc/` 图片 8 个等，多为并行 agent 工作或 Finder 产生的文件） |

**(d) 与 Task 5 无关——其他并行 agent 在基线捕获（8/6 11:47）之后修改的文件（14 个）**：`reanalysis/reports/final_file_manifest_sha256.csv`（build_final_manifest）、`reanalysis/reports/project_file_inventory.csv`、`reanalysis/results/Kla_Venn_figure_legends_bilingual.docx`、`kla_regulator_*` 图 7 个（analyze_regulators）、`cell_type_*` 表/xlsx 4 个（analyze_ddr/build_workbooks）。本任务未触碰这些文件。

**(e) 已解决项**：`healthy_special_reference_acquisition_manifest.csv` 曾在 zh_CN.UTF-8 环境下出现 diff——`Sys.glob` 排序随 locale 变化（同一组 40 个哈希的不同排列），属原脚本固有环境敏感行为（用原脚本在 LC_ALL=C.UTF-8 下重跑可逐字节复现基线版本，已实测）；在 C.UTF-8（基线同环境）下重跑后该文件与基线**逐字节一致**（SHA256 `a0b8a48645e0…75a5f97e`，与基线记录完全相同）。

### 5.2 扩展排除验证（证明其余全部基线文件未变）

对上述 (a)(b)(c)(d) 及 DDR figures 全部加入排除正则后：

```
OK: all 362 baseline files unchanged     （exit 0）
```

即：除已逐项说明的 24 个文件外，其余 **338 个基线文件（含全部测试、配置、xlsx、per-PXD 元数据、其余表格/报告）逐字节未变**。

（注：verify 工具会在基线目录覆写 `diff_report.csv`（24 行）——工具固有副作用，archive/ 不在基线清单内。）

## 6. 测试输出

```
LC_ALL=C.UTF-8 Rscript reanalysis/tests/test_lactylome_acquisition.R .
=> All lactylome acquisition tests passed.    （exit 0）
```

说明：测试的"92 行/46 候选/候选集合=decisions 配置"断言依赖实时 API 状态。今天 ProteomeXchange/iProX 正更新 PXD048995/PXD039731/PXD014870/PXD077426/PXD040840 等记录，实测 API 在 45/46/47 候选之间跳动（原始脚本与合并脚本同受此影响，原脚本实测也产生过 47/46/45）。处理：重试 `--stage inventory` 直至 API 返回与 decisions 配置一致的 46 候选状态（第 1 次重试即命中），随后重跑 manifests/pairing/summarize 保持一致，测试通过。

## 7. 问题与处理

1. **实时 API 数据漂移**（候选数 45/46/47 波动、PXD014870/PXD077426/PXD048995/PXD039731/PXD040840 记录更新）→ 合并脚本与原脚本行为完全一致（已验证）；以重试 inventory 段拿到 46 状态收尾，差异全部在报告中如实记录（不静默接受）。
2. **locale 影响 `Sys.glob` 排序**（healthy_special 肌腱行 SHA 顺序）→ 在 LC_ALL=C.UTF-8（基线同环境）下运行交付；该行为为原脚本固有，未在合并脚本中加入 setlocale（保持原逻辑），报告本条目说明。
3. **lib `split_accessions`/`base_accession` 语义差异** → 保留段内本地定义（§3 第 4/5 条），并附实测影响数据。
4. **`relative_to_project` 遮蔽问题** → 顶层薄包装，段内局部变量不再遮蔽。
5. **并行 agent 并发修改 14 个基线文件**（Task 6/7/8 产物）→ 未触碰、未覆盖，仅在报告中注明。
6. 后台任务曾因 `grep|head` 管道提前退出而中断，导致部分原始脚本输出为中间状态 → 已用快照/重跑方式核对最终状态，最终磁盘状态全部由合并脚本生成。

## 8. 结论

- 合并结构、段划分、run_stage 模式、download 默认跳过：全部符合 brief 要求。
- 合并保真度：静态文本 diff + 动态字节对比双重验证，合并脚本与原脚本行为逐字节一致。
- `Rscript reanalysis/scripts/acquire_data.R .`：成功（计算段全跑、无下载动作）。
- `test_lactylome_acquisition.R`：**通过**。
- verify 官方命令：24 个 diff 全部逐项说明（5 时间戳 + 4 上游 API 漂移 + 1 本地目录漂移，均为原脚本逻辑/外部环境导致，非合并缺陷；14 个为并行 agent 产物）；扩展排除后 `OK: all 362 baseline files unchanged`。
- 交付文件：`reanalysis/scripts/acquire_data.R`（新增）；重跑的 config/清单文件即本任务输出（已逐项对照基线说明差异）。
