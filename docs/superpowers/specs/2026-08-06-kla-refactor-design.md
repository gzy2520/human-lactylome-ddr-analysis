# Kla 项目重构设计 v2（全 R 合并重构）

日期：2026-08-06（v1 获批同日修订）
状态：修订版，待用户审阅

## 目标

在不改变任何分析功能的前提下，重构 kla 项目：

- **目录结构完全不动**（`data/`、`reanalysis/`、`archive/`、`previous_umap/` 原样）。
- **显著减少代码文件数量**：scripts/ 从 ~30 个文件合并为 **6 个 R 文件 + 1 个 lib/**。
- **尽可能全部用 R 完成分析**：所有 Python 分析/工具与 Node 工作簿工具迁移到 R；
  迁移后 scripts/ 内不再有 .py/.mjs/.js 文件。
- 不改变输出文件路径、列名、数值与分析规则；不删除数据（被替换文件移 `archive/`）。

## 一、目标文件结构（scripts/ 内，位置不动）

```
reanalysis/scripts/
├── lib/                          # 新增共享模块（R）
│   ├── accession_utils.R         # base_accession, is_uniprot, split_accessions,
│   │                             #   split_protein_identifiers, accession_feature,
│   │                             #   match_target_accession
│   ├── io_utils.R                # read_delimited, safe_numeric, relative_path,
│   │                             #   valid_maxquant_rows, write_csv_std, save_figure
│   └── extractors.R              # extract_maxquant_sites/proteins,
│                                 #   extract_spectronaut_*, extract_pd_*,
│                                 #   extract_huvec_xml, 及 R 版 Kla/GO/Venn 提取器
├── acquire_data.R                # 数据获取与清单（合并下列 R 脚本）
├── analyze_ddr.R                 # Kla/GO/Venn + 参照 DDR 占比 + 目标蛋白审计
├── analyze_regulators.R          # 调控因子 3 图（Kla 热图 / 全蛋白热图 / 景观图）
├── build_workbooks.R             # 交付工作簿/文档（原 .mjs/.js + metadata）
├── build_final_manifest.R        # sha256 清单（原 build_final_manifest.py）
└── __pycache__/、reanalysis/ 空目录 → 清理
```

## 二、合并映射

### 1. `acquire_data.R`（原 10 个 R 脚本 → 1 个）

`build_lactylome_reference_pairing.R`（1239 行）、`build_human_lactylome_inventory.R`（496）、
`build_lactylome_acquisition_manifests.R`（218）、`build_ensembl_uniprot_mapping.R`（161）、
`download_healthy_tissue_references.R`（157）、`download_lactylome_pair_files.R`（175）、
`extract_lactylome_pair_archives.R`（70）、`probe_lactylome_pair_files.R`（104）、
`register_additional_lactylome_pair_files.R`（73）、`summarize_acquired_lactylome_data.R`（214）、
`build_healthy_special_reference_manifest.R`（66）。
用参数/分段函数组织（`--stage download|pair|manifest` 或按函数调用）。

### 2. `analyze_ddr.R`（Python 主流程 + R DDR 分析 → 1 个）

- `run_pipeline.py`（718）＋ `extractors.py`（922）＋ `common.py`（235）→ R 重写
  （Kla 位点/蛋白提取、GO 附注、Venn 区域、肿瘤特异表、目标蛋白追踪、回归表）。
- `analyze_expanded_ddr_fraction_by_accession.R`（1474）→ 并入（去重复工具函数）。
- `analyze_reference_proteome_ddr.py`（728）→ R 重写并入。
- `audit_target_sources.py`（186）→ R 重写并入。

### 3. `analyze_regulators.R`（3 个 R 热图脚本 → 1 个）

`analyze_kla_regulator_intensity.R`（1197）、
`analyze_kla_regulator_whole_proteome_intensity.R`（1130）、
`plot_kla_regulator_landscape.R`（1019）→ 合并，共享 lib，保留各段独立函数
（`--analysis heatmap|whole_proteome|landscape` 或全跑）。

### 4. `build_workbooks.R`（Node 工具 + metadata → 1 个）

`build_cell_type_statistics_workbook.mjs`、`build_hippocampus_review_md.mjs`、
`build_reference_ddr_comparison_workbook.mjs`、`build_reference_proteome_selection_workbook.mjs`、
`build_venn_combined_workbook.mjs`、`create_bilingual_figure_legends_docx.js`、
`build_project_metadata.py` → 合并重写。
xlsx 用 writexl（需安装）；docx 用 officer 或等效 R 方案（需安装）。

### 5. `build_final_manifest.R`（原 .py → R）

sha256 用 `openssl::sha256` 或 `digest::digest(algo = "sha256")`；
输出路径、表头、行数与原版一致。

## 三、迁移关键点

| 原实现 | R 方案 | 注意 |
|---|---|---|
| matplotlib Venn（run_pipeline.py） | ggVennDiagram（已装） | 区域数字/标签/输出路径一致；渲染风格允许不同（软件差异） |
| pandas 大表（158MB 全蛋白 long） | data.table::fread/fwrite | 内存与速度可控 |
| PEAKS scan 归一化/marker-156 解析 | R 逐函数翻译 | 由 test_pipeline.R 逐点验证 |
| sha256 清单 | openssl::sha256 | 表头/路径格式不变 |
| exceljs 工作簿 | writexl（需安装） | 内容一致，格式从简可接受 |
| docx 双语图注 | officer（需安装） | 内容一致 |

需要安装的 R 包：`writexl`、`officer`（如 docx 方案需要）。

## 四、测试迁移

- `tests/test_pipeline.py`（14KB）→ `tests/test_pipeline.R`（纯 R 断言风格，
  不用 testthat；与现有 5 个 R 测试风格一致），覆盖原 14 项测试点。
- 现有 5 个 R 测试适配新脚本入口（调用方式与输出路径不变）。

## 五、功能不变验证

1. **R 测试全过**（重写后的 test_pipeline.R + 适配的 5 个 R 测试）。
2. **输出回归**：重构前后重跑主流程，全部输出文件 sha256 diff；
   统计表/中间表/清单必须字节一致；图因渲染库不同允许视觉差异但数字与标签一致。
3. **manifest 重建**：`Rscript reanalysis/scripts/build_final_manifest.R`，
   清单与磁盘一致。
4. **统计数字复核**：3112 Kla 蛋白、275 DDR、37 样本组等不变。

## 六、边界（不做）

- 不移动任何目录/文件位置；输出文件路径与列名不变。
- 不改 PXD 提取规则、阈值、GO 注释处理、ID 映射逻辑。
- 不删除任何数据；被合并/替代的 .py/.mjs/.js/.R 原文件移 `archive/`。
- 每 PXD 显式分支保留，不做数据驱动重构。
- `data/` 原始文件与 `previous_umap/` 不动。
