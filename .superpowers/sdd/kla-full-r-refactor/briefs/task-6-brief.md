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

