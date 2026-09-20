# Kla 31组蛋白组与 bulk RNA-seq 工作区

当前分支：`research/kla31-rnaseq-20260912`。蛋白组材料 31 行（包含 ESCC），RNA 对应 28 份不重复参考矩阵、1,898 个样本标识；脚本与输出分开保存。

**当前 RNA 发布入口：[config/rna_current_release.csv](config/rna_current_release.csv)**。不要从旧的 31/28 混版目录任意选择文件。

```sh
Rscript --vanilla workflow/preflight_kla31.R
```

- [当前工作交接](AI_HANDOFF_CURRENT.md)
- [2026-09-19 修复与验收](audit/20260919_rna_repair/REPORT.md)
- [修复前审计及原始证据](audit/20260919_project_review/REPORT.md)
- [当前 RNA 方法](docs/RNA_DATA_PROCESSING_AND_ANALYSIS_WORKFLOW.md)
- [RNA 简报](docs/RNA_WORKFLOW_BRIEF.md)
- [来源选择台账](audit/20260916_full_31_rna_status/rna_31_group_status.csv)

本轮按用户要求保留显著性分析方法，修复源表达换算、运行安全、图形和版本一致性。验收通过表示数据链和交付完整性通过，不表示已证明来源批次、供体独立性或所有统计假设。

原始 RNA 下载及 Stage 2 对象在隔离服务器目录，主要表达矩阵/元数据和当前交付同步本地。当前发布 SHA-256 清单记录本地文件；大原始数据不进入 Git。

蛋白组旧冻结包 `corrected_final_result_20260905_sample_only/` 保留，`final_result` 仅为该包快捷入口，不代表当前 RNA 版本。旧说明和方法保存在 `docs/history/`。旧 `workflow/run_pipeline.R` 为历史流程，不能代替当前 RNA 重建入口。
