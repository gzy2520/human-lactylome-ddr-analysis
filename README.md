# Kla 31组蛋白组与 bulk RNA-seq 工作区

当前分支：`research/kla31-rnaseq-20260912`，基于远端main `bc02d1d`，包含sample-only修正与强度审计。31组包含ESCC；分类9/3/12/7。脚本与新输出分别保存。

先运行只读检查：

```sh
Rscript --vanilla workflow/preflight_kla31.R
```

- [当前工作约定](AI_HANDOFF_CURRENT.md)
- [RNA候选清单](docs/MATCHED_RNASEQ_REFERENCE_CATALOG_31.md)
- [2026-09-14组织匹配补充与普通全蛋白复审](outputs/20260914_reference_material_audit/REPORT.md)：新增RNA候选；MCF10A蛋白FDR问题、HK-2样本别名证据缺口单列。
- [此前add_data工作区审计](audit/20260912_workspace_rnaseq/REPORT.md)（历史快照，其分支状态不描述当前工作树）
- 31组冻结注册表：`audit/20260912_workspace_rnaseq/groups_31.csv`
- 已有sample-only基线：`corrected_final_result_20260905_sample_only/`；`final_result`仅为其本地快捷入口，不代表最终批准的新版本。
- 新脚本放`R/`、`workflow/`；新输出放`outputs/<日期_任务>/`，不写回已有结果目录。
- 原30组说明保存在`docs/history/`，旧`workflow/run_pipeline.R`须显式启用历史模式才能运行。

本地已独立复制31组扩展输入（约26 MB）和R库目录；未复制45 GB原始蛋白数据。输入哈希清单见`config/kla31_local_inputs_sha256.csv`。外部视觉参考仍为`/Users/gzy2520/Desktop/renew/kla`；当前未跟踪Word副本在被忽略的`local_reference/manuscript/`，不替换main内的文稿。

RNA仍为候选来源，尚未完成表达矩阵验收。2026-09-14已新增PC-3M身份匹配候选；NSC具体模型等剩余限制见补充报告。全量源重建所需原始数据仍在原工作区，不能仅凭本检查PASS声称完整重建已通过。
