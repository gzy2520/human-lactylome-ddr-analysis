# SDD ledger — plan: docs/superpowers/plans/2026-08-06-kla-full-r-refactor.md

适配说明：本项目非 git 仓库，无 commit/worktree。检查点机制 =
archive/refactor_baseline_2026-08-06/output_sha256_baseline.csv（Task 1 建立）
+ 每任务"任务前文件快照"（BASE）与 verify_outputs diff。
裁决记录：
- Task 13 归档被合并/迁移原文件后，同步从基线中剔除这些条目并记录（否则 Task 14
  verify 会把已归档文件报 MISSING）；剔除清单随归档清单走。
Task 2: complete (writexl/officer installed, requireNamespace OK — controller 直跑，因 Agent 服务暂不可用)
Task 1: complete (基线 362 文件建立；verify_outputs.R 用 digest::digest(file=) 修复 openssl raw 返回问题；6 测试文件全绿)
Task 3: complete (lib 4 文件创建；358/362 数据文件字节不变；controller 逐函数核对等价)
  - 裁决1: 4 个 cell_type_kla_vs_reference_ddr_fraction* 图 CHANGED = 渲染环境差异（PNG 重跑稳定、
    PDF 含时间戳非确定）；数据表不变→图数字绑定不变；后续 verify 排除该前缀图，Task 14 目检
  - 裁决2: extract_maxquant_proteins 依赖 map_ensembl_proteins（expanded 私有）→ Task 7 在
    analyze_ddr.R 内定义后再用
  - 裁决3: match_target_accession 闭包→显式 target_accessions 参数；relative_path 闭包→显式 root
    参数 → Task 4 合并热图脚本时改调用点传参（target_accessions <- unique(regulators$BaseAccession)
    保留为脚本级变量）
Task 4: in_progress
用户指示: 并行派发多个 subagent 加速。约束: 仅文件不重叠的任务并行。
Task 5/6/10/11 并行中 (独立新脚本, 独立输出文件, 均已核实无重叠)。
会话进程退出事件: 5 个并行代理(4/5/6/10/11)被停止且无落地产出(磁盘核查: 无目标新文件、
无任务报告)；lib/ 完好。决定: 重新派发 fresh agents 重做 Task 4/5/6/10/11(并行)。
Task 11: complete (build_final_manifest.R 与 Python 逐字节一致, 498 文件, MANIFEST_IDENTICAL)
Task 5: complete (acquire_data.R 3063 行合并 10 脚本, 原脚本未动, 测试过; 24 diff 全部可解释:
  5 Sys.Date 时间戳 + 4 上游 API 实时漂移 + 1 data/ 目录新增 + 14 并行代理产物)
  - 裁决: acquire_data 段内保留与 lib 语义不同的本地 split_accessions/base_accession
    (实测替换破坏基线计数), 属设计内"私有变体保留"边界; Task 14 verify 排除正则需涵盖
    时间戳/API 漂移类文件
Task 4/6/10: in_progress (运行中)
[2026-08-06] 用户暂停重构，转做新进展。3 个代理(Task 4/6/10)已停。
暂停时磁盘状态:
  - Task 4: analyze_regulators.R 存在(103KB, 12:58) 未验证/无报告
  - Task 6: analyze_ddr.R 存在(122KB, 13:16) 未验证/无报告; test_pipeline.R 未创建
  - Task 10: build_workbooks.R 存在(61KB, 13:18) 未验证/无报告
恢复方法: 对 3 个部分产出先跑语法检查+测试+verify 基线 diff；若验证不过，
          重新派发代理从零执行或修复；ledger 与 briefs 均在 .superpowers/sdd/kla-full-r-refactor/。
