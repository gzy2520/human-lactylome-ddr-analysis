# B-R v3 执行日志与运行记录（2026-09-21）

## 环境信息
- 机器：macOS (Darwin aarch64), Apple Silicon M-series
- 统计环境：R version 4.4.3 (2025-02-28), renv 隔离运行
- 工作路径：`/Users/gzy2520/Desktop/Research/kla/.worktrees/kla31-rnaseq-20260912`
- 分支：`research/kla31-rnaseq-20260912`，基于提交 `5fb5ec6`
- 随机种子：25

## 执行步骤与耗时
1. **独立复算旧预测 (`05b_v2_independent_recalculation.R`)**：
   - 耗时：~1.8 秒
   - 结果：完全复现提示词中给定的参考值，输出 `v2_errata_comparison_table.csv`。
2. **冻结配置与规范 (`config/lactylation_ml_review/model_config_v3.md`)**：
   - 完成 8 个模型比较定义、优化目标、严格内层隔离与加权标准化公式冻结。
3. **核心工具与训练引擎构建 (`05c_metrics_and_diagnostics_v3.R`, `05d_training_engine_v3.R`)**：
   - 编写阈值聚合 AP、梯形 PR-AUC、ROC、IRLS 岭回归及严格内层隔离逻辑。
   - 求解器与 `stats::optim` 比对测试通过（偏差 $< 6e-7$）。
   - 三重隔离单元测试通过（扰动外层标签、扰动内层验证标签、扰动 MP 交叉拟合保留组标签，差异均为 0）。
4. **主外层 19 折验证重跑 (`08_train_BR_v3.R`)**：
   - 耗时：25.26 秒
   - 结果：6,744 行外层预测全部无缺失生成，19 折全部数值收敛，输出主性能汇总表与材料校准表。
5. **两类敏感性分析 (`09_sensitivity_v3.R`)**：
   - S1 评价汇总敏感性完成（19 折逐折剔除后均值重算）。
   - S2 真实训练敏感性重跑完成：19 个连通组各剔除一次，每次完整重跑 18 折嵌套验证，耗时 521.73 秒。
   - 19/19 运行检查清单全部记录为 COMPLETED。
6. **自动化验收测试套件 (`10_acceptance_tests_v3.R`)**：
   - 17 项验收测试全部 PASS（0 失败）。
7. **图表与源数据自动生成 (`11_figures_v3.R`)**：
   - 生成 4 套核心图件（PNG 与 PDF，均通过渲染与文件体积校验），并输出对应源数据 CSV。
8. **审计报告与哈希清单 (`12_generate_audit_report_v3.R`)**：
   - 自动聚合生成最终 REPORT.md、RUN_LOG.md 与 hashes.csv。

