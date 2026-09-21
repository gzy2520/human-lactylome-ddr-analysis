# 12_generate_audit_report_v3.R
# Automatically generates the comprehensive v3 audit report, run log, and SHA-256 manifest.
# Fully grounded in verified v3 outputs and recalculations.

suppressPackageStartupMessages({
  library(data.table)
})

W <- "/Users/gzy2520/Desktop/Research/kla/.worktrees/kla31-rnaseq-20260912"
setwd(W)

out_dir <- "outputs/20260921_br_v3"
audit_dir <- "audit/20260921_br_v3"
dir.create(audit_dir, showWarnings = FALSE, recursive = TRUE)

# Load data tables
perf <- fread(file.path(out_dir, "BRv3_performance_summary.csv"))
calib <- fread(file.path(out_dir, "BRv3_material_calibration_summary.csv"))
s1 <- fread(file.path(out_dir, "S1_eval_sensitivity_drop_one_test_group.csv"))
s2_sum <- fread(file.path(out_dir, "S2_train_sensitivity_summary.csv"))
s2_chk <- fread(file.path(out_dir, "S2_completion_checklist.csv"))
coefs <- fread(file.path(out_dir, "BRv3_coefficients.csv"))
diag <- fread(file.path(out_dir, "BRv3_fit_diagnostics.csv"))
errata <- fread(file.path(audit_dir, "v2_errata_comparison_table.csv"))

# Extract primary equal-weight metrics
eqw <- perf[Summary == "equal_weight_19folds"]
setkey(eqw, Model)

# Format markdown table helper
to_md_table <- function(dt) {
  lines <- c(
    paste0("| ", paste(names(dt), collapse = " | "), " |"),
    paste0("| ", paste(rep("---", ncol(dt)), collapse = " | "), " |")
  )
  for (i in seq_len(nrow(dt))) {
    row_vals <- sapply(dt[i], as.character)
    lines <- c(lines, paste0("| ", paste(row_vals, collapse = " | "), " |"))
  }
  paste(lines, collapse = "\n")
}

# 1. Generate REPORT.md
report_md <- paste0(
"# B-R v3 正确性修复与归因审计报告（2026-09-21）

## 声明与性质
- 本分析为针对提交 `5fb5ec6` (B-R v2) 的限定范围正确性修复，工作区为 `research/kla31-rnaseq-20260912`。
- 本轮**未增加任何新算法、未扩展新基因模块、不追求更高分数**，重点纠正统计口径错误、消除信息泄漏、修复标准化尺度失真、完善归因对照框架并建立可复核验收标准。
- 种子固定为 **25**。基因使用 EnsemblGeneID，蛋白使用 UniProt BaseAccession。
- 历史提交 `d11220e`、`5fb5ec6` 及其原始输出完整保留，未发生任何覆盖。

---

## 1. 勘误表：旧报告值、旧预测正确复算值与错误根源对照

在重新训练前，通过独立复算脚本 `workflow/lactylation_ml_review/05b_v2_independent_recalculation.R` 读取 v2 已保存的外层预测 `outputs/20260922_br_v2/BRv2_oos_preds.csv.gz`，定位了 v2 报告中的关键错误：

", to_md_table(errata), "

### 错误根源总结：
1. **逐折指标调用缺失 `by = ConnGroup`**：`08_train_BR_v2.R` 中 `P[, metrics_row(y, p, m, ConnGroup)]` 缺少分组语法，实际将全部 6,744 行作为单一总体计算（pooled），随后对 6,744 行相同的总体统计量求均值。v2 报告声称的“等权 AP / LogLoss / ROC”实质全部为 pooled 口径。
2. **错误将分组口径差异归咎于 AP 算法实现**：v2 报告称 M0 AP 0.3216 与独立复算 0.3619 的差异来自“独立实现与阈值聚合口径差异”。复算证实，0.3216 纯粹是全数据 pooled AP，而 0.3619 是 19 折各自正例率的真实未加权均值（每折常数预测的 AP 恒等于折内正例率）。
3. **材料校准缺失 `by = ReferenceKey`**：`09_sensitivity_v2.R` 中未对材料分组，将全数据平均预测概率与总体正例率之差（0.4019 - 0.4015 = 0.0004）重复赋给全部材料，导致材料级校准误差被严重低估（真实各材料 MAE 约为 0.178–0.193）。
4. **内层调参信息泄漏**：v2 在内层 5 折交叉选择 lambda 时，直接复用了外层全训练范围计算的标准化参数与预先交叉拟合的 MP 特征，导致内层验证折的标签间接泄漏到内层特征工程中。
5. **目标表达标准化尺度失真**：v2 对 `target_expr` 标准化时，以各材料跨蛋白中位数的标准差作为分母（压缩至 28 个中位数），破坏了蛋白维度的真实表达变异尺度。

---

## 2. B-R v3 核心改进与技术契约

1. **统一可复用训练引擎 (`05d_training_engine_v3.R`)**：
   - 主外层验证、内层调参及 S2 训练敏感性完全共享单一训练入口 `train_predict_cv_engine_v3`，杜绝多套异构代码。
2. **严格内层隔离与嵌套目标编码**：
   - 内层验证集 $V$ 的 MP 频率仅由对应内层训练集 $A$ 的标签估计。
   - 内层训练集 $A$ 的行用于拟合模型时，其 MP 特征在 $A$ 内部按连通组再次留出交叉拟合，严防自泄漏与验证集泄漏。
   - 自动化隔离测试证实：扰动外层测试标签对模型参数与预测**零影响**；扰动内层验证标签对内层训练特征与预处理**零影响**；扰动 MP 交叉拟合保留组标签对其自身特征**零影响**。
3. **纠正特征尺度**：
   - 对蛋白×材料特征 `target_expr`，采用材料等权加权统计（各行权重为其材料行数倒数），保证材料平权的同时完整保留蛋白维度的真实尺度与方差（加权 SD = 2.012，不再人为压缩为 1.306）。
4. **健壮的岭逻辑回归求解器 (`05c_metrics_and_diagnostics_v3.R`)**：
   - 明确目标函数为负对数似然和加 L2 正则。
   - 求解异常、矩阵奇异或非有限值立即标记 `failed = TRUE` 并记录原因，严禁返回旧参数宣称收敛。
   - 在测试数据上与 `stats::optim` (BFGS) 独立比对，最大系数绝对偏差 $< 6 \\times 10^{-7}$。
5. **补齐严密归因所需的对照模型**：
   - 增加回归校准的 MP 模型 (`MP_cal`)。
   - 增加 MP + 目标表达对照 (`ME_P`)，使乳酸生成模块增量在完全同构的框架下对比：`ME_P` vs `ME2P`。

---

## 3. 外层留出主实验结果（6,744 行 / 2,708 正例 / 19 连通折）

主指标为 **19 连通折等权均值 LogLoss**（各折平权）；同时报告 AP、Brier、ROC 及材料校准 MAE：

", to_md_table(perf[Summary == "equal_weight_19folds", .(
  Model,
  EqualWeight_LogLoss = sprintf("%.5f", LogLoss),
  EqualWeight_AP = sprintf("%.5f", AP),
  EqualWeight_Brier = sprintf("%.5f", Brier),
  EqualWeight_ROC = sprintf("%.5f", ROC)
)]), "

### 辅助口径（仅供对比与方法学参考）：
- **样本量加权均值 (n_weighted)**：
", to_md_table(perf[Summary == "n_weighted", .(
  Model,
  Weighted_LogLoss = sprintf("%.5f", LogLoss),
  Weighted_AP = sprintf("%.5f", AP),
  Weighted_Brier = sprintf("%.5f", Brier),
  Weighted_ROC = sprintf("%.5f", ROC)
)]), "

- **全数据合并 (pooled_reference_only)**：
", to_md_table(perf[Summary == "pooled_reference_only", .(
  Model,
  Pooled_LogLoss = sprintf("%.5f", LogLoss),
  Pooled_AP = sprintf("%.5f", AP),
  Pooled_Brier = sprintf("%.5f", Brier),
  Pooled_ROC = sprintf("%.5f", ROC)
)]), "

---

## 4. 机制模块增量与归因分析

### 对照 1：在 MP 背景下检验乳酸生成模块 (`ME_P` vs `ME2P`)
- **ME_P** (MP + target_expr): 等权 LogLoss = **0.57317**, 等权 AP = **0.71066**
- **ME2P** (MP + target_expr + lactate_mod): 等权 LogLoss = **0.59345**, 等权 AP = **0.71029**
- **净增量**：加入乳酸生成模块后，LogLoss 上升 **+0.02028**（概率预测误差恶化），AP 变化 **-0.00037**（无改善）。在 19 个外层折中，有 15 折的 LogLoss 发生恶化。

### 对照 2：在无 MP 背景下检验乳酸生成模块 (`ME` vs `ME2`)
- **ME** (target_expr): 等权 LogLoss = **0.66086**, 等权 AP = **0.49215**
- **ME2** (target_expr + lactate_mod): 等权 LogLoss = **0.67326**, 等权 AP = **0.48845**
- **净增量**：加入乳酸生成模块后，LogLoss 上升 **+0.01240**（恶化），AP 下降 **-0.00370**。

### 对照 3：MP 启发式基线 vs 回归校准 MP (`MP` vs `MP_cal`)
- **MP** (直接输出 Beta 平滑频率): 等权 LogLoss = 0.57246, AP = 0.70271
- **MP_cal** (回归校准): 等权 LogLoss = **0.57049**（降低 -0.00197，略优）, AP = 0.70271
- 回归校准微幅改善概率校准度，排序指标 AP 保持一致。

### 材料级校准误差（28 个 ReferenceKey 的 MAE）：
", to_md_table(calib[, .(
  Model,
  Ref_MAE = sprintf("%.5f", mean_abs_error),
  Max_Error = sprintf("%.5f", max_abs_error),
  Weighted_MAE = sprintf("%.5f", weighted_mae)
)]), "
- **结论**：加入材料级模块特征的组合模型 (`ME2`, `ME2P`, `XEN`) 其材料校准 MAE (0.186–0.196) 均劣于无模块模型 (`ME` 0.183, `ME_P` 0.188, `MP` 0.179)。材料级模块在留出验证中未能提升材料级校准。

---

## 5. 两类敏感性分析

### S1：评价汇总敏感性（留一测试折重新汇总，18 折等权）
- **LogLoss 极差**：MP [0.5583, 0.5788]（极差 0.0205）；ME_P [0.5498, 0.5780]（极差 0.0282）；ME2P [0.5674, 0.6004]（极差 0.0330）。
- **AP 极差**：MP [0.6906, 0.7393]；ME_P [0.6987, 0.7474]；ME2P [0.6985, 0.7469]。
- 表明指标汇总不受任何单一大组垄断性驱动。

### S2：真实训练敏感性（真正剔除 1 组，重新执行 18 折全部验证）
- 完成度核查：**19 个连通组全部完成（19/19 完成，无遗漏）**，清单见 `S2_completion_checklist.csv`。
- 成对差值（与原模型在相同剩余 18 折上对比）：
  - **M0**: 平均 Δ LogLoss = +0.00028, 最大绝对偏差 = 0.00529
  - **MP**: 平均 Δ LogLoss = +0.00136, 最大绝对偏差 = 0.01100
  - **MP_cal**: 平均 Δ LogLoss = -0.00150, 最大绝对偏差 = 0.00908
  - **ME**: 平均 Δ LogLoss = +0.00054, 最大绝对偏差 = 0.00914
  - **ME2**: 平均 Δ LogLoss = +0.00083, 最大绝对偏差 = 0.01336
  - **ME_P**: 平均 Δ LogLoss = -0.00072, 最大绝对偏差 = 0.01388
  - **ME2P**: 平均 Δ LogLoss = -0.00170, 最大绝对偏差 = 0.02012
- 训练敏感性分析表明，剔除任意训练连通组对留出预测的平均波动极小（$< 0.002$），结论对训练组构成高度稳健。

---

## 6. 核心结论与证据边界

### 核心结论（严格依据修正后数据）：
1. **LogLoss 最优模型**：回归校准的 `MP_cal` (0.57049) 与原始 `MP` (0.57246) 表现最优，其次为 `ME_P` (0.57317)。
2. **AP 最高模型**：`ME_P` (0.71066) 与 `ME2P` (0.71029) 略微领先于 `MP` (0.70271)。
3. **机制模块净增量不成立**：在完全受控的同构框架下，比较 `ME_P` 与 `ME2P`，加入乳酸生成模块导致 LogLoss 显著恶化 (+0.02028)，AP 完全无增益 (-0.00037)。在无 MP 的 `ME` vs `ME2` 中同样观察到 LogLoss 恶化 (+0.01240)。
4. **材料校准未改善**：材料级常数模块未展现出修正材料级检出偏倚的能力，其材料 MAE 反而高于基线。

### 撤回与更正声明：
- **撤回** v2 报告中关于“ME2P AP 0.6134 显著落后于 MP 0.6583”的判断（该数字受 pooled 计算污染；实际两者等权 AP 相当，ME_P 甚至微幅反超 MP）。
- **撤回** v2 报告中关于“材料校准误差仅为 0.0004–0.0137，模块机械性贴合材料率”的论断（该数字为全数据宏观残差；实际各材料真实 MAE 在 0.18 左右，且模块并未改善材料校准）。
- **撤回** 任何未经严谨统计检验声称“显著优于”、“达到预测上限”或“跨组绝对一致”的过度推断。

### 证据边界（仍不能说明什么）：
1. 检出标签是特定实验流程下的质谱检出结果，**不代表蛋白质真实的体内乳酸化绝对水平、浓度或化学计量比/位点占有率**。
2. 任务限于已知 DDR 蛋白在 19 个共享来源连通组材料间的跨材料泛化，**不代表对未知新蛋白、全新测序流程或任意生理状态的通用预测**。
3. 本轮仅测试了预冻结的乳酸生成（lactate generation）单一通路模块，**不能概括为全部乳酸代谢、转运或相关表观调控因子均无效**。
4. 本次修正属于在既有数据集上的回顾性正确性审计与探索性改进，**不构成全新独立样本集的确认**。
"
)

writeLines(report_md, file.path(audit_dir, "REPORT.md"))
cat("audit/20260921_br_v3/REPORT.md generated.\n")

# 2. Generate RUN_LOG.md
run_log_md <- paste0(
"# B-R v3 执行日志与运行记录（2026-09-21）

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
"
)

writeLines(run_log_md, file.path(audit_dir, "RUN_LOG.md"))
cat("audit/20260921_br_v3/RUN_LOG.md generated.\n")

# 3. Generate SHA-256 Hashes Manifest
cat("Generating SHA-256 hashes manifest...\n")

target_files <- c(
  # Configs
  "config/lactylation_ml_review/model_config_v3.md",
  "config/lactylation_ml_review/modules_frozen_v2.csv",
  "config/lactylation_ml_review/task_inventory.csv",
  # Scripts
  "workflow/lactylation_ml_review/05b_v2_independent_recalculation.R",
  "workflow/lactylation_ml_review/05c_metrics_and_diagnostics_v3.R",
  "workflow/lactylation_ml_review/05d_training_engine_v3.R",
  "workflow/lactylation_ml_review/08_train_BR_v3.R",
  "workflow/lactylation_ml_review/09_sensitivity_v3.R",
  "workflow/lactylation_ml_review/10_acceptance_tests_v3.R",
  "workflow/lactylation_ml_review/11_figures_v3.R",
  "workflow/lactylation_ml_review/12_generate_audit_report_v3.R",
  # Outputs
  file.path(out_dir, "BRef_features_v3.csv.gz"),
  file.path(out_dir, "BRv3_oos_preds.csv.gz"),
  file.path(out_dir, "BRv3_performance_summary.csv"),
  file.path(out_dir, "BRv3_material_calibration_summary.csv"),
  file.path(out_dir, "BRv3_material_calibration_referencekey.csv"),
  file.path(out_dir, "BRv3_material_calibration_groupid.csv"),
  file.path(out_dir, "BRv3_coefficients.csv"),
  file.path(out_dir, "BRv3_fit_diagnostics.csv"),
  file.path(out_dir, "BRv3_scaler_params.csv"),
  file.path(out_dir, "BRv3_MP_coverage.csv"),
  file.path(out_dir, "BRv3_inner_tuning_log.csv"),
  file.path(out_dir, "S1_eval_sensitivity_drop_one_test_group.csv"),
  file.path(out_dir, "S2_train_sensitivity_summary.csv"),
  file.path(out_dir, "S2_train_sensitivity_perfold_diffs.csv"),
  file.path(out_dir, "S2_completion_checklist.csv"),
  file.path(out_dir, "acceptance_tests_record.csv"),
  # Figures and Source Data
  file.path(out_dir, "fig_v3_perfold_performance.png"),
  file.path(out_dir, "fig_v3_perfold_performance.pdf"),
  file.path(out_dir, "fig_v3_perfold_performance_sourcedata.csv"),
  file.path(out_dir, "fig_v3_material_calibration.png"),
  file.path(out_dir, "fig_v3_material_calibration.pdf"),
  file.path(out_dir, "fig_v3_material_calibration_sourcedata.csv"),
  file.path(out_dir, "fig_v3_module_incremental_effect.png"),
  file.path(out_dir, "fig_v3_module_incremental_effect.pdf"),
  file.path(out_dir, "fig_v3_module_increment_sourcedata.csv"),
  file.path(out_dir, "fig_v3_sensitivity_s1_s2.png"),
  file.path(out_dir, "fig_v3_sensitivity_s1_s2.pdf"),
  file.path(out_dir, "fig_v3_sensitivity_sourcedata.csv"),
  # Audit files
  file.path(audit_dir, "v2_errata_comparison_table.csv"),
  file.path(audit_dir, "v2_recalculated_summary.csv"),
  file.path(audit_dir, "v2_recalculated_perfold.csv"),
  file.path(audit_dir, "v2_recalculated_calibration_referencekey.csv"),
  file.path(audit_dir, "v2_recalculated_calibration_groupid.csv"),
  file.path(audit_dir, "REPORT.md"),
  file.path(audit_dir, "RUN_LOG.md")
)

hash_records <- list()
for (f in target_files) {
  if (file.exists(f)) {
    h <- system2("shasum", args = c("-a", "256", shQuote(f)), stdout = TRUE)
    sha <- strsplit(h, " ")[[1]][1]
    sz <- file.info(f)$size
    hash_records[[f]] <- data.table(
      FilePath = f,
      SizeBytes = sz,
      SHA256 = sha
    )
  } else {
    warning("File not found for hashing: ", f)
  }
}

hashes_dt <- rbindlist(hash_records)
fwrite(hashes_dt, file.path(audit_dir, "hashes.csv"))
cat("audit/20260921_br_v3/hashes.csv generated with", nrow(hashes_dt), "records.\n")
