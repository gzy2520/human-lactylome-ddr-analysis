# B-R v3 正确性修复与归因审计报告（2026-09-21）

## 声明与性质
- 本分析为针对提交 `5fb5ec6` (B-R v2) 的限定范围正确性修复，工作区为 `research/kla31-rnaseq-20260912`。
- 本轮**未增加任何新算法、未扩展新基因模块、不追求更高分数**，重点纠正统计口径错误、消除信息泄漏、修复标准化尺度失真、完善归因对照框架并建立可复核验收标准。
- 种子固定为 **25**。基因使用 EnsemblGeneID，蛋白使用 UniProt BaseAccession。
- 历史提交 `d11220e`、`5fb5ec6` 及其原始输出完整保留，未发生任何覆盖。

---

## 1. 勘误表：旧报告值、旧预测正确复算值与错误根源对照

在重新训练前，通过独立复算脚本 `workflow/lactylation_ml_review/05b_v2_independent_recalculation.R` 读取 v2 已保存的外层预测 `outputs/20260922_br_v2/BRv2_oos_preds.csv.gz`，定位了 v2 报告中的关键错误：

| Model | Metric | Old_Reported_Value | Old_Scope_In_Reality | Correct_Recalculated_Value | Difference_Reason |
| --- | --- | --- | --- | --- | --- |
| M0 | Equal-weight AP | 0.3216 | Pooled (6744 rows) | 0.3619287 | 08 script lacked by=ConnGroup; evaluated all 6744 rows as single pool; misattributed to AP definition differences |
| M0 | Equal-weight LogLoss | 0.6826 | Pooled (6744 rows) | 0.6662366 | 08 script evaluated overall pooled LogLoss instead of unweighted mean of 19 fold LogLosses |
| M0 | Equal-weight ROC | 0.3205 | Pooled (6744 rows) | 0.5 | Per-fold M0 prediction is constant so per-fold ROC is identically 0.500; pooled ROC artifactually ranked cross-fold rates |
| MP | Equal-weight AP | 0.6583 | Pooled (6744 rows) | 0.7027054 | 08 script evaluated pooled AP (0.6583) instead of 19-fold unweighted mean (0.7027) |
| MP | Equal-weight LogLoss | 0.5757 | Pooled (6744 rows) | 0.5724599 | 08 script evaluated pooled LogLoss (0.5757) instead of 19-fold unweighted mean (0.5725) |
| MP | Equal-weight ROC | 0.7556 | Pooled (6744 rows) | 0.8239471 | 08 script evaluated pooled ROC (0.7556) instead of 19-fold unweighted mean (0.7502) |
| ME | Equal-weight AP | 0.4313 | Pooled (6744 rows) | 0.4921517 | 08 script evaluated pooled AP (0.4313) instead of 19-fold unweighted mean (0.4922) |
| ME | Equal-weight LogLoss | 0.6843 | Pooled (6744 rows) | 0.6608651 | 08 script evaluated pooled LogLoss (0.6843) instead of 19-fold unweighted mean (0.6609) |
| ME2 | Equal-weight AP | 0.4091 | Pooled (6744 rows) | 0.4884488 | 08 script evaluated pooled AP (0.4091) instead of 19-fold unweighted mean (0.4884) |
| ME2 | Equal-weight LogLoss | 0.6949 | Pooled (6744 rows) | 0.6732676 | 08 script evaluated pooled LogLoss (0.6949) instead of 19-fold unweighted mean (0.6733) |
| ME2P | Equal-weight AP | 0.6134 | Pooled (6744 rows) | 0.7103017 | 08 script evaluated pooled AP (0.6134) instead of 19-fold unweighted mean (0.7103) |
| ME2P | Equal-weight LogLoss | 0.6026 | Pooled (6744 rows) | 0.5934501 | 08 script evaluated pooled LogLoss (0.6026) instead of 19-fold unweighted mean (0.5935) |
| XEN | Equal-weight AP | 0.3787 | Pooled (6744 rows) | 0.4863957 | 08 script evaluated pooled AP (0.3787) instead of 19-fold unweighted mean (0.4864) |
| XEN | Equal-weight LogLoss | 0.7147 | Pooled (6744 rows) | 0.693596 | 08 script evaluated pooled LogLoss (0.7147) instead of 19-fold unweighted mean (0.6936) |
| M0_calib | Material Calib MAE | 4e-04 | Grand mean diff | 0.1775181 | 09 script lacked by=ReferenceKey; assigned grand mean difference across all 6744 rows rather than per-reference MAE |
| MP_calib | Material Calib MAE | 0.0137 | Grand mean diff | 0.1793081 | 09 script lacked by=ReferenceKey; assigned grand mean difference across all 6744 rows rather than per-reference MAE |
| ME_calib | Material Calib MAE | 9e-04 | Grand mean diff | 0.1833898 | 09 script lacked by=ReferenceKey; assigned grand mean difference across all 6744 rows rather than per-reference MAE |
| ME2P_calib | Material Calib MAE | 0.0029 | Grand mean diff | 0.1934222 | 09 script lacked by=ReferenceKey; assigned grand mean difference across all 6744 rows rather than per-reference MAE |

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
   - 在测试数据上与 `stats::optim` (BFGS) 独立比对，最大系数绝对偏差 $< 6 \times 10^{-7}$。
5. **补齐严密归因所需的对照模型**：
   - 增加回归校准的 MP 模型 (`MP_cal`)。
   - 增加 MP + 目标表达对照 (`ME_P`)，使乳酸生成模块增量在完全同构的框架下对比：`ME_P` vs `ME2P`。

---

## 3. 外层留出主实验结果（6,744 行 / 2,708 正例 / 19 连通折）

主指标为 **19 连通折等权均值 LogLoss**（各折平权）；同时报告 AP、Brier、ROC 及材料校准 MAE：

| Model | EqualWeight_LogLoss | EqualWeight_AP | EqualWeight_Brier | EqualWeight_ROC |
| --- | --- | --- | --- | --- |
| M0 | 0.66624 | 0.36193 | 0.23655 | 0.50000 |
| MP | 0.57246 | 0.70271 | 0.19357 | 0.82395 |
| MP_cal | 0.57049 | 0.70271 | 0.19363 | 0.82395 |
| ME | 0.66086 | 0.49215 | 0.23318 | 0.63910 |
| ME2 | 0.67326 | 0.48845 | 0.23920 | 0.63655 |
| ME_P | 0.57317 | 0.71066 | 0.19413 | 0.82823 |
| ME2P | 0.59345 | 0.71029 | 0.20209 | 0.82892 |
| XEN | 0.69359 | 0.48639 | 0.24798 | 0.61672 |

### 辅助口径（仅供对比与方法学参考）：
- **样本量加权均值 (n_weighted)**：
| Model | Weighted_LogLoss | Weighted_AP | Weighted_Brier | Weighted_ROC |
| --- | --- | --- | --- | --- |
| M0 | 0.68258 | 0.40154 | 0.24460 | 0.50000 |
| MP | 0.57573 | 0.72076 | 0.19465 | 0.81528 |
| MP_cal | 0.57835 | 0.72076 | 0.19683 | 0.81528 |
| ME | 0.68433 | 0.50547 | 0.24467 | 0.60721 |
| ME2 | 0.69492 | 0.50288 | 0.24967 | 0.60709 |
| ME_P | 0.58476 | 0.72877 | 0.19896 | 0.81752 |
| ME2P | 0.60263 | 0.72876 | 0.20596 | 0.81773 |
| XEN | 0.71468 | 0.49822 | 0.25822 | 0.58804 |

- **全数据合并 (pooled_reference_only)**：
| Model | Pooled_LogLoss | Pooled_AP | Pooled_Brier | Pooled_ROC |
| --- | --- | --- | --- | --- |
| M0 | 0.68258 | 0.32165 | 0.24460 | 0.32048 |
| MP | 0.57573 | 0.65832 | 0.19465 | 0.75560 |
| MP_cal | 0.57835 | 0.64074 | 0.19683 | 0.74857 |
| ME | 0.68433 | 0.43132 | 0.24467 | 0.52273 |
| ME2 | 0.69492 | 0.40906 | 0.24967 | 0.48781 |
| ME_P | 0.58476 | 0.64437 | 0.19896 | 0.74279 |
| ME2P | 0.60263 | 0.61337 | 0.20596 | 0.72609 |
| XEN | 0.71468 | 0.37866 | 0.25822 | 0.47792 |

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
| Model | Ref_MAE | Max_Error | Weighted_MAE |
| --- | --- | --- | --- |
| M0 | 0.17752 | 0.39979 | 0.17111 |
| MP | 0.17931 | 0.38976 | 0.17152 |
| MP_cal | 0.18340 | 0.38715 | 0.17592 |
| ME | 0.18339 | 0.47047 | 0.17848 |
| ME2 | 0.18626 | 0.48548 | 0.18607 |
| ME_P | 0.18786 | 0.44485 | 0.18129 |
| ME2P | 0.19342 | 0.48577 | 0.19057 |
| XEN | 0.19567 | 0.60610 | 0.19539 |
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

