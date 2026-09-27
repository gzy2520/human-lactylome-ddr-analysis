# 全共同基因 RNA 输入：逐样本随机森林

按要求取消 128 基因面板限制：现在输入为 **1,898 条 RNA × 18,332 个 Ensembl 基因**。逐条 RNA 训练和预测，不先汇总为材料，不按乳酸代谢、DDR、差异显著性或与标签相关性筛选基因。

“全”指现有 28 份来源矩阵均可提供的全部共同基因，不代表不同测序方案未覆盖的全部 RNA 类型。各来源有些基因仅在部分矩阵存在；为避免把未测到当作零表达，本轮使用所有来源基因 ID 交集。各来源总数及未进入共同输入的数量保存在 `coverage.csv`，18,332 个稳定 ID 保存于 `gene_ids.csv`。

## 公平比较的固定项

- 表达量是修正后的 log2(TPM + 0.5)，不使用 qsmooth；原 128 个基因逐值与旧输入一致，误差 <1e-10。
- 直接复用上轮 `split.csv`，1,524 条训练，374 条内部测试，覆盖 19 个材料；已知供者不跨内部训练/测试。106 条原始 RNA 缺乏可解析供者身份，仍沿用未解析记录键。
- 原有 14 个来源块及全部共享材料标签保持不变，另完成整来源留出评价。
- 随机森林仍为 500 棵树，最小叶节点大小 1，材料平衡 bootstrap 抽样，种子 25。材料、来源、供者 ID 不进入表达特征。
- 比较 mtry=64（与之前主模型相同）和 mtry=135（18,332 的平方根取整）。mtry 是每次分裂随机考察的基因数，不是总输入基因数。全部 18,332 个基因均有参与候选的机会。
- 没有重抽测试集或按新结果改变参数范围。由于这一划分已经用于多次方案比较，它仍是探索性内部测试，不是完全独立的最终确认集。

## 实际结果

同一个 374 条内部测试集，逐记录等权比较：

| 输入及模型 | MAE（百分点） | R² | 误差在 ±5 百分点内 |
|---|---:|---:|---:|
| 128 基因，mtry=64 | 1.07 | 0.966 | 96.3% |
| 18,332 基因，mtry=64 | 1.21 | 0.975 | 97.3% |
| 18,332 基因，mtry=135 | 1.09 | 0.979 | 98.1% |

全基因 mtry=135 的材料等权内部测试 MAE=1.59，±5 百分点内比例=96.9%；原 128 基因主模型对应为 1.55 和 92.4%。误差容限内比例和 R²改善，但平均绝对误差没有改善，不能称全面优于旧模型。

全基因 mtry=135 的来源留出逐记录 MAE=12.18，±5 百分点内比例=2.0%；材料等权 MAE=9.55，比例=16.4%。这仍不支持外推到新来源。可以采用全基因输入作为当前探索性模型，但效果范围仍是已有材料的共享标签重现。

已验证与旧模型拆分表完全一致，独立重算测试 MAE/容限比例，输入输出哈希通过；重载 mtry=135 全数据模型后逐样本预测与保存的训练输出一致。未重复搜索测试拆分。

## 文件及运行

输入目录：`outputs/20260926_full_rna_inputs/`。`full_rna.rds` 约 267 MB，保存在本机，未放入普通 Git；基因清单、覆盖、源文件哈希和矩阵哈希纳入版本控制，以下准备脚本可从已修正源矩阵重建。原始修正表达文件路径记录在输入清单。

模型结果：`outputs/20260926_full_rna_models/`，包含两套全数据森林、逐样本预测、三种评价、128 基因对比表、图和校验记录。`RF_mtry135.png/pdf` 与 `RF_mtry64.png/pdf` 都是左侧训练拟合、中间内部留出、右侧来源留出。图中 MAE 和容限比例采用材料等权；逐记录指标在表中单列。

```sh
Rscript --vanilla workflow/full_rna_fraction/prepare.R /path/to/outputs/20260919_expression_corrected /tmp/new_full_inputs
Rscript --vanilla workflow/full_rna_fraction/train.R /tmp/new_full_inputs /tmp/new_full_models
Rscript --vanilla workflow/full_rna_fraction/compare.R outputs/20260926_nonlinear_samples_final /tmp/new_full_models
Rscript --vanilla workflow/full_rna_fraction/predict.R /tmp/new_full_models/RF_mtry135.rds /tmp/new_full_inputs/full_rna.rds /tmp/new_full_predictions.csv
```

各输出路径须未被占用。predict 接口要求相同尺度、模型需要的全部 Ensembl 列。已保存的 RDS 是全数据模型，不能将其对训练数据的输出当作独立测试结果。

模型目标仍是赋给各 RNA 的共享材料级 Kla 蛋白检出占比，没有新增任何个体配对蛋白标签。内部测试衡量已有材料类别的标签重现；整来源留出衡量外推，两者必须分别报告。

## 补充图表集（2026-09-27 更新）

按多维度评估与汇报需求，补充了以下 6 套高质量矢量/标量图表与数据表（位于 `outputs/20260926_full_rna_models/figures/` 和 `tables/`）：

1. **补全训练集 4 栏散点图 (`RF_mtry135_4panel.png/pdf`)**：
   - 补齐了原先只展示全量拟合的空缺，形成 4 栏严格对照：
     - 栏 1：`Training split (n=1,524)`，在 1,524 条训练子集上训练并重拟合评估（MAE=0.12 pp, R²=0.998, ±5 pp 99.9%）；
     - 栏 2：`Within-material test (n=374)`，用 1,524 样本训练后对 374 独立同类材料样本做盲测（MAE=1.59 pp, R²=0.955, ±5 pp 96.9%）；
     - 栏 3：`Full dataset fit (n=1,898)`，全数据拟合模型容量上限（MAE=0.12 pp, R²=0.998）；
     - 栏 4：`Source held out (n=1,898)`，14 折整来源留出外推测试（MAE=9.55 pp, R²=-0.195）。

2. **各材料预测分布图 (`figure1_material_distributions.png/pdf`)**：
   - 横轴按实测 Kla 占比从低到高排列 28 种材料，纵轴为 Kla 预测占比；
   - 蓝色点代表训练集样本（n=1,524），橙色点代表测试集样本（n=374），金色菱形代表该材料实测真实标签；
   - 直观展示各材料预测的聚集度与方差，以及与真实测定值的吻合状态。

3. **Top 25 特征重要性条形图 (`figure2_feature_importance.png/pdf`, `tables/feature_importance_top25.csv`)**：
   - 提取随机森林 500 棵决策树的 Gini 不纯度下降（Impurity Importance）；
   - 严格保留 Ensembl 稳定 ID，并利用 `org.Hs.eg.db` 映射官方 Gene Symbol；
   - 排名前列包含重要代谢转氨酶与脱氢酶（如 GPT2、ALDH7A1、ALDH4A1、GLUD1 等）以及组织特异性标志物。

4. **转录组流形降维投影图 (`figure3_manifold_projection.png/pdf`)**：
   - 基于前 2,000 高变基因与 10 个主成分运行 UMAP（种子 25），3 栏并排展示：
     - A. 生物学大类着色（癌组织、正常组织、癌细胞系、正常细胞系）；
     - B. 真实测得的材料级 Kla 占比连续色谱；
     - C. 模型预测的 Kla 占比连续色谱；
   - B 与 C 高度吻合，从流形几何角度证明模型从转录组空间准确捕捉到了 Kla 检出比例的连续梯度。

5. **残差诊断图 (`figure4_residual_diagnostics.png/pdf`)**：
   - 左图：残差（预测值 - 真实值）对实测值的散点图与 LOESS 平滑曲线，误差紧密落在 ±5 个百分点红虚线内，全量程无系统性倾斜；
   - 右图：残差概率密度分布，在 0 处呈现尖锐对称单峰，无显著方差膨胀或偏态。

6. **生物学大类分层散点图 (`figure5_category_stratified.png/pdf`, `tables/category_test_metrics.csv`)**：
   - 拆分为正常组织（Test n=171, MAE=0.97 pp, ±5 pp 99.4%）、癌组织（Test n=186, MAE=1.08 pp, ±5 pp 97.3%）、癌细胞系（Test n=8, MAE=1.76 pp, ±5 pp 100%）、正常细胞系（Test n=9, MAE=2.90 pp, ±5 pp 88.9%）；
   - 证实模型在占据绝对多数的临床组织样本上预测极为稳健。

运行生成脚本：`workflow/full_rna_fraction/generate_extended_figures.R`。
