# 全 RNA 输入预测 DDR / DNA repair 的 Kla 检出比例

2026-09-27。旧总蛋白占比模型保留，本目录单独构建老师指定的两个新目标。统计与建模使用 R，稳定 Ensembl / UniProt ID，seed=25。

## 目标定义

DDR = GO:0006974 及其 is_a/part_of 后代；DNA repair = GO:0006281 及其后代。使用项目冻结的人类 UniProt GO 注释，排除 NOT；后代列表保存于本轮 QuickGO 快照。它不是一次重新下载全人类最新蛋白注释的更新。UniProt isoform 后缀折叠为 BaseAccession，蛋白组映射歧义按既有映射展开。GO 定义不以 Kla 阳性蛋白反向筛选。

每一蛋白材料组的目标为：

`100 × |普通蛋白组 ∩ Kla 蛋白组 ∩ GO 蛋白集合| / |普通蛋白组 ∩ GO 蛋白集合|`

计数的是不同蛋白 accession 的检出比例，不是位点比例、修饰强度或占有率。普通蛋白组未检出的 Kla 蛋白另列 KlaOutsideReference，不挤入分子。31 个蛋白材料组映射到 28 个 RNA ReferenceKey；同一个 ReferenceKey 对应多组时沿用旧框架取组比例算术平均，保留各组范围。全部 62 个目标/组组合分母非零，无需排除材料。

DDR 注释集合 2660 个 BaseAccession，每组分母 83–535；DNA repair 2135 个，每组分母 54–325。后者蛋白集合是前者的子集。集合大小不是表达矩阵的特征数量。

## 输入、切分及模型

使用冻结的 1898 × 18332 全 RNA 交集矩阵（log2 TPM 口径），没有先限于 DDR 或 DNA repair RNA 基因。输入来自材料匹配的 RNA 参考，标签来自蛋白组材料；并非 1898 名供者均实测了该占比。

沿用既有 1524 / 374 划分，已知供者不跨训练和留出；106 条未知供者记录沿用伪 ID，不能据此声称所有人的身份都被确认隔离。内部留出覆盖 19 种材料，其余小样本材料未进入此内部验证。另运行 14 个连接来源块留出；连接块用于合并相关来源/共享材料或供者，不等同于 14 个独立实验室。

两个目标各使用同一固定随机森林：ranger、500 棵树、mtry=135、min.node.size=1、seed=25。case.weights 使材料等权、材料内供者等权、供者内记录等权。本轮没有根据新目标内部测试结果挑选超参数，但划分和方法已有开发历史，所以不把内部留出称为全新独立盲测。GO 集合仅用于定义输出标签。

每个目标计算：1524 条训练内拟合、374 条内部留出、1898 条全量再拟合、14 折整来源留出。落盘 RF_mtry135.rds 是全量模型，不能用于复算内部留出性能；该性能来自只训练 1524 条的独立拟合。训练集 Top25 使用 1524 条模型的重要性；覆盖全部 18332 基因后排名，未以 Top25 重新训练，也不是 Gini 或因果贡献。

## 本次性能

逐 RNA 记录结果：

|目标|评价|MAE（百分点）|R²|±5 个百分点内|
|---|---|---:|---:|---:|
|DDR|内部留出 374|1.274|0.982|96.79%|
|DNA repair|内部留出 374|1.473|0.983|94.65%|
|DDR|整来源留出 1898|15.414|-0.508|1.32%|
|DNA repair|整来源留出 1898|17.807|-0.527|1.21%|

考虑不同材料 RNA 数量悬殊，主图采用材料等权、材料内供者等权评价：

|目标|评价|MAE（百分点）|R²|±5 个百分点内|
|---|---|---:|---:|---:|
|DDR|内部留出 374 / 19 材料|2.462|0.950|84.44%|
|DNA repair|内部留出 374 / 19 材料|2.794|0.950|74.42%|
|DDR|整来源留出 / 28 材料|13.617|-0.090|14.58%|
|DNA repair|整来源留出 / 28 材料|15.724|-0.125|13.46%|

无 RNA 的训练材料标签中位数基线，材料等权 MAE：内部 DDR 11.380、DNA repair 12.742；整来源 DDR 12.640、DNA repair 14.075。模型在内部留出明显优于基线，在整来源留出没有优于基线。误差阈值比例不是分类“正确率”，两个加权口径不可混用。

因此，本次回答老师的目标已经变为 DDR / DNA repair 范围内的乳酸化蛋白检出占比。已有材料内共享标签可较准确重现；新来源外推尚未得到支持。较差外推可能涉及来源、材料、检测深度及共享标签结构，不能单凭这些结果归因于批次效应。两目标高度相关，不能当成两份独立复制证据。

## 验收文件

- `outputs/20260927_ddr_fraction_labels/`：31组原始计数、28材料标签、逐蛋白分母/分子成员、GO集合和输入校验。
- `outputs/20260927_ddr_fraction_models/comparison_metrics.csv`：两目标、四评价、两加权口径共16行指标。
- 各模型 `figures/RF_mtry135.pdf` / `.png`：四联散点，横轴赋予RNA记录的材料占比、纵轴预测占比；每点一条RNA。依次训练内、内部留出、全量拟合、整来源留出。
- 各模型 `figures/top25.pdf` / `.png` 与 `tables/training_importance.csv`：Top25及全基因排序。
- 各模型 `tables/predictions.csv.gz`：逐条预测；`internal_test_by_material.csv`：按材料核查。
- `audit/20260927_ddr_fraction/external_screening.md`：新独立数据检索、候选与缺口。本次尚未完成新增独立外部测试。
- `outputs/20260927_ddr_fraction_models/validation.txt`：成员重计数、指标重算、冻结切分一致、模型重载、校验和均通过。4张 PNG 已目视核验，无文字裁切。

## 复现

在仓库根目录运行。保留冻结输入和本次快照。脚本要求新的输出目录，避免缓存/旧结果混入；以下正式路径已有结果，复现时改成新的目录。train.R 明确读取本次冻结 labels 路径；如重新计算标签，先显式更新该输入路径及其 SHA 记录。

```sh
Rscript --vanilla workflow/ddr_fraction/labels.R outputs/20260927_ddr_fraction_labels
Rscript --vanilla workflow/ddr_fraction/train.R outputs/20260926_full_rna_inputs outputs/20260927_ddr_fraction_models/DDR DDR
Rscript --vanilla workflow/ddr_fraction/train.R outputs/20260926_full_rna_inputs outputs/20260927_ddr_fraction_models/DNA_repair DNA_repair
Rscript --vanilla workflow/ddr_fraction/validate.R
python3 workflow/ddr_fraction/catalog_external.py
```

267 MB full_rna.rds 使用原本本地输入及已存 SHA256，不重复提交大矩阵。新模型、预测、图表和脚本单独版本化；没有改动旧模型、原始蛋白数据或用户 Word 文档。
