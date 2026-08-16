# 当前普通全蛋白参照审计

## 当前结论

当前30个Kla样本组均具有材料身份可接受的非Kla富集普通全蛋白参照，并映射到
28条唯一热图展示行。参照关系并非全部为同一批样本或相同处理状态：

- 同一批实际样本；
- 同研究匹配样本；
- 同材料基线；
- 独立队列的相同材料或表型。

因此Methods和图注必须区分材料身份匹配与实验状态匹配，不能统一称为
“condition-matched controls”。

## 当前范围变化

原33个可配对组中，`PXD055230`、`PXD057709`和`PXD014870/MCF7`按老师要求
排除，形成当前30组。另有4个Kla定量组因缺少可接受普通全蛋白参照而未进入原33组；
3个PXD037371临床组因TMT通道无法可靠拆分而不属于37个可定量组。

## 可审计输出

- 英文材料身份审计：
  `results/tables/strict_reference_material_identity_audit.csv`
- 中文材料身份审计：
  `results/tables/strict_reference_material_identity_audit_zh.csv`
- 30组老师审阅表：
  `results/tables/kla_and_reference_teacher_review_zh.csv`
- 28条普通全蛋白展示行：
  `results/tables/kla_regulator_whole_proteome_heatmap_rows.csv`

每行保留Kla材料、参照材料、取材粒度、实际使用文件/列、共享参照关系、
`ExperimentalStateMatch`和`PairingCaveat`。
