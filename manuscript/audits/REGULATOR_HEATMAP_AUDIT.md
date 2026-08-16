# 当前乳酸化调控蛋白热图审计

## 两份热图

1. Kla/lactylome热图：30个当前Kla样本组。
2. 普通全蛋白热图：与30组配对的28条唯一参照展示行。

两图按non-tumor tissues、tumor tissues、cancer cell lines、
normal cell lines排列。列按Writer、Eraser、Writer-Eraser和Reader分面。
英文和中文图使用相同数据顺序，字体为Arial Unicode MS。

## 数值定义

跨研究原始强度不可直接比较。每个样本内部先对有效信号执行`log2(signal+1)`
（来源已为log2时不重复转换），再在该样本全部有效特征中计算0–100平均秩百分位。
同组多个样本取中位数。

可定量样本中未检出的调控蛋白记为0，仅表示绘图缺失约定，不表示真实生物学零丰度。
Kla信号从不用于补普通全蛋白信号。

## 标签和对齐

- Kla行名描述Kla数据自身的材料和处理状态。
- 普通全蛋白行名独立描述参照数据自身的材料/状态，并附参照PXD。
- 共享同一普通全蛋白矩阵和相同子集的Kla组只显示一条参照行。
- 两图色块尺寸一致；分类间距用于容纳标签，不改变单个色块大小。

## 当前输出

- Kla热图输入：
  `results/tables/kla_regulator_normalized_intensity_long.csv`
- 普通全蛋白热图行：
  `results/tables/kla_regulator_whole_proteome_heatmap_rows.csv`
- 普通全蛋白算法审计：
  `results/tables/kla_regulator_whole_proteome_algorithm_audit.csv`
- 定量可用性审计：
  `results/tables/kla_regulator_intensity_availability_audit.csv`
- 图：
  `results/figures/kla_regulator_*_relative_intensity_heatmap_*.png`

热图表示各样本内部的相对位置，不是跨研究绝对丰度、fold change或差异分析。
