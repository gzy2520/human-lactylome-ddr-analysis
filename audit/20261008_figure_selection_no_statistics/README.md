# 无统计标注选图包

按最终选图前的展示需求，另行输出无生物学显著性标注的 111 对当前独立 PNG/PDF，不修改此前统计版本，不重新计算任何检验或模型。

13 张分类图从保存的 ggplot 对象中移除旧比较图层及总体检验副标题，去掉统计说明，收回比较横线占用的留白。两张已美化的 RNA 图使用 20261007 保存的无横线版本，其冻结数据与当前一致。其余 96 张当前图原样复制；旧组合图及历史模型目标不收录。

原始点、均值、中位数、箱体、误差线和样本量保留。机器学习的 MAE、R² 等模型评价保留。这里的“无统计标注”指移除生物学显著性比较展示，不代表源数据从未进行过分析，也不代表删去全部描述性统计量。

核对 13 张重绘图的原始数据层与汇总坐标一致、比较图层数为零、统计副标题与图注为空；完整包中 67 个非机器学习 PDF 检查没有检验 P/q、ANOVA、ns/NE 或 CR2 标注，乳酸代谢图中的大写 NS 材料缩写保留。完整包校验 111 对图片及 ZIP 内容 SHA256。

复现：

```sh
Rscript --vanilla workflow/figure_delivery/render_selection_no_statistics.R
python3 workflow/figure_delivery/package_selection_no_statistics.py
```

输入对象、原画布尺寸记录在 `render_plan.csv`；输出 `outputs/20261008_figure_selection_no_statistics`，桌面交付“乳酸化选图包_无统计标注_20261008”。
