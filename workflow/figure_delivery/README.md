# 乳酸化图片汇总与机器学习零起点重绘

运行 `Rscript --vanilla workflow/figure_delivery/render_ml_zero.R`，再运行 `python3 workflow/figure_delivery/package_figures.py`（Python 只使用标准库）。输出位于独立的 outputs/20261004_all_figures，交付目录与 ZIP 位于桌面。

R 读取已保存的预测表，重绘三版主要机器学习的 7 个模型，不训练模型，不修改旧结果。依赖 data.table、ggplot2、ragg、uwot。Ensembl 为分析键，seed=25。现行模型为组学内比例四目标，历史模型为跨组学 DDR/DNA repair 两目标和全蛋白 Kla 占比。

非负数值轴通过 limits=c(0,NA)、下端 expansion=0 强制从零起点显示。分类轴保留分类标签。绝对误差图与计数直方图从零开始。带符号残差单独作为诊断附录保留负数，避免丢失低估信息。UMAP 在相同 RNA 矩阵上计算后分别减去两轴最小值，图轴注明平移坐标，点间距离不变；全数据拟合着色不作为独立验证。

正式蛋白组图片来自 corrected_final_result_20260905_sample_only 下 Figure_1–3 与 Supplementary_Figure_S1–S4；RNA 来源 outputs/20260919_repair_release；乳酸代谢来源 outputs/20260920_lactate_metabolism/figures。不收录旧试跑目录、重复生成器目录或错误外部 RNA 转换结果。图片索引记录原始来源与 SHA256；PNG、PDF 均保留。
