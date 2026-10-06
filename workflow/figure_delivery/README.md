# 乳酸化图片汇总与机器学习零起点重绘

当前入口（2026-10-06）：运行 `Rscript --vanilla workflow/figure_delivery/render_latest_ml_zero.R`，再运行 `python3 workflow/figure_delivery/package_latest_figures.py`（Python 只使用标准库）。输出位于独立的 outputs/20261006_figures_latest；交付目录为桌面“乳酸化全部图片_当前版本_20261006”，同时生成 ZIP。按选定面板打包，每张图有当前/历史标记、来源和 SHA256；版本选择依据表记录替换理由。

2026-10-04 的 `render_ml_zero.R` / `package_figures.py` 是历史交付入口，保留用于复现。其重绘丢失了最新图的基因符号、材料参考菱形及类别 UMAP 面板，也漏了当前四目标概览；本次恢复并补齐。现行四目标放主目录，跨组学两目标与全蛋白 Kla 占比移入明确标记的历史目录。蛋白组只保留标准面板编号，不重复收录别名；旧配套图单独归档。

R 读取已保存的预测表，重绘三版主要机器学习的 7 个模型，不训练模型，不修改旧结果。依赖 data.table、ggplot2、ragg、patchwork、uwot。Ensembl 为分析键，Symbol 仅为显示注释，直接读取冻结注释表并核对重要性 ID/数值；seed=25。现行模型为组学内比例四目标，历史模型为跨组学 DDR/DNA repair 两目标和全蛋白 Kla 占比。

非负数值轴通过 limits=c(0,NA)、下端 expansion=0 强制从零起点显示，并以 ggplot_build 核对面板下界。分类轴保留分类标签。绝对误差图与计数直方图从零开始，直方图核对全部记录均进入分箱。带符号残差单独作为诊断附录保留负数，避免丢失低估信息。UMAP 在相同 RNA 矩阵上按最新参数计算后分别减去两轴最小值，保存原始及平移坐标，核对差分不变；全数据拟合着色不作为独立验证。原图未保存坐标，因此不保证新嵌入与旧图片逐点一致。

正式蛋白组图片来自 corrected_final_result_20260905_sample_only 下 Figure_1–3 与 Supplementary_Figure_S1–S4，已与原 fork、RNA 工作树及 20260917 蛋白组复制版逐图比较；RNA 以 config/rna_current_release.csv 指定的 outputs/20260919_repair_release 为准；乳酸代谢来源 outputs/20260920_lactate_metabolism/figures。不收录旧试跑目录、重复生成器目录或错误外部 RNA 转换结果。机器学习展示内容依据 workflow/assay_composition/render_all_extended_figures.R、summarize.R、workflow/report_revision/render_extended.R / render_top25.R 与 workflow/full_rna_fraction/generate_extended_figures.R；回归重要性术语沿用科学修订版的“方差减少”。
