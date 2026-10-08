# 乳酸化图片汇总与机器学习零起点重绘

2026-10-08 无统计标注选图版：`render_selection_no_statistics.R` 与 `package_selection_no_statistics.py` 输出 111 对当前独立 PNG/PDF，去除生物学显著性横线、星号、ns/NE、P/q 值及检验副标题；保留数据点、误差线和机器学习评价指标。两张 RNA 主图采用最新无横线美化版，不覆盖此前统计版本。详见 `audit/20261008_figure_selection_no_statistics/README.md`。

2026-10-08 全套分类统计补齐：`render_all_category_statistics.R` 在此前两张 RNA 单图之外重绘其余 13 张四类别 barplot/boxplot，包括蛋白组 DDR 占比、七个 Kla 通路、三个 MKI67 比值、RNA DDR 表达与检出基因占比。`package_all_category_statistics.py` 交付 111 对当前单图和 15 对统计图便捷包，逐图注明处理范围；材料明细图保留描述性展示，未新增材料级检验。无组合图，旧版本保留。详见 `audit/20261008_all_category_statistics/README.md`。

2026-10-08 图二统计标注修订：两个 RNA 单图均以肿瘤组织为参照，分别与其余三类比较；左图分别在 DDR 和 Kla 蛋白对应基因集内检验，共 6 条彩色比较横线，右图为 3 条横线。沿用来源平衡、来源聚类 CR2 检验，BH 分别校正左图 6 次和右图 3 次检验，全部如实标为 ns。运行 `render_rna_figure2_tumor_reference.R` 和 `package_rna_figure2_tumor_reference.py`；新输出为 `outputs/20261008_rna_figure2_tumor_reference`，详见 `audit/20261008_rna_figure2_tumor_reference/README.md`。无网格、不组合，旧版保留。

2026-10-07 图二单图修订：本次“图二”指用户截图中的左右箱线面板，最终确认两个面板均使用 RNA。左图为 DDR 与乳酸化蛋白对应基因的转录表达，右图为 Hallmark G2M RNA 增殖评分；使用 28 份去重材料。运行 `render_rna_figure2_refined.R` 后运行 `package_rna_figure2_refined.py`，输出至 `outputs/20261007_rna_figure2_refined` 和桌面“图二_RNA单图美化_20261007”。仅交付两个独立单图，均无网格、无上方比较括号或 ns/星号；详见 `audit/20261007_rna_figure2_refined/README.md`。此编号来自当前论文排版，冻结 RNA 源文件仍叫 Figure_1a/1b，不要与此前编号为 Figure_2 的 UpSet/七通路图混淆。

2026-10-06 教师修订交付：`outputs/20261006_teacher_figures` 与桌面“乳酸化全部图片_无网格与统计标注_20261006”在以上版本基础上去除全部背景网格，补充来源聚类的组间比较，生成正文图 1–3 及可编辑 ML 表。方法和版本说明见 `audit/20261006_teacher_figures/REPORT.md`，图注单列保存。旧输出保留。

复现依次运行：

```sh
R_LIBS_USER=/tmp/kla_teacher_fig_Rlib Rscript --vanilla workflow/figure_delivery/teacher_pairwise.R
KLA_ML_FIGURE_OUTPUT=outputs/20261006_teacher_figures/ML KLA_UMAP_CACHE=outputs/20261006_figures_latest/ML/umap_coordinates.csv.gz KLA_FIGURE_NO_GRID=1 KLA_SAVE_PLOTS=1 Rscript --vanilla workflow/figure_delivery/render_latest_ml_zero.R
Rscript --vanilla workflow/figure_delivery/render_teacher_nogrid.R
Rscript --vanilla workflow/figure_delivery/assemble_teacher_main.R
Rscript --vanilla workflow/figure_delivery/qa_teacher_release.R
python3 workflow/figure_delivery/package_teacher_figures.py
```

新增统计依赖 clubSandwich（本次 0.7.0）、sandwich（3.1-3），本次安装在临时 R 库 `/tmp/kla_teacher_fig_Rlib`，不修改系统库；库不存在时需先安装。矩阵默认只读使用原 RNA 工作树，`KLA_RNA_MATRIX_ROOT` 可指向另一份已核对的冻结矩阵目录。`KLA_TEACHER_OUTPUT` 可指定新输出目录；若完整重绘旧目录，应先移除其 `render_completed` 标记和已有 RNA/lactate 渲染输出，优先使用新目录。保存的绘图对象仅为本地中间文件，不进入桌面包或 Git。

正文 Word 表用已加载的 Codex 文档运行时执行 `build_ml_body_table.py`；该脚本只读取 `Table_1_ML_performance.csv`。渲染使用文档技能的 `render_docx.py`，在 macOS 需让 FONTCONFIG_FILE/FONTCONFIG_PATH 指向可见的系统字体配置；本次为 `/opt/homebrew/etc/fonts/fonts.conf` 与 `/opt/homebrew/etc/fonts`。表为黑白、可编辑格式，已逐页检查。

本次无网格修订前的汇总入口：运行 `Rscript --vanilla workflow/figure_delivery/render_latest_ml_zero.R`，再运行 `python3 workflow/figure_delivery/package_latest_figures.py`（Python 只使用标准库）。输出位于独立的 outputs/20261006_figures_latest；交付目录为桌面“乳酸化全部图片_当前版本_20261006”，同时生成 ZIP。按选定面板打包，每张图有当前/历史标记、来源和 SHA256；版本选择依据表记录替换理由。

2026-10-04 的 `render_ml_zero.R` / `package_figures.py` 是历史交付入口，保留用于复现。其重绘丢失了最新图的基因符号、材料参考菱形及类别 UMAP 面板，也漏了当前四目标概览；本次恢复并补齐。现行四目标放主目录，跨组学两目标与全蛋白 Kla 占比移入明确标记的历史目录。蛋白组只保留标准面板编号，不重复收录别名；旧配套图单独归档。

R 读取已保存的预测表，重绘三版主要机器学习的 7 个模型，不训练模型，不修改旧结果。依赖 data.table、ggplot2、ragg、patchwork、uwot。Ensembl 为分析键，Symbol 仅为显示注释，直接读取冻结注释表并核对重要性 ID/数值；seed=25。现行模型为组学内比例四目标，历史模型为跨组学 DDR/DNA repair 两目标和全蛋白 Kla 占比。

非负数值轴通过 limits=c(0,NA)、下端 expansion=0 强制从零起点显示，并以 ggplot_build 核对面板下界。分类轴保留分类标签。绝对误差图与计数直方图从零开始，直方图核对全部记录均进入分箱。带符号残差单独作为诊断附录保留负数，避免丢失低估信息。UMAP 在相同 RNA 矩阵上按最新参数计算后分别减去两轴最小值，保存原始及平移坐标，核对差分不变；全数据拟合着色不作为独立验证。原图未保存坐标，因此不保证新嵌入与旧图片逐点一致。

正式蛋白组图片来自 corrected_final_result_20260905_sample_only 下 Figure_1–3 与 Supplementary_Figure_S1–S4，已与原 fork、RNA 工作树及 20260917 蛋白组复制版逐图比较；RNA 以 config/rna_current_release.csv 指定的 outputs/20260919_repair_release 为准；乳酸代谢来源 outputs/20260920_lactate_metabolism/figures。不收录旧试跑目录、重复生成器目录或错误外部 RNA 转换结果。机器学习展示内容依据 workflow/assay_composition/render_all_extended_figures.R、summarize.R、workflow/report_revision/render_extended.R / render_top25.R 与 workflow/full_rna_fraction/generate_extended_figures.R；回归重要性术语沿用科学修订版的“方差减少”。
