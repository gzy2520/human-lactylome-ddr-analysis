# 图二 RNA 单图修订

本次按确认后的 RNA 指标提供两个独立面板，组合排版由论文作者完成。原截图左图是质谱 DDR 蛋白占比，本次左图改为已有 RNA 基因集表达结果；不能继续使用原截图的百分比纵轴或蛋白组图例。原蛋白组、RNA 冻结结果和 20261006 图片包均保留。

“图二”是本次论文排版编号，所用冻结 RNA 源图和表仍沿用 Figure_1a/1b 名称；与旧项目中编号为 Figure_2 的 UpSet/七通路图不是同一组面板。

## 两个面板

- `Figure_2a_RNA_DDR_and_Kla_gene_expression`：固定 DDR 基因集和已检出乳酸化蛋白对应基因集的 RNA 表达，蓝色为 DDR，橙色为 Kla 蛋白对应基因。读取冻结的 `figure1a_rna_ddr_and_lactylated_28materials.csv`，共 56 个基因集与材料组合点、8 个箱体；两类基因集均有 28 份 RNA 材料。
- `Figure_2b_RNA_proliferation`：196 个冻结 Hallmark G2M Checkpoint 基因的 RNA 增殖相关转录评分，读取 `figure1b_hallmark_proliferation_28materials.csv`，共 28 个材料点、4 个箱体。评分不等同于实测增殖速率。

四类依次为非肿瘤组织、肿瘤组织、正常细胞系、癌细胞系，材料数为 9、3、6、10。横轴 n 是 RNA 材料数，不能理解为独立患者数。相同 HCT116、HK-2 RNA 基线已去重复用，未回退到截图中的 31 行展开版本。

## 图形说明

每个点是一份材料的汇总值。左图先在每条 RNA 记录内取对应基因集的变换后表达中位数，再取材料内 RNA 记录的中位数。右图先在每条 RNA 记录内取 196 个基因的变换后表达均值，再取材料内 RNA 记录的中位数。变换为 qsmooth 后的 log2(TPM + 0.5)，与冻结处理一致。

箱体为四分位间距，黑线为中位数，须线按 1.5 倍四分位距规则；红线为类别内材料点的未加权均值。所有材料点均保留，包括箱体之外的点。两个基因集可能重叠，不是互斥集合。Kla 蛋白对应基因的转录表达不是乳酸化蛋白占比、修饰强度或占据率。

保持原有颜色和分类顺序，去除网格，细化点与箱体边线，缩小图例，把 n 合并到横轴标签。按最后确认要求，两个面板均不显示上方的组间比较括号、ns 或星号，并收回相应留白；箱内中位数和红色均值线保留。每张图为 3.8 × 3.45 英寸；PNG 为 600 dpi，另提供 PDF 和 SVG 矢量单图。未生成组合图。

## 比较与核对

左图不能复用蛋白组 P 值。本次对两个 RNA 基因集分别采用来源平衡类别均值模型、来源块聚类的 CR2 标准误及 Satterthwaite 自由度。四个预设比较为非肿瘤组织分别与其他三类比较，以及正常细胞系与癌细胞系比较；两个基因集共 8 个比较作 BH 校正。统计结果仅随 CSV 提供，图中不标注。

右图统计表沿用 20261006 已保存的 4 个 RNA 增殖评分来源聚类比较，P 和 q 未改动。两图的预设比较均未通过校正阈值（q ≥ 0.05）。来源平衡检验与未加权箱体描述的均值口径不同；来源聚类不能消除处理条件或跨研究混杂。

核对包括材料去重及数量、绘制点数、点值与冻结数据的一致性、箱体数量、无网格、原始表文件一致性、RNA 增殖评分比较未变以及 PDF/PNG/SVG 文件。全部检查通过；两个矢量 PDF 已渲染检查，标题、刻度与全部材料点可读。

## 复现

```sh
R_LIBS_USER=/tmp/kla_teacher_fig_Rlib Rscript --vanilla workflow/figure_delivery/render_rna_figure2_refined.R
python3 workflow/figure_delivery/package_rna_figure2_refined.py
```

R 依赖 data.table、ggplot2、ragg、clubSandwich；本次 clubSandwich 0.7.0 和 sandwich 3.1-3 使用既有临时 R 库。输入版本为 `outputs/20260919_repair_release/core`，来源块元数据为 `outputs/20260925_sample_fraction_inputs_final/sample_metadata.csv`。输出为 `outputs/20261007_rna_figure2_refined`，未重跑 RNA 矩阵或机器学习。
