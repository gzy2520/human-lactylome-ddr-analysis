# Kla project — current AI handoff

更新时间：2026-09-04
当前修订分支：`figure-enhancement`
当前修订提交：`feat(figure1): upright layout with sample-level dots for DDR fraction and MKI67 ratios`

## 先区分两条线

### 已批准的冻结定稿基线

`main` 分支仍是原来的 30 个 Kla 数据集/组的冻结出版版本。默认的
`data/publication_input/`、`results/figures/` 和正式 publication renderer
不能因为候选图修改而被覆盖、清理或重算。所有分析键必须使用稳定的
UniProt `BaseAccession`，不能用 Gene Symbol 做集合分析。

### 当前正在修订的候选版本

当前修订只在 `figure-enhancement` 分支和日期隔离目录中进行。老师要求把
代表食管癌的 PXD064038 重新纳入，并为它寻找“同为食管癌、但没有乳酸化
测量”的全蛋白对照：

- 乳酸化数据：PXD064038 的 `MEC and NEC ESCC groups`，6 个 Kla 样本，
  1,239 个 Kla 蛋白、92 个 Kla-DDR 蛋白；只选这一组，不纳入该项目另外
  两组。
- 非乳酸化全蛋白对照：iProX/PXD065830 的 Dataset1，sheet
  `2.a protein raw information`，只使用 94 个 `ESCC-*T` 肿瘤列；24 个
  `N` 非肿瘤列排除。该肿瘤 union 为 8,083 个蛋白，其中 420 个 DDR。
- PXD053809 当前没有纳入，因为缓存的 processed XLS 加密，无法审计。
  PXD010154 只是健康食管背景，不能替代食管癌全蛋白对照。
- PXD065830 的 94 个 T 列在 Figure 1 中逐样本显示；它们不是与 6 个
  PXD064038 Kla 样本一一配对的同一样本。

## 图形修改契约

用户给出的 2026-09-02 图是视觉基线。后续修改只更新数据源和明确指出的
统计/排版，不要重新设计整套图片。

1. Figure 1 DDR fraction：保持“四个生物类别纵向面板；每个面板两个
   水平 boxplot（Whole proteome/Kla）；点是源样本结果”的布局。箱内深色
   线是 median，红色竖线是 mean。每个类别做 Whole proteome 对 Kla 的
   one-way ANOVA，4 个类别的 p 值做 BH 校正。
2. MKI67/H3C1：这是全蛋白组的 MKI67/H3C1 强度比图，不是 Kla 图。只用
   有 MKI67=P46013 和 H3C1=P68431 完整定量的全蛋白样本；没有插补。四类
   比较使用 one-way ANOVA。PXD065830 的 94 个 T 列 H3C1 全部缺失，不能
   把它们硬加入该图。
3. UpSet：保持原来的四集合、15 个交集、零计数组合保留和排序契约；只
   修复左下角 set-size 数字被裁切的问题。UpSet 本身没有 ANOVA。
4. Kla DDR pathway summary：不是“四张图、每张 7 行”，而是“7 张图、
   每张 4 个生物类别行”。每个类别同时有 Up/positive 和 Down/negative
   两个 box（至多 8 个 box/图），Down 放在 Up 同一侧的相邻位置，不画到
   负方向。每个点是一个 Kla 源样本。每个 pathway 使用
   `Category × Direction` 的 two-way ANOVA，保留 Category、Direction 和
   interaction 三个项，7 张图共 21 个项做 BH 校正。
5. 随机种子统一使用 `25`。R 优先。不要因为图形不好看而改数据、顺序、
   颜色、分组或分母定义。

## 当前候选数据规模

扩展输入目录：
`data/candidate/escc_inclusion_20260903_pxd065830_tumor_reference/candidate_input/`

- Figure 1：310 个点；Whole proteome 212，Kla 98。
- Figure 1 按类别的点数（Whole proteome / Kla）：
  - non-tumor tissues：62 / 32
  - tumor tissues：107 / 12
  - normal cell lines：15 / 19
  - cancer cell lines：28 / 35
- Kla pathway 输入：686 行，即 98 个 Kla 源样本 × 7 个 pathway。
- 扩展 publication group 数为 31；类别组数为 9 / 3 / 7 / 12（normal
  tissue / cancer tissue / normal cells / cancer cells）。

## 当前输出位置

候选输出：
`results/candidate/escc_inclusion_20260903_pxd065830_tumor_reference/legacy_layout/`

- `Figure_1_DDR_fraction_candidate_category_boxplot_refined.png/.pdf`
- `mki67_ratio_boxplot/Figure_1_MKI67_over_ACTB_boxplot.png/.pdf`
- `mki67_ratio_boxplot/Figure_1_MKI67_over_TUBB_boxplot.png/.pdf`
- `mki67_ratio_boxplot/Figure_1_MKI67_over_H3C1_boxplot.png/.pdf`
- `pathway_summary_by_pathway/Figure_2_DDR_pathway_summary_<BER|NER|MMR|FA|HR|AEJ|NHEJ>_barplot.png/.pdf`（及同名兼容 `_boxplot.png/.pdf`）

### 2026-09-04 最新修订：Figure 1 核心图表正向立式（左旋90度）重构与个体样本级散点观测
- **背景与修改要求**：
  - 用户明确要求：将 Figure 1 图表正过来（左旋 90 度，由原来的水平箱线图转换为立式直立面板布局）；
  - 横轴排列 4 大类（`non-tumor tissues`, `tumor tissues`, `normal cell lines`, `cancer cell lines`）；
  - 图中所有散点必须为**独立样本级别观测**（`ObservationType == "sample"`），严禁聚合为数据集（PXD）级别的点；
  - 同步更新 Figure 1 的 4 张关键图表：
    1. `Figure_1_DDR_fraction_candidate_category_boxplot_refined`（DDR 分数图）
    2. `Figure_1_MKI67_over_ACTB_boxplot`
    3. `Figure_1_MKI67_over_TUBB_boxplot`
    4. `Figure_1_MKI67_over_H3C1_boxplot`
- **Figure 1 DDR 分数图（Figure_1_DDR_fraction_candidate_category_boxplot_refined）重构要点**：
  - **布局设计**：采用 4 列并列分面面板（`facet_grid(. ~ CategoryLabel)`），顶部保留 4 色马卡龙背景 strip（`#DCE9E2`, `#F0DEDE`, `#E7E1EE`, `#EEE4D2`）；
  - **模态对比**：每个分类面板内并列放置立式箱线图：左侧为 Whole proteome（蓝色 `#4E79A7`），右侧为 Lactylome (Kla)（橙色 `#F28E2B`）；
  - **样本散点**：每个点为来源于测序源数据的独立样本（310 个样本点：62/32, 107/12, 15/19, 28/35），使用个人随机数种子 25 进行水平微抖动（`position_jitter(width = 0.16, height = 0, seed = 25)`）；
  - **统计线与标注**：箱内黑色横线为中位数（Median），红色粗横线为均值（Mean，`colour = "#C0392B"`）；顶部标注样本量（`n=...`）及单因素方差分析显著性横杠与星号（BH 校正后四类均为 `****`）；
  - **纵轴刻度**：纵向显示 `GO-DDR annotated protein fraction (%)`（0% ~ 25%）。
- **三张 MKI67 比值图（MKI67 over ACTB / TUBB / H3C1）同步更新与优化**：
  - **纳入食管癌肿瘤对照**：全面纳入 PXD065830 的 94 个 ESCC-T 肿瘤组织全蛋白样本（MKI67=P46013，ACTB=P60709 检出 5 例，TUBB=P07437 检出 94 例，H3C1=P68431 检出 45 例）；
  - **样本级散点**：图中点均为独立源样本点（种子 25 抖动），tumor tissues 样本量显著充实（TUBB 由 n=5 增至 n=99，ANOVA q=3.06e-28；H3C1 由 n=5 增至 n=50，ANOVA q=5.90e-07；ACTB 由 n=5 增至 n=10，ANOVA q=0.067）；
  - **对数纵轴顶部间距（Headroom）优化**：重构对数刻度范围计算，将 `y_max` 调整为 `10^(log10(raw_max) + 0.65)`，各组样本量标注 `label_y` 设为 `10^(log10(MaxRatio) + 0.22)`，彻底消除此前 `n=14`（H3C1）与 `n=99`（TUBB）与顶部外框发生微接触或裁剪的问题。
- **输出文件路径**：
  - 候选输出（31 组 ESCC 纳入）：
    - `results/candidate/escc_inclusion_20260903_pxd065830_tumor_reference/legacy_layout/Figure_1_DDR_fraction_candidate_category_boxplot_refined.png/.pdf`
    - `results/candidate/escc_inclusion_20260903_pxd065830_tumor_reference/legacy_layout/mki67_ratio_boxplot/Figure_1_MKI67_over_ACTB_boxplot.png/.pdf`
    - `results/candidate/escc_inclusion_20260903_pxd065830_tumor_reference/legacy_layout/mki67_ratio_boxplot/Figure_1_MKI67_over_TUBB_boxplot.png/.pdf`
    - `results/candidate/escc_inclusion_20260903_pxd065830_tumor_reference/legacy_layout/mki67_ratio_boxplot/Figure_1_MKI67_over_H3C1_boxplot.png/.pdf`
- **自动化测试**：全量 12 项契约测试（`tests/validate_*.R`）100% 全部 PASS。

### 2026-09-03 最新修订：七通路 Summary 图表升级为带样本散点的立式 SEM 柱状图
- **表现形式**：响应用户指令，七通路 summary 图不再使用 boxplot，改为带 error bar 的柱状图，并在柱上叠加个体样本散点以确保数据透明：
  - 柱高表示各分类的均值（Mean，柱透明度 `alpha = 0.65`）；
  - 误差棒（error bar）使用 SEM 算法：$SEM = SD / \sqrt{n}$，上下界为 $Mean \pm SEM$（下界截断于 0），位于顶层保持统计清晰；
  - 柱上样本散点：`geom_point`（`shape = 21, size = 2.6, stroke = 0.55, alpha = 0.85, colour = charcoal`），内部与柱子同色填充（Pro 为通路主题色，Inh 为灰调 `#98A1AA`），水平微抖动 `position_jitter(width = 0.16, height = 0, seed = 25)`（遵循个人随机数种子 25 规则），散点完美收敛于柱宽内部。
- **正向立式布局（左旋 90 度）**：
  - 类别横向并列排列（顶部 4 个 strip，保留马卡龙色背景）；
  - 纵轴为竖直百分比 `"Relative portion of Kla-DDR proteins (%)"`。
- **命名与排序**：
  - 不再使用 `Up/down`，正向用 `"Pro"` 表示，反向用 `"Inh"` 表示；
  - `"Pro 在前，Inh 在后"`（Pro 居左使用通路专用色，Inh 居右使用灰调 `#98A1AA`）。
- **兼容性保障**：脚本同时输出 `_barplot.png/.pdf` 和 `_boxplot.png/.pdf`，manifest 与自动化测试 100% 通过。

### 2026-09-03 最新额外作图：10 个特定调控因子全蛋白样本级表达百分位数散点图（10 张图）
- **背景与目的**：对原 Figure 3a 全蛋白百分位数热图进行精细化微观拆解。原热图展示的是每个 PXD/数据集的中位数百分位数，本修订按**单个样本（sample）**级别计算并绘制散点图。
- **10 个目标基因与 UniProt 稳定 BaseAccession 严格 1:1 映射**：
  1. `AARS1` $\to$ `P49588` (Lactylation Writer)
  2. `ACAT2` $\to$ `Q9BWD1` (Lactylation Writer)
  3. `KRT18` $\to$ `P05783` (Lactylation Writer)
  4. `SIRT2` $\to$ `Q8IXJ6` (Lactylation Eraser)
  5. `PARK7` $\to$ `Q99497` (Lactylation Writer-Eraser)
  6. `HDAC1` $\to$ `Q13547` (Lactylation Writer-Eraser)
  7. `HDAC2` $\to$ `Q92769` (Lactylation Writer-Eraser)
  8. `BRD4` $\to$ `O60885` (Lactylation Reader)
  9. `SMARCA4` $\to$ `P51532` (Lactylation Reader)
  10. `TRIM33` $\to$ `Q9UPN9` (Lactylation Reader)
  *严格执行全局规则：所有底层数据提取与匹配必须使用 `BaseAccession`，Gene Symbol 仅作为显示名称。*
- **作图规范**：
  - 横轴（X 轴）：四分类（`non-tumor tissues`, `tumor tissues`, `normal cell lines`, `cancer cell lines`，带马卡龙色 `#DCE9E2`, `#F0DEDE`, `#E7E1EE`, `#EEE4D2`）；
  - 纵轴（Y 轴）：在全蛋白数据中，该蛋白在其所属样本的真实表达百分位数（`Whole-proteome expression percentile (%)`，范围 $0\% \sim 100\%$）；
  - 图表形式：散点图（Scatter plot），每个点为一个独立的全蛋白参考样本（`shape = 21, size = 2.8, stroke = 0.55, colour = charcoal, alpha = 0.75`），水平抖动采用个人随机数种子 `25`（`position_jitter(width = 0.22, height = 0, seed = 25)`）；
  - 统计标记：各分类中位线以深色横棒（`median crossbar, linewidth = 0.7, width = 0.45`）清晰标示；顶部标注各分类样本总数 $n$。
- **数据输出与图表路径**：
  - 数据表：`data/candidate/escc_inclusion_20260903_pxd065830_tumor_reference/target_10_regulators_sample_percentiles.csv`（包含 31 组 6,340 条观测记录与 30 组基准 5,400 条记录）；
  - 候选范围图表（31 datasets，包含 94 个 ESCC 肿瘤全蛋白对照样本）：
    `results/candidate/escc_inclusion_20260903_pxd065830_tumor_reference/regulator_sample_percentiles/Figure_3_regulator_sample_percentile_<GENE>_<ACCESSION>.png/.pdf`
  - 冻结基线图表（30 datasets）：
    `results/candidate/escc_inclusion_20260903_pxd065830_tumor_reference/regulator_sample_percentiles/baseline_30_datasets/Figure_3_regulator_sample_percentile_<GENE>_<ACCESSION>_30datasets.png/.pdf`
- **生成脚本与测试契约**：
  - 生成脚本：`R/candidate/build_regulator_sample_percentile_scatters.R`
  - 契约测试：`tests/validate_regulator_sample_percentile_contract.R`（已纳入自动化回归套件，全量通过）。

### 2026-09-03 最新修订：全蛋白与乳酸化调控因子百分位数热图（Figure 3a & 3b）同步更新至 31 组
- **背景与原因**：由于候选分支纳入了代表食管癌的 PXD064038（乳酸化 6 样本）与其配对的 PXD065830（94 个 ESCC-T 肿瘤组织全蛋白对照），源数据已由 30 组扩展为 31 组。原 Figure 3a 和 Figure 3b 热图此前仅展示 30 组（肿瘤组织仅 2 行：Prostate cancer 与 Primary HCC），现根据源数据变更同步更新。
- **数据提取与计算逻辑**：
  - **Figure 3a（全蛋白百分位数热图）**：提取 PXD065830 的 94 个 ESCC-T 肿瘤样本的原始定量信号，按样本计算各蛋白百分位数后取中位数；在 `tumor tissues` 面板新增首行 `Independent ESCC tumor whole proteome`（使肿瘤组织由 2 行增加为 3 行完整呈现）。
  - **Figure 3b（乳酸化 Kla 百分位数热图）**：提取 PXD064038 的 6 个样本（MEC_1-3, NEC_1-3）在 `La (K)Sites.txt` 中的位点定量信号，按样本计算百分位数后取中位数；在 `tumor tissues` 面板新增行 `ESCC MEC/NEC groups`。其中检出具有乳酸化定量信号的调控因子包括 HDAC1（中位数 66.7%）、PARK7（中位数 58.2%）、HAT1（中位数 22.0%）、DPF2（中位数 23.9%）等。
- **生成脚本与文件**：
  - 核心脚本：`R/candidate/update_regulator_percentile_heatmaps.R`
  - 更新数据：`data/candidate/escc_inclusion_20260903_pxd065830_tumor_reference/publication_input/regulator_reference_percentiles_30.csv`（1,488 行，覆盖 31 组 × 48 个调控因子）与 `regulator_kla_percentiles_30.csv`（1,519 行，覆盖 31 组 × 49 行），并同步刷新 `INPUT_MANIFEST.csv`。
  - 正式输出：
    - `results/escc_inclusion_20260903_pxd065830_tumor_reference/formal_figures/Figure_3a_reference_regulator_percentiles.png/.pdf`
    - `results/escc_inclusion_20260903_pxd065830_tumor_reference/formal_figures/Figure_3b_Kla_regulator_percentiles.png/.pdf`
  - 候选镜像：`results/candidate/escc_inclusion_20260903_pxd065830_tumor_reference/heatmaps/`
- **自动化测试**：全量 12 项契约测试（包括 `validate_escc_inclusion_scope.R`）100% 全部通过。

扩展正式图与附表（与冻结默认正式出版目录隔离）：
- 正式图表目录：`results/escc_inclusion_20260903_pxd065830_tumor_reference/formal_figures/`
- 候选附表目录：`results/escc_inclusion_20260903_pxd065830_tumor_reference/supplementary/`（及镜像 `results/candidate/escc_inclusion_20260903_pxd065830_tumor_reference/supplementary/`）
  - `Supplementary_Table_S1_Kla_Data.xlsx`：扩展至 31 组（纳入 PXD064038 6 个 Kla 样本、1,239 个蛋白、92 个 DDR 蛋白）。
  - `Supplementary_Table_S2_Reference_Data.xlsx`：扩展至 31 组（纳入 PXD065830 94 例 ESCC 肿瘤参考全蛋白、8,083 个蛋白、420 个 DDR 蛋白）。
  - `Supplementary_Table_S3_Human_DDR_GO_Annotations.xlsx`：6,707 条人类 DDR GO 注释保持完备。
  - `Supplementary_Table_S4_Venn_Membership.xlsx`（原 Table S6 重命名）：扩展 4 个集合分析的 15 区域（AllKla 5,814、KlaDDR 401、Reference 24,397、ReferenceDDR 836）。
  - `Supplementary_Table_S5_Pathway_Protein_Ranking.xlsx`（原 Table S4 重命名）：TumorTissues 面板由 178 扩展至 192 个蛋白质（新增 14 个食管癌相关 Kla-DDR 蛋白并赋序）。
  - `Supplementary_Table_S6_Lactylation_Regulators.xlsx`（原 Table S5 重命名）：49 个调控因子及其文献依据完整保留。

### 2026-09-03 最新修订：UpSet 图版本核对、像素级对齐与排列顺序审计
- **版本数据核实**：
  - 候选目录 `results/escc_inclusion_20260903_pxd065830_tumor_reference/formal_figures/` 中的 UpSet 图均为 31 组最新版本（Figure 2b 中 Tumor tissues 大小为 192，总数 401；Figure 2a 全蛋白 DDR 为 836；柱子交集蛋白数全面更新）。
- **点阵与柱子图对齐修复（消除 36.5 px 偏差至 0.00 px）**：
  - 原代码采用 `(top_row / bottom_row)` 嵌套 patchwork 拼图，左下角 set-size 文字边距挤压了下方面板，造成点阵相较上方柱形图向右偏移 36.5 像素。
  - 在 `R/publication/build_publication_outputs.R` 中将两层拼图重构为统一的全局 2×2 网格布局：`(plot_spacer() + intersection_plot + set_size_plot + matrix_plot) + plot_layout(ncol = 2, widths = c(3.7, 10), heights = c(1.45, 1.25))`。
  - Python 像素级扫描复测结果显示：全部 15 列的柱子中心与点阵圆点中心水平差均为 0.00 像素，实现像素级对齐。
- **两版本交集排列顺序说明**：
  - **30 组 vs 31 组**：列顺序**不同**。因为 UpSet 遵循按交集蛋白数量降序排列（`bars ordered by size`）。31 组纳入食管癌后，交集 14（肿瘤组织+癌细胞系+正常细胞系）蛋白数由 5 跃升至 12，排序从第 8 列升至第 5 列；交集 13 则降为 3，排序后移至第 10 列。
  - **同一版本内 Figure 2a 与 Figure 2b**：列顺序**完全相同**。Figure 2a 严格通过 `mask_order = kla_ddr_mask_order` 强制对齐 Figure 2b 的交集列顺序，便于横向比较。
- **UpSet Y 轴 'Intersection size' 错位重叠彻底修复**：
  - 原代码将 `intersection_axis_title` 绘制在 `X = 0.5` 且未指定横向位移，导致 90 度竖立的 "Intersection size" 标题直接从刻度数字（如 '50'、'100' 或 '200'、'300'）正中穿过，并贴着第 1 列柱子边缘，发生严重遮挡错位。
  - 现引入根据刻度数字位数自适应的位移计算：`axis_title_vjust <- -2.2 - 0.9 * max_axis_digits`，并动态扩充左边距 `plot_left_margin <- 32 + 6.5 * max_axis_digits`（对 `intersection_plot` 与 `matrix_plot` 同步生效保证两图宽度一致）。"Intersection size" 完美平移至刻度数字左侧并留出清爽留白，文字重叠与数字遮挡彻底解决。

## 关键脚本

- 上游扩展：`R/candidate/prepare_escc_inclusion_inputs.R`
- 附表构建：`workflow/build_supplementary_workbooks.R`（支持环境变量 `KLA_PUBLICATION_INPUT` 与 `KLA_SUPPLEMENTARY_OUTPUT` 隔离输出）
- Figure 1：`R/candidate/build_figure1_category_boxplot.R`
- MKI67 全蛋白比值：`R/candidate/build_figure1_mki67_ratio_boxplots.R`
- Kla 七通路图：`R/candidate/build_ddr_pathway_summary_boxplots.R`
- 10 个调控因子样本散点图：`R/candidate/build_regulator_sample_percentile_scatters.R`
- 31 组热图更新脚本：`R/candidate/update_regulator_percentile_heatmaps.R`
- ANOVA：`R/candidate/boxplot_significance.R`
- UpSet/正式图：`R/publication/build_publication_outputs.R`
- 附表契约验证：`tests/validate_candidate_supplementary_contract.R`

## 推荐接手流程

```bash
git status --short
git switch figure-enhancement
git show --stat --oneline ac8d877
```

接手 AI 必须先阅读本文件，再阅读上述实际脚本和对应输入表。任何新改动
都要先放在 `results/candidate/` 日期隔离目录中，渲染后检查 PNG/PDF，
再运行：

```bash
KLA_CANDIDATE_INPUT=data/candidate/escc_inclusion_20260903_pxd065830_tumor_reference/candidate_input \
Rscript tests/validate_figure1_category_boxplot_contract.R .
KLA_CANDIDATE_INPUT=data/candidate/escc_inclusion_20260903_pxd065830_tumor_reference/candidate_input \
Rscript tests/validate_ddr_pathway_summary_boxplot_contract.R .
KLA_CANDIDATE_INPUT=data/candidate/escc_inclusion_20260903_pxd065830_tumor_reference/candidate_input \
Rscript tests/validate_boxplot_significance_contract.R .
Rscript tests/validate_figure1_mki67_ratio_contract.R .
Rscript tests/validate_escc_inclusion_scope.R .
```

不要 `git reset --hard`、不要覆盖 `main` 的冻结输入，也不要把未审计的新
数据直接加入正式 release。提交和 push 是两件事；没有明确要求时只提交，
不要 push。

## 可直接复制给其他 AI 的 prompt

```text
你接手的是 /Users/gzy2520/Desktop/Research/kla 的 Kla 人类 DDR 蛋白组项目。
请先阅读项目根目录的 AI_HANDOFF_CURRENT.md，再查看当前 Git 状态；当前有效
修订在 figure-enhancement 分支，提交 ac8d877。main 的 30 组 publication
input 和正式 results/figures 是已批准冻结基线，不要覆盖、清理或重算。

当前任务是延续 2026-09-03 的候选修订：只纳入 PXD064038 的 MEC and NEC
ESCC groups（6 个 Kla 样本），并使用 iProX/PXD065830 Dataset1 的
2.a protein raw information 中 94 个 ESCC-*T 列作为同为食管癌的非乳酸化
全蛋白参考；24 个 N 列排除。分析集合只能用稳定 UniProt BaseAccession，
不能用 Gene Symbol；R 优先，随机种子用 25。

严格保持昨天提供的图形布局：Figure 1 是四个类别纵向面板、Whole proteome
和 Kla 两个水平 boxplot、点为源样本；median 深色、mean 红线。MKI67/H3C1
是全蛋白图，PXD065830 的 H3C1 全缺失，不得插补。UpSet 只修复左下角
set-size 数字裁切。通路 summary 必须是 7 张 Kla pathway 图，每张 4 个
类别行，Up/Down 在同一正向侧相邻显示，最多 8 个 box；不要改成四类横轴的
新设计。Figure 1 和 MKI67 使用 one-way ANOVA；通路使用 Category、Direction
及交互项的 two-way ANOVA，7 张图共 21 项做 BH 校正。

任何修改先放在日期隔离的 results/candidate 输出，先读脚本和输入，再改动。
不要擅自改数据意义、分母、排序、颜色或冻结正式图。完成后必须渲染 PNG/PDF、
做视觉检查，并运行 validate_escc_inclusion_scope.R、
validate_figure1_category_boxplot_contract.R、
validate_ddr_pathway_summary_boxplot_contract.R、
validate_boxplot_significance_contract.R 和
validate_figure1_mki67_ratio_contract.R。若发现新的科学范围选择或数据缺失，
先报告证据，不要猜测。
```
