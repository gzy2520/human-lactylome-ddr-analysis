# Kla project — current AI handoff

更新时间：2026-09-03
当前修订分支：`figure-enhancement`
当前修订提交：`27ec757 feat(regulators): generate sample-level whole-proteome percentile scatter plots for 10 target genes`

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

扩展正式图与附表（与冻结默认正式出版目录隔离）：
- 正式图表目录：`results/escc_inclusion_20260903_pxd065830_tumor_reference/formal_figures/`
- 候选附表目录：`results/escc_inclusion_20260903_pxd065830_tumor_reference/supplementary/`（及镜像 `results/candidate/escc_inclusion_20260903_pxd065830_tumor_reference/supplementary/`）
  - `Supplementary_Table_S1_Kla_Data.xlsx`：扩展至 31 组（纳入 PXD064038 6 个 Kla 样本、1,239 个蛋白、92 个 DDR 蛋白）。
  - `Supplementary_Table_S2_Reference_Data.xlsx`：扩展至 31 组（纳入 PXD065830 94 例 ESCC 肿瘤参考全蛋白、8,083 个蛋白、420 个 DDR 蛋白）。
  - `Supplementary_Table_S3_Human_DDR_GO_Annotations.xlsx`：6,707 条人类 DDR GO 注释保持完备。
  - `Supplementary_Table_S4_Venn_Membership.xlsx`（原 Table S6 重命名）：扩展 4 个集合分析的 15 区域（AllKla 5,814、KlaDDR 401、Reference 24,397、ReferenceDDR 836）。
  - `Supplementary_Table_S5_Pathway_Protein_Ranking.xlsx`（原 Table S4 重命名）：TumorTissues 面板由 178 扩展至 192 个蛋白质（新增 14 个食管癌相关 Kla-DDR 蛋白并赋序）。
  - `Supplementary_Table_S6_Lactylation_Regulators.xlsx`（原 Table S5 重命名）：49 个调控因子及其文献依据完整保留。

## 关键脚本

- 上游扩展：`R/candidate/prepare_escc_inclusion_inputs.R`
- 附表构建：`workflow/build_supplementary_workbooks.R`（支持环境变量 `KLA_PUBLICATION_INPUT` 与 `KLA_SUPPLEMENTARY_OUTPUT` 隔离输出）
- Figure 1：`R/candidate/build_figure1_category_boxplot.R`
- MKI67 全蛋白比值：`R/candidate/build_figure1_mki67_ratio_boxplots.R`
- Kla 七通路图：`R/candidate/build_ddr_pathway_summary_boxplots.R`
- 10 个调控因子样本散点图：`R/candidate/build_regulator_sample_percentile_scatters.R`
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
