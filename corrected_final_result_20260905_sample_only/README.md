# Kla DDR Project — Final Publication Figures & Supplementary Tables

本目录汇总了本项目（31 组数据集，纳入食管癌 PXD064038 乳酸化及配对 PXD065830 肿瘤全蛋白对照）最新精细化修订版本的所有正式主图（Figure 1~3）、附图（Figure S1~S4）及附表（Tables S1~S6、Table 14 与统计检验表）。

---

## 目录结构速览

```text
results/final_figures_and_tables/
├── Figure_1/
│   ├── Figure_1a_DDR_fraction_boxplot.pdf / .png             # Fig 1a: 4大类全蛋白与乳酸化 DDR 比例箱线图（精简副标题，去除长说明文字；mean/median线细化；四大类在X轴底部相邻）
│   ├── Figure_1b_MKI67_over_H3C1_boxplot.pdf / .png         # Fig 1b: 全蛋白 MKI67/H3C1 增殖比率箱线图（四大类采用统一色彩体系无缺失，mean/median线细化，画幅与1a一致）
│   └── companion_controls/                                  # 增殖对照内参（ACTB, TUBB）
│       ├── Figure_1_MKI67_over_ACTB_boxplot.pdf / .png
│       └── Figure_1_MKI67_over_TUBB_boxplot.pdf / .png
├── Figure_2/
│   ├── Figure_2a_whole_proteome_DDR_UpSet.pdf / .png        # Fig 2a: 全蛋白 DDR 4集合交集 UpSet 图 (836 蛋白)
│   ├── Figure_2b_Kla_DDR_UpSet.pdf / .png                   # Fig 2b: 乳酸化 Kla-DDR 4集合交集 UpSet 图 (401 蛋白)
│   ├── Figure_2c_DDR_pathway_summary_BER_barplot.pdf / .png # Fig 2c: BER 通路立式 SEM 柱状图 + 样本散点 (1:1 正方形画幅，精简副标题与图例，去长文本)
│   ├── Figure_2d_DDR_pathway_summary_NER_barplot.pdf / .png # Fig 2d: NER 通路立式 SEM 柱状图 + 样本散点 (1:1 画幅)
│   ├── Figure_2e_DDR_pathway_summary_MMR_barplot.pdf / .png # Fig 2e: MMR 通路立式 SEM 柱状图 + 样本散点 (1:1 画幅)
│   ├── Figure_2f_DDR_pathway_summary_FA_barplot.pdf / .png  # Fig 2f: FA 通路立式 SEM 柱状图 + 样本散点 (1:1 画幅)
│   ├── Figure_2g_DDR_pathway_summary_HR_barplot.pdf / .png  # Fig 2g: HR 通路立式 SEM 柱状图 + 样本散点 (1:1 画幅)
│   ├── Figure_2h_DDR_pathway_summary_AEJ_barplot.pdf / .png # Fig 2h: AEJ 通路立式 SEM 柱状图 + 样本散点 (1:1 画幅)
│   └── Figure_2i_DDR_pathway_summary_NHEJ_barplot.pdf / .png# Fig 2i: NHEJ 通路立式 SEM 柱状图 + 样本散点 (1:1 画幅)
├── Figure_3/
│   ├── Figure_3a_reference_regulator_percentiles.pdf / .png            # Fig 3a: 全蛋白 48 个调控因子百分位数热图 (肿瘤组织扩充为3组，右侧标签竖排左旋90°，图例bar移至右侧，红框高亮)
│   ├── Figure_3a_reference_regulator_percentiles_no_frame.pdf / .png   # Fig 3a: 无高亮框纯净版热图
│   ├── Figure_3b_Kla_regulator_percentiles.pdf / .png                  # Fig 3b: 乳酸化调控因子百分位数热图 (肿瘤组织扩充为3组，右侧标签竖排左旋90°，图例bar移至右侧，蓝框高亮版)
│   ├── Figure_3b_Kla_regulator_percentiles_with_blue_frame.pdf / .png  # Fig 3b: 蓝框高亮版明确别名
│   └── Figure_3b_Kla_regulator_percentiles_no_frame.pdf / .png         # Fig 3b: 无框纯净版热图
├── Supplementary_Figure_S1/
│   ├── Figure_S1a_DDR_fraction_by_PXD.pdf / .png            # Fig S1a: 31个具体组织/细胞系数据集全蛋白与乳酸化DDR比例图（横轴采用直观组织名，顶部四大类彩色条标示）
│   ├── Figure_S1b_MKI67_over_H3C1_by_PXD.pdf / .png         # Fig S1b: 31个数据集全蛋白MKI67/H3C1比值图（未检出标记ND，顶部四大类条，与S1a严格对齐1:1）
│   └── Figure_S1b_MKI67_over_H3C1_11_detected_only.pdf / .png # Fig S1b备选: 仅含11个定量检出数据集的紧凑版
├── Supplementary_Figure_S2/
│   ├── Figure_S2a_whole_proteome_UpSet.pdf / .png           # Fig S2a: 全蛋白组 4集合全量蛋白 UpSet 图 (原S1a自动加一，24,397 蛋白)
│   └── Figure_S2b_Kla_proteome_UpSet.pdf / .png             # Fig S2b: 乳酸化组 4集合全量蛋白 UpSet 图 (原S1b自动加一，5,814 蛋白)
├── Supplementary_Figure_S3/
│   ├── Figure_S3a_DDR_pathway_matrix_tumor_tissues.pdf / .png     # Fig S3a: 肿瘤组织 DDR 通路矩阵独立图 (原2c拆分一，192 蛋白)
│   ├── Figure_S3b_DDR_pathway_matrix_non_tumor_tissues.pdf / .png # Fig S3b: 非肿瘤组织 DDR 通路矩阵独立图 (原2c拆分二，183 蛋白)
│   └── legacy_tissue_summaries/                                   # 组织横向总结备选图
│       ├── Figure_2d_DDR_pathway_summary_tumor_tissue.pdf / .png
│       └── Figure_2e_DDR_pathway_summary_non_tumor_tissue.pdf / .png
├── Supplementary_Figure_S4/
│   ├── Figure_S4a_DDR_pathway_matrix_cancer_cell_lines.pdf / .png # Fig S4a: 癌细胞系 DDR 通路矩阵独立图 (原S2a拆分一，381 蛋白)
│   ├── Figure_S4b_DDR_pathway_matrix_normal_cell_lines.pdf / .png # Fig S4b: 正常细胞系 DDR 通路矩阵独立图 (原S2a拆分二，292 蛋白)
│   └── companion_summaries/                                       # 细胞系横向总结备选图
│       ├── Supplementary_Figure_S2b_DDR_pathway_summary_cancer_cell_lines.pdf / .png
│       └── Supplementary_Figure_S2c_DDR_pathway_summary_normal_cell_lines.pdf / .png
└── Tables/
    ├── Table_S1_Kla_Data.xlsx                               # Table S1: 31组乳酸化蛋白数据及位点详情
    ├── Table_S2_Reference_Data.xlsx                         # Table S2: 31组全蛋白参考数据
    ├── Table_S3_Human_DDR_GO_Annotations.xlsx               # Table S3: 人类 6,707 条 DDR GO 注释基准
    ├── Table_S4_Venn_Membership.xlsx                        # Table S4: 4集合分析 15 区域归属清单 (含401 Kla-DDR)
    ├── Table_S5_Pathway_Protein_Ranking.xlsx                # Table S5: DDR 通路蛋白打分排序表 (含食管癌 192 蛋白)
    ├── Table_S6_Lactylation_Regulators.xlsx                 # Table S6: 49 个乳酸化调控因子及文献依据
    ├── Table_14_new_ESCC_Kla_DDR_proteins.xlsx / .csv       # 14 个新增食管癌 Kla-DDR 蛋白核对分类明细
    └── statistical_tests/                                   # 统计检验表 (ANOVA 等)
        ├── figure1_category_omnibus_anova.csv
        ├── figure1_category_one_way_anova.csv
        ├── figure1_category_boxplot_mean_median.csv
        └── pathway_summary_two_way_anova.csv
```

---

## 本次修订核心要点回顾

1. **画幅文字精简（去冗余副标题与长说明）**：
   - 因画布收窄与 1:1 比例调整，移除了图底可能发生截断的多行文字说明（caption），保留简洁明了的统计显著性标注（ANOVA p/q 值），文章正文图例（Figure Legend）可完整承载背景信息。
2. **Mean 和 Median 线条细化**：
   - 全面收细各箱线图中的 Median（中位数，linewidth 由 1.35~1.40 pt 降至 0.65~0.75 pt）与红色 Mean（均值，linewidth 由 1.35~1.55 pt 降至 0.70~0.80 pt）线宽，视觉更加精细利落。
3. **四大类色彩系统完备统一**：
   - 四大生物分类（Non-tumor tissues: `#0072B2`, Tumor tissues: `#D55E00`, Normal cell lines: `#009E73`, Cancer cell lines: `#CC79A7`）全面赋色，彻底消除了“有的有颜色有的为灰色/无色”的不一致现象。
4. **Figure 3 热图排版、图例与数据重大修正**：
   - **癌组织扩充为 3 组**：成功补齐食管癌（ESCC，PXD064038 乳酸化 6 样本 / PXD065830 全蛋白 94 样本）的调控因子百分位数，使 `tumor tissues` 面板完整呈现 3 组（ESCC、Prostate cancer、Primary HCC），全图为完整的 31 组。
   - **右侧分类标签竖向排列（左旋 90 度）**：由下至上竖向排布（`angle = 90`，头部左倾 90° 可自然横向阅读）。
   - **图例移动至右侧**：渐变色条（colorbar）由底部移至画布右侧垂直放置。
   - **蓝色版本双版输出**：提供一版带蓝色高亮框（`#08519C`，10 个关键调控因子列）与一版完全无框纯净版（`_no_frame`）。
