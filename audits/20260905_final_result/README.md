# final_result 审计：需要修订后才能作为统一定稿

审计日期：2026-09-05（Asia/Shanghai）。检查对象为 `/Users/gzy2520/orca/workspaces/kla/figure_refine/final_result`，实际指向 `/Users/gzy2520/Desktop/Research/kla/results/final_figures_and_tables`。当前 fork 没有自己的 final_result，因此不能用 fork 的目录缺失来推断最终产物不存在。

已检查全部 40 个 PNG、40 个 PDF 的文件清单和来源；PNG 按 SHA-256 去重后为 32 个不同图像。所有 PDF 已提取文字，所有不同 PNG 已通过联系表检查，对主要问题另行打开原图；S3a PDF 另行栅格化确认与 PNG 的差异。审计前后核对 93 个文件，哈希未变。本次只新增审计文件，没有修改作图脚本、输入、final_result，也没有提交或推送。

结论：最终包不是统一的最新 31 组版本，也不满足“所有散点均为 sample”。存在旧图、新表混用、非 sample 混入、重复使用参考样本、重复测量统计模型错误及可见渲染问题。数据质量技能要求的观测粒度核查与验证技能要求的独立重算，使这些问题能够与单纯排版问题区分。

## 1. 每张图究竟使用什么版本

这里的“31 对”应准确理解为 31 个 Kla 数据集/样本组及其参考对应关系，而不是 31 对一一配对的生物样本。新增 PXD064038 的 6 个 Kla 样本与 PXD065830 的 94 个 ESCC 肿瘤全蛋白样本来自独立队列。

| 最终图/文件族 | 实际来源与覆盖 | 结论 |
|---|---|---|
| Figure 1a DDR fraction | 默认 30 组 candidate，210 观测（118 WP + 92 Kla）；肿瘤为 13 WP / 6 Kla | 旧数据；31 组应为 310 观测，肿瘤 107 / 12 |
| Figure 1b MKI67/H3C1 | 30 组基线中的可计算子集：35 观测、10 个组 | 旧数据；新增 ESCC 的 45 个有效比值全部缺失，31 组子集应为 80 观测、11 个组 |
| companion_controls/Figure_1_MKI67_over_ACTB_boxplot | 基线子集 44 观测 | 旧数据 |
| companion_controls/Figure_1_MKI67_over_TUBB_boxplot | 基线子集 51 观测 | 旧数据 |
| companion_controls/Figure_1_companion_MKI67_over_ACTB_boxplot | 扩展子集 49 观测，含 ESCC 5 例 | 31 组来源的有效子集；与上面的旧图并存 |
| companion_controls/Figure_1_companion_MKI67_over_TUBB_boxplot | 扩展子集 145 观测，含 ESCC 94 例 | 31 组来源的有效子集；与上面的旧图并存 |
| Figure 2a / 2b UpSet | 扩展正式输出；836 whole-proteome DDR / 401 Kla-DDR 蛋白 | 31 组来源核对通过 |
| Figure 2c BER、2d NER、2e MMR、2f FA、2g HR、2h AEJ、2i NHEJ | 全部匹配默认 30 组 candidate 的 PNG；每图 92 个观测 × 2 个方向；肿瘤 n=6 | 七张均为旧数据；最新应为 98 个观测、肿瘤 n=12 |
| Figure_2_DDR_pathway_summary_*_barplot 七个无字母别名 | 与 Figure 2c–i 分别完全相同 | 七个别名同样旧，不是另外一套新结果 |
| Figure 3a / 3b，含 no_frame、with_blue_frame 所有版本 | 扩展 heatmaps；全蛋白输入 1488 行，Kla 输入 1519 行，均覆盖 31 组 | ESCC 行存在，范围正确；缺失值语义另见第 3 节 |
| S1a DDR fraction by PXD | 扩展输入 310 观测、31 个组 | 31 组来源正确；混合观测粒度仍有问题 |
| S1b MKI67/H3C1 by PXD 与 11_detected_only | 扩展输入 80 个有效比值、11 个组；主布局列出 31 组，20 组标 ND | 新来源的有效子集；不是全部 31 组都有可用比值 |
| S2a / S2b UpSet | 24397 全蛋白 / 5814 Kla 蛋白 | 扩展来源核对通过 |
| S3a 肿瘤 / S3b 非肿瘤矩阵 | 扩展评分工作簿，192 / 183 个稳定蛋白 ID | 成员和状态计数核对通过；PNG 标题异常 |
| S4a 癌细胞 / S4b 正常细胞矩阵 | 扩展评分工作簿，381 / 292 个稳定蛋白 ID | 成员和状态计数核对通过；PNG 标题异常 |
| S3 legacy_tissue_summaries 两图、S4 companion_summaries 两图 | 匹配扩展正式输出，分母为各类别 Kla-DDR 蛋白并集 | 来源正确；是并集比例图，不是样本统计图 |

完整 40 行逐 PNG 清单：[`per_figure_audit.csv`](per_figure_audit.csv)。每个 PDF/PNG 与候选来源的精确匹配：[`figure_hash_provenance.csv`](figure_hash_provenance.csv)。CSV 中无来源哈希匹配的 S1 和独立矩阵采用脚本、输入、PDF 标签及图面交叉核查，不声称它们具有不存在的复制来源证明。

## 2. barplot / boxplot 的点是不是 sample

下表以输入中的 `ObservationType` 为证据。标记为 sample 是最低必要条件，不能仅凭这个字段证明生物学独立性；condition/model 中有些可能可进一步解析为独立观测，但目前不能直接当作通过独立样本审计。

| 图/输入 | 观测总数 | sample 标记 | 其他类型 |
|---|---:|---:|---|
| 当前 Figure 1a | 210 | 172 | 38：WP aggregate 12；Kla pool 2、dataset_union 4、condition 12、model 8 |
| 最新 S1a / 31 组 Figure 1 输入 | 310 | 272 | 相同的 38 个非 sample 观测仍存在 |
| 当前七通路图，每个方向 | 92 | 66 | 26：pool 2、dataset_union 4、condition 12、model 8 |
| 31 组七通路候选，每个方向 | 98 | 72 | 相同的 26 个非 sample 观测仍存在 |
| 当前 Figure 1b H3C1 | 35 | 34 | aggregate 1 |
| S1b H3C1，两种布局 | 80 | 79 | aggregate 1 |
| 旧 ACTB / 新 ACTB | 44 / 49 | 43 / 48 | 各有 aggregate 1 |
| 旧 TUBB / 新 TUBB | 51 / 145 | 43 / 137 | 各有 aggregate 8 |

七通路图中同一观测各贡献一个 Pro 点和一个 Inh 点，所以当前每张实际绘制 184 行点数据，但只有 92 个来源观测；不能把它说成 184 个独立样本。

具体例子：PXD075377 的 `HCC_pool` / `Control_pool` 是混样；PXD028488 的 `HCT116_dataset_union` / `HEK293T_dataset_union` 是组并集；MKI67 图的 PC-3M 观测来自 `PXD022005_aggregate`；TUBB 还包含 ProCan averaged protein matrix 的汇总谱。现有绘图脚本没有用 `ObservationType == "sample"` 过滤或分层展示。

额外发现真实重复：PXD058534 与 PXD078736 的两组 HK-2 都使用 PXD072220 的 `amostra1 / amostra3 / amostra4`。源文件、定量列和 H3C1 比值完全相同（约 0.4850662、0.5292882、0.9559819）。在 S1 按组展示时可注明参考复用关系；在 Figure 1 分类合并及其 ANOVA 中不能重复当作六个独立样本。只检查 `(PXD, SampleGroup, SampleID)` 唯一性抓不住这种跨组复用。

并集 summary 和 UpSet 的柱/点不应强行改称 sample：UpSet 的圆点表示集合归属，柱表示蛋白数；S3/S4 附带 summary 的柱表示类别蛋白并集比例，没有样本散点。它们本身属于不同分析单位，需要准确标注。

证据：[`observation_type_counts.csv`](observation_type_counts.csv)、[`non_sample_observations.csv`](non_sample_observations.csv)、[`baseline30_reused_reference_identifiers.csv`](baseline30_reused_reference_identifiers.csv)。参考 ID 重复表也包括可能同名但来自不同组织的观测，不能把整张表的每行都自动判为真实重复；上述 HK-2 是逐定量列确认的例子。

## 3. 统计和作图用法的问题

### 高：图与统计表对不上

Figure 1a 图上 F=142.81，是 210 行旧输入的结果；随包 `figure1_category_omnibus_anova.csv` 却为 N=310、F=217.39573。`figure1_category_boxplot_mean_median.csv` 的肿瘤均值为 WP 5.490037%、Kla 6.887262%，来自扩展数据，不是当前图的旧输入。

七通路图显示肿瘤 n=6，总来源观测 92；随包 `pathway_summary_two_way_anova.csv` 却是 NPoint=98、N=196、残差 df=188。独立重算当前图的模型是 N=184、残差 df=176。HR 的交互项当前图是 **（旧 q≈0.00282），扩展表为 *（q≈0.01138）。即使某些星号碰巧相同，也不能认作图表一致。

ACTB 两个版本甚至给出不同的显著性结论：旧图 q≈0.01347，新 companion 图 q≈0.06748。以上都是重现现有模型的数值，**不代表认可该模型有效性**。

### 高：七通路 ANOVA 漏掉同一样本内的配对

`R/candidate/boxplot_significance.R::compute_pathway_sample_two_way_anova()` 创建了 `SourceSampleID`，但拟合的是 `ValuePercent ~ CategoryFactor * DirectionFactor`，没有样本误差项或随机效应。Pro/Inh 来自同一来源观测，属于重复测量。现在把两列拉长后按独立误差拟合，p/q 值不能直接作为可靠推断使用。

后续需要按实际独立样本定义设计重复测量/混合模型，并处理研究来源/组的相关性。BH 校正只能调整多重检验，不能修复独立性错误。R 官方文档说明了多误差层与 `Error` 的作用：[aov 文档](https://www.stat.ethz.ch/R-manual/R-devel/library/stats/html/aov.html)。

### 高：Figure 1a 类别效应的表述不符合实际检验

代码写的是 `DdrFractionPercentage ~ CategoryFactor * DatasetFactor`，用默认顺序 ANOVA 取第一项类别 F，却将其描述成 `controlling for modality`。在这里明显不平衡的样本配置下，先进入模型的类别项并不等于“调整模态后的类别效应”。

只调换项的进入顺序、不改任何数值：30 组 F 从 142.8078 变为 103.7668；31 组从 217.3957 变为 151.1639。该诊断说明现有标签不成立，不把“调换顺序”本身作为最终统计修复。需要先明确要检验的类别效应及交互，并处理来源相关性。[R 顺序 ANOVA 文档](https://www.stat.ethz.ch/R-manual/R-devel/library/stats/html/anova.lm.html)。见 [`anova_term_order_sensitivity.csv`](anova_term_order_sensitivity.csv)。

### 高：跨研究合并不能忽略来源与真实独立单位

Figure 1 和七通路的推断把不同研究、组织/细胞系、实验条件、样本数及来源表处理方法直接合并，没有建模研究层级或参考复用。31 组方案里肿瘤 WP 的 107 个观测有 94 个来自一个 ESCC 队列；H3C1 非肿瘤组全部 5 个值来自 BPH 一个组。图可以描述当前观测分布，但四类之间的总体差异不能仅凭现有 pooled ANOVA 归因于生物类别。

SEM=`sd/sqrt(N)` 的算式本身无误，但 N 混合了 sample、pool、union、condition、model，因此当前 SEM 不能无条件解释为独立生物样本均值的不确定性。

### 中：热图把未检出记为 0，且与已检出后的低百分位混在一起

新增 ESCC 脚本在某样本未检出时填 `pcts[j] <- 0`，然后对全部样本取 median。它表示“零编码后的总体百分位中位数”，不是“检出样本中的丰度百分位中位数”。在整个 31 组输入中，全蛋白有 58 行、Kla 有 31 行满足 `Value == 0 && DetectedSampleCount > 0`（Kla 按角色行计数，不等于独立蛋白数）。例如 ESCC 全蛋白 TRIM33 检出 46/94 例，热图值仍为 0；Kla BRD4 检出 2/6，仍为 0。

因此白色不能直接解释为“无表达/完全未检出”。需要在图注明确零编码定义，或区分无定量、未检出与真实百分位。保留现有定义也可以，但不能把它包装成单纯表达量热图。证据：[`heatmap_detected_but_zero_percentile.csv`](heatmap_detected_but_zero_percentile.csv)。

### 中：S1b 的 ND 不能唯一表示“未检出”

脚本先把无有效 H3C1 比值的组补空，再统一标 ND。一个比值缺失可能来自 MKI67 缺失、H3C1 缺失、定量列不可用或数据仅为组汇总；也不能由“没有比值”推出两个蛋白都未检出。应区分 unavailable ratio、MKI67 not detected、H3C1 not detected 等原因。该图的 31 个组槽位不是 31 个有可用比值的数据集。

### 中：可见渲染/说明缺陷

- S3a、S3b、S4a、S4b 的 PNG 标题均出现空数字，例如 `tumor tissues (    proteins)`；对应 PDF 文本中数字存在。S3a PDF 重新渲染可见完整 `192 proteins`。原因尚未定位，不能直接归因于人为涂改或特定字体。证据图 [`S3a_pdf_render.png`](S3a_pdf_render.png)。
- `S1b_11_detected_only` 左上方 Non-tumor tissues 被裁成不完整文字。
- S4 正常细胞 companion summary 的 HR 标注 `124 (42.5%)` 右侧被裁切。
- 当前七通路主图移除了说明文字，图面未解释误差棒是 SEM；Pro/Inh 与具体分母应写入正式图例。
- 四张矩阵及并集 summary 需要图例说明彩色=促进、黑色=抑制、浅灰=未归到相应有符号状态。矩阵中的 0 不应被解读为实验测得“无作用”，且蛋白的多通路标注不是互斥构成比例。

## 4. 本轮核对通过的部分与边界

- 31 组注册表确有 31 个唯一 `(PXD, SampleGroup)`，Table S1/S2 的组汇总页也均为 31 行。
- 四套 UpSet 的 BaseAccession 无重复，15 个互斥交集重新求和与并集一致；四个矩阵成员与同类别 Kla-DDR 集合一致，状态只包含 -1/0/+1。
- 肿瘤 192 蛋白的评分摘要重算为 HR 促进74/抑制8、NHEJ 46/8、BER 27/0、NER 23/2、AEJ 16/6、MMR 13/0、FA 11/0，与当前并集 summary 一致。
- 当前及扩展 sample 输入的 DDR 比例、通路比例重新按对应计数/分母计算，未发现算术差错；不能把该结论延伸成全部源样本均独立或上游鉴定全部正确。
- 直接读取 PXD065830 原始工作簿，按 P46013/P60709/P07437/P68431 稳定 accession 核查：94 个 ESCC-T 列；ACTB/TUBB/H3C1 分别有 5/94/45 个有效比值，存储值与原表重算最大误差小于 6e-15。旧 handoff 中“H3C1 全部缺失”的文字已过时。
- 集合/比值/通路的核心分析采用稳定 UniProt accession；Gene Symbol 用作展示。未发现本轮检查路径中用 symbol 替代核心集合键的错误。随机抖动种子为25。
- 本轮是最终产物与当前冻结/扩展源表的审计，不是从全部原始质谱重新鉴定、不是对每个研究所有生物样本来源的完整人工复核，也不验证所有通路人工注释的文献真实性。未运行会覆盖结果的 publication pipeline。

## 5. 建议修订顺序

1. 固定唯一 31 组输入入口及逐图 manifest，阻止默认 30 组路径进入新定稿；由同一次执行导出图片和统计表。
2. 建立来源层面的独立样本键，审查 reference 复用；sample、pool、model、condition、union 分类处理。不要将汇总行复制成伪重复来满足“全31组都有sample”。
3. 保留全31组的覆盖/并集视图；样本推断图使用真正可解析的样本，并透明列出无法提供独立样本的组及原因。
4. 按研究/组层级及同一样本内 Pro/Inh 配对重设推断模型，再算p/q和误差棒。合理性诊断应针对拟合模型，不能仅要求生成有限p值。
5. 重新生成 Figure 1 与七通路图、同步表格；整理新旧别名；修复标题、裁切及缺失值图例，并对 PNG/PDF 成对检查。

复查命令（只写本审计目录）：

```sh
Rscript --vanilla audits/20260905_final_result/audit.R
```

脚本：[`audit.R`](audit.R)。计算环境：[`sessionInfo.txt`](sessionInfo.txt)。保护检查：[`preservation_check.csv`](preservation_check.csv)。数据和代码指纹：[`inputs_and_scripts_sha256.csv`](inputs_and_scripts_sha256.csv)。
