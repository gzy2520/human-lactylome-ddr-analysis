# 最早DDR柱状图算法与TALL历史

## 已核实的最早版本

主仓库不是浅克隆（`git rev-parse --is-shallow-repository`返回false）。所有现存分支只有一个根提交：`e5a89ac59bf37c5424b3f23cfaa414f1fb409bf6`，提交时间为2026-08-07。根提交保存了更早的Python分析与其R迁移/扩展版本。`reanalysis/reports/REFERENCE_PROTEOME_DDR_ANALYSIS.md`文内注明分析日期2026-08-05，且给出明确TALL数值。

这份最早可验证的10模型柱状图记录为：

|系列|材料|DDR数|分母|比例|
|---|---|---:|---:|---:|
|Kla|PXD028488 TALL-104，旧标签T-ALL|37|368|10.0543478261%|
|Whole proteome|PXD030304 TALL-1主替代背景|267|3383|7.8924031924%|

对应数值直接见快照`e5a89ac_REFERENCE_PROTEOME_DDR_ANALYSIS.md`第48行。最早图名为`cell_type_kla_vs_reference_ddr_fraction`，另有`_v2`别名；生成脚本为根提交的`reanalysis/scripts/analyze_reference_proteome_ddr.py`（快照497–508行绘制两系列barh，538–562行以分子/分母标注）。

## Kla柱的实际算法

`e5a89ac_run_pipeline.py`第146–158行：按CellOrTissueType选择primary证据，建立`set(subset["BaseAccession"])`，与GO accession集合求交，以交集大小除以集合大小。每个细胞材料只有一个去重蛋白集合。重复PSM、不同位点、不同采集和所纳入两臂中的相同蛋白只计一次。

TALL的入选范围由`e5a89ac_extractors.py`第256行决定：同时遍历`Enrichment-Search files`和`Nonenrichment-Search files`。第270–272行只有HEK全HCD聚合目录被标为validation_only；第350行设置`PrimaryIncluded=not validation_only`。因此TALL的两个目录均进入primary集合：

- `data/PXD028488/search_results/Enrichment-Search files/TALL-NALAC-Search files/DB search psm.csv`
- `data/PXD028488/search_results/Nonenrichment-Search files/TALL-Nonenrichment-Search files/DB search psm.csv`

同目录的`protein-peptides.csv`提供肽起点和accession匹配；`proteins.csv`提供注释；`PSM ions.csv`提供156.1025诊断离子记录。PSM处理排除反向/污染accession，要求`-10lgP >= 20`，识别Kla修饰位置，并核对肽到蛋白位置映射。诊断离子用于证据注释，不作为所有PSM的强制纳入门槛（extractors第292–350行）。

因此37/368来自纳入上述两臂后、按BaseAccession去重的组集合。它不是三个富集样本比例的均值，也不包含样本间SD/SEM。后续只取富集臂三个raw file再平均，改变了输入范围和聚合方法，数值应分别解释。

## 分母与DDR注释

Kla分母为该细胞材料全部通过筛选的去重BaseAccession；isoform后缀去除。GO文件为`data/annotations/GO-repair+damage(human).tsv`。`e5a89ac_common.py`第209–214行将GENE PRODUCT ID规范为BaseAccession，排除NOT，保留taxon 9606或缺失taxon，按accession去重后做交集。此Kla统计不通过基因符号替代匹配。

最早Python全蛋白组函数以每模型提取的去重蛋白/特征表行数作为分母（`e5a89ac_analyze_reference_proteome_ddr.py`第365–386行）。先按BaseAccession匹配GO，未命中时可用已有GeneSymbol辅助（第318–360行）。TALL-1读取ProCan averaged矩阵的`SIDM00370;TALL-1`非空蛋白特征行；6个DIA运行未分别计算误差线。

同根提交的扩展R脚本`e5a89ac_analyze_expanded_ddr_fraction_by_accession.R`已明确accession-only：第623–647行将TALL-104映射到primary表的T-ALL并去重；第1047–1075行求DDR交集，直接计算`length(kla_ddr)/length(kla)`及reference对应比例；第1546行用geom_col绘图。对于ENSP背景，这一历史扩展算法保留源蛋白ID作分母，显式映射后判定源ID是否命中DDR；这一点与当前按映射UniProt去重的恢复值需要分开记录。

## 后续TALL背景替换

`372651b7fa454b6f4c450c96d38dc497db37c476`（2026-08-08）的原始diff已保存为`372651b_reference_change.diff`：删除TALL-104→T-ALL旧背景别名，新增TALL-104调用`extract_peaks_reference()`读取同研究背景文件。该函数选取任何Area列大于0的蛋白并去重，仍然输出组检测集合；不是样本比例平均。

当前冻结`group_summary_31.csv`的TALL Kla 37/368延续上述历史Kla集合。后续whole背景换成同研究TALL non-enrichment，与最早TALL-1的267/3383是不同来源；不能混用。

## 更早归档记录

根提交的run_pipeline.py第491行引用`archive/reanalysis_v1_2026-07-21/outputs/01_kla_extraction/kla_proteins_by_dataset_cell_type_group.csv`，第560–563行还记录OldWorkflowDirectory只含TALL-NALAC等富集目录，证明当时存在更早流程。当前主目录不存在该archive目录，Git根提交未跟踪该归档数据；因此无法仅凭该路径还原其柱状图数值。当前能由脚本、输入范围和正式数值表三者共同证实的最早柱状图就是上述37/368版本。

所有调查采用git show/read-only读取。当前源文件、图和已冻结结果均未更改。快照原始路径、提交号和SHA256记录在`source_manifest.csv`。


## 隐藏及忽略文件补查

对主项目执行`rg --files --hidden --no-ignore`，排除`.git`、`.worktrees`与`node_modules`，再查找archive、reanalysis_v1、2026-07-21、20260721及kla_proteins_by_dataset_cell_type_group关键字，均未找到对应本地文件。另以精确glob查询两种目标名称，同样无结果。完整枚举文件数与命令保存在`ignored_archive_search.txt`。据此，当前主项目本地和Git可交叉验证的最早数值仍为分析日期2026-08-05的37/368。

## 1.79%的真正首次引入提交

提交`4239155006e4c37740f276a0bba7362c552d898b`时间为2026-09-01 21:57:31 +08:00，标题为Add sample-level boxplot candidate for Figure 1。此时脚本尚无`read_tall_by_sample`；TALL仍在single_groups内，以冻结`TALL-104_dataset_union`输入。该提交实际跟踪的`sample_boxplot_values.csv`记录368、37、10.054347826087%。因此1.79%不是在4239155引入。

真正引入发生在`6735dc8`（2026-09-02 14:29:48 +08:00，Add four-category Figure 1 boxplot candidate）：

1. `6735dc8_prepare_sample_boxplot_inputs.R`第402–416行新增`read_tall_by_sample`，仅按`Area Sample 1/2/3 > 0`从proteins.csv提取所有蛋白；没有Kla修饰肽或位点筛选。
2. 第717–720行对富集目录proteins.csv调用该函数，将全蛋白检测集放入Kla sample_records。
3. 第759–766行将这些蛋白left join至冻结Kla membership的IsDdr，并将未命中赋FALSE。
4. 第768–778行用全蛋白检测数作为分母，用冻结Kla集合中命中的DDR蛋白数作为分子。

该提交实际跟踪的`figure1_sample_boxplot_values.csv`已保存：

|样本|错误分母|冻结Kla-DDR匹配分子|比例|
|---|---:|---:|---:|
|Sample 1|1662|29|1.744885679904%|
|Sample 2|1686|30|1.779359430605%|
|Sample 3|1678|31|1.847437425507%|

三个值的算术平均为1.790560845338%。这条链将富集蛋白表全部检测蛋白作为Kla分母，再只给冻结Kla成员保留DDR标签，导致分母和分子采用不同纳入规则。冻结join与NA赋FALSE在4239155已存在，但TALL当时输入是正确的368蛋白冻结集合；到6735dc8改变TALL输入后，才产生上述错误。

两提交的脚本与实际数据表均已保存快照并登记SHA256，时间定位由Git内容差异和提交数据表共同确认。
