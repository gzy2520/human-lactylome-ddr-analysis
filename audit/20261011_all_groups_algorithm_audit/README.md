# 31组原始算法与当前算法的全量核查

## 结论

本次从项目原表审查了 **31个材料×两种检测类型＝62个位置**，覆盖前版325条来源记录：组织Kla 44条、细胞Kla 54条、全蛋白组227条。逐组比较原始组级算法、旧样本算法与当前计算，检查原始蛋白ID字段、分子分母集合及样本/技术采集拆分。

发现并修正两类问题：**4个细胞Kla材料的7个样本误计反向序列；4个位置把作者完整处理矩阵当作逐样本检出集合。** 另将HK2未核定的条件解释恢复为源标签ctr/man。当前版本为 `deliverables/20261011_lactylation_ddr_all_groups_audited/`。

总表 `all_62_algorithm_comparison.csv` 为每个位置列出历史组比例、审计前后均值、ID解析、分子分母、观察单位和具体证据。Historical指冻结31组流程的组并集值；TALL最早Git版本的专项追溯另见`audit/20261011_tall_algorithm_history/`，没有将不同年代的数据范围当作同一统计量。

## 蛋白ID解析错误：四材料、七样本

|材料|受影响样本|分母修正|DDR分子|
|---|---|---|---|
|HepG2|SIRT3KO_rep1/rep2|1831→1830；1596→1595|不变|
|HMC3|H0/H24|1277→1276；1087→1086|不变|
|pretreated HK2|HK2|2164→2163|不变|
|HK2 ctr/man|ctr_1/ctr_2|962→961；952→951|不变|

原表代表Protein中存在`REV__O14732`或`REV__Q14093`，但Reverse与污染标志为空。历史解析有的先去掉REV__前缀，有的读入已去前缀的Proteins字段，均把反向序列当作正常UniProt。修正时先按原始代表字段剔除反向行，再重建样本集合；7个样本都没有对应accession的其他合格forward支持，因此分母各减1。全部DDR分子不变。

|材料|前版均值%|纠正后均值%|
|---|---:|---:|
|HepG2|11.696698|11.698228|
|HMC3|11.341986|11.351657|
|pretreated HK2|10.166359|10.171059|
|HK2 ctr/man|10.596034|10.599823|

新数值变化均小于0.01个百分点。TALL的完整多accession解析重新核对通过，仍为11.421069%。本次还扫描原表全部Protein、Proteins、Majority及Leading字段；不能只查实际用于计数的那一列。全蛋白组中的混合候选字段含REV/CON时，正常forward成员仍保留，未把“附属候选含reverse”等同于整条正常代表均为reverse。

证据：`kla_cells/confirmed_all_reverse_rows.csv`、`affected_id_forward_support.csv`、`corrected_affected_sample_counts.csv`及独立重建membership。组织12位置和全蛋白31位置均未发现这类新增计数错误。

## 样本拆分与误差线：四个完整处理矩阵

|位置|原展示来源列数|处理矩阵的DDR/蛋白|当前处理|
|---|---:|---:|---|
|Hippocampus Kla|3|29/853|1个作者处理panel，无样本SEM|
|Hippocampus whole|3|219/6082|同上|
|Adjacent liver whole|8|194/4251|同上|
|HCC whole|8|271/5188|同上|

海马Kla S3是筛选后的quantifiable Class I位点表。三列均全正值，相同最小值分别出现73、70、1279次；论文另记原始3650个位点仅1167个位点在三样本均鉴定。Whole矩阵三列相同最小值分别出现84、69、594次。当前量化表没有保留完整原始检出/缺失模式，不能用>0恢复原样本检出集合。

邻癌肝与HCC来源ANTs/CISs是相对比值矩阵，8列全正，Non-rec1整列为1。其重复相同蛋白集合得到的SEM=0也不表达原始样本检出率的分散程度。四位置现均按一个作者处理矩阵集合展示，图上标`p`，保留原始来源列ID及列数，不再把重复集合用作样本误差估计。没有推定作者具体采用哪种缺失填补方法。

这四处柱高不变。分析记录由325变为307，62个柱位保持完整，可画SEM的位置由40变为36。Tendon的源表确有0/正值检测分界，三个RCT通道去重集合相同是原表事实，保留其来源单位。

## 有依据的方法变化，单独记录

- **Tendon Kla**：原全部鉴定并集571/30；当前reporter正值并集464/23。被删107蛋白、7DDR在全部六通道都无正reporter，差异不是排除正常组织通道造成。`___1/2/3`是位点multiplicity列，不能算三个技术样本；当前OR合并不增加样本数。
- **Whole的Tendon、HepG2、sperm**：全部鉴定集合→所选列正值检出集合，分别从3984/154→3574/140、6772/450→6713/448、10680/461→9492/385。每个减少项均有原表无正值证据。
- **Lung、placenta whole**：原分母是ENSP等SourceProteinID，当前分母是映射后去重UniProt。保留一对多和未映射记录，每个UniProt直接判DDR，避免将某源ID的any-DDR标签传播给全部映射对象。
- **MCF10A whole**：旧6字符正则漏81个合法10字符UniProt，当前已保留；使用Leading Razor Protein，避免把共享肽全部候选当作额外检测蛋白。
- **ProCan背景**：旧averaged表达feature矩阵→实际run肽检出集合；每个技术制备两次仪器采集合并为一个单位，三个制备的SEM标`t`。HCT两位置共用参考，复用已记录。
- **HMC3 whole**：当前是LFQ可定量集合；与Intensity>0集合不同。H0/H24和HCT116的NC116/R116各为两个不同条件，每条件一列，不计算跨条件重复SEM。
- **HEK、HCT亚细胞分级和MCF10A**：HCD能量、fraction、DDA/DIA及重复采集按原材料合并，不增加生物样本数。NSC/GBM不同模型保留模型均值，不画模型内样本SEM。

这些选择改变了检测集合或统计对象，逐位置均在总表和来源报告中注明。当前图中的样本均值不要求与历史组并集比例逐一相等。

## 仍需原始来源补足的样本注释

PXD078736的本地metadata仅确认ctrl与man各三个重复，man的生物学意义注明待发表。历史mannitol与hypoxia-reoxygenation标签冲突，显示标签已改`HK-2 (ctr/man)`，当前Condition采用ctr/man。SampleGroup字段保留冻结注册键用于历史关联，其旧字符串不作为当前条件释义。

Tendon的1–3→RCT映射沿用既有样本映射，当前源表仅给通道号；本次未获得新增临床通道对应证据。海马若需原始每样本检出率，需要原始site/PSM sample-identification表或作者S2，目前本地处理表不足以恢复该mask。

## 复现、核验和文件

- `kla_tissues/`：12个组织Kla位置的原表重建、全部ID字段扫描、tendon正值筛选分解、海马量化表证据。
- `kla_cells/`：19个细胞Kla位置的独立解析、PEAKS/RKO等原解析器集合对照、4组反向序列修正。
- `whole_proteome/`：31个全蛋白位置、227来源记录的重建与ID/原始列扫描、技术制备登记。
- `cross_group/`：组并集与均值对照、来源复用、条件、完整panel核查；最终数据独立28项检查全部通过。
- `all_62_algorithm_comparison.csv`：62位置总表。

最终图从上版冻结membership和本次原表独立修正集合生成：

```sh
Rscript workflow/lactylation_ddr/apply_all_group_audit.R
Rscript workflow/lactylation_ddr/render.R
```

各子目录保留原表重建脚本，Python使用含pandas/openpyxl的环境，R使用data.table/readxl。所有原始文件保持不变；历史算法与输出保留用于比较。新图仅4个Kla柱高微调、4个panel取消样本SEM，以及HK2条件显示修正；RNA输入和表达值全部保持不变。

图中增加31个RNA样本数标注和62个DDR来源单位数标注。全蛋白组与乳酸化组分别统计；处理矩阵标为1 panel，技术制备、模型和QC单位保留类型标记。
