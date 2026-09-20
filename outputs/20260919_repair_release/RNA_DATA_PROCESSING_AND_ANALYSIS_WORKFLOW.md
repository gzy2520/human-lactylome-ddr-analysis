# RNA 数据处理与分析方法（2026-09-19 修复版）

当前发布路径以 [`config/rna_current_release.csv`](../config/rna_current_release.csv) 为准；逐文件 SHA-256、输入指纹及验收表随发布保存。旧文档已归档至 `docs/history/20260919_pre_repair_*`，其中历史数值不代表本次重建。

本次修复表达量换算、流程隔离、元数据、图形排版和版本管理；按用户要求，保留既有单因素及双因素 ANOVA 方法，不修订显著性分析。修正输入后所有数值按既有算法重新计算。

## 1. 范围与观察单位

蛋白组/Kla 有 31 条 GroupID；RNA 对应 28 个不重复 ReferenceKey、1,898 个唯一 SampleID。HCT116 的 KLA31_14/15/16 共用 RNA，HK-2 的 KLA31_27/28 共用 RNA。RNA 核心组图每 ReferenceKey 只使用一次（保留 KLA31_14 和 KLA31_27）；样本点图包含全部 1,898 个唯一样本。

31 行映射用于蛋白组联合状态和对应展示，不能解释成 31 份独立 RNA。1,898 个 SampleID 是矩阵样本标识数，不能据此声称所有样本都是独立供体或独立生物学重复。RNA 为外部同类材料基线，不能当作与质谱同一样本、同供体、同处理的配对测量。

范围与选样的机器可读入口为：

- `config/rnaseq_expression_extraction_contract_20260919.csv`：GroupID、ReferenceKey、SourceKind、SourcePath、SelectorValue、ExpectedN。
- `audit/20260916_full_31_rna_status/rna_31_group_status.csv`：来源选择和材料对应台账。
- 当前 expression 目录：`group_index.csv`、`group_manifest.csv`、`group_sample_qc.csv`。

## 2. 稳定 ID 与来源内提取

分析主键为 Ensembl gene ID；蛋白组主键为 UniProt BaseAccession，通过稳定 ID 映射衔接。Symbol 仅展示，或在原作者只提供符号的矩阵中用于一次性转为稳定 ID。

- Ensembl 版本后缀去除，保留 `_PAR_Y`。
- Entrez 经人类 NCBI gene2ensembl 映射；一对多歧义按既有规则排除，聚合规则记录在 IDRule。
- 胎盘原始 ENST 以 GRCh37 转录本注释优先映射，NCBI 表为后备。
- A549 与 T-47D 的符号源使用官方 Symbol→Entrez→Ensembl，并补入 HGNC `prev_symbol`；不将可歧义的 `alias_symbol` 全量合并。
- 原有源样本选择、缺失基因处理、基因聚合算法和组分类保持不变。HGNC、映射和基因长度文件指纹保存在提取输入记录中。

## 3. 表达量换算与 QC

### 3.1 Counts → TPM

对保留的有长度基因，以 Ensembl 111 外显子并集长度计算：

`rate[g,s] = counts[g,s] / (length_bp[g] / 1000)`

`TPM[g,s] = rate[g,s] / sum_g(rate[g,s]) * 1e6`

RSEM 等分数 counts 以原始精度计算 TPM；归档 counts 块沿用既有四舍五入规则。GSE181540 若有 counts 列，实际走 counts 路由；不能因为 handler 同时返回 FPKM 就把 provenance 标成 FPKM。

### 3.2 FPKM/RPKM → TPM

FPKM/RPKM 已经除过长度，正确换算是：

`TPM[g,s] = FPKM[g,s] / sum_g(FPKM[g,s]) * 1e6`

BPH 与 TALL-104 使用此路由。为保持既有纳入范围，仍取长度注释可用基因集合，但长度仅用于定义纳入集合，**不再除以长度**。旧结果重复做了长度归一化，本次从源文件重建；旧矩阵不能继续作为当前结果。

### 3.3 其他来源及限制

recount3 的 coverage 按既有方法换算；coverage 与 read count 仅在有效覆盖长度因子可视为一致时成比例，不能无条件宣称与所有 counts 定量完全等价。不同来源的建库、注释、基因空间和定量方式仍可能造成差异；统一为 TPM 不等于消除了来源批次。

### 3.4 QC

主表达值为 `log2(TPM+0.5)`。检查稳定 ID、样本列唯一性、有限非负 TPM、每列 TPM 总和、counts/log2 TPM 对应、样本数和合同完整性；报告 LibrarySize、GenesDetected、ZeroFraction 及组内 Spearman 范围。

相关性和零值比例仅为诊断，不能自动排除降解或样本误标。本流程不根据一个相关性阈值自动剔除样本。n=1 的 TALL-104 无组内相关性，报告 NA。现有质量指标也不能单独证明生物学重复独立性。

## 4. qsmooth：保留算法，准确解释诊断

28 个源矩阵的 Ensembl 交集为 18,332 基因。算法和随机种子保持原样（seed=25）：

- **A**：每材料先取样本中位数，18,332×28，以 ReferenceKey 分组。
- **B**：18,332×1,898 全样本，以 ReferenceKey 分组，保留组内样本变化。

权重对应**排序后的分位数位置**，因此边车列名为 `QuantileRank`，不再把它标成某个固定基因 ID。A 每组只有一列，组内变异无法从重复中估计，输出接近原始材料中位数；权重中位数为 0 不能作为“已辨认真实生物学差异”的证据。B 的权重需单独报告。qsmooth 本身不证明已排除来源批次混杂。

A/B 对照为 A 与 B 按材料汇总后的逐基因比较。`MeanAbsDelta` 是每个基因跨材料的平均绝对差，`MaxAbsDelta` 是最大绝对差，二者不能互换。当前 `summary.csv` 分别记录超过 0.5 的基因数和逐元素最大差；实际数值以本次产物为准，不再从旧报告手抄。

## 5. 蛋白组映射与联合状态

DDR 映射和蛋白组输入沿用既有稳定 ID 合同，不重新选择面板。`overlay_ddr_panel_20260916.py` 将当前 A 投影到 Kla∩DDR、reference DDR、GO DDR 三个面板，保留 31 行对应关系及共享基线标记。

`explore_rna_assisted_ddr_20260917.R` 重建材料×蛋白联合表。`plot_rna_reference_31group_20260917.R` 生成 `RNA_1` 与 `RNA_2*` 描述图。RNA 表达调用为平滑尺度 `log2(TPM+0.5) >= log2(1.5)`；没有 RNA 信号、全蛋白未检出或 Kla 未检出，均不能直接证明蛋白不存在或不存在乳酰化。

联合状态的数字由当前 reference 目录 CSV 提供；不沿用旧的 11,067 对及对应百分比。31 行描述性展开不是增加独立 RNA 样本数。

## 6. 核心图的实际计算口径

| 图/表 | 当前输入和计算 |
|---|---|
| Hallmark 增殖分数 | 冻结的 196 个 Ensembl G2M Checkpoint 基因；B 中每样本基因均值，再取材料内样本中位数；不是实测增殖速率 |
| 固定 DDR 表达 | 既有 Kla∩DDR 注释中的 371 个可用 Ensembl 基因；B 中每样本基因中位数，再取材料中位数 |
| Kla 靶基因表达 | 既有映射的 5,224 个基因；B 中每样本基因中位数，再取材料中位数 |
| DDR 表达占比 | B 中固定 DDR 面板达到阈值的基因数 / 18,332 共有基因空间中达到阈值的基因数 ×100%；材料值为样本占比中位数 |
| DDR/Kla 双面板 | 同一套 B 来源的 28 材料分数；保留当前 `Category * GeneSet` ANOVA |
| 调控因子热图 | A 中 48 个调控因子、49 条角色条目 ×28 材料；分别显示表达和跨材料行标准化 Z 值，无缺失填补 |

占比是受分子和分母共同影响的指标，不能直接表述为“DDR 基因表达更高”。组图的单因素 ANOVA 及双面板双因素 ANOVA 本次不变；该保留不等于本次审计证明了检验独立性等假设。

Hallmark 输入冻结为 `config/rna_hallmark_g2m_checkpoint_frozen_20260919.csv`，来自此前已采用的官方 Hallmark 交付映射，保存 Ensembl/Entrez 及展示标签。渲染时不再在线查询可能变化的基因集。

## 7. 重跑、验收与发布

1. Stage 2 只接受空输出目录；共享 ReferenceKey 只在同一次运行内复用。任何 handler 失败都会非零退出。部分提取不会写成完整 manifest/QC。
2. Stage 3 先确认合同中全部 28 个对象、对应矩阵及样本数齐备，再生成完整 manifest、QC 和 31 行 index。
3. `verify_31_expression_matrices_20260916.R` 检查全部 28 个对象；独立回归测试验证 FPKM 不受长度除法影响。
4. 本地运行 `rebuild_rna_release_20260919.sh <expression_dir> <new_label>`，重建 qsmooth、DDR overlay、联合表、reference 图、core 图及 validation。
5. 数值验收确认仅 BPH/TALL-104 源矩阵发生预期改变，其余 26 份不变；归一化的下游数值可能跨材料变化。输出图需另外渲染检查。
6. 渲染器不再自动复制到旧交付或桌面。通过 QA 后更新单一 release 指针、SHA-256 清单和桌面副本；旧版本保留为历史。
7. `Rscript --vanilla workflow/preflight_kla31.R` 检查当前发布与输入指纹，同时核验旧冻结包；旧历史 preflight 单独归档，不再要求当前 RNA 维持 AnalysisReady=FALSE。

## 8. 来源选择的既有特殊决定

以下为历史选样决定，本次不改变；进一步核实供体、建库及模型匹配属于额外的来源审核，不能由数值修复替代。
### 8.1 源数据替换与人工判定（必须披露）

| 组 | 问题 | 处理 |
|---|---|---|
| **KLA31_27/28 HK-2** | 原候选 GSE269418 的"溶剂对照"实际含**过夜血清饥饿 + DMSO + HCl + BSA**，非未处理细胞 | 换用 GSE240748 未处理臂，逐样本核验 SOFT 字段；且因作者表行名为 `gene-MIR6859-1` 式非稳定 ID，改用 NCBI 重定量表，**未把缺失补零** |
| **KLA31_30 NSC** | 蛋白组侧使用 ENSA 与 HMP1 两个模型；GSE119834 的 9 个文库中**仅 1 个与 ENSA 精确匹配**，HMP1 完全缺失 | 经人工确认 9 个作为"类别匹配外部模型"纳入；计划中的原始 FASTQ 重定量**最终未执行**，改用 NCBI 表；`GSM3384848` 因不在 NCBI 表中被剔除 → **n=8** |
| **KLA31_07 精子** | 计划 7 例 | NCBI 表中仅存 4 例 → **n=4**（报告中已声明不能报 n=7） |
| **KLA31_31 HUVEC** | 原候选 ENCODE `ENCSR000COZ` 不足以证明 n≥3 | 换用 GSE203551（4 例未处理） |
| **KLA31_01 肌腱 / KLA31_22 PC-3M** | 原计划 FASTQ 重定量 | 发现 NCBI 已有计数矩阵 → 改用 NCBI 表，FASTQ 仅留档 |
| **KLA31_14/15/16 HCT116、KLA31_24 RKO** | DepMap 每个模型仅 1 行汇总，无法估计重复间误差 | DepMap 降级为背景，改用作者计数矩阵 |
| **KLA31_08 BPH** | — | **全部 18 例供体均在服用 5α-还原酶抑制剂**，多数并用 α 受体阻滞剂。属背景用药而非实验臂，**披露而不校正** |
| **KLA31_06 胎盘** | GEO 无孕周/单胎字段 | 计划中的足月单胎筛选**未执行、不声称**；仅用源文件的 control 组（已排除 PE/IUGR） |
| **KLA31_17 TALL-104** | 全球唯一公开未扰动 bulk 文库 | **n=1**，无法估计组内方差、不能画误差棒 |

### 8.2 人工改标签与人工映射（必须披露）

1. **T-47D 列名与条件不一致**：作者矩阵把所选 4 列写作 `rep*_0h`（零时），而 GEO 记录为 `DMSO (24 hrs)`。按 GEO 判定为溶剂对照臂，4 个 C48 药物列排除。
2. **`GCN5 (KAT2A)`**：热图轴把 UniProt `Q92830` 人工显示为 `GCN5 (KAT2A)`（Q92830 官方名为 KAT2A，GCN5 为历史名），仅为展示。
3. **HCT116 / HK-2 的 (a)(b)(c) 后缀**：`plot_rna_reference_31group_20260917.R` 中 31 条手写标签给共享同一矩阵的行加上 (a)/(b)/(c)，**该后缀在数据中不对应任何生物学差异**。
4. **ESCC 组织学筛选**：TCGA-ESCA 中对自由文本 `primary_diagnosis` 做 `grepl("squamous")` 子串匹配以剔除腺癌，**仅此一个队列做了该筛选**。
5. **TCGA-ESCC 组织学筛选、GTEx 注释样本丢弃**：GTEx 肺 867 例注释中仅 578 例有表达行（海马 243 → 197），差额为**该 release 无表达行的注释子切片**，契约中"All GTEx v8 Lung samples"并非字面属实。
