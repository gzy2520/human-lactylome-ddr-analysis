# 当前项目审计：2026-09-19

审计基线：`b05872d2331f82edb3457e808435b270cec2a389`，研究分支 `research/kla31-rnaseq-20260912`。开始时 Git 工作区干净；本地比本地记录的远端分支领先 1 个提交，未刷新远端引用。此次仅新增审计脚本、证据和报告，没有修改原始数据、分析脚本、正式图表或服务器文件。

## 结论与范围

当前 RNA 交付不能认定为已验证的最终版。已确认数值换算错误仍在产物中、重复计数影响显著性、双因素模型未处理材料内配对，以及输出和文档版本混用。并非所有历史结果都失效：旧冻结包哈希完整，当前 RNA 矩阵的稳定 ID、样本唯一性和局部数据传递检查通过。

本轮覆盖当前 RNA Stage 2 → qsmooth → 教师问题图表主链，服务器两个 FPKM/RPKM 来源的原文件复算，最新交付目录，方法文档，以及蛋白组冻结完整性和新样本表观察单位。未重新核验全部公共数据的供体/建库信息，未从原始质谱重建所有蛋白组图，也未逐图完成全项目视觉验收；因此不能把局部通过表述为整个项目无问题。历史笔记只用于定位，结论依据当前文件和实际计算。

## 已确认的问题

### P1：FPKM/RPKM 错误仍存在于当前 RNA 数值

- 本地 `workflow/lib_kla31_expression_20260916.R:203` 已改为 FPKM 按列重标至一百万；服务器同名 helper 仍再次除以基因长度。
- 从服务器原始 BPH、TALL-104 文件重新提取，已有 RDS TPM 与旧错误公式逐元素最大误差均为 **0**。
- 正确公式与存量结果的绝对 log2(TPM+0.5) 差：BPH 中位数 **0.5911485**，最大 **6.501166**；TALL-104 中位数 **0.09945899**，最大 **5.219513**。这些是各自源矩阵中的逐元素差异，不是最终统计量变化。
- 两个源矩阵服务器与本地 MD5 一致；本地 28 个源矩阵到 qsmooth 保存的输入逐元素差异全为 0。因此错误确实进入当前 qsmooth 输入，不是仅服务器遗留备份有错。
- 正确 FPKM→TPM 为 `FPKM / colSums(FPKM) * 1e6`，无需再次除长度。本地修复仍保留长度表交集作为基因纳入过滤，审计复算保持相同基因集合。
- qsmooth_B 使用所有材料共同估计，修复两个输入后应重建整条归一化链，不能只替换两张终端图。

证据：`remote_check.txt`、`check.txt`、`statistics.txt`。原理：[COEX-Seq 原始论文](https://pmc.ncbi.nlm.nih.gov/articles/PMC6440655/)。

### P1：旧 31 行重复计数会改变显著性，且正式目录仍保留旧表

- `results/rna_teacher_questions/figure1a_fixed_ddr_31groups.csv`：31 行，28 个 ReferenceKey。
- 同目录样本表：**1,907 行 / 1,898 个唯一 SampleID**，多出的 9 行来自 HCT116 三次展开和 HK-2 两次展开。
- 交付目录已有 28 材料表及 1,898 行样本表，但还保留旧 31 行表和图。历史终端“只在 /tmp 试跑”的描述不代表当前交付状态。
- 对当前存量表达数值复算：DDR 四类 ANOVA，31 行 **F=3.925，p≈0.019**；28 行 **F=2.937，p≈0.0537**。重复计数改变了 0.05 阈值下的判定。
- DDR 占比对应 p 为约 `2.18e-7` 与 `4.09e-7`。以上均基于尚未修正 FPKM 的输入，不应作为最终结果。

证据：`check.txt`、`statistics.txt`；当前 renderer 的 `group_exp_use` 已去重，但不能据此认定所有产物已更新。

### P1：双因素 ANOVA 将同一材料的 DDR/Kla 分数当作独立观测

`workflow/render_teacher_questions_rna_exact.R:680` 使用 `MedianExpression ~ Category * GeneSet`，未建模同一 GroupID 的两个分数之间的相关性。28 材料表中两项分数 Pearson r=**0.9693074**。

对同一存量数据作模型敏感性检查：

| Category 主效应 | F | p |
|---|---:|---:|
| 当前普通双因素 ANOVA | 3.937 | 0.0137 |
| 增加 `Error(GroupID/GeneSet)` 配对结构 | 1.984 | 0.143 |

这表明误差结构有实质影响。配对模型在这里仅用于诊断，不是批准的新正式分析；还需结合研究问题、来源依赖及模型假设决定最终模型。不能以当前 p=0.0137 支撑可靠显著性结论。

证据：`statistics.txt`，复算脚本 `workflow/audit_statistics_20260919.R`。

### P1：现有重跑方式可能跳过修复或破坏元数据完整性

`workflow/build_31_expression_matrices_20260916.R:394`：只要 RDS 已存在就 `next`，没有源文件或代码哈希的缓存失效判据。仅同步新 helper 后原目录重跑会继续使用旧矩阵。

同脚本末尾仅将本次 `manifest_rows`/`qc_rows` 写入 `group_manifest.csv`/`group_sample_qc.csv`。若在原目录只重建 BPH、TALL-104，汇总元数据可能被覆盖成两组，而旧其余矩阵还在。旧 `group_failures.csv` 也不自动清理，handler 失败后可继续并最终正常退出。需要独立输出目录、完整元数据汇总及失败验收，不能直接照历史命令链无条件重跑。

### P2：qsmooth 的实际行为与文档结论不一致

- 当前 A 为 18,332×28，每 ReferenceKey 一列，组内没有重复。A 与未平滑的组中位数输入最大绝对差仅 **0.0008745567**；A 基本保留原谱。
- A 权重中位数为 0；B 权重中位数约 **0.2794**。不能将 A 的退化行为解释为“算法已证明差异是真实生物学而非批次”。研究来源与材料类型混杂仍存在，28 个 ReferenceKey 去重也不等于建立了完全独立的随机材料抽样。
- `docs/RNA_DATA_PROCESSING_AND_ANALYSIS_WORKFLOW.md:279` 声称无基因绝对差大于 0.5，但 `AB_comparison.csv` 中 **14,775 个基因的 MaxAbsDelta >0.5**，最大 **1.708813**。实际为所有基因的 **MeanAbsDelta** ≤0.5；文档把跨材料平均误差写成了逐元素上界。
- 当前 renderer 的 DDR/Kla 组统计由 B 样本分数再取材料中位数；调控因子热图使用 A。文档 §7 将部分组图归为 A 也不准确。

证据：`check.txt`、`statistics.txt`、`workflow/qsmooth_31group_20260916.R:114`。算法含义参照 [Bioconductor 官方手册](https://www.bioconductor.org/packages/release/bioc/manuals/qsmooth/man/qsmooth.pdf)。本轮不以此擅自更换归一化算法。

### P2：实际交付图仍有 NA、标题截断和多版并存

已渲染检查交付目录的 `Figure_1a_RNA_DDR_and_lactylated_gene_expression_28materials_boxplot.pdf`，其副标题仍为 **p = NA (F = NA)**，主标题左右截断。说明脚本中的 `trimws` 修复尚未体现在这张交付 PDF。截图见 `dual_28.png`。

交付目录同时存在旧 20 基因增殖评分、Hallmark 版本、31 行版、28 材料版以及旧 percentile 图。`QSMOOTH_README.md` 仍写 17,340 基因，当前实际为 18,332。`results/` 与 `outputs/` 均被 Git 忽略，Git 干净不能证明产物版本一致。需要唯一发布清单和文件/输入指纹。

### P2：方法文档和项目入口检查已过期

- 方法 §3.4.2 仍明确写错误的 FPKM 再除长度公式；§7 仍写旧 20 基因增殖面板与 31 点。
- `AI_HANDOFF_CURRENT.md` 仍说 results 为空、全部 RNA AnalysisReady=FALSE，已不能描述当前现场。
- `preflight_kla31.R` 当前失败。具体原因是旧本地输入清单中 **6 个历史 candidate 文件路径缺失、1 个 INPUT_MANIFEST 哈希不符**；并非冻结包损坏。
- 该 preflight 还要求旧 registry 的所有 AnalysisReady 为 FALSE，只验证旧合同，不能充当当前 RNA 全链验收。
- 文档把所有样本称为独立生物学样本，并将低相关性确定归因于生物学异质性，证据强度不足；唯一 SampleID、TPM 总量及相关系数不能单独证明生物学独立性或排除质量问题。
- GSE181540 handler 同时返回 counts 和 FPKM，运行分支优先使用 counts；manifest 的 ValueRoute 却优先判断 FPKM，因此这两个来源的路线说明也有误导性。

## 已通过的检查

1. 原冻结包 163 文件全部存在且 SHA-256 全部匹配；没有证据表明此前冻结文件被覆盖。
2. 当前 B：18,332 Ensembl 基因 × 1,898 唯一样本，无重复基因 ID、重复列名或非有限值，样本与 QC 元数据完整对应；A：18,332×28。
3. 28 个本地源矩阵在共有基因上的数值与 qsmooth 保存输入完全一致。此项证明传递一致，不证明源转换正确。
4. 当前 renderer 已实施 ReferenceKey 去重，并使用稳定 ID 定量；symbol 用于显示及部分原始符号源的一次性映射。基因集不能仅因标签写有 symbol 就判为违规。
5. 新蛋白组包三张样本输入实查为 272/264/504 行，全部 ObservationType=sample；其 7 项已有 validation 全为 TRUE。这不是对所有质谱源数据、统计模型及图形的全面背书。

## 建议的修复顺序

1. 保存本次基线及当前发布文件清单，明确唯一候选发布目录。
2. 同步正确的源转换函数，在独立输出目录重建受影响来源，完整恢复 28 材料元数据；补充转换数值验收及缓存失效规则。
3. 重建 qsmooth 及所有相关下游产物；固定 28 材料推断单位，保留蛋白组 31 行的对应映射。
4. 修复材料内配对统计结构，报告批次/来源混杂限制，重新评估结论。
5. 同步 Methods、入口说明、图名及版本清单；逐图验收后发布。历史产物移入明确归档，不覆盖以丢失审计线索。

## 复核方式

在研究工作树运行，脚本只读已有分析数据：

```sh
Rscript --vanilla workflow/audit_project_state_20260919.R
Rscript --vanilla workflow/audit_statistics_20260919.R
ssh -T -o BatchMode=yes -o ConnectTimeout=8 -o HostKeyAlias=192.168.3.45 user@100.121.229.123 'Rscript --vanilla -' < workflow/audit_remote_fpkm_20260919.R
```

远端复核脚本显式计算旧错误公式，避免服务器函数修复后对照定义发生变化。未来重建后的存量结果应与正确公式一致，而不再与旧公式一致。
