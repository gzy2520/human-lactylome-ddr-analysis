# RNA 分析简报（2026-09-19）

当前文件入口为 [`config/rna_current_release.csv`](../config/rna_current_release.csv)，方法见 [完整流程](RNA_DATA_PROCESSING_AND_ANALYSIS_WORKFLOW.md)。

- 31 条蛋白组材料记录对应 28 份不重复 RNA 参考矩阵、1,898 个唯一样本标识。HCT116/HK-2 共用参考在核心 RNA 组图中各只计一次。
- 以 Ensembl/Entrez 做 RNA 分析键，UniProt BaseAccession 做蛋白组键；Symbol 用于展示和原始符号源的一次性映射。
- Counts 先除基因长度再按列重标为 TPM；FPKM/RPKM 只按列重标，不再除长度。本次修复 BPH/TALL-104 并重建下游。
- qsmooth 仍使用原 A/B 两套算法：A 为材料中位数谱，B 为全部样本。A 每组一列时接近保留原谱，不能把权重 0 当作已排除批次的证据。
- 增殖图采用冻结的 196 基因 Hallmark G2M Checkpoint 面板；DDR 表达采用既有 371 基因固定面板；Kla 靶基因采用既有 5,224 基因集合。样本表达图及组比较来自 B，48 调控因子热图来自 A。
- 联合状态描述 RNA、全蛋白检出及 Kla 检出的对应关系；外部 RNA 不是与质谱样本配对的测量。未检出不等于蛋白或修饰不存在。
- 按用户要求保留既有显著性分析方法。本次只修复数值、运行安全、文档和交付一致性；统计假设未重新确证。

所有当前数值均在本次 release 的表格中；旧结果及旧比例表保留在历史目录，不用于当前汇报。验收包括完整范围、源矩阵变化范围、输入输出指纹和 PDF/PNG 检查。
