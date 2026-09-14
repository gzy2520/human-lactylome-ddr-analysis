# 31组 RNA 参照：按解剖匹配与可整合性重新筛选

日期：2026-09-14。此文件是数据选路候选，不是表达分析结果。完整机器可读路由在 [rnaseq_source_routing_candidate_31.csv](../outputs/20260914_reference_material_audit/rnaseq_source_routing_candidate_31.csv)，候选注册表在 [rnaseq_reference_candidate_31.csv](../outputs/20260914_reference_material_audit/rnaseq_reference_candidate_31.csv)。两者全部保持 `AnalysisReady=FALSE`。

## 结论

应当按这一要求重新筛选，但“更好”不等于寻找一个覆盖所有组织和模型的单一数据库。先固定解剖组织、取材类别和细胞模型；其次才选择在该范围内单位、稳定 ID 和样本信息最可审计的数据。将 GTEx、GDC、DepMap、ENCODE、GEO 的 counts/TPM/FPKM/RPKM 直接拼成一个矩阵，会把来源与建库差异当成生物学差异，不能作为主分析输入。

GSE181540 的 GEO 简表只写 normal skin；原始论文的方法说明三例增生性瘢痕与其**邻近全层正常皮肤**一起采集。因此 KLA31_03 和 KLA31_04 满足瘢痕/邻近皮肤的严格解剖学要求，证据应固定为论文方法而非只引用 GEO 简表。[论文](https://doi.org/10.3389/fcell.2021.748703) · [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE181540)

邻近肝仍固定为 TCGA-LIHC `Solid Tissue Normal`，而不是 GTEx liver；HCC 固定为同项目 `Primary Tumor`。这种取材类别优先于文件格式统一。

## 数据路由

|路由|组别|输入|稳定分析主键|处理原则|
|---|---|---|---|---|
|严格原始重定量|KLA31_01、KLA31_22|GSE236746、GSE235595 选定 FASTQ|Ensembl gene ID|前者仅发布 ENST 转录本 FPKM；后者没有选定 PC-3M DMSO 对照的矩阵。两者用同一声明的基因层流程重定量。|
|统一的公开重处理计数|KLA31_06、07、29、30|recount3 G026：SRP148556、014899、117021、161553|Ensembl gene ID|逐一核对所选 SRR 均在对应项目的同一基因计数矩阵内。NSC 的具体模型不匹配仍保留。|
|来源内基因矩阵|KLA31_03、04、08、17、23、25、26、27、28|GEO 已发布的 counts/RPKM/RSEM/STAR matrix|Ensembl 或 Entrez ID|保留其真实单位；在来源内 QC 和汇总，不和其他来源原始量级合并。|
|固定发布版量化|KLA31_02、05、09–16、18–21、24、31|GTEx v8、GDC STAR counts、DepMap Public 24Q4、ENCODE|Ensembl 或 DepMap 的 Entrez ID|每种来源固定版本、样本筛选和文件哈希。|

KLA31_13–16、18–21、24 是 DepMap 外部亲本基线；同一 HCT116 参照也不应重复充当三个独立转录组。它们不能被解释为对应 KO、感染或共培养条件。

## 这次新增核实

- **HEK293T**：GSE203529 的 `RAW.tar` 已确认三个 WT untreated 文件：GSM6175912、GSM6175918、GSM6175922。它们是 RSEM 基因层结果，含 versioned Ensembl gene ID、`expected_count`、TPM 和 FPKM；不需要再下载 FASTQ。
- **BPH**：GSE132714 的 `RAW.tar` 已确认含 GSM3890639–GSM3890656 的 18 个 BPH 每样本 RPKM 文件，行名为 versioned Ensembl ID。使用此来源内矩阵，并在临床注释审核后锁最终纳入样本；不为了获得整数 counts 而改用不相同的前列腺组织。
- **瘢痕/邻近皮肤**：GSE181540 输入臂的同一 XLSX 具有 Ensembl gene ID、三例 HS 和三例 NS 的 count 与 FPKM 列；仅使用 input，排除 IP/甲基化峰数据。
- **MES28、HK-2、HMC3**：分别已有可定位的 `read_count`、报告 raw-count 和 STAR readcount 矩阵；保留原始稳定 ID，核对具体对照列后使用。

服务器上的 FASTQ 队列已从 137 个文件收缩为 12 个真正的基因层缺口文件：肩袖肌腱 6 个文件和 PC-3M DMSO 6 个文件。已完成的文件与 `.part` 均保留，队列使用 MD5 和期望字节数验收。通过 `./workflow/check_server_rnaseq_download_20260914.sh` 查看当前范围和进度。

## 检索后不采用的替代项

- **GSE180836**：31 例肩关节肌腱撕裂的样本量较大，但混合创伤性与退行性撕裂。它可作为病理肩袖肌腱的敏感性队列，不能无条件替换当前明确的非糖尿病撕裂组。
- **GSE199484**：材料是撕裂冈上肌和未损伤肩胛下肌，属于肌肉，不可替代肩袖**腱**。
- **GSE189134**：明确是瘢痕与邻近组织配对，但检测为 miRNA qPCR array，不是普通全转录组。
- **DepMap**：24Q4 `Model.csv` 有 MCF10A 和 TALL-1，却没有可用于替代 MCF10A 或 TALL-104 的相应表达行；TALL-1 也不是 TALL-104，因此不替换现有精确模型来源。

## 跨来源分析契约

1. 不建立未经论证的“31组共同 counts/TPM/FPKM 矩阵”。
2. 在每一来源内完成样本 QC、重复汇总和单位记录；跨来源只使用预先声明的可比摘要，例如每样本的基因秩、稳定 ID 映射后的基因集分数，或各来源内得到的效应量。
3. Ensembl/Entrez 是分析键，Symbol 只用于显示。先锁注释版本，再剥离 Ensembl version；DepMap 重复 Entrez 列按可追溯规则处理。
4. 所有来源继续记录样本、材料、条件、文库、单位、文件 SHA256 和选择理由。随机步骤固定 `set.seed(25)`。

recount3 原始文件与 R/Bioconductor 访问方式见其[官方快速访问说明](https://rna.recount.bio/docs/quick-access.html)；其 G026 基因计数文件使用 GENCODE v26，应作为固定版本记录，而不是与其他注释版本默认为同一版本。
