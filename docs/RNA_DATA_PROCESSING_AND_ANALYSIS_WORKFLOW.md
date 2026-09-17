# 31 组 Kla ∩ DDR 转录组（RNA-seq）数据处理与分析全流程规范

**文档版本**：V1.0 (2026-09-17)  
**分析范围**：31 类生物学材料（9 非肿瘤组织、3 肿瘤组织、12 癌细胞系、7 正常细胞/模型）  
**主要脚本**：
- 矩阵构建与归一化：[`workflow/qsmooth_31group_20260916.R`](../workflow/qsmooth_31group_20260916.R)
- DDR 映射与覆盖：[`workflow/map_ddr_uniprot_to_ensembl_20260916.py`](../workflow/map_ddr_uniprot_to_ensembl_20260916.py)、[`workflow/overlay_ddr_panel_20260916.py`](../workflow/overlay_ddr_panel_20260916.py)
- 转录组辅助生物学解析：[`workflow/explore_rna_assisted_ddr_20260917.R`](../workflow/explore_rna_assisted_ddr_20260917.R)
- 出版级可视化：[`workflow/plot_rna_reference_31group_20260917.R`](../workflow/plot_rna_reference_31group_20260917.R)

---

## 1. 背景与核心设计原则

### 1.1 为什么需要在 Kla 蛋白质组中引入转录组？
蛋白质质谱分析（LC-MS/MS）基于随机抽样（Data-Dependent Acquisition / DDA 或 DIA）以及检测动态范围的物理限制，且修饰蛋白质组（Lactylome）本质上仅提供二元或半定量的“检出/未检出”状态。  
当某一 DNA 损伤修复（DDR）蛋白质在特定生物材料中**未被检出乳酰化**时，存在两种机制假设：
1. **转录/翻译沉默**：该材料中该基因未转录或未翻译成蛋白质；
2. **生物学修饰特异性（真阴性）**：该基因在该材料中正常转录且蛋白质高丰度存在，但赖氨酸残基特异性未发生乳酰化。

质谱本身在结构上无法区分上述两种情况。引入全转录组（RNA-seq）作为**连续定量表达基准**，是解析蛋白质乳酰化选择性、剥离转录丰度混杂的决定性证据链。

### 1.2 核心分析约束与生信规范
- **严禁使用 Gene Symbol 进行分析**：Gene Symbol 存在别名混淆、版本退役和命名漂移风险；全流程以不可变的 **Ensembl Gene ID**（`ENSG...`，GRCh38 / Ensembl v111）和 UniProt **BaseAccession** 为唯一主键，Symbol 仅用于最终展示。
- **固定随机种子**：涉及所有随机化、抖动采样及数值计算，严格固定随机种子为 `25`（项目统一标识）。
- **材料类别参考谱属性（Reference Profile）**：31 组 RNA-seq 数据来自公开发表的权威高质量队列（GTEx、TCGA、DepMap、GEO、ENCODE），与项目中的 Kla 质谱样品为**材料类别级匹配**，而非同受试者严格配对。因此，所有跨组学比较均在类别尺度（Material Class）进行，绝不伪造配对检验。
- **版本控制与目录隔离**：脚本受 Git 严密追踪；代码与输出目录完全隔离，避免历史错误数据污染。

---

## 2. 31 组材料来源与多级审计架构

31 组生物材料涵盖人类主要组织与细胞模型，划分为四大生物学类别：

```
31 组材料结构
├── 非肿瘤组织 (Non-tumor tissues, n=9)
│   ├── KLA31_01: 肩袖病理性撕裂肌腱 (GSE236746)
│   ├── KLA31_02: 正常人肺组织 (GTEx v8)
│   ├── KLA31_03: 增生性瘢痕组织 (GSE178411)
│   ├── KLA31_04: 瘢痕邻近未受损皮肤 (GSE178411)
│   ├── KLA31_05: 人海马脑区 (GTEx v8)
│   ├── KLA31_06: 正常足月妊娠胎盘 (GSE114691)
│   ├── KLA31_07: 人精子 (GSE40181)
│   ├── KLA31_08: 良性前列腺增生 BPH (GSE132714)
│   └── KLA31_09: 癌旁正常肝组织 (TCGA-LIHC Solid Tissue Normal)
├── 肿瘤组织 (Tumor tissues, n=3)
│   ├── KLA31_10: 食管鳞状细胞癌 ESCC (TCGA-ESCA 确诊鳞癌，排除腺癌)
│   ├── KLA31_11: 前列腺癌 (TCGA-PRAD 原发肿瘤)
│   └── KLA31_12: 肝细胞癌 HCC (TCGA-LIHC 原发肿瘤)
├── 癌细胞系 (Cancer cell lines, n=12)
│   ├── KLA31_13: MCF7 (DepMap ACH-000019)
│   ├── KLA31_14/15/16: HCT116 共享基线 (DepMap ACH-000971)
│   ├── KLA31_17: TALL-104 (GSE163787，Plan A 严格 bulk 库)
│   ├── KLA31_18: HepG2 亲本基线 (DepMap ACH-000739)
│   ├── KLA31_19: A549 (DepMap ACH-000681)
│   ├── KLA31_20: MDA-MB-468 (DepMap ACH-000849)
│   ├── KLA31_21: T-47D (DepMap ACH-000147)
│   ├── KLA31_22: PC-3M 细胞对照 (GSE235595)
│   ├── KLA31_23: 胶质母细胞瘤干细胞 MES28 载体对照 (GSE266884)
│   └── KLA31_24: RKO 亲本基线 (DepMap ACH-000943)
└── 正常细胞系与培养模型 (Normal cell lines/models, n=7)
    ├── KLA31_25: HEK293T 未处理野生型 (GSE203529)
    ├── KLA31_26: HMC3 溶剂对照 (GSE275256)
    ├── KLA31_27/28: HK-2 共享未处理对照 (GSE240748)
    ├── KLA31_29: MCF10A 对照 (GSE103520)
    ├── KLA31_30: 神经干细胞 NSC 对照模型 (GSE119834)
    └── KLA31_31: HUVEC 全细胞基线 (ENCODE ENCSR000COZ / GSE203551)
```

> [!IMPORTANT]
> **共享参照显式追踪**：  
> 三个 HCT116 组（KLA31_14/15/16）以及两个 HK-2 组（KLA31_27/28）分别共享同一份转录组基线，因此实际包含 **28 个独立参考矩阵**。在下游归一化及膨胀至 31 组时，通过 [`outputs/20260916_qsmooth_31group/group_expansion_31.csv`](../outputs/20260916_qsmooth_31group/group_expansion_31.csv) 显式记录 `ReferenceKey`，杜绝重复计算伪自由度。

---

## 3. 表达量统一化与平滑分位数归一化（qsmooth）

### 3.1 表达尺度统一
所有数据源经服务器端清洗、去重复、ID 映射后，统一换算为 $\log_2(\text{TPM} + 0.5)$ 尺度。在 28 个独立矩阵中，取全集交集，共锁定 **17,340 个共同 Ensembl 基因**（交集瓶颈为 T-47D 的 18,815 个基因）。

### 3.2 组织感知归一化（qsmooth / YARN）
- **传统分位数归一化（Global QN）的局限**：强制假设所有样本的经验分位数分布一致。然而，脑组织（海马）、生殖组织（精子）、上皮肿瘤与培养细胞系之间存在巨大的内源性转录组差异，Global QN 会人为抹杀真实生物学差异。
- **平滑分位数归一化算法**：使用 Bioconductor 包 `qsmooth`（v1.22.0），根据组内变异与组间变异的比率，为每个基因的分位数计算连续权重 $w_g \in [0, 1]$：
  - $w_g \to 1$：偏向全局分位数归一化（技术噪音为主）；
  - $w_g \to 0$：保留组内特有分位数分布（生物学真实差异）。
- **归一化诊断结果**：
  - 全基因权重中位数严格为 **0**（`summary.csv`），证实大多数基因保留了其材料类别的固有分位数，符合组织特异性转录组特征。
  - **方案 A 与方案 B 交叉验证**：
    - 方案 A（每组材料取中位数向量，17,340 × 28）
    - 方案 B（全样本矩阵输入，17,340 × 1,898）
    两者基因间表达相关系数中位数达到 **0.9898**，绝对差值中位数仅 0.21 $\log_2$ 单位，无任何基因差值大于 0.5。因此，方案 A 坍缩参考谱作为跨材料基准高度可靠。

---

## 4. DDR 修复通路与 UniProt-Ensembl 双向映射

### 4.1 严格无 Symbol 映射策略
1. **DDR 靶蛋白集合**：401 个 Kla ∩ DDR 靶蛋白。
2. **双重官方映射**：
   - 途径一：UniProt 官方 `idmapping.dat`（UniProtKB-AC $\to$ Ensembl Gene ID）；
   - 途径二：NCBI `gene2ensembl` 通过 Entrez GeneID 间接映射。
   - 交叉比对：在两种途径皆覆盖的 882 个 accession 中，879 个完全一致（99.66%）；3 个分歧记录在案并显式解析。
3. **覆盖度**：
   - 401 个 Kla ∩ DDR 蛋白质中，**377 个**成功覆盖到 17,340 表达谱空间（对应 **357 个非冗余 Ensembl Gene ID**）。未覆盖的 24 个蛋白主要由于基因在某一稀有组织矩阵中缺失。

### 4.2 DDR 修复通路体系（8 大分面）
基于 S4 DDR 注释库，将 357 个基因划分至 8 大修复通路：
- 同源重组（**HR**, 106 基因）
- 碱基切除修复（**BER**, 54 基因）
- 核苷酸切除修复（**NER**, 28 基因）
- 范可尼贫血通路（**FA**, 22 基因）
- 非同源末端连接（**NHEJ**, 17 基因）
- 错配修复（**MMR**, 12 基因）
- 替代末端连接（**AEJ**, 3 基因）
- 未归类 DDR 相关（**unassigned**, 115 基因）

> [!NOTE]
> **关键逻辑修正**：S4 注释表中的方向状态包括 `+1`（上调促活）、`-1`（下调抑制）和 `0`（未判定方向性）。早期代码错误将 `0` 判定为非成员，导致基因错误坍缩至单一通路。现已严格修正为包含全部有效成员。

---

## 5. 转录组辅助的核心生物学发现与阴性控制

### 5.1 乳酰化未检出状态的四分类解析（11,067 对）
对 31 个材料类别 × 357 个 DDR 基因构成的 **11,067 个 (Material $\times$ Protein) 观测对**进行全面状态分解：

| 状态分类 | 组对数 (Pairs) | 占比 (%) | 生物学与技术解释 |
|---|---|---|---|
| **Captured as Kla** | 3,085 | **27.9%** | 转录本活跃、全蛋白存在且质谱明确捕获到乳酰化修饰肽段 |
| **Present, expressed, not Kla** | 2,319 | **21.0%** | **转录本活跃（>组内中位数）、全蛋白稳定存在，但质谱特异性无乳酰化** |
| **Expressed, not detected** | 1,299 | **11.7%** | 转录本活跃，但受限于质谱检测深度，全蛋白与修饰谱均未捕获 |
| **Low transcript** | 4,364 | **39.4%** | 转录本处于低丰度或沉默状态（$\le$ 组内中位数），修饰缺失源于转录本低下 |

> [!TIP]
> **21.0% 的生物学决定性意义**：  
> 蛋白质组学仅能回答“有没有捕获到”，无法证实未检出蛋白是否“真正缺乏修饰”。**这 21.0%（2,319 对）是转录组给予蛋白质组最关键的赋能**——它证明了大量 DDR 蛋白在细胞中高丰度表达且翻译成蛋白，但细胞对其乳酰化进行了高度特异性的修饰抑制，直接确立了 Kla 对 DDR 通路的生物学调控选择性。

### 5.2 表达水平梯度验证
三类检测状态在转录组中的表达水平呈现阶梯式分布：
- **Captured as Kla**：$\log_2(\text{TPM} + 0.5)$ 中位数为 **4.94**；
- **In reference proteome only**：中位数为 **4.56**；
- **Neither（均未检出）**：中位数为 **3.53**。  
这证实全蛋白检出与乳酰化检出均依赖于基础转录丰度，但 Kla 检出蛋白对高表达具有更强偏好。

### 5.3 关键阴性对照：过度乳酰化率受质谱深度主导
在探索性分析中，计算了各材料类别的“过度乳酰化率（excess rate）”。然而混杂因素分析表明：
- 组级 excess rate 与蛋白质组的**质谱检测深度（Proteome Depth）呈现高达 80% 的决定系数（$R^2 \approx 0.80$）**。
- 一旦通过多元回归移除深度影响，四大生物学类别之间的 excess 差异完全消失。
- **决定**：严格记录该阴性结果（见 [`outputs/20260917_rna_assisted_ddr/validity_vs_depth.csv`](../outputs/20260917_rna_assisted_ddr/validity_vs_depth.csv)），不在论文正文中将 excess 作为生物学机制主张，坚决采用经受住深度校正的稳健比值指标。

---

## 6. 出版级图表成果汇总

全部图表均采用顶刊（Nature/Cell 风格）标准重构，支持系统矢量字体（Arial Unicode MS），并输出 300 DPI PNG 与 cairo-PDF：

| 图表编号 | 图像文件 | 尺寸 (in) | 核心内容与视觉重构亮点 |
|---|---|---|---|
| **RNA_1** | [`RNA_1_DDR_panel_expression_31groups.png`](../results/rna_reference_31group/RNA_1_DDR_panel_expression_31groups.png) | 11.0 × 13.5 | **357 DDR 基因 × 31 组材料表达热图**。<br>• 顶部集成四大生物类别彩色 Banner（n=9, n=3, n=12, n=7）；<br>• 左侧 8 大 DDR 通路采用独立科研色标分面；<br>• 31 列标注生物学材料名称（如 Tendon, Lung, ESCC, HCC）；<br>• 9 阶高对比度 RdBu 渐变。 |
| **RNA_2a** | [`RNA_2a_non_detection_meaning.png`](../results/rna_reference_31group/RNA_2a_non_detection_meaning.png) | 6.8 × 5.0 | **乳酰化未检出状态分解柱状图**。<br>• 柱顶标明精确对数，柱内嵌白色百分比；<br>• 突出标注 21.0%“表达翻译但未乳酰化”的核心生物学特异性。 |
| **RNA_2b** | [`RNA_2b_expression_by_detection_state.png`](../results/rna_reference_31group/RNA_2b_expression_by_detection_state.png) | 6.8 × 5.0 | **检出状态对应的表达密度曲线**。<br>• 柔和半透明填充配合清晰轮廓线；<br>• 虚线标示三状态真实中位数（4.94 vs 4.56 vs 3.53）；<br>• 右上角卡片化图例。 |
| **RNA_2** | [`RNA_2_lactylome_transcriptome_coupling.png`](../results/rna_reference_31group/RNA_2_lactylome_transcriptome_coupling.png) | 12.8 × 5.2 | **主刊级 2a + 2b 组合复合图**。<br>• Nature 风格加粗 A/B 标签，直接用于论文正文或扩展数据展示。 |

---

## 7. 目录文件结构对照

```
.
├── audit/
│   ├── 20260916_full_31_rna_status/
│   │   └── rna_31_group_status.csv               # 31 组 RNA 样本元数据与审计底册
│   └── 20260917_rna_assisted_ddr/
│       ├── absence_decomposition.csv             # 四状态频数与比例
│       ├── expression_threshold_sweep.csv        # TPM 绝对阈值扫描
│       ├── within_group_percentile_sweep.csv     # 组内分位数扫描
│       └── validity_vs_depth.csv                 # 质谱深度混杂阴性对照
├── outputs/
│   ├── 20260916_qsmooth_31group/
│   │   ├── matrices/
│   │   │   └── qsmooth_A_collapsed_log2tpm.tsv.gz # 归一化后 17,340 × 28 参考表达矩阵
│   │   └── group_expansion_31.csv                # 28 矩阵向 31 组的显式扩展映射
│   ├── 20260916_ddr_panel_31group/
│   │   ├── kla_ddr_annotation.csv                # DDR 八通路与 GO 注释表
│   │   └── kla_ddr_expression_31groups.tsv.gz    # 377 个 DDR 蛋白对应表达子矩阵
│   └── 20260917_rna_assisted_ddr/
│       └── group_by_protein_pairs.tsv.gz         # 11,067 组对详细判定底表
├── results/
│   └── rna_reference_31group/                    # 出版级最终图表与数据包
│       ├── RNA_1_DDR_panel_expression_31groups.png / .pdf
│       ├── RNA_2a_non_detection_meaning.png / .pdf
│       ├── RNA_2b_expression_by_detection_state.png / .pdf
│       ├── RNA_2_lactylome_transcriptome_coupling.png / .pdf
│       ├── RNA_1_ddr_panel_genes.csv
│       ├── RNA_2_non_detection_counts.csv
│       ├── README.md
│       └── sessionInfo.txt
└── workflow/
    ├── qsmooth_31group_20260916.R                # qsmooth 归一化主脚本
    ├── explore_rna_assisted_ddr_20260917.R       # 转录组辅助生物学解析脚本
    └── plot_rna_reference_31group_20260917.R     # 出版级图表绘制脚本
```
