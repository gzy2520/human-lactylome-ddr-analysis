# 31 组 Kla ∩ DDR 转录组（RNA-seq）数据处理与分析全流程规范

**文档版本**：V1.2 (2026-09-18)  
**分析范围**：31 类生物学材料（9 非肿瘤组织、3 肿瘤组织、12 癌细胞系、7 正常细胞/模型），涵盖 28 个独立参考矩阵，共计 1,898 个高质量生物学测序样本；共同基因空间 **18,332 个 Ensembl 基因**  
**本次修订（V1.2）**：补入 §7 逐图方法学索引与 §8 特殊操作清单；修正因 HGNC 改名回退而失效的全部旧数值（17,340 → 18,332、357 → 371 DDR 基因、11,067 → 11,501 配对），详见 §3.3.2 与 §8.4。
**主要脚本体系**：
- **Stage 0 注释标尺与映射表构建**：
  - 外显子并集长度标尺：[`workflow/server_prepare_gene_annotation_20260916.sh`](../workflow/server_prepare_gene_annotation_20260916.sh)
  - 跨组装转录本映射：[`workflow/server_prepare_grch37_transcript_map_20260916.sh`](../workflow/server_prepare_grch37_transcript_map_20260916.sh)
- **Stage 1 原始数据获取与路由**：
  - DepMap 24Q4 基线获取与校验：[`workflow/server_download_depmap_20260914.R`](../workflow/server_download_depmap_20260914.R)、[`workflow/validate_depmap_expression_20260914.py`](../workflow/validate_depmap_expression_20260914.py)
  - TCGA GDC STAR Counts 获取：[`workflow/build_gdc_star_counts_manifest_20260914.R`](../workflow/build_gdc_star_counts_manifest_20260914.R)、[`workflow/server_download_gdc_star_counts_20260914.R`](../workflow/server_download_gdc_star_counts_20260914.R)
- **Stage 2 & 3 多源矩阵标准化提取与汇总**：
  - 核心提取与标准化库：[`workflow/lib_kla31_expression_20260916.R`](../workflow/lib_kla31_expression_20260916.R)
  - 31 组多源矩阵流式提取：[`workflow/build_31_expression_matrices_20260916.R`](../workflow/build_31_expression_matrices_20260916.R)
  - 审计底册与元数据汇总：[`workflow/finalize_31_expression_matrices_20260916.R`](../workflow/finalize_31_expression_matrices_20260916.R)
- **Stage 4 跨组织组织感知归一化**：
  - 平滑分位数归一化：[`workflow/qsmooth_31group_20260916.R`](../workflow/qsmooth_31group_20260916.R)
- **Stage 5 DDR 映射与覆盖**：
  - 双向映射与通路覆盖：[`workflow/map_ddr_uniprot_to_ensembl_20260916.py`](../workflow/map_ddr_uniprot_to_ensembl_20260916.py)、[`workflow/overlay_ddr_panel_20260916.py`](../workflow/overlay_ddr_panel_20260916.py)
- **Stage 6 & 7 辅助生物学解析与出版级可视化**：
  - 转录组辅助生物学解析：[`workflow/explore_rna_assisted_ddr_20260917.R`](../workflow/explore_rna_assisted_ddr_20260917.R)
  - 出版级可视化：[`workflow/plot_rna_reference_31group_20260917.R`](../workflow/plot_rna_reference_31group_20260917.R)

---

## 1. 背景与核心设计原则

### 1.1 为什么需要在 Kla 蛋白质组中引入转录组？
蛋白质质谱分析（LC-MS/MS）基于随机抽样及检测动态范围的物理限制，修饰蛋白质组（Lactylome）仅提供"检出/未检出"状态。

本项目的蛋白质组数据**已有配对的未富集全蛋白参照组（Ref）**。因此对每个（材料 × 蛋白）组合可以得到 2×2 的检出状态，其中"**蛋白在全蛋白组检出、但不在乳酰化组（Kla−/Ref+）**"这类组合，**仅凭蛋白质组自身即可判定**，不需要转录组。

转录组能补充的只有一种情形：**两套质谱均未检出（Kla−/Ref−）**。此时无法从蛋白质组判断该蛋白是"确实不存在"还是"存在但两套方法都没捕获"。引入全转录组（RNA-seq）作为连续表达基准，用于区分这两种来源。

> 需要说明的是：本项目蛋白质组数据的核心指标是**检出集合的比值**（`KlaDdrFraction` 及与参照组的差值），并非蛋白丰度。转录组在本流程中定位于**解释检出缺失的来源**，不与蛋白丰度做定量比较。

### 1.2 核心分析约束与生信规范
- **严禁使用 Gene Symbol 进行分析**：Gene Symbol 存在别名混淆、版本退役和命名漂移风险；全流程以不可变的 **Ensembl Gene ID**（`ENSG...`，GRCh38 / Ensembl v111）和 UniProt **BaseAccession** 为唯一主键，Symbol 仅用于最终展示。
- **固定随机种子**：涉及所有随机化、抖动采样及数值计算，严格固定随机种子为 `25`（项目统一标识）。
- **材料类别参考谱属性（Reference Profile）**：31 组 RNA-seq 数据来自公开发表的权威高质量队列（GEO 18 组、GTEx v8 2 组、TCGA/GDC 4 组、recount3 1 组），与项目中的 Kla 质谱样品为**材料类别级匹配**，而非同受试者严格配对。因此，所有跨组学比较均在类别尺度（Material Class）进行，绝不伪造配对检验。
- **版本控制与目录隔离**：脚本受 Git 严密追踪；代码与输出目录完全隔离，避免历史错误数据污染。

---

## 2. 31 组材料来源与 28 个参考矩阵全景架构

31 组生物材料涵盖人类主要组织与细胞模型，划分为四大生物学类别。其中，三个 HCT116 组（KLA31_14/15/16）以及两个 HK-2 组（KLA31_27/28）分别共享同一份高质量转录组基线，因此全套材料由 **28 个独立参考矩阵**（包含 1,898 个独立样本）严格映射生成。

```
31 组材料结构
├── 非肿瘤组织 (Non-tumor tissues, n=9)
│   ├── KLA31_01: 肩袖病理性撕裂肌腱 (GSE236746)
│   ├── KLA31_02: 正常人肺组织 (GTEx v8)
│   ├── KLA31_03: 增生性瘢痕组织 (GSE181540)
│   ├── KLA31_04: 瘢痕邻近未受损皮肤 (GSE181540)
│   ├── KLA31_05: 人海马脑区 (GTEx v8)
│   ├── KLA31_06: 正常妊娠胎盘 (GSE114691；GEO 无孕周字段，未做足月筛选)
│   ├── KLA31_07: 人精子 (GSE65683)
│   ├── KLA31_08: 良性前列腺增生 BPH (GSE132714)
│   └── KLA31_09: 癌旁正常肝组织 (TCGA-LIHC Solid Tissue Normal)
├── 肿瘤组织 (Tumor tissues, n=3)
│   ├── KLA31_10: 食管鳞状细胞癌 ESCC (TCGA-ESCA 确诊鳞癌，排除腺癌)
│   ├── KLA31_11: 前列腺癌 (TCGA-PRAD 原发肿瘤)
│   └── KLA31_12: 肝细胞癌 HCC (TCGA-LIHC 原发肿瘤)
├── 癌细胞系 (Cancer cell lines, n=12)
│   ├── KLA31_13: MCF7 (GSE157383 / DepMap ACH-000019)
│   ├── KLA31_14/15/16: HCT116 共享基线 (GSE253699 / DepMap ACH-000971)
│   ├── KLA31_17: TALL-104 (GSE163787，Plan A 严格 bulk 库)
│   ├── KLA31_18: HepG2 亲本基线 (GSE158552 / DepMap ACH-000739)
│   ├── KLA31_19: A549 (GSE171750 / DepMap ACH-000681)
│   ├── KLA31_20: MDA-MB-468 (GSE157383 / DepMap ACH-000849)
│   ├── KLA31_21: T-47D (GSE283812 / DepMap ACH-000147)
│   ├── KLA31_22: PC-3M 细胞对照 (GSE235595)
│   ├── KLA31_23: 胶质母细胞瘤干细胞 MES28 载体对照 (GSE266884)
│   └── KLA31_24: RKO 亲本基线 (GSE318640 / DepMap ACH-000943)
└── 正常细胞系与培养模型 (Normal cell lines/models, n=7)
    ├── KLA31_25: HEK293T 未处理野生型 (GSE203529)
    ├── KLA31_26: HMC3 溶剂对照 (GSE275256)
    ├── KLA31_27/28: HK-2 共享未处理对照 (GSE240748)
    ├── KLA31_29: MCF10A 对照 (GSE103520 / recount3)
    ├── KLA31_30: 神经干细胞 NSC 对照模型 (GSE119834)
    └── KLA31_31: HUVEC 全细胞基线 (GSE203551 / ENCODE)
```

> [!IMPORTANT]
> **关于括号中的 `DepMap ACH-*`**：这是**蛋白组侧**用于锁定细胞系身份的模型编号，**不是 RNA 来源**；本流程的 RNA 全部来自 GSE / GTEx / TCGA / recount3。
>
> **共享参照显式追踪**：  
> 三个 HCT116 组（KLA31_14/15/16）以及两个 HK-2 组（KLA31_27/28）分别共享同一份转录组基线，因此实际包含 **28 个独立参考矩阵**。在下游归一化及膨胀至 31 组时，通过 [`outputs/20260918_qsmooth_31group_hgnc/group_expansion_31.csv`](../outputs/20260918_qsmooth_31group_hgnc/group_expansion_31.csv) 显式记录 `ReferenceKey`，杜绝重复计算伪自由度。（2026-09-18 起改用这一版：`20260916_qsmooth_31group` 由符号映射尚未补入 HGNC 旧名的矩阵算出，缺失全部改名基因。）

---

## 3. 上游数据获取、注释标尺构建与矩阵标准化处理 (Upstream Processing)

本章节系统阐述从公开发表的原始数据检索、服务器端流式提取、外显子并集基因长度构建、跨版本转录本映射、到表达量代数换算与组内质量控制（QC）的完整上游技术细节。

```mermaid
flowchart TD
    subgraph "Stage 0: 长度与 ID 基准构建"
        A1["Ensembl Release 111 GRCh38 GTF"] -->|"awk 提取外显子 + 一维区间并集融合"| A2["human_gene_lengths_ensembl111.tsv<br>(唯一外显子并集非重叠长度标尺)"]
        A3["Ensembl Release 87 GRCh37 GTF"] -->|"awk 解析 transcript_id & gene_id"| A4["human_grch37_transcript_to_gene.tsv<br>(跨组装版本 ENST 映射表)"]
        A5["NCBI Homo_sapiens.gene_info.gz<br>(taxid 9606)"] -->|"官方单射过滤 (剔除1对多歧义)"| A6["human_symbol_to_entrez.tsv<br>(仅用于无法避免 Symbol 的作者矩阵)"]
    end

    subgraph "Stage 1: 多源数据获取与解析 (1,898 样本)"
        B1["GTEx v8 (578 肺, 197 海马)"] -->|"awk 流式提取 GCT 目标列"| M1["GTEx Counts"]
        B2["TCGA / GDC API (ESCC 95, PRAD 501, LIHC 371/50)"] -->|"多核并行提取 unstranded counts, 过滤 ^N_"| M2["TCGA STAR Counts"]
        B3["GEO / recount3 / SRA (MCF10A)"] -->|"sra.gene_sums 覆盖度累加 (4 Runs/样本)"| M3["Recount3 Coverage"]
        B4["GEO 胎盘 GSE114691 (21 例)"] -->|"ENST 转录本计数通过 A4 映射并 rowsum 累加"| M4["Placenta Gene Counts"]
        B5["GEO 作者 RSEM / Counts / FPKM / RPKM (BPH, TALL-104 等)"] -->|"tar / openpyxl / 官方 Entrez 映射"| M5["GEO 各类表达矩阵"]
    end

    subgraph "Stage 2: 矩阵标准化与代数统一"
        M1 & M2 & M3 & M4 & M5 --> C1{"数值物理属性判定"}
        C1 -->|"加性计数 (Counts / Coverage / RSEM)"| C2["基于 A2 长度标尺换算 Rate:<br>Rate = Counts / (Length_bp / 1000)"]
        C1 -->|"分数型 (FPKM / RPKM)"| C3["基于 A2 长度标尺再归一化 Rate"]
        C2 & C3 --> C4["缩放为每百万转录本 (TPM):<br>TPM = Rate / Σ(Rate) × 1e6"]
        C4 --> C5["稳健对数变换 (假计数 0.5):<br>Y = log2(TPM + 0.5)"]
    end

    subgraph "Stage 3: 质控审计与标准化资产生成"
        C5 --> D1["组内样本 QC: Spearman 秩相关 (>0.90), 零值率, 库深度"]
        D1 --> D2["28 个独立标准化矩阵 (.tsv.gz & .rds)"]
        D1 --> D3["审计底册 (group_manifest, group_sample_qc, group_gene_summary)"]
    end
```

### 3.1 原始多源数据检索、获取与路由策略

针对 31 类生物学材料对应的 28 个参考组，项目建立了严格的分类接入路由协议（`handlers` 调度架构，见 [`workflow/build_31_expression_matrices_20260916.R`](../workflow/build_31_expression_matrices_20260916.R)）：

| 数据来源与接入 Handler | 代表性材料组 | 样本数 (n) | 原始数据形态 | 提取与处理技术细节 |
|---|---|---|---|---|
| **GTEx v8** (`h_gtex_v8`) | KLA31_02 (正常肺)<br>KLA31_05 (海马脑区) | 肺: 578<br>海马: 197 | `GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz` (1.2 GB) | 1. 读取样本属性表 `SampleAttributesDS.txt`，依 `SMTSD` 提取目标组织条目；<br>2. 剔除仅有属性登记但该 release 无测序行的子切片（海马剔除 46 个，肺剔除 289 个）；<br>3. 动态生成 `awk` 流式过滤脚本，单遍流式读取 GCT，提取目标列，避免耗尽宿主机内存。 |
| **TCGA / GDC API** (`h_gdc_star`) | KLA31_09 (癌旁正常肝)<br>KLA31_10 (食管鳞癌 ESCC)<br>KLA31_11 (前列腺癌 PRAD)<br>KLA31_12 (肝细胞癌 HCC) | 正常肝: 50<br>ESCC: 95<br>PRAD: 501<br>HCC: 371 | GDC augmented STAR gene counts (GENCODE v36, TSV per sample) | 1. 严格锁定组织学分类：TCGA-ESCA 中仅提取病理确诊的 **ESCC（鳞状细胞癌）**，坚决剔除腺癌（Adenocarcinoma），保证材料真实对齐；<br>2. 排除 `Solid Tissue Normal` 以外的重复 aliquot；<br>3. 使用 `mclapply` 多核并行读取单样本 TSV，利用 `colClasses` 跳过冗余列仅加载 `unstranded` 计数；<br>4. 过滤 `^N_` 质控行（`N_unmapped`, `N_multimapping`, `N_noFeature`, `N_ambiguous`）。 |
| **recount3 / SRA** (`h_recount3_runs`) | KLA31_29 (MCF10A 对照) | 3 | `sra.gene_sums.SRP117021.G026.gz` | 1. 跳过开头的 `##annotation=` 元数据注释头；<br>2. recount3 的 `sra.gene_sums` 记录的是单碱基覆盖度总和（base-level coverage sums），每个样本包含 4 条独立的测序 Run（SRR），在样本内对 Runs 进行精确行累加。 |
| **GEO 跨版本转录本计数** (`h_transcript_counts`) | KLA31_06 (正常胎盘) | 21 | `GSE114691_MasterCount_ControlONLY.txt.gz` | 原作者基于 GRCh37 的转录本定量（ENST）。调用 Ensembl GRCh37.87 转录本映射表汇聚为 ENSG（详见 3.3 节）。 |
| **GEO / NCBI 稳定 Counts** (`h_ncbi_counts`) | KLA31_01 (肌腱)<br>KLA31_07 (精子)<br>KLA31_13 (MCF7)<br>KLA31_18 (HepG2)<br>KLA31_20 (MDA-MB-468)<br>KLA31_22 (PC-3M)<br>KLA31_23 (MES28)<br>KLA31_27/28 (HK-2)<br>KLA31_30 (NSC)<br>KLA31_31 (HUVEC) | 3 ~ 8 每组 | NCBI 官方统一重定量矩阵（`raw_counts_GRCh38.p13_NCBI.tsv.gz`） | 1. 基于 NCBI 官方 `gene2ensembl` 表（taxid 9606）将数字 Entrez GeneID 映射至 Ensembl ENSG；<br>2. 丢弃 1 对多歧义 Entrez ID；多对一由 `rowsum` 累加。 |
| **GEO 单一 Ensembl Counts** (`h_single_ensembl_counts`) | KLA31_14/15/16 (HCT116)<br>KLA31_24 (RKO) | 各 3 | 作者 raw counts 矩阵 | 去除版本后缀（如 `.15`），通过 `rowsum` 对相同无版本 ENSG 进行折叠合并。 |
| **GEO 作者 Symbol Counts** (`h_single_symbol_counts`) | KLA31_19 (A549)<br>KLA31_21 (T-47D) | 各 3 ~ 4 | GSE171750 RSEM counts<br>GSE283812 raw counts | 仅在原作者未提供稳定 ID 时调用。经 NCBI 官方 `gene_info` 严格单射过滤为 Entrez，再映射至 ENSG，杜绝直接使用 Symbol（详见 3.3 节）。 |
| **GEO RSEM 分数型计数** (`h_rsem_per_sample`) | KLA31_25 (HEK293T 未处理) | 3 | `GSE203529_RAW.tar` 内单样本 RSEM 文件 | 从 tar 归档提取解压，读取 `expected_count`（EM 算法估计的分数型计数值），基于全集外显子并集长度换算 TPM。 |
| **GEO FPKM / RPKM** (`h_per_sample_rpkm`, `h_xlsx_fpkm`, `h_single_ensembl_fpkm`) | KLA31_03/04 (瘢痕及皮肤)<br>KLA31_08 (BPH)<br>KLA31_17 (TALL-104) | 1 ~ 18 | 压缩包单样本 RPKM (BPH)<br>Excel 工作簿 (GSE181540)<br>作者 CSV (GSE163787) | 1. 对 GSE181540，因服务器无 `openpyxl`/`pandas`，由 Python 标准库 `zipfile` + `ElementTree` 直接解析 xlsx 导出 TSV；<br>2. 对 BPH，解压 tar.gz 提取 18 例 transition-zone 样本 RPKM 并取共有基因；<br>3. 全部基于统一 Ensembl 111 长度标尺再归一化为 TPM。 |

---

### 3.2 统一基因外显子并集长度标尺构建（Stage 0：Merged-Exon Length Reference）

在 RNA-seq 中，从原始 Counts 或 FPKM/RPKM 换算 TPM 时，**基因长度（Gene Length）的选择直接决定了表达定量的准确性**：
- **致命缺陷警示（为什么不能直接累加转录本外显子？）**：一个基因往往存在多条转录异构体（Isoforms），这些异构体高度共享外显子。若简单将各转录本的长度直接相加，会导致共享外显子被重复计算数次，造成基因长度虚高数倍，严重低估高转录异构体基因的 TPM。
- **外显子区间并集融合算法（Exon Union / Merged Exon）**：  
  在 [`workflow/server_prepare_gene_annotation_20260916.sh`](../workflow/server_prepare_gene_annotation_20260916.sh) 中，以官方 Ensembl Release 111 GRCh38 GTF（`Homo_sapiens.GRCh38.111.gtf.gz`）为唯一起点：
  1. 提取第 3 列为 `exon` 的特征行，提取第 9 列属性中的 `gene_id`，剔除版本号与修饰后缀；
  2. 按照 `gene_id`（字典序）与染色体起始坐标（数值升序）排序：`sort -k1,1 -k2,2n`；
  3. 执行一维区间并集融合算法（1D Interval Merging）：
     ```awk
     $1 != g {
       if (g != "") print g"\t"len
       g = $1; cs = $2; ce = $3; len = $3 - $2 + 1; next
     }
     {
       if ($2 <= ce) { if ($3 > ce) { len += $3 - ce; ce = $3 } }
       else          { len += $3 - $2 + 1; ce = $3 }
     }
     END { if (g != "") print g"\t"len }
     ```
  4. 生成 [`metadata/annotation/human_gene_lengths_ensembl111.tsv`](../metadata/annotation/human_gene_lengths_ensembl111.tsv)（包含 63,241 个 Ensembl 基因的物理非重叠外显子总长度，单位 bp）。
- **唯一物理基准地位**：全套 31 组数据中所有涉及“从 Counts 算 TPM”或“从 FPKM/RPKM 重算 TPM”的步骤，**全部且唯一使用该表**，彻底抹平了不同作者、不同历史版本 GTF 带来的基因长度系统偏差。

---

### 3.3 跨版本转录本映射与跨命名空间防歧义映射机制

#### 3.3.1 GSE114691 胎盘跨组装版本映射（GRCh37 ENST $\to$ GRCh38 ENSG）
- **问题**：KLA31_06（胎盘组织，GSE114691）原作者仅发布了基于 GRCh37 的转录本水平计数（ENST）。若直接使用当前 NCBI `gene2ensembl` 映射，因当前表仅保留含有 RefSeq RNA 对应关系的转录本，**仅能解析该矩阵中约 23% 的旧版 ENST 编号**。
- **解决方案**：在 [`workflow/server_prepare_grch37_transcript_map_20260916.sh`](../workflow/server_prepare_grch37_transcript_map_20260916.sh) 中，下载官方 Ensembl GRCh37 release 87 GTF（`Homo_sapiens.GRCh37.87.gtf.gz`），提取完整的 `transcript_id` $\to$ `gene_id` 映射表（`human_grch37_transcript_to_gene.tsv`）。
- **映射效能**：
  - 源矩阵共有 204,940 个转录本；
  - 成功映射 **196,501 个转录本（覆盖率高达 95.9%）**，顺利汇集至 **57,905 个 Ensembl 基因**；
  - 丢失的未映射转录本仅占总测序计数的 **3.08%**。
  - 由于 Ensembl Gene ID（`ENSG...`）在 GRCh37 与 GRCh38 组装版本间保持主键稳定，因此通过基因级汇聚，实现了跨基因组版本的无损对齐。

#### 3.3.2 严格隔离 Gene Symbol 的单射路由控制
- **问题**：在全部 28 个参考组中，仅有两个矩阵的作者原始文件以 Gene Symbol 作为行名：GSE171750 (A549, KLA31_19) 与 GSE283812 (T-47D, KLA31_21)。这两个矩阵写成于 2020 年前后，仍在使用**当时的**符号。
- **严苛规则**：根据生信分析红线，绝不允许直接基于 Symbol 进行任何下游分析或表连接。为此设计了逐级收紧的单射通道（`build_symbol_lookup()`，见 [`workflow/lib_kla31_expression_20260916.R`](../workflow/lib_kla31_expression_20260916.R)）：
  1. **第一级**：下载 NCBI 官方 `Homo_sapiens.gene_info.gz`（taxid 9606），构建官方 Symbol $\to$ Entrez 字典（`human_symbol_to_entrez.tsv`）；统计频数并**歧义剔除**——若一个 Symbol 对应多于 1 个官方 Entrez GeneID，**整行直接丢弃**，仅保留严格一对一单射。
  2. **第二级（2026-09-18 补入）**：官方符号只含**今天的**命名，因此凡被改名过的基因（`AARS`$\to$`AARS1`、`CSRP2BP`$\to$`KAT14`、`ADCK3`$\to$`COQ8A`、`ACPP`$\to$`ACP3`）在第一级一律解析失败并被静默丢弃。补入 HGNC 官方改名史（`hgnc_complete_set.tsv.gz` 的 `prev_symbol`）作为回退。**刻意不使用 `alias_symbol`**：alias 只是习惯简称，历史上可能被多个位点共用，按 alias 合并会把不同基因错误融合（实测会融合 `HNRNPU`+`HNRNPU-AS1`、`CT45A3`+`CT45A4`）。
  3. **收敛保护**：别名回退行若与同一矩阵内其它行落到**同一个 Ensembl 基因**，则该行**丢弃而不合并**。保护必须跑在 `rowsum()` 实际分组所用的 **Ensembl 键**上而非 GeneID——多个 Entrez 可共享同一 ENSG，只在 GeneID 层设防仍会让别名行并进无关基因、虚增其计数（实测 `ENSG00000267858` 即为此情形）。
  4. 映射结果：GSE171750 中 **23,830 / 26,475** 行成功映射（23,771 基因）；GSE283812 中 **19,159 / 19,308** 行（19,154 基因）。
- **为何这一条至关重要（2026-09-18 修正）**：共同基因空间是 28 个矩阵的**交集**，因此**一个矩阵丢基因 $=$ 全部下游分析丢基因**。未经第二级回退时，仅 A549 的符号丢弃就使 **764 个基因**被排除在全部跨材料比较之外，T-47D 再排除 104 个；`AARS1`（最主要的乳酸化 Writer）即在其中。修正后共同基因空间由 **17,340 增至 18,332**（$+992$）。

#### 3.3.3 歧义折叠的代数与统计学法则
在将源数据映射至无版本 Ensembl ID（`strip_ensembl_version`）时，项目严格遵守数值的代数物理意义（见 [`workflow/lib_kla31_expression_20260916.R`](../workflow/lib_kla31_expression_20260916.R)）：
- **加性物理量（Counts / Recount3 Coverage Sums）**：
  当多个源转录本、外显子或历史 ID 汇聚至同一个 ENSG 时，其物理含义为落入同一基因的测序读数片段累加。因此采用 `rowsum()` 进行精确线性求和：
  $$\text{Counts}_{\text{ENSG}} = \sum_{t \in \text{ENSG}} \text{Counts}_t$$
- **比例与分数量（FPKM / RPKM / TPM）**：
  分数量的分母已包含了文库大小和基因长度，**在数学上不可相加，也不可简单求平均**。流水线规定：凡是遇到多个条目映射至同一个 ENSG 的分数量，**直接剔除该歧义 ENSG**（`n_dropped_ambiguous`），绝不在分数值上做非法代数加和。

---

### 3.4 统一表达量换算数学标准

为确保来自不同定量流水线的表达谱在同一物理尺度下可比，所有矩阵在 Stage 2 均标准化为统一的 $\log_2(\text{TPM} + 0.5)$ 空间。

#### 3.4.1 原始 Counts $\to$ TPM 算法
对任意基因 $g$ 与样本 $s$：
1. 长度标尺换算：
   $$\text{Length\_kb}_g = \frac{\text{Length\_bp}_g}{1000}$$
2. 每千碱基读数速率（Rate）：
   $$\text{Rate}_{g,s} = \frac{\text{Counts}_{g,s}}{\text{Length\_kb}_g}$$
3. 相对丰度缩放（TPM）：
   $$\text{TPM}_{g,s} = \frac{\text{Rate}_{g,s}}{\sum_{i} \text{Rate}_{i,s}} \times 10^6$$
   （保证每个样本的 $\sum_{g} \text{TPM}_{g,s} \equiv 1,000,000$）。

#### 3.4.2 FPKM/RPKM $\to$ TPM 算法
对于仅提供 FPKM/RPKM 的矩阵（如 GSE132714 BPH, GSE163787 TALL-104），由于 FPKM 本质上已除以长度，流水线使用相同的 Ensembl 111 外显子并集长度进行再标定：
$$\text{Rate}_{g,s}' = \frac{\text{FPKM}_{g,s}}{\text{Length\_kb}_g}, \quad \text{TPM}_{g,s} = \frac{\text{Rate}_{g,s}'}{\sum_{i} \text{Rate}_{i,s}'} \times 10^6$$

#### 3.4.3 特殊物理量（RSEM 与 recount3）的代数自洽性证明
- **RSEM expected_count**：RSEM 使用 Expectation-Maximization（EM）算法处理多重比对，其输出的 expected_count 本质为连续非负实数（带有小数）。算法直接以高精度浮点数带入 Counts $\to$ TPM 公式计算，归档时再记录四舍五入整数，保证数学期望无偏。
- **recount3 coverage sum**：`sra.gene_sums` 中的数值是覆盖到每个外显子碱基上的 Base Coverage 累加，数值规模约为 $10^9$（等于 $\text{ReadCounts} \times \text{ReadLength}$）。当将其代入 TPM 公式时：
  $$\text{Rate}_{g,s} = \frac{\text{Coverage}_{g,s}}{L_g} = \frac{\text{Counts}_{g,s} \times \text{ReadLength}}{L_g} = \text{ReadLength} \times \left( \frac{\text{Counts}_{g,s}}{L_g} \right)$$
  在计算 $\text{TPM} = \frac{\text{Rate}_{g,s}}{\sum_i \text{Rate}_{i,s}} \times 10^6$ 时，恒定的读长因子 $\text{ReadLength}$ 在分子与分母中**精确对消**！因此从 recount3 覆盖度直接换算的 TPM 与原始 Counts 换算的 TPM 完全等价，具备严格的数理自洽性。

#### 3.4.4 对数变换与伪计数设定
表达矩阵统一执行稳健对数变换：
$$Y_{g,s} = \log_2(\text{TPM}_{g,s} + 0.5)$$
- **伪计数取 0.5 的理论依据**：
  - 当基因未表达（$\text{TPM} = 0$）时，$Y = \log_2(0.5) = -1.0$，数值稳定且有明确下界，避免产生负无穷大（$-\infty$）；
  - 相比于传统加 1.0（产生 0），加 0.5 能更好地分离“真零表达”与“微弱背景噪音”，且能有效压缩极低表达区间的泊松随机波动。

---

### 3.5 样本级质量控制（QC）与组内生物学重复一致性检验

为杜绝任何离群、降解或标注错误的测序样本污染下游分析，在提取过程中针对全量 **1,898 个生物学样本** 执行了自动化质控流水线（`group_qc`）：

1. **核心 QC 指标体系**：
   - **`LibrarySize`**：文库测序深度总读数（GTEx 与 TCGA 中位深度 $> 5 \times 10^7$ reads）；
   - **`GenesDetected`**：有效检出基因数（实测跨度较大，13,564 ~ 55,681：GTEx/TCGA 深层队列偏高，单细胞系与小样本队列偏低）；
   - **`ZeroFraction`**：零表达基因比例（实测 8.1% ~ 64.9%）；
   - **`TPMTotal`**：总 TPM 校验（严格等于 $1.0 \times 10^6$，确保代数正确）。
2. **组内生物学重复一致性检验（Spearman Correlation）**：
   - 为避免全基因集中数万个零表达基因虚假拉高相关系数，QC 算法先筛选出在组内至少半数样本中 $\text{TPM} \ge 1.0$ 的稳健表达基因集（通常为 12,000 ~ 15,000 基因）；
   - 在该子集上计算组内所有生物学重复样本的两两 Spearman 秩相关系数矩阵（$r_s$），提取极小值 `PairwiseSpearmanMin` 与极大值 `PairwiseSpearmanMax`；
   - **检验结果（实测）**：组内重复一致性呈明显的来源依赖性。技术重复型细胞系极高（MCF10A **0.970 ~ 0.972**，HMC3 **0.986 ~ 0.988**）；而供体异质的临床队列自然偏低——**28 组中有 14 组极小值低于 0.95**，最低为 BPH 0.34、TCGA-PRAD 0.35、TCGA-LIHC 原发 0.44、NSC 0.52、ESCC 0.58、海马 0.59、精子 0.60、肺 0.70。**该指标反映材料本身的生物学异质性（不同供体、不同个体），不是数据质量缺陷**；样本可用性应以 `LibrarySize`、`ZeroFraction` 与 `TPMTotal` 为准，相关系数仅作来源特性描述。
3. **输出的标准化审计资产**（存放于 [`outputs/20260916_expression_extraction/`](../outputs/20260916_expression_extraction/)）：
   - `group_manifest.csv`：28 个独立参考组的完整元数据、源文件 MD5、提取路由与基因/样本统计；
   - `group_sample_qc.csv`：全量 1,898 个独立测序样本的单样本 QC 明细；
   - `group_gene_summary.csv`：每个基因在各参考组内的平均表达量与检出率；
   - `group_index.csv`：31 个材料类别与 28 个物理矩阵的对应索引字典。

---

## 4. 表达量统一化与跨组织平滑分位数归一化（qsmooth / YARN 算法体系）

### 4.1 共同基因全集交集锁定
不同公共数据源由于所用参考基因组注释（GENCODE v26 / v36 / Ensembl 111 / NCBI RefSeq）版本不同，各自覆盖的基因数从 19,154（T-47D GSE283812）到 61,365（HMC3 GSE275256）不等。  
在 28 个独立参考矩阵中取 Ensembl ENSG 全集交集，**最终锁定 18,332 个全局共有 Ensembl 基因**。所有 31 组材料的下游分析均在该封闭的 18,332 空间中进行。

> [!IMPORTANT]
> **交集瓶颈的准确定位（2026-09-18 修正）**：交集由两个覆盖面最窄的矩阵限定——T-47D 矩阵（19,154 基因）独占地排除 164 个基因，A549 矩阵再排除 74 个，其余 26 个矩阵各自排除 0 个。
>
> 在 HGNC 改名修正之前，本节曾把 17,340 归因于"T-47D 矩阵只有 18,815 基因"。那个归因**是错的**：T-47D 与 A549 偏低并非它们的作者测得少，而是其符号路由丢掉了全部改名基因。修正后两矩阵分别回升 1,411 与 339 个基因，共同基因空间随之由 17,340 增至 18,332。详见 3.3.2 节。

### 4.2 跨组织场景下传统分位数归一化的失效与平滑分位数算法
- **传统全局分位数归一化（Global Quantile Normalization, Global QN）的致命局限**：  
  Global QN 假设所有样本的整体基因表达分布（Empirical Quantile Distribution）在数学上完全一致，仅存在技术批次扰动。然而，本项目涵盖了**人类大脑（海马）、生殖组织（精子）、间叶组织（肌腱）、上皮性癌症与培养永生化细胞系**。不同组织类型之间存在极其显著的内源性转录组异质性（例如脑组织中特定转录本超高丰度表达、精子组织中广泛的转录本截短与静默）。使用 Global QN 会强制拉平这些真实的生物学峰度与偏度，人为抹杀组织特异性转录组特征。
- **平滑分位数归一化算法（qsmooth）**：  
  引入 Bioconductor 核心算法 `qsmooth`（v1.22.0，基于 YARN 框架）。算法将样本按生物学材料类别分组，对每个分位数级别 $q$，分别计算组内变异 $S_{\text{within}}(q)$ 与组间变异 $S_{\text{between}}(q)$，并拟合平滑加权函数 $w(q) \in [0, 1]$：
  - 当变异主要由实验技术噪音驱动时（组间变异小），$w(q) \to 1$，偏向全局分位数归一化，消除技术漂移；
  - 当变异主要由生物学真实组织差异驱动时（组间变异大），$w(q) \to 0$，保留各组织材料特有的固有分位数分布。
- **归一化诊断结果**：
  - 对 18,332 基因的权重诊断显示，权重中位数严格为 **0**（详见 `summary.csv`），证实大多数基因保留了其材料类别的固有分位数分布，完全符合预期。
  - **方案 A 与方案 B 跨尺度交叉验证**：
    - **方案 A（Collapsed Reference）**：每个参考材料组取样本中位数向量，构成 $18{,}332 \times 28$ 矩阵进行归一化；
    - **方案 B（Full Sample-Aware）**：全量 $18{,}332 \times 1{,}898$ 样本完整输入归一化后再按材料求中位数。
    - **结果**：方案 A 与方案 B 计算得到的跨组织基因表达向量，两两相关系数中位数高达 **0.9899**；每个基因表达值的绝对差值中位数仅为 **0.213 $\log_2$ 单位**，且全基因集中无任何基因的绝对差值大于 0.5。这无可辩驳地证实：方案 A 坍缩参考谱作为跨材料基准高度可靠且代数稳健。

---

## 5. DDR 修复通路与 UniProt-Ensembl 双向映射

### 5.1 严格无 Symbol 映射策略
1. **DDR 靶蛋白集合**：401 个 Kla ∩ DDR 靶蛋白。
2. **双重官方映射**：
   - 途径一：UniProt 官方 `idmapping.dat`（UniProtKB-AC $\to$ Ensembl Gene ID）；
   - 途径二：NCBI `gene2ensembl` 通过 Entrez GeneID 间接映射。
   - 交叉比对：在两种途径皆覆盖的 882 个 accession 中，879 个完全一致（99.66%）；3 个分歧记录在案并显式解析。
3. **覆盖度**：
   - 401 个 Kla ∩ DDR 蛋白质中，**394 个**成功覆盖到 18,332 表达谱空间（对应 **371 个非冗余 Ensembl Gene ID**）。未覆盖的 7 个蛋白，其基因在至少一个参考矩阵中缺失。
   - （2026-09-18 修正前为 377 个蛋白 / 357 个基因；差额来自 3.3.2 节的 HGNC 改名回退。）

### 5.2 DDR 修复通路体系（8 大分面）
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
> **关键逻辑修正**：S4 注释表中的方向状态包括 `+1`（上调促活）、`-1`（下调抑制）和 `0`（**未判定方向性，不代表该蛋白不属于此通路**）。早期代码的成员判据只排除了空字符串，因而**把 `0` 也当成了成员**，导致每个基因同时匹配全部 8 条通路、按顺序坍缩到第一条（BER），热图只剩单个分面。现已修正为要求状态为 `+1` 或 `-1` 才算该通路成员。

---

## 6. 分析内容与结果

### 6.1 三个来源在每一个（材料 × DDR 蛋白）对上的联合状态（11,501 对）

每个观测对同时有三个来源的测量：**乳酰化组（Kla）**、**配对的未富集全蛋白参照组（Ref）**、**该基因的表达水平（RNA）**。三者联合分布如下：

| Kla | Ref | 表达（TPM） | 对数 | 占比 |
|---|---|---|---|---|
| + | + | ≥ 1 | 2,651 | 23.1% |
| + | + | < 1 | 57 | 0.5% |
| + | − | ≥ 1 | 480 | 4.2% |
| + | − | < 1 | 28 | 0.2% |
| − | + | ≥ 1 | 3,935 | 34.2% |
| − | + | < 1 | 101 | 0.9% |
| − | − | ≥ 1 | 3,827 | 33.3% |
| − | − | < 1 | 422 | 3.7% |
| | | **合计** | **11,501** | **100%** |

### 6.2 以表达水平为主线看两套质谱的检出率

| 表达（TPM） | 对数 | 全蛋白组（Ref）检出率 | 乳酰化组（Kla）检出率 |
|---|---|---|---|
| ≥ 1 | 10,893 | **60.5%** | **28.7%** |
| < 1 | 608 | **26.0%** | **14.0%** |

表达阈值取 **TPM ≥ 1**（项目质控所用的同一判据），而非组内排名——组内排名会强制每组得到相同的高/低比例，抹平材料间差异。

各联合状态对应的转录本水平中位数：

| 状态 | 对数 | $\log_2(\text{TPM}+0.5)$ 中位数 |
|---|---|---|
| Kla+ / Ref+ | 2,708 | 5.14 |
| Kla− / Ref+ | 4,036 | 4.56 |
| Kla+ / Ref− | 508 | 4.12 |
| Kla− / Ref− | 4,249 | 3.53 |

> [!NOTE]
> 以上为**检出状态与表达水平的联合描述**，未作统计检验，也不对乳酰化选择性作机制推断。
> 关于"哪个来源提供了什么"的说明见 1.1 节：`Kla−/Ref+`（蛋白在全蛋白组检出、未在乳酰化组检出）这一状态由蛋白质组自身即可判定，转录组补充的是 `Kla−/Ref−`（两套质谱均未检出）内部"转录本本身偏低"与"转录本不低但未被捕获"的区分。

### 6.3 阴性对照：基于计数的"超额乳酰化率"受质谱深度主导

探索阶段曾计算各材料类别的"超额乳酰化率"（实际检出数 − 按表达预期检出数）。混杂因素分析结果：

- 该指标与蛋白质组的**质谱检测深度**决定系数 $R^2 \approx 0.80$（Spearman 0.930）；
- 回归移除深度后，四大生物学类别之间的差异消失（正常组织 −0.027 vs 癌细胞系 +0.033，Wilcoxon $p = 0.31$）。

**记录**：该指标不成立，不作为分析结论使用（见 [`outputs/20260917_rna_assisted_ddr/validity_vs_depth.csv`](../outputs/20260917_rna_assisted_ddr/validity_vs_depth.csv)）。项目原有的**比值型**指标（`KlaDdrFraction`、`DdrFractionPercentagePointDifference`）经同一检验不受深度显著影响（$\rho$ 0.19 ~ 0.22，$p$ 0.24 ~ 0.31），比率自带归一化。


---

## 7. 逐图方法学索引（每张图各自用了什么）

前六节描述的是从原始数据到封闭基因空间的**共用路径**；本节列出每张交付图**特有的**计算方法、统计检验与呈现约定。所有图均读取同一份 `qsmooth_B`（样本级，18,332 × 1,898）或 `qsmooth_A`（组级，18,332 × 28 → 膨胀至 31）。

### 7.1 图 1a 系列：DDR 基因与乳酸化靶基因（脚本 3A–3F）

| 图 | 数据来源 | 每点/每组的统计量 | 检验 |
|---|---|---|---|
| `Figure_1a_RNA_DDR_expression_by_tissue_boxplot` | `qsmooth_B` 样本级 | 每个**样本**取 371 个 DDR 基因表达的**中位数**（非均值） | 无（描述性） |
| `..._faceted_boxplot` | 同上 | 同上，按四大类分面 | 无 |
| `Figure_1a_RNA_DDR_and_lactylated_gene_expression_31groups_boxplot` | `qsmooth_A` 组级 | 每材料的组内中位数（Scheme B，31 点）；DDR（371）与 Kla 靶基因（5,224）并列 | **双因素 ANOVA**：`MedianExpression ~ Category * GeneSet`，报 Category 主效应 |
| `Figure_1a_RNA_DDR_fraction_31groups_boxplot` | `qsmooth_B` → 组内中位数 | DDR 表达基因数 ÷ 该材料表达基因总数 × 100% | **单因素 ANOVA** `~ Category`（31 组） |
| `Figure_1a_companion_*` | `qsmooth_B` 样本级 | 全 1,898 样本 | 单因素 / 双因素 ANOVA（样本级） |
| `Figure_1a_companion_RNA_DDR_fraction_by_tissue_boxplot` | `qsmooth_B` | 分组织展开 | 无 |

### 7.2 图 1b 系列：细胞增殖（脚本 1A–1F）

| 图 | 数据来源 | 每点/每组的统计量 | 检验 |
|---|---|---|---|
| `Figure_1b_RNA_proliferation_by_tissue_sample_boxplot` | `qsmooth_B` 样本级 | 每个**样本**取 20 个增殖标志基因表达的**均值**（非中位数） | 无（描述性） |
| `..._faceted_boxplot` | 同上 | 同上，四大类分面 | 无 |
| `Figure_1b_RNA_proliferation_score_31groups_boxplot` | `qsmooth_B` → 组内**中位数** | Scheme B，31 点 | **单因素 ANOVA** `~ Category` |
| `Figure_1b_RNA_proliferation_ranking_31groups` | 同上 | 按组中位数排序的柱状图 | 无 |
| `Figure_1b_companion_RNA_proliferation_score_sample_boxplot` | `qsmooth_B` | 全 1,898 样本 | 单因素 ANOVA（样本级） |
| `Figure_1b_companion_RNA_MKI67_over_{H3C1,ACTB,TUBB}` | `qsmooth_B` | 样本级 `MKI67 / 参照基因` 比值，**log10** 纵轴 | 单因素 ANOVA on `log10(Ratio)` |

> [!WARNING]
> **增殖与 DDR 的汇总统计量不一致**：增殖打分用 20 基因的**均值**（`colMeans`），DDR 表达用 371 基因的**中位数**（`colMedians`）。两者不是同一种中心趋势度量，跨这两张图直接比较基因集"表达高低"并不严格可比。这是刻意保留还是历史遗留，需要作者确认。
>
> **MKI67 比值图是旧算法**：增殖已改用 20 基因标志物集打分（见 1B/1C），比值图仅作备查保留，不作为主结论。

### 7.3 图 3c 系列：乳酸化调控因子热图（脚本 2A–2C）

| 图 | 数值 | 色彩 |
|---|---|---|
| `Figure_3c_RNA_regulator_qsmooth_heatmap` | `qsmooth_A` 原始 `log2(TPM+0.5)` | 墨绿/青绿连续渐变 |
| `..._zscore_heatmap` | 每个基因在 31 组内标准化 | 蓝-白-红发散，**显示时截断至 ±2.5** |
| `Figure_3c_RNA_regulator_percentiles`（备查） | 组内相对百分位（**已废弃**） | 渐变 |

- 48 个唯一调控因子；`HDAC8` 在源注释表中同时属于 Eraser 与 Writer-Eraser，故被画成 **2 列（共 49 列）**——这是忠实于源表，不是重复计数错误。
- 10 个重点基因（AARS1, ACAT2, KRT18, SIRT2, PARK7, HDAC1, HDAC2, BRD4, SMARCA4, TRIM33）加连续框选。

### 7.4 图 2 系列：RNA 辅助解析（`plot_rna_reference_31group_20260917.R`）

| 图 | 方法 |
|---|---|
| `RNA_1_DDR_panel_expression_31groups` | 371 基因 × 31 材料热图，按 DDR 通路分面；**行内标准化**（`t(scale(t(...)))`），故展示的是相对模式而非绝对丰度 |
| `RNA_2a_three_source_joint` | 11,501 个（材料 × DDR 蛋白）对的四象限联合状态，按 qsmooth TPM ≥ 1 拆分 |
| `RNA_2b_detection_by_transcript_level` | 同一批对，从 RNA 侧看两套质谱在"表达 / 低表达"上的检出率 |
| `RNA_2` | 2a + 2b 合成图 |

> [!NOTE]
> **"表达"判据是平滑后的**：`RNA_2*` 与所有 `Figure_1a` 占比图的表达调用均为 `qsmooth log2(TPM+0.5) ≥ log2(1.5)`，即在**跨组织平滑后的矩阵**上判定，而非原始单样本 TPM。该判据与项目 QC 一致，但严格说是"平滑 TPM ≥ 1"，图注已据此措辞。

---

## 8. 特殊操作与已知脆弱点清单

以下为流程中**非标准、人工干预或易被误读**的环节。前四类为方法学上必须披露的；后两类为实现层面的工程债。

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

### 8.3 下游管线中的隐式替换（已全部清除）

> [!IMPORTANT]
> **Figure 3c 热图曾写入伪造值。** 2026-09-18 之前，`AARS1` 与 `CSRP2BP` 因符号映射缺陷不在平滑矩阵中，脚本依次回退到：①未平滑的组均值；②**三个无关细胞系（MCF7/HCT116/HepG2）的均值**。A549 列因此携带一个凭空的 `6.297704`（实测值 `8.443779`），且这两行与相邻基因**不同标度**。该回退已删除，改为**硬报错**。详见 3.3.2 节与提交 `63ba6c9`。

其余已清除或需注意的实现层操作：

- **重建两个矩阵后所有既有基因的 TPM 被重新标定**：分母纳入了新找回的基因，A549 全体 ×0.960286、T-47D ×0.988064（在所有表达基因上恒定到 6 位小数），**非**逐基因扰动。
- **零方差基因在两张热图中被画成 z = 0**：`t(scale(t(...)))` 产生的 `NaN` 被 `scaled[!is.finite(scaled)] <- 0` 填为 0，与真实的"均值"不可区分。影响 `RNA_1` 与 `teacher_questions_rna_20260917.R` 的调控因子热图（`Figure_3c` 不受影响）。
- **`h_per_sample_rpkm` 忽略契约样本清单**：BPH 组按 tar 成员名中的字面量 `"BPH_"` 选样，契约里列出的 18 个 GSM 实际未参与筛选。
- **构建缓存不校验源文件哈希**：`if (file.exists(obj_path)) cached; next`。源文件变更后必须手工删除 `.rds` 才会重建——这正是 2026-09-18 需要备份并手工清理两个对象的原因。
- **`build_symbol_lookup()` 在 HGNC 文件缺失时仅告警并退回官方符号表**，即重新回到会丢失全部改名基因的行为。已加强为 `warning()`，但**仍是软失败**；若要绝对安全应改为 `stop()`。

### 8.4 已修正的正确性缺陷（2026-09-18）

| 缺陷 | 影响 | 修正 |
|---|---|---|
| 两个符号路由矩阵丢弃全部改名基因 | 共同基因空间少 992 个基因（含 AARS1） | HGNC `prev_symbol` 回退 + Ensembl 键收敛保护 |
| 热图伪造回退值 | 图 3c 中 A549 列一个凭空数值 | 删除回退，改硬报错 |
| 别名合并会融合不同位点 | 若用 `alias_symbol` 会融合 `HNRNPU`+`HNRNPU-AS1` 等 | 只用 `prev_symbol` |
| `overlay_summary.csv` 的 `shared_reference_rows` 恒为 0 | 已发布边车里一个错误数字（应为 5） | Python 端 `"True"` vs R 端 `TRUE` 大小写修正 |
| DDR 基因数硬编码为 357 | 图注与数据不符（实为 371） | 改为 `length(ddr_in_qb)` 动态计算 |
| 31 行底册与 `verification_summary.csv` 仍是修正前的基因数 | 与 `group_manifest.csv` 自相矛盾 | 重新生成 |
