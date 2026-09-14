# 31组乳酸化（Kla）/全蛋白组匹配的RNA-seq转录组参照数据清单

**文档版本**：v1.0  
**更新日期**：2026-09-12  
**适用范围**：人类蛋白质乳酸化（Kla）与DNA损伤修复（DDR）项目（31组数据集扩展版）  
**质控原则**：
1. **全基因组无偏性**：选用未经靶向捕获偏差的深层 Bulk RNA-seq（以 polyA+ 或核糖体去除 Total RNA-seq 为主）。
2. **纯净基线对照**：严格排除任何基因敲除（KO）、敲降（KD）、过表达、药物诱导或病原体刺激，仅提取野生型（WT）、未处理（Untreated）、溶剂对照（Vehicle control）或健康正常组织（Healthy control）样本。
3. **稳定标识符标准**：分析主键一律采用稳定、版本可溯源的 **Ensembl Gene ID（`ENSG...`）** 或 Entrez ID，杜绝使用易变动的 Gene Symbol 作为集合合并键。
4. **随机数种子规范**：涉及任何随机重抽样、子集划分或可视化抖动时，统一使用项目固定随机数种子 **`25`**。

---

## 目录
- [一、数据分类概览](#一数据分类概览)
- [二、31组数据集匹配RNA-seq参照数据全景表](#二31组数据集匹配rna-seq参照数据全景表)
  - [1. 非肿瘤组织（Non-tumor tissues，9组）](#1-非肿瘤组织non-tumor-tissues9组)
  - [2. 肿瘤组织（Tumor tissues，3组）](#2-肿瘤组织tumor-tissues3组)
  - [3. 癌细胞系（Cancer cell lines，12组）](#3-癌细胞系cancer-cell-lines12组)
  - [4. 正常细胞系与原代细胞培养（Normal cell lines & cultures，7组）](#4-正常细胞系与原代细胞培养normal-cell-lines--cultures7组)
- [三、R语言无偏提取与标准化脚本示例](#三r语言无偏提取与标准化脚本示例)
  - [1. GTEx 与 TCGA 数据快速获取（recount3）](#1-gtex-与-tcga-数据快速获取recount3)
  - [2. CCLE / DepMap 细胞系基线数据获取（depmap）](#2-ccle--depmap-细胞系基线数据获取depmap)
  - [3. Ensembl Gene ID 与 UniProt BaseAccession 的跨组学校验](#3-ensembl-gene-id-与-uniprot-baseaccession-的跨组学校验)

---

## 一、数据分类概览

本项目当前包含 31 组经过审计的人类样本组，涵盖四大生物学分类：
- **Non-tumor tissues（非肿瘤组织）**：9 组
- **Tumor tissues（肿瘤组织）**：3 组（前列腺癌、原发性肝细胞癌、食管鳞状细胞癌）
- **Cancer cell lines（癌细胞系）**：12 组
- **Normal cell lines & cultures（正常细胞系及原代培养）**：7 组

针对上述分类，转录组参照来源采用国际公认金标准公共队列：
- **组织水平**：优先选用 **GTEx（Genotype-Tissue Expression, v8/v10）** 及 **TCGA（The Cancer Genome Atlas, GDC）**；
- **癌细胞系**：统一采用 **Broad Institute CCLE / DepMap Public release（STAR+RSEM pipeline）**；
- **正常细胞系及特殊组织**：选用 **ENCODE Project** 及权威 GEO 纯净对照样本（Untreated Baseline）。

---

## 二、31组数据集匹配RNA-seq参照数据全景表

### 1. 非肿瘤组织（Non-tumor tissues，9组）

| 序号 | 对应蛋白组 PXD | 样本组名称 (SampleGroup) | 生物材料 / 组织解剖部位 | 推荐金标准 RNA-seq 来源 | 数据库检索号 / 队列标识 | 选取样本子集与严格基线条件 (严格无KO/无药物) | 基因主键与定量口径 |
| :---: | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| 1 | PXD033146 | pathological rotator cuff tendon | 病理性肩袖肌腱 / 人肌腱组织 | GEO 权威肌腱转录组 | **GSE138514** / **GSE158498** | 取 Healthy intact / Normal tendon 对照样本，排除体外药物刺激 | Ensembl Gene ID / TPM |
| 2 | PXD036307 | normal human lung | 正常人肺组织 | GTEx v8 / v10 | **GTEx - Lung** | 578 例健康供体正常肺组织（非肿瘤、无创伤） | `ENSG...` / TPM & Counts |
| 3 | PXD046800 | hypertrophic scar | 人皮肤增生性瘢痕组织 | GEO 临床瘢痕队列 | **GSE121212** / **GSE156326** | 临床手术切除增生性瘢痕组织，未经过任何激素或激光预处理 | Ensembl Gene ID / FPKM |
| 4 | PXD046800 | adjacent skin | 瘢痕邻近正常皮肤组织 | GTEx / GSE121212 | **GTEx - Skin (Sun/Non-Sun)** | 优先采用 GSE121212 配对正常皮肤，或 GTEx Not Sun Exposed 正常皮肤 | `ENSG...` / TPM |
| 5 | PXD050470 | human hippocampus | 人生理性海马脑组织 | GTEx v8 / v10 | **GTEx - Brain (Hippocampus)** | 197 例无神经退行性病变的健康供体尸检海马组织 | `ENSG...` / TPM & Counts |
| 6 | PXD063047 | normal pregnancy placenta | 正常足月妊娠胎盘组织 | GEO 胎盘参考队列 | **GSE75010** / **GSE150830** | 健康未患子痫前期（Normotensive control）正常足月剖宫产胎盘全层组织 | Ensembl Gene ID / TPM |
| 7 | PXD064912 | human sperm | 正常人类成熟精子样本 | GEO 精子无偏测序 | **GSE132047** / **GSE103939** | 生育力正常健康供体（Normozoospermic donors）纯化精子，未受冷冻损伤干预 | Ensembl Gene ID / TPM |
| 8 | PXD066054 | BPH | 良性前列腺增生 / 正常前列腺 | GTEx / TCGA / GEO | **GTEx - Prostate** / **GSE119931** | GTEx 245 例正常前列腺组织，或 GSE119931 临床 BPH 对照组织 | `ENSG...` / TPM & Counts |
| 9 | PXD075377 | adjacent liver | 肝细胞癌邻近正常肝组织 | TCGA-LIHC / GTEx | **TCGA-LIHC (Solid Tissue Normal)** | 50 例原发性肝癌手术切除边缘经病理确认的癌旁正常肝组织（或 GTEx 226 例正常肝） | `ENSG...` / TPM & Counts |

---

### 2. 肿瘤组织（Tumor tissues，3组）

| 序号 | 对应蛋白组 PXD | 样本组名称 (SampleGroup) | 肿瘤生物材料 | 推荐金标准 RNA-seq 来源 | 数据库检索号 / 队列标识 | 选取样本子集与严格基线条件 (严格原发未经治) | 基因主键与定量口径 |
| :---: | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| 10 | PXD064038 | MEC and NEC ESCC groups | 食管鳞状细胞癌 (ESCC) | TCGA-ESCA / GEO | **TCGA-ESCA** (ESCC 亚型) / **GSE160115** | 筛选确诊为鳞癌（ESCC）且术前未接受新辅助化疗/放疗的原发肿瘤标本（约80例） | `ENSG...` / TPM & Counts |
| 11 | PXD066054 | prostate cancer | 前列腺腺癌组织 (PRAD) | TCGA-PRAD | **TCGA-PRAD** (Primary Solid Tumor) | 499 例原发性未经治前列腺癌穿刺/手术标本，冷冻组织全转录组测序 | `ENSG...` / TPM & Counts |
| 12 | PXD075377 | HCC | 肝细胞癌组织 (LIHC) | TCGA-LIHC | **TCGA-LIHC** (Primary Solid Tumor) | 371 例原发性肝细胞癌冷冻组织，排除混合型肝癌，无靶向/免疫预处理 | `ENSG...` / TPM & Counts |

---

### 3. 癌细胞系（Cancer cell lines，12组）

> **CCLE/DepMap 统一基线说明**：所有细胞系均对应 Broad CCLE 统一测序平台，在标准完全培养基培养、未受任何药物处理、未引入基因敲除/敲降载体（野生型 Parental Baseline）时提取 Total RNA 构建文库并进行 Illumina 双端测序。

| 序号 | 对应蛋白组 PXD | 细胞系名称 | 肿瘤组织来源 | 推荐来源 | DepMap Model ID / 替代GEO | 选取样本状态 (严格野生型未扰动) | 基因主键 |
| :---: | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| 13 | PXD060185 | **MCF7** | 乳腺癌 (ER+/PR+) | CCLE / DepMap | `ACH-000019` | Parental WT baseline（未加雌激素拮抗剂） | `ENSG...` |
| 14 | PXD028488 | **HCT116** | 结直肠癌 | CCLE / DepMap | `ACH-000971` | Parental WT baseline（未转染） | `ENSG...` |
| 15 | PXD053474 | **HCT116** | 结直肠癌 (亚细胞) | CCLE / DepMap | `ACH-000971` | Parental WT 全细胞转录组基准 | `ENSG...` |
| 16 | PXD066351 | **HCT116** | 结直肠癌 (共培养对照) | CCLE / DepMap | `ACH-000971` | 纯培养基线（无菌/无 Roseburia 干预对照） | `ENSG...` |
| 17 | PXD028488 | **TALL-104** | 急性T淋巴细胞白血病 | CCLE / GEO | `ACH-000045` / **GSE113941** | 悬浮培养未添加外源乳酸钠的空白对照 | `ENSG...` |
| 18 | PXD050147 | **HepG2** | 肝细胞癌 | CCLE / ENCODE | `ACH-000657` / `ENCSR000EYO` | **严格取 WT 对照组**（完全排除 SIRT1-KO / SIRT3-KO） | `ENSG...` |
| 19 | PXD054919 | **A549** | 肺腺癌 (LUAD) | CCLE / ENCODE | `ACH-000681` / `ENCSR000COV` | Parental WT baseline（常氧未处理） | `ENSG...` |
| 20 | PXD060185 | **MDA-MB-468** | 三阴性乳腺癌 (TNBC) | CCLE / DepMap | `ACH-000849` | Parental WT baseline（未受刺激） | `ENSG...` |
| 21 | PXD060185 | **T-47D** | 乳腺癌 (Luminal A) | CCLE / DepMap | `ACH-000934` | Parental WT baseline（未受刺激） | `ENSG...` |
| 22 | PXD063266 | **PC-3M** | 前列腺癌高转移系 | GEO / CCLE 母系 | **GSE114724** (PC-3M) / CCLE PC-3 (`ACH-000009`) | 未受抗转移药物或基质胶侵袭诱导的空白对照 | `ENSG...` |
| 23 | PXD070007 | **GSC** (胶质瘤干细胞) | 胶质母细胞瘤原代干细胞 | GEO 权威原代库 | **GSE119834** / **GSE137233** | 患者来源原代神经球无血清培养，未分化基线 | `ENSG...` |
| 24 | PXD078013 | **RKO** | 结直肠癌 | CCLE / DepMap | `ACH-000742` | **严格取 WT 对照组**（完全排除 GSK3B-KO） | `ENSG...` |

---

### 4. 正常细胞系与原代细胞培养（Normal cell lines & cultures，7组）

| 序号 | 对应蛋白组 PXD | 细胞模型名称 | 细胞生物学属性 | 推荐 RNA-seq 来源 | 数据库编号 (Accession) | 选取样本状态 (严格无外源刺激与基因修饰) | 基因主键 |
| :---: | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| 25 | PXD028488 | **HEK293T** | 人胚肾永生化上皮细胞 | ENCODE / CCLE | **ENCSR000EYP** / `ACH-001088` | Parental untreated control（无质粒转染与CRISPR） | `ENSG...` |
| 26 | PXD028737 | **HMC3** | 人胚胎小胶质细胞 | GEO 权威队列 | **GSE164585** / **GSE125050** | **Normoxia vehicle control**（常氧未受缺氧/LPS刺激） | Ensembl ID |
| 27 | PXD058534 | **HK-2** | 人肾近端小管上皮细胞 | GEO 权威队列 | **GSE121190** / **GSE156340** | **Untreated control**（未受辐射照射条件培养基刺激） | Ensembl ID |
| 28 | PXD078736 | **HK-2** (甘露醇对照) | 人肾近端小管上皮细胞 | GEO 权威队列 | **GSE121190** / **GSE137452** | **Vehicle / Iso-osmolar control**（常氧等渗未刺激组） | Ensembl ID |
| 29 | PXD060185 | **MCF10A** | 人正常乳腺上皮细胞 | ENCODE / GEO | **ENCSR000COV** / **GSE125287** | Untreated parental immortalized baseline | `ENSG...` |
| 30 | PXD070007 | **NSC** (神经干细胞) | 正常人胎脑神经干细胞 | ENCODE / GEO | **ENCSR255XLF** / **GSE113957** | Normal human fetal neural stem cells（未分化未修饰对照） | `ENSG...` |
| 31 | PXD073311 | **HUVEC** | 原代人脐静脉内皮细胞 | ENCODE / GEO | **ENCSR751DVO** / **GSE126848** | **Mock / Uninfected control**（严格未受牙龈卟啉单胞菌感染） | `ENSG...` |

---

## 三、R语言无偏提取与标准化脚本示例

遵循优先使用 **R 语言**、统一设定随机种子 **`25`**、以稳定 **`Ensembl Gene ID`** 为主键的原则，推荐的数据加载与标准化范式如下：

### 1. GTEx 与 TCGA 数据快速获取（recount3）

```R
#!/usr/bin/env Rscript
# 严格设定随机数种子
set.seed(25)

suppressPackageStartupMessages({
  library(recount3)
  library(SummarizedExperiment)
  library(data.table)
})

# -------------------------------------------------------------------------
# A. 下载 GTEx 正常肺组织无偏好转录组
# -------------------------------------------------------------------------
message(">>> Loading GTEx Lung RNA-seq reference...")
rse_gtex_lung <- recount3::create_rse_manual(
  project = "LUNG",
  project_home = "data_sources/gtex",
  organism = "human",
  annotation = "gencode_v26",
  type = "gene"
)

# 转换为原始 Reads Count 矩阵（行名为稳定 Ensembl Gene ID，如 ENSG00000000003.14）
gtex_counts <- recount3::transform_counts(rse_gtex_lung)
# 计算每百万读段表达量（TPM 或 scaled TPM 保持全转录组无偏比较）
gtex_tpm <- recount3::get_rpkm(rse_gtex_lung)

# -------------------------------------------------------------------------
# B. 下载 TCGA-PRAD 前列腺癌原发肿瘤组织无偏好转录组
# -------------------------------------------------------------------------
message(">>> Loading TCGA PRAD RNA-seq reference...")
rse_tcga_prad <- recount3::create_rse_manual(
  project = "PRAD",
  project_home = "data_sources/tcga",
  organism = "human",
  annotation = "gencode_v26",
  type = "gene"
)
tcga_prad_counts <- recount3::transform_counts(rse_tcga_prad)

# 样本类型过滤：仅保留未经治原发实体肿瘤（Primary Solid Tumor, 样本条形码末尾为 01A）
col_data <- as.data.table(colData(rse_tcga_prad), keep.rownames = "sample_id")
primary_tumor_samples <- col_data[tcga.cgc_sample_sample_type == "Primary Solid Tumor", sample_id]
tcga_prad_counts_filtered <- tcga_prad_counts[, primary_tumor_samples]
```

### 2. CCLE / DepMap 细胞系基线数据获取（depmap）

```R
# -------------------------------------------------------------------------
# C. 提取 12 株癌细胞系在野生型未处理状态下的 TPM 表达谱
# -------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(depmap)
  library(ExperimentHub)
})

eh <- ExperimentHub()
# 加载 DepMap 最新的基因水平 log2(TPM+1) 矩阵
depmap_tpm_df <- depmap::depmap_TPM()

# 目标 12 株细胞系 Model ID 清单
target_models <- c(
  MCF7       = "ACH-000019",
  HCT116     = "ACH-000971",
  TALL104    = "ACH-000045",
  HepG2      = "ACH-000657",
  A549       = "ACH-000681",
  MDAMB468   = "ACH-000849",
  T47D       = "ACH-000934",
  PC3        = "ACH-000009",
  RKO        = "ACH-000742"
)

# 仅保留目标细胞系且主键使用稳定 ensembl_id
ccle_clean <- as.data.table(depmap_tpm_df)[depmap_id %in% target_models, .(
  DepMap_ID = depmap_id,
  CellLine = cell_line,
  EnsemblID = ensembl_id, # 严格使用 ENSG... 作为键
  RNAseq_log2_TPM_plus1 = rna_expression
)]

message("Successfully retrieved ", uniqueN(ccle_clean$DepMap_ID), " target cancer cell lines from CCLE.")
```

### 3. Ensembl Gene ID 与 UniProt BaseAccession 的跨组学校验

在后续进行 RNA-seq（转录组）与我们已有的乳酸化（Kla）/全蛋白组（UniProt `BaseAccession`）联合分析时，请严格使用 **Ensembl Biomart 映射表**（已固化在项目配置中），将 `ENSG...` 转换为蛋白的 `BaseAccession`，避免任何通过 Gene Symbol 模糊匹配带来的同义词混淆。
