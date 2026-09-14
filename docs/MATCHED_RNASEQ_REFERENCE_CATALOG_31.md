# 31组 Kla / 全蛋白组的 bulk RNA-seq 参照候选（已审计）

**2026-09-14补充**：按同一解剖组织与取材类别重新核查，新增瘢痕/邻近皮肤GSE181540、PC-3M GSE235595及同研究MES28 GSE266884候选；见[补充报告](../outputs/20260914_reference_material_audit/REPORT.md)和[新版31组候选](../outputs/20260914_reference_material_audit/rnaseq_reference_candidate_31.csv)。以下为2026-09-12历史候选，新增版仍全部AnalysisReady=FALSE。

**2026-09-14数据选路补充**：按“解剖取材优先、统一处理其次”重新核实可用矩阵和公共重处理资源。31组的实际来源路由、明确排除的替代数据集和跨来源分析边界见[RNA来源筛选说明](RNASEQ_SOURCE_SELECTION_20260914.md)及[机器可读路由表](../outputs/20260914_reference_material_audit/rnaseq_source_routing_candidate_31.csv)。

更新：2026-09-12。**这是来源核查与候选注册表，不是已下载验收的表达矩阵。** 原Gemini版本存在物种、组织、细胞系和技术错误，已归档至 [原版](../audit/20260912_workspace_rnaseq/Gemini_catalog_original.md)，不能使用其编号或R代码直接分析。

当前生物学范围为31组（9非肿瘤组织、3肿瘤组织、12癌细胞组、7非癌细胞/培养模型）。31组对应22个Kla PXD与18个全蛋白参考PXD，不代表31组同受试者配对RNA。

[完整工作区审计](../audit/20260912_workspace_rnaseq/REPORT.md) · [原37个编号逐项核查](../audit/20260912_workspace_rnaseq/original_accession_audit.csv) · [31组机器可读注册表](../audit/20260912_workspace_rnaseq/rnaseq_reference_registry_31.csv) · [冻结组清单](../audit/20260912_workspace_rnaseq/groups_31.csv)

证据分级：material_verified＝已核实人类bulk及材料/样本；model_id_verified＝模型ID正确但尚未锁表达release/profile；cohort_candidate＝大队列方向合理但具体样本待筛；approximate＝仅近似材料/模型；exact_match_gap＝精确匹配仍缺。所有记录仍需表达矩阵及样本QC，AnalysisReady均为FALSE。polyA与rRNA去除totalRNA均可作为普通bulk候选，但需要分别记录；单细胞、微阵列、smallRNA、靶向测序和不符需求的亚细胞分画不混入。

|序号|Kla PXD / 原始组名|RNA来源与选样|证据级别|限制|
|---|---|---|---|---|
|1|PXD033146 / pathological rotator cuff tendon|[GSE236746](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE236746)；GSM7574837;GSM7574838;GSM7574839|material_verified|非糖尿病撕裂肩袖肌腱；不是健康肌腱。排除糖尿病组；核对慢性病理和原蛋白样本背景。|
|2|PXD036307 / normal human lung|[GTEx_v8](https://gtexportal.org/home/downloads/adult-gtex/bulk_tissue_expression)；Lung; sample IDs pending|cohort_candidate|固定v8及样本属性；供体排除标准和最终n待核对，不能声称全部健康无创伤。|
|3|PXD046800 / hypertrophic scar|[GSE178411](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE178411)；hypertrophic scar; GSM list pending|material_verified|只取增生性瘢痕；排除普通伤口及瘢痕成熟阶段不符者；供体与治疗史待锁定；已过滤counts不能直接充当全基因分母。|
|4|PXD046800 / adjacent skin|[GSE178411](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE178411)；uninjured skin; GSM list pending|approximate|未损伤皮肤并不等于同患者瘢痕邻近皮肤；必须查取材部位/配对信息，未证实前只作外部背景。|
|5|PXD050470 / human hippocampus|[GTEx_v8](https://gtexportal.org/home/downloads/adult-gtex/bulk_tissue_expression)；Brain - Hippocampus; sample IDs pending|cohort_candidate|固定v8；核对年龄、死亡和神经疾病排除条件；不沿用197例声明。|
|6|PXD063047 / normal pregnancy placenta|[GSE114691](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE114691)；GSM3147325-GSM3147345; gestational-age filter pending|material_verified|21个control为初始候选，不等于21个足月；核对孕周、产式和胎盘取材层；排除PE/IUGR。|
|7|PXD064912 / human sperm|[GSE40181](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE40181)；GSM987942;GSM987943|material_verified|仅两个供体mRNA库；排除smallRNA GSM987944/45；处理表有RPKM阈值，建议原始reads重定量并检查精子纯度。|
|8|PXD066054 / BPH|[GSE132714](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132714)；GSM3890639-GSM3890656; clinical filter pending|material_verified|18个BPH候选；核对治疗和移行区信息。处理表基因名不能直接分析，优先raw reads以Ensembl重定量；不同建库需分层。|
|9|PXD075377 / adjacent liver|[TCGA-LIHC](https://portal.gdc.cancer.gov/projects/TCGA-LIHC)；Solid Tissue Normal; donor/file UUID pending|cohort_candidate|癌旁独立队列；不是GTEx健康肝；与HCC参考可配对需以caseID确认；去除重复aliquot。|
|10|PXD064038 / MEC and NEC ESCC groups|[TCGA-ESCA](https://portal.gdc.cancer.gov/projects/TCGA-ESCA)；Primary Tumor AND squamous histology; UUID pending|cohort_candidate|只保留ESCC而非EAC；MEC/NEC状态并未匹配；需实际临床治疗字段，不用01A推断未经治疗。|
|11|PXD066054 / prostate cancer|[TCGA-PRAD](https://portal.gdc.cancer.gov/projects/TCGA-PRAD)；Primary Tumor; UUID pending|cohort_candidate|核对组织学、术前治疗和独立病例；不预填499例。|
|12|PXD075377 / HCC|[TCGA-LIHC](https://portal.gdc.cancer.gov/projects/TCGA-LIHC)；Primary Tumor; UUID pending|cohort_candidate|核对HCC组织学及治疗；不预填371例；与本项目蛋白组不是同受试者。|
|13|PXD060185 / MCF7|[DepMap](https://depmap.org/portal/cell_line/ACH-000019)；ACH-000019|model_id_verified|MCF7身份核实；固定同一release、Model及OmicsProfiles后才能锁定表达行和是否log2(TPM+1)。|
|14|PXD028488 / HCT116|[DepMap](https://depmap.org/portal/cell_line/ACH-000971)；ACH-000971|model_id_verified|HCT116身份核实；与另外两组共享RNA参照，不构成三份独立转录组。|
|15|PXD053474 / HCT116|[DepMap](https://depmap.org/portal/cell_line/ACH-000971)；ACH-000971|model_id_verified|同一HCT116外部基线；全细胞RNA不能视为亚细胞分画匹配。|
|16|PXD066351 / HCT116 control and Roseburia co-culture|[DepMap](https://depmap.org/portal/cell_line/ACH-000971)；ACH-000971|model_id_verified|只代表外部基线；不能代表Roseburia共培养，蛋白组混合条件保持原状。|
|17|PXD028488 / TALL-104|[GSE163787](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE163787)；GSM4987488; SRX9725233|material_verified|真正bulk子系列，不是混合细胞scRNA伪bulk；仅一个TALL104库，不能估计组内生物变异；混合物种参考需只取人类Ensembl。|
|18|PXD050147 / HepG2 WT and SIRT1 or SIRT3 KO|[DepMap](https://depmap.org/portal/cell_line/ACH-000739)；ACH-000739|model_id_verified|修正原ACH-000657；外部亲本基线不匹配SIRT1/3 KO。|
|19|PXD054919 / A549|[DepMap](https://depmap.org/portal/cell_line/ACH-000681)；ACH-000681|model_id_verified|模型ID正确；原备用ENCODE指向H1应移除。|
|20|PXD060185 / MDA-MB-468|[DepMap](https://depmap.org/portal/cell_line/ACH-000849)；ACH-000849|model_id_verified|MDA-MB-468模型ID核实；固定release及表达profile。|
|21|PXD060185 / T-47D|[DepMap](https://depmap.org/portal/cell_line/ACH-000147)；ACH-000147|model_id_verified|修正原ACH-000934；亲本癌细胞并不意味着所有基因野生型。|
|22|PXD063266 / PC-3M|[GSE143451](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE143451)；GSM4259672;GSM4259673;GSM4259674|exact_match_gap|仅近似候选：PC-3M-1E8且scrambled载体对照，不是未经转染PC-3M。严格精确匹配仍缺；PC3不能顶替。|
|23|PXD070007 / glioblastoma stem cells|[GSE119834](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE119834)；44 GSC RNA libraries; GSM list pending|approximate|GSC群体独立模型面板；需核对G2907/G3028/G3264/GSC23/MES28/RKI；系列GSC23序号不能证明同一细胞模型。|
|24|PXD078013 / RKO WT and GSK3B KO|[DepMap](https://depmap.org/portal/cell_line/ACH-000943)；ACH-000943|model_id_verified|修正原ACH-000742；外部基线不能代表GSK3B KO。|
|25|PXD028488 / HEK293T|[GSE203529](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE203529)；GSM6175912; additional replicates pending|material_verified|该样本为HEK293T WT untreated；排除IFN及基因扰动；HEKTE或HEK293不可自动当293T。|
|26|PXD028737 / HMC3|[GSE275256](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE275256)；GSM8473662;GSM8473663;GSM8473664|material_verified|vehicle基线；排除SLO。是HMC3永生化模型而非原代小胶质。|
|27|PXD058534 / pretreated HK-2|[GSE269418](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE269418)；GSM8315283;GSM8315284;GSM8315285;GSM8315286;GSM8315297|material_verified|5个常氧溶剂24h对照候选；排除BAY/缺氧/TGF；不能代表原pretreatment条件；与下一组共享参照。|
|28|PXD078736 / HK-2 control and mannitol|[GSE269418](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE269418)；GSM8315283;GSM8315284;GSM8315285;GSM8315286;GSM8315297|material_verified|共享HK2基线；不是mannitol条件匹配，甘露醇也不自动等于无扰动。|
|29|PXD060185 / MCF10A|[GSE103520](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103520)；GSM2773021;GSM2773022;GSM2773023|material_verified|3个未处理MCF10A；排除其他细胞及转染。处理diff文件不是原始counts，优先获取完整表达或原始reads。|
|30|PXD070007 / neural stem cells|[GSE119834](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE119834)；GSM3384847-GSM3384855|approximate|9个NSC库；对应本项目两个具体NSC模型和发育年龄尚未证明，作为外部NSC面板。|
|31|PXD073311 / HUVEC control and Pg infection|[ENCSR000COZ](https://www.encodeproject.org/experiments/ENCSR000COZ/)；released whole-cell long polyA RNA; file IDs pending|material_verified|官方元数据确认HUVEC全细胞；继续锁定donor/replicates/file assembly；未感染基线不代表Pg感染条件。|

## 下游约束

1. 保留疾病和原蛋白组处理范围；非肿瘤不等于健康，亲本癌细胞不等于所有基因WT。外部基线只能叫参照，不能声称KO/感染/共培养条件匹配。
2. 三组HCT116及两组HK2分别共享RNA参照；用ReferenceKey显式追踪，不能重复计作独立转录组。GSC/NSC需额外核对具体模型身份。
3. 固定GTEx单一版本、GDC文件UUID与筛选快照、DepMap单一release及Model/OmicsProfiles；不能沿用旧清单的固定样本数或用Primary Tumor推断未治疗。
4. 冻结sample/donor/condition、文库、原始/处理文件、基因注释版本和SHA。以Ensembl/Entrez分析，Symbol仅展示；UniProt映射使用可追踪的BaseAccession规则。随机种子25。
5. 保留counts/TPM/FPKM/log转换的真实单位。DEG子表、已按表达过滤的矩阵及smallRNA不能定义全转录组分母。不同来源先独立QC和汇总；外部队列间直接比较表达强度需另行论证。
6. 本次未下载完整表达矩阵或执行跨组学分析；明确候选不等于承诺全部31组已经找到精确配对数据。严格PC3M仍是缺口，PC3和PC3M-1E8不自动替代。
