# 31组 RNA候选补充与普通全蛋白参照复审

日期：2026-09-14。审核基准：同一解剖组织和取材类别。癌旁肝不等于普通/远端正常肝，瘢痕邻近皮肤不等于任意正常皮肤，BPH不等于普通前列腺。条件、基因型和供体差异单独披露，不自动否决组织层面的匹配。细胞培养物核对细胞身份，不用整个器官替代培养细胞。

本次检查覆盖冻结注册表31组、18个普通全蛋白PXD、20个实际参考文件（全部存在并记录SHA256）、取数代码、文件表头与关键原始参数；结合公开数据库和论文核验。并未重新计算全部蛋白集合、重搜质谱或重建统计结果；没有修改冻结分析包。下表的“组织一致”不等于供体、处理、分画或统计独立性全部一致。

## RNA：找到的改进来源

|目标|新候选及选样|结论与限制|
|---|---|---|
|瘢痕邻近皮肤|GSE181540：GSM5505085、5505087、5505089（NS input）|论文明确采集瘢痕邻近全层正常皮肤，取材类别比旧GSE178411的uninjured skin更准确。|
|增生性瘢痕|同一GSE181540：GSM5505079、5505081、5505083（HS input）|可与邻近皮肤统一到同一研究；仍是外部RNA队列，不是与Kla同供体。|
|PC-3M|GSE235595：GSM7506018、GSM7506019、GSM7506020，DMSO对照|确认为PC-3M，解决旧PC-3M-1E8近似模型问题；保留DMSO培养4天的条件说明，排除KMI169及KMI169Ctrl。|
|GSC|GSE266884：优先GSM8255530、GSM8255531，MES28 shMCT1Control|与PXD070007/PXD069969同一论文，且MES28就是六个原模型之一。仅覆盖一个模型，不能声称六模型RNA齐全；另两份shKU70对照不自动合并成四份独立重复。|
|NSC|保留GSE119834的9个NSC RNA样本候选|GEO确认NSC培养物而非整块脑组织；符合NSC类别的外部背景。未找到ENSA/HMP1的精确RNA，继续保留具体模型未匹配的标记。|

DepMap 的9个细胞模型参照已在服务器固定到官方 `DepMap Public 24Q4`（2024-12-16）：蛋白编码基因 TPM 矩阵、`Model.csv` 和 `OmicsProfiles.csv` 均已下载并通过目录公布的 MD5。当前 `DepMap Public 26Q1` 目录对相应文件只返回元数据而不给可用 URL，因此 24Q4 是明确记录的可复现回退版本；这些外部亲本表达仍不能代表各自的 KO、感染或共培养条件，且分析前必须用 ACH 模型 ID 和稳定基因 ID 完成行列及尺度核对。

GSE181540虽然GEO统一标为RIP-Seq，但NS input样本明确注明抗体为none，采用rRNA去除total RNA、未做免疫富集；只有input臂可作为普通bulk候选。IP臂、甲基化峰表不能纳入。论文有Gallus gallus参考基因组文字，而GEO写人类hg19，存在来源内矛盾：取材与input建库证据支持候选资格，真正分析前应核验reads与参考基因组，不能直接信任全部处理结果。PC-3M处理表为RefSeq NM转录本主键，需稳定版本映射至Ensembl/Entrez；不能走Symbol合并。

证据：[皮肤论文取材与方法](https://doi.org/10.3389/fcell.2021.748703)、[皮肤input样本](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5505085)、[PC-3M系列](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE235595)、[GSC系列](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE266884)、[同源蛋白/RNA论文](https://doi.org/10.1038/s41556-025-01839-y)、[NSC系列](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE119834)、[DepMap公开目录接口说明](https://forum.depmap.org/t/provide-an-open-endpoint-for-latest-version-retrieval/4652)。

新版31组候选见`rnaseq_reference_candidate_31.csv`，保留原GroupID、PriorSource和全部AnalysisReady=FALSE。替换候选不表示完成矩阵验收。旧20260912审计注册表保留作为历史快照。

## 普通全蛋白：组织取材与实际选择

|组别|实际普通全蛋白来源与选择|复审意见|
|---|---|---|
|肩袖病理肌腱|PXD033146，HA119TQ普通蛋白文件|与病理肩袖肌腱一致。|
|正常肺、正常胎盘|PXD010154，分别独立lung/placenta文件|器官一致；未锁定更细的肺亚区或胎盘层次。|
|瘢痕、邻近皮肤|PXD046800，QB002普通蛋白，分别HSP1–4、NSP1–4|官方明确4个瘢痕和4个邻近正常皮肤；代码按HSP/NSP分列，无普通皮肤替代。|
|海马|PXD050470，原研究tables4|当前是同研究海马来源；不是旧建议配置中的PXD043880。|
|精子|PXD066517，20240275.tsv|材料是精子；表内同时有L/N样本，L/N供体表型未独立核实，不能统称健康供体。|
|BPH、前列腺癌|PXD066054，DA/Protein_Quant.tsv，分别NAT和PCa列|官方明确5 BPH和5 PCa。此项目的NAT命名对应BPH，不应仅凭缩写解释成癌旁。|
|邻近肝、HCC|PXD065775，分别ANTs和CISs工作表|官方区分ANTs/PTTs/DNTs；本地表为ANTs/CISs/DNTs。代码读取ANTs用于邻近肝、CISs用于HCC，没有误选DNTs。|
|ESCC|PXD065830，Dataset1的2.a普通蛋白表，94个T列|使用ESCC肿瘤，排除24个N列；不是健康食管或磷酸化表。|
|MCF7、两组HCT116、A549、MDA-MB-468、T47D、RKO|PXD030304，指定SIDM/细胞行|实际行核对通过；HCT-116与HCT116、T47D与T-47D是显式别名，不是错模型。|
|HCT116共培养|PXD066351，DA普通蛋白文件|细胞模型一致。|
|TALL-104|PXD028488，TALL-Nonenrichment/proteins.csv|当前使用TALL-104非富集臂；没有调用旧TALL-1/Jurkat备选。|
|HepG2|PXD050147，SIRT_proteinGroups.txt|9列pro普通蛋白重复，WT/SIRT1KO/SIRT3KO各3列；论文以蛋白量归一化修饰比值支持普通蛋白臂。不能仅因该研究含Kla/Kac就判整项是富集结果。原始run到pro文件的全链未重新重建。|
|PC-3M|PXD022005，txt_proteomics.zip，Intensity H|论文证实重标为PC-3M、轻标为PC-3；当前选重通道，未用PC-3替代。普通蛋白参数为Protein FDR=0.01。|
|GSC、NSC|PXD069969，SA206LQB1，六GSC与ENSA/HMP1列|与Kla六GSC/两NSC模型对应，表头支持分组；没有用胶质瘤整组织或普通脑替代。|
|HEK293T|PXD030304，Control_HEK293T_lys|身份一致；属于401条过程质控裂解物run，不能视为401个生物重复。|
|HMC3|PXD028737，普通proteinGroups，H0/H24|模型一致。|
|两组HK-2|PXD072220，前三个Log2Quantity列|身份是HK-2；内部列amostra1/3/4与公开Control_1/3/4的显式重命名映射尚未找到，暂保留样本别名证据缺口。|
|MCF10A|PXD002400，evidence.txt中10A前缀|身份一致；来源参数为Protein FDR=1、PSM FDR=0.01，按本项目要求保留来源处理方案并如实披露。|

## 需要处理的问题

1. **来源参数披露：MCF10A。** 直接读取本地msms.zip/parameters.txt得到`Protein FDR=1`、`PSM FDR=0.01`，对应论文也说明关闭蛋白FDR。按老师确认的项目规则，不对Kla或普通全蛋白另加概率/蛋白FDR阈值，也不因此排除或重处理该来源；当前`load_pxd002400()`保留既有10A raw前缀、非reverse/contaminant和合法UniProt提取规则。此参数必须随结果披露，且不得写成蛋白层1% FDR。证据：[原论文](https://doi.org/10.1186/s13059-015-0742-x)及本目录`PXD002400_parameter_evidence.txt`。
2. **中优先级：HK-2样本别名证据。** 现有脚本按前三列而非明确样本名选择，实际是amostra1/3/4。与公开Control_1/3/4编号相符，但没有显式rename manifest，不能把推断写成逐run证明。现阶段属于证据不足，不是已经证实选错。
3. **用途限制：技术重复/汇总参照。** ProCan细胞系为技术重复/均值，HEK293T为过程控制；PC-3M有合并SILAC通道。它们可支持检测背景，但不能因此成为独立生物重复。sample-only包已排除aggregate，这次未改变其策略。
4. **标签与配置漂移。** 旧`reference_proteome_selection.json`仍推荐海马PXD043880，当前31组是PXD050470；代码也残留非当前使用的TALL-1/Jurkat和PC-3备选分支。重建必须以冻结31组注册表和实际选择器为准。HK-2旧SampleGroup仍写mannitol，而当前显示标签/命名依据已写缺氧-复氧；不能用旧名字反推实验条件。本次不改稳定组键。
5. **普通全蛋白不应称“无乳酸化蛋白组”。** 它是未经Kla富集的全蛋白测量，并不排除其中蛋白发生乳酸化。

关键一手证据：[ANTs/PTTs/DNTs](https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD065775)、[瘢痕/邻近皮肤](https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD046800)、[BPH/PCa](https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD066054)、[ProCan技术重复](https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD030304)、[PC-3M重标](https://pmc.ncbi.nlm.nih.gov/articles/PMC10372916/)、[HepG2普通蛋白归一化](https://pmc.ncbi.nlm.nih.gov/articles/PMC11440250/)。

## 复现与限制

从当前工作树执行`Rscript --vanilla workflow/audit_reference_materials_20260914.R /Users/gzy2520/Desktop/Research/kla`。脚本读取原工作区源文件但仅向本目录输出，使用R和seed 25；不使用Symbol进行蛋白或基因集合分析。`reference_source_files.csv`记录20文件SHA，`source_headers.csv`记录表头/压缩包成员，`targeted_checks.csv`记录检查与未定项，`proteome_material_review_31.csv`逐组记录结论。

初次检查中HCT116/MDA-MB-468别名和HK-2多行表头造成的检查脚本假阴性，已依照实际源格式纠正；它们不是生产流程选错数据的证据。当前所有布尔结构检查通过，NA行明确表示仅描述或证据待补。此通过不消除HK-2别名证据缺口；MCF10A的FDR设置作为保留来源方案的参数披露。
