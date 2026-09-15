# 老师意见处理与重复数复核（2026-09-15）

本轮优先于 2026-09-14 来源清单。老师原始截图保存在 `teacher_markup.png`；逐组结果为 `teacher_review_31.csv`。不把老师颜色改写为我们的判断：绿色13行、黄色12行、红色2行、未标4行。绿色代表老师认可来源，仍需独立完成重复数、样本资格和数据QC；未标代表老师未打开网站，不代表否决。

## 生效规则

1. 组织以独立供体、细胞以独立培养/实验重复计数。GSM数量只作初筛；SRR、lane、R1/R2、重复定量文件不增加生物学n。供体pool按pool统计，不能把pool内人数当可估计的独立样本。
2. 尽量获得同一来源、同一材料和同一条件的至少3个独立样本。n<3暂列描述性备选，不进入本项目主要带误差线的比较。n=2数学上能计算SD，但估计不稳；n=1不能估计组内变异。不能通过复制数据、对基因/reads bootstrap或混合不同处理组来凑n。
3. 优先同一研究的完整表达矩阵；来源内先QC，保留Ensembl/Entrez。跨研究样本不能直接当同批重复；确需跨研究扩充须保留研究/批次层级。
4. 分开记录组织身份、处理条件、重复数和矩阵可用性。绿色、n>=3或文件已下载，均不单独构成AnalysisReady。

## 红色 HK-2：旧参照停用于主分析，锁定新候选

旧 GSE269418 所选5例是常氧溶剂对照，但方案包含过夜血清饥饿、0.0001% DMSO、4 μM HCl及0.0001% BSA。不是简单的未处理HK-2；也没有证据显示所选对照本身接受BAY/TGF或缺氧。旧清单仅写“常氧溶剂24h”遗漏了重要条件。KLA31_27/28本轮均停止沿用该来源作为主参照。

新候选 **GSE240748：GSM7708662、GSM7708663、GSM7708664**。三个样本逐项确认：人HK2、WT、Control、24h、rep1–3，处理说明为未处理或TGF-β1处理，本次只取未处理臂；培养为DMEM/F12+10%FBS。每例提供 `NC1/2/3.gene.txt.gz`，说明为基因counts/表达表。独立培养来源、实际表头/ID/完整性及原蛋白处理范围仍要验收，故暂不提升AnalysisReady，也不能把此未处理基线称为mannitol或其他预处理条件的匹配实验。

证据：[旧对照](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8315283)、[新来源](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE240748)、[rep1](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM7708662)、[rep2](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM7708663)、[rep3](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM7708664)。原始SOFT已保存在evidence。

文件复核补充：3个gzip文件已下载并通过完整性检查，每份42,112行，提供Counts/TPM/FPKM，三份行标一致且无重复。每份这三个数值列合计72个缺失值、无负值，缺失不得自动当0。但`Gene_id`形如`gene-MIR6859-1`，并非Ensembl/Entrez稳定ID。不能仅去掉`gene-`后按Symbol连接；需通过确切注释版本/坐标核对建立稳定ID映射。来源内元数据较旧HK-2合适，尚不等于可直接纳入分析。详见`hk2_replacement_file_qc.csv`。

## 黄色条件复核

- **瘢痕/邻近皮肤（03/04）**：所选6个GSM为input，样本说明`rip antibody: none`，提供单独的去rRNA input建库/表达分析流程。因此研究m6A不等于这些input受到甲基化药物处理，更不能把IP峰表当表达矩阵。老师黄色保留，另记“input属性有证据”；原论文物种/参考注释矛盾仍需核查。每臂名义3例，需锁配对供体。[GSM5505079](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5505079)
- **DepMap（13–16、18–21、24，共9行/7个模型）**：每个模型当前只有一行汇总表达，不能据此估计该细胞系重复间误差。全部降为外部描述性背景，继续找相同模型的>=3独立RNA样本。HCT116三组共享同一参照；HepG2/RKO亲本不代表KO，HCT116亲本不代表共培养。老师对加药的疑问尚不能由ACH身份本身消除。
- **TALL-104（17）**：原GSE163787仅一个bulk库，达不到目标。新检索到GSM6934151为DMSO 24h的第3个生物学重复，但平台是Human Transcriptome Array 2.0芯片，不能悄悄替换bulk RNA-seq。当前仍保留缺口。[样本](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM6934151)

## 即使标绿也需补重复的组

|组别|当前可观察单位|本轮处理|
|---|---|---|
|精子07|2供体|寻找>=3个独立供体；GSE65683是后续线索，需锁生育状态和完整基因定量，不能直接用RNA element百分位矩阵|
|MES28 23|2个shMCT1对照库|继续找同模型同条件；不能把shKU70对照直接拼成4例|
|HUVEC31|2个ENCODE文件|供体/重复关系仍待锁定；无论如何现有文件不足以证明n>=3|

连同9行DepMap与TALL-104，共**13/31行当前所选数据不足以支持n>=3**。这不是13种独立模型。另有6行大队列尚未锁最终样本；NSC的9个库不是已经核实的9个独立生物学重复。

## 新检索候选及排除项

|目标|来源|本轮判断|
|---|---|---|
|HCT116/RKO|GSE81194|论文/series声明各3个培养重复，但GEO每细胞系仅一个GSM；必须追溯SRA BioSample/文库层级和是否混池，不能直接宣称补齐。附表主要是变异鉴定结果，不是完整表达矩阵。|
|MDA-MB-468|GSE157383，GSM4763919为DMSO第3例|存在逐基因readcounts线索；7天DMSO须记录，其他两个GSM和独立性尚待逐项核查。保持条件候选。|
|T47D/MCF7|GSE310291|有RNA-seq DMSO对照的时间点设计；需锁同一时间点三个重复、激素剥夺等培养条件，排除CUT&Tag/药物臂；尚未通过。|
|A549|GSE171750|series声明triplicate；需核对所选untreated臂有无siRNA，暂作待核查线索。|
|MCF7|GSE149241/GSM4494786|虽然标题Control/untreated，但growth protocol涉及移植瘤及雌激素补充；不采用作普通培养MCF7替代。|
|RKO/A549|GSE145108|每条件只有两个独立重复，不解决n>=3问题。|
|HK-2|GSE193192|shCtrl/AGEs设计，存在转染条件，不优先替代。|

GSE81194：[GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE81194)、[原论文](https://pmc.ncbi.nlm.nih.gov/articles/PMC5305277/)。其他来源原始记录见evidence及`new_source_evidence.csv`。本轮只是上述候选的证据等级，不声称所有缺口已补齐。

## 当前交付和未完成项

已交付31行老师意见/重复数审核和新的HK-2三个候选样本；已下载的DepMap与ENCODE QC结果继续作为技术整理产物，不能凭QC通过升级为主分析样本。其他未解决来源继续保留明确缺口。四个未标色GDC行继续按癌旁/肿瘤/组织学、供体去重及临床条件审核，不能标为老师已通过。
