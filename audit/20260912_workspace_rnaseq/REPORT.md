# Kla 工作区与 RNA-seq 参照审计（2026-09-12）

本次目标是版本/来源审计与RNA候选筛选；不重新分析蛋白组、不改已批准图、不把外部RNA称为同样本配对组学。检查覆盖Git状态、默认运行器、候选统计入口、输入注册表、验证声明、两份Word正文、renew中24个交付文件及原RNA清单全部37个唯一编号。没有逐一重新计算45 GB源蛋白数据。

## 1. 当前分支不能直接视为干净的最新分析基线

- `add_data` @ `aa195c6`。开始审计时已跟踪文件无修改，但两份Word未跟踪：`manuscript/DDR_Kla_manuscript_V3.docx`、`DDR_Kla_manuscript_V3_backup.docx`。因此“已跟踪文件干净”不等于整个工作树干净。
- `add_data`没有upstream、远端同名分支不存在。live `origin/main`为`bc02d1d`；两边分叉：本分支独有1提交，远端main独有10提交，包括`c349c7f`（31组sample-only修正）及`92d37c2`（强度审计）。本地main仍在`037645d`。
- 不应在本分支直接继续生成“最终图”。后续整合应从核实后的远端main创建隔离工作树，移入本次审计资料，再逐项核对31组来源；不要用旧分支整体覆盖新修正。本次未merge/reset/checkout，也未提交或push，未触碰两份Word。

## 2. 已确认的31组合同

`groups_31.csv`来自本地`data/candidate/escc_inclusion_20260903_pxd065830_tumor_reference/publication_input/group_summary_30.csv`。源文件名仍含30，但实有31行；这是一个误导性的历史文件名，不是删除ESCC的理由。

- 9非肿瘤组织、3肿瘤组织、12癌细胞、7非癌细胞/培养模型。
- 22个不同Kla PXD、18个不同全蛋白参考PXD；31表示匹配分析组，不是31个独立队列或同样本配对研究。
- ESCC：Kla `PXD064038`，全蛋白参考`PXD065830`。
- 原RowOrder跳过13、26、27，保留不改；另加`KLA31_01`—`KLA31_31`供本次注册表关联。组键为PXD+完整SampleGroup，不能用PXD单独join。
- 蛋白组继续使用BaseAccession；RNA使用带版本来源的Ensembl/Entrez。跨组学映射应保存版本、下载日期、SHA、1对多处理规则及未匹配清单。当前未找到支持原清单“配置里已经冻结BioMart映射”的证据。

## 3. 会干扰下一轮的旧文件与声明

| 位置 | 查到的问题 | 本轮处理/后续约束 |
|---|---|---|
| README、AI_HANDOFF_CURRENT.md | 仍声称30组/9-2-12-7、399；handoff还写旧分支和旧H3C1缺失结论 | 增加显著历史警示和本审计入口；历史内容保留供追踪 |
| workflow/run_pipeline.R | 默认硬编码30组输入；覆盖输入环境设置；publication/source模式先清空results/figures与results/supplementary | 未运行；仅可作历史30组重现，不得当31组入口 |
| tests/validate_publication_contract.R及source准备脚本 | 硬性断言30组、399蛋白等 | 通过这些检查不能证明当前31组正确；本轮未机械改成31掩盖输入问题 |
| R/candidate/build_figure1_category_boxplot.R | 读取输入后没有ObservationType==sample过滤；默认旧candidate | 当前310条输入中只有272 sample；其余38为聚合/条件/模型/池等 |
| R/candidate/build_ddr_pathway_summary_boxplots.R | 同样缺sample过滤 | 当前686条通路记录中504是sample，剩182不是；7通路记录不是独立样本 |
| R/candidate/boxplot_significance.R | Pro/Inh同源样本构成配对，但普通双因素aov未建模此依赖；另有参考源复用 | 后续必须采用配对/重复测量或适当混合模型；本轮不擅自重算 |
| manuscript主文与backup | 主文Methods为31；backup Methods仍30但Results31。主文“29”指全蛋白展示行，不能写成29个独立PXD | 两份未跟踪文档保留，backup不得作为最新Methods来源；主文自动验证等声明尚需与最新代码对齐 |
| .gitignore | 31组扩展输入目录、results等未纳入追踪 | 保存本轮31组合同及关键SHA；Git干净不能证明这些实际输入/图一致 |
| PXD010154本地fullproteome_RNAseq_txt.zip | 检查胎盘压缩包为MaxQuant蛋白结果，未发现RNA表达矩阵 | 文件名含RNAseq不能作为已获得RNA的证据 |

`observation_unit_counts.csv`和`proteome_reference_reuse.csv`给出具体分母/复用。新审计脚本只读历史输入，不写results。

## 4. renew图不是已验证的同一套可复现发布

对`/Users/gzy2520/Desktop/renew/kla`的18 PNG和6 XLSX，按同文件名搜索当前results及data/candidate，再比SHA-256：6张PNG和1份XLSX找到完全相同文件，其余未找到同名同SHA副本。详见`renew_file_hash_matches.csv`。这不证明其余文件错误，也未证明像素不同（PNG元数据也可能改变SHA）；只是不能把整个renew目录自动绑定到当前某个renderer/commit。

后续视觉风格以用户给定renew为参考，但应先逐图绑定来源脚本、输入、提交、参数和SHA。本轮没有重绘或覆盖参考图。

## 5. Gemini清单结论

桌面文件与仓库原`docs/MATCHED_RNASEQ_REFERENCE_CATALOG_31.md`逐字节相同。原版存入`Gemini_catalog_original.md`；全部37个唯一GSE/ENCSR/ACH编号逐项见`original_accession_audit.csv`，不再将原版代码示例当生产流程。

典型严重错误：

- [GSE160115](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE160115)是拟南芥，不能补ESCC；[GSE103939](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103939)是水蚤，不能补人精子。
- [GSE75010](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75010)是胎盘芯片；[GSE114724](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE114724)是乳腺肿瘤免疫单细胞，均不满足对应bulk需求。
- ENCODE官方JSON表明[ENCSR000EYO](https://www.encodeproject.org/experiments/ENCSR000EYO/)为K562，[ENCSR000EYP](https://www.encodeproject.org/experiments/ENCSR000EYP/)为H1；[ENCSR000COV](https://www.encodeproject.org/experiments/ENCSR000COV/)为H1胞质RNA，不能同时补A549与MCF10A。
- HepG2应查ACH-000739、T47D为ACH-000147、RKO为ACH-000943；原编号不能照用。PC3纠正为ACH-000090后仍不是PC3M。
- GSE119834可保留为GSC/NSC bulk候选，但是否为本项目具体模型仍未证明。其余编号的保留/排除/未解决状态见逐项表，未把404解释成确定不存在。

原清单还把疾病组织替换成健康对照，把癌细胞parental称为严格WT，把不同GTEx版本样本数混用，并将recount3 RPKM示例命名为TPM。这些是设计问题，不只是错别字。用户本轮只要求bulk RNA补充，没有要求把所有疾病/KO/处理组改成健康未处理组；本审计保留原蛋白组范围，并明确哪些RNA只是外部基线。

## 6. 新参照清单如何用

`rnaseq_reference_registry_31.csv`覆盖31组，包含源编号、具体样本或待筛选规则、证据级别、建库/单位、复用键、限制和一手URL。可读版在仓库docs同名清单。EvidenceStatus不等于AnalysisReady；**31行当前均为AnalysisReady=FALSE**，因为本轮完成来源/身份核查，没有完成所有表达矩阵下载、样本筛选与基因ID/单位QC。

当前分级为12组材料/样本已核实、9组模型ID已核实、6组队列待筛、3组近似匹配、1组精确缺口。

较明确的新来源包括肌腱GSE236746、胎盘GSE114691、精子GSE40181、BPH GSE132714、TALL104 GSE163787、HMC3 GSE275256、HK2 GSE269418、MCF10A GSE103520、HUVEC ENCSR000COZ。TALL104的GSM4987488是实际bulk单库，不能用同SuperSeries的混合单细胞库替代；一个库也不能估计生物学重复方差。[原始样本记录](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4987488)

精确缺口是PC3M（目前只有PC3M-1E8+scrambled载体近似候选）；瘢痕邻近皮肤、GSC和NSC具体模型也只达到近似匹配。HCT116共培养、HepG2/RKO KO、HK2处理及HUVEC感染均不能由亲本基线RNA证明条件匹配。

下一阶段每条可纳入来源应冻结：sample/donor/condition、raw run或processed file ID、建库及分画、版本及SHA、全基因矩阵单位和主键、过滤造成的基因缺失、参考复用关系。先在来源内汇总，不能把跨平台TPM/FPKM/logTPM/counts直接拼接做强度比较，也不能把未检测/未报告当0。禁止仅用DDR/DEG筛选后的表定义全转录组分母。

## 7. 验证与边界

从仓库根运行：`Rscript --vanilla audit/20260912_workspace_rnaseq/scripts/build_audit.R`。它验证31组唯一键/分类/ESCC、RNA注册表逐行匹配，导出观察单位和来源复用表；不联网、不重画、不验证RNA表达值。

sources保留可访问的一手ENCODE JSON、ENA run metadata、论文XML以及检索证据快照。部分NCBI下载端点SSL失败、DepMap下载页为CAPTCHA；`depmap_download.html`是失败页面，不是数据。未下载到的矩阵绝不声称已经验收。检索索引快照只用于可追踪筛选，正式取数仍须冻结真实文件。

补充取数证据：已从ENA获取7个BioProject的run元数据（含已有肌腱项目），新增6个项目合计247条run记录；这些包含未纳入条件，并非247个合格样本。TALL104核实到SRR13296314。FASTQ地址/MD5在sources的TSV中，尚未下载FASTQ。工作区99个已跟踪文件及关键本地输入的SHA见workspace_file_manifest.csv（其中跟踪文件为审计完成时版本）。
