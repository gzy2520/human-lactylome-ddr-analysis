# 首次材料层级外部验证：食管鳞癌与癌旁（2026-09-28）

已实际下载并处理新蛋白组和新 RNA，使用 2026-09-27 冻结的 DDR / DNA repair 随机森林预测，没有重训、调参或按外部结果挑选模型。结果不支持当前模型跨数据来源准确预测占比。它是组织材料层级外部验证，不是同一患者多组学配对验证。

## 数据对应

|用途|来源|实际使用|
|---|---|---|
|普通蛋白组与乳酸化组|PXD063945，PMC13195773 补充材料 mmc3 / mmc4|分别取 H/L 组的 T（肿瘤）与 N（癌旁），排除所有 NAC（新辅助化疗）列|
|外部 RNA|GSE130078，韩国 YSH 队列|23 条肿瘤 + 23 条癌旁；去 rRNA total RNA bulk RNA-seq；全部为新的 GSM ID|
|旧 RNA 对照|既有 TCGA_ESCA_ESCC|95 条既有 RNA，配新蛋白标签；单列为对照，不称独立 RNA 测试|
|模型|outputs/20260927_ddr_fraction_models|DDR 与 DNA_repair 的 RF_mtry135.rds 原文件，校验和不变|

[质谱论文](https://pmc.ncbi.nlm.nih.gov/articles/PMC13195773/)对应 Huashan Hospital 临床材料；[旧 ESCC 项目记录](https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD064038)来自 West China Hospital、不同提交者。新项目没有出现在旧训练来源映射中。没有把不同 PXD 自动解释为逐患者层面已经证明完全不重叠；本次依据可核实的项目与机构来源作材料级来源区分。

[GSE130078](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE130078) 是另一个研究来源的患者组织，不要求与质谱同人配对。论文另有细胞 HERES 敲低实验，但本次未使用那些细胞数据。临床补充表没有明确治疗前史，因此 RNA 的“治疗前、未经任何临床治疗”状态未得到确认；这项外部验证保留此限制。目录中的 ESCC_untreated 表示用于定义标签的未治疗质谱组，不代表已确认 RNA 的治疗前史。

## 蛋白计数

使用既有 GO 集合与稳定 UniProt BaseAccession，取正强度检出，分组做蛋白并集：

`占比 = 100 × |目标集合 ∩ 普通蛋白组 ∩ Kla组| / |目标集合 ∩ 普通蛋白组|`

未治疗肿瘤普通蛋白列 25 个，Kla 列 22 个；癌旁也为 25 / 22。原补充表两种检测的编号不完全一致，未臆造患者对应表。按用户确认的材料级口径分别汇总。H/L 分化组分别计数作为补充，主结果合并未经治疗的肿瘤。不能将其解释为25人逐个实测的占比。

原表共有 9260 个蛋白行和1975个位点行；筛选至少一个样本正强度后分别为8281和1836，与论文报告一致。没有把全缺失行当成检出，也没有把多个 Kla 位点重复计为多个蛋白。

|材料|DDR 分子/分母|DDR 占比|DNA repair 分子/分母|DNA repair 占比|
|---|---:|---:|---:|---:|
|未经 NACT 的肿瘤组|47/399|11.7794%|33/233|14.1631%|
|相应癌旁组|41/366|11.2022%|27/209|12.9187%|

癌旁有1个目标Kla蛋白未在普通蛋白组中检出，保留于审计计数，不纳入分子。详见 external_labels.csv、protein_members.csv、selected_columns.csv。

## RNA 处理

作者提供每样本 57783 个 Ensembl ID 的 featureCounts 计数。按作者方法使用 GRCh37 / GENCODE v19，以每基因合并外显子区间长度（不重复累加重叠外显子）计算 rate=count/length，再按全57783基因归一到 TPM，最后 log2(TPM+1)，提取冻结模型的18332基因并固定顺序。

46条记录均满足：18332模型基因全部存在；57783个计数基因全部有长度；TPM列和为1000000；无补零、缺失值插补、symbol转ID或测试集联合归一化。样本ID与1898条原始训练RNA零重叠。这是使用作者计数的处理后验证，未从FASTQ重比对；不同研究的定量流程仍可能产生差异。

## 结果

|RNA与材料|目标|RNA条数|新蛋白参考占比|平均预测占比|MAE（百分点）|±5个百分点内|
|---|---|---:|---:|---:|---:|---:|
|新RNA，肿瘤|DDR|23|11.78%|20.85%|9.07|0%|
|新RNA，肿瘤|DNA repair|23|14.16%|23.98%|9.81|0%|
|新RNA，癌旁|DDR|23|11.20%|21.31%|10.11|0%|
|新RNA，癌旁|DNA repair|23|12.92%|24.78%|11.86|0%|
|旧TCGA RNA，新肿瘤标签|DDR|95|11.78%|21.66%|9.88|0%|
|旧TCGA RNA，新肿瘤标签|DNA repair|95|14.16%|26.59%|12.42|0%|

只有两种材料共享标签，且来自同一个蛋白项目，因此没有用46条记录包装成46份独立蛋白真值，也没有报告无意义的单材料R²或把两个材料的R²当泛化证据。肿瘤是与既有训练材料相同类别的新来源验证；癌旁为附加探索，不能与原本健康食管或其他器官正常组织混用。

旧训练 ESCC 材料标签为 DDR 21.19%、DNA repair 26.12%；新蛋白参考为11.78%、14.16%。新RNA的预测仍较接近原训练标签。这与模型依赖已有材料标签的解释相容，但不能仅据此断定它只识别组织，也不能把差异全部归因于批次。真实生物差异、取材、检测深度及分组并集口径都可能参与。

本结果支持保留“已有材料体系内重现较好”的原结论，同时明确：当前版本对这个新来源没有达到准确预测。没有用这次外部标签做常数校准再宣称外部验证成功。

## 验收与复现

- 主要图：outputs/20260928_external_escc_evaluation/newRNA/figures/predictions.png（及PDF）。红菱形是新蛋白组材料占比，蓝点是逐RNA预测。
- 旧RNA对照：同目录下 reused_TCGA/。
- 逐样本预测与汇总：各子目录 predictions.csv / summary.csv。
- 质谱标签与成员：outputs/20260928_external_escc_labels/。
- RNA矩阵、元数据、QC、长度表：outputs/20260928_external_escc_rna/。
- validation.txt：标签复算、NAC排除、样本ID、TPM、冻结模型hash、预测重载与指标复算通过。
- audit/20260927_external_material/retrieval_manifest.csv：成功下载文件的URL/SHA256。大型原始下载、完整论文附件和最初失败文件仅保存在本地，不混入结果目录或Git；关键蛋白表、处理后RNA和可复现脚本入库。

命令从仓库根运行；输出必须为新目录。原下载目录保留，fetch.sh可重新下载到新目录并核验归档。首次下载的PMC13195773_supp.zip不完整；本次仅使用通过ZIP校验的PMC13195773_supp_complete.zip。NCBI FTP的GSE130078_RAW.tar下载超时；使用GEO HTTPS下载入口完整取回的GSE130078_download.tar（46成员，逐成员gzip通过）。

```sh
bash workflow/external_material/fetch.sh NEW_DOWNLOAD_DIRECTORY
# 使用审计manifest核对下载，并解压RNA与两篇论文补充文件到脚本指定目录。
Rscript --vanilla workflow/external_material/escc_labels.R NEW_LABEL_OUTPUT
Rscript --vanilla workflow/external_material/prepare_rna.R NEW_RNA_OUTPUT
# predict.R固定读取本次已冻结的外部标签路径。要更换标签须显式改路径并留存输入hash。
Rscript --vanilla workflow/external_material/predict.R outputs/20260928_external_escc_rna/full_rna.rds NEW_PREDICTION_OUTPUT GSE130078_new_RNA
Rscript --vanilla workflow/external_material/validate.R
```
