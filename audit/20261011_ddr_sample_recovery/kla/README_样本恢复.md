# Kla 样本检出蛋白恢复审计

2026-10-11。原始数据来自 `/Users/gzy2520/Desktop/Research/kla/data`。本目录仅恢复缺失位置的乳酸化组样本信息，未修改共同作图输入、图件或冻结注释。

## 恢复结果

恢复 11 个材料组的 26 个观测单位，保留原组全部实验条件。具体数值见 `recovered_kla_sample_values.csv`，逐样本稳定 UniProt ID 见 `recovered_kla_membership.csv.gz`。

| 材料 | 恢复单位 | 原表证据与误差线依据 |
|---|---:|---|
| adjacent liver / HCC，PXD075377 | 每种组织各 1 个样本 | 元数据明确来自同一患者的一对组织；原表为 Control 和 HCC 两列。旧脚本的 pool 标签需要修正。每组 n=1。 |
| HCT116，PXD028488 | 1 个样本 | PSM 中一个 enrichment 原始文件；从 PSM 与蛋白肽段映射重建。 |
| HEK293T，PXD028488 | 1 个技术合并单位 | HCD24/28/32/36/40 是碰撞能量，原始文件含 mix。各能量数据合并，旧设计文件的 timepoints 描述错误。 |
| HCT116，PXD053474 | 1 个技术合并单位 | 重新解析 DDA/DIA 和作者 S3，合并亚细胞组分、分级及采集模式；保留 S3 位点或 DDA/DIA 双重支持位点。590 蛋白，65 DDR。 |
| pretreated HK-2，PXD058534 | 1 个样本 | summary.txt 只有一个 HK2 原始文件。 |
| HCT116 control / Roseburia，PXD066351 | 每条件各 1 个样本 | MSstats 表 R.FileName 仅 NC116 与 R116，R.Replicate=0，R.Fraction=NA；DIA_result 同样只有两次采集。 |
| NSC / GBM stem，PXD070007 | 2 / 6 个模型 | 原表有 ENSA、HMP1 及六个 GBM 模型列；没有模型内样本或原始文件列。可呈现模型间变异；UnitType=model。 |
| HUVEC，PXD073311 | 6 个样本 | A0h_1/2/3、A6h_1/2/3 六个 raw 强度列，分别为两个条件各三次样本采集。 |
| RKO，PXD078013 | 4 个样本 | sg_1/sg_2 与 sg737_1/sg737_2；元数据明确 WT/KO 每组两个重复。SDRF 生物学重复字段空；样本命名、原始文件与元数据共同支持分组。 |

HUVEC 和 RKO 旧脚本将样本分布归入 condition，导致后续只接收 sample 的图丢失这些位置。模型数据仍保留 model 分类。所有条件名写入 Condition；没有只选择 control 或 WT。

## 计算与验证

每个单位先合并检出位点对应蛋白，并按去除 isoform 后缀的 UniProt accession 去重。乳酸化组 DDR 比例为 `DDR 检出蛋白数 / 全部 Kla 检出蛋白数`。复用历史解析器的污染物、反向库、位点与蛋白组代表选择规则；PXD053474 复用其作者表和双采集模式联合支持规则。

DDR 使用冻结 `data/publication_input/human_ddr_go_annotations.tsv`，筛选条件与 `prepare_escc_inclusion_inputs.R` 一致：TAXON ID 为 9606 或缺失，QUALIFIER 不含 NOT，GENE PRODUCT DB 为 UniProtKB。没有将原始数据中缺失的 membership 注释默认为非 DDR。

验证结果：所有恢复 accession 都存在于对应组冻结 membership，DDR 标记差异为 0；无重复的 PXD × SampleGroup × SampleID × BaseAccession；各单位分母大于零，DDR 分子不大于分母。差异表仅有表头。

`condition_summary.csv` 提供按条件的均值、SD、SEM。HUVEC/RKO 的条件内样本组计算 SEM；n=1 与模型单位保留 SEM 空值。模型组同时提供 SD，供明确标注模型间差异时使用。

## 文件与复现

- `source_ledger.csv`：每项恢复结论、列名、无法继续拆分的具体原因。
- `source_file_inventory.csv`：这些 PXD 的全部本地文件清单，包括可替代处理结果和压缩包。
- `source_headers.json`：核心原表列名。
- `pxd028488_source_sites.csv`、`pxd053474_source_sites.csv`：从原表重建的位点证据。
- `pxd053474_site_reconciliation.csv`：作者表与 DDA/DIA 位点支持比对。
- `membership_annotation_differences.csv`：恢复 accession 与冻结 membership 的差异，目前为空。

先在本工作区根目录运行 `extract_legacy_sources.py`，再运行 `Rscript audit/20261011_ddr_sample_recovery/kla/recover_kla_samples.R`，最后运行 `audit_sources.py`。Python 需要 pandas 与 openpyxl，R 需要 data.table 与 readxl。恢复脚本依赖本工作区原有 `workflow/source_legacy/python/utils/extractors.py` 与 `R/candidate/prepare_sample_boxplot_inputs.R` 的只读解析函数；所有新输出限定在本目录。
