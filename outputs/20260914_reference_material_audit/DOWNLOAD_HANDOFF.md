# RNA 下载交接（2026-09-14）

RNA 候选的下载已转移到服务器 `192.168.3.45` 的独立目录：

`/home/user/gzy/kla31-rnaseq-20260914`

本目录不触碰服务器上已有项目，也不把原始数据写入 Git。下载合同是当前工作树的 `config/rnaseq_server_download_20260914.tsv`，服务器副本在 `metadata/download_specification.tsv`。

RNA 下载任务分别写入服务器 `logs/`：

- `download_20260914.log`：GEO 的 SOFT、series metadata、表达补充文件；GTEx V8 counts/TPM 和样本/供体注释；ENCODE HUVEC 两个 GRCh38 RSEM gene-quantification 文件。状态和 SHA-256 在 `metadata/download_status.tsv`。
- `gdc_download_20260914.log`：GDC 清单中的 1,017 个 STAR-count UUID，按 `config/gdc_star_counts_manifest_20260914.tsv` 的官方 MD5 下载。已全部完成；状态和实际 MD5 在 `metadata/gdc_star_counts_download_status.tsv`。
- `geo_raw_resolution_20260914.log`：将清单中的 GSM/SRX 解析为 SRA run 和 ENA FASTQ。已解析全部137个 FASTQ；下载队列仍在运行，未解析项会保留为显式状态，不会被猜测替代。

GSE163787 的 GEO 矩阵目录实际使用 GPL20301/GPL21103 两个文件名，已在服务器用平台文件补取；选定的 GSM4987488/SRX9725233 对应 GPL20301。历史默认文件名失败记录保留在状态历史中，最新状态快照已标记为被平台文件替代。

DepMap 另有一条独立下载队列，合同为 `config/depmap_expression_20260914.tsv`，服务器副本在 `metadata/depmap_expression_20260914.tsv`，日志为 `logs/depmap_download_20260914.log`。当前 26Q1 官方无验证码目录能返回文件名和 MD5，但这些文件的 URL 为空；因此按官方目录中可复现的固定链接锁定 `DepMap Public 24Q4`（2024-12-16），并在服务器保存目录快照 `metadata/depmap_no_captcha_download_catalog_20260914.csv`。这属于明确标注的版本回退，不与 26Q1 混用。

原始数据目录为服务器 `raw/`。PC-3M 的 `KI169`/`KI169Ctrl` 处理臂以及 m6A 位点表放在 `quarantine/`，不进入候选表达数据。服务器任务和所有下载完成后，仍需按 GSM/病例/细胞模型核对样本身份、建库与物种，统一到 Ensembl/Entrez 稳定 ID并核对 counts/FPKM/TPM 单位；因此 31 组注册表和补充候选表继续全部保持 `AnalysisReady=FALSE`。
