# 当前工作交接：31组 + bulk RNA-seq

完成配置：2026-09-14。工作树创建于2026-09-12。
分支：research/kla31-rnaseq-20260912；起点main bc02d1d。

1. 保留31组（9非肿瘤组织、3肿瘤组织、12癌细胞组、7非癌细胞/培养模型），ESCC不得再次遗漏。31组不是31个独立队列。
2. 蛋白分析以BaseAccession，RNA以Ensembl/Entrez；Symbol只展示。优先R，随机种子25。
3. 修正包原为冻结；2026-09-16经用户明确授权解除，因为冻结的是已被拒的初稿版本，现按31组重新分析。新脚本放R或workflow，新结果放outputs下独立日期任务目录。data/publication_input已由30组更新为31组（见该目录SCOPE.md），原30组版本备份在work/20260916_publication_input_30group_backup/；由旧30组输入派生的图表（corrected_final_result_20260905_sample_only/、audits/20260905_final_result/）尚未重算，重新使用前必须先重建。
4. results现在是独立空目录；renv/library为独立目录（包缓存链接保留），不再通过整个目录链接到旧工作树。勿直接修改包缓存内部文件，依赖按renv管理。
5. 本地扩展输入为独立副本。原始大数据未迁入；需全量重建时先明确源根与独立输出位置。旧纠正脚本仍有KLA_SOURCE_ROOT等参数，不能假定无参数就是当前新分析入口。
6. 旧add_data审计及Gemini原文是历史记录；RNA候选入口看docs/MATCHED_RNASEQ_REFERENCE_CATALOG_31.md，2026-09-14补充版为outputs/20260914_reference_material_audit/rnaseq_reference_candidate_31.csv。全部AnalysisReady=FALSE。按同组织/取材类别纳入，条件差异单独记录；本轮新增瘢痕邻近皮肤、PC-3M及MES28候选。普通全蛋白MCF10A来源参数为Protein FDR=1、PSM FDR=0.01，按老师确认不另加概率/FDR阈值但必须披露；HK-2内部/公开run别名未锁，见同目录REPORT.md；未改冻结结果。
7. RNA下载已转到服务器`192.168.3.45:/home/user/gzy/kla31-rnaseq-20260914`，合同为`config/rnaseq_server_download_20260914.tsv`；GEO/GTEx/ENCODE、GDC STAR counts、GEO raw-read解析和DepMap表达队列均有独立日志与SHA/MD5记录。DepMap因26Q1官方无验证码目录对这些文件不提供URL，固定使用官方目录提供的DepMap Public 24Q4（2024-12-16）Figshare链接，并保存目录快照；不与26Q1混用。下载后仍须锁定样本/文件版本/建库、做稳定ID和单位QC；处理Pro/Inh配对和研究来源复用。外部亲本RNA不能当KO/感染条件匹配。
8. 提交应明确描述变更，只提交代码、合同和应发布的结果；原始下载、临时输出、私人草稿不混入。推送当前研究分支，不推main，不自动合并。

启动检查：Rscript --vanilla workflow/preflight_kla31.R
历史main说明已归档至docs/history。
