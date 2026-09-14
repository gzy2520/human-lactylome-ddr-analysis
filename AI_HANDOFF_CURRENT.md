# 当前工作交接：31组 + bulk RNA-seq

完成配置：2026-09-14。工作树创建于2026-09-12。
分支：research/kla31-rnaseq-20260912；起点main bc02d1d。

1. 保留31组（9非肿瘤组织、3肿瘤组织、12癌细胞组、7非癌细胞/培养模型），ESCC不得再次遗漏。31组不是31个独立队列。
2. 蛋白分析以BaseAccession，RNA以Ensembl/Entrez；Symbol只展示。优先R，随机种子25。
3. 已有修正包冻结；新脚本放R或workflow，新结果放outputs下独立日期任务目录。不得运行旧30组流程覆盖已有结果。
4. results现在是独立空目录；renv/library为独立目录（包缓存链接保留），不再通过整个目录链接到旧工作树。勿直接修改包缓存内部文件，依赖按renv管理。
5. 本地扩展输入为独立副本。原始大数据未迁入；需全量重建时先明确源根与独立输出位置。旧纠正脚本仍有KLA_SOURCE_ROOT等参数，不能假定无参数就是当前新分析入口。
6. 旧add_data审计及Gemini原文是历史记录；RNA当前候选看docs/MATCHED_RNASEQ_REFERENCE_CATALOG_31.md。全部AnalysisReady=FALSE，不能自动进入跨组学推断。
7. 下一步锁定RNA样本/文件版本/建库、下载并做稳定ID和单位QC；处理Pro/Inh配对和研究来源复用。外部亲本RNA不能当KO/感染条件匹配。
8. 提交应明确描述变更，只提交代码、合同和应发布的结果；原始下载、临时输出、私人草稿不混入。推送当前研究分支，不推main，不自动合并。

启动检查：Rscript --vanilla workflow/preflight_kla31.R
历史main说明已归档至docs/history。
