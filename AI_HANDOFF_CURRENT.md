# 当前交接：31 组蛋白组 + 28 份 RNA 参考

更新：2026-09-20（修复启动于 2026-09-19）。分支：`research/kla31-rnaseq-20260912`。

1. 当前 RNA 发布入口是 `config/rna_current_release.csv`，修复验收见 `audit/20260919_rna_repair/REPORT.md`。旧审计见 `audit/20260919_project_review/REPORT.md`；旧目录不是当前交付。
2. RNA 是 28 份不重复 ReferenceKey、1,898 个唯一样本，与 31 条蛋白组 GroupID 对应。HCT116 三行、HK-2 两行共用 RNA，在核心组图中只计一次。保留 31 行映射用于材料对照；SampleID 唯一不等同于供体独立。
3. 蛋白组使用 BaseAccession，RNA 使用 Ensembl/Entrez；Symbol 只展示或用于原始符号源的一次性 ID 转换。优先 R，seed=25。
4. 用户于本轮明确要求“显著性不用管”：不更换单因素/双因素 ANOVA，不修复配对模型、不重新选择检验。输入重建后原方法统计量随输入重新计算；此保留不表示统计假设已经验证。
5. 修复 FPKM/RPKM 重复长度归一化、陈旧缓存、部分元数据覆盖、图标题和 NA 显示、文档错误及混版交付。Stage 2/渲染都要求独立输出目录；不要覆盖历史对象，不要把旧目录的文件整体复制进新发布。
6. 大型原始 RNA 数据在 `user@100.121.229.123:/home/user/gzy/kla31-rnaseq-20260914`，SSH 使用 `HostKeyAlias=192.168.3.45`。初次提取脚本保存在 `workflow/releases/20260919_repair/`，后续重跑使用完整验证版 `workflow/releases/20260920_validated/`；本轮 expression 为 `outputs/20260919_expression_corrected`。本地同步主要矩阵和元数据，原始下载和对象保留服务器。
7. 本轮没有修改蛋白组分析或旧冻结包。`final_result` 仍是旧冻结包快捷入口，不能当作当前 RNA 发布；原冻结 163 文件哈希检查继续保留。
8. `workflow/preflight_kla31.R` 已改为检查当前发布和当前输入。旧检查在 `workflow/preflight_kla31_legacy_20260912.R`，其中候选路径与 AnalysisReady=FALSE 要求只是历史合同。
9. 来源局限、NSC 模型匹配、胎盘取样限制、TALL-104 n=1 等没有被数值修复消除，详见方法文档。不得将外部 RNA 当作蛋白组同一样本的配对定量。
10. 脚本放 workflow/R，结果放 dated outputs；发布以 SHA-256 清单固定。正式绘图器不再自动写桌面。提交和推送仅当前研究分支，不改 main、不合并。

运行：`Rscript --vanilla workflow/preflight_kla31.R`。

本地全链重建：`bash workflow/rebuild_rna_release_20260919.sh <已验收的expression目录> <新日期标签>`。输出完成仍需视觉 QA，再更新当前发布清单。

## 乳酸代谢探索（2026-09-21）

已按老师要求完成第一阶段来源/去向表达分析。新包：`outputs/20260920_lactate_metabolism/`；中文汇报在包内 `REPORT_老师汇报.md`，方法为 `docs/LACTATE_METABOLISM_METHODS_20260920.md`，审计为 `audit/20260920_lactate_metabolism/REPORT.md`。128实测基因、6 GO集合、10机制模块、23套图；17张表独立重跑一致。现有RNA发布入口不变。新图使用实际RNA对照条件标签，保留31条蛋白映射；不把Kla检出数当强度，不拟合相关/机器学习。BPH背景用药与跨研究比较须披露。

## B-R 机器学习当前入口（v4，2026-09-21）

以 `outputs/20260921_br_v4/` 和 `audit/20260921_br_v4/REPORT.md` 为当前修复版；冻结参数为 `config/lactylation_ml_review/model_config_v4.md`，入口 `workflow/lactylation_ml_review/run_v4.sh`（输出已存在会拒绝重建，应另用隔离副本）。v2/v3保留为历史，不沿用v3手写的“15/19折变差”（实际9/19）。v4内层按连通组等权选参；主模块为4基因L-乳酸/丙酮酸转换表达，LDHD的D-乳酸氧化分开展示。备用LDHA ID在28源矩阵均缺失，覆盖不能按NA占位行计数。`outputs/20260921_br_v4_weight_only/` 仅作权重修复对照。此任务是DDR蛋白流程检出预测，非乳酸浓度或修饰占有率；不改变正式RNA发布入口及冻结蛋白组结果。
