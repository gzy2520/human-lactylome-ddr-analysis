# 乳酸代谢探索交付

开始：2026-09-20；完成：2026-09-21。

先读 [中文汇报](REPORT_老师汇报.md)，方法见 [METHODS.md](METHODS.md)。

23套图在 figures/；逐图源数据位置见 tables/figure_sources.csv；全部表在 tables/；冻结证据在 evidence/。源数据库响应在仓库 config/lactate_metabolism/raw/，由 input_sha256.csv 追溯。

本包是表达探索，不代表乳酸浓度、代谢通量、Kla修饰强度或因果机制。显著性方法未调整；未训练预测模型。

数值复现入口：bash workflow/lactate_metabolism/run.sh <新空输出目录>。
