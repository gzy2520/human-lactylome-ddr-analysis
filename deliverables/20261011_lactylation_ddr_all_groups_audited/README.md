# 乳酸化调控基因与DDR蛋白比例：全组核查版

打开 `Combined_RNA_DDR.pdf` 查看四页合并图。每页上方为RNA箱线图，中间为组织/细胞材料，下方为全蛋白组与乳酸化组DDR比例柱状图。

RNA面板逐列标注独立SampleID数；下方两类DDR柱分别标注各自来源单位数。t为技术制备，m为细胞模型，q为混合QC来源；p位置直接标“1 panel”，表示作者处理后的蛋白集合。来源记录、条件和来源列数见 `tables/ddr_source_unit_values.csv`；共享RNA参考见 `tables/display_mapping_31.csv`。RNA和DDR来源分别统计。

本版保留31个材料、62根DDR柱、307个分析单位，36个柱位具备SEM。4个完整处理矩阵按panel展示并取消样本SEM。4个细胞乳酸化材料的7个样本去除无正常序列支持的反向蛋白，分母各减少1，DDR分子不变。RNA表达值保持原版不变。

- `算法核查.md`：中文核查报告。
- `tables/all_62_algorithm_comparison.csv`：所有62位置的历史算法、当前算法、分子分母与样本拆分对照。
- `tables/rna_sample_counts_31.csv`、`tables/ddr_count_labels_62.csv`：图中数量标注。
- `tables/final_release_validation.csv`：28项独立数据核验。
- `SHA256SUMS.tsv`：交付文件校验值。

复现：在项目工作树执行 `Rscript workflow/lactylation_ddr/apply_all_group_audit.R`，再执行 `Rscript workflow/lactylation_ddr/render.R`。原始重建脚本与证据位于 `audit/20261011_all_groups_algorithm_audit/`。基因文献依据沿用上一版 `deliverables/20261011_lactylation_ddr_tall_verified/` 的记录。
