# Tables audit — 2026-09-05

审计对象是冻结目录
`/Users/gzy2520/Desktop/Research/kla/results/final_figures_and_tables/Tables`。
本轮只读检查，没有修改或覆盖任何发表表格；脚本与结果写在当前 fork 的本目录。

## 结论

不是所有表都有问题。核心集合表的稳定 ID、行数和算术大多通过，但当前 Tables 包仍不能直接作为与图片同步的定稿，主要有四类问题：

1. **统计检验表不可直接采信**：Figure 1 和七通路统计表把 aggregate、pool、dataset_union、condition、model 混入 sample；七通路 Pro/Inh 是同一来源观测的配对测量，却按独立行拟合。`figure1_category_omnibus_anova.csv` 的 N=310（其中 38 条非 sample），`pathway_summary_two_way_anova.csv` 的 NPoint=98（其中 26 个来源观测非 sample，长表 N=196）。
2. **图表版本不一致**：冻结 Figure 1a 仍是旧 30 组输入（210 条观测、肿瘤 WP/Kla 为 13/6），但随包统计表是扩展 31 组（N=310、肿瘤为 107/12）；七通路图仍是旧 92 个来源观测/方向，而统计表是 NPoint=98。表本身与图不是同一次输入导出。
3. **Table S2 有未映射的稳定 ID**：PXD036307 normal human lung 有 595 条、PXD063047 normal pregnancy placenta 有 546 条 `ENSEMBLPROT` 的 `MappedBaseAccessions` 为空。摘要仍计入这些源蛋白：lung 为 14,022（可映射 accession 12,358；DDR 668 vs 606），placenta 为 12,395（可映射 accession 10,847；DDR 556 vs 502）。如果后续按 BaseAccession 做交集，这些行会静默丢失；必须显式保留为 unmapped 或在摘要中给出可映射分母。
4. **Table S6 调控因子表有记录/文献问题**：`Regulator_Annotations` 有 50 行、49 个唯一 `(Role, GeneSymbol, BaseAccession)` 记录、48 个唯一 accession；HDAC5/Q9UQL6（记录 23、24）同角色重复，其中一行没有文献。另有 12 行 `References` 为空，PARK7 记录 42 把两个 DOI 直接粘连、没有分隔符。

另外，S1/S2 的 31 行是 **31 个分组**，不是 31 个独立数据集对：当前表里 Kla PXD 只有 22 个、参考 PXD 只有 18 个；例如 PXD030304 被 8 个细胞系分组共用。若“31 对”在正文中被理解为 31 个独立配对，会造成分母和独立性误读。

## 通过的检查

- Table S1 和 S2 各有 31 个唯一 `(PXD, SampleGroup)`；S1 类别为 normal tissue 9、tumor tissue 3、cancer cell 12、normal cell 7。
- S1 的 `Group_Summary` 与两个 membership sheet 的唯一 `BaseAccession` 计数完全一致；S1/S2 的参考组摘要字段完全一致；membership 复合键没有重复。
- Table S3 与当前 `data/publication_input/human_ddr_go_annotations.tsv` 逐单元格（NA/空白/空格归一化后）完全一致，实际是 **5,111 行**；因此 README 中“6,707 条”是过时说明，不是当前 workbook 的行数。
- Table S4 四套 membership 的 `BaseAccession` 无重复，Region 与四个布尔集合列一致，`Set_Counts`/`Region_Counts` 可由明细重算。
- Table S5 四个 sheet 的 accession 集合与 Table S4 的 Kla-DDR 类别完全一致，状态只含 -1/0/+1，`SignedScore` 加权重算无误。
- Table 14 的 CSV/XLSX 完全一致，14 个 accession 唯一，`SignedScore` 算术无误。

S1 的 `RowOrder` 为 1–12、14–25、28–34，存在 13、26、27 的空档；这不改变计数，但应说明它是旧显示顺序而不是连续行号。

## 结果文件

- [`audit_tables.py`](audit_tables.py)：只读审计脚本
- [`workbook_inventory.csv`](workbook_inventory.csv)：工作簿、sheet、缺失值、重复行、公式清单
- [`s1_checks.csv`](s1_checks.csv)、[`s1_membership_reconciliation.csv`](s1_membership_reconciliation.csv)
- [`s2_checks.csv`](s2_checks.csv)、[`s2_mapping_reconciliation.csv`](s2_mapping_reconciliation.csv)
- [`s3_source_reconciliation.csv`](s3_source_reconciliation.csv)
- [`s4_checks.csv`](s4_checks.csv)、[`s5_checks.csv`](s5_checks.csv)
- [`s6_checks.csv`](s6_checks.csv)、[`s6_duplicate_records.csv`](s6_duplicate_records.csv)、[`s6_missing_references.csv`](s6_missing_references.csv)
- [`statistics_table_checks.csv`](statistics_table_checks.csv)：统计表的分母、观测粒度和模型风险
- [`cross_artifact_mismatch.csv`](cross_artifact_mismatch.csv)：统计表与冻结图片的版本错配
- [`table_hashes.csv`](table_hashes.csv)：本轮读取的表文件 SHA-256

与图片层的完整证据见 [`../20260905_final_result/README.md`](../20260905_final_result/README.md)。
