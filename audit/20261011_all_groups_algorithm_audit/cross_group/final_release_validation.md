# 全组审计版最终数据核验

核验对象：`deliverables/20261011_lactylation_ddr_all_groups_audited`。独立执行28项数据检查，**28通过、0失败**；逐项结果见 `final_release_validation.csv`。

- 307个分析观察单位、62个材料×组学位置、36处绘图SEM，实际表一致。
- 每样本/面板membership唯一；所有DDR标签直接匹配原冻结2660个GO accession。307个分母、分子、比例以及62个N、均值、SD、SEM独立回算一致。
- 4个processed_matrix均为N=1，无SD/SEM/PlotSEM，并带p标记。SourceObservationN为3、3、8、8；SourceObservationIDs完整保留原列名；每个panel集合与被合并的各旧列集合完全一致，柱高不变。
- 4个修正细胞Kla组的最终membership与同伴从原始记录独立重建的完整集合逐键一致。恰有7个来源样本各减少1个非DDR蛋白，DDR分子不变；被删ID的合格正向支持行均为0。
- 其余保留源样本的membership及分子分母全部与tall_verified版一致，没有额外静默改动。
- HK2的当前Condition正确为ctr/man，样本名前缀与Condition一致；显示为`HK-2 (ctr/man)*`。原SampleGroup仍含mannitol作为历史关联键，未将其用作当前处理条件释义。
- 5个RNA输入文件与tall_verified版字节完全一致，实际绘图RNA表内容完全一致：12×1898记录、28个参考来源，稳定ID及表达值均未改变。

已阅读 `workflow/lactylation_ddr/apply_all_group_audit.R`，其针对原记录修正版membership的替换、4个面板折叠、历史key保留及当前显示标签更新与源证据一致。本报告核验数据文件与计算结果；图正在生成，尚未包括最终渲染页面检查。
