# DDR逐样本恢复与原图核对

2026-10-11。当前合并图为`deliverables/20261011_lactylation_ddr_sample_rebuilt/`。本目录保存恢复脚本、原表证据、稳定UniProt membership以及核对结果。

## 核对结论

原272条记录中269条的蛋白分母与DDR分子完全一致，TALL乳酸化3条因原分母含未乳酸化蛋白而修正。新增53条来源单位显示记录后，共325条、62个材料×检测类型位置；其中40处绘制均值±SEM。数据来自当前原表，旧文档用于寻找来源，再依据真实列、metadata与采集信息判断单位。

- `base_rebuild/`：30组原解析函数逐蛋白重建，另外恢复ESCC Kla六样本。`all_membership`包含用于追溯的condition/model/union等记录；合并入口仅选择所需sample记录并替换TALL。
- `existing_samples/`：ESCC全蛋白组94例、PC3M Kla、TALL的独立复算；TALL最终采用`TALL_Kla_legacy_mapping`。`final_delivery_validation.csv`中20项检查全部通过。
- `kla/`：缺失11个Kla显示位置的26条来源单位；包括可拆分样本、单例、细胞模型及技术合并单位。
- `whole_proteome/`：缺失12个全蛋白组显示位置的27条主分析记录；使用`UseRole=primary`。另外两条HEK标准品只用于来源比较。

所有最终IsDdr均以冻结GO文件直接逐UniProt判定。可绘图误差由单位性质与条件内重复情况确定；图中的技术SEM、细胞模型与QC已单独标识。

## 从原表重新恢复

在工作区根运行以下入口。本地原始数据目录为`/Users/gzy2520/Desktop/Research/kla/data`；具体环境路径在各脚本开头设置。

```sh
Rscript audit/20261011_ddr_sample_recovery/base_rebuild/rebuild.R . /Users/gzy2520/Desktop/Research/kla/data
Rscript audit/20261011_ddr_sample_recovery/base_rebuild/rebuild_escc.R
Rscript audit/20261011_ddr_sample_recovery/existing_samples/audit_existing.R
python3 audit/20261011_ddr_sample_recovery/kla/extract_legacy_sources.py
Rscript audit/20261011_ddr_sample_recovery/kla/recover_kla_samples.R
python3 audit/20261011_ddr_sample_recovery/whole_proteome/recover.py
Rscript workflow/lactylation_ddr/rebuild_sample_inputs.R
Rscript workflow/lactylation_ddr/render.R
Rscript audit/20261011_ddr_sample_recovery/existing_samples/validate_final.R
```

R恢复依赖data.table、readxl；Python恢复依赖pandas、openpyxl。大型源位点中间表和原论文全文不纳入Git，位点表可由上述入口重建。压缩membership保留，可直接重新生成最终计数与图。

各子目录README记录逐来源的筛选条件、参考论文、分级或重复采集合并，以及定位阈值与标识映射的具体影响。
