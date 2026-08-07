# Kla 重新分析方法

## 数据范围

纳入 PXD014870、PXD028488、PXD050470、PXD053474、PXD060185、
PXD078013 和 PXD078736。PXD038880/PXD050906 的下载文件和元数据保留，
但不进入 Kla、GO 或 Venn 分析。

## 证据读取

分析不把 `.raw` 或 `.d.zip` 直接当作 Kla 证据。不同软件和数据集使用独立解析器：

- PXD014870：四组 MCF7 的 `Lactyl (K)Sites.txt`，并用
  `modificationSpecificPeptides.txt` 中的 Kla site ID 验证。
- PXD028488：PEAKS `DB search psm.csv` 中的 Lactylation/PTM/AScore 或
  `K(+72.02)`、`K(Lactyl)`，再通过 `protein-peptides.csv` 定位蛋白位点。
  BV2/RAW 小鼠目录排除；HCT116 non-enrichment 因缺少 PSM 表排除；all-HCD
  聚合目录只用于验证。marker 156 只作辅助支持，不能单独定义 Kla。
- PXD050470：作者 Supplementary Tables S3/S12，不重新检索缺失的 27 个 raw。
- PXD053474：分别读取 enriched/unenriched DDA 与 DIA，并与作者 S3 对照。
  主结果保留 S3 全部位点，再加入同时获得 DDA 和 DIA 支持的 search-only 位点。
- PXD060185：`RESULT/combined/txt/La (K)Sites.txt`，按 A/B/C/D 保留
  MCF7、MDA-MB-468、MCF10A、T-47D 的样本来源。
- PXD078013：联合 `evidence.txt` 和 `proteinGroups.txt`；要求 evidence 的
  `La (K)>0`、存在 site ID，并能映射到 proteinGroups 的 site position。
- PXD078736：`La(K)Sites.txt`，按 ctrl/man 和六个重复保留样本来源。

主分析定位概率阈值为 0，即不额外删除作者已报告的 Kla 位点；仍排除 reverse、
contaminant、非人源和无可追溯 Kla 位点的记录。PXD014870 另存阈值 0.75 的敏感性结果。

## 蛋白与 GO

读取阶段保留 sample-level 长表，不提前合并样本、实验组或重复。蛋白集合在最后一步
按 UniProt `BaseAccession` 去重，isoform 后缀去除。GO repair/damage 先按
BaseAccession 匹配，未匹配时才以 GeneSymbol 辅助；`NOT` 注释排除，GO evidence code
保留但主分析不按 evidence code 缩减。

## 分组与 Venn

最终确认的三个集合为：

- `hippocampus_tissue`：人海马体组织。
- `normal_immortalized_cell_lines`：HEK293T、HK-2、MCF10A。
- `tumor_cell_lines`：MCF7、HCT116、T-ALL、MDA-MB-468、T-47D、RKO。

HEK293T/HK-2 是转化或永生化模型，不描述为正常组织；MCF10A 是第三个非肿瘤/永生化
模型，MCF7 归肿瘤细胞系。Venn 的计数单位是唯一 BaseAccession，每张图的七个精确区域
均有对应 CSV。

## 可复现命令

```bash
python3 reanalysis/scripts/build_project_metadata.py --project-root /Users/gzy2520/Desktop/Research/kla
PYTHONPATH=reanalysis/scripts python3 reanalysis/scripts/run_pipeline.py --project-root /Users/gzy2520/Desktop/Research/kla
PYTHONPATH=reanalysis/scripts python3 -m unittest discover -s reanalysis/tests -p 'test_*.py' -v
```
