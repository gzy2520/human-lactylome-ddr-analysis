# Human Kla proteomics publication repository

更新时间：2026-08-12

本仓库以老师确认的30组配对分析和399个Kla∩DDR蛋白为当前发表版主线。旧分析、
探索图和V1/V2版本不再位于当前代码树，必要时通过Git历史追溯。

## 目录

| 目录 | 用途 |
|---|---|
| `data/` | 本地原始数据、检索结果、补充表和人工评分表；大型文件不进入Git |
| `config/` | 固定样本范围、参照配对、ID映射和UniProt注释缓存 |
| `python/` | 7个核心PXD异构原始表解析，仅负责基础Kla证据层 |
| `R/data_preparation/` | 40/37/30组目录和通路功能输入准备 |
| `R/analysis/` | 强度、DDR占比和五集合降维参数分析 |
| `R/figures/` | 论文图和补充图 |
| `workflow/` | 统一运行、环境记录和SHA256清单 |
| `tests/` | 发表版结果合同 |
| `manuscript/` | 中英文Methods与审计说明 |
| `docs/` | 数据来源和分析决策记录 |
| `results/` | 可重新生成的表、图和报告；Git忽略 |
| `work/` | 可重新生成的中间表与日志；Git忽略 |

## 固定分析口径

- 40个来源样本组进入可审计目录，其中37组具有Kla定量。
- 30组具有可用普通全蛋白参照，四分类显示顺序为非肿瘤组织9、肿瘤组织2、癌细胞系12、正常细胞系7。
- 30个Kla组对应28条唯一普通全蛋白参照展示行。
- 配对分析排除PXD062720、PXD063047重度子痫前期胎盘、PXD064038和
  PXD075014。
- PXD037371三个临床组因TMT通道无法可靠映射，不进入Kla定量范围。
- 老师要求排除PXD055230、PXD057709和PXD014870。
- 当前399个Kla∩DDR蛋白；五集合按显示顺序为183/178/381/292/399。
- 分析键为去isoform的UniProt `BaseAccession`；GeneSymbol仅展示和审计。
- 所有随机步骤使用种子25。

## 一键运行

```bash
cd "/Users/gzy2520/Desktop/Research/kla"
python3 -m pip install -r requirements.txt
Rscript workflow/install_r_dependencies.R
Rscript workflow/run_pipeline.R selected_figures
Rscript workflow/record_environment.R .
Rscript workflow/build_manifest.R .
```

分阶段运行：

```bash
Rscript workflow/run_pipeline.R core
Rscript workflow/run_pipeline.R validate
```

## 发表版入口

- 中文Methods：`manuscript/methods/METHODS_ZH.md`
- 英文Methods：`manuscript/methods/METHODS_EN.md`
- 数据来源：`docs/DATA_PROVENANCE.md`
- 结果合同：`tests/validate_publication_contract.R`
- 完整交接：`NEW_CHAT_PROJECT_PROMPT.md`
