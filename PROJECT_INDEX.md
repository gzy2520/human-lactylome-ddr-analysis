# Human Kla蛋白组发表项目：当前项目索引

更新时间：2026-08-16

本文件是当前项目状态的主要入口。判定分析范围时，优先级依次为：

1. `config/`中的范围和配对配置；
2. 当前执行脚本；
3. `tests/validate_publication_contract.R`；
4. 本文件、双语Methods和当前审计文档；
5. 可重新生成的`results/`。

`archive/`、Git历史、旧`results/reports/`以及文件名含33group、507、UMAP、
t-SNE、PCA、Cytoscape、V1/V2/V3或manual score的材料仅供追溯，不是当前依据。

## 当前固定范围

- 来源目录：40个候选`PXD+SampleGroup`。
- 可进行Kla逐蛋白定量：37组。
- 在老师最后范围调整前具有参照的配对组：33组。
- 当前主分析：30组、28条唯一普通全蛋白参照展示行。
- 四分类顺序：non-tumor tissues 9组、tumor tissues 2组、
  cancer cell lines 12组、normal cell lines 7组。
- 当前Kla∩DDR并集：399个唯一`BaseAccession`。
- 4+1集合顺序与大小：183、178、381、292、399。

老师明确从原33组中排除：

- `PXD055230 / human fibroblasts mock and HCMV or HSV-1`；
- `PXD057709 / human fibroblasts mock and HCMV`；
- `PXD014870 / MCF7`。

另外，40组目录中有3个PXD037371临床组无法可靠拆分TMT通道，4个定量组缺少
可接受的普通全蛋白参照。这些状态由配置和审计表保存，不应与老师最后排除的3组
混为一类。

## 当前分析规则

- 匹配、去重、GO交集和集合运算只使用去isoform的UniProt
  `BaseAccession`。
- GeneSymbol和ProteinName仅用于显示和审计。
- GO-DDR来自`data/annotations/GO-repair+damage(human).tsv`，排除`NOT`，
  主分析不按evidence code进一步缩减。
- 普通全蛋白必须来自非Kla富集数据。材料身份、同批样本、实验状态和独立队列
  关系分别记录，不能把所有参照统称为同条件对照。
- Kla和普通全蛋白DDR比例使用各自独立分母。
- Venn图使用固定四集合几何，面积不表达数量；数字来自精确membership。
- 通路评分使用399蛋白的全部直接BP、CC、MF GO term。term可同时属于多条通路；
  不命中七通路的term归`Others`。旧人工±1评分表和系数不参与当前结果。

## 当前代码入口

| 路径 | 当前用途 |
|---|---|
| `python/data_preparation/build_core_kla_inputs.py` | 构建可审计的核心Kla证据层 |
| `R/data_preparation/build_kla_regulator_landscape.R` | 建立40/37组目录和Kla热图输入 |
| `R/analysis/analyze_regulator_*_intensity.R` | 生成30行Kla和28行普通全蛋白热图数据 |
| `R/analysis/analyze_ddr_fraction.R` | 计算30组Kla与参照DDR比例 |
| `R/figures/plot_four_class_venn.R` | 生成四套精确membership Venn图与审计表 |
| `R/data_preparation/build_go_term_pathway_scores.R` | 生成399蛋白的直接GO term通路计数 |
| `R/figures/plot_five_set_pathway_matrix_go_term.R` | 生成4+1线性图和summary |
| `workflow/run_pipeline.R` | 统一运行入口 |
| `tests/validate_publication_contract.R` | 当前结果合同 |

## 当前目录

| 目录 | 用途 |
|---|---|
| `config/` | 样本范围、参照关系、ID映射、GO规则和显示配置 |
| `data/` | 本地原始/处理数据及少量Git跟踪注释输入 |
| `python/` | 异构Kla证据解析 |
| `R/` | 数据准备、统计和论文图 |
| `workflow/` | 执行、依赖、环境记录和SHA256清单 |
| `tests/` | 发表版结果合同 |
| `manuscript/` | 当前双语Methods和审计说明 |
| `docs/` | 当前来源、评分和复现说明 |
| `results/` | 可重新生成结果，Git忽略 |
| `work/` | 可重新生成中间表，Git忽略 |
| `archive/` | 历史材料，Git忽略，不是当前分析入口 |

## 运行和验证

```bash
cd "/Users/gzy2520/Desktop/Research/kla"
python3 -m pip install -r requirements.txt
Rscript workflow/install_r_dependencies.R
Rscript workflow/run_pipeline.R selected_figures
Rscript workflow/record_environment.R .
Rscript workflow/build_manifest.R .
```

当前发表范围不运行UMAP、t-SNE、PCA或Cytoscape。
