# Kla项目交接Prompt（发表版结构）

继续处理项目：

`/Users/gzy2520/Desktop/Research/kla`

先阅读 `PROJECT_INDEX.md`、`manuscript/methods/METHODS_ZH.md`、
`docs/DATA_PROVENANCE.md` 和 `tests/validate_publication_contract.R`。

## 当前范围

项目整合人源蛋白乳酸化质谱数据，并以普通全蛋白组作参照评估DDR相关蛋白。
来源目录含40个样本组，37组具有Kla定量，30组进入当前配对分析；四分类显示顺序为
非肿瘤组织9、肿瘤组织2、癌细胞系12、正常细胞系7。普通全蛋白热图为28条唯一参照行。

排除的4个无可用配对参照组为PXD062720、PXD063047重度子痫前期胎盘、
PXD064038和PXD075014。PXD037371三个临床组因TMT通道无法可靠映射，只保留
来源审计。老师另要求排除PXD055230、PXD057709和PXD014870。当前Kla∩DDR
蛋白并集为399，五集合按显示顺序为183/178/381/292/399。

## 不可改变的规则

- 匹配、去重、GO交集和集合运算只使用去isoform的UniProt `BaseAccession`。
- GeneSymbol和ProteinName仅用于显示或审计，不能作为命中回退。
- Ensembl protein ID必须经显式表转换至UniProt。
- GO-DDR输入为 `data/annotations/GO-repair+damage(human).tsv`，排除NOT；
  主分析不按evidence code缩减。
- 普通全蛋白信号不得由Kla/PTM富集信号替代。
- 随机种子固定25。
- 不移动或删除 `data/` 下大型原始数据。
- 旧V1/V2、碰撞比较、圆形矩阵和Cytoscape试验不属于当前发表流程。

## 代码层次

1. `python/data_preparation/build_core_kla_inputs.py` 解析7个核心PXD，生成
   `work/intermediate/kla_by_dataset/all_primary_sample_level_kla_sites.csv`。
2. `R/data_preparation/build_kla_regulator_landscape.R` 建立40/37组Kla目录。
3. 两个 `R/analysis/analyze_regulator_*_intensity.R` 生成28行普通全蛋白轴和
   30行Kla轴。
4. `R/data_preparation/build_reference_material_audit.R` 固定30/7参照状态。
5. `R/analysis/analyze_ddr_fraction.R` 与
   `R/figures/plot_four_class_venn.R` 生成DDR统计和四类集合。
6. `R/data_preparation/prepare_protein_function_inputs.R` 保留人工评分、GO和颜色底表。
7. `R/figures/plot_five_set_pathway_matrix.R` 从当前membership筛出399蛋白，
   生成最终4+1线性通路展示。
8. 本轮不运行UMAP、t-SNE、PCA和Cytoscape。

统一入口：

```bash
cd "/Users/gzy2520/Desktop/Research/kla"
python3 -m pip install -r requirements.txt
Rscript workflow/install_r_dependencies.R
Rscript workflow/run_pipeline.R selected_figures
```

完成后必须运行：

```bash
Rscript workflow/record_environment.R .
Rscript workflow/run_pipeline.R validate
Rscript workflow/build_manifest.R .
```

发表前必须确认37/30/28、9/2/12/7、399和183/178/381/292/399均通过结果合同，
并同时人工查看中英文柱状图、热图、Venn图和线性图是否截断。
