# Kla项目新对话交接Prompt

继续处理项目：

`/Users/gzy2520/Desktop/Research/kla`

首先阅读：

1. `PROJECT_INDEX.md`
2. `manuscript/methods/METHODS_ZH.md`
3. `docs/DATA_PROVENANCE.md`
4. `docs/GO_TERM_PATHWAY_SCORING_30GROUPS.md`
5. `tests/validate_publication_contract.R`

不要把`archive/`、旧`results/reports/`或含33group、507、UMAP、t-SNE、PCA、
Cytoscape、V1/V2/V3、manual score字样的文件当成当前项目状态。

## 当前不可擅自改变的范围

- 40个候选样本组，37组有Kla定量，原有33个可配对组，当前主分析30组。
- 当前30组由non-tumor tissues 9组、tumor tissues 2组、
  cancer cell lines 12组、normal cell lines 7组构成。
- 普通全蛋白热图显示28条唯一参照行。
- 当前Kla∩DDR并集为399个唯一`BaseAccession`。
- 4+1集合按上述四类加全集排列，大小为183/178/381/292/399。
- 老师最后排除`PXD055230`、`PXD057709`和`PXD014870/MCF7`。

## 当前分析合同

- `BaseAccession`是唯一分析键；GeneSymbol和ProteinName只用于显示/审计。
- GO-DDR输入为`data/annotations/GO-repair+damage(human).tsv`，排除`NOT`。
- 普通全蛋白和Kla使用各自分母；Kla信号不能替代普通全蛋白信号。
- 参照的材料身份、实验状态、同批样本或独立队列关系必须分别表述。
- 四套Venn分别为全部Kla、Kla-DDR、普通全蛋白、普通全蛋白-DDR。
  最终图为固定几何、非面积比例图；图中数字来自精确membership和region表。
- 通路评分不再读取旧人工评分Excel。对399蛋白的全部直接BP、CC、MF GO term
  分配BER、NER、MMR、FA、HR、NHEJ、AEJ或`Others`；一个term允许多通路。
- 单蛋白的通路得分是该通路的不同直接GO term数；七通路term数之和只用于排序。
- 当前图形范围仅为柱状图、两份热图、四套Venn和4+1线性图/summary。
  不重跑UMAP、t-SNE、PCA或Cytoscape。
- `乳酸化DDR基因评分表_Revised_20260816.xlsx`对应的4+1图只是老师评分的
  独立试绘版；不能用它覆盖当前直接GO-term线性图、Methods或发表合同。

## 当前数据流

1. `python/data_preparation/build_core_kla_inputs.py`
2. `R/data_preparation/build_kla_regulator_landscape.R`
3. `R/analysis/analyze_regulator_reference_intensity.R`
4. `R/analysis/analyze_regulator_kla_intensity.R`
5. `R/data_preparation/build_reference_material_audit.R`
6. `R/analysis/analyze_ddr_fraction.R`
7. `R/figures/plot_four_class_venn.R`
8. `R/data_preparation/build_teacher_review_table.R`
9. `R/data_preparation/build_go_term_pathway_scores.R`
10. `R/figures/plot_five_set_pathway_matrix_go_term.R`
11. `R/analysis/summarize_four_class_venn_counts.R`
12. `tests/validate_publication_contract.R`

独立评分表试绘入口：

```bash
Rscript workflow/run_pipeline.R revised_score_preview
```

其说明为`docs/REVISED_SCORE_WORKBOOK_PREVIEW_20260816.md`。

统一命令：

```bash
Rscript workflow/run_pipeline.R selected_figures
Rscript workflow/record_environment.R .
Rscript workflow/build_manifest.R .
```

完成修改后必须核对：40/37/30/28、9/2/12/7、399、
183/178/381/292/399、10,605个蛋白–GO配对、2,785个唯一term、
103个七通路term、8个多通路term以及35行summary。对PNG进行人工目视检查，
再执行`git diff --check`、精确暂存、提交并推送。
