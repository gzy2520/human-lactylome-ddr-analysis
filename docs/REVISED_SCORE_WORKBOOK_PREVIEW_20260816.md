# 老师修订评分表的独立线性图预览

## 定位

输入为：

`data/identifier/乳酸化DDR基因评分表_Revised_20260816.xlsx`

本版本只用于观察老师修订评分后的4+1线性图效果。当前正式线性图仍来自399蛋白的
直接BP、CC、MF GO term分类；GO版本的脚本、表格、图和Methods均未被覆盖。

## 输入审计

- `评分表`含507个唯一UniProt`BaseAccession`，无重复。
- 当前30组Kla∩DDR membership的399个蛋白全部被评分表覆盖。
- 评分表中其余108个蛋白不进入本次试绘。
- BER、NER、MMR、FA、HR、NHEJ和AEJ仅含`-1`、`0`或`+1`，无缺失。
- 工作簿`score`与以下七通路权重逐行一致：
  BER 1、NER 2、MMR 3、FA 4、HR 5、AEJ 6、NHEJ 7。
- `Chromatin Interaction`和`Others`不进入七行主图或加权排序分数。

与上一版人工评分表相比，当前399蛋白中有12个蛋白发生变化，共涉及19个
蛋白–通路状态；12个蛋白的加权排序分数随之改变。逐蛋白差异表随结果一并输出。

## 作图规则

- 仍使用当前四分类membership和4+1顺序：
  non-tumor tissues 183、tumor tissues 178、cancer cell lines 381、
  normal cell lines 292和全集399。
- 每个集合独立按工作簿加权`score`升序排列；同分以`BaseAccession`稳定排序。
- `+1`使用相应通路实色，`-1`使用深炭灰，`0`使用浅灰。
- 行顺序和颜色与当前GO-term图一致：
  BER、NER、MMR、FA、HR、NHEJ、AEJ。
- 线性矩阵和summary分开生成；中英文均导出PNG、PDF和SVG。

## 复现

```bash
Rscript workflow/run_pipeline.R revised_score_preview
```

输出：

- 图：`results/figures/five_set_pathway_matrix_revised_excel_20260816/`
- 表：`results/tables/five_set_pathway_matrix_revised_excel_20260816/`
- 完整性检查：`tests/validate_revised_score_preview.R`

该预览不进入`tests/validate_publication_contract.R`或
`results/reports/publication_manifest_sha256.csv`。
