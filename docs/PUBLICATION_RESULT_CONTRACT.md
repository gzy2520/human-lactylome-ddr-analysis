# 当前发表结果合同

项目只有在`tests/validate_publication_contract.R`通过时，才可描述为当前可复现版本。

合同固定：

- 40个候选组、37个Kla定量组、30个当前配对组和28条唯一参照展示行；
- 四分类显示顺序为non-tumor tissues 9组、tumor tissues 2组、
  cancer cell lines 12组、normal cell lines 7组；
- 当前Kla∩DDR并集为399个唯一`BaseAccession`；
- 4+1集合大小依次为183、178、381、292和399；
- 四套Venn分析均可由membership和region表精确重建；
- Venn最终几何不按面积拟合，15个逻辑区域均保留，数字是唯一数量编码；
- 399蛋白含10,605条直接蛋白–GO配对和2,785个唯一直接term；
- 103个term命中至少一条七通路规则，8个term同时命中两条通路；
- `Others`与七通路互斥且共同覆盖全部2,785个term；
- 4+1乘7条通路产生35行summary；
- 中英文柱状图、两份热图、四套Venn、4+1线性图和summary均存在。

合同不比较PNG字节，因为字体渲染可随操作系统变化。仍需目视检查双语图是否存在
文字截断、字号异常、标签错位或意外标题。

当前合同不要求UMAP、t-SNE、PCA或Cytoscape；这些历史结果不进入当前SHA256
发表清单。

老师修订人工评分表生成的
`five_set_pathway_matrix_revised_excel_20260816`属于独立预览，同样不进入本合同
或发表SHA256清单。其完整性由`tests/validate_revised_score_preview.R`单独验证。
