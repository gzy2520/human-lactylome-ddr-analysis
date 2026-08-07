# 当前版 KLA-DDR 分析流程说明

## 1. 目的

本流程用于整理乳酸化相关蛋白中的 DNA damage repair / DDR 候选基因，
并将其按修复机制做保守分类，最终生成 UMAP + 饼图可视化和对应表格。

本目录是与旧流程隔离的独立版本，不依赖根目录下的历史脚本和历史输出。

## 2. 输入

### 2.1 乳酸化蛋白来源表

路径：

- `data/article_tables/PXD014870.csv`
- `data/article_tables/PXD028488.csv`
- `data/article_tables/PXD050470.csv`
- `data/article_tables/PXD053474.csv`

每个表提供乳酸化相关蛋白/基因信息，至少包含以下列中的一组：

- `GeneSymbol` 或 `Symbol`
- `UniProtKB` 或 `UniProt`

同时保留：

- `ProteinName`
- `KlaSites`

### 2.2 GO repair/damage 注释表

路径：

- `data/GO-repair+damage(human).tsv`

这是 human 的 GO repair / DNA damage 相关表，用作主证据来源。

### 2.3 旧布局对照表

路径：

- `data/legacy_layout/umap_8type_data.csv`

用途仅限于：

- 保留旧图空间位置作为对照
- 让新增流程可与旧图比较

### 2.4 MSigDB gene set 文件

路径：

- `data/genesets.tsv`

默认不参与主分析。
仅在显式运行 `--include-msigdb` 时作为敏感性分析补充证据。
因为 MSigDB gene set 较宽，加入后会显著增加多机制归类，不建议作为默认主图依据。

## 3. 分析流程

### 3.1 候选基因整理

脚本 `run_workflow.py` 会：

1. 读取 4 个 PXD 表
2. 统一基因名为大写 `Symbol`
3. 统一 UniProt accession，并去掉 isoform 后缀，生成 `BaseAcc`
4. 按 `BaseAcc` 合并重复项

输出：

- `results/kla_ddr_source_provenance.csv`
- `results/kla_ddr_unique_go_repair_damage.csv`
- `results/gene_symbol_conflicts.csv`

其中 `gene_symbol_conflicts.csv` 用于记录同一 UniProt accession 对应多个来源基因名的情况，
便于后续人工核查。

### 3.2 GO 证据筛选

主证据来源是 `GO-repair+damage(human).tsv`。

筛选规则：

- 只保留 human 条目
- 不按 evidence code 过滤
- 排除 `QUALIFIER` 中含 `NOT` 的记录
- 用 GO term / GO name 将证据分到不同修复类别

### 3.3 分类口径

当前采用 9 类：

- `HR`
- `NHEJ`
- `BER`
- `NER`
- `MMR`
- `TLS`
- `DRR`
- `CP`
- `Other`

原则是尽量保守、可解释：

- 只有明确指向具体机制的 GO term 才进入机制类
- `CP` 只放 DDR signaling、checkpoint、泛 DNA repair
- `Other` 放确实和 DNA damage/repair 相关，但不应硬塞进机制类的 term

`Other` 主要包括：

- interstrand cross-link repair
- DNA alkylation repair
- protein-DNA covalent cross-linking repair
- mitochondrial DNA repair
- telomere maintenance in response to DNA damage
- DNA repair-dependent chromatin remodeling

分类规则会写入：

- `results/classification_rules.csv`

### 3.4 证据矩阵和主分类

分类后会生成：

- `results/gene_repair_evidence_long.csv`
- `results/gene_repair_category_matrix.csv`
- `results/gene_primary_category.csv`
- `results/repair_category_counts.csv`
- `results/evidence_code_counts.csv`

其中：

- `gene_repair_evidence_long.csv` 是长表证据记录
- `gene_repair_category_matrix.csv` 是基因 x 类别 0/1 矩阵
- `gene_primary_category.csv` 给每个基因一个主分类
- `repair_category_counts.csv` 给出每一类的基因数

主分类并不是简单按固定类别顺序选择。当前规则为：

1. 先按具体机制证据条数选择 primary mechanism
2. 若机制证据条数相同，优先选择更具体的机制 GO term
3. 若仍相同，再按 evidence code 权重排序
4. 若仍无法区分，则保留 `Primary_Is_Tie` 和 `Primary_Tie_Categories`

例如 `XRCC5/Ku80` 同时有 `GO:0000725 recombinational repair` 和
`GO:0006303 double-strand break repair via nonhomologous end joining`。
后者是更具体的 NHEJ 机制 term，因此 primary 归为 `NHEJ`，而不是依赖
固定类别顺序归入 `HR`。

### 3.5 UMAP 重算

这版 UMAP 不再沿用旧图的类别比例坐标，而是基于更细的 GO 证据特征重算：

1. 先构建 GO term 二值矩阵
2. 再叠加各类别的 log-count 特征
3. 用 UMAP 重新拟合空间布局

对应输出：

- `results/umap_go_term_feature_matrix.csv`
- `results/umap_pie_data.csv`

说明：

- `UMAP_1` / `UMAP_2` 为本版重算坐标
- `Legacy_UMAP_1` / `Legacy_UMAP_2` 为旧布局对照坐标
- `Has_Legacy_Coordinates` 标记是否能对上旧图

## 4. 作图

脚本：

- `plot_umap.R`

输出两张主图：

- `figures/umap_repair_pies.png`
- `figures/umap_repair_primary.png`

同时也会输出 PDF：

- `figures/umap_repair_pies.pdf`
- `figures/umap_repair_primary.pdf`

另外会生成旧布局对照图：

- `figures/umap_repair_pies_legacy_layout.png`
- `figures/umap_repair_pies_legacy_layout.pdf`

### 4.1 饼图

每个点对应一个基因，饼图扇区显示各类别证据比例。

### 4.2 Primary 图

主分类图中：

- 颜色表示 primary category
- 形状区分是否存在多个具体机制证据

## 5. 产出文件

### 5.1 结果表

- `results/kla_ddr_source_provenance.csv`
- `results/kla_ddr_unique_go_repair_damage.csv`
- `results/gene_repair_evidence_long.csv`
- `results/gene_repair_category_matrix.csv`
- `results/gene_primary_category.csv`
- `results/repair_category_counts.csv`
- `results/repair_category_evidence_source_counts.csv`
- `results/evidence_code_counts.csv`
- `results/gene_symbol_conflicts.csv`
- `results/classification_rules.csv`
- `results/umap_go_term_feature_matrix.csv`
- `results/umap_pie_data.csv`
- `results/umap_plot_data.csv`
- `results/genes_without_legacy_umap_coordinates.csv`

### 5.2 图

- `figures/umap_repair_pies.png`
- `figures/umap_repair_pies.pdf`
- `figures/umap_repair_primary.png`
- `figures/umap_repair_primary.pdf`
- `figures/umap_repair_pies_legacy_layout.png`
- `figures/umap_repair_pies_legacy_layout.pdf`

## 6. 重跑命令

默认主分析：

```bash
python3 run_workflow.py
Rscript plot_umap.R
```

加入 MSigDB 敏感性分析：

```bash
python3 run_workflow.py --include-msigdb
Rscript plot_umap.R
```

## 7. 解释注意点

- 这版主图的坐标是重算的，不是旧图硬继承的
- 旧布局只用于对照
- `Other` 不是“垃圾桶”，而是无法科学归入 7 个具体机制时的保守分类
- 主图的标签和点大小是展示层参数，分类结果不受影响
- GO evidence code 未过滤，因此结果适合做注释概览，不应表述为高置信功能验证
- MSigDB 版是敏感性分析，不能和 GO-only 主分析混为同一证据强度
- UMAP 只表示当前证据矩阵下的相似性，不代表真实生物距离或通路层级
