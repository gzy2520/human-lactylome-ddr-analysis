# 30组 Kla∩DDR 蛋白的 GO term–通路重建

## 当前分析对象

- 当前主分析范围为 30 个成对样本组。
- 分析键为去除 isoform 后缀的 UniProt `BaseAccession`，不使用 Gene Symbol 进行匹配或计算。
- 四分类并集包含 399 个 Kla∩DDR 蛋白；4+1 集合大小依次为 183、178、381、292 和 399。
- GO 注释来自 UniProt release 2026_02（Taxon ID 9606），使用每个蛋白的直接 BP、CC 和 MF term。

## GO term 分类算法

1. 从 399 个当前蛋白提取所有直接 GO 注释，共 10,605 个蛋白–term 配对和 2,785 个唯一 term。
2. 为 BER、NER、MMR、FA、HR、NHEJ 和 AEJ 建立人工审阅的 GO 种子规则。
3. 仅在种子具有明确通路特异性时扩展其 GO 后代；本体层级不足以表达通路特异性时使用精确 term 规则。
4. 一个 GO term 可以同时归入多条通路。例如：
   - recombinational interstrand cross-link repair：FA 和 HR；
   - nucleotide-excision repair involved in interstrand cross-link repair：NER 和 FA；
   - generic nonhomologous end joining 及其调控：NHEJ 和 AEJ。
5. 未命中七条通路规则的 term 归为 `Others`。宽泛的 DNA repair、DNA damage response、DNA binding、DNA ligase activity 等不会仅因名称相关而被强行分配。
6. 正调控和负调控都表示蛋白“涉及”对应通路，不再解释为促进或抑制，也不再使用 ±1。

规则表为 `config/go_term_pathway_seed_rules.csv`。每个非 Others term 均保存命中的种子 term、精确/后代关系、GO 名称和定义，支持逐条复核。

## 计分与作图

- 单个蛋白在某条通路的得分，为其直接 GO 注释中归入该通路的不同 term 数。
- 七通路总分为七个通路 term 数之和；多通路 term 会在其涉及的每条通路各计一次。
- 线性图在每个集合内按七通路总分升序排列；同分时仅用 `BaseAccession` 稳定排序，避免人为制造通路色块。
- 图中只编码是否至少有一个 term 命中该通路：命中为通路实色，未命中为浅灰。
- `Others` 用于完整性审计，不参与七通路总分，也不作为主图第八行。
- summary 显示每个集合中至少有一个对应直接 GO term 的蛋白数和比例。

## 审计结果

- BP：5,246 个蛋白–term 配对、1,713 个唯一 term、399 个蛋白。
- CC：2,757 个蛋白–term 配对、413 个唯一 term、398 个蛋白。
- MF：2,602 个蛋白–term 配对、659 个唯一 term、380 个蛋白。
- 103 个唯一 term 命中至少一条七通路规则。
- 8 个 term 同时命中两条通路。
- GO.db 3.20.0 缺少 20 个较新的 UniProt term；逐条检查后均不属于七条 DNA 修复通路，并保存在缺失 term 审计表中。

## 主要可复现输出

- `results/tables/go_term_pathway_scoring_30groups/direct_go_annotations_399proteins.csv`
- `results/tables/go_term_pathway_scoring_30groups/go_term_decision_audit_2785.csv`
- `results/tables/go_term_pathway_scoring_30groups/go_term_to_pathway_long.csv`
- `results/tables/go_term_pathway_scoring_30groups/multi_pathway_go_terms.csv`
- `results/tables/go_term_pathway_scoring_30groups/protein_pathway_direct_term_count_matrix.csv`
- `results/tables/go_term_pathway_scoring_30groups/protein_seven_pathway_binary_matrix.csv`
- `results/figures/five_set_pathway_matrix_go_term_30groups/`

旧人工评分 Excel、±1 状态和通路系数不参与本版本的计分或作图。
