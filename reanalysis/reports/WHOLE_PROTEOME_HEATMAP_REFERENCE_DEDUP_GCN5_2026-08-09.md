# 普通全蛋白调控因子热图去重与 GCN5 修正

## 老师意见

- HK-2、MCF7、HCT116 共享的普通全蛋白参照在热图中只显示一次。
- 其他相同细胞系的不同处理或研究相邻显示，但不得错误合并。
- Writer 中补入 GCN5。

## 修正结果

严格配对审计仍保留 33 个 Kla 样本组；显示层按细胞材料、参照 PXD、
定量来源文件、普通全蛋白特征数和参照蛋白数生成唯一参照键，最终显示
30 个普通全蛋白参照行：

- HK-2：PXD072220 对应 PXD058534 和 PXD078736。
- MCF7：PXD030304 对应 PXD014870 和 PXD060185。
- HCT116：PXD030304 对应 PXD028488 和 PXD053474。

HCT116/Roseburia 共培养使用 PXD066351 的独立普通全蛋白参照，因此保留为
同一 HCT116 簇中的另一行。两个人成纤维细胞研究同样相邻显示，但因参照文件
不同而不合并。

## GCN5

原 identifier 表中已经存在官方基因名 `KAT2A`，对应 UniProt BaseAccession
`Q92830`；GCN5 是同一蛋白的常用名。因此本次没有新增第二条蛋白记录，而是将
热图显示名改为 `GCN5 (KAT2A)`。定量匹配仍只使用 `Q92830`，GeneSymbol
不参与命中或缺失回退。

30 个唯一参照中有 7 个检出 Q92830，组内百分位范围为 0--68.996。

## 审计与验证

- `reanalysis/results/tables/kla_regulator_whole_proteome_heatmap_rows.csv`
- `reanalysis/results/tables/kla_regulator_whole_proteome_heatmap_display_long.csv`
- `reanalysis/results/figures/kla_regulator_whole_proteome_relative_intensity_heatmap.png`
- `reanalysis/results/figures/kla_regulator_whole_proteome_relative_intensity_heatmap.pdf`

自动测试验证 33 条严格配对未删除、30 个显示键不重复、只有上述三组共享参照、
Q92830 在 Writer 中每个参照只出现一次，以及所有命中均为
`UniProt_BaseAccession_only`。
