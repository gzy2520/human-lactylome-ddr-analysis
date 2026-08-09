# DDR BaseAccession 分母修复

日期：2026-08-09

## 修复范围

普通全蛋白参照侧的 DDR 总蛋白数和 DDR 蛋白数现已统一使用映射后、去除
isoform 后缀并去重的人源 UniProt `BaseAccession`。

旧实现对大多数参照使用 BaseAccession，但 PXD010154 的正常肺和正常胎盘参照
仍以源 Ensembl protein ID 计数。源 ID 现仅用于来源追溯；一对多映射被展开，
未映射 ID 保留在审计中但不进入 accession-only 主分析分母。

## 实际变化

| Kla组 | 普通全蛋白参照 | 修复前总数 | 修复后总数 | 修复前DDR数 | 修复后DDR数 | 修复前占比 | 修复后占比 |
|---|---|---:|---:|---:|---:|---:|---:|
| PXD036307 / normal human lung | PXD010154 healthy lung | 14,022 | 12,358 | 668 | 535 | 4.7639% | 4.3292% |
| PXD063047 / normal pregnancy placenta | PXD010154 healthy placenta | 12,395 | 10,847 | 556 | 434 | 4.4857% | 4.0011% |

其余31条严格参照关系的 DDR 计数未改变；全部 Kla 蛋白数和 Kla-DDR 蛋白数
未改变。33组范围、9/2/9/13四分类、30条唯一普通全蛋白显示轴和两套调控因子
热图均未改变。

## 重建结果

- DDR中英文 PNG/PDF 正式图；内容完全相同的无语言后缀和
  `accession_only` 图像别名已归档；
- 37组和严格配对33组DDR统计表；
- DDR来源、映射及纳入排除审计；
- 四套四分类Euler/Venn图及精确membership/region表；
- 老师审阅CSV和XLSX；
- 最终SHA256清单。

四套Euler/Venn的精确membership和region counts与修复前一致，因为该流程原本
已经使用 `MappedBaseAccessions`。

旧版结果已复制到：

`archive/2026-08-09_pre_baseaccession_ddr_fix/`

归档文件及SHA256见该目录的 `MIGRATION_MANIFEST.csv`。

重复输出别名另行归档到：

`archive/2026-08-09_redundant_output_aliases/`
