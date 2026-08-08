# Kla 与普通全蛋白参照老师审阅表

## 交付内容

- `reanalysis/results/tables/kla_and_reference_teacher_review_zh.xlsx`
- `reanalysis/results/tables/kla_and_reference_teacher_review_zh.csv`
- `reanalysis/results/tables/kla_and_reference_teacher_review_focus_zh.csv`

XLSX 包含“说明与汇总”“老师审阅总表”“需重点确认”三页。每行对应
一个最终 Kla 样本组，并同时列出 Kla 文件、普通全蛋白参照文件、实际
读取样本/工作表/列、两侧蛋白与 DDR 数量、匹配限制和纳入决定。

## 最终核对

- 最终范围为 37 个唯一的 `PXD + 样本组`。
- 33 组具有严格生物材料身份/取材粒度匹配且可定量的普通全蛋白参照。
- 4 组因无完全匹配且有逐蛋白定量强度的参照而排除普通全蛋白侧：
  PXD062720 膀胱癌/EPI、PXD063047 重度子痫前期胎盘、
  PXD064038 MEC/NEC 食管鳞癌、PXD075014 AC16 对照/低氧。
- 37 组 Kla 证据路径全部存在；33 组已纳入参照路径全部存在。
- Kla 与普通全蛋白的 DDR 匹配均使用 UniProt BaseAccession；
  GeneSymbol 回退数为 0。
- HCC 与癌旁肝虽然共用 PXD065775 文件，但分别读取 CISs 和 ANTs。
- 增生性瘢痕与癌旁皮肤分别读取 HSP1-HSP4 和 NSP1-NSP4。
- BPH 与前列腺癌分别读取 NAT1-NAT5 和 PCa1-PCa5。
- 海马使用 PXD050470 同研究 Table S4 的 H072/H081/H0187，
  为 6,082 个普通全蛋白、219 个 DDR 蛋白。
- 两个 HK-2 组均使用 PXD072220 的 amostra1/amostra3/amostra4
  `PG.Log2Quantity`；已把原先过于笼统的配置文字修正为实际读取列。

## 尚需老师决定的限制

当前没有发现会推翻已有计数的代码错误，但“材料匹配”不能解释为所有
实验状态完全一致：

- 13 组只使用同一材料的基线参照，不代表 Kla 处理或时间点。
- 5 组使用独立队列的同一组织/疾病表型，不是同一供体。
- 1 组为同一细胞系的独立实验。
- 1 组同研究存在对应状态，但普通全蛋白定量方法仍需确认。
- 18 行的配对元数据尚未记录或核实论文 DOI；PXD 链接均已保留，
  不影响当前蛋白计数。

如果老师进一步要求“同一供体且同一处理/时间点”，应对上述基线参照
和独立队列另做更严格的排除或敏感性分析，不能把它们静默称为完全匹配。

## 验证

8 个 R 测试和 19 个 Python 测试全部通过。工作簿三页均完成渲染检查，
未发现公式错误、截断的关键表头或空白工作表。
