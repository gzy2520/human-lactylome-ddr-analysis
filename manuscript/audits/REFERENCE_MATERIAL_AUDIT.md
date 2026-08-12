# 普通全蛋白参照材料粒度严格审计

更新日期：2026-08-08

## 结论

- 最终 Kla 定量范围仍为 37 个样本组。
- 其中 33 组具有材料身份和取材粒度严格匹配的普通全蛋白参照。
- 4 组仍因没有完全匹配且含逐蛋白强度的普通全蛋白数据而排除。
- “材料身份匹配”和“实验处理状态匹配”已拆开记录。材料身份严格匹配是主分析纳入条件；只提供基线参照的处理状态不会再被写成完全匹配。

逐行证据见：

- `results/tables/strict_reference_material_identity_audit.csv`
- `results/tables/strict_reference_material_identity_audit_zh.csv`

## HCC 与癌旁肝

`PXD075377 / HCC` 和 `PXD075377 / adjacent liver` 的确共用
`PXD065775` 工作簿，但没有使用同一批蛋白强度：

- HCC 只读取 `CISs` 工作表的 `Non-rec1-4` 和 `Rec1-4`。
- 癌旁肝只读取 `ANTs` 工作表的 `Non-rec1-4` 和 `Rec1-4`。
- `DNTs` 远端正常组织不进入上述任一行。

因此这是“一份工作簿内的两个独立组织子集”，不是用统一肝组织同时替代 HCC 和癌旁肝。局限是普通全蛋白参照来自独立患者队列，并非 Kla 研究的同一患者。

## 海马修正

旧流程把未指定亚区的整个人海马 Kla 数据配到独立研究的 CA1 亚区
`PXD043880`。按照“相近不等于相同”的标准，这不合格。

现已改用 Kla 论文同研究的 Table S4：

- 文件：`data/PXD050470/supplementary/prca2331-sup-0006-tables4.xlsx`
- 样本：`H072`、`H081`、`H0187`
- 蛋白：6,082 个唯一有效 UniProt BaseAccession
- DDR：219 个
- 强度：三份同一海马样本的普通全蛋白相对强度
- 匹配：直接按 BaseAccession，不使用 GeneSymbol 转换

旧 CA1 图表和映射审计已保存到
`archive/2026-08-08_pre_material_granularity_fix`。

## 其他共用参照

- `PXD046800`：瘢痕只读 `HSP1-HSP4`；邻近皮肤只读 `NSP1-NSP4`。
- `PXD066054`：BPH 只读 `NAT1-NAT5`；前列腺癌只读 `PCa1-PCa5`。
- `PXD069969`：GSC 只读六个命名 GSC 模型；NSC 只读 `ENSA/HMP1`。
- `PXD030304`：共享的是细胞系图谱文件，但每个细胞系读取独立命名行；同一 MCF7 或 HCT116 被不同 Kla 研究复用时，属于同材料基线复用。
- `PXD072220`：两个 HK-2 Kla 研究都使用同一组未处理 HK-2 基线。细胞身份匹配，但预处理或甘露醇状态未被参照覆盖，已标记为 `baseline_reference_only`。
- `PXD073311`：HUVEC 参照只用 A0h 基线；材料匹配，但 Pg 感染状态未进入普通全蛋白参照。

## 图表影响

- Kla Venn：8,719 个蛋白，修改前后集合完全相同。
- Kla-DDR Venn：510 个蛋白，修改前后集合完全相同。
- 普通全蛋白 Venn 并集：24,566 变为 24,735，新增 175，移除 6。
- 普通全蛋白 DDR Venn 并集：764 变为 765，新增 1，移除 0。
- 正常/非肿瘤组织普通全蛋白集合：18,137 变为 18,468。
- 正常/非肿瘤组织普通全蛋白 DDR 集合：646 变为 649。
- 海马 Kla 蛋白和 Kla-DDR 数保持 853 和 29，不受参照更换影响。

## 判读限制

当前主分析要求材料身份和取材粒度完全一致，但允许独立队列或基线参照，并将实验状态差异单独披露。若老师进一步要求“同一供者、同一处理、同一时间点”也必须匹配，则审计表中 `baseline_reference_only` 和
`independent_cohort_same_phenotype` 行还需要再次排除或另作敏感性分析。
