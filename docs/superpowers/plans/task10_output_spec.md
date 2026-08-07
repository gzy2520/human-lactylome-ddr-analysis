# Task 10 输出规格（task10_output_spec.md）

依据：7 个被合并的 Node/Python 交付工具（只读，未修改）：
`build_cell_type_statistics_workbook.mjs`、`build_hippocampus_review_md.mjs`、
`build_reference_ddr_comparison_workbook.mjs`、`build_reference_proteome_selection_workbook.mjs`、
`build_venn_combined_workbook.mjs`、`create_bilingual_figure_legends_docx.js`、
`build_project_metadata.py`。本文件是 R 实现（`reanalysis/scripts/build_workbooks.R`）与内容一致性验证的依据。
所有项目根路径 = `<project_root>`（默认 `/Users/gzy2520/Desktop/Research/kla`）。

## 0. 总体实现约定

- xlsx 用 `writexl::write_xlsx(list(Sheet=df,...), path)`：**sheet 名与列名照抄**；数据从第 1 行表头、第 2 行起。
  .mjs 版有合并标题行、颜色、条件格式、冻结窗格、Excel 表格与公式（writexl 不支持），R 版以**直接计算的单元格值**代替公式，
  列内容一致（sheet 名/列名/行数/单元格值），版式细节有意省略（见报告）。
- 数字→字符串格式：JS `Number(s)→String()` 等价于"最短回读表示"。已实测当前全部输入字符串均为 JS 规范形
  （Python repr 验证 0 差异），R 版直接透传原字符串；整数串保持整数形式。
- csv 输出（.mjs 生成）均为：UTF-8 BOM + LF 行尾 + `QUOTE_MINIMAL` 引号规则 + **无结尾换行**。
- metadata 输出（.py 生成）均为：UTF-8 无 BOM + LF + `QUOTE_MINIMAL` + **每行（含最后一行）以 \n 结尾**。
- docx 用 officer 生成，正文文本（标题/段落/图注）与 .js 一致；页脚页码字段以 zip 后处理注入。

---

## 1. build_cell_type_statistics_workbook.mjs → cell_type_kla_ddr_statistics.xlsx

- 输入：`reanalysis/results/tables/cell_type_kla_ddr_statistics.csv`（BOM；4 列
  CellOrTissueType, TotalKlaProteins, KlaGoDdrProteins, KlaGoDdrFraction；10 数据行）。
- 输出：`reanalysis/results/tables/cell_type_kla_ddr_statistics.xlsx`（现存在，7/28 .mjs 版，5806B）。
- sheet「按细胞类型统计」（唯一 sheet）：
  - 表头（原第 4 行，R 版第 1 行）：`细胞系或组织 | 总 Kla 蛋白数 | 与 DDR 交集蛋白数 | DDR 交集 / 总数`
  - 数据（原第 5–14 行，10 行）：A=CellOrTissueType；B=Number(TotalKlaProteins)；C=Number(KlaGoDdrProteins)；
    D=原公式 `=IFERROR(C/B,0)` 的计算值（B==0 → 0），即 `ifelse(B==0, 0, C/B)`（IEEE 双精度，与 JS 缓存值逐位一致，已实测）。
  - 原版第 1–2 行为合并标题/说明（呈现层，writexl 不重现）。

## 2. build_hippocampus_review_md.mjs → HUMAN_HIPPOCAMPUS_25_DATASET_REVIEW_TABLE.md

- 输入：`reanalysis/config/human_hippocampus_25_datasets.json`（25 条、PXD 全唯一，断言）。
- 输出：`reanalysis/reports/HUMAN_HIPPOCAMPUS_25_DATASET_REVIEW_TABLE.md`（现存在，7/30 .mjs 版，32743B，**基线字节一致要求**）。
- 内容结构（逐行照抄）：
  1. `# ProteomeCentral 人海马 25 份数据审核大表`
  2. 空行；`生成日期：2026-07-30`（.mjs 硬编码，R 版同样硬编码以保字节一致）；空行
  3. `检索条件：[ProteomeCentral：hippocampus + Homo sapiens](https://proteomecentral.proteomexchange.org/ui?view=datasets&search=hippocampus&species=Homo%20sapiens)`；空行
  4. `> 审核提醒：数据库物种字段不能直接等同于真实样本。当前 25 条中，14 条含原生人海马，1 条为含人尸检海马的混合模型，10 条并非原生人海马。普通蛋白组鉴定也不能直接作为 Kla 位点证据。`；空行
  5. 表头行（17 列）：`PXD | 数据集标题 | 年份/发表 | 真实材料与模型 | 是否原生人海马 | 人供者数 | 性别 | 年龄 | 疾病/对照 | 研究者额外操作 | 样本处理与MS方法 | 论文 | 提交者/PI | 重复或关联关系 | Kla适用性 | 等级与纳入建议 | 证据来源/备注`
  6. 分隔行：17 个 `---` 以 ` | ` 连接
  7. 25 行数据行；每格按 .mjs 规则组装（`cell()`：`String(v ?? "未报告")` 后 `|`→`\|`、`\r\n`/`\n`→`<br>`；doiLinks/pmidLinks 按 `;` 拆分、trim、链接 `<br>` 连接；`article` = [doi, pmid, 第一作者] 过滤空后 `<br>` 连接；model/methods/grade/evidence 各段 `<br>` 连接）
  8. 尾部：空行、`## 等级说明`、空行、6 条 `- A：…/B：…/C：…/D：…/E：…/R：…`、空行、`## 重点复核项`、空行、
     4 条 `- \`PXD010543\` 与 \`PXD010544\` …`（逐字照抄 .mjs）、末尾元素为空串 → 文件以单个 `\n` 结尾、无 BOM。

## 3. build_reference_ddr_comparison_workbook.mjs

- 输入：`reanalysis/results/tables/cell_type_kla_vs_reference_ddr_statistics.csv`（10 行，断言）、
  `tall104_surrogate_ddr_sensitivity.csv`（4 行）、`cell_type_reference_control_information.csv`（10 行，断言，23 列）。
- 硬编码 `chineseInfo`：10 个键（Human hippocampus/HEK293T/HK-2/MCF10A/MCF7/HCT116/T-ALL/MDA-MB-468/T-47D/RKO），
  每键 13 个中文字段（name/matchType/subset/baseline/sampleCount/acquisition/search/ptm/rationale/caveat/countDetail/completeness/statistics）——逐字照抄 .mjs。
- 输出 A：`reanalysis/results/tables/cell_type_kla_vs_reference_ddr_statistics.xlsx`（现存在，7/30 版，18082B），4 sheets：
  - 「DDR占比对照」16 列：`细胞系或组织 | 参考PXD | 参考样本子集 | 参考蛋白数 | 参考DDR蛋白数 | 参考DDR占比 | Kla蛋白数 | Kla∩DDR蛋白数 | Kla DDR占比 | Kla-参考百分点差 | Kla/参考占比倍数 | 参考计数单位 | 登录号匹配数 | 基因符号辅助匹配数 | 蛋白数计数口径 | 描述性结论`
    数值：D,E,G,H,M,N=CSV 数字；F=`E/D`（D==0→0）；I=`H/G`（G==0→0）；J=`I-F`；K=`I/F`（F==0→0）；
    L=`UniProt基础登录号`（ReferenceCountUnit=="BaseAccession"）否则`论文发布的基因/蛋白特征`；
    P=`ifelse(abs(J)<0.002,"基本相当",ifelse(J>0,"Kla中DDR占比更高","参考组中DDR占比更高"))`。
  - 「对照选择信息」23 列（表头照抄 .mjs；值按 chineseInfo 映射 + controlInfo CSV 数字列；K=DDR占比 原样透传）。
  - 「TALL替代敏感性」4 列：`参考集合 | 蛋白数 | DDR蛋白数 | DDR占比`；D=`C/B`（B==0→0）。
  - 「方法与解释」2 列：`项目 | 说明`，7 行硬编码（逐字照抄 .mjs）。
- 输出 B/C（**基线字节一致要求**）：
  - `cell_type_kla_vs_reference_ddr_statistics_zh.csv`：BOM+LF、13 列
    `细胞或组织 | 参考PXD | 参考样本子集 | 参考蛋白数 | 参考DDR蛋白数 | 参考DDR占比 | Kla蛋白数 | Kla与DDR交集蛋白数 | Kla DDR占比 | Kla减参考百分点差 | Kla与参考占比倍数 | 选择理由 | 主要限制`
    （数值列 = CSV 原串透传——已实测全部为 JS 规范形）；无结尾换行。
  - `cell_type_reference_control_information_zh.csv`：BOM+LF、23 列（表头与值照抄 .mjs 的 controlZh）；无结尾换行。

## 4. build_reference_proteome_selection_workbook.mjs

- 输入：`reanalysis/config/reference_proteome_selection.json`（10 条，断言）、
  `reanalysis/results/tables/cell_type_kla_ddr_statistics.csv`（取 TotalKlaProteins/KlaGoDdrProteins/KlaGoDdrFraction）。
- 行构造：JSON 行 ∪ {total_kla_proteins=Number, kla_go_ddr_proteins=Number, kla_go_ddr_fraction=Number, decision_date="2026-07-30"}。
- 输出 A：`reanalysis/results/tables/cell_type_reference_proteome_selection.csv`（**基线字节一致**；BOM+LF+35 列；表头见 .mjs mainColumns 中文列名；无结尾换行）。
- 输出 B：`reference_proteome_pxd030304_sample_audit.csv`（**基线字节一致**；BOM+LF；6 列
  `当前类型 | PXD030304标签 | SIDM/项目标识 | 运行数 | 用途 | 说明`；9 行硬编码；无结尾换行）。
- 输出 C：`reanalysis/results/tables/cell_type_reference_proteome_selection.xlsx`（现存在，7/30 版，18754B），3 sheets：
  - 「10类对照选择」36 列（35 主列 + 辅助列 `QC辅助：PXD首次出现`）：
    表头=mainColumns 中文名；F 列=`Kla∩DDR/Kla`=E/D 计算值（D==0→0）；辅助列=参考PXD 首现行 →1，否则 0。
  - 「PXD030304样本核验」6 列，9 行硬编码（auditRows 逐字照抄）。
  - 「数据集与比较规则」9 列、29 行布局：
    行1 表头 `选择摘要 | 数量`（C–I 空）；行2–6 A 列五标签 + B 列计算值（B2=10、B3=9、B4=1、B5=7、B6=4；
    B3=Exact cell-line match(8)+Exact tissue match(1)，B6=参考PXD 去重首现数）；
    行7–8 空；行9 数据集表头（9 列）；行10–15 六条数据集行（硬编码 datasetRows）；
    行16–17 空；行18 `规则 | 执行说明` 表头；行19–29 十一条规则行（硬编码 ruleRows）。
    空行以 NA 单元格表示（readxl 读回为空）。

## 5. build_venn_combined_workbook.mjs → venn_combined_tables.xlsx

- 输出：`reanalysis/results/tables/venn_combined_tables.xlsx`（**当前不存在**，以本规格为准）。
- 4 个数据 sheet（内容 = CSV 原样，表头为 CSV 表头）：
  - `Kla_unique` ← `all_kla_three_groups_combined.csv`（21 列，3112 数据行）
  - `Kla_non_dedup` ← `all_kla_three_groups_combined_non_deduplicated.csv`（22 列，含 SourceCategory，4940 行）
  - `Kla_DDR_unique` ← `kla_go_ddr_three_groups_combined.csv`（21 列，275 行）
  - `Kla_DDR_non_dedup` ← `kla_go_ddr_three_groups_combined_non_deduplicated.csv`（22 列，452 行）
- `README` sheet（3 列，14 行布局）：
  行1 A=`Kla Venn combined protein tables`；行2 空；行3 表头 `Sheet | Rows | Meaning`；
  行4–7 四 sheet 名 + 行数（原 COUNTA 计算值：3112/4940/275/452）+ 含义文字（照抄 .mjs）；
  行8 空；行9 空；行10 A=`Source CSV files`；行11–14 四个 CSV 文件名；B/C 空。

## 6. create_bilingual_figure_legends_docx.js → Kla_Venn_figure_legends_bilingual.docx

- 输出：`reanalysis/results/Kla_Venn_figure_legends_bilingual.docx`（现存在，7/29 版，214064B，基线记录；docx 允许字节差异、文本一致）。
- 输入图片：`reanalysis/results/figures/all_kla_three_group_venn.png`（500×449 px）、`kla_go_ddr_three_group_venn.png`（478×466 px）。
- 正文结构（officer）：
  1. 标题段：`Bilingual Figure Legends / 双语图注`（居中、加粗、16pt、Arial/SimHei）
  2. 注释段：`注：Figure X/Y 与图 X/Y 为占位符，使用时请替换为稿件中的实际图号。`（居中、斜体、9pt、灰 #666666）
  3. 小节：`Figure X / 图 X`（加粗 12pt）
  4. 图 1（居中，5.2083×4.6771 in，alt="Venn diagram of all Kla proteins"）
  5. 图注 EN：`Figure X. Distribution of lysine-lactylated proteins among hippocampal tissue, immortalized models, and tumor cell lines. ` + 正文（10.5pt，标签加粗）
  6. 图注 ZH：`图 X. 人海马组织、永生化模型与肿瘤细胞系中赖氨酸乳酰化蛋白的分布。` + 正文
  7. 小节：`Figure Y / 图 Y`（**段前分页**）
  8. 图 2（居中，4.9791×4.8542 in，alt="Venn diagram of Kla proteins associated with DNA repair and damage responses"）
  9. 图注 EN：`Figure Y. Distribution of Kla proteins associated with DNA repair and DNA damage responses among hippocampal tissue, immortalized models, and tumor cell lines. ` + 正文
  10. 图注 ZH：`图 Y. 人海马组织、永生化模型与肿瘤细胞系中 DNA 修复和 DNA 损伤应答相关 Kla 蛋白的分布。` + 正文
- 图注正文（EN/ZH 各两段）逐字照抄 .js 常量（figure1EnglishBody/figure1ChineseBody/figure2EnglishBody/figure2ChineseBody）。
- 页脚：居中页码字段（`PAGE`，9pt，灰 #777777，Times New Roman/SimSun）；officer 无页脚创建 API →
  在 `print()` 后对 zip 做后处理（写 word/footer1.xml、补 [Content_Types].xml、document.xml.rels、sectPr footerReference）。
- 页面：Letter（11906×16838 twips）；边距 top/bottom 1134、left/right 1418 twips。
- 文档属性：creator=`Codex`，title=`Bilingual figure legends for Kla Venn diagrams`，
  description=`Publication-ready English and Chinese legends for two Kla Venn diagrams.`。

## 7. build_project_metadata.py（metadata 段）

- 输入：`reanalysis/config/datasets.csv`（utf-8-sig，10 行含表头，断言 PXD 集合 == EXPECTED_PXDS 九个）；
  `archive/migration_manifest_2026-07-21.csv`（utf-8-sig）。
- EXPECTED_PXDS = {PXD014870, PXD028488, PXD038880, PXD050470, PXD050906, PXD053474, PXD060185, PXD078013, PXD078736}。
- 输出（**均字节一致要求**；无 BOM，LF，QUOTE_MINIMAL，末行带 \n）：
  - `data/<PXD>/metadata/dataset_metadata.csv`：表头 = datasets.csv 全部列名（顺序照抄），单行 = 该 PXD 行（值原样透传，含引号转义）。
  - `data/<PXD>/metadata/file_inventory.csv`：7 列
    `PXD | Area | RelativePath | CurrentPath | FileName | Extension | SizeBytes`；
    递归扫描 `data/<PXD>/`（排除 .DS_Store、.pyc 后缀、该目录下 dataset_metadata.csv/file_inventory.csv 两个生成文件）；
    排序键 = posix 相对路径（Python `rglob` 后按 as_posix 排序，UTF-8 字节序 == 码点序，R radix 排序等价）；
    Area = 相对路径第一段（无嵌套 → `root`）；Extension = Python Path.suffix 语义（小写）。
  - `reanalysis/reports/project_file_inventory.csv`：9 个 PXD（按 PXD 升序）的 file_inventory 行拼接。
  - `archive/migration_manifest_reconstructed_2026-07-22.csv`：9 列
    `current_path | size_bytes | reconstructed_original_path | recorded_timestamp_utc | recorded_action | recorded_reason | provenance_evidence | reconstruction_confidence | notes`；
    扫描整个项目根（排除 .DS_Store、.pyc、输出自身）；provenance 链解析逻辑照译（destination→source 回溯、环检测、
    archive/2026-07-21_pre_restructure/ 与 archive/reanalysis_v1_2026-07-21/ 前缀推断、默认 unknown 记录）。
- 运行顺序：9 个 PXD 逐个（dataset_metadata → file_inventory），再项目 inventory，最后 reconstructed manifest。

---

## 8. 验证锚点（现磁盘状态）

- 基线：`archive/refactor_baseline_2026-08-06/output_sha256_baseline.csv`（362 文件）。
- 需字节一致：md、3 个 zh/selection csv、audit csv、project_file_inventory.csv、per-PXD metadata csv 与 file_inventory.csv。
- 内容一致（readxl 文本模式逐单元格比较）：3 个 xlsx。
- docx 文本一致：`unzip -p ... word/document.xml` 文本抽取对比。
- 已知漂移：`data/PXD050470/metadata/mqpar.xml`（2026-08-05 14:17 创建，已录入 8/5 的 per-PXD file_inventory.csv）
  未录入 7/22 生成的 project_file_inventory.csv（312 行）→ R 重跑后为 313 行，与基线 SHA 不同（预存在漂移，非 R 实现差异）。
- venn_combined_tables.xlsx 不在基线中（当前不存在），以本规格为准。
