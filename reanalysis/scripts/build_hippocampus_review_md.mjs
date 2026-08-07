import fs from "node:fs/promises";

const sourcePath =
  "/Users/gzy2520/Desktop/Research/kla/reanalysis/config/human_hippocampus_25_datasets.json";
const outputPath =
  "/Users/gzy2520/Desktop/Research/kla/reanalysis/reports/HUMAN_HIPPOCAMPUS_25_DATASET_REVIEW_TABLE.md";

const rows = JSON.parse(await fs.readFile(sourcePath, "utf8"));

if (rows.length !== 25 || new Set(rows.map((row) => row.pxd)).size !== 25) {
  throw new Error("Expected exactly 25 unique PXD records.");
}

function cell(value) {
  return String(value ?? "未报告")
    .replaceAll("|", "\\|")
    .replaceAll("\r\n", "<br>")
    .replaceAll("\n", "<br>");
}

function pxdLink(pxd) {
  return `[${pxd}](https://proteomecentral.proteomexchange.org/?pxid=${pxd})`;
}

function doiLinks(value) {
  if (!value || value === "未报告") return "未报告";
  return value
    .split(";")
    .map((doi) => doi.trim())
    .filter(Boolean)
    .map((doi) => `[${doi}](https://doi.org/${doi})`)
    .join("<br>");
}

function pmidLinks(value) {
  if (!value || value === "未报告") return "";
  return value
    .split(";")
    .map((pmid) => pmid.trim())
    .filter(Boolean)
    .map(
      (pmid) =>
        `[PMID ${pmid}](https://pubmed.ncbi.nlm.nih.gov/${pmid}/)`,
    )
    .join("<br>");
}

const headers = [
  "PXD",
  "数据集标题",
  "年份/发表",
  "真实材料与模型",
  "是否原生人海马",
  "人供者数",
  "性别",
  "年龄",
  "疾病/对照",
  "研究者额外操作",
  "样本处理与MS方法",
  "论文",
  "提交者/PI",
  "重复或关联关系",
  "Kla适用性",
  "等级与纳入建议",
  "证据来源/备注",
];

const lines = [
  "# ProteomeCentral 人海马 25 份数据审核大表",
  "",
  "生成日期：2026-07-30",
  "",
  "检索条件：[ProteomeCentral：hippocampus + Homo sapiens](https://proteomecentral.proteomexchange.org/ui?view=datasets&search=hippocampus&species=Homo%20sapiens)",
  "",
  "> 审核提醒：数据库物种字段不能直接等同于真实样本。当前 25 条中，14 条含原生人海马，1 条为含人尸检海马的混合模型，10 条并非原生人海马。普通蛋白组鉴定也不能直接作为 Kla 位点证据。",
  "",
  `| ${headers.join(" | ")} |`,
  `| ${headers.map(() => "---").join(" | ")} |`,
];

for (const row of rows) {
  const article = [
    doiLinks(row.doi),
    pmidLinks(row.pmid),
    row.first_author && row.first_author !== "未报告"
      ? `第一作者：${cell(row.first_author)}`
      : "",
  ]
    .filter(Boolean)
    .join("<br>");

  const model = [
    cell(row.actual_material),
    `模型：${cell(row.model_type)}`,
  ].join("<br>");

  const methods = [
    `处理：${cell(row.sample_processing)}`,
    `MS：${cell(row.ms_method)}`,
  ].join("<br>");

  const grade = [
    `等级：${cell(row.recommendation_grade)}`,
    cell(row.inclusion_advice),
  ].join("<br>");

  const evidence = [
    cell(row.evidence_level),
    `人口学来源：${cell(row.demographic_source)}`,
    cell(row.notes),
  ].join("<br>");

  const values = [
    pxdLink(row.pxd),
    cell(row.title),
    `${cell(row.paper_year)}<br>${cell(row.publication_status)}<br>数据库发布：${cell(row.announce_date)}`,
    model,
    cell(row.native_human_hippocampus),
    cell(row.human_donor_n),
    cell(row.sex),
    cell(row.age),
    cell(row.disease_group),
    cell(row.extra_operation),
    methods,
    article,
    cell(row.submitter_pi),
    cell(row.duplicate_related || "无已知重复登记"),
    cell(row.kla_suitability),
    grade,
    evidence,
  ];

  lines.push(`| ${values.join(" | ")} |`);
}

lines.push(
  "",
  "## 等级说明",
  "",
  "- A：正常或基线原生人海马。",
  "- B：原生人海马，但包含疾病、用药、保存方法或其他处理。",
  "- C：混合模型，包含人尸检海马但还混有体外或动物模型。",
  "- D：人源体外模型，如 iPSC、类器官或海马球体。",
  "- E：小鼠海马、其他人细胞或关键词误命中。",
  "- R：旧数据重分析或与其他 PXD 明确重叠。",
  "",
  "## 重点复核项",
  "",
  "- `PXD010543` 与 `PXD010544` 是同一批 4 TLE + 4 control 生物样本的两个 TMT 技术批次。",
  "- `PXD000395` 与 `PXD000950` 存在 control PE 原始文件重叠。",
  "- `PXD062981` 有相关论文，但论文 Data availability 未明确列出该 PXD，关联仍需复核。",
  "- `PXD050470` 是本批数据中唯一直接以正常人海马 Kla 图谱为主要目标的数据。",
  "",
);

await fs.writeFile(outputPath, lines.join("\n"), "utf8");
console.log(JSON.stringify({ outputPath, rows: rows.length }));
