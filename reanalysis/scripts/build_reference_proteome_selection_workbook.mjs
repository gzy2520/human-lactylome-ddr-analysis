import fs from "node:fs/promises";
import { SpreadsheetFile, Workbook } from "@oai/artifact-tool";

const projectRoot = "/Users/gzy2520/Desktop/Research/kla";
const selectionPath =
  `${projectRoot}/reanalysis/config/reference_proteome_selection.json`;
const statsPath =
  `${projectRoot}/reanalysis/results/tables/cell_type_kla_ddr_statistics.csv`;
const outputDir = `${projectRoot}/reanalysis/results/tables`;
const previewDir = "/tmp/kla-reference-proteomes-work/previews";

const selection = JSON.parse(await fs.readFile(selectionPath, "utf8"));
const statsText = (await fs.readFile(statsPath, "utf8")).replace(/^\uFEFF/, "");
const statsLines = statsText.trim().split(/\r?\n/);
const statsHeaders = statsLines[0].split(",");
const stats = new Map(
  statsLines.slice(1).map((line) => {
    const values = line.split(",");
    const row = Object.fromEntries(
      statsHeaders.map((header, index) => [header, values[index]]),
    );
    return [row.CellOrTissueType, row];
  }),
);

if (selection.length !== 10) {
  throw new Error(`Expected 10 selected cell/tissue types, found ${selection.length}`);
}

const rows = selection.map((row) => {
  const current = stats.get(row.cell_type);
  if (!current) {
    throw new Error(`Missing current Kla statistics for ${row.cell_type}`);
  }
  return {
    ...row,
    total_kla_proteins: Number(current.TotalKlaProteins),
    kla_go_ddr_proteins: Number(current.KlaGoDdrProteins),
    kla_go_ddr_fraction: Number(current.KlaGoDdrFraction),
    decision_date: "2026-07-30",
  };
});

const mainColumns = [
  ["cell_type", "当前细胞/组织类型", 20],
  ["exact_identity", "准确模型身份", 38],
  ["kla_pxd", "当前Kla来源PXD", 26],
  ["total_kla_proteins", "Kla蛋白数", 13],
  ["kla_go_ddr_proteins", "Kla∩DDR蛋白数", 16],
  ["kla_go_ddr_fraction", "Kla∩DDR/Kla", 15],
  ["reference_pxd", "推荐常规蛋白组PXD", 19],
  ["reference_subset", "应使用的样本子集", 48],
  ["match_type", "匹配类型", 26],
  ["reference_year", "参考年份", 12],
  ["baseline_definition", "正常/基线定义", 44],
  ["selection_rationale", "选择理由", 60],
  ["reference_sample_count", "参考样本/采集数", 34],
  ["reference_protein_count_main", "对照蛋白数（主计数）", 18],
  ["reference_protein_count_detail", "蛋白数计数口径", 58],
  ["reference_protein_count_source", "蛋白数来源文件", 52],
  ["proteome_depth", "蛋白组深度", 42],
  ["acquisition", "采集方式/仪器", 34],
  ["search_quantification", "检索与定量", 40],
  ["ptm_enrichment", "是否PTM富集", 15],
  ["raw_data_status", "原始数据情况", 42],
  ["processed_data_status", "处理结果情况", 48],
  ["repository_completeness", "仓库与分析完整度", 48],
  ["recommended_file", "优先读取文件", 52],
  ["download_priority", "下载优先级", 44],
  ["suitability_grade", "适用等级", 12],
  ["use_for_detection_background", "可作检出背景", 18],
  ["use_for_statistical_differential", "可作统计差异", 34],
  ["main_caveat", "主要限制", 56],
  ["backup_reference", "备用数据", 48],
  ["dataset_url", "数据集URL", 48],
  ["paper_url", "论文URL", 44],
  ["processed_data_url", "处理数据URL", 48],
  ["selection_status", "选择状态", 34],
  ["decision_date", "核验日期", 14],
];
const helperColumn = columnLetter(mainColumns.length);

const auditRows = [
  ["HEK293T", "Control_HEK293T_lys", "Control_HEK293T", 401, "主要检出背景", "未处理裂解物过程控制；不能当作401个生物重复"],
  ["HEK293T", "Control_HEK293T_std_H002", "Control_HEK293T", 663, "仅用于技术稳定性检查", "标准QC，不进入主要蛋白集合"],
  ["MCF7", "MCF7", "SIDM00148", 6, "主要对照", "准确细胞系"],
  ["HCT116", "HCT-116", "SIDM00783", 6, "主要对照", "名称需映射 HCT116 -> HCT-116"],
  ["MDA-MB-468", "MDA-MB-468", "SIDM00628", 6, "主要对照", "准确细胞系"],
  ["T-47D", "T47D", "SIDM00097", 6, "主要对照", "名称需映射 T-47D -> T47D"],
  ["RKO", "RKO", "SIDM01090", 6, "主要对照", "准确细胞系"],
  ["T-ALL", "TALL-1", "SIDM00370", 6, "主要替代对照", "与TALL-104不是同一细胞系"],
  ["T-ALL", "Jurkat", "SIDM01016", 6, "敏感性替代对照", "与TALL-104不是同一细胞系"],
];

const datasetRows = [
  [
    "PXD030304",
    "HEK293T、MCF7、HCT116、T-ALL替代、MDA-MB-468、T-47D、RKO",
    "PXD030304为PARTIAL；论文配套Figshare提供最终矩阵和映射",
    "6,981次采集（仓库说明）；分析矩阵README列出6,864次运行",
    "6,692蛋白高置信矩阵；8,498蛋白敏感性矩阵",
    "mapping_file_averaged / mapping_file_replicates / protein_matrix_6692 / protein_matrix_8498",
    "先下载处理矩阵，不下载全部原始DIA",
    "https://proteomecentral.proteomexchange.org/?pxid=PXD030304",
    "https://doi.org/10.6084/m9.figshare.19345397",
  ],
  [
    "PXD043880",
    "Human hippocampus",
    "PARTIAL，但74个raw和74x2,092蛋白处理矩阵可用",
    "74名神经学正常供者",
    "2,092个蛋白特征",
    "13024_2023_650_MOESM1_ESM.xlsx / Source Data Proteins",
    "先使用补充矩阵；统一重检索时再下载raw",
    "https://proteomecentral.proteomexchange.org/?pxid=PXD043880",
    "https://doi.org/10.1186/s13024-023-00650-3",
  ],
  [
    "PXD072220",
    "HK-2",
    "PARTIAL，但9个raw和Spectronaut关键报告可用",
    "3 control + 3 copper + 3 TTM",
    "4,933蛋白；89-90%完整度",
    "HK-2_Spectronaut-report_PG_Quantity.txt",
    "只读取3个Control样本",
    "https://proteomecentral.proteomexchange.org/?pxid=PXD072220",
    "https://doi.org/10.1152/ajpcell.00311.2026",
  ],
  [
    "PXD002400",
    "MCF10A；MCF7备用",
    "PARTIAL；106个raw和一个261.6MB MaxQuant结果包",
    "MCF10A 10个IEF组分x2次进样=20个raw",
    "深度分级DDA蛋白组",
    "msms.zip",
    "先下载MaxQuant结果包",
    "https://proteomecentral.proteomexchange.org/?pxid=PXD002400",
    "https://www.ebi.ac.uk/pride/archive/projects/PXD002400",
  ],
  [
    "PXD027472",
    "HEK293T备用",
    "jPOST数据；whole-cell lysate和crude membrane两组",
    "样本数需在下载前复核",
    "全细胞与膜蛋白组",
    "whole-cell-lysate arm",
    "只使用全细胞裂解物，排除膜富集组",
    "https://proteomecentral.proteomexchange.org/?pxid=PXD027472",
    "https://doi.org/10.1016/j.mcpro.2022.100206",
  ],
  [
    "PXD028488",
    "TALL-104同细胞系技术敏感性检查",
    "现有Kla研究本身，不是独立外部基线",
    "TALL-104 enriched/non-enriched",
    "乳酸暴露后的非富集蛋白组",
    "TALL-Nonenrichment",
    "不能作为正常主对照；只作同细胞系敏感性检查",
    "https://proteomecentral.proteomexchange.org/?pxid=PXD028488",
    "https://doi.org/10.1038/s41592-022-01523-1",
  ],
];

const ruleRows = [
  ["正常的定义", "同一细胞系/组织、未处理或使用明确的 untreated/control 子集；癌细胞系的“正常”不代表非癌，只代表该细胞系的基线状态。"],
  ["无偏的定义", "不做Kla、磷酸化、乙酰化、泛素化、免疫沉淀、BioID或膜蛋白富集。允许为了提高深度进行普通肽段分级。"],
  ["主比较层级", "仅在蛋白层面比较。常规蛋白组没有Kla位点信息，不能用于判断某个赖氨酸位点是否未乳酸化。"],
  ["主背景集合", "每种模型的参考蛋白组中通过质量控制并被检出的BaseAccession集合。"],
  ["PXD030304主规则", "使用6,692蛋白矩阵作为主结果；8,498蛋白矩阵作为敏感性结果。"],
  ["Kla未检出解释", "常规蛋白组检测到、Kla表未出现的蛋白只能称为“未在当前Kla实验中鉴定”，不能称为“不乳酸化”。"],
  ["参考组未检出解释", "Kla蛋白不在参考蛋白组中可能来自低丰度、批次、仪器、搜索库或富集增敏，不应直接删除。"],
  ["UniProt规则", "继续使用去isoform后缀的BaseAccession优先匹配，GeneSymbol仅作辅助。"],
  ["TALL-104规则", "TALL-1和Jurkat结果必须分别输出；两者共同支持时才标记为T-ALL替代背景稳定检出。"],
  ["禁止的统计", "不能把技术进样、IEF组分或HEK293T QC运行数当作独立生物重复进行显著性检验。"],
  ["建议输出", "每个模型保存reference_detected、Kla∩reference、Kla-reference、DDR归一化比例和映射失败表。"],
];

function csvEscape(value) {
  const text = String(value ?? "");
  return /[",\r\n]/.test(text) ? `"${text.replaceAll('"', '""')}"` : text;
}

function columnLetter(index) {
  let value = index + 1;
  let result = "";
  while (value > 0) {
    const remainder = (value - 1) % 26;
    result = String.fromCharCode(65 + remainder) + result;
    value = Math.floor((value - 1) / 26);
  }
  return result;
}

const headers = mainColumns.map(([, header]) => header);
const csvRows = rows.map((row) =>
  mainColumns.map(([key]) => row[key] ?? ""),
);
const csvText = [headers, ...csvRows]
  .map((row) => row.map(csvEscape).join(","))
  .join("\n");

await fs.mkdir(outputDir, { recursive: true });
await fs.writeFile(
  `${outputDir}/cell_type_reference_proteome_selection.csv`,
  `\uFEFF${csvText}`,
  "utf8",
);

const auditCsv = [
  ["当前类型", "PXD030304标签", "SIDM/项目标识", "运行数", "用途", "说明"],
  ...auditRows,
]
  .map((row) => row.map(csvEscape).join(","))
  .join("\n");
await fs.writeFile(
  `${outputDir}/reference_proteome_pxd030304_sample_audit.csv`,
  `\uFEFF${auditCsv}`,
  "utf8",
);

const workbook = Workbook.create();

function addMainSheet() {
  const sheet = workbook.worksheets.add("10类对照选择");
  const values = [
    headers,
    ...rows.map((row) =>
      mainColumns.map(([key]) =>
        key === "kla_go_ddr_fraction" ? null : row[key] ?? "",
      ),
    ),
  ];
  const lastColumn = columnLetter(mainColumns.length - 1);
  sheet.getRange(`A1:${lastColumn}${values.length}`).values = values;
  sheet.getRange("F2").formulas = [["=E2/D2"]];
  sheet.getRange(`F2:F${values.length}`).fillDown();
  sheet.getRange(`F2:F${values.length}`).format.numberFormat = "0.0%";

  const used = sheet.getRange(`A1:${lastColumn}${values.length}`);
  used.format.wrapText = true;
  used.format.verticalAlignment = "top";
  used.format.font = { size: 10 };
  used.format.borders = {
    insideHorizontal: { style: "thin", color: "#D9E2F3" },
    bottom: { style: "thin", color: "#B4C6E7" },
  };

  const header = sheet.getRange(`A1:${lastColumn}1`);
  header.format.fill = "#1F4E78";
  header.format.font = { bold: true, color: "#FFFFFF", size: 10 };
  header.format.rowHeight = 34;
  header.format.verticalAlignment = "center";
  sheet.getRange(`A2:${lastColumn}${values.length}`).format.rowHeight = 72;

  mainColumns.forEach(([, , width], index) => {
    sheet
      .getRangeByIndexes(0, index, values.length, 1)
      .format.columnWidth = width;
  });

  const gradeIndex = mainColumns.findIndex(([key]) => key === "suitability_grade");
  const statusIndex = mainColumns.findIndex(([key]) => key === "selection_status");
  rows.forEach((row, index) => {
    const gradeFill =
      row.suitability_grade === "A"
        ? "#D9EAD3"
        : row.suitability_grade === "B"
          ? "#FFF2CC"
          : "#F4CCCC";
    sheet.getCell(index + 1, gradeIndex).format.fill = gradeFill;
    sheet.getCell(index + 1, gradeIndex).format.font = { bold: true };
    sheet.getCell(index + 1, gradeIndex).format.horizontalAlignment = "center";
    sheet.getCell(index + 1, statusIndex).format.fill =
      row.selection_status.startsWith("Recommended") ? "#E2F0D9" : "#FCE4D6";
  });

  const table = sheet.tables.add(
    `A1:${lastColumn}${values.length}`,
    true,
    "ReferenceProteomeSelectionTable",
  );
  table.style = "TableStyleMedium2";
  table.showFilterButton = true;
  table.showBandedRows = true;
  sheet.freezePanes.freezeRows(1);
  sheet.freezePanes.freezeColumns(3);
  sheet.showGridLines = false;

  sheet.getRange(`${helperColumn}1`).values = [["QC辅助：PXD首次出现"]];
  sheet
    .getRange(`${helperColumn}2`)
    .formulas = [["=IF(COUNTIF($G$2:G2,G2)=1,1,0)"]];
  sheet.getRange(`${helperColumn}2:${helperColumn}${values.length}`).fillDown();
  sheet.getRange(`${helperColumn}1:${helperColumn}${values.length}`).format = {
    fill: "#E7E6E6",
    font: { color: "#666666", size: 9 },
    numberFormat: "0",
  };
  sheet.getRange(`${helperColumn}1`).format.font = {
    bold: true,
    color: "#666666",
    size: 9,
  };
  sheet
    .getRange(`${helperColumn}1:${helperColumn}${values.length}`)
    .format.columnWidth = 18;
}

function addAuditSheet() {
  const sheet = workbook.worksheets.add("PXD030304样本核验");
  const values = [
    ["当前类型", "PXD030304标签", "SIDM/项目标识", "运行数", "用途", "说明"],
    ...auditRows,
  ];
  sheet.getRange(`A1:F${values.length}`).values = values;
  sheet.getRange(`A1:F${values.length}`).format.wrapText = true;
  sheet.getRange(`A1:F${values.length}`).format.verticalAlignment = "top";
  sheet.getRange("A1:F1").format.fill = "#1F4E78";
  sheet.getRange("A1:F1").format.font = { bold: true, color: "#FFFFFF" };
  sheet.getRange("A2:F10").format.rowHeight = 44;
  [18, 30, 22, 12, 24, 52].forEach((width, index) => {
    sheet.getRangeByIndexes(0, index, values.length, 1).format.columnWidth = width;
  });
  sheet.getRange("D2:D10").format.numberFormat = "0";
  const table = sheet.tables.add(
    `A1:F${values.length}`,
    true,
    "PXD030304AuditTable",
  );
  table.style = "TableStyleMedium2";
  sheet.freezePanes.freezeRows(1);
  sheet.showGridLines = false;
}

function addRulesSheet() {
  const sheet = workbook.worksheets.add("数据集与比较规则");
  sheet.getRange("A1:B1").values = [["选择摘要", "数量"]];
  sheet.getRange("A1:B1").format.fill = "#1F4E78";
  sheet.getRange("A1:B1").format.font = { bold: true, color: "#FFFFFF" };
  sheet.getRange("A2:A6").values = [
    ["当前细胞/组织类型"],
    ["精确匹配"],
    ["替代匹配"],
    ["A级"],
    ["主要参考PXD"],
  ];
  sheet.getRange("B2:B6").formulas = [
    ["=COUNTA('10类对照选择'!A2:A11)"],
    ["=COUNTIF('10类对照选择'!I2:I11,\"Exact cell-line match\")+COUNTIF('10类对照选择'!I2:I11,\"Exact tissue match\")"],
    ["=COUNTIF('10类对照选择'!I2:I11,\"Disease-matched surrogate, not an exact cell-line match\")"],
    ["=COUNTIF('10类对照选择'!U2:U11,\"A\")"],
    [`=SUM('10类对照选择'!${helperColumn}2:${helperColumn}11)`],
  ];
  sheet.getRange("B2:B6").format.numberFormat = "0";

  const datasetHeaderRow = 9;
  const datasetHeaders = [
    "PXD",
    "覆盖模型",
    "完整度说明",
    "样本/采集规模",
    "蛋白组深度",
    "优先文件",
    "下载策略",
    "数据集URL",
    "处理数据/论文URL",
  ];
  sheet.getRange(`A${datasetHeaderRow}:I${datasetHeaderRow}`).values = [
    datasetHeaders,
  ];
  sheet
    .getRange(`A${datasetHeaderRow}:I${datasetHeaderRow}`)
    .format.fill = "#5B9BD5";
  sheet
    .getRange(`A${datasetHeaderRow}:I${datasetHeaderRow}`)
    .format.font = { bold: true, color: "#FFFFFF" };
  sheet.getRange(`A${datasetHeaderRow + 1}:I${datasetHeaderRow + datasetRows.length}`).values =
    datasetRows;

  const ruleHeaderRow = datasetHeaderRow + datasetRows.length + 3;
  sheet.getRange(`A${ruleHeaderRow}:B${ruleHeaderRow}`).values = [
    ["规则", "执行说明"],
  ];
  sheet.getRange(`A${ruleHeaderRow}:B${ruleHeaderRow}`).format.fill = "#70AD47";
  sheet.getRange(`A${ruleHeaderRow}:B${ruleHeaderRow}`).format.font = {
    bold: true,
    color: "#FFFFFF",
  };
  sheet.getRange(
    `A${ruleHeaderRow + 1}:B${ruleHeaderRow + ruleRows.length}`,
  ).values = ruleRows;

  const lastRow = ruleHeaderRow + ruleRows.length;
  sheet.getRange(`A1:I${lastRow}`).format.wrapText = true;
  sheet.getRange(`A1:I${lastRow}`).format.verticalAlignment = "top";
  sheet.getRange(`A2:I${lastRow}`).format.borders = {
    insideHorizontal: { style: "thin", color: "#D9E2F3" },
    bottom: { style: "thin", color: "#B4C6E7" },
  };
  [20, 48, 46, 34, 38, 50, 42, 48, 48].forEach((width, index) => {
    sheet.getRangeByIndexes(0, index, lastRow, 1).format.columnWidth = width;
  });
  sheet
    .getRange(`A${datasetHeaderRow + 1}:I${datasetHeaderRow + datasetRows.length}`)
    .format.rowHeight = 64;
  sheet
    .getRange(`A${ruleHeaderRow + 1}:B${lastRow}`)
    .format.rowHeight = 54;
  sheet.freezePanes.freezeRows(1);
  sheet.showGridLines = false;
}

addMainSheet();
addAuditSheet();
addRulesSheet();

const formulaErrors = await workbook.inspect({
  kind: "match",
  searchTerm: "#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A",
  options: { useRegex: true, maxResults: 100 },
  summary: "final formula error scan",
});

const checks = [];
for (const [sheetName, range] of [
  ["10类对照选择", "A1:J11"],
  ["10类对照选择", "K1:T11"],
  ["10类对照选择", "U1:AD11"],
  ["10类对照选择", "AE1:AI11"],
  ["PXD030304样本核验", "A1:F10"],
  ["数据集与比较规则", "A1:I16"],
  ["数据集与比较规则", "A18:B29"],
]) {
  const inspected = await workbook.inspect({
    kind: "table",
    range: `${sheetName}!${range}`,
    include: "values,formulas",
    tableMaxRows: 16,
    tableMaxCols: 10,
    maxChars: 6000,
  });
  checks.push({ sheetName, range, ndjson: inspected.ndjson });
}

await fs.mkdir(previewDir, { recursive: true });
for (const [fileName, sheetName, range] of [
  ["selection_left.png", "10类对照选择", "A1:J11"],
  ["selection_middle.png", "10类对照选择", "K1:T11"],
  ["selection_right.png", "10类对照选择", "U1:AD11"],
  ["selection_end.png", "10类对照选择", "AE1:AI11"],
  ["pxd030304_audit.png", "PXD030304样本核验", "A1:F10"],
  ["datasets.png", "数据集与比较规则", "A1:I16"],
  ["rules.png", "数据集与比较规则", "A18:B29"],
]) {
  const preview = await workbook.render({
    sheetName,
    range,
    scale: 1,
    format: "png",
  });
  await fs.writeFile(
    `${previewDir}/${fileName}`,
    new Uint8Array(await preview.arrayBuffer()),
  );
}

const exported = await SpreadsheetFile.exportXlsx(workbook);
await exported.save(
  `${outputDir}/cell_type_reference_proteome_selection.xlsx`,
);

await fs.writeFile(
  "/tmp/kla-reference-proteomes-work/verification.json",
  JSON.stringify(
    {
      rows: rows.length,
      cellTypes: rows.map((row) => row.cell_type),
      referencePxd: [...new Set(rows.map((row) => row.reference_pxd))],
      formulaErrors: formulaErrors.ndjson,
      checks,
    },
    null,
    2,
  ),
);

console.log(
  JSON.stringify({
    rows: rows.length,
    exactMatches: rows.filter((row) => row.match_type.startsWith("Exact")).length,
    surrogateMatches: rows.filter((row) =>
      row.match_type.startsWith("Disease-matched"),
    ).length,
    output: `${outputDir}/cell_type_reference_proteome_selection.xlsx`,
  }),
);
