import fs from "node:fs/promises";
import path from "node:path";
import { SpreadsheetFile, Workbook } from "@oai/artifact-tool";

const projectRoot = process.argv[2] || "/Users/gzy2520/Desktop/Research/kla";
const tableRoot = path.join(projectRoot, "reanalysis/results/tables");
const previewRoot = "/tmp/kla-37group-workbook-previews";

function parseCsv(text) {
  const rows = [];
  let row = [];
  let field = "";
  let quoted = false;
  const source = text.replace(/^\uFEFF/, "");
  for (let i = 0; i < source.length; i += 1) {
    const char = source[i];
    if (quoted) {
      if (char === '"' && source[i + 1] === '"') {
        field += '"';
        i += 1;
      } else if (char === '"') {
        quoted = false;
      } else {
        field += char;
      }
    } else if (char === '"') {
      quoted = true;
    } else if (char === ",") {
      row.push(field);
      field = "";
    } else if (char === "\n") {
      row.push(field.replace(/\r$/, ""));
      rows.push(row);
      row = [];
      field = "";
    } else {
      field += char;
    }
  }
  if (field.length || row.length) {
    row.push(field.replace(/\r$/, ""));
    rows.push(row);
  }
  return rows.filter((values) => values.some((value) => value !== ""));
}

async function readCsv(relativePath) {
  return parseCsv(await fs.readFile(path.join(projectRoot, relativePath), "utf8"));
}

function columnName(index) {
  let value = index + 1;
  let result = "";
  while (value > 0) {
    const remainder = (value - 1) % 26;
    result = String.fromCharCode(65 + remainder) + result;
    value = Math.floor((value - 1) / 26);
  }
  return result;
}

function typedValue(value) {
  if (value === "") return null;
  if (/^-?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?$/.test(value)) {
    const number = Number(value);
    if (Number.isFinite(number)) return number;
  }
  if (value === "TRUE") return true;
  if (value === "FALSE") return false;
  return value;
}

function addDataSheet(workbook, name, csvRows, options = {}) {
  const sheet = workbook.worksheets.add(name);
  sheet.showGridLines = false;
  const rows = csvRows.map((row, rowIndex) =>
    row.map((value) => (rowIndex === 0 ? value : typedValue(value))),
  );
  const rowCount = rows.length;
  const columnCount = Math.max(...rows.map((row) => row.length));
  const normalizedRows = rows.map((row) => [
    ...row,
    ...Array(columnCount - row.length).fill(null),
  ]);
  const lastColumn = columnName(columnCount - 1);
  sheet.getRange(`A1:${lastColumn}${rowCount}`).values = normalizedRows;
  sheet.getRange(`A1:${lastColumn}1`).format = {
    fill: options.headerFill || "#A64B2A",
    font: { bold: true, color: "#FFFFFF", size: 10 },
    horizontalAlignment: "center",
    verticalAlignment: "center",
    wrapText: true,
    borders: { preset: "outside", style: "thin", color: "#7A351E" },
  };
  sheet.getRange(`A2:${lastColumn}${rowCount}`).format = {
    font: { color: "#2B2B2B", size: 9 },
    verticalAlignment: "top",
    borders: {
      insideHorizontal: { style: "thin", color: "#E5E1DD" },
    },
  };
  sheet.getRange(`A1:${lastColumn}${Math.min(rowCount, 60)}`).format.wrapText = true;
  sheet.freezePanes.freezeRows(1);
  if (options.freezeColumns) sheet.freezePanes.freezeColumns(options.freezeColumns);

  for (let col = 0; col < columnCount; col += 1) {
    const values = normalizedRows
      .slice(0, Math.min(rowCount, 200))
      .map((row) => String(row[col] ?? ""));
    const maxLength = Math.max(...values.map((value) => value.length), 8);
    const width = Math.min(Math.max(maxLength + 2, 10), options.maxWidth || 32);
    sheet.getRange(`${columnName(col)}:${columnName(col)}`).format.columnWidth = width;
  }
  sheet.getRange(`A1:${lastColumn}1`).format.rowHeight = 34;
  return { sheet, rowCount, columnCount, lastColumn };
}

function addReadme(workbook, title, lines) {
  const sheet = workbook.worksheets.add("说明");
  sheet.showGridLines = false;
  sheet.getRange("A1:H1").merge();
  sheet.getRange("A1").values = [[title]];
  sheet.getRange("A1:H1").format = {
    fill: "#8F2D1D",
    font: { bold: true, color: "#FFFFFF", size: 16 },
    horizontalAlignment: "left",
    verticalAlignment: "center",
  };
  sheet.getRange("A1:H1").format.rowHeight = 32;
  const rows = lines.map((line) => [line]);
  sheet.getRange(`A3:H${rows.length + 2}`).merge(true);
  sheet.getRange(`A3:A${rows.length + 2}`).values = rows;
  sheet.getRange(`A3:H${rows.length + 2}`).format = {
    fill: "#FFF8F0",
    font: { color: "#3A2B25", size: 10 },
    wrapText: true,
    verticalAlignment: "top",
  };
  sheet.getRange("A:H").format.columnWidth = 16;
  sheet.getRange(`A3:H${rows.length + 2}`).format.rowHeight = 28;
  return sheet;
}

async function verifyAndExport(workbook, outputPath, previewPrefix) {
  const sheets = workbook.worksheets.items;
  const overview = await workbook.inspect({
    kind: "sheet",
    include: "id,name",
    maxChars: 5000,
  });
  console.log(overview.ndjson);
  const errors = await workbook.inspect({
    kind: "match",
    searchTerm: "#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A",
    options: { useRegex: true, maxResults: 100 },
    summary: `${previewPrefix} formula error scan`,
  });
  console.log(errors.ndjson);

  await fs.mkdir(previewRoot, { recursive: true });
  for (const sheet of sheets) {
    const used = sheet.getUsedRange(true);
    if (!used) continue;
    const rowCount = sheet.name === "说明"
      ? 8
      : Math.min(used.rowCount, 24);
    const columnCount = sheet.name === "说明"
      ? 8
      : Math.min(used.columnCount, 14);
    const range = `A1:${columnName(columnCount - 1)}${rowCount}`;
    const preview = await workbook.render({
      sheetName: sheet.name,
      range,
      scale: 1,
      format: "png",
    });
    const safeName = sheet.name.replace(/[^\p{L}\p{N}_-]+/gu, "_");
    await fs.writeFile(
      path.join(previewRoot, `${previewPrefix}_${safeName}.png`),
      new Uint8Array(await preview.arrayBuffer()),
    );
  }

  const output = await SpreadsheetFile.exportXlsx(workbook);
  await output.save(outputPath);
  await fs.rm(`${outputPath}.inspect.ndjson`, { force: true });
}

async function buildDdrWorkbook() {
  const stats = await readCsv(
    "reanalysis/results/tables/cell_type_kla_vs_reference_ddr_statistics_accession_only.csv",
  );
  const statsZh = await readCsv(
    "reanalysis/results/tables/cell_type_kla_vs_reference_ddr_statistics_accession_only_zh.csv",
  );
  const audit = await readCsv(
    "reanalysis/results/tables/cell_type_kla_vs_reference_ddr_accession_only_audit.csv",
  );
  const materialAudit = await readCsv(
    "reanalysis/results/tables/strict_reference_material_identity_audit_zh.csv",
  );
  const plotRows = await readCsv(
    "reanalysis/results/tables/cell_type_kla_vs_reference_ddr_plot_rows.csv",
  );
  const workbook = Workbook.create();
  addReadme(workbook, "37组 Kla 与33组严格普通全蛋白参照 DDR 占比统计", [
    "Kla最终范围为定量可用的37个 PXD+样本组；PXD037371 三组已排除。",
    "其中33组具有生物材料匹配且可审计的逐蛋白强度参照，进入柱状图和普通全蛋白分析；其余4组只保留Kla结果和排除原因。",
    "Kla 与普通全蛋白的 DDR 交集均按去除 isoform 后缀的 UniProt BaseAccession 计算；GeneSymbol 回退数为0。",
    "HUVEC 普通全蛋白参照为 PXD073311 同研究 A0h_1/A0h_2/A0h_3 PG矩阵；A6h不进入基线。",
    "两个 HK-2 组均使用 PXD072220 的 amostra1/amostra3/amostra4 未处理对照。",
    "人海马使用 PXD050470 同研究 Table S4 的 H072/H081/H0187，不再使用独立 CA1 参照。",
  ]);
  const zhSheet = addDataSheet(workbook, "统计_中文", statsZh, {
    freezeColumns: 2,
    maxWidth: 42,
    headerFill: "#A64B2A",
  });
  const zhLastRow = zhSheet.rowCount;
  zhSheet.sheet.getRange(`I2:I${zhLastRow}`).formulas = Array.from(
    { length: zhLastRow - 1 },
    (_, index) => [`=IFERROR(H${index + 2}/G${index + 2},0)`],
  );
  zhSheet.sheet.getRange(`N2:N${zhLastRow}`).formulas = Array.from(
    { length: zhLastRow - 1 },
    (_, index) => [`=IFERROR(M${index + 2}/L${index + 2},0)`],
  );
  zhSheet.sheet.getRange(`O2:O${zhLastRow}`).formulas = Array.from(
    { length: zhLastRow - 1 },
    (_, index) => [`=I${index + 2}-N${index + 2}`],
  );
  zhSheet.sheet.getRange(`I2:I${zhLastRow}`).format.numberFormat = "0.0%";
  zhSheet.sheet.getRange(`N2:O${zhLastRow}`).format.numberFormat = "0.0%";

  const enSheet = addDataSheet(workbook, "Statistics_EN", stats, {
    freezeColumns: 2,
    maxWidth: 42,
    headerFill: "#B5651D",
  });
  const enLastRow = enSheet.rowCount;
  enSheet.sheet.getRange(`H2:H${enLastRow}`).formulas = Array.from(
    { length: enLastRow - 1 },
    (_, index) => [`=IFERROR(G${index + 2}/F${index + 2},0)`],
  );
  enSheet.sheet.getRange(`M2:M${enLastRow}`).formulas = Array.from(
    { length: enLastRow - 1 },
    (_, index) => [`=IFERROR(L${index + 2}/K${index + 2},0)`],
  );
  enSheet.sheet.getRange(`N2:N${enLastRow}`).formulas = Array.from(
    { length: enLastRow - 1 },
    (_, index) => [`=H${index + 2}-M${index + 2}`],
  );
  enSheet.sheet.getRange(`H2:H${enLastRow}`).format.numberFormat = "0.0%";
  enSheet.sheet.getRange(`M2:N${enLastRow}`).format.numberFormat = "0.0%";
  addDataSheet(workbook, "纳入审计", audit, {
    freezeColumns: 2,
    maxWidth: 38,
    headerFill: "#8F2D1D",
  });
  addDataSheet(workbook, "材料粒度审计", materialAudit, {
    freezeColumns: 2,
    maxWidth: 46,
    headerFill: "#6F3D2E",
  });
  addDataSheet(workbook, "绘图显示去重", plotRows, {
    freezeColumns: 4,
    maxWidth: 46,
    headerFill: "#7A4A35",
  });
  await verifyAndExport(
    workbook,
    path.join(tableRoot, "cell_type_kla_vs_reference_ddr_statistics.xlsx"),
    "ddr37",
  );
}

async function buildPairingWorkbook() {
  const pairing = await readCsv(
    "reanalysis/results/tables/lactylome_and_reference_proteome_pairing_zh.csv",
  );
  const compact = await readCsv(
    "reanalysis/results/tables/lactylome_group_two_reference_columns_complete_zh.csv",
  );
  const summary = await readCsv(
    "reanalysis/results/tables/lactylome_reference_pairing_summary_zh.csv",
  );
  const materialAudit = await readCsv(
    "reanalysis/results/tables/strict_reference_material_identity_audit_zh.csv",
  );
  const workbook = Workbook.create();
  addReadme(workbook, "乳酸化数据与普通全蛋白参照配对", [
    "本工作簿记录乳酸化证据、普通全蛋白参照、健康组织基线、蛋白数和匹配限制。",
    "普通全蛋白参照必须是非Kla富集数据；不能用乳酸化富集强度替代。",
    "最终37组Kla中，33组进入严格普通全蛋白配对分析；4组因没有完全匹配且可定量的参照而排除普通全蛋白分析。",
    "PXD073311 HUVEC 已更换为同研究 A0h 普通全蛋白矩阵，唯一 BaseAccession 数为7,794。",
    "PXD050470 海马使用同研究同三份样本的 Table S4，唯一 BaseAccession 数为6,082。",
    "共享同一PXD或文件的组织/细胞已逐行记录实际读取子集；材料身份与实验状态分别审计。",
    "分析身份键为稳定蛋白ID；GeneSymbol只用于显示和人工审计。",
  ]);
  addDataSheet(workbook, "配对明细", pairing, {
    freezeColumns: 3,
    maxWidth: 46,
    headerFill: "#A64B2A",
  });
  addDataSheet(workbook, "两列精简表", compact, {
    freezeColumns: 3,
    maxWidth: 46,
    headerFill: "#B5651D",
  });
  addDataSheet(workbook, "配对汇总", summary, {
    freezeColumns: 1,
    maxWidth: 42,
    headerFill: "#8F2D1D",
  });
  addDataSheet(workbook, "材料粒度审计", materialAudit, {
    freezeColumns: 2,
    maxWidth: 46,
    headerFill: "#6F3D2E",
  });
  await verifyAndExport(
    workbook,
    path.join(tableRoot, "lactylome_and_reference_proteome_pairing_zh.xlsx"),
    "pairing37",
  );
}

async function buildVennWorkbook() {
  const workbook = Workbook.create();
  addReadme(workbook, "四分类 Kla 与普通全蛋白 Venn/Euler 审计表", [
    "四分类顺序：正常/非肿瘤组织、癌症组织、正常/非肿瘤细胞、癌症细胞。",
    "Kla样本组数量依次为10、3、10、14，总计37组；普通全蛋白严格参照组数量依次为9、2、9、13，总计33组。",
    "Kla集合与修复前完全一致；普通全蛋白集合仅使用完全匹配且具有逐蛋白定量强度的33组参照。",
    "membership表以UniProt BaseAccession为分析键，GeneSymbol和ProteinName仅用于人工审计。",
  ]);
  const analyses = [
    ["all_kla_four_class_venn", "Kla全部"],
    ["kla_ddr_four_class_venn", "Kla_DDR"],
    ["reference_proteome_four_class_venn", "普通全蛋白"],
    ["reference_proteome_ddr_four_class_venn", "普通_DDR"],
  ];
  for (const [analysis, label] of analyses) {
    const base = `reanalysis/results/tables/four_class_venn/${analysis}`;
    addDataSheet(workbook, `${label}_集合数`, await readCsv(`${base}/set_counts.csv`), {
      maxWidth: 30,
      headerFill: "#A64B2A",
    });
    addDataSheet(workbook, `${label}_区域数`, await readCsv(`${base}/region_counts.csv`), {
      maxWidth: 34,
      headerFill: "#B5651D",
    });
    addDataSheet(workbook, `${label}_成员`, await readCsv(`${base}/membership.csv`), {
      freezeColumns: 3,
      maxWidth: 38,
      headerFill: "#8F2D1D",
    });
  }
  await verifyAndExport(
    workbook,
    path.join(tableRoot, "four_class_venn_tables.xlsx"),
    "venn37",
  );
}

await fs.mkdir(tableRoot, { recursive: true });
await buildDdrWorkbook();
await buildPairingWorkbook();
await buildVennWorkbook();
console.log("Built 37-group deliverable workbooks.");
