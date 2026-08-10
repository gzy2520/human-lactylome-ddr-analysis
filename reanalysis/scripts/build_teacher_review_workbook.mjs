import fs from "node:fs/promises";
import path from "node:path";
import { SpreadsheetFile, Workbook } from "@oai/artifact-tool";

const projectRoot = process.argv[2] || "/Users/gzy2520/Desktop/Research/kla";
const tableRoot = path.join(projectRoot, "reanalysis/results/tables");
const previewRoot = "/tmp/kla-teacher-review-previews";
const outputPath = path.join(tableRoot, "kla_and_reference_teacher_review_zh.xlsx");

function parseCsv(text) {
  const rows = [];
  let row = [];
  let field = "";
  let quoted = false;
  const source = text.replace(/^\uFEFF/, "");
  for (let index = 0; index < source.length; index += 1) {
    const character = source[index];
    if (quoted) {
      if (character === '"' && source[index + 1] === '"') {
        field += '"';
        index += 1;
      } else if (character === '"') {
        quoted = false;
      } else {
        field += character;
      }
    } else if (character === '"') {
      quoted = true;
    } else if (character === ",") {
      row.push(field);
      field = "";
    } else if (character === "\n") {
      row.push(field.replace(/\r$/, ""));
      rows.push(row);
      row = [];
      field = "";
    } else {
      field += character;
    }
  }
  if (field.length || row.length) {
    row.push(field.replace(/\r$/, ""));
    rows.push(row);
  }
  return rows.filter((values) => values.some((value) => value !== ""));
}

async function readCsv(name) {
  const content = await fs.readFile(path.join(tableRoot, name), "utf8");
  return parseCsv(content);
}

function typedValue(value) {
  if (value === "") return null;
  if (/^-?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?$/.test(value)) {
    const parsed = Number(value);
    if (Number.isFinite(parsed)) return parsed;
  }
  return value;
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

function writeDataSheet(workbook, name, csvRows, options = {}) {
  const sheet = workbook.worksheets.add(name);
  sheet.showGridLines = false;
  const rows = csvRows.map((row, rowIndex) =>
    row.map((value) => (rowIndex === 0 ? value : typedValue(value))),
  );
  const columnCount = Math.max(...rows.map((row) => row.length));
  const normalized = rows.map((row) => [
    ...row,
    ...Array(columnCount - row.length).fill(null),
  ]);
  const lastColumn = columnName(columnCount - 1);
  const lastRow = normalized.length;
  const usedRange = sheet.getRange(`A1:${lastColumn}${lastRow}`);
  usedRange.values = normalized;
  usedRange.format = {
    font: { color: "#272321", size: 9 },
    verticalAlignment: "top",
    wrapText: true,
  };
  sheet.getRange(`A1:${lastColumn}1`).format = {
    fill: options.headerFill || "#8F2D1D",
    font: { bold: true, color: "#FFFFFF", size: 10 },
    horizontalAlignment: "center",
    verticalAlignment: "center",
    wrapText: true,
    borders: { preset: "outside", style: "thin", color: "#6E2015" },
  };
  sheet.getRange(`A2:${lastColumn}${lastRow}`).format.borders = {
    insideHorizontal: { style: "thin", color: "#E6DED8" },
  };
  sheet.getRange(`A1:${lastColumn}1`).format.rowHeight = 38;
  sheet.getRange(`A2:${lastColumn}${lastRow}`).format.rowHeight =
    options.rowHeight || 54;
  sheet.freezePanes.freezeRows(1);
  sheet.freezePanes.freezeColumns(options.freezeColumns || 4);
  const table = sheet.tables.add(
    `A1:${lastColumn}${lastRow}`,
    true,
    options.tableName,
  );
  table.style = "TableStyleMedium2";
  table.showFilterButton = true;
  return { sheet, lastRow, lastColumn, columnCount };
}

function setColumnWidths(sheet, widths) {
  widths.forEach((width, index) => {
    const column = columnName(index);
    sheet.getRange(`${column}:${column}`).format.columnWidth = width;
  });
}

function buildSummarySheet(workbook) {
  const sheet = workbook.worksheets.getItem("说明与汇总");
  sheet.showGridLines = false;
  sheet.getRange("A1:H1").merge();
  sheet.getRange("A1").values = [["Kla与普通全蛋白参照：老师审阅表"]];
  sheet.getRange("A1:H1").format = {
    fill: "#8F2D1D",
    font: { bold: true, color: "#FFFFFF", size: 16 },
    verticalAlignment: "center",
  };
  sheet.getRange("A1:H1").format.rowHeight = 34;

  sheet.getRange("A3:A7").values = [
    ["最终Kla样本组"],
    ["纳入严格配对分析"],
    ["未纳入普通全蛋白配对"],
    ["分析匹配键"],
    ["GeneSymbol回退"],
  ];
  sheet.getRange("B3:B7").formulas = [
    ["=COUNTA('老师审阅总表'!$C$2:$C$34)"],
    ["=COUNTIF('老师审阅总表'!$Y$2:$Y$34,\"是\")"],
    ["=COUNTIF('老师审阅总表'!$Y$2:$Y$34,\"否\")"],
    ['="UniProt BaseAccession"'],
    ['="0（禁止参与分析命中）"'],
  ];
  sheet.getRange("A3:B7").format = {
    fill: "#FFF6EC",
    font: { size: 11 },
    borders: { preset: "outside", style: "thin", color: "#D6B7A3" },
    verticalAlignment: "center",
  };
  sheet.getRange("A3:A7").format.font = { bold: true, color: "#6D2C1B" };
  sheet.getRange("A3:B7").format.rowHeight = 25;

  sheet.getRange("D3:F3").values = [[
    "四分类",
    "Kla组数",
    "有严格参照组数",
  ]];
  sheet.getRange("D4:D7").values = [
    ["正常/非肿瘤组织"],
    ["癌症组织"],
    ["正常/非肿瘤细胞"],
    ["癌症细胞"],
  ];
  sheet.getRange("E4:E7").formulas = [4, 5, 6, 7].map((row) => [
    `=COUNTIF('老师审阅总表'!$B$2:$B$34,D${row})`,
  ]);
  sheet.getRange("F4:F7").formulas = [4, 5, 6, 7].map((row) => [
    `=COUNTIFS('老师审阅总表'!$B$2:$B$34,D${row},'老师审阅总表'!$Y$2:$Y$34,"是")`,
  ]);
  sheet.getRange("D3:F7").format = {
    borders: { preset: "outside", style: "thin", color: "#D6B7A3" },
    verticalAlignment: "center",
  };
  sheet.getRange("D3:F3").format = {
    fill: "#B45F35",
    font: { bold: true, color: "#FFFFFF" },
    horizontalAlignment: "center",
  };
  sheet.getRange("D4:F7").format.fill = "#FFF9F3";
  sheet.getRange("D3:F7").format.rowHeight = 25;

  const notes = [
    "阅读口径：每行是一组最终Kla样本。Kla与普通全蛋白均按去除isoform后缀的UniProt BaseAccession计数；GeneSymbol只用于显示或人工核对。",
    "“材料身份严格匹配=是”表示细胞/组织身份及取材粒度可对应，不等同于同一供体、同一处理或同一时间点。实验状态限制已单独列出。",
    "最终正式分析固定为33组，且33组均进入普通全蛋白配对分析；已删除的4组不再保留在活动Kla数据中。",
    "同一PXD或同一文件被多个组使用时，必须查看“实际读取样本/工作表/列”。例如HCC读取CISs，邻近肝读取ANTs；两组没有混用。",
    "两个HK-2组读取PXD072220的amostra1、amostra3、amostra4 PG.Log2Quantity，文件表头前两行跳过。海马读取PXD050470同研究Table S4的H072/H081/H0187。",
    "需老师重点确认的情况集中在“需重点确认”页，包括基线参照、独立队列和独立实验；已删除组不再显示。",
  ];
  sheet.getRange("A10:H10").merge();
  sheet.getRange("A10").values = [["审阅说明"]];
  sheet.getRange("A10:H10").format = {
    fill: "#E8B04A",
    font: { bold: true, color: "#3C2816", size: 11 },
  };
  notes.forEach((note, index) => {
    const row = 11 + index;
    sheet.getRange(`A${row}:H${row}`).merge();
    sheet.getRange(`A${row}`).values = [[note]];
    sheet.getRange(`A${row}:H${row}`).format = {
      fill: index % 2 === 0 ? "#FFF9F0" : "#FFFFFF",
      font: { color: "#3A302B", size: 10 },
      wrapText: true,
      verticalAlignment: "center",
    };
    sheet.getRange(`A${row}:H${row}`).format.rowHeight = 34;
  });
  sheet.getRange("A:H").format.columnWidth = 15;
  sheet.getRange("A:A").format.columnWidth = 24;
  sheet.getRange("B:B").format.columnWidth = 22;
  sheet.getRange("D:D").format.columnWidth = 24;
  sheet.freezePanes.freezeRows(1);
  return sheet;
}

async function main() {
  const reviewRows = await readCsv("kla_and_reference_teacher_review_zh.csv");
  const focusRows = await readCsv("kla_and_reference_teacher_review_focus_zh.csv");
  const workbook = Workbook.create();
  workbook.worksheets.add("说明与汇总");

  const review = writeDataSheet(workbook, "老师审阅总表", reviewRows, {
    tableName: "TeacherReviewTable",
    freezeColumns: 4,
    rowHeight: 66,
    headerFill: "#8F2D1D",
  });
  setColumnWidths(review.sheet, [
    7, 18, 13, 28, 30, 20, 48, 24, 22, 12,
    14, 13, 18, 30, 20, 48, 40, 28, 22, 15,
    15, 15, 13, 38, 13, 38, 56, 28, 42,
  ]);
  review.sheet.getRange(`J2:L${review.lastRow}`).format.numberFormat = "#,##0";
  review.sheet.getRange(`L2:L${review.lastRow}`).format.numberFormat = "0.0%";
  review.sheet.getRange(`T2:V${review.lastRow}`).format.numberFormat = "#,##0";
  review.sheet.getRange(`V2:V${review.lastRow}`).format.numberFormat = "0.0%";
  review.sheet.getRange(`W2:W${review.lastRow}`).conditionalFormats.add(
    "containsText",
    { text: "否", format: { fill: "#FDE8E6", font: { color: "#A5261B", bold: true } } },
  );
  review.sheet.getRange(`Y2:Y${review.lastRow}`).conditionalFormats.add(
    "containsText",
    { text: "是", format: { fill: "#E8F2E8", font: { color: "#2E6B3B", bold: true } } },
  );
  review.sheet.getRange(`Y2:Y${review.lastRow}`).conditionalFormats.add(
    "containsText",
    { text: "否", format: { fill: "#FDE8E6", font: { color: "#A5261B", bold: true } } },
  );
  review.sheet.getRange(`X2:X${review.lastRow}`).conditionalFormats.add(
    "containsText",
    { text: "独立队列", format: { fill: "#FFF0C9", font: { color: "#7A4A00" } } },
  );
  review.sheet.getRange(`X2:X${review.lastRow}`).conditionalFormats.add(
    "containsText",
    { text: "基线参照", format: { fill: "#FFF7DE", font: { color: "#6F4A00" } } },
  );

  const focus = writeDataSheet(workbook, "需重点确认", focusRows, {
    tableName: "TeacherFocusTable",
    freezeColumns: 4,
    rowHeight: 68,
    headerFill: "#B45F35",
  });
  setColumnWidths(focus.sheet, [
    7, 18, 13, 30, 32, 18, 32, 42, 40, 14, 64,
  ]);
  focus.sheet.getRange(`J2:J${focus.lastRow}`).conditionalFormats.add(
    "containsText",
    { text: "否", format: { fill: "#FDE8E6", font: { color: "#A5261B", bold: true } } },
  );
  buildSummarySheet(workbook);

  const inspection = await workbook.inspect({
    kind: "table",
    range: "老师审阅总表!A1:AC8",
    include: "values,formulas",
    tableMaxRows: 8,
    tableMaxCols: 29,
    maxChars: 9000,
  });
  console.log(inspection.ndjson);
  const errors = await workbook.inspect({
    kind: "match",
    searchTerm: "#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A",
    options: { useRegex: true, maxResults: 100 },
    summary: "teacher review workbook formula error scan",
  });
  console.log(errors.ndjson);

  await fs.rm(previewRoot, { recursive: true, force: true });
  await fs.mkdir(previewRoot, { recursive: true });
  const previews = [
    ["说明与汇总", "A1:H17", "summary"],
    ["老师审阅总表", "A1:O14", "review_left"],
    ["老师审阅总表", "P1:AC14", "review_right"],
    ["需重点确认", `A1:K${Math.min(focus.lastRow, 16)}`, "focus"],
  ];
  for (const [sheetName, range, fileName] of previews) {
    const image = await workbook.render({
      sheetName,
      range,
      scale: 1.25,
      format: "png",
    });
    await fs.writeFile(
      path.join(previewRoot, `${fileName}.png`),
      new Uint8Array(await image.arrayBuffer()),
    );
  }

  await fs.mkdir(tableRoot, { recursive: true });
  const output = await SpreadsheetFile.exportXlsx(workbook);
  await output.save(outputPath);
  await fs.rm(`${outputPath}.inspect.ndjson`, { force: true });
  console.log(`Saved ${outputPath}`);
}

await main();
