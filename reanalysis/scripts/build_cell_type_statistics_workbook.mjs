import fs from "node:fs/promises";
import path from "node:path";
import { Workbook, SpreadsheetFile } from "@oai/artifact-tool";

const projectRoot = "/Users/gzy2520/Desktop/Research/kla";
const inputPath = path.join(
  projectRoot,
  "reanalysis/results/tables/cell_type_kla_ddr_statistics.csv",
);
const outputPath = path.join(
  projectRoot,
  "reanalysis/results/tables/cell_type_kla_ddr_statistics.xlsx",
);
const previewPath = "/tmp/kla_cell_type_statistics_preview.png";

function parseCsv(text) {
  const lines = text.trim().split(/\r?\n/);
  const headers = lines[0].replace(/^\uFEFF/, "").split(",");
  return lines.slice(1).map((line) => {
    const values = line.split(",");
    return Object.fromEntries(headers.map((header, index) => [header, values[index] ?? ""]));
  });
}

const rows = parseCsv(await fs.readFile(inputPath, "utf8"));
const workbook = Workbook.create();
const sheet = workbook.worksheets.add("按细胞类型统计");
sheet.showGridLines = false;

sheet.getRange("A1:D1").merge();
sheet.getRange("A1").values = [["Kla 蛋白与 DDR 交集统计"]];
sheet.getRange("A1:D1").format = {
  fill: "#1F4E78",
  font: { bold: true, color: "#FFFFFF", size: 16 },
  horizontalAlignment: "center",
  verticalAlignment: "center",
};
sheet.getRange("A1:D1").format.rowHeight = 32;

sheet.getRange("A2:D2").merge();
sheet.getRange("A2").values = [[
  "统计口径：每个细胞系或组织内按 UniProt BaseAccession 去重；海马体为组织样本。",
]];
sheet.getRange("A2:D2").format = {
  fill: "#D9EAF7",
  font: { color: "#17365D", size: 10 },
  wrapText: true,
  verticalAlignment: "center",
};
sheet.getRange("A2:D2").format.rowHeight = 28;

sheet.getRange("A4:D4").values = [[
  "细胞系或组织",
  "总 Kla 蛋白数",
  "与 DDR 交集蛋白数",
  "DDR 交集 / 总数",
]];
sheet.getRange("A4:D4").format = {
  fill: "#5B9BD5",
  font: { bold: true, color: "#FFFFFF" },
  horizontalAlignment: "center",
  verticalAlignment: "center",
  wrapText: true,
  borders: { preset: "outside", style: "thin", color: "#2F5597" },
};
sheet.getRange("A4:D4").format.rowHeight = 28;

const startRow = 5;
const endRow = startRow + rows.length - 1;
const values = rows.map((row) => [
  row.CellOrTissueType,
  Number(row.TotalKlaProteins),
  Number(row.KlaGoDdrProteins),
]);
sheet.getRange(`A${startRow}:C${endRow}`).values = values;
sheet.getRange(`D${startRow}`).formulas = [[`=IFERROR(C${startRow}/B${startRow},0)`]];
sheet.getRange(`D${startRow}:D${endRow}`).fillDown();

sheet.getRange(`A${startRow}:D${endRow}`).format = {
  borders: {
    insideHorizontal: { style: "thin", color: "#D9E2F3" },
    bottom: { style: "thin", color: "#A6A6A6" },
  },
  verticalAlignment: "center",
};
sheet.getRange(`B${startRow}:C${endRow}`).format = {
  numberFormat: "#,##0",
  horizontalAlignment: "right",
};
sheet.getRange(`D${startRow}:D${endRow}`).format = {
  numberFormat: "0.00%",
  horizontalAlignment: "right",
};
sheet.getRange(`A${startRow}:A${endRow}`).format.horizontalAlignment = "left";

sheet.getRange(`A${startRow}:D${endRow}`).conditionalFormats.add("expression", {
  formula: `=MOD(ROW(),2)=1`,
  format: { fill: "#F4F8FC" },
});
sheet.tables.add(`A4:D${endRow}`, true, "CellTypeKlaDdrStatistics");
sheet.freezePanes.freezeRows(4);

sheet.getRange("A:A").format.columnWidth = 25;
sheet.getRange("B:B").format.columnWidth = 18;
sheet.getRange("C:C").format.columnWidth = 22;
sheet.getRange("D:D").format.columnWidth = 20;
sheet.getRange(`A${startRow}:D${endRow}`).format.rowHeight = 22;

const inspected = await workbook.inspect({
  kind: "table",
  range: `按细胞类型统计!A1:D${endRow}`,
  include: "values,formulas",
  tableMaxRows: 20,
  tableMaxCols: 6,
});
console.log(inspected.ndjson);

const errors = await workbook.inspect({
  kind: "match",
  searchTerm: "#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A",
  options: { useRegex: true, maxResults: 100 },
  summary: "final formula error scan",
});
console.log(errors.ndjson);

const preview = await workbook.render({
  sheetName: "按细胞类型统计",
  range: `A1:D${endRow}`,
  scale: 2,
});
await fs.writeFile(previewPath, new Uint8Array(await preview.arrayBuffer()));

const output = await SpreadsheetFile.exportXlsx(workbook);
await output.save(outputPath);
console.log(JSON.stringify({ outputPath, previewPath, rowCount: rows.length }));
