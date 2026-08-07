#!/usr/bin/env node

import fs from "node:fs/promises";
import path from "node:path";
import { SpreadsheetFile, Workbook } from "@oai/artifact-tool";


const projectRoot = "/Users/gzy2520/Desktop/Research/kla";
const tableRoot = path.join(projectRoot, "reanalysis/results/tables");
const outputPath = path.join(tableRoot, "venn_combined_tables.xlsx");
const previewRoot = "/tmp/kla_venn_workbook_20260723/previews";

const sources = [
  {
    sheetName: "Kla_unique",
    tableName: "KlaUniqueTable",
    fileName: "all_kla_three_groups_combined.csv",
  },
  {
    sheetName: "Kla_non_dedup",
    tableName: "KlaNonDedupTable",
    fileName: "all_kla_three_groups_combined_non_deduplicated.csv",
  },
  {
    sheetName: "Kla_DDR_unique",
    tableName: "KlaDdrUniqueTable",
    fileName: "kla_go_ddr_three_groups_combined.csv",
  },
  {
    sheetName: "Kla_DDR_non_dedup",
    tableName: "KlaDdrNonDedupTable",
    fileName: "kla_go_ddr_three_groups_combined_non_deduplicated.csv",
  },
];

function columnName(index) {
  let value = index;
  let name = "";
  while (value > 0) {
    value -= 1;
    name = String.fromCharCode(65 + (value % 26)) + name;
    value = Math.floor(value / 26);
  }
  return name;
}

function styleDataSheet(sheet, tableName) {
  const used = sheet.getUsedRange(true);
  const rowCount = used.rowCount;
  const columnCount = used.columnCount;
  const headers = used.getRow(0).values[0].map(String);

  sheet.showGridLines = false;
  sheet.freezePanes.freezeRows(1);
  sheet.freezePanes.freezeColumns(4);
  used.format.font = { size: 10, color: "#1F2937" };
  used.format.verticalAlignment = "center";
  used.getRow(0).format = {
    fill: "#1F4E78",
    font: { bold: true, color: "#FFFFFF", size: 10 },
    verticalAlignment: "center",
    wrapText: true,
    rowHeight: 34,
    borders: { preset: "outside", style: "thin", color: "#17365D" },
  };

  const widths = {
    BaseAccession: 16,
    GeneSymbol: 18,
    ProteinName: 42,
    KlaSites: 30,
    SourceCategory: 27,
    InHippocampusTissue: 20,
    InNormalImmortalizedCellLines: 28,
    InTumorCellLines: 18,
    DetectedGroupCount: 18,
    VennRegion: 31,
    AnalysisSet: 16,
    PXD: 18,
    Sample: 32,
    CellType: 24,
    ExperimentalGroup: 32,
    Category: 29,
    SourceFile: 60,
    EvidenceMode: 38,
    LocalizationProb: 18,
    SourceConfidence: 34,
    EvidenceRows: 14,
  };
  headers.forEach((header, index) => {
    sheet.getRangeByIndexes(0, index, rowCount, 1).format.columnWidth = widths[header] ?? 18;
  });

  const lastColumn = columnName(columnCount);
  sheet.tables.add(`A1:${lastColumn}${rowCount}`, true, tableName);
  return { rowCount, columnCount, lastColumn };
}

const firstCsv = await fs.readFile(path.join(tableRoot, sources[0].fileName), "utf8");
const workbook = await Workbook.fromCSV(firstCsv, { sheetName: sources[0].sheetName });
for (const source of sources.slice(1)) {
  const csvText = await fs.readFile(path.join(tableRoot, source.fileName), "utf8");
  await workbook.fromCSV(csvText, { sheetName: source.sheetName });
}

const dimensions = {};
for (const source of sources) {
  const sheet = workbook.worksheets.getItem(source.sheetName);
  dimensions[source.sheetName] = styleDataSheet(sheet, source.tableName);
}

const readme = workbook.worksheets.add("README");
readme.showGridLines = false;
readme.getRange("A1:C1").merge();
readme.getRange("A1").values = [["Kla Venn combined protein tables"]];
readme.getRange("A1:C1").format = {
  fill: "#1F4E78",
  font: { bold: true, color: "#FFFFFF", size: 16 },
  verticalAlignment: "center",
  rowHeight: 30,
};
readme.getRange("A3:C8").values = [
  ["Sheet", "Rows", "Meaning"],
  ["Kla_unique", null, "All Kla proteins; one row per BaseAccession"],
  ["Kla_non_dedup", null, "Three Kla groups appended; intersection proteins repeat by SourceCategory"],
  ["Kla_DDR_unique", null, "Kla intersected with GO repair/damage; one row per BaseAccession"],
  ["Kla_DDR_non_dedup", null, "Three Kla-DDR groups appended; intersection proteins repeat by SourceCategory"],
  ["Membership", null, "Yes/No columns and VennRegion reproduce the exact Venn inputs"],
];
readme.getRange("B4:B7").formulas = [
  [`=COUNTA('Kla_unique'!A2:A${dimensions.Kla_unique.rowCount})`],
  [`=COUNTA('Kla_non_dedup'!A2:A${dimensions.Kla_non_dedup.rowCount})`],
  [`=COUNTA('Kla_DDR_unique'!A2:A${dimensions.Kla_DDR_unique.rowCount})`],
  [`=COUNTA('Kla_DDR_non_dedup'!A2:A${dimensions.Kla_DDR_non_dedup.rowCount})`],
];
readme.getRange("A3:C3").format = {
  fill: "#D9EAF7",
  font: { bold: true, color: "#17365D" },
  borders: { preset: "outside", style: "thin", color: "#9FBAD0" },
};
readme.getRange("B4:B7").format.numberFormat = "#,##0";
readme.getRange("A10:C14").values = [
  ["Source CSV files", null, null],
  [sources[0].fileName, null, null],
  [sources[1].fileName, null, null],
  [sources[2].fileName, null, null],
  [sources[3].fileName, null, null],
];
readme.getRange("A10:C10").format = {
  fill: "#E2F0D9",
  font: { bold: true, color: "#375623" },
};
readme.getRange("A1:A14").format.columnWidth = 38;
readme.getRange("B1:B14").format.columnWidth = 14;
readme.getRange("C1:C14").format.columnWidth = 76;
readme.freezePanes.freezeRows(3);

await fs.mkdir(previewRoot, { recursive: true });
for (const sheetName of ["README", ...sources.map((source) => source.sheetName)]) {
  const range = sheetName === "README" ? "A1:C14" : "A1:U14";
  const preview = await workbook.render({ sheetName, range, scale: 1, format: "png" });
  await fs.writeFile(
    path.join(previewRoot, `${sheetName}.png`),
    new Uint8Array(await preview.arrayBuffer()),
  );
}

const inspection = await workbook.inspect({
  kind: "sheet,table",
  maxChars: 8000,
  tableMaxRows: 4,
  tableMaxCols: 6,
});
console.log(inspection.ndjson);
const errors = await workbook.inspect({
  kind: "match",
  searchTerm: "#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A",
  options: { useRegex: true, maxResults: 100 },
  summary: "final formula error scan",
});
console.log(errors.ndjson);

const output = await SpreadsheetFile.exportXlsx(workbook);
await output.save(outputPath);
console.log(`Saved ${outputPath}`);
