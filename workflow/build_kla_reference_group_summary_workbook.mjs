#!/usr/bin/env node

import crypto from "node:crypto";
import fs from "node:fs/promises";
import path from "node:path";
import process from "node:process";
import {
  FileBlob,
  SpreadsheetFile,
  Workbook,
} from "@oai/artifact-tool";

const scriptPath = path.resolve(process.argv[1]);
const projectRoot = process.argv[2]
  ? path.resolve(process.argv[2])
  : path.resolve(path.dirname(scriptPath), "..");
const outputDir = process.argv[3]
  ? path.resolve(process.argv[3])
  : path.join(
      "/Users/gzy2520/.codex/visualizations/2026/08/18",
      "01a01339-e298-7cb3-8f94-4905c15ffd09",
      "outputs",
    );
const projectOutputDir = path.join(projectRoot, "results", "supplementary");
const outputFilename =
  "Supplementary_Data_S1_Kla_and_Reference_Group_Summary_20260818.xlsx";
const inputWorkbookPath = path.join(
  projectRoot,
  "results/supplementary/Supplementary_Data_S1_Kla_Evidence_Draft_20260817.xlsx",
);
const pairedPath = path.join(
  projectRoot,
  "results/tables/cell_type_kla_vs_reference_ddr_statistics_accession_only_paired_30.csv",
);
const previewDir = path.join(outputDir, "previews");

const palette = {
  navy: "#1F4E78",
  text: "#243447",
  white: "#FFFFFF",
};

function assert(condition, message) {
  if (!condition) throw new Error(message);
}

function cleanString(value) {
  if (value === null || value === undefined) return "";
  return String(value).replace(/^\uFEFF/, "").trim();
}

function normalizeMatrix(matrix) {
  return matrix.map((row) =>
    row.map((value) => {
      if (value === undefined) return null;
      if (typeof value === "number" && !Number.isFinite(value)) return null;
      return value;
    }),
  );
}

async function readCsvMatrix(relativePath) {
  const text = (await fs.readFile(path.join(projectRoot, relativePath), "utf8"))
    .replace(/^\uFEFF/, "");
  const workbook = await Workbook.fromCSV(text, { sheetName: "Data" });
  const used = workbook.worksheets.getItemAt(0).getUsedRange(true);
  assert(used, `No used range found in ${relativePath}.`);
  return normalizeMatrix(used.values);
}

async function readXlsxSheetMatrix(absolutePath, sheetName) {
  const workbook = await SpreadsheetFile.importXlsx(
    await FileBlob.load(absolutePath),
  );
  const sheet = workbook.worksheets.getItem(sheetName);
  const used = sheet.getUsedRange(true);
  assert(used, `No used range found in ${sheetName}.`);
  return normalizeMatrix(used.values);
}

function objectsFromMatrix(matrix) {
  const headers = matrix[0].map(cleanString);
  return matrix.slice(1).map((row) => {
    const object = {};
    headers.forEach((header, index) => {
      object[header] = row[index] ?? null;
    });
    return object;
  });
}

function matrixFromObjects(objects, headers) {
  return [
    headers,
    ...objects.map((object) =>
      headers.map((header) =>
        object[header] === undefined ? null : object[header],
      ),
    ),
  ];
}

function excelColumn(indexZeroBased) {
  let value = indexZeroBased + 1;
  let result = "";
  while (value > 0) {
    value -= 1;
    result = String.fromCharCode(65 + (value % 26)) + result;
    value = Math.floor(value / 26);
  }
  return result;
}

function chooseColumnWidth(header, values) {
  const headerLength = cleanString(header).length;
  const valueLength = Math.max(
    0,
    ...values.map((value) => cleanString(value).length),
  );
  const longText = /Evidence|Material|SampleGroup|Biological|Label|File/i.test(
    cleanString(header),
  );
  if (longText) return Math.min(Math.max(Math.max(headerLength, valueLength) + 2, 18), 48);
  return Math.min(Math.max(Math.max(headerLength, valueLength) + 2, 10), 28);
}

function applyNumberFormat(sheet, header, columnIndex, rowCount) {
  if (rowCount <= 1) return;
  const range = sheet.getRangeByIndexes(1, columnIndex, rowCount - 1, 1);
  if (/Fraction/i.test(header)) range.format.numberFormat = "0.0%";
  if (/Count|RecordNo/i.test(header)) range.format.numberFormat = "#,##0";
}

function styleDataSheet(sheet, matrix, tableName) {
  const rows = matrix.length;
  const cols = matrix[0].length;
  sheet.showGridLines = false;
  sheet.freezePanes.freezeRows(1);
  const used = sheet.getRangeByIndexes(0, 0, rows, cols);
  used.format.font = { name: "Arial Unicode MS", size: 10, color: palette.text };
  used.format.verticalAlignment = "center";

  const header = sheet.getRangeByIndexes(0, 0, 1, cols);
  header.format = {
    fill: palette.navy,
    font: {
      name: "Arial Unicode MS",
      size: 10,
      bold: true,
      color: palette.white,
    },
    horizontalAlignment: "center",
    verticalAlignment: "center",
    wrapText: true,
    borders: {
      bottom: { style: "medium", color: palette.navy },
    },
  };
  header.format.rowHeight = 30;

  for (let index = 0; index < cols; index += 1) {
    const samples = matrix.slice(1, Math.min(rows, 100)).map((row) => row[index]);
    sheet
      .getRangeByIndexes(0, index, rows, 1)
      .format.columnWidth = chooseColumnWidth(matrix[0][index], samples);
    applyNumberFormat(sheet, cleanString(matrix[0][index]), index, rows);
  }
  const lastCell = `${excelColumn(cols - 1)}${rows}`;
  const table = sheet.tables.add(`A1:${lastCell}`, true, tableName);
  table.style = "TableStyleMedium2";
  table.showFilterButton = true;
  table.showBandedRows = true;
}

async function sha256(filePath) {
  return crypto
    .createHash("sha256")
    .update(await fs.readFile(filePath))
    .digest("hex");
}

async function main() {
  await fs.mkdir(outputDir, { recursive: true });
  await fs.mkdir(projectOutputDir, { recursive: true });
  await fs.mkdir(previewDir, { recursive: true });

  const sourceSheet1 = await readXlsxSheetMatrix(inputWorkbookPath, "Group_Summary");
  assert(sourceSheet1.length === 31, "The current Kla Group_Summary must contain 30 rows.");

  const paired = objectsFromMatrix(
    await readCsvMatrix(
      "results/tables/cell_type_kla_vs_reference_ddr_statistics_accession_only_paired_30.csv",
    ),
  ).sort((left, right) => Number(left.RowOrder) - Number(right.RowOrder));
  assert(paired.length === 30, "The reference Group_Summary must contain 30 rows.");

  const referenceHeaders = [
    "RecordNo",
    "Category",
    "ReferencePXD",
    "SampleGroup",
    "BiologicalMaterial",
    "ReferenceLabelEn",
    "ReferenceProteinCount",
    "ReferenceDdrProteinCount",
    "ReferenceDdrFraction",
    "ReferenceEvidenceFile",
  ];
  const referenceRows = paired.map((row, index) => ({
    RecordNo: index + 1,
    Category: cleanString(row.CategoryEn),
    ReferencePXD: cleanString(row.ReferencePXD),
    SampleGroup: cleanString(row.SampleGroup),
    BiologicalMaterial: cleanString(row.BiologicalMaterial),
    ReferenceLabelEn: cleanString(row.ReferenceLabelEn),
    ReferenceProteinCount: Number(row.ReferenceProteinCount),
    ReferenceDdrProteinCount: Number(row.ReferenceDdrProteinCount),
    ReferenceDdrFraction: Number(row.ReferenceDdrFraction),
    ReferenceEvidenceFile: cleanString(row.ReferenceEvidenceFile),
  }));
  assert(
    referenceRows.every((row) => row.ReferenceProteinCount > 0),
    "Every current reference group must have a positive full-proteome count.",
  );

  const workbook = Workbook.create();
  const klaSheet = workbook.worksheets.add("Group_Summary");
  klaSheet.getRangeByIndexes(
    0,
    0,
    sourceSheet1.length,
    sourceSheet1[0].length,
  ).writeValues(sourceSheet1);
  styleDataSheet(klaSheet, sourceSheet1, "Current30KlaGroupSummaryTable");

  const referenceMatrix = matrixFromObjects(referenceRows, referenceHeaders);
  const referenceSheet = workbook.worksheets.add("Reference_Group_Summary");
  referenceSheet
    .getRangeByIndexes(0, 0, referenceMatrix.length, referenceMatrix[0].length)
    .writeValues(referenceMatrix);
  styleDataSheet(
    referenceSheet,
    referenceMatrix,
    "Current30ReferenceGroupSummaryTable",
  );

  for (const [sheetName, range] of [
    ["Group_Summary", "A1:I14"],
    ["Reference_Group_Summary", "A1:J14"],
  ]) {
    const preview = await workbook.render({
      sheetName,
      range,
      scale: 1,
      format: "png",
    });
    await fs.writeFile(
      path.join(previewDir, `${sheetName}.png`),
      new Uint8Array(await preview.arrayBuffer()),
    );
  }

  const errors = await workbook.inspect({
    kind: "match",
    searchTerm: "#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A",
    options: { useRegex: true, maxResults: 100 },
    summary: "Kla/reference Group_Summary formula scan",
    maxChars: 3000,
  });
  assert(
    !/#REF!|#DIV\/0!|#VALUE!|#NAME\?|#N\/A/.test(errors.ndjson),
    "The new workbook contains formula errors.",
  );

  const output = await SpreadsheetFile.exportXlsx(workbook);
  const externalPath = path.join(outputDir, outputFilename);
  const projectPath = path.join(projectOutputDir, outputFilename);
  await output.save(externalPath);
  await fs.copyFile(externalPath, projectPath);
  const manifest = {
    file: outputFilename,
    projectPath,
    externalPath,
    sheet1Rows: sourceSheet1.length - 1,
    sheet2Rows: referenceRows.length,
    sheet2ReferenceProteinTotal: referenceRows.reduce(
      (sum, row) => sum + row.ReferenceProteinCount,
      0,
    ),
    sha256: await sha256(projectPath),
  };
  await fs.writeFile(
    path.join(projectOutputDir, "kla_reference_group_summary_manifest.json"),
    `${JSON.stringify(manifest, null, 2)}\n`,
    "utf8",
  );
  console.log(JSON.stringify(manifest, null, 2));
}

await main();
