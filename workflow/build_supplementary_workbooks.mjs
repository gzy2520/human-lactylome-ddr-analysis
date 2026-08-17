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

process.on("uncaughtException", (error) => {
  console.error("KLA_SUPPLEMENT_BUILD_ERROR");
  console.error(error?.message ?? String(error));
  console.error(error?.stack ?? "");
  process.exit(1);
});

const scriptPath = path.resolve(process.argv[1]);
const projectRoot = process.argv[2]
  ? path.resolve(process.argv[2])
  : path.resolve(path.dirname(scriptPath), "..");
const externalOutputDir = process.argv[3]
  ? path.resolve(process.argv[3])
  : path.join(
      "/tmp",
      "kla-supplementary-workbooks",
      "outputs",
      "01a00a28-237a-75e2-8bf2-91c62d0c0110",
    );
const projectOutputDir = path.join(projectRoot, "results", "supplementary");
const previewDir = path.join(externalOutputDir, "previews");
const releaseDate = "2026-08-17";

const outputFiles = {
  tables: `Supplementary_Tables_S1-S6_Draft_${releaseDate.replaceAll("-", "")}.xlsx`,
  dataS1: `Supplementary_Data_S1_Kla_Evidence_Draft_${releaseDate.replaceAll("-", "")}.xlsx`,
  dataS2: `Supplementary_Data_S2_Human_GO_DDR_Draft_${releaseDate.replaceAll("-", "")}.xlsx`,
  dataS3: `Supplementary_Data_S3_GO_Pathway_Draft_${releaseDate.replaceAll("-", "")}.xlsx`,
};

const palette = {
  navy: "#1F4E78",
  blue: "#D9EAF7",
  pale: "#F4F7FA",
  border: "#B8C4CE",
  text: "#243447",
  white: "#FFFFFF",
  warning: "#FFF2CC",
  success: "#E2F0D9",
};

function assert(condition, message) {
  if (!condition) throw new Error(message);
}

function cleanString(value) {
  if (value === null || value === undefined) return "";
  return String(value).replace(/^\uFEFF/, "").trim();
}

function isTrue(value) {
  if (value === true || value === 1) return true;
  return ["TRUE", "T", "YES", "Y", "1"].includes(cleanString(value).toUpperCase());
}

function keyOf(pxd, sampleGroup) {
  return `${cleanString(pxd)}__${cleanString(sampleGroup)}`;
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
  const absolutePath = path.join(projectRoot, relativePath);
  const csvText = (await fs.readFile(absolutePath, "utf8")).replace(/^\uFEFF/, "");
  const temporary = await Workbook.fromCSV(csvText, { sheetName: "Data" });
  const sheet = temporary.worksheets.getItemAt(0);
  const used = sheet.getUsedRange(true);
  assert(used, `No used range found in ${relativePath}`);
  return normalizeMatrix(used.values);
}

async function readTsvMatrix(relativePath) {
  const absolutePath = path.join(projectRoot, relativePath);
  const text = (await fs.readFile(absolutePath, "utf8")).replace(/^\uFEFF/, "");
  const lines = text.split(/\r?\n/).filter((line) => line.length > 0);
  return lines.map((line) => line.split("\t"));
}

async function readXlsxSheetMatrix(relativePath, sheetName) {
  const absolutePath = path.join(projectRoot, relativePath);
  const blob = await FileBlob.load(absolutePath);
  const sourceWorkbook = await SpreadsheetFile.importXlsx(blob);
  const sheet = sourceWorkbook.worksheets.getItem(sheetName);
  const used = sheet.getUsedRange(true);
  assert(used, `No used range found in ${relativePath}:${sheetName}`);
  return normalizeMatrix(used.values);
}

function objectsFromMatrix(matrix) {
  assert(matrix.length >= 1, "Cannot convert an empty matrix to objects.");
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

function chooseColumnWidth(header, sampleValues) {
  const headerLength = cleanString(header).length;
  let maxLength = headerLength;
  for (const value of sampleValues) {
    const text = value === null || value === undefined ? "" : String(value);
    const longestLine = Math.max(...text.split(/\r?\n/).map((part) => part.length));
    maxLength = Math.max(maxLength, longestLine);
  }
  const longText = /Reason|Caveat|Evidence|Definition|Reference|Source|File|URL|DOI|Note|Title|Material|Subset|Peptide|ProteinName|GOName|Rationale/i.test(
    cleanString(header),
  );
  const identifier = /Accession|PXD|GO_ID|GeneSymbol|Taxon|ID$/i.test(cleanString(header));
  if (longText) return Math.min(Math.max(maxLength + 2, 20), 44);
  if (identifier) return Math.min(Math.max(maxLength + 2, 12), 24);
  return Math.min(Math.max(maxLength + 2, 10), 28);
}

function applyColumnNumberFormat(sheet, header, columnIndex, rowCount) {
  if (rowCount <= 1) return;
  const range = sheet.getRangeByIndexes(1, columnIndex, rowCount - 1, 1);
  const name = cleanString(header);
  if (/Fraction|Proportion|Percent/i.test(name)) {
    range.format.numberFormat = "0.0%";
  } else if (/Count|Order|Rank|Year|ID$|State|Score$|PathwayCount|ProteinCount/i.test(name)) {
    range.format.numberFormat = "#,##0";
  } else if (/Probability|Prob|PEP|Intensity|Delta|Score/i.test(name)) {
    range.format.numberFormat = "0.000";
  }
}

function styleDataSheet(sheet, matrix, tableName) {
  const rowCount = matrix.length;
  const columnCount = matrix[0]?.length ?? 0;
  assert(rowCount >= 1 && columnCount >= 1, `Sheet ${sheet.name} is empty.`);

  sheet.showGridLines = false;
  sheet.freezePanes.freezeRows(1);
  const used = sheet.getRangeByIndexes(0, 0, rowCount, columnCount);
  used.format.font = { name: "Arial Unicode MS", size: 10, color: palette.text };
  used.format.verticalAlignment = "center";

  const header = sheet.getRangeByIndexes(0, 0, 1, columnCount);
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

  const sampleLimit = Math.min(rowCount, 250);
  for (let columnIndex = 0; columnIndex < columnCount; columnIndex += 1) {
    const samples = matrix
      .slice(1, sampleLimit)
      .map((row) => row[columnIndex]);
    const columnRange = sheet.getRangeByIndexes(0, columnIndex, rowCount, 1);
    columnRange.format.columnWidth = chooseColumnWidth(
      matrix[0][columnIndex],
      samples,
    );
    applyColumnNumberFormat(
      sheet,
      matrix[0][columnIndex],
      columnIndex,
      rowCount,
    );
  }

  if (rowCount > 1) {
    const lastCell = `${excelColumn(columnCount - 1)}${rowCount}`;
    const table = sheet.tables.add(`A1:${lastCell}`, true, tableName);
    table.style = "TableStyleMedium2";
    table.showFilterButton = true;
    table.showBandedRows = true;
  }
}

function writeMatrixInChunks(sheet, matrix, chunkSize = 2000) {
  const normalized = normalizeMatrix(matrix);
  const columnCount = normalized[0]?.length ?? 0;
  assert(columnCount > 0, `Cannot write empty matrix to ${sheet.name}.`);
  for (let start = 0; start < normalized.length; start += chunkSize) {
    const chunk = normalized.slice(start, Math.min(start + chunkSize, normalized.length));
    sheet
      .getRangeByIndexes(start, 0, chunk.length, columnCount)
      .writeValues(chunk);
  }
}

function addMatrixSheet(workbook, sheetName, matrix, tableName) {
  const sheet = workbook.worksheets.add(sheetName);
  writeMatrixInChunks(sheet, matrix);
  styleDataSheet(sheet, matrix, tableName);
  return {
    sheet,
    meta: {
      name: sheetName,
      rows: matrix.length - 1,
      columns: matrix[0].length,
      headers: matrix[0].map(cleanString),
    },
  };
}

async function addCsvSheet(workbook, sheetName, relativePath, tableName) {
  const matrix = await readCsvMatrix(relativePath);
  return addMatrixSheet(workbook, sheetName, matrix, tableName);
}

function addReadmeSheet(workbook, title, descriptionRows, contentsRows) {
  const sheet = workbook.worksheets.add("README");
  sheet.showGridLines = false;
  sheet.mergeCells("A1:D1");
  sheet.getRange("A1:D1").values = [[title]];
  sheet.getRange("A1:D1").format = {
    fill: palette.navy,
    font: {
      name: "Arial Unicode MS",
      size: 16,
      bold: true,
      color: palette.white,
    },
    verticalAlignment: "center",
  };
  sheet.getRange("A1:D1").format.rowHeight = 28;

  const metadata = [
    ["Release date", releaseDate],
    ["Status", "Draft supplementary material for manuscript preparation"],
    ...descriptionRows,
  ];
  sheet.getRangeByIndexes(2, 0, metadata.length, 2).writeValues(metadata);
  for (let index = 0; index < metadata.length; index += 1) {
    sheet.mergeCells(`B${index + 3}:D${index + 3}`);
  }
  sheet.getRangeByIndexes(2, 0, metadata.length, 1).format = {
    fill: palette.blue,
    font: { name: "Arial Unicode MS", bold: true, color: palette.text },
  };
  sheet.getRangeByIndexes(2, 1, metadata.length, 1).format = {
    font: { name: "Arial Unicode MS", color: palette.text },
    wrapText: true,
  };

  const contentsStart = metadata.length + 4;
  sheet.getRangeByIndexes(contentsStart, 0, 1, 4).values = [
    ["Worksheet", "Rows", "Purpose", "Primary source"],
  ];
  const contentsHeader = sheet.getRangeByIndexes(contentsStart, 0, 1, 4);
  contentsHeader.format = {
    fill: palette.navy,
    font: {
      name: "Arial Unicode MS",
      bold: true,
      color: palette.white,
    },
    wrapText: true,
  };
  sheet
    .getRangeByIndexes(contentsStart + 1, 0, contentsRows.length, 4)
    .writeValues(contentsRows);
  sheet
    .getRangeByIndexes(contentsStart + 1, 0, contentsRows.length, 4)
    .format.font = { name: "Arial Unicode MS", size: 10, color: palette.text };

  sheet.getRange("A:A").format.columnWidth = 25;
  sheet.getRange("B:B").format.columnWidth = 14;
  sheet.getRange("C:C").format.columnWidth = 55;
  sheet.getRange("D:D").format.columnWidth = 60;
  sheet.getRange(`A3:D${contentsStart + contentsRows.length + 1}`).format.verticalAlignment =
    "top";
  sheet.getRange(`B3:D${contentsStart + contentsRows.length + 1}`).format.wrapText =
    true;
  sheet.freezePanes.freezeRows(1);
  return sheet;
}

function reorderObject(source, outputHeaders) {
  const result = {};
  for (const header of outputHeaders) result[header] = source[header] ?? null;
  return result;
}

async function buildSupplementaryTablesWorkbook() {
  const catalog = objectsFromMatrix(
    await readCsvMatrix("config/sample_group_catalog.csv"),
  );
  const paired30 = objectsFromMatrix(
    await readCsvMatrix(
      "results/tables/cell_type_kla_vs_reference_ddr_statistics_accession_only_paired_30.csv",
    ),
  );
  const referenceAudit = objectsFromMatrix(
    await readCsvMatrix("results/tables/strict_reference_material_identity_audit.csv"),
  );
  const orderedPaired30 = [...paired30].sort(
    (left, right) => Number(left.RowOrder) - Number(right.RowOrder),
  );

  const pairedByKey = new Map(
    paired30.map((row) => [keyOf(row.PXD, row.SampleGroup), row]),
  );
  const catalogByKey = new Map(
    catalog.map((row) => [
      keyOf(row["乳酸化PXD"], row["样本组"]),
      row,
    ]),
  );

  const s1Headers = [
    "RecordNo",
    "Category",
    "LactylomePXD",
    "SampleGroup",
    "BiologicalMaterial",
    "KlaLabelEn",
    "KlaEvidenceFile",
    "KlaProteinCount",
    "KlaDdrProteinCount",
    "KlaDdrFraction",
    "ReferencePXD",
    "ReferenceLabelEn",
    "ReferenceEvidenceFile",
    "ReferenceProteinCount",
    "ReferenceDdrProteinCount",
    "ReferenceDdrFraction",
    "DdrFractionPercentagePointDifference",
    "ReferenceMatchNote",
    "NamingBasis",
    "PublicationYear",
    "DOI",
    "DatasetTitle",
    "ProteomeXchangeURL",
  ];
  const s1Rows = orderedPaired30.map((paired, index) => {
    const pxd = cleanString(paired.PXD);
    const sampleGroup = cleanString(paired.SampleGroup);
    const catalogRow = catalogByKey.get(keyOf(pxd, sampleGroup));
    assert(catalogRow, `No catalogue row found for current group ${pxd}/${sampleGroup}.`);
    return {
      RecordNo: index + 1,
      Category: cleanString(paired.CategoryEn),
      LactylomePXD: pxd,
      SampleGroup: sampleGroup,
      BiologicalMaterial: cleanString(paired.BiologicalMaterial),
      KlaLabelEn: cleanString(paired.KlaLabelEn),
      KlaEvidenceFile: cleanString(paired.KlaEvidenceFile),
      KlaProteinCount: Number(paired.KlaProteinCount),
      KlaDdrProteinCount: Number(paired.KlaDdrProteinCount),
      KlaDdrFraction: Number(paired.KlaDdrFraction),
      ReferencePXD: cleanString(paired.ReferencePXD),
      ReferenceLabelEn: cleanString(paired.ReferenceLabelEn),
      ReferenceEvidenceFile: cleanString(paired.ReferenceEvidenceFile),
      ReferenceProteinCount: Number(paired.ReferenceProteinCount),
      ReferenceDdrProteinCount: Number(paired.ReferenceDdrProteinCount),
      ReferenceDdrFraction: Number(paired.ReferenceDdrFraction),
      DdrFractionPercentagePointDifference: Number(
        paired.DdrFractionPercentagePointDifference,
      ),
      ReferenceMatchNote: cleanString(paired.ReferenceMatchNote),
      NamingBasis: cleanString(paired.NamingBasis),
      PublicationYear: Number(catalogRow["发表年份"]) || null,
      DOI: cleanString(catalogRow["论文DOI"]),
      DatasetTitle: cleanString(catalogRow["数据集标题"]),
      ProteomeXchangeURL: cleanString(catalogRow["数据集链接"]),
    };
  });
  assert(s1Rows.length === 30, `Expected 30 S1 rows, found ${s1Rows.length}.`);

  const auditByKey = new Map(
    referenceAudit.map((row) => [keyOf(row.KlaPXD, row.SampleGroup), row]),
  );
  const s2Headers = [
    "RecordNo",
    "Category",
    ...Object.keys(referenceAudit[0]),
  ];
  const s2Rows = orderedPaired30.map((paired, index) => {
    const key = keyOf(paired.PXD, paired.SampleGroup);
    const row = auditByKey.get(key);
    assert(row, `No reference audit row found for current group ${key}.`);
    return {
      RecordNo: index + 1,
      Category: cleanString(paired.CategoryEn),
      ...row,
    };
  });
  assert(s2Rows.length === 30, `Expected 30 S2 rows, found ${s2Rows.length}.`);
  assert(
    new Set(s2Rows.map((row) => keyOf(row.KlaPXD, row.SampleGroup))).size === 30,
    "The current reference rows are not unique.",
  );

  const regulatorRaw = objectsFromMatrix(
    await readXlsxSheetMatrix(
      "data/identifier/乳酸化调控因子_Writer-Eraser-Reader.xlsx",
      "乳酸化调控因子",
    ),
  ).filter((row) => cleanString(row["基因"]));
  const regulatorMap = objectsFromMatrix(
    await readCsvMatrix("config/lactylation_regulator_uniprot_mapping.csv"),
  );
  const regulatorByGene = new Map(
    regulatorMap.map((row) => [cleanString(row.GeneSymbol), row]),
  );
  const s3Headers = [
    "RecordNo",
    "Role",
    "GeneSymbol",
    "BaseAccession",
    "References",
    "MappingSource",
    "MappingRetrievedDate",
  ];
  const s3Rows = regulatorRaw.map((row, index) => {
    const gene = cleanString(row["基因"]);
    const mapping = regulatorByGene.get(gene);
    assert(mapping, `No BaseAccession mapping found for regulator ${gene}.`);
    return {
      RecordNo: index + 1,
      Role: cleanString(row["分组"]),
      GeneSymbol: gene,
      BaseAccession: cleanString(mapping.BaseAccession),
      References: cleanString(row["参考文献"]),
      MappingSource: cleanString(mapping.MappingSource),
      MappingRetrievedDate: cleanString(mapping.RetrievedDate),
    };
  });
  assert(s3Rows.length === 50, `Expected 50 S3 role rows, found ${s3Rows.length}.`);

  const workbook = Workbook.create();
  const sheetMeta = [];

  const s1 = addMatrixSheet(
    workbook,
    "S1_Datasets",
    matrixFromObjects(s1Rows, s1Headers),
    "S1DatasetsTable",
  );
  sheetMeta.push(s1.meta);
  const s2 = addMatrixSheet(
    workbook,
    "S2_References",
    matrixFromObjects(s2Rows, s2Headers),
    "S2ReferencesTable",
  );
  sheetMeta.push(s2.meta);
  const s3 = addMatrixSheet(
    workbook,
    "S3_Regulators",
    matrixFromObjects(s3Rows, s3Headers),
    "S3RegulatorsTable",
  );
  sheetMeta.push(s3.meta);

  const vennAnalyses = [
    {
      key: "AllKla",
      label: "All Kla proteins",
      dir: "all_kla_four_class_venn",
    },
    {
      key: "KlaDDR",
      label: "Kla-DDR proteins",
      dir: "kla_ddr_four_class_venn",
    },
    {
      key: "RefProt",
      label: "Whole-proteome proteins",
      dir: "reference_proteome_four_class_venn",
    },
    {
      key: "RefDDR",
      label: "Whole-proteome DDR proteins",
      dir: "reference_proteome_ddr_four_class_venn",
    },
  ];
  const combinedSetCounts = [];
  const combinedRegionCounts = [];
  for (const analysis of vennAnalyses) {
    const root = `results/tables/four_class_venn/${analysis.dir}`;
    const setObjects = objectsFromMatrix(await readCsvMatrix(`${root}/set_counts.csv`));
    const regionObjects = objectsFromMatrix(
      await readCsvMatrix(`${root}/region_counts.csv`),
    );
    combinedSetCounts.push(
      ...setObjects.map((row) => ({
        Analysis: analysis.label,
        ...row,
      })),
    );
    combinedRegionCounts.push(
      ...regionObjects.map((row) => ({
        Analysis: analysis.label,
        ...row,
      })),
    );
    const membership = await addCsvSheet(
      workbook,
      `S4_${analysis.key}_Members`,
      `${root}/membership.csv`,
      `S4${analysis.key}MembersTable`,
    );
    sheetMeta.push(membership.meta);
  }
  assert(combinedSetCounts.length === 16, "Expected 16 combined Venn set counts.");
  assert(combinedRegionCounts.length === 60, "Expected 60 combined Venn regions.");
  const s4Set = addMatrixSheet(
    workbook,
    "S4_SetCounts",
    matrixFromObjects(combinedSetCounts, [
      "Analysis",
      ...Object.keys(combinedSetCounts[0]).filter((key) => key !== "Analysis"),
    ]),
    "S4SetCountsTable",
  );
  sheetMeta.push(s4Set.meta);
  const s4Regions = addMatrixSheet(
    workbook,
    "S4_RegionCounts",
    matrixFromObjects(combinedRegionCounts, [
      "Analysis",
      ...Object.keys(combinedRegionCounts[0]).filter((key) => key !== "Analysis"),
    ]),
    "S4RegionCountsTable",
  );
  sheetMeta.push(s4Regions.meta);

  for (const [sheetName, relativePath, tableName] of [
    [
      "S5_GO_Pathways",
      "results/tables/go_term_pathway_scoring_30groups/go_term_to_pathway_long.csv",
      "S5GOPathwaysTable",
    ],
    [
      "S5_DecisionAudit",
      "results/tables/go_term_pathway_scoring_30groups/go_term_decision_audit_2785.csv",
      "S5DecisionAuditTable",
    ],
    [
      "S5_SeedRules",
      "config/go_term_pathway_seed_rules.csv",
      "S5SeedRulesTable",
    ],
  ]) {
    const added = await addCsvSheet(workbook, sheetName, relativePath, tableName);
    sheetMeta.push(added.meta);
  }

  const revisedScoreRaw = objectsFromMatrix(
    await readXlsxSheetMatrix(
      "data/identifier/乳酸化DDR基因评分表_Revised_20260816.xlsx",
      "评分表",
    ),
  ).filter((row) => cleanString(row.BaseAccession));
  const currentMembership = objectsFromMatrix(
    await readCsvMatrix(
      "results/tables/four_class_venn/kla_ddr_four_class_venn/membership.csv",
    ),
  );
  const current399 = new Set(
    currentMembership.map((row) => cleanString(row.BaseAccession)),
  );
  const scoreByAccession = new Map(
    revisedScoreRaw.map((row) => [cleanString(row.BaseAccession), row]),
  );
  assert(
    [...current399].every((accession) => scoreByAccession.has(accession)),
    "The revised score table does not cover every protein in the current Kla-DDR union.",
  );

  const s6Headers = [
    "RecordNo",
    "BaseAccession",
    "GeneSymbol",
    "ProteinName",
    "BER",
    "NER",
    "MMR",
    "FA",
    "HR",
    "AEJ",
    "NHEJ",
    "ChromatinInteraction",
    "OtherSupport",
    "Note",
    "WorkbookScore",
    "RecalculatedScore",
    "ScoreMatches",
  ];
  const s6CurrentRows = [...current399]
    .sort((left, right) => left.localeCompare(right))
    .map((accession, index) => {
    const row = scoreByAccession.get(accession);
    const recalc =
      Number(row.BER) +
      2 * Number(row.NER) +
      3 * Number(row.MMR) +
      4 * Number(row.FA) +
      5 * Number(row.HR) +
      6 * Number(row.AEJ) +
      7 * Number(row.NHEJ);
    assert(
      recalc === Number(row.score),
      `Weighted score mismatch for ${cleanString(row.BaseAccession)}.`,
    );
    return {
      RecordNo: index + 1,
      BaseAccession: cleanString(row.BaseAccession),
      GeneSymbol: cleanString(row.GeneSymbol),
      ProteinName: cleanString(row.ProteinName),
      BER: Number(row.BER),
      NER: Number(row.NER),
      MMR: Number(row.MMR),
      FA: Number(row.FA),
      HR: Number(row.HR),
      AEJ: Number(row.AEJ),
      NHEJ: Number(row.NHEJ),
      ChromatinInteraction: Number(row["Chromatin Interaction"]),
      OtherSupport: Number(
        row["Others (Transcription, RNA processing and proteostasis)"],
      ),
      Note: cleanString(row.Note),
      WorkbookScore: Number(row.score),
      RecalculatedScore: null,
      ScoreMatches: null,
    };
  });
  assert(s6CurrentRows.length === 399, "Expected 399 current signed-score rows.");
  const s6Current = addMatrixSheet(
    workbook,
    "S6_Signed399",
    matrixFromObjects(s6CurrentRows, s6Headers),
    "S6Signed399Table",
  );
  const s6CurrentEnd = s6CurrentRows.length + 1;
  s6Current.sheet.getRange("P2").formulas = [
    ["=E2+2*F2+3*G2+4*H2+5*I2+6*J2+7*K2"],
  ];
  s6Current.sheet.getRange(`P2:P${s6CurrentEnd}`).fillDown();
  s6Current.sheet.getRange("Q2").formulas = [["=O2=P2"]];
  s6Current.sheet.getRange(`Q2:Q${s6CurrentEnd}`).fillDown();
  sheetMeta.push(s6Current.meta);

  const contents = [
    ["S1_Datasets", 30, "Current dataset catalogue and quantitative summary", "results/tables/cell_type_kla_vs_reference_ddr_statistics_accession_only_paired_30.csv"],
    ["S2_References", 30, "Current whole-proteome reference pairing and state-matching audit", "results/tables/strict_reference_material_identity_audit.csv filtered to the current 30 groups"],
    ["S3_Regulators", 50, "Writer, Eraser and Reader role entries with BaseAccession mappings", "data/identifier/乳酸化调控因子_Writer-Eraser-Reader.xlsx"],
    ["S4_SetCounts", 16, "Four biological categories across four Venn analyses", "results/tables/four_class_venn/*/set_counts.csv"],
    ["S4_RegionCounts", 60, "All 15 exact logical regions across four Venn analyses", "results/tables/four_class_venn/*/region_counts.csv"],
    ["S4_*_Members", "variable", "Exact BaseAccession membership for each Venn analysis", "results/tables/four_class_venn/*/membership.csv"],
    ["S5_GO_Pathways", 2793, "Long GO-term-to-pathway decisions; multi-pathway terms retained", "results/tables/go_term_pathway_scoring_30groups/go_term_to_pathway_long.csv"],
    ["S5_DecisionAudit", 2785, "Exhaustive one-row-per-direct-term audit", "results/tables/go_term_pathway_scoring_30groups/go_term_decision_audit_2785.csv"],
    ["S5_SeedRules", "current", "Curated exact/descendant seed rules", "config/go_term_pathway_seed_rules.csv"],
    ["S6_Signed399", 399, "Revised manually curated signed pathway annotations for the current Kla-DDR union", "data/identifier/乳酸化DDR基因评分表_Revised_20260816.xlsx filtered by current Kla-DDR membership"],
  ];
  addReadmeSheet(
    workbook,
    "Draft Supplementary Tables S1-S6",
    [
      ["Analysis scope", "30 current Kla-reference groups"],
      ["Biological classes", "9 non-tumor tissues; 2 tumor tissues; 12 cancer cell lines; 7 normal cell lines"],
      ["Current Kla-DDR union", "399 unique UniProt BaseAccessions"],
      ["Analytical identifier", "Isoform-stripped UniProt BaseAccession; gene symbols and protein names are display/audit fields only"],
      ["Version note", "S1-S5 support the direct-GO-term primary analysis. S6 contains the revised signed pathway annotations described as an alternative method."],
    ],
    contents,
  );

  return { workbook, sheetMeta };
}

async function buildDataS1Workbook() {
  const paired30 = objectsFromMatrix(
    await readCsvMatrix(
      "results/tables/cell_type_kla_vs_reference_ddr_statistics_accession_only_paired_30.csv",
    ),
  );
  const currentKeys = new Set(
    paired30.map((row) => keyOf(row.PXD, row.SampleGroup)),
  );
  const pairedByKey = new Map(
    paired30.map((row) => [keyOf(row.PXD, row.SampleGroup), row]),
  );
  const currentRecordNoByKey = new Map(
    [...paired30]
      .sort((left, right) => Number(left.RowOrder) - Number(right.RowOrder))
      .map((row, index) => [keyOf(row.PXD, row.SampleGroup), index + 1]),
  );
  const allMembership = objectsFromMatrix(
    await readCsvMatrix(
      "work/intermediate/expanded_ddr_by_accession/kla_proteins_by_sample_group.csv",
    ),
  );
  const currentMembership = allMembership
    .filter((row) => currentKeys.has(keyOf(row.PXD, row.SampleGroup)))
    .map((row) => {
      const paired = pairedByKey.get(keyOf(row.PXD, row.SampleGroup));
      return {
        RecordNo: currentRecordNoByKey.get(keyOf(row.PXD, row.SampleGroup)),
        Category: cleanString(paired.CategoryEn),
        PXD: cleanString(row.PXD),
        SampleGroup: cleanString(row.SampleGroup),
        BaseAccession: cleanString(row.BaseAccession),
        IsDdr: isTrue(row.IsDdr),
      };
    })
    .sort((left, right) =>
      left.RecordNo - right.RecordNo ||
      left.BaseAccession.localeCompare(right.BaseAccession),
    );
  const observedKeys = new Set(
    currentMembership.map((row) => keyOf(row.PXD, row.SampleGroup)),
  );
  assert(currentKeys.size === 30, "The current dataset scope is not 30 groups.");
  assert(
    observedKeys.size === 30 &&
      [...currentKeys].every((key) => observedKeys.has(key)),
    "Current Kla protein membership does not cover all 30 groups.",
  );
  const currentDdrMembership = currentMembership.filter((row) => row.IsDdr);
  assert(
    new Set(currentDdrMembership.map((row) => row.BaseAccession)).size === 399,
    "Current Kla-DDR membership does not contain 399 unique proteins.",
  );

  const summaryHeaders = [
    "RecordNo",
    "Category",
    "PXD",
    "SampleGroup",
    "BiologicalMaterial",
    "KlaProteinCount",
    "KlaDdrProteinCount",
    "KlaDdrFraction",
    "KlaEvidenceFile",
  ];
  const summaryRows = [...paired30]
    .sort((left, right) => Number(left.RowOrder) - Number(right.RowOrder))
    .map((row, index) => ({
      RecordNo: index + 1,
      Category: cleanString(row.CategoryEn),
      PXD: cleanString(row.PXD),
      SampleGroup: cleanString(row.SampleGroup),
      BiologicalMaterial: cleanString(row.BiologicalMaterial),
      KlaProteinCount: Number(row.KlaProteinCount),
      KlaDdrProteinCount: Number(row.KlaDdrProteinCount),
      KlaDdrFraction: Number(row.KlaDdrFraction),
      KlaEvidenceFile: cleanString(row.KlaEvidenceFile),
    }));

  const workbook = Workbook.create();
  const sheetMeta = [];
  const groupSummary = addMatrixSheet(
    workbook,
    "Group_Summary",
    matrixFromObjects(summaryRows, summaryHeaders),
    "Current30GroupSummaryTable",
  );
  sheetMeta.push(groupSummary.meta);
  const allProteins = addMatrixSheet(
    workbook,
    "Kla_Protein_Membership",
    matrixFromObjects(currentMembership, [
      "RecordNo",
      "Category",
      "PXD",
      "SampleGroup",
      "BaseAccession",
      "IsDdr",
    ]),
    "Current30KlaProteinMembershipTable",
  );
  sheetMeta.push(allProteins.meta);
  const ddrProteins = addMatrixSheet(
    workbook,
    "Kla_DDR_Membership",
    matrixFromObjects(currentDdrMembership, [
      "RecordNo",
      "Category",
      "PXD",
      "SampleGroup",
      "BaseAccession",
      "IsDdr",
    ]),
    "Current30KlaDDRMembershipTable",
  );
  sheetMeta.push(ddrProteins.meta);
  addReadmeSheet(
    workbook,
    "Draft Supplementary Data S1: current Kla protein membership",
    [
      ["Analysis scope", "30 current Kla-reference groups only"],
      ["Analytical identifier", "Isoform-stripped UniProt BaseAccession"],
      ["Current Kla-DDR union", "399 unique BaseAccessions"],
      ["Interpretation", "Each membership row records one protein in one current biological sample group. IsDdr indicates intersection with the human GO-DDR set."],
    ],
    [
      ["Group_Summary", groupSummary.meta.rows, "Current Kla and Kla-DDR counts for each of the 30 groups", "results/tables/cell_type_kla_vs_reference_ddr_statistics_accession_only_paired_30.csv"],
      ["Kla_Protein_Membership", allProteins.meta.rows, "Current group-by-BaseAccession Kla membership", "work/intermediate/expanded_ddr_by_accession/kla_proteins_by_sample_group.csv filtered to current groups"],
      ["Kla_DDR_Membership", ddrProteins.meta.rows, "Current group-by-BaseAccession Kla-DDR membership", "Kla_Protein_Membership filtered by IsDdr"],
    ],
  );
  return { workbook, sheetMeta };
}

async function buildDataS2Workbook() {
  const workbook = Workbook.create();
  const matrix = await readTsvMatrix("data/annotations/GO-repair+damage(human).tsv");
  const added = addMatrixSheet(
    workbook,
    "Human_GO_DDR",
    matrix,
    "HumanGODDRTable",
  );
  assert(added.meta.rows === 5111, `Expected 5,111 GO-DDR rows, found ${added.meta.rows}.`);
  addReadmeSheet(
    workbook,
    "Draft Supplementary Data S2: human GO-DDR annotations",
    [
      ["Taxon", "Homo sapiens (taxon 9606)"],
      ["Primary filtering", "Records with the GO qualifier NOT are excluded during the analysis; evidence codes are retained without primary evidence-code restriction."],
      ["Analytical identifier", "UniProtKB accessions normalized to isoform-stripped BaseAccession during intersection"],
    ],
    [
      ["Human_GO_DDR", added.meta.rows, "Human DNA-damage-response GO annotations used to define the DDR set", "data/annotations/GO-repair+damage(human).tsv"],
    ],
  );
  return { workbook, sheetMeta: [added.meta] };
}

async function buildDataS3Workbook() {
  const workbook = Workbook.create();
  const sheetMeta = [];
  const sources = [
    [
      "Direct_GO_Annotations",
      "results/tables/go_term_pathway_scoring_30groups/direct_go_annotations_399proteins.csv",
      "DirectGOAnnotationsTable",
    ],
    [
      "Protein_GO_Pathway",
      "results/tables/go_term_pathway_scoring_30groups/protein_go_term_pathway_long.csv",
      "ProteinGOPathwayTable",
    ],
    [
      "Pathway_Counts_Long",
      "results/tables/go_term_pathway_scoring_30groups/protein_pathway_direct_term_counts_long.csv",
      "PathwayCountsLongTable",
    ],
    [
      "Pathway_Count_Matrix",
      "results/tables/go_term_pathway_scoring_30groups/protein_pathway_direct_term_count_matrix.csv",
      "PathwayCountMatrixTable",
    ],
    [
      "Pathway_Binary_Matrix",
      "results/tables/go_term_pathway_scoring_30groups/protein_seven_pathway_binary_matrix.csv",
      "PathwayBinaryMatrixTable",
    ],
  ];
  for (const [sheetName, source, tableName] of sources) {
    const added = await addCsvSheet(workbook, sheetName, source, tableName);
    sheetMeta.push(added.meta);
  }
  assert(sheetMeta[0].rows === 10605, "Direct GO annotation rows are not 10,605.");
  assert(sheetMeta[3].rows === 399, "Pathway count matrix is not 399 proteins.");
  assert(sheetMeta[4].rows === 399, "Pathway binary matrix is not 399 proteins.");
  addReadmeSheet(
    workbook,
    "Draft Supplementary Data S3: direct GO annotations and pathway matrices",
    [
      ["Analysis scope", "399 unique BaseAccessions in the current 30-group Kla-DDR union"],
      ["Direct annotations", "10,605 unique protein-GO pairs and 2,785 unique direct GO terms from UniProt release 2026_02"],
      ["Pathways", "BER, NER, MMR, FA, HR, NHEJ and AEJ; unmatched direct terms are retained as Others"],
      ["Analytical identifier", "Isoform-stripped UniProt BaseAccession"],
    ],
    sources.map(([sheetName, source], index) => [
      sheetName,
      sheetMeta[index].rows,
      [
        "Direct protein-GO pairs",
        "Protein-GO-pathway long table",
        "Direct term counts by protein and pathway",
        "Wide direct-term count matrix",
        "Wide seven-pathway presence/absence matrix",
      ][index],
      source,
    ]),
  );
  return { workbook, sheetMeta };
}

async function verifyWorkbook(workbook, workbookLabel, sheetMeta) {
  const sheetInspection = await workbook.inspect({
    kind: "sheet",
    include: "id,name",
    maxChars: 10000,
  });
  const sheetNamesText = sheetInspection.ndjson;
  for (const meta of sheetMeta) {
    assert(
      sheetNamesText.includes(meta.name),
      `${workbookLabel} is missing worksheet ${meta.name}.`,
    );
  }

  const errors = await workbook.inspect({
    kind: "match",
    searchTerm: "#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A",
    options: { useRegex: true, maxResults: 300 },
    summary: `${workbookLabel} formula error scan`,
    maxChars: 5000,
  });
  assert(
    !/#REF!|#DIV\/0!|#VALUE!|#NAME\?|#N\/A/.test(errors.ndjson),
    `${workbookLabel} contains formula errors.`,
  );

  for (const meta of [{ name: "README", rows: 15, columns: 4 }, ...sheetMeta]) {
    const previewRows = Math.max(2, Math.min(meta.rows + 1, 14));
    const previewColumns = Math.max(2, Math.min(meta.columns, 12));
    const range = `A1:${excelColumn(previewColumns - 1)}${previewRows}`;
    const preview = await workbook.render({
      sheetName: meta.name,
      range,
      scale: 1,
      format: "png",
    });
    const safeLabel = workbookLabel.replace(/[^A-Za-z0-9_-]+/g, "_");
    const safeSheet = meta.name.replace(/[^A-Za-z0-9_-]+/g, "_");
    await fs.writeFile(
      path.join(previewDir, `${safeLabel}__${safeSheet}.png`),
      new Uint8Array(await preview.arrayBuffer()),
    );
  }
}

async function exportWorkbook(workbook, outputFilename) {
  const externalPath = path.join(externalOutputDir, outputFilename);
  const projectPath = path.join(projectOutputDir, outputFilename);
  const output = await SpreadsheetFile.exportXlsx(workbook);
  await output.save(externalPath);
  await fs.copyFile(externalPath, projectPath);
  return { externalPath, projectPath };
}

async function verifyExportedWorkbook(filePath, workbookLabel, sheetMeta) {
  const imported = await SpreadsheetFile.importXlsx(await FileBlob.load(filePath));
  const sheets = await imported.inspect({
    kind: "sheet",
    include: "id,name",
    maxChars: 10000,
  });
  for (const meta of sheetMeta) {
    assert(
      sheets.ndjson.includes(meta.name),
      `${workbookLabel} export is missing worksheet ${meta.name}.`,
    );
  }
  const errors = await imported.inspect({
    kind: "match",
    searchTerm: "#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A",
    options: { useRegex: true, maxResults: 300 },
    summary: `${workbookLabel} exported formula error scan`,
    maxChars: 5000,
  });
  assert(
    !/#REF!|#DIV\/0!|#VALUE!|#NAME\?|#N\/A/.test(errors.ndjson),
    `${workbookLabel} export contains formula errors.`,
  );
}

async function sha256(filePath) {
  const bytes = await fs.readFile(filePath);
  return crypto.createHash("sha256").update(bytes).digest("hex");
}

async function main() {
  await fs.mkdir(externalOutputDir, { recursive: true });
  await fs.mkdir(projectOutputDir, { recursive: true });
  await fs.mkdir(previewDir, { recursive: true });

  const builds = [
    ["tables", "Supplementary_Tables", await buildSupplementaryTablesWorkbook()],
    ["dataS1", "Supplementary_Data_S1", await buildDataS1Workbook()],
    ["dataS2", "Supplementary_Data_S2", await buildDataS2Workbook()],
    ["dataS3", "Supplementary_Data_S3", await buildDataS3Workbook()],
  ];

  const manifest = [];
  for (const [outputKey, label, built] of builds) {
    await verifyWorkbook(built.workbook, label, built.sheetMeta);
    const exported = await exportWorkbook(built.workbook, outputFiles[outputKey]);
    await verifyExportedWorkbook(exported.projectPath, label, built.sheetMeta);
    manifest.push({
      File: outputFiles[outputKey],
      ProjectPath: path.relative(projectRoot, exported.projectPath),
      ExternalOutputPath: exported.externalPath,
      SHA256: await sha256(exported.projectPath),
      WorksheetCount: built.sheetMeta.length + 1,
      DataRows: built.sheetMeta.reduce((sum, item) => sum + item.rows, 0),
      ReleaseDate: releaseDate,
    });
  }

  const manifestHeaders = [
    "File",
    "ProjectPath",
    "ExternalOutputPath",
    "SHA256",
    "WorksheetCount",
    "DataRows",
    "ReleaseDate",
  ];
  const csv = [
    manifestHeaders.join(","),
    ...manifest.map((row) =>
      manifestHeaders
        .map((header) => {
          const value = row[header] ?? "";
          const text = String(value);
          return `"${text.replaceAll('"', '""')}"`;
        })
        .join(","),
    ),
  ].join("\n");
  await fs.writeFile(
    path.join(projectOutputDir, "supplementary_workbook_manifest.csv"),
    `${csv}\n`,
    "utf8",
  );

  console.log(JSON.stringify({ projectRoot, externalOutputDir, manifest }, null, 2));
}

await main();
