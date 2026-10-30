import fs from 'node:fs/promises';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
import { Workbook, SpreadsheetFile } from '@oai/artifact-tool';

const root = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');
const outputDir = process.argv[2] || path.join(root, 'outputs/01a18b74-d839-7a53-b3c3-3c45f301073e');
const payload = JSON.parse(await fs.readFile(path.join(outputDir, 'prism_points.json'), 'utf8'));
const workbook = Workbook.create();
const previews = path.join(outputDir, 'preview_plain');
await fs.mkdir(previews, { recursive: true });

function columnName(n) {
  let result = '';
  while (n) { n--; result = String.fromCharCode(65 + n % 26) + result; n = Math.floor(n / 26); }
  return result;
}

for (const panel of payload.panels) {
  const sheet = workbook.worksheets.add(panel.name);
  const columns = panel.matrix[0].length;
  const rows = panel.matrix.length;
  const end = columnName(columns);
  sheet.getRange(`A1:${end}${rows}`).values = panel.matrix;
}

workbook.recalculate();
for (const detail of payload.details) {
  const value = workbook.worksheets.getItem(detail.Sheet).getRange(detail.Cell).values[0][0];
  if (value !== detail.Value) throw new Error(`Point mismatch: ${detail.Sheet}!${detail.Cell}`);
}
const summary = await workbook.inspect({ kind: 'table', range: "'a_DDR比例'!A1:F3", include: 'values', tableMaxRows: 3, tableMaxCols: 6, maxChars: 1500 });
console.log(summary.ndjson);
const errors = await workbook.inspect({ kind: 'match', searchTerm: '#REF!|#DIV/0!|#VALUE!|#NAME\\?|#NUM!|#NULL!|#SPILL!|#CALC!', options: { useRegex: true, maxResults: 20 }, maxChars: 1500 });
console.log(errors.ndjson);
for (const panel of payload.panels) {
  const range = panel.type === 'Column' ? 'A1:D13' : 'A1:H3';
  const preview = await workbook.render({ sheetName: panel.name, range, scale: 1.5, format: 'png' });
  await fs.writeFile(path.join(previews, `${panel.name}.png`), new Uint8Array(await preview.arrayBuffer()));
}
const output = await SpreadsheetFile.exportXlsx(workbook);
const outputPath = path.join(outputDir, 'GraphPad_图a_b及七通路_逐点数据_20261030.xlsx');
await output.save(outputPath);
console.log(JSON.stringify({ outputPath, sheets: payload.panels.length, pointsVerified: payload.details.length }));
