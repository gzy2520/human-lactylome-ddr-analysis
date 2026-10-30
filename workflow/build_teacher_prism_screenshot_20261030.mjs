import fs from 'node:fs/promises';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
import { Workbook, SpreadsheetFile } from '@oai/artifact-tool';

const root = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');
const outputDir = process.argv[2] || path.join(root, 'outputs/01a18b74-d839-7a53-b3c3-3c45f301073e');
const payload = JSON.parse(await fs.readFile(path.join(outputDir, 'prism_points.json'), 'utf8'));
const workbook = Workbook.create();
const previews = path.join(outputDir, 'preview_four_categories');
await fs.mkdir(previews, { recursive: true });
const categoryLabels = ['非肿瘤组织', '肿瘤组织', '正常细胞系', '癌细胞系'];
const names = new Map(payload.panels.map(panel => [panel.name,
  panel.name.startsWith('a_') ? '示例a_DDR比例' : panel.name.startsWith('b_') ? '示例b_增殖得分' : panel.name]));

function columnName(n) {
  let result = '';
  while (n) { n--; result = String.fromCharCode(65 + n % 26) + result; n = Math.floor(n / 26); }
  return result;
}

for (const panel of payload.panels) {
  const sheet = workbook.worksheets.add(names.get(panel.name));
  const matrix = [panel.type === 'Column' ? categoryLabels :
    [null, ...categoryLabels.flatMap(label => Array(panel.replicates).fill(label))],
    ...panel.matrix.slice(1)];
  const end = columnName(matrix[0].length);
  sheet.getRange(`A1:${end}${matrix.length}`).values = matrix;
  sheet.getRange(`A1:${end}${matrix.length}`).format.columnWidth = 24;
}

workbook.recalculate();
for (const detail of payload.details) {
  const value = workbook.worksheets.getItem(names.get(detail.Sheet)).getRange(detail.Cell).values[0][0];
  if (value !== detail.Value) throw new Error(`Point mismatch: ${detail.Sheet}!${detail.Cell}`);
}
const summary = await workbook.inspect({ kind: 'table', range: "'示例a_DDR比例'!A1:E3", include: 'values', tableMaxRows: 3, tableMaxCols: 5, maxChars: 1500 });
console.log(summary.ndjson);
const errors = await workbook.inspect({ kind: 'match', searchTerm: '#REF!|#DIV/0!|#VALUE!|#NAME\\?|#NUM!|#NULL!|#SPILL!|#CALC!', options: { useRegex: true, maxResults: 20 }, maxChars: 1500 });
console.log(errors.ndjson);
for (const panel of payload.panels) {
  const range = panel.type === 'Column' ? 'A1:D13' : 'A1:E3';
  const preview = await workbook.render({ sheetName: names.get(panel.name), range, scale: 1.2, format: 'png' });
  await fs.writeFile(path.join(previews, `${panel.name}.png`), new Uint8Array(await preview.arrayBuffer()));
}
const output = await SpreadsheetFile.exportXlsx(workbook);
const outputPath = path.join(outputDir, 'GraphPad_四类样品_纯数据_20261030.xlsx');
await output.save(outputPath);
console.log(JSON.stringify({ outputPath, sheets: payload.panels.length, pointsVerified: payload.details.length }));
