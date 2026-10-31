import fs from 'node:fs/promises';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
import { Workbook, SpreadsheetFile } from '@oai/artifact-tool';

const root = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');
const outputDir = process.argv[2] || path.join(root, 'outputs/01a18b74-d839-7a53-b3c3-3c45f301073e');
const payload = JSON.parse(await fs.readFile(path.join(outputDir, 'prism_points.json'), 'utf8'));
const workbook = Workbook.create();
const previews = path.join(outputDir, 'preview_adjacent_four_columns');
await fs.mkdir(previews, { recursive: true });
const categoryLabels = ['非肿瘤组织', '肿瘤组织', '正常细胞系', '癌细胞系'];
const names = new Map(payload.panels.map(panel => [panel.name,
  panel.name.startsWith('a_') ? '示例a_DDR比例' : panel.name.startsWith('b_') ? '示例b_增殖得分' : panel.name]));
const targets = new Map();

function columnName(n) {
  let result = '';
  while (n) { n--; result = String.fromCharCode(65 + n % 26) + result; n = Math.floor(n / 26); }
  return result;
}

for (const panel of payload.panels) {
  const sheet = workbook.worksheets.add(names.get(panel.name));
  let matrix;
  if (panel.type === 'Column') {
    matrix = [categoryLabels, ...panel.matrix.slice(1)];
    for (const detail of payload.details.filter(value => value.Sheet === panel.name)) {
      targets.set(`${detail.Sheet}!${detail.Cell}`, { sheet: sheet.name, cell: detail.Cell });
    }
  } else {
    matrix = [[panel.name.startsWith('a_') ? '检测类型' : '方向', ...categoryLabels]];
    for (let seriesIndex = 1; seriesIndex < panel.matrix.length; seriesIndex++) {
      const series = panel.matrix[seriesIndex][0];
      const count = Math.max(...panel.counts[series]);
      const firstRow = matrix.length + 1;
      for (let replicate = 0; replicate < count; replicate++) {
        matrix.push([series, ...categoryLabels.map((_, categoryIndex) =>
          panel.matrix[seriesIndex][1 + categoryIndex * panel.replicates + replicate])]);
      }
      for (const detail of payload.details.filter(value => value.Sheet === panel.name && value.Series === series)) {
        const sourceColumn = detail.Cell.match(/[A-Z]+/)[0].split('').reduce((n, c) => n * 26 + c.charCodeAt(0) - 64, 0);
        const categoryIndex = Math.floor((sourceColumn - 2) / panel.replicates);
        targets.set(`${detail.Sheet}!${detail.Cell}`, { sheet: sheet.name,
          cell: `${columnName(categoryIndex + 2)}${firstRow + detail.Replicate - 1}` });
      }
    }
  }
  const end = columnName(matrix[0].length);
  sheet.getRange(`A1:${end}${matrix.length}`).values = matrix;
  sheet.getRange(`A1:${end}${matrix.length}`).format.columnWidth = 24;
}

workbook.recalculate();
for (const detail of payload.details) {
  const target = targets.get(`${detail.Sheet}!${detail.Cell}`);
  const value = workbook.worksheets.getItem(target.sheet).getRange(target.cell).values[0][0];
  if (value !== detail.Value) throw new Error(`Point mismatch: ${detail.Sheet}!${detail.Cell}`);
}
const summary = await workbook.inspect({ kind: 'table', range: "'示例a_DDR比例'!A1:E5", include: 'values', tableMaxRows: 5, tableMaxCols: 5, maxChars: 1500 });
console.log(summary.ndjson);
const errors = await workbook.inspect({ kind: 'match', searchTerm: '#REF!|#DIV/0!|#VALUE!|#NAME\\?|#NUM!|#NULL!|#SPILL!|#CALC!', options: { useRegex: true, maxResults: 20 }, maxChars: 1500 });
console.log(errors.ndjson);
for (const panel of payload.panels) {
  const range = panel.type === 'Column' ? 'A1:D13' : 'A1:E8';
  const preview = await workbook.render({ sheetName: names.get(panel.name), range, scale: 1.2, format: 'png' });
  await fs.writeFile(path.join(previews, `${panel.name}.png`), new Uint8Array(await preview.arrayBuffer()));
}
const output = await SpreadsheetFile.exportXlsx(workbook);
const outputPath = path.join(outputDir, '四类样品_相邻四列_纯数据_20261031.xlsx');
await output.save(outputPath);
console.log(JSON.stringify({ outputPath, sheets: payload.panels.length, pointsVerified: payload.details.length }));
