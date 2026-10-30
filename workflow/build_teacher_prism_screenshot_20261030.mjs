import fs from 'node:fs/promises';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
import { Workbook, SpreadsheetFile } from '@oai/artifact-tool';

const root = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');
const outputDir = process.argv[2] || path.join(root, 'outputs/01a18b74-d839-7a53-b3c3-3c45f301073e');
const payload = JSON.parse(await fs.readFile(path.join(outputDir, 'prism_points.json'), 'utf8'));
const workbook = Workbook.create();
const previews = path.join(outputDir, 'preview');
await fs.mkdir(previews, { recursive: true });

function columnName(n) {
  let result = '';
  while (n) { n--; result = String.fromCharCode(65 + n % 26) + result; n = Math.floor(n / 26); }
  return result;
}

function plain(sheet, range) {
  sheet.showGridLines = true;
  sheet.getRange(range).format.font = { name: 'Arial', size: 11, color: '#000000' };
  sheet.getRange(range).format.verticalAlignment = 'center';
  sheet.getRange(range).format.rowHeight = 22;
}

for (const panel of payload.panels) {
  const sheet = workbook.worksheets.add(panel.name);
  const columns = panel.matrix[0].length;
  const rows = panel.matrix.length;
  const end = columnName(columns);
  sheet.getRange(`A1:${end}${rows}`).values = panel.matrix;
  plain(sheet, `A1:${end}${rows}`);
  sheet.getRange(`A1:${end}1`).format.font.bold = true;
  sheet.getRange(`A1:${end}1`).format.wrapText = true;
  sheet.getRange(`A1:${end}1`).format.horizontalAlignment = 'center';
  sheet.getRange(`A1:${end}1`).format.rowHeight = 44;
  sheet.getRange(`A1:${end}1`).format.borders = { bottom: { style: 'thin', color: '#000000' } };
  if (panel.type === 'Grouped') {
    sheet.getRange(`B1:${end}${rows}`).format.columnWidth = 16;
    sheet.getRange(`A1:A${rows}`).format.columnWidth = 23;
    sheet.getRange(`A2:A${rows}`).format.font.bold = true;
    sheet.getRange(`B2:${end}${rows}`).setNumberFormat('0.000000');
    sheet.freezePanes.freezeColumns(1);
    for (const categoryIndex of [1, 2, 3]) {
      const start = columnName(2 + categoryIndex * panel.replicates);
      sheet.getRange(`${start}1:${start}${rows}`).format.borders = { left: { style: 'thin', color: '#000000' } };
    }
  } else {
    sheet.getRange(`A1:${end}${rows}`).format.columnWidth = 24;
    sheet.getRange(`A2:${end}${rows}`).setNumberFormat('0.000000');
    sheet.freezePanes.freezeRows(1);
  }
}

const detailSheet = workbook.worksheets.add('逐点明细');
const detailKeys = Object.keys(payload.details[0]);
const detailMatrix = [detailKeys, ...payload.details.map(row => detailKeys.map(key => row[key] ?? null))];
const detailEnd = columnName(detailKeys.length);
detailSheet.getRange(`A1:${detailEnd}${detailMatrix.length}`).values = detailMatrix;
plain(detailSheet, `A1:${detailEnd}${detailMatrix.length}`);
detailSheet.getRange(`A1:${detailEnd}1`).format.font.bold = true;
detailSheet.getRange(`A1:${detailEnd}1`).format.rowHeight = 36;
detailSheet.getRange(`A1:${detailEnd}1`).format.wrapText = true;
detailSheet.getRange(`A1:${detailEnd}${detailMatrix.length}`).format.columnWidth = 20;
detailSheet.getRange(`C1:D${detailMatrix.length}`).format.columnWidth = 23;
detailSheet.getRange(`G1:L${detailMatrix.length}`).format.columnWidth = 32;
detailSheet.getRange(`G1:L${detailMatrix.length}`).format.wrapText = true;
detailSheet.getRange(`T1:T${detailMatrix.length}`).format.columnWidth = 85;
detailSheet.getRange(`T1:T${detailMatrix.length}`).format.wrapText = true;
detailSheet.getRange(`A2:${detailEnd}${detailMatrix.length}`).format.autofitRows();
detailSheet.getRange(`P2:P${detailMatrix.length}`).setNumberFormat('0.000000');
detailSheet.freezePanes.freezeRows(1);
detailSheet.freezePanes.freezeColumns(2);

const guide = workbook.worksheets.add('使用说明');
const guideRows = [
  ['项目', '说明', '各组点数（非肿瘤组织 / 肿瘤组织 / 正常细胞系 / 癌细胞系）'],
  ['数据版本', '对应本次提供的两张组合图，使用历史作图数值。', null],
  ['Grouped 表', '四类样品为列组；列组内连续子列为 Y1、Y2……；第一列为行名。', null],
  ['Prism 导入', '新建 Grouped 表，选择将重复观测值录入并排子列；设4个列组。', null],
  ['复制范围', 'Grouped表第1行为列标题，第1列为行标题，第2–3行为数值；空白保持为空。', null],
  ['百分比', '图a及通路表均为0–100的百分数，例如4.967603表示4.967603%。', null],
  ['a_DDR比例', 'Grouped；每个列组107个子列；Whole与Kla分别排序，编号分别对应各自明细。', 'Whole: 60 / 107 / 13 / 20; Kla: 31 / 11 / 9 / 21'],
  ['b_增殖得分', 'Column；四列分别输入四类样品；直接使用已计算的log2得分。', '9 / 3 / 7 / 12'],
  ['d_BER至j_NHEJ', 'Grouped；每个列组31个子列；同一子列的Pro和Inh对应同一样本。', '每个方向：31 / 11 / 9 / 21'],
  ['图a计算', 'DDR蛋白数 / 该来源观测的全部检出蛋白数 × 100。', null],
  ['通路计算', '对应通路Pro或Inh蛋白数 / 该样本Kla-DDR蛋白数 × 100。', null],
  ['图b计算', '每RNA样本取196个G2M_CHECKPOINT基因的log2(qsmooth TPM+0.5)均值，再按历史GroupID取样本中位数。', null],
  ['图b条目', '保留原图31个GroupID；共享ReferenceKey的条目可在逐点明细中核查。', null],
  ['样本查找', '逐点明细的Sheet及Cell定位数值；Replicate为对应类别内的子列序号。', null],
  ['精度', '单元格显示6位小数，文件保留源数值精度；0为源表实际值。', null],
  ['数据来源', '以下列出历史作图源表及SHA256。', null],
  ...payload.sources.map(source => [source.key, source.file, source.sha256]),
];
guide.getRange(`A1:C${guideRows.length}`).values = guideRows;
plain(guide, `A1:C${guideRows.length}`);
guide.getRange(`A1:A${guideRows.length}`).format.columnWidth = 22;
guide.getRange(`B1:B${guideRows.length}`).format.columnWidth = 85;
guide.getRange(`C1:C${guideRows.length}`).format.columnWidth = 68;
guide.getRange(`A1:C${guideRows.length}`).format.wrapText = true;
guide.getRange(`A1:C${guideRows.length}`).format.autofitRows();
guide.getRange('A1:C1').format.font.bold = true;

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
for (const [name, range] of [['逐点明细', 'A1:F8'], ['使用说明', `A1:C${guideRows.length}`]]) {
  const preview = await workbook.render({ sheetName: name, range, scale: 1.2, format: 'png' });
  await fs.writeFile(path.join(previews, `${name}.png`), new Uint8Array(await preview.arrayBuffer()));
}
const output = await SpreadsheetFile.exportXlsx(workbook);
const outputPath = path.join(outputDir, 'GraphPad_图a_b及七通路_逐点数据_20261030.xlsx');
await output.save(outputPath);
console.log(JSON.stringify({ outputPath, sheets: payload.panels.length + 2, pointsVerified: payload.details.length }));
