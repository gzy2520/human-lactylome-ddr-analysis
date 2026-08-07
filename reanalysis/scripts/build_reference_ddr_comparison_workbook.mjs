import fs from "node:fs/promises";
import path from "node:path";
import { SpreadsheetFile, Workbook } from "@oai/artifact-tool";

const projectRoot = "/Users/gzy2520/Desktop/Research/kla";
const tablesDir = path.join(projectRoot, "reanalysis/results/tables");
const comparisonPath = path.join(
  tablesDir,
  "cell_type_kla_vs_reference_ddr_statistics.csv",
);
const tallPath = path.join(tablesDir, "tall104_surrogate_ddr_sensitivity.csv");
const controlInfoPath = path.join(
  tablesDir,
  "cell_type_reference_control_information.csv",
);
const outputPath = path.join(
  tablesDir,
  "cell_type_kla_vs_reference_ddr_statistics.xlsx",
);
const previewDir = "/tmp/kla-reference-ddr-previews";

function parseCsv(text) {
  const rows = [];
  let row = [];
  let field = "";
  let quoted = false;
  const source = text.replace(/^\uFEFF/, "");
  for (let index = 0; index < source.length; index += 1) {
    const char = source[index];
    if (quoted) {
      if (char === '"' && source[index + 1] === '"') {
        field += '"';
        index += 1;
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
  if (field || row.length) {
    row.push(field.replace(/\r$/, ""));
    rows.push(row);
  }
  const headers = rows[0];
  return rows.slice(1).filter((item) => item.some(Boolean)).map((values) =>
    Object.fromEntries(headers.map((header, index) => [header, values[index] ?? ""])),
  );
}

const comparison = parseCsv(await fs.readFile(comparisonPath, "utf8"));
const tall = parseCsv(await fs.readFile(tallPath, "utf8"));
const controlInfo = parseCsv(await fs.readFile(controlInfoPath, "utf8"));
if (comparison.length !== 10) {
  throw new Error(`Expected 10 comparison rows, found ${comparison.length}`);
}
if (controlInfo.length !== 10) {
  throw new Error(`Expected 10 control-information rows, found ${controlInfo.length}`);
}

const chineseInfo = {
  "Human hippocampus": {
    name: "人海马组织",
    matchType: "组织精确匹配",
    subset: "74名神经学正常的CA1海马供者",
    baseline: "神经学正常的死后CA1海马组织；未进行PTM富集",
    sampleCount: "74名供者，74个raw文件",
    acquisition: "非标记LC-MS/MS；Q Exactive HF",
    search: "Spectronaut Pulsar X；1% FDR；至少3条肽段用于定量",
    ptm: "否",
    rationale: "人CA1海马组织精确匹配；74名神经学正常供者；年龄范围接近Kla海马样本；无PTM富集且有论文发布的处理矩阵。",
    caveat: "供者年龄为66-104岁，材料为CA1；适合作为老年海马背景，但不能代表年轻海马或整个海马组织。",
    countDetail: "74名供者共2,092个论文发布的基因/蛋白特征；单供者2,018-2,092个，中位数2,090个。",
    completeness: "ProteomeXchange为PARTIAL，但74个raw文件和论文处理矩阵均可获得。",
    statistics: "可在控制年龄、性别和死后间隔后进行供者层面统计。",
  },
  HEK293T: {
    name: "HEK293T",
    matchType: "细胞系精确匹配",
    subset: "Control_HEK293T_lys；主分析排除Control_HEK293T_std_H002标准QC",
    baseline: "未处理HEK293T裂解物过程控制；未进行PTM富集",
    sampleCount: "401个裂解物运行；另有663个标准QC运行",
    acquisition: "DIA/SWATH；TripleTOF 6600",
    search: "DIA-NN检索，随后使用maxLFQ定量",
    ptm: "否",
    rationale: "细胞系精确匹配；未处理裂解物；标准化DIA平台且运行数多；标准QC材料被明确排除。",
    caveat: "这些样本是过程控制而非独立生物学重复，只能作为蛋白检出背景，不能把401个运行当作401个生物重复。",
    countDetail: "401个Control_HEK293T_lys运行中共检出6,441个高置信BaseAccession；单运行1,987-4,759个，中位数3,512个。",
    completeness: "PRIDE为PARTIAL；Figshare提供最终矩阵、映射和肽段计数文件。",
    statistics: "不可将过程控制运行作为独立生物重复进行差异检验。",
  },
  "HK-2": {
    name: "HK-2",
    matchType: "细胞系精确匹配",
    subset: "HK-2_Control_1、HK-2_Control_3和HK-2_Control_4",
    baseline: "KSFM加2% FBS培养24小时的未处理对照；未进行PTM富集",
    sampleCount: "3个未处理对照；项目共9个raw文件",
    acquisition: "DirectDIA非标记LC-MS/MS；Orbitrap Exploris 240",
    search: "Spectronaut 19 directDIA；跨运行归一化；未插补",
    ptm: "否",
    rationale: "细胞系精确匹配；有3个明确的未处理对照；蛋白组完整度高；提供Spectronaut蛋白组报告。",
    caveat: "只纳入Control_1、Control_3和Control_4，不能混入Cu或TTM处理组；保留仓库原始编号。",
    countDetail: "3个对照的4,791个定量蛋白组展开后得到4,897个唯一BaseAccession；单对照蛋白组数为4,514-4,670。",
    completeness: "ProteomeXchange为PARTIAL，但9个raw文件和关键Spectronaut报告可获得。",
    statistics: "可对3个未处理对照进行组内描述，但不能混入Cu或TTM样本。",
  },
  MCF10A: {
    name: "MCF10A",
    matchType: "细胞系精确匹配",
    subset: "未处理MCF-10A基线细胞裂解物",
    baseline: "汇合状态的未处理MCF-10A细胞；未进行PTM富集",
    sampleCount: "10个肽段IEF组分×2次进样，共20个raw文件",
    acquisition: "DDA；LTQ Orbitrap Velos；240分钟梯度",
    search: "MaxQuant 1.4.1.2 / Andromeda",
    ptm: "否；肽段等电聚焦属于分析分级，不是PTM富集",
    rationale: "未处理MCF-10A精确匹配；分级深度较高；仓库提供MaxQuant evidence，适合作为检出背景。",
    caveat: "平台较旧且进行了深度分级；组分和重复进样属于技术深度，不能作为生物重复，也不宜与现代DIA绝对强度直接比较。",
    countDetail: "20个MCF-10A组分/进样的evidence中共有4,839个有效leading-razor BaseAccession；仓库未提供proteinGroups.txt。",
    completeness: "ProteomeXchange为PARTIAL；包含大量raw文件和一个合并MaxQuant检索结果包。",
    statistics: "不可将IEF组分和重复进样作为独立生物重复。",
  },
  MCF7: {
    name: "MCF7",
    matchType: "细胞系精确匹配",
    subset: "MCF7 / SIDM00148",
    baseline: "未处理基线培养；未进行PTM富集",
    sampleCount: "6次DIA采集",
    acquisition: "DIA/SWATH；TripleTOF 6600",
    search: "DIA-NN检索，随后使用maxLFQ定量",
    ptm: "否",
    rationale: "未处理MCF7精确匹配；来自统一的泛癌DIA图谱；具有高置信蛋白矩阵和明确的细胞系映射。",
    caveat: "Kla MCF7数据包含DCA、缺氧、鱼藤酮和同位素示踪条件；该对照只能作为未处理检出背景，不是配对实验对照。",
    countDetail: "MCF7/SIDM00148平均细胞系行中有3,099个非缺失高置信BaseAccession。",
    completeness: "PRIDE为PARTIAL；Figshare提供完整分析矩阵。",
    statistics: "6次采集主要反映技术采集，不能直接作为独立生物重复。",
  },
  HCT116: {
    name: "HCT116",
    matchType: "细胞系精确匹配",
    subset: "HCT-116 / SIDM00783",
    baseline: "未处理基线培养；未进行PTM富集",
    sampleCount: "6次DIA采集",
    acquisition: "DIA/SWATH；TripleTOF 6600",
    search: "DIA-NN检索，随后使用maxLFQ定量",
    ptm: "否",
    rationale: "未处理HCT-116精确匹配；统一DIA平台；可作为独立于Kla亚细胞分级实验的全细胞检出背景。",
    caveat: "PXD053474含亚细胞组分；全细胞参考不能用于推断某蛋白在特定亚细胞组分中缺失。",
    countDetail: "HCT-116/SIDM00783平均细胞系行中有4,010个非缺失高置信BaseAccession。",
    completeness: "PRIDE为PARTIAL；Figshare提供完整分析矩阵。",
    statistics: "6次采集主要反映技术采集，不能直接作为独立生物重复。",
  },
  "T-ALL": {
    name: "T-ALL（TALL-1替代）",
    matchType: "疾病类别替代，非精确细胞系匹配",
    subset: "主替代：TALL-1 / SIDM00370；敏感性替代：Jurkat / SIDM01016",
    baseline: "未处理T淋巴母细胞白血病细胞系基线；未进行PTM富集",
    sampleCount: "TALL-1为6次DIA采集；Jurkat为6次DIA采集",
    acquisition: "DIA/SWATH；TripleTOF 6600",
    search: "DIA-NN检索，随后使用maxLFQ定量",
    ptm: "否",
    rationale: "未找到独立、未处理的TALL-104全蛋白组；因此TALL-1作为主疾病类别替代，Jurkat单独用于敏感性分析。",
    caveat: "TALL-1和Jurkat均不是TALL-104。PXD028488的TALL-104非富集组经历乳酸暴露，也不是独立正常基线。",
    countDetail: "TALL-1主分析3,383个蛋白；Jurkat 3,363个；并集3,835个；交集2,911个。",
    completeness: "PRIDE为PARTIAL；Figshare提供完整分析矩阵。",
    statistics: "只能作为替代背景，不可写成TALL-104精确对照。",
  },
  "MDA-MB-468": {
    name: "MDA-MB-468",
    matchType: "细胞系精确匹配",
    subset: "MDA-MB-468 / SIDM00628",
    baseline: "未处理基线培养；未进行PTM富集",
    sampleCount: "6次DIA采集",
    acquisition: "DIA/SWATH；TripleTOF 6600",
    search: "DIA-NN检索，随后使用maxLFQ定量",
    ptm: "否",
    rationale: "未处理MDA-MB-468精确匹配；来自统一的泛癌DIA图谱和高置信蛋白矩阵。",
    caveat: "只用于同细胞系的蛋白检出背景，不与MaxQuant Kla富集数据直接比较绝对强度。",
    countDetail: "MDA-MB-468/SIDM00628平均细胞系行中有3,760个非缺失高置信BaseAccession。",
    completeness: "PRIDE为PARTIAL；Figshare提供完整分析矩阵。",
    statistics: "6次采集主要反映技术采集，不能直接作为独立生物重复。",
  },
  "T-47D": {
    name: "T-47D",
    matchType: "细胞系精确匹配",
    subset: "T47D / SIDM00097",
    baseline: "未处理基线培养；未进行PTM富集",
    sampleCount: "6次DIA采集",
    acquisition: "DIA/SWATH；TripleTOF 6600",
    search: "DIA-NN检索，随后使用maxLFQ定量",
    ptm: "否",
    rationale: "未处理T47D精确匹配；来自统一泛癌DIA图谱；已将来源标签T47D明确映射到项目标签T-47D。",
    caveat: "来源中的名称为T47D而非T-47D；只作检出背景。",
    countDetail: "T47D/SIDM00097平均细胞系行中有3,516个非缺失高置信BaseAccession。",
    completeness: "PRIDE为PARTIAL；Figshare提供完整分析矩阵。",
    statistics: "6次采集主要反映技术采集，不能直接作为独立生物重复。",
  },
  RKO: {
    name: "RKO",
    matchType: "细胞系精确匹配",
    subset: "RKO / SIDM01090",
    baseline: "未处理基线培养；未进行PTM富集",
    sampleCount: "6次DIA采集",
    acquisition: "DIA/SWATH；TripleTOF 6600",
    search: "DIA-NN检索，随后使用maxLFQ定量",
    ptm: "否",
    rationale: "未处理RKO精确匹配；来自统一泛癌DIA图谱；作为Kla WT/KO实验的共同检出背景。",
    caveat: "外部RKO蛋白组不能替代GSK3B WT与KO的组内比较。",
    countDetail: "RKO/SIDM01090平均细胞系行中有3,615个非缺失高置信BaseAccession。",
    completeness: "PRIDE为PARTIAL；Figshare提供完整分析矩阵。",
    statistics: "6次采集主要反映技术采集，不能直接作为独立生物重复。",
  },
};

function csvEscape(value) {
  const text = String(value ?? "");
  return /[",\r\n]/.test(text) ? `"${text.replaceAll('"', '""')}"` : text;
}

function toCsv(headers, rows) {
  return [headers, ...rows]
    .map((row) => row.map(csvEscape).join(","))
    .join("\n");
}

const workbook = Workbook.create();

function addComparisonSheet() {
  const sheet = workbook.worksheets.add("DDR占比对照");
  sheet.showGridLines = false;
  sheet.getRange("A1:P1").merge();
  sheet.getRange("A1").values = [["常规蛋白组与 Kla 蛋白组的 DDR 占比对照"]];
  sheet.getRange("A1:P1").format = {
    fill: "#1F4E78",
    font: { bold: true, color: "#FFFFFF", size: 16 },
    horizontalAlignment: "center",
    verticalAlignment: "center",
  };
  sheet.getRange("A1:P1").format.rowHeight = 32;

  sheet.getRange("A2:P2").merge();
  sheet.getRange("A2").values = [[
    "DDR占比 = GO repair/damage 相交蛋白数 / 对应蛋白集合总数。T-ALL 主参考为 TALL-1；海马体参考分母为论文发布的 gene/protein feature。",
  ]];
  sheet.getRange("A2:P2").format = {
    fill: "#D9EAF7",
    font: { color: "#17365D", size: 10 },
    wrapText: true,
    verticalAlignment: "center",
  };
  sheet.getRange("A2:P2").format.rowHeight = 32;

  const headers = [
    "细胞系或组织",
    "参考PXD",
    "参考样本子集",
    "参考蛋白数",
    "参考DDR蛋白数",
    "参考DDR占比",
    "Kla蛋白数",
    "Kla∩DDR蛋白数",
    "Kla DDR占比",
    "Kla-参考百分点差",
    "Kla/参考占比倍数",
    "参考计数单位",
    "登录号匹配数",
    "基因符号辅助匹配数",
    "蛋白数计数口径",
    "描述性结论",
  ];
  sheet.getRange("A4:P4").values = [headers];
  sheet.getRange("A4:P4").format = {
    fill: "#5B9BD5",
    font: { bold: true, color: "#FFFFFF" },
    horizontalAlignment: "center",
    verticalAlignment: "center",
    wrapText: true,
  };
  sheet.getRange("A4:P4").format.rowHeight = 32;

  const startRow = 5;
  const endRow = startRow + comparison.length - 1;
  const values = comparison.map((row) => {
    const zh = chineseInfo[row.CellOrTissueType];
    return [
    zh.name,
    row.ReferencePXD,
    zh.subset,
    Number(row.ReferenceProteinCount),
    Number(row.ReferenceDdrProteinCount),
    null,
    Number(row.TotalKlaProteins),
    Number(row.KlaGoDdrProteins),
    null,
    null,
    null,
    row.ReferenceCountUnit === "BaseAccession"
      ? "UniProt基础登录号"
      : "论文发布的基因/蛋白特征",
    Number(row.BaseAccessionMatches),
    Number(row.GeneSymbolFallbackMatches),
    zh.countDetail,
    null,
  ];
  });
  sheet.getRange(`A${startRow}:P${endRow}`).values = values;
  sheet.getRange(`F${startRow}`).formulas = [[`=IFERROR(E${startRow}/D${startRow},0)`]];
  sheet.getRange(`F${startRow}:F${endRow}`).fillDown();
  sheet.getRange(`I${startRow}`).formulas = [[`=IFERROR(H${startRow}/G${startRow},0)`]];
  sheet.getRange(`I${startRow}:I${endRow}`).fillDown();
  sheet.getRange(`J${startRow}`).formulas = [[`=I${startRow}-F${startRow}`]];
  sheet.getRange(`J${startRow}:J${endRow}`).fillDown();
  sheet.getRange(`K${startRow}`).formulas = [[`=IFERROR(I${startRow}/F${startRow},0)`]];
  sheet.getRange(`K${startRow}:K${endRow}`).fillDown();
  sheet.getRange(`P${startRow}`).formulas = [[
    `=IF(ABS(J${startRow})<0.002,"基本相当",IF(J${startRow}>0,"Kla中DDR占比更高","参考组中DDR占比更高"))`,
  ]];
  sheet.getRange(`P${startRow}:P${endRow}`).fillDown();

  sheet.getRange(`A${startRow}:P${endRow}`).format = {
    wrapText: true,
    verticalAlignment: "top",
    borders: {
      insideHorizontal: { style: "thin", color: "#D9E2F3" },
      bottom: { style: "thin", color: "#A6A6A6" },
    },
  };
  sheet.getRange(`D${startRow}:E${endRow}`).format.numberFormat = "#,##0";
  sheet.getRange(`G${startRow}:H${endRow}`).format.numberFormat = "#,##0";
  sheet.getRange(`M${startRow}:N${endRow}`).format.numberFormat = "#,##0";
  sheet.getRange(`F${startRow}:F${endRow}`).format.numberFormat = "0.00%";
  sheet.getRange(`I${startRow}:J${endRow}`).format.numberFormat = "0.00%";
  sheet.getRange(`K${startRow}:K${endRow}`).format.numberFormat = "0.00x";
  sheet.getRange(`A${startRow}:P${endRow}`).format.rowHeight = 58;
  [
    22, 15, 40, 14, 17, 16, 14, 18, 16, 18, 18, 28, 16, 16, 58, 24,
  ].forEach((width, index) => {
    sheet.getRangeByIndexes(0, index, endRow, 1).format.columnWidth = width;
  });
  sheet.getRange(`J${startRow}:J${endRow}`).conditionalFormats.add("cellIs", {
    operator: "greaterThan",
    formula: 0,
    format: { fill: "#FCE4D6", font: { color: "#9C0006" } },
  });
  sheet.getRange(`J${startRow}:J${endRow}`).conditionalFormats.add("cellIs", {
    operator: "lessThan",
    formula: 0,
    format: { fill: "#E2F0D9", font: { color: "#006100" } },
  });
  sheet.tables.add(`A4:P${endRow}`, true, "ReferenceVsKlaDdrTable");
  sheet.freezePanes.freezeRows(4);
  sheet.freezePanes.freezeColumns(1);
}

function addTallSheet() {
  const sheet = workbook.worksheets.add("TALL替代敏感性");
  sheet.showGridLines = false;
  const headers = ["参考集合", "蛋白数", "DDR蛋白数", "DDR占比"];
  sheet.getRange("A1:D1").values = [headers];
  sheet.getRange("A1:D1").format = {
    fill: "#1F4E78",
    font: { bold: true, color: "#FFFFFF" },
    horizontalAlignment: "center",
  };
  const values = tall.map((row) => [
    row.TAllReferenceSet,
    Number(row.ProteinCount),
    Number(row.DdrProteinCount),
    null,
  ]);
  sheet.getRange(`A2:D${values.length + 1}`).values = values;
  sheet.getRange("D2").formulas = [["=IFERROR(C2/B2,0)"]];
  sheet.getRange(`D2:D${values.length + 1}`).fillDown();
  sheet.getRange(`B2:C${values.length + 1}`).format.numberFormat = "#,##0";
  sheet.getRange(`D2:D${values.length + 1}`).format.numberFormat = "0.00%";
  sheet.getRange(`A1:D${values.length + 1}`).format.borders = {
    insideHorizontal: { style: "thin", color: "#D9E2F3" },
    bottom: { style: "thin", color: "#A6A6A6" },
  };
  [34, 16, 18, 16].forEach((width, index) => {
    sheet.getRangeByIndexes(0, index, values.length + 1, 1).format.columnWidth = width;
  });
  sheet.tables.add(`A1:D${values.length + 1}`, true, "TallSensitivityTable");
  sheet.freezePanes.freezeRows(1);
}

function addControlInformationSheet() {
  const sheet = workbook.worksheets.add("对照选择信息");
  sheet.showGridLines = false;
  const headers = [
    "细胞或组织",
    "参考PXD",
    "年份",
    "匹配类型",
    "适用等级",
    "选用样本或子集",
    "正常或基线条件",
    "样本或采集数",
    "对照蛋白数",
    "DDR蛋白数",
    "DDR占比",
    "采集方式或仪器",
    "检索与定量",
    "是否PTM富集",
    "选择理由",
    "主要限制",
    "蛋白数计数口径",
    "实际分析来源文件",
    "仓库与分析完整度",
    "可作统计差异",
    "数据集URL",
    "论文URL",
    "处理数据URL",
  ];
  sheet.getRange("A1:W1").values = [headers];
  sheet.getRange("A1:W1").format = {
    fill: "#1F4E78",
    font: { bold: true, color: "#FFFFFF" },
    horizontalAlignment: "center",
    verticalAlignment: "center",
    wrapText: true,
  };
  sheet.getRange("A1:W1").format.rowHeight = 34;
  const values = controlInfo.map((row) => {
    const zh = chineseInfo[row["细胞或组织"]];
    return [
    zh.name,
    row["参考PXD"],
    Number(row["年份"]),
    zh.matchType,
    row["适用等级"],
    zh.subset,
    zh.baseline,
    zh.sampleCount,
    Number(row["对照蛋白数"]),
    Number(row["DDR蛋白数"]),
    Number(row["DDR占比"]),
    zh.acquisition,
    zh.search,
    zh.ptm,
    zh.rationale,
    zh.caveat,
    zh.countDetail,
    row["实际分析来源文件"],
    zh.completeness,
    zh.statistics,
    row["数据集URL"],
    row["论文URL"],
    row["处理数据URL"],
  ];
  });
  sheet.getRange(`A2:W${values.length + 1}`).values = values;
  sheet.getRange(`A2:W${values.length + 1}`).format = {
    wrapText: true,
    verticalAlignment: "top",
    borders: {
      insideHorizontal: { style: "thin", color: "#D9E2F3" },
      bottom: { style: "thin", color: "#A6A6A6" },
    },
  };
  sheet.getRange(`C2:C${values.length + 1}`).format.numberFormat = "0";
  sheet.getRange(`I2:J${values.length + 1}`).format.numberFormat = "#,##0";
  sheet.getRange(`K2:K${values.length + 1}`).format.numberFormat = "0.00%";
  sheet.getRange(`A2:W${values.length + 1}`).format.rowHeight = 76;
  [
    22, 15, 10, 28, 12, 46, 44, 34, 16, 16, 14, 34, 40, 18, 58, 58,
    58, 54, 50, 36, 46, 42, 46,
  ].forEach((width, index) => {
    sheet.getRangeByIndexes(0, index, values.length + 1, 1).format.columnWidth = width;
  });
  controlInfo.forEach((row, index) => {
    const gradeFill =
      row["适用等级"] === "A"
        ? "#D9EAD3"
        : row["适用等级"] === "B"
          ? "#FFF2CC"
          : "#F4CCCC";
    sheet.getCell(index + 1, 4).format.fill = gradeFill;
    sheet.getCell(index + 1, 4).format.font = { bold: true };
    sheet.getCell(index + 1, 4).format.horizontalAlignment = "center";
  });
  sheet.tables.add(`A1:W${values.length + 1}`, true, "ReferenceControlInformationTable");
  sheet.freezePanes.freezeRows(1);
  sheet.freezePanes.freezeColumns(2);
}

function addMethodsSheet() {
  const sheet = workbook.worksheets.add("方法与解释");
  sheet.showGridLines = false;
  const rows = [
    ["项目", "说明"],
    ["DDR定义", "使用 data/annotations/GO-repair+damage(human).tsv；仅保留人源且排除 qualifier 含 NOT 的注释。"],
    ["匹配顺序", "UniProt BaseAccession 优先，来源数据具有 GeneSymbol 时再作辅助匹配。"],
    ["主分母", "每个参考对照实际纳入的唯一 BaseAccession；海马体为论文发布的 2,092 个 gene/protein feature。"],
    ["PXD030304", "使用 6,692 蛋白高置信矩阵。HEK293T 由 401 个 Control_HEK293T_lys 运行的 peptide-count 并集计算，标准QC排除。"],
    ["T-ALL", "10类主表使用 TALL-1 作为 TALL-104 的疾病类别替代；Jurkat、并集和交集结果单独保留。"],
    ["统计解释", "本结果是描述性蛋白集合比例，不把技术进样当作生物重复，不进行跨平台丰度显著性检验。"],
    ["未检出解释", "常规蛋白组或 Kla 数据未检出均不能直接解释为蛋白不存在或不发生乳酸化。"],
  ];
  sheet.getRange(`A1:B${rows.length}`).values = rows;
  sheet.getRange("A1:B1").format = {
    fill: "#70AD47",
    font: { bold: true, color: "#FFFFFF" },
  };
  sheet.getRange(`A1:B${rows.length}`).format.wrapText = true;
  sheet.getRange(`A1:B${rows.length}`).format.verticalAlignment = "top";
  sheet.getRange(`A2:B${rows.length}`).format.rowHeight = 50;
  sheet.getRange("A:A").format.columnWidth = 22;
  sheet.getRange("B:B").format.columnWidth = 86;
  sheet.getRange(`A1:B${rows.length}`).format.borders = {
    insideHorizontal: { style: "thin", color: "#D9E2F3" },
    bottom: { style: "thin", color: "#A6A6A6" },
  };
}

addComparisonSheet();
addControlInformationSheet();
addTallSheet();
addMethodsSheet();

const comparisonZhHeaders = [
  "细胞或组织",
  "参考PXD",
  "参考样本子集",
  "参考蛋白数",
  "参考DDR蛋白数",
  "参考DDR占比",
  "Kla蛋白数",
  "Kla与DDR交集蛋白数",
  "Kla DDR占比",
  "Kla减参考百分点差",
  "Kla与参考占比倍数",
  "选择理由",
  "主要限制",
];
const comparisonZhRows = comparison.map((row) => {
  const zh = chineseInfo[row.CellOrTissueType];
  return [
    zh.name,
    row.ReferencePXD,
    zh.subset,
    Number(row.ReferenceProteinCount),
    Number(row.ReferenceDdrProteinCount),
    Number(row.ReferenceDdrFraction),
    Number(row.TotalKlaProteins),
    Number(row.KlaGoDdrProteins),
    Number(row.KlaGoDdrFraction),
    Number(row.DdrFractionPercentagePointDifference),
    Number(row.KlaToReferenceDdrFractionRatio),
    zh.rationale,
    zh.caveat,
  ];
});
const controlZhHeaders = [
  "细胞或组织",
  "参考PXD",
  "年份",
  "匹配类型",
  "适用等级",
  "选用样本或子集",
  "正常或基线条件",
  "样本或采集数",
  "对照蛋白数",
  "DDR蛋白数",
  "DDR占比",
  "采集方式或仪器",
  "检索与定量",
  "是否PTM富集",
  "选择理由",
  "主要限制",
  "蛋白数计数口径",
  "实际分析来源文件",
  "仓库与分析完整度",
  "可作统计差异",
  "数据集URL",
  "论文URL",
  "处理数据URL",
];
const controlZhRows = controlInfo.map((row) => {
  const zh = chineseInfo[row["细胞或组织"]];
  return [
    zh.name,
    row["参考PXD"],
    Number(row["年份"]),
    zh.matchType,
    row["适用等级"],
    zh.subset,
    zh.baseline,
    zh.sampleCount,
    Number(row["对照蛋白数"]),
    Number(row["DDR蛋白数"]),
    Number(row["DDR占比"]),
    zh.acquisition,
    zh.search,
    zh.ptm,
    zh.rationale,
    zh.caveat,
    zh.countDetail,
    row["实际分析来源文件"],
    zh.completeness,
    zh.statistics,
    row["数据集URL"],
    row["论文URL"],
    row["处理数据URL"],
  ];
});
await fs.writeFile(
  path.join(tablesDir, "cell_type_kla_vs_reference_ddr_statistics_zh.csv"),
  `\uFEFF${toCsv(comparisonZhHeaders, comparisonZhRows)}`,
  "utf8",
);
await fs.writeFile(
  path.join(tablesDir, "cell_type_reference_control_information_zh.csv"),
  `\uFEFF${toCsv(controlZhHeaders, controlZhRows)}`,
  "utf8",
);

const errors = await workbook.inspect({
  kind: "match",
  searchTerm: "#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A",
  options: { useRegex: true, maxResults: 100 },
  summary: "final formula error scan",
});

const checks = [];
for (const [sheetName, range] of [
  ["DDR占比对照", "A1:P14"],
  ["对照选择信息", "A1:K11"],
  ["对照选择信息", "L1:W11"],
  ["TALL替代敏感性", "A1:D5"],
  ["方法与解释", "A1:B8"],
]) {
  const inspected = await workbook.inspect({
    kind: "table",
    range: `${sheetName}!${range}`,
    include: "values,formulas",
    tableMaxRows: 20,
    tableMaxCols: 18,
    maxChars: 12000,
  });
  checks.push({ sheetName, range, ndjson: inspected.ndjson });
}

await fs.mkdir(previewDir, { recursive: true });
for (const [fileName, sheetName, range] of [
  ["comparison.png", "DDR占比对照", "A1:P14"],
  ["control_info_left.png", "对照选择信息", "A1:K11"],
  ["control_info_right.png", "对照选择信息", "L1:W11"],
  ["tall.png", "TALL替代敏感性", "A1:D5"],
  ["methods.png", "方法与解释", "A1:B8"],
]) {
  const preview = await workbook.render({ sheetName, range, scale: 1, format: "png" });
  await fs.writeFile(
    path.join(previewDir, fileName),
    new Uint8Array(await preview.arrayBuffer()),
  );
}

const output = await SpreadsheetFile.exportXlsx(workbook);
await output.save(outputPath);
await fs.writeFile(
  path.join(previewDir, "verification.json"),
  JSON.stringify({ formulaErrors: errors.ndjson, checks }, null, 2),
);
console.log(JSON.stringify({ outputPath, rows: comparison.length, formulaErrors: errors.ndjson }));
