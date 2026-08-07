const fs = require("fs");
const path = require("path");
const {
  AlignmentType,
  Document,
  Footer,
  ImageRun,
  LineRuleType,
  PageBreak,
  PageNumber,
  Packer,
  Paragraph,
  TextRun,
} = require("docx");

const projectRoot = path.resolve(__dirname, "../..");
const outputPath = path.join(
  projectRoot,
  "reanalysis/results/Kla_Venn_figure_legends_bilingual.docx",
);

const englishFont = {
  ascii: "Times New Roman",
  hAnsi: "Times New Roman",
  eastAsia: "SimSun",
};
const chineseFont = {
  ascii: "Times New Roman",
  hAnsi: "Times New Roman",
  eastAsia: "SimSun",
};
const headingFont = {
  ascii: "Arial",
  hAnsi: "Arial",
  eastAsia: "SimHei",
};

function title(text) {
  return new Paragraph({
    alignment: AlignmentType.CENTER,
    spacing: { after: 120, line: 288, lineRule: LineRuleType.AUTO },
    children: [
      new TextRun({
        text,
        bold: true,
        size: 32,
        font: headingFont,
      }),
    ],
  });
}

function note(text) {
  return new Paragraph({
    alignment: AlignmentType.CENTER,
    spacing: { after: 180, line: 240, lineRule: LineRuleType.AUTO },
    children: [
      new TextRun({
        text,
        italics: true,
        color: "666666",
        size: 18,
        font: chineseFont,
      }),
    ],
  });
}

function sectionLabel(text, pageBreakBefore = false) {
  return new Paragraph({
    pageBreakBefore,
    spacing: { before: 60, after: 80, line: 264, lineRule: LineRuleType.AUTO },
    children: [
      new TextRun({
        text,
        bold: true,
        size: 24,
        font: headingFont,
      }),
    ],
  });
}

function figure(imagePath, width, height, description) {
  return new Paragraph({
    alignment: AlignmentType.CENTER,
    spacing: { after: 100 },
    children: [
      new ImageRun({
        type: "png",
        data: fs.readFileSync(imagePath),
        transformation: { width, height },
        altText: {
          title: description,
          description,
          name: path.basename(imagePath),
        },
      }),
    ],
  });
}

function legend(label, body, language) {
  const font = language === "zh" ? chineseFont : englishFont;
  return new Paragraph({
    alignment: AlignmentType.JUSTIFIED,
    spacing: {
      before: language === "zh" ? 70 : 30,
      after: language === "zh" ? 100 : 50,
      line: 276,
      lineRule: LineRuleType.AUTO,
    },
    widowControl: true,
    children: [
      new TextRun({
        text: label,
        bold: true,
        size: 21,
        font,
      }),
      new TextRun({
        text: body,
        size: 21,
        font,
      }),
    ],
  });
}

const figure1EnglishTitle =
  "Figure X. Distribution of lysine-lactylated proteins among hippocampal tissue, immortalized models, and tumor cell lines. ";
const figure1EnglishBody =
  "Venn diagram showing the overlap of 3,112 lysine-lactylated (Kla) proteins identified across human hippocampal tissue, immortalized/non-tumor cell models (HEK293T, HK-2, and MCF10A), and tumor cell lines (MCF7, HCT116, T-ALL, MDA-MB-468, T-47D, and RKO). Proteins were collapsed to unique UniProt base accessions after removal of isoform suffixes. Numbers indicate the counts of unique proteins in each mutually exclusive region.";
const figure1ChineseTitle =
  "图 X. 人海马组织、永生化模型与肿瘤细胞系中赖氨酸乳酰化蛋白的分布。";
const figure1ChineseBody =
  "维恩图展示在人海马组织、永生化/非肿瘤细胞模型（HEK293T、HK-2 和 MCF10A）及肿瘤细胞系（MCF7、HCT116、T-ALL、MDA-MB-468、T-47D 和 RKO）中鉴定到的 3,112 个赖氨酸乳酰化（Kla）蛋白的重叠关系。去除异构体后缀后，以唯一 UniProt 基础登录号（BaseAccession）作为蛋白计数单位。数字表示各互斥区域内的唯一蛋白数。";

const figure2EnglishTitle =
  "Figure Y. Distribution of Kla proteins associated with DNA repair and DNA damage responses among hippocampal tissue, immortalized models, and tumor cell lines. ";
const figure2EnglishBody =
  "Venn diagram showing the overlap of 275 Kla proteins annotated to DNA repair- or DNA damage response-related Gene Ontology (GO) biological processes across human hippocampal tissue, immortalized/non-tumor cell models (HEK293T, HK-2, and MCF10A), and tumor cell lines (MCF7, HCT116, T-ALL, MDA-MB-468, T-47D, and RKO). Kla proteins were matched to human GO annotations primarily by UniProt base accession, with gene-symbol matching used only when accession-level matching was unavailable; annotations carrying the qualifier “NOT” were excluded. Numbers indicate the counts of unique proteins in each mutually exclusive region.";
const figure2ChineseTitle =
  "图 Y. 人海马组织、永生化模型与肿瘤细胞系中 DNA 修复和 DNA 损伤应答相关 Kla 蛋白的分布。";
const figure2ChineseBody =
  "维恩图展示在人海马组织、永生化/非肿瘤细胞模型（HEK293T、HK-2 和 MCF10A）及肿瘤细胞系（MCF7、HCT116、T-ALL、MDA-MB-468、T-47D 和 RKO）中鉴定到的 275 个 Kla 蛋白的重叠关系；这些蛋白被注释至 DNA 修复或 DNA 损伤应答相关的基因本体（GO）生物学过程。Kla 蛋白优先通过 UniProt 基础登录号与人源 GO 注释匹配，仅在登录号无法匹配时使用基因符号进行辅助匹配；带有“NOT”限定词的注释被排除。数字表示各互斥区域内的唯一蛋白数。";

const document = new Document({
  creator: "Codex",
  title: "Bilingual figure legends for Kla Venn diagrams",
  description:
    "Publication-ready English and Chinese legends for two Kla Venn diagrams.",
  styles: {
    default: {
      document: {
        run: {
          font: englishFont,
          size: 21,
        },
        paragraph: {
          spacing: { line: 276, lineRule: LineRuleType.AUTO },
        },
      },
    },
  },
  sections: [
    {
      properties: {
        page: {
          size: { width: 11906, height: 16838 },
          margin: {
            top: 1134,
            right: 1418,
            bottom: 1134,
            left: 1418,
          },
        },
      },
      footers: {
        default: new Footer({
          children: [
            new Paragraph({
              alignment: AlignmentType.CENTER,
              children: [
                new TextRun({
                  children: [PageNumber.CURRENT],
                  size: 18,
                  color: "777777",
                  font: englishFont,
                }),
              ],
            }),
          ],
        }),
      },
      children: [
        title("Bilingual Figure Legends / 双语图注"),
        note("注：Figure X/Y 与图 X/Y 为占位符，使用时请替换为稿件中的实际图号。"),
        sectionLabel("Figure X / 图 X"),
        figure(
          path.join(
            projectRoot,
            "reanalysis/results/figures/all_kla_three_group_venn.png",
          ),
          500,
          449,
          "Venn diagram of all Kla proteins",
        ),
        legend(figure1EnglishTitle, figure1EnglishBody, "en"),
        legend(figure1ChineseTitle, figure1ChineseBody, "zh"),
        sectionLabel("Figure Y / 图 Y", true),
        figure(
          path.join(
            projectRoot,
            "reanalysis/results/figures/kla_go_ddr_three_group_venn.png",
          ),
          478,
          466,
          "Venn diagram of Kla proteins associated with DNA repair and damage responses",
        ),
        legend(figure2EnglishTitle, figure2EnglishBody, "en"),
        legend(figure2ChineseTitle, figure2ChineseBody, "zh"),
      ],
    },
  ],
});

Packer.toBuffer(document).then((buffer) => {
  fs.writeFileSync(outputPath, buffer);
  process.stdout.write(`${outputPath}\n`);
});
