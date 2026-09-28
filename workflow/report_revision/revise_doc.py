#!/usr/bin/env python3
"""Revise the supplied report while preserving its layout and original file."""
from pathlib import Path
from copy import deepcopy
import csv, hashlib, json
from docx import Document
from docx.oxml.ns import qn
ROOT=Path(__file__).resolve().parents[2]
SOURCE=Path('/Users/gzy2520/Desktop/ML+DDR&repair.docx')
DEST=Path('/Users/gzy2520/Desktop/ML+DDR&repair_科学修订版.docx')
OUT=ROOT/'outputs/20260928_report_revision'
doc=Document(SOURCE)
assert len(doc.paragraphs)==74, 'Source document changed; review paragraph anchors.'
changes={
2:'另引入 2026 年新公开的 ESCC 普通蛋白组与乳酸化组数据（PXD063945 / PMC13195773），并使用独立 YSH 队列的 23 位患者、共 46 份配对肿瘤及癌旁 RNA 样本（GSE130078），开展外部材料级一致性评估。RNA 与质谱来自不同队列，不构成患者级配对验证。',
4:'表 1 DDR 与 DNA repair 模型在内部留出、来源关联块留出及外部材料级评估中的性能',
6:'内部留出的 374 条 RNA 记录覆盖 19 种材料，DDR 与 DNA repair 的逐记录 MAE 分别为 1.27 和 1.47 个百分点；材料等权 MAE 分别为 2.46 和 2.79。整个数据集包含 28 种材料，同一材料的 RNA 记录共享蛋白组参考标签。整来源关联块留出的逐记录 MAE 分别为 15.41 和 17.81；外部材料级比较的 MAE 为 12.30～18.53 个百分点。结果支持已知材料体系内的标签重现，尚不支持跨来源或个体乳酸化水平的准确预测。',
14:'▪ DNA repair 集合：GO:0006281（DNA repair）及沿 is_a / part_of 关系纳入的后代，在冻结注释中包含 2,135 个人类 UniProt BaseAccession；当前集合属于 DDR 集合的子集，各材料中的检出分母为 54～325 个蛋白。',
15:'▪ 分母非零：所有材料组的分母均大于零，无需因零分母剔除材料。普通蛋白组未检出的 Kla 蛋白记入 KlaOutsideReference，不纳入分子，确保分子集合包含于分母集合，比例范围为 0%～100%。该指标是不同蛋白种类的检出占比，不代表修饰强度或位点占有率。',
17:'模型输入 18,332 个各来源共有的 Ensembl 基因，训练与外部输入统一为 log2(TPM + 0.5)。DDR 和 DNA repair 限定预测目标 Y，不限制输入特征 X；模型探索全转录组表达与材料级检出占比的关联。随机森林采用 500 棵回归树，mtry = 135，最小节点数 1，随机种子 25，并沿用材料及供者平衡的训练抽样权重。31 个蛋白组材料组按 ReferenceKey 对应为 28 个 RNA 材料标签，同一 ReferenceKey 的多个组占比取算术平均。增加 RNA 记录不等于增加独立蛋白组标签。',
19:'外部数据以组织类型和疾病状态进行材料级匹配，并核查干预与治疗信息。质谱明确排除新辅助化疗组；RNA 为临床组织测序，未见针对这些组织 RNA 的基因敲低等实验干预，但其临床治疗前史尚未确认，因此不能宣称两队列均已证实未经治疗。',
20:'▪ 外部质谱参考：PXD063945 排除新辅助化疗组后，普通蛋白组各有 25 列肿瘤和癌旁数据，Kla 组各有 22 列；不将其视作 25 对完整的两组学配对样本。按材料合并检出列表，肿瘤检出 8,186 个普通蛋白及 688 个 Kla 蛋白，癌旁分别为 7,649 和 472 个。肿瘤 DDR 占比为 47/399 = 11.78%，DNA repair 为 33/233 = 14.16%；癌旁分别为 41/366 = 11.20% 和 27/209 = 12.92%。论文全队列的 8,281 个蛋白与 1,836 个 Kla 位点不作为本筛选子集的统计，位点数也不用于蛋白种类占比的分子。',
21:'▪ 独立 RNA 队列：GSE130078 为 YSH 队列 23 位患者的 23 份 ESCC 肿瘤和 23 份配对癌旁组织，共 46 份 bulk RNA-seq（去 rRNA total RNA-seq）。原始计数按 GENCODE v19 合并外显子长度换算 TPM，再作 log2(TPM + 0.5) 转换，提取模型所需的全部 18,332 个基因。样本标识与训练集无重叠，但肿瘤与癌旁来自同一患者，不能视为 46 位独立受试者。',
24:'图 1 依次展示训练子集回代（n=1,524）、内部留出验证（n=374）、全数据回代拟合（n=1,898）及 14 折来源关联块留出。内部划分沿用供者信息进行隔离；该划分已用于模型开发比较，不作为另行封存的最终盲测集。14 折按来源关联块划分，不能直接解释为 14 个独立实验室。图中指标采用材料等权口径，与表 1 的逐记录口径分开报告。',
27:'3.2 材料级参考占比与逐 RNA 记录预测分布',
28:'图 2 展示 28 种材料的预测分布。金色菱形是材料级共享参考占比，蓝点为训练 RNA 记录（n=1,524），橙点为内部留出记录（n=374，覆盖其中 19 种材料）。这些参考值不是各 RNA 受试者分别测得的蛋白组占比。',
30:'图 2 Kla-DDR / DDR 模型逐 RNA 记录预测与共享材料标签的比较',
31:'多数内部留出记录接近对应材料参考值，但不同材料与类别的误差不一致，仍存在偏差较大的记录。未包含内部留出记录的材料不能据此认定其验证性能。该图用于描述已知材料体系内的标签重现程度。',
32:'3.3 Top 25 预测特征',
33:'图 3 按随机森林回归的方差减少重要性，对全部 18,332 个输入基因排序并展示前 25 个特征。分析使用 Ensembl ID，Gene Symbol 仅作为辅助显示。该排名反映模型内的预测贡献，可能受到相关特征和材料差异影响，不能解释为乳酸化调控因果效应，也不是 Gini 分类重要性。',
37:'图 4 使用 1,898 条 RNA 记录的前 2,000 个高变基因、前 10 个主成分构建 UMAP，依次按生物学类别、共享材料参考占比及全数据训练模型的回代拟合值着色。',
39:'图 4 Kla-DDR / DDR 的转录组 UMAP：类别、材料参考值与全数据拟合值',
40:'参考标签与拟合值的颜色分布相近，与训练内拟合较好一致。由于 C 图使用回代拟合值，该图仅用于描述转录组结构与标签分布，不构成独立预测能力或个体修饰连续梯度的验证证据。',
42:'图 5 展示训练与内部留出记录的残差、LOESS 趋势和误差密度；图 6 按四种生物学类别分别展示训练与内部留出结果。',
46:'图 6 Kla-DDR / DDR 按正常／癌组织及正常／癌细胞系分层的训练与内部留出结果',
47:'内部留出中，正常组织 MAE 为 1.20 个百分点、误差在 ±5 个百分点内的比例为 98.2%；癌组织分别为 1.10 和 97.8%。残差在零附近集中，但仍需结合趋势与尾部偏差判断，不能据此声称无偏或临床预测稳健。正常细胞系与癌细胞系的留出记录分别仅 9 和 8 条，类别间表现不能一概而论。',
53:'图 8 Kla-DNA repair / DNA repair 的材料级参考标签与逐 RNA 记录预测分布',
57:'图 10 Kla-DNA repair / DNA repair 的 UMAP：材料标签与全数据回代拟合，仅作描述性展示',
59:'图 11 Kla-DNA repair / DNA repair 的训练与内部留出残差诊断',
61:'图 12 Kla-DNA repair / DNA repair 按生物学类别分层的训练与内部留出结果',
62:'五、外部队列的材料级一致性评估',
63:'保持 DDR 与 DNA repair 模型冻结，输入 GSE130078 的 46 份组织 RNA，与 PXD063945 对应的肿瘤及癌旁材料参考占比比较。每个预测目标只有两个材料级参考标签；MAE 衡量逐 RNA 预测值对这些参考标签的偏离，并非针对 46 个个体实测占比计算的误差。外部数据未用于重新训练或调参。本次修订将外部表达转换统一为训练端的 log2(TPM + 0.5)，以下表格和图均采用更正后的结果。',
65:'图 13 冻结模型的外部材料级比较：独立 RNA 队列预测值与另一质谱队列的材料参考占比',
66:'表 2 更正输入转换后的外部材料级评估结果',
68:'六、外部差异的解释与适用范围',
69:'外部 RNA 的平均预测值约为 24%～31%，对应质谱材料参考占比约为 11%～14%；四项比较的 MAE 为 12.30～18.53 个百分点，误差在 ±5 个百分点内的比例均为 0%。相对于此次外部材料参考值，预测整体偏高。下列因素是可能解释，当前分析尚不能区分其贡献。',
70:'▪ 材料标签与来源差异：既有 ESCC 材料的训练参考标签为 DDR 21.19%、DNA repair 26.12%，高于新质谱参考值。模型输出可能受到训练材料与来源特征影响，但回归结果不能证明模型完成了 ESCC 分类，也不能据此断言其只学习了材料身份。',
71:'▪ 质谱检测差异：普通蛋白组覆盖范围、Kla 富集与检出效率可能影响分子和分母。具体抗体、仪器或色谱差异的作用未在本分析中分离验证，不能据此作因果归因；尤其不能将 Kla 位点数与蛋白种类数直接相比来解释通路蛋白占比。',
72:'▪ 队列与配对限制：RNA 与质谱来自不同人群，仅在组织类型和疾病类别上匹配；RNA 治疗前史未确认。该比较用于考察外部材料标签的一致性，不能估计患者级预测误差。综上，模型可在已知材料体系内较好地重现共享参考占比，但跨来源表现不足，尚不支持个体乳酸化水平预测或调控机制推断。'
}
def put(p,text):
    if p.runs:
        p.runs[0].text=text
        for r in p.runs[1:]:r.text=''
    else:p.add_run(text)
for i,t in changes.items():put(doc.paragraphs[i],t)
# Table text and figures are updated together; no model refitting.
for tab in doc.tables:
 for row in tab.rows:
  for cell in row.cells:
   for p in cell.paragraphs:
    t=p.text.replace('内部盲测','内部留出').replace('外部新盲测','外部材料比较').replace('整来源留出','来源块留出').replace('外部质谱真值','质谱材料参考值').replace('盲测 MAE','材料比较 MAE').replace('样本量 (n)','RNA 记录数').replace('ESCC 未治疗肿瘤切片','ESCC 肿瘤组织').replace('ESCC 配对癌旁切片','ESCC 配对癌旁组织')
    if t!=p.text:put(p,t)
rows=list(csv.DictReader((OUT/'external/summary.csv').open()))
metrics={(r['Target'],r['Material']):r for r in rows}
pairs=[('DDR','ESCC_untreated'),('DDR','ESCC_adjacent_untreated'),('DNA_repair','ESCC_untreated'),('DNA_repair','ESCC_adjacent_untreated')]
def cellset(c,t):
 put(c.paragraphs[0],str(t))
 for p in c.paragraphs[1:]:put(p,'')
for i,key in enumerate(pairs,1):
 r=metrics[key];cells=doc.tables[1].rows[i].cells
 for j,value in {3:f"{float(r['Label']):.2f}%",4:f"{float(r['MeanPrediction']):.2f}%",5:f"{float(r['MAE']):.2f}",6:f"{float(r['Within5']):.1f}%"}.items():cellset(cells[j],value)
t=doc.tables[0]
new=deepcopy(t.rows[13]._tr);t._tbl.append(new)
for i,key in zip([6,7,13,14],pairs):
 r=metrics[key];cellset(t.rows[i].cells[4],f"{float(r['MAE']):.2f}")
cellset(t.rows[14].cells[1],'外部材料比较 (GSE130078 癌旁 vs PXD063945)')
# Clarify paired RNA and shared-label weighting below tables using existing paragraph 5.
put(doc.paragraphs[5],'注：pp 为百分点；±5 pp 命中率指误差在该范围内的记录比例。材料等权指标先赋予各材料相同总权重，再计算记录级误差。外部肿瘤与癌旁来自同一批 23 位患者，各组只有一个共享参考标签，故不报告组内 R²。内部验证不等于独立封存的最终测试集。')
imgs={1:OUT/'DDR/figures/figure1_material_distributions.png',2:OUT/'DDR/figures/top25.png',3:OUT/'DDR/figures/figure3_manifold_projection.png',4:OUT/'DDR/figures/figure4_residual_diagnostics.png',5:OUT/'DDR/figures/figure5_category_stratified.png',7:OUT/'DNA_repair/figures/figure1_material_distributions.png',8:OUT/'DNA_repair/figures/top25.png',9:OUT/'DNA_repair/figures/figure3_manifold_projection.png',10:OUT/'DNA_repair/figures/figure4_residual_diagnostics.png',11:OUT/'DNA_repair/figures/figure5_category_stratified.png',12:OUT/'external/figures/predictions.png'}
for i,path in imgs.items():
 rid=doc.inline_shapes[i]._inline.graphic.graphicData.pic.blipFill.blip.embed
 doc.part.related_parts[rid]._blob=path.read_bytes()
# Keep table headers visible after page breaks; do not split a single table row.
from docx.oxml import OxmlElement
for table in doc.tables:
 for row in table.rows:
  prop=row._tr.get_or_add_trPr(); prop.append(OxmlElement('w:cantSplit'))
 table.rows[0]._tr.get_or_add_trPr().append(OxmlElement('w:tblHeader'))
# Normalize inherited Word grid spacing so Chinese text remains compact and readable.
from docx.shared import Pt
from docx.enum.text import WD_ALIGN_PARAGRAPH
for p in doc.paragraphs:
 f=p.paragraph_format
 f.line_spacing=1.18; f.space_before=Pt(0); f.space_after=Pt(6)
 grid=OxmlElement('w:snapToGrid');grid.set(qn('w:val'),'0');p._p.get_or_add_pPr().append(grid)
 if p.text:
  isheading=p.style.name.startswith('Heading')
  size=16 if p.style.name=='Heading 1' else 13 if isheading else 11
  if isheading:f.space_before=Pt(10);f.keep_with_next=True
  for r in p.runs:r.font.size=Pt(size)
for table in doc.tables:
 for row in table.rows:
  for cell in row.cells:
   for p in cell.paragraphs:
    p.alignment=WD_ALIGN_PARAGRAPH.LEFT
    p.paragraph_format.line_spacing=1.05
    p.paragraph_format.space_before=Pt(2);p.paragraph_format.space_after=Pt(2)
    grid=OxmlElement('w:snapToGrid');grid.set(qn('w:val'),'0');p._p.get_or_add_pPr().append(grid)
    for r in p.runs:r.font.size=Pt(9)
# Use locally available Chinese font families for reliable PDF rendering.
for element in [doc.element, doc.styles.element] + [s.header.part.element for s in doc.sections] + [s.footer.part.element for s in doc.sections]:
 for fonts in element.iter(qn('w:rFonts')):
  for attr in ('ascii','hAnsi','eastAsia','cs'):
   key=qn('w:'+attr); val=fonts.get(key)
   if val in ('SimSun','宋体'):fonts.set(key,'Songti SC')
   elif val in ('SimHei','黑体'):fonts.set(key,'Heiti SC')
doc.save(DEST)
(OUT/'revision_manifest.json').write_text(json.dumps({'source':str(SOURCE),'source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),'output':str(DEST),'output_sha256':hashlib.sha256(DEST.read_bytes()).hexdigest(),'paragraphs_revised':sorted(changes),'images_replaced':list(imgs)},ensure_ascii=False,indent=2))
print(DEST)
