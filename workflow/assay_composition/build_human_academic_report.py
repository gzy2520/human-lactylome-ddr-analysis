#!/usr/bin/env python3
"""
build_human_academic_report.py
Generates a natural, rigorous, professional Chinese academic research report (.docx)
for the four within-assay composition targets (Figure 1a definition).
Completely removes AI meta-commentary, coaching language, and exaggerated rhetoric.
Strictly adheres to Chinese academic thesis/report standards and OOXML validation.
"""

from pathlib import Path
import csv, json, os
import docx
from docx.shared import Pt, RGBColor, Cm
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.oxml import parse_xml, OxmlElement
from docx.oxml.ns import nsdecls, qn

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'outputs/20260929_assay_composition'

def read_csv(name):
    return list(csv.DictReader((OUT / name).open()))

def create_report():
    doc = docx.Document()

    # 1. Page setup: A4 with 2.5 cm margins
    section = doc.sections[0]
    section.page_width = Cm(21.0)
    section.page_height = Cm(29.7)
    section.top_margin = Cm(2.5)
    section.bottom_margin = Cm(2.5)
    section.left_margin = Cm(2.5)
    section.right_margin = Cm(2.5)

    # Fix settings.xml zoom percent to satisfy OOXML XSD validation
    settings = doc.settings.element
    zoom = settings.find(qn('w:zoom'))
    if zoom is not None:
        zoom.set(qn('w:percent'), '100')
    else:
        zoom_elem = parse_xml(r'<w:zoom %s w:percent="100"/>' % nsdecls('w'))
        settings.append(zoom_elem)

    # 2. Header and Footer
    header = section.header
    hp = header.paragraphs[0]
    hp.alignment = WD_ALIGN_PARAGRAPH.RIGHT
    hrun = hp.add_run("人类乳酸化修饰组学与 DDR 机制分析 · 实验研究报告")
    hrun.font.name = "SimSun"
    hrun.font.size = Pt(8.5)
    hrun.font.color.rgb = RGBColor(128, 128, 128)

    footer = section.footer
    fp = footer.paragraphs[0]
    fp.alignment = WD_ALIGN_PARAGRAPH.CENTER
    f_run = fp.add_run("第 ")
    f_run.font.name = "SimSun"
    f_run.font.size = Pt(9)
    f_run.font.color.rgb = RGBColor(100, 100, 100)

    fldSimple1 = OxmlElement('w:fldSimple')
    fldSimple1.set(qn('w:instr'), 'PAGE')
    fp._p.append(fldSimple1)

    f_run2 = fp.add_run(" 页")
    f_run2.font.name = "SimSun"
    f_run2.font.size = Pt(9)
    f_run2.font.color.rgb = RGBColor(100, 100, 100)

    # 3. Helper functions
    def set_run_font(run, font_chinese="SimSun", font_ascii="Arial", size_pt=11, bold=False, italic=False, color_rgb=None):
        run.bold = bold
        run.italic = italic
        run.font.size = Pt(size_pt)
        if color_rgb:
            run.font.color.rgb = color_rgb
        run.font.name = font_ascii
        rPr = run._r.get_or_add_rPr()
        rFonts = parse_xml(r'<w:rFonts %s w:ascii="%s" w:hAnsi="%s" w:eastAsia="%s"/>' %
                           (nsdecls('w'), font_ascii, font_ascii, font_chinese))
        rPr.append(rFonts)

    def add_title(text):
        p = doc.add_paragraph()
        p.alignment = WD_ALIGN_PARAGRAPH.CENTER
        p.paragraph_format.space_before = Pt(6)
        p.paragraph_format.space_after = Pt(4)
        p.paragraph_format.line_spacing = 1.2
        run = p.add_run(text)
        set_run_font(run, font_chinese="SimHei", font_ascii="Arial", size_pt=16, bold=True, color_rgb=RGBColor(20, 30, 50))
        return p

    def add_subtitle(text):
        p = doc.add_paragraph()
        p.alignment = WD_ALIGN_PARAGRAPH.CENTER
        p.paragraph_format.space_before = Pt(0)
        p.paragraph_format.space_after = Pt(12)
        run = p.add_run(text)
        set_run_font(run, font_chinese="SimSun", font_ascii="Arial", size_pt=10.5, color_rgb=RGBColor(100, 100, 100))
        return p

    def add_meta_box(meta_dict):
        table = doc.add_table(rows=len(meta_dict), cols=2)
        table.alignment = WD_TABLE_ALIGNMENT.CENTER
        tblPr = table._tbl.tblPr
        borders_xml = parse_xml(r'''
            <w:tblBorders %s>
                <w:top w:val="single" w:sz="4" w:space="0" w:color="CCCCCC"/>
                <w:left w:val="none"/>
                <w:bottom w:val="single" w:sz="4" w:space="0" w:color="CCCCCC"/>
                <w:right w:val="none"/>
                <w:insideH w:val="none"/>
                <w:insideV w:val="none"/>
            </w:tblBorders>
        ''' % nsdecls('w'))
        tblLook = tblPr.find(qn('w:tblLook'))
        if tblLook is not None:
            tblLook.addprevious(borders_xml)
        else:
            tblPr.append(borders_xml)

        for i, (k, v) in enumerate(meta_dict.items()):
            row = table.rows[i]
            row.cells[0].width = Cm(3.2)
            row.cells[1].width = Cm(12.8)
            p0 = row.cells[0].paragraphs[0]
            p0.alignment = WD_ALIGN_PARAGRAPH.RIGHT
            p0.paragraph_format.space_before = Pt(2)
            p0.paragraph_format.space_after = Pt(2)
            r0 = p0.add_run(k + "：")
            set_run_font(r0, font_chinese="SimHei", font_ascii="Arial", size_pt=9.5, bold=True, color_rgb=RGBColor(70, 70, 70))

            p1 = row.cells[1].paragraphs[0]
            p1.alignment = WD_ALIGN_PARAGRAPH.LEFT
            p1.paragraph_format.space_before = Pt(2)
            p1.paragraph_format.space_after = Pt(2)
            r1 = p1.add_run(v)
            set_run_font(r1, font_chinese="SimSun", font_ascii="Arial", size_pt=9.5, color_rgb=RGBColor(40, 40, 40))

        sp = doc.add_paragraph()
        sp.paragraph_format.space_before = Pt(0)
        sp.paragraph_format.space_after = Pt(6)

    def add_h1(text):
        p = doc.add_paragraph()
        p.paragraph_format.space_before = Pt(14)
        p.paragraph_format.space_after = Pt(5)
        p.paragraph_format.keep_with_next = True
        run = p.add_run(text)
        set_run_font(run, font_chinese="SimHei", font_ascii="Arial", size_pt=13, bold=True, color_rgb=RGBColor(25, 55, 95))
        return p

    def add_h2(text):
        p = doc.add_paragraph()
        p.paragraph_format.space_before = Pt(9)
        p.paragraph_format.space_after = Pt(3)
        p.paragraph_format.keep_with_next = True
        run = p.add_run(text)
        set_run_font(run, font_chinese="SimHei", font_ascii="Arial", size_pt=11.5, bold=True, color_rgb=RGBColor(40, 75, 120))
        return p

    def add_body(text, bold_prefix=""):
        p = doc.add_paragraph()
        p.paragraph_format.first_line_indent = Cm(0.74)
        p.paragraph_format.line_spacing = 1.25
        p.paragraph_format.space_before = Pt(0)
        p.paragraph_format.space_after = Pt(3)
        if bold_prefix:
            r_pre = p.add_run(bold_prefix)
            set_run_font(r_pre, font_chinese="SimSun", font_ascii="Arial", size_pt=10.5, bold=True)
        r_text = p.add_run(text)
        set_run_font(r_text, font_chinese="SimSun", font_ascii="Arial", size_pt=10.5)
        return p

    def add_bullet(text, bold_prefix=""):
        p = doc.add_paragraph()
        p.paragraph_format.left_indent = Cm(0.74)
        p.paragraph_format.first_line_indent = Cm(-0.4)
        p.paragraph_format.line_spacing = 1.25
        p.paragraph_format.space_before = Pt(1)
        p.paragraph_format.space_after = Pt(2)
        r_bullet = p.add_run("▪ ")
        set_run_font(r_bullet, font_chinese="SimSun", font_ascii="Arial", size_pt=9.5, bold=True, color_rgb=RGBColor(31, 78, 121))
        if bold_prefix:
            r_pre = p.add_run(bold_prefix)
            set_run_font(r_pre, font_chinese="SimSun", font_ascii="Arial", size_pt=10.5, bold=True)
        r_text = p.add_run(text)
        set_run_font(r_text, font_chinese="SimSun", font_ascii="Arial", size_pt=10.5)
        return p

    def add_caption(text, is_table=False):
        p = doc.add_paragraph()
        p.alignment = WD_ALIGN_PARAGRAPH.CENTER
        p.paragraph_format.space_before = Pt(6) if is_table else Pt(4)
        p.paragraph_format.space_after = Pt(4) if is_table else Pt(8)
        p.paragraph_format.first_line_indent = Cm(0)
        p.paragraph_format.keep_with_next = is_table
        run = p.add_run(text)
        set_run_font(run, font_chinese="SimSun", font_ascii="Arial", size_pt=10, bold=True, color_rgb=RGBColor(50, 50, 50))
        return p

    def add_figure(img_path, caption_text, width_cm=15.5):
        if not os.path.exists(img_path):
            print(f"[Warning] Image not found: {img_path}")
            return
        p_img = doc.add_paragraph()
        p_img.alignment = WD_ALIGN_PARAGRAPH.CENTER
        p_img.paragraph_format.space_before = Pt(6)
        p_img.paragraph_format.space_after = Pt(2)
        p_img.paragraph_format.first_line_indent = Cm(0)
        p_img.paragraph_format.keep_with_next = True
        run_img = p_img.add_run()
        run_img.add_picture(str(img_path), width=Cm(width_cm))
        add_caption(caption_text, is_table=False)

    def add_three_line_table(headers, rows, col_widths_cm):
        n_cols = len(headers)
        table = doc.add_table(rows=len(rows) + 1, cols=n_cols)
        table.alignment = WD_TABLE_ALIGNMENT.CENTER
        tblPr = table._tbl.tblPr

        borders_xml = parse_xml(r'''
            <w:tblBorders %s>
                <w:top w:val="single" w:sz="12" w:space="0" w:color="000000"/>
                <w:left w:val="none"/>
                <w:bottom w:val="single" w:sz="12" w:space="0" w:color="000000"/>
                <w:right w:val="none"/>
                <w:insideH w:val="none"/>
                <w:insideV w:val="none"/>
            </w:tblBorders>
        ''' % nsdecls('w'))
        tblLook = tblPr.find(qn('w:tblLook'))
        if tblLook is not None:
            tblLook.addprevious(borders_xml)
        else:
            tblPr.append(borders_xml)

        hdr_row = table.rows[0]
        hdr_row._tr.get_or_add_trPr().append(parse_xml(r'<w:tblHeader %s/>' % nsdecls('w')))
        for c_idx, h_text in enumerate(headers):
            cell = hdr_row.cells[c_idx]
            cell.width = Cm(col_widths_cm[c_idx])
            tcPr = cell._tc.get_or_add_tcPr()
            tcBorders = parse_xml(r'''
                <w:tcBorders %s>
                    <w:bottom w:val="single" w:sz="6" w:space="0" w:color="000000"/>
                </w:tcBorders>
            ''' % nsdecls('w'))
            tcPr.append(tcBorders)
            shd = parse_xml(r'<w:shd %s w:val="clear" w:color="auto" w:fill="F5F7FA"/>' % nsdecls('w'))
            tcPr.append(shd)

            p = cell.paragraphs[0]
            p.alignment = WD_ALIGN_PARAGRAPH.CENTER
            p.paragraph_format.space_before = Pt(3)
            p.paragraph_format.space_after = Pt(3)
            p.paragraph_format.first_line_indent = Cm(0)
            run = p.add_run(h_text)
            set_run_font(run, font_chinese="SimHei", font_ascii="Arial", size_pt=9.5, bold=True)

        for r_idx, row_data in enumerate(rows):
            row = table.rows[r_idx + 1]
            for c_idx, val in enumerate(row_data):
                cell = row.cells[c_idx]
                cell.width = Cm(col_widths_cm[c_idx])
                p = cell.paragraphs[0]
                p.alignment = WD_ALIGN_PARAGRAPH.CENTER if c_idx > 0 else WD_ALIGN_PARAGRAPH.LEFT
                p.paragraph_format.space_before = Pt(2.5)
                p.paragraph_format.space_after = Pt(2.5)
                p.paragraph_format.first_line_indent = Cm(0)
                run = p.add_run(str(val))
                set_run_font(run, font_chinese="SimSun", font_ascii="Arial", size_pt=9.5)

        sp = doc.add_paragraph()
        sp.paragraph_format.space_before = Pt(0)
        sp.paragraph_format.space_after = Pt(6)

    # Read data
    m = read_csv('comparison_metrics.csv')
    ext = read_csv('external_metrics.csv')
    targets = ['Proteome_DDR', 'Kla_DDR', 'Proteome_DNA_repair', 'Kla_DNA_repair']
    target_names = {
        'Proteome_DDR': '普通组 · DDR 占比',
        'Kla_DDR': 'Kla 组 · DDR 占比',
        'Proteome_DNA_repair': '普通组 · DNA repair 占比',
        'Kla_DNA_repair': 'Kla 组 · DNA repair 占比'
    }

    # Document Header
    add_title("RNA 表达预测组学内 DDR 与 DNA 修复蛋白占比模型分析报告")
    add_subtitle("图 1a 计数口径下的四目标非线性建模与独立队列材料级对照验证")

    add_meta_box({
        "研究课题": "人类乳酸化修饰组学与 DDR 机制交互分析",
        "预测目标": "组学内部构成比 (普通蛋白组与 Kla 蛋白组各自内部的 DDR / DNA repair 检出占比)",
        "计数口径": "严格对齐论文图 1a 定义 (分子与分母同属单一组学，不求跨组学交集)",
        "输入特征": "全转录组 18,332 个共有 Ensembl 稳定基因 (log2(TPM + 0.5) 表达矩阵)",
        "内部验证设计": "1,898 条 RNA (1,524 训练 / 374 内部留出，供者隔离，种子 25)",
        "外部对照队列": "质谱: PXD063945 (华山医院 ESCC)；RNA: GSE130078 (韩国 YSH 队列 46 例配对组织)",
        "报告日期": "2026 年 9 月 30 日"
    })

    # 一、研究背景与主要结论
    add_h1("一、研究背景与主要结论")
    add_body(
        "在前期探索中，预测目标曾被设定为‘普通蛋白组检出的通路蛋白中同时发生乳酸化的比例’。"
        "然而，该跨组学交集比值的分子依赖抗体富集质谱、分母依赖常规全蛋白组质谱，"
        "受制于两套实验检测深度的交乘影响，在跨实验室验证时容易出现系统性基线偏移，且与论文图 1a 的统计口径存在出入。"
    )
    add_body(
        "为此，本轮分析将预测目标调整为与论文图 1a 完全一致的【组学内部组成比例（Within-Assay Composition）】："
        "分别在普通蛋白质组和乳酸化修饰组内部，计算 DDR 及 DNA 修复相关蛋白占各自总检出蛋白的百分比。"
        "四个目标各自独立训练随机森林回归模型，输入特征保持为全转录组 18,332 个共有 Ensembl 基因。"
    )
    add_body("四个预测目标的定义如下：")
    add_bullet("普通组 DDR 占比 (Proteome_DDR)：普通蛋白组中 DDR 检出蛋白数 / 全部普通组检出蛋白数 × 100；", bold_prefix="1. ")
    add_bullet("Kla 组 DDR 占比 (Kla_DDR)：Kla 蛋白组中 DDR 检出蛋白数 / 全部 Kla 组检出蛋白数 × 100；", bold_prefix="2. ")
    add_bullet("普通组 DNA repair 占比 (Proteome_DNA_repair)：普通蛋白组中 DNA 修复蛋白数 / 全部普通组检出蛋白数 × 100；", bold_prefix="3. ")
    add_bullet("Kla 组 DNA repair 占比 (Kla_DNA_repair)：Kla 蛋白组中 DNA 修复蛋白数 / 全部 Kla 组检出蛋白数 × 100。", bold_prefix="4. ")

    add_body(
        "核心评估结果表明：",
        bold_prefix="主要结论概括："
    )
    add_bullet("内部留出表现优异：在 374 条未参与训练的同类材料 RNA 记录中，四个目标的逐记录 MAE 均在 0.03～0.20 个百分点，R² 为 0.94～0.98，超过 96% 的记录误差在 ±1 个百分点以内；", bold_prefix="内部留出一致性高：")
    add_bullet("外部肿瘤材料对照高度吻合：在外部独立的食管鳞癌未治疗肿瘤组织中，四个目标的平均预测值与新质谱真值的偏差均在 0.05～0.19 个百分点，显著优于无 RNA 基线；", bold_prefix="外部肿瘤组织验证良好：")
    add_bullet("外部癌旁预测存在一定偏离：在外部配对癌旁切片中，Kla 组的预测误差为 1.0～1.8 个百分点，反映出癌旁正常组织的修饰基线与模型在肿瘤中学习到的规律存在一定差异；", bold_prefix="组织类型差异客观存在：")
    add_bullet("模型定位于组织水平参考重现：每条 RNA 记录对应的是所属材料的质谱共享参考标签，当前结果支持在材料层面重现通路占比，尚不能证明可用于单人水平的个体差异预测。", bold_prefix="适用边界明确：")

    # 表 1: 四目标综合性能
    rows_tbl1 = []
    for t in targets:
        r_tr = next(r for r in m if r['Target'] == t and r['Evaluation'] == 'Training split' and r['Weighting'] == 'Record')
        r_in = next(r for r in m if r['Target'] == t and r['Evaluation'] == 'Within-material test' and r['Weighting'] == 'Record')
        r_in_m = next(r for r in m if r['Target'] == t and r['Evaluation'] == 'Within-material test' and r['Weighting'] == 'Material')
        r_so_m = next(r for r in m if r['Target'] == t and r['Evaluation'] == 'Source held out' and r['Weighting'] == 'Material')
        rows_tbl1.append([
            target_names[t],
            f"{float(r_tr['MAE']):.3f}",
            f"{float(r_in['MAE']):.3f}",
            f"{float(r_in['R2']):.3f}",
            f"{float(r_in_m['MAE']):.3f}",
            f"{float(r_in_m['Within1']):.1f}%",
            f"{float(r_so_m['MAE']):.3f}",
            f"{float(r_so_m['BaselineMAE']):.3f}"
        ])

    add_caption("表 1 四大组学内部通路占比模型在不同验证阶段的性能指标汇总表", is_table=True)
    add_three_line_table(
        headers=["预测目标", "训练 MAE", "内部记录 MAE", "内部 R²", "内部材料 MAE", "内部 ±1pp 率", "跨来源 MAE", "无 RNA 基线"],
        rows=rows_tbl1,
        col_widths_cm=[3.8, 1.6, 1.8, 1.6, 1.8, 1.8, 1.8, 1.8]
    )

    # 二、数据架构与计数口径
    add_h1("二、数据架构与计数口径")
    add_h2("2.1 图 1a 计数口径的严格复现")
    add_body(
        "在既有数据集的统计中，普通蛋白组保留原始 SourceProteinID 计数单位（包含部分 Ensembl protein），"
        "当一个源 ID 映射到多个 UniProt 时，任一映射命中目标 GO 集合即判定为命中，但源 ID 只计一次，"
        "以保证与冻结的 group_summary_31.csv 基准一致。Kla 组统一按去重 BaseAccession 计数。"
        "所有指标的分子均为分母的真子集，数值严格有界于 0%～100%。"
    )
    add_body(
        "DDR 通路采用 GO:0006974 及其全部 is_a / part_of 后代定义（包含 2,660 个 UniProt BaseAccession）；"
        "DNA repair 采用 GO:0006281 及其后代定义（包含 2,135 个 UniProt BaseAccession）。"
        "GO 注释仅用于界定输出标签，不作为输入基因特征的筛选条件。"
    )

    add_h2("2.2 建模参数与样本划分")
    add_body(
        "输入特征为各数据来源共有的全部 18,332 个 Ensembl 基因，统一采用修正后的 log2(TPM + 0.5) 尺度，"
        "不使用经验分位数平滑（qsmooth）。回归模型统一使用随机森林（ranger 包实现），"
        "设置 500 棵回归树，mtry = 135，最小叶节点大小为 1，随机种子固定为 25。"
        "样本权重综合平衡材料与供者数量，防止样本规模较大的组织主导特征分裂。"
    )
    add_body(
        "验证方案分为两部分：一是供者隔离的内部留出验证（1,524 训练 / 374 留出，来自同一供者的样本不跨两侧）；"
        "二是 14 折跨来源关联块留出验证（每次留出一个来源块作为盲测，循环 14 次）。"
    )

    # 三、模型训练与内部留出结果
    add_h1("三、模型训练与内部留出结果")
    add_body(
        "图 1 展示了四个模型在 374 条独立内部测试样本上的预测散点分布；图 2 展示了 14 折来源留出下的表现。"
    )
    add_figure(
        "outputs/20260929_assay_composition/figures/internal_test.png",
        "图 1 四大组学内部通路占比模型在 374 条独立内部留出样本上的预测散点图",
        width_cm=15.5
    )
    add_figure(
        "outputs/20260929_assay_composition/figures/source_heldout.png",
        "图 2 四大模型在 14 折来源关联块留出 (Leave-One-Source-Out) 下的预测对照图",
        width_cm=15.5
    )
    add_body(
        "结果显示：在内部留出测试中，所有样本点紧密贴近对角虚线（y = x），"
        "逐记录平均绝对误差在 0.03～0.20 个百分点之间，判定系数 R² 达到 0.94～0.98，"
        "表明模型在已见过的材料类型中能够极好地重现该材料对应的通路占比标签。"
    )
    add_body(
        "在跨来源关联块留出测试中，Kla 组的两个目标在材料等权口径下的 MAE 分别为 1.61 和 1.05 个百分点，"
        "均显著优于不使用 RNA 的简单中位数基线（3.33 和 2.45 个百分点），误差减少幅度为 51.7%～57.1%。"
        "这说明 RNA 表达特征即使在来源留出条件下也能提供一定的有效信息，但整体误差相较内部测试有所上升。"
    )

    # 四、外部独立队列的材料级对照验证
    add_h1("四、外部独立队列的材料级对照验证")
    add_body(
        "为评估模型跨实验队列的泛化表现，我们检索并引入了 2026 年新公布的华山医院食管鳞癌蛋白质组与乳酸化组数据（PXD063945 / PMC13195773），"
        "以及韩国 YSH 队列 46 例配对食管癌切片组织 RNA 数据（GSE130078）。"
        "在完全冻结四个模型的前提下，直接输入外部 RNA 并与新质谱计算出的材料参考值进行比对。"
    )
    add_body(
        "预测散点与分布如图 3 所示，具体数值详见表 2。"
    )
    add_figure(
        "outputs/20260929_assay_composition/figures/external_newRNA.png",
        "图 3 冻结模型在外部独立食管鳞癌队列 (GSE130078 组织 RNA vs PXD063945 质谱参考值) 中的预测分布图",
        width_cm=15.5
    )

    # 表 2: 外部验证详细指标
    rows_tbl2 = []
    for t in targets:
        for mat, cn in [('ESCC_untreated', '肿瘤组织'), ('ESCC_adjacent_untreated', '配对癌旁')]:
            r_ext = next(r for r in ext if r['Scope'] == 'newRNA' and r['Target'] == t and r['Material'] == mat)
            rows_tbl2.append([
                target_names[t],
                cn,
                f"{float(r_ext['Reference']):.3f}%",
                f"{float(r_ext['MeanPrediction']):.3f}%",
                f"{float(r_ext['MAE']):.3f}",
                f"{float(r_ext['BaselineMAE']):.3f}",
                f"{float(r_ext['Within1']):.1f}%"
            ])

    add_caption("表 2 外部独立队列 (GSE130078 组织 RNA vs PXD063945 质谱参考值) 性能对照表", is_table=True)
    add_three_line_table(
        headers=["预测目标", "组织类型", "质谱材料参考值", "模型预测均值", "材料比较 MAE", "无 RNA 基线", "±1 pp 命中率"],
        rows=rows_tbl2,
        col_widths_cm=[4.0, 2.0, 2.2, 2.2, 2.0, 2.0, 1.6]
    )

    add_body(
        "从外部验证结果中可以观察到两个鲜明现象：",
        bold_prefix="对比观察："
    )
    add_bullet(
        "在 ESCC 肿瘤组织中，四个目标的模型预测均值与新质谱参考值高度吻合："
        "Kla-DDR 参考值为 6.831%，预测均值为 6.750%（误差 0.190 个百分点）；"
        "Kla-Repair 参考值为 4.797%，预测均值为 4.730%（误差 0.159 个百分点）；"
        "普通组两目标的误差仅为 0.143 和 0.052 个百分点，所有 23 条肿瘤 RNA 样本的误差均在 ±1 个百分点以内，且显著优于基线；",
        bold_prefix="肿瘤预测一致性高：")
    add_bullet(
        "在配对癌旁组织中，Kla 组两个目标的预测均值出现了一定程度的偏低："
        "Kla-DDR 参考值为 8.898%，模型预测均值为 7.017%（误差 1.881 个百分点）；"
        "Kla-Repair 参考值为 5.932%，模型预测均值为 4.899%（误差 1.033 个百分点），未优于无 RNA 基线。",
        bold_prefix="癌旁预测偏低：")

    # 五、结果分析与适用边界
    add_h1("五、结果分析与适用边界")
    add_h2("5.1 为什么组学内占比在外部肿瘤中吻合度更好？")
    add_body(
        "此前基于两组学交集比值的预测模型在外部验证中出现约 9 个百分点的偏移，"
        "主要原因在于抗体富集深度差异改变了两套质谱检出列表的交集大小。"
        "而本次采用的【组学内部组成比例】仅依赖单一组学自身的检出总量作为分母，"
        "反映的是恶性肿瘤细胞总蛋白池或总修饰池中分配给特定修复通路的固有比例。"
        "这类基础组分构成在同类肿瘤细胞中表现出较高的稳态特性，因而减少了跨实验室富集深度波动带来的影响。"
    )

    add_h2("5.2 癌旁组织偏离的可能原因与客观限制")
    add_body(
        "癌旁组织中 Kla 占比预测偏低，可能与以下因素有关：",
        bold_prefix="潜在限制："
    )
    add_bullet(
        "癌旁组织细胞异质性较高：癌旁切片常混杂间质细胞、浸润免疫细胞及不同程度的慢性炎症改变，"
        "其 RNA 表达谱与纯肿瘤组织存在显著差异，且各实验室取材距肿瘤边缘的距离并不统一；",
        bold_prefix="取材与异质性：")
    add_bullet(
        "训练材料覆盖度不均衡：现有训练集中正常食管组织样本数量极少，模型对食管鳞状上皮在非癌状态下的修饰特征学习不够充分；",
        bold_prefix="正常组织样本偏少：")
    add_bullet(
        "队列非同一患者多组学配对：RNA 来自韩国 YSH 队列，质谱来自华山医院，两者仅在组织病理层面匹配，"
        "无法排除人群遗传背景、生活方式及未知临床治疗史对基线修饰水平的影响。",
        bold_prefix="非个体配对限制：")

    add_h2("5.3 结论与建议")
    add_body(
        "综上所述，当前模型在肿瘤组织内的表现支持其作为‘已知材料体系内基于全转录组特征重现组学内通路占比’的有效工具。"
        "但面对正常/癌旁组织或更广泛的非肿瘤材料时，模型尚不能做到无条件准确预测。"
        "若要在未来开展真正的盲测应用，仍需在固定模型的前提下，引入小批量同批次配对基准数据进行局部标定。"
    )

    # 六、交付文件与数据索引
    add_h1("六、交付文件与数据索引")
    add_caption("表 3 四目标模型交付产出文件路径索引表", is_table=True)
    add_three_line_table(
        headers=["文件类别", "相对路径 / 文件名", "规格说明"],
        rows=[
            ["综合指标表", "outputs/20260929_assay_composition/comparison_metrics.csv", "四目标全景评价指标表 (含基线对照)"],
            ["外部验证指标", "outputs/20260929_assay_composition/external_metrics.csv", "新 RNA 食管癌肿瘤与癌旁逐项指标明细表"],
            ["内部测试图", "outputs/20260929_assay_composition/figures/internal_test.png/.pdf", "四目标内部留出预测 4 联散点图"],
            ["跨来源留出图", "outputs/20260929_assay_composition/figures/source_heldout.png/.pdf", "14 折来源留出交叉验证对比散点图"],
            ["外部盲测图", "outputs/20260929_assay_composition/figures/external_newRNA.png/.pdf", "GSE130078 组织 RNA 盲测外部质谱参考值分布图"],
            ["主报告文件", "reports/图1a口径四目标RNA预测模型及外部盲测实验报告_20260929.docx", "本项目的完整正式学术 Word 报告"],
            ["桌面同步副本", "/Users/gzy2520/Desktop/ML_图1a四目标_20260929.docx", "直接供审阅与交付的桌面端 Word 副本"]
        ],
        col_widths_cm=[2.8, 8.2, 5.0]
    )

    # Save to reports/
    out_dir = Path("reports")
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "图1a口径四目标RNA预测模型及外部盲测实验报告_20260929.docx"
    doc.save(str(out_path))

    # Also save to Desktop
    desktop_path = Path("/Users/gzy2520/Desktop/ML_图1a四目标_20260929.docx")
    doc.save(str(desktop_path))

    # Also save to outputs/ report directory
    out_report_dir = OUT / 'report'
    out_report_dir.mkdir(parents=True, exist_ok=True)
    out_report_docx = out_report_dir / 'ML_图1a四目标_20260929.docx'
    doc.save(str(out_report_docx))

    print(f"Human-academic report successfully saved to:\n  - {out_path}\n  - {desktop_path}\n  - {out_report_docx}")
    return out_path

if __name__ == "__main__":
    create_report()
