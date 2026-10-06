#!/usr/bin/env python3
"""Package explicitly selected current versions; put superseded targets in history."""
from pathlib import Path
import csv
import hashlib
import json
import shutil

root = Path(__file__).resolve().parents[2]
out = root / 'outputs/20261006_figures_latest'
dest = Path('/Users/gzy2520/Desktop/乳酸化全部图片_当前版本_20261006')
dest.mkdir(exist_ok=True)
manifest, decisions = [], []

def csv_write(path, headers, rows):
    with path.open('w', encoding='utf-8-sig', newline='') as f:
        writer = csv.writer(f, lineterminator='\n')
        writer.writerow(headers)
        writer.writerows(rows)

def collect_file(source, relative, status, evidence):
    target = dest / relative
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, target)
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    manifest.append([str(relative), status, str(source), digest, evidence])

def collect_pair(source_root, stem, label, status, evidence):
    for ext in ('png', 'pdf'):
        source = source_root / (stem + '.' + ext)
        assert source.is_file(), source
        collect_file(source, Path(label) / (stem + '.' + ext), status, evidence)

def collect_tree(source, label, status, evidence):
    for p in sorted(source.rglob('*')):
        if p.suffix.lower() in ('.png', '.pdf'):
            collect_file(p, Path(label) / p.relative_to(source), status, evidence)

formal = root / 'corrected_final_result_20260905_sample_only'
evidence = 'sample-only 修复发布 c349c7f；与原 figure_refine-fork、RNA 工作树及 20260917 蛋白组复制版逐图一致'
main = [
 'Figure_1/Figure_1a_DDR_fraction_boxplot',
 'Figure_1/Figure_1b_MKI67_over_H3C1_boxplot',
 'Figure_2/Figure_2a_whole_proteome_DDR_UpSet',
 'Figure_2/Figure_2b_Kla_DDR_UpSet',
 'Figure_3/Figure_3a_reference_regulator_percentiles',
 'Figure_3/Figure_3b_Kla_regulator_percentiles',
 'Supplementary_Figure_S1/Figure_S1a_DDR_fraction_by_PXD',
 'Supplementary_Figure_S1/Figure_S1b_MKI67_over_H3C1_by_PXD',
 'Supplementary_Figure_S2/Figure_S2a_whole_proteome_UpSet',
 'Supplementary_Figure_S2/Figure_S2b_Kla_proteome_UpSet',
 'Supplementary_Figure_S3/Figure_S3a_DDR_pathway_matrix_tumor_tissues',
 'Supplementary_Figure_S3/Figure_S3b_DDR_pathway_matrix_non_tumor_tissues',
 'Supplementary_Figure_S4/Figure_S4a_DDR_pathway_matrix_cancer_cell_lines',
 'Supplementary_Figure_S4/Figure_S4b_DDR_pathway_matrix_normal_cell_lines',
]
for letter, pathway in zip('cdefghi', ['BER', 'NER', 'MMR', 'FA', 'HR', 'AEJ', 'NHEJ']):
    main.append(f'Figure_2/Figure_2{letter}_DDR_pathway_summary_{pathway}_barplot')
for stem in main:
    collect_pair(formal, stem, '01_当前蛋白组主图与附图', '当前', evidence)
alternatives = [
 'Figure_1/companion_controls/Figure_1_companion_MKI67_over_ACTB_boxplot',
 'Figure_1/companion_controls/Figure_1_companion_MKI67_over_TUBB_boxplot',
 'Figure_3/Figure_3a_reference_regulator_percentiles_no_frame',
 'Figure_3/Figure_3b_Kla_regulator_percentiles_no_frame',
 'Supplementary_Figure_S1/Figure_S1b_MKI67_over_H3C1_11_detected_only',
]
for stem in alternatives:
    collect_pair(formal, stem, '01_当前蛋白组主图与附图/有效备选与对照', '有效备选', evidence)
historical = [
 'Supplementary_Figure_S3/legacy_tissue_summaries/Figure_2d_DDR_pathway_summary_tumor_tissue',
 'Supplementary_Figure_S3/legacy_tissue_summaries/Figure_2e_DDR_pathway_summary_non_tumor_tissue',
 'Supplementary_Figure_S4/companion_summaries/Supplementary_Figure_S2b_DDR_pathway_summary_cancer_cell_lines',
 'Supplementary_Figure_S4/companion_summaries/Supplementary_Figure_S2c_DDR_pathway_summary_normal_cell_lines',
]
for stem in historical:
    collect_pair(formal, stem, '90_历史版本_仅供追溯/蛋白组旧配套图', '历史', '原包保留的配套视图；当前主目录采用更新的 sample-only 七通路图，不作为最新主图')
decisions.append(['蛋白组主附图', str(formal), '当前', evidence,
                  '不重复收录七通路无面板编号别名、ACTB/TUBB 重复别名和蓝框重复别名；旧配套图归档'])

with (root / 'config/rna_current_release.csv').open() as f:
    release = dict(csv.reader(f))
rna = root / release['delivery']
assert release['delivery'] == 'outputs/20260919_repair_release'
collect_tree(rna / 'core', '02_当前RNA/核心图', '当前', 'config/rna_current_release.csv 指定的修复发布；28 个独立 RNA 材料、1,898 条记录')
collect_tree(rna / 'reference', '02_当前RNA/蛋白组映射参考图', '当前参考视图', '同一当前修复发布；31 行是蛋白组映射，部分共享 RNA，不是 31 份独立 RNA 材料')
decisions.append(['RNA', str(rna), '当前', 'config/rna_current_release.csv', '排除修复前 31/28 混版目录；31 行映射参考视图单独标注'])

lactate = root / 'outputs/20260920_lactate_metabolism/figures'
collect_tree(lactate, '03_乳酸代谢', '当前探索结果', '已冻结的第一阶段乳酸来源/去向表达分析；release_sha256.csv')
decisions.append(['乳酸代谢', str(lactate), '当前探索结果', 'outputs/20260920_lactate_metabolism/release_sha256.csv', '不收录 pre_labels、pre_mapping_qa、pre_visual 试跑目录'])

ml = out / 'ML'
collect_tree(ml / 'current_four_targets', '04_当前机器学习_四个组学内目标', '当前',
             '20260929 四目标训练结果；最新版 render_all_extended_figures.R 与 summarize.R 的图内容；20261006 零起点重绘')
collect_tree(ml / 'history_cross_assay', '90_历史版本_仅供追溯/机器学习_跨组学两目标', '历史目标的有效修订版',
             '内部预测沿用 20260927；展示术语和外部结果以 20260928_report_revision 为准，外部 RNA 使用 log2(TPM+0.5)')
collect_tree(ml / 'history_total_Kla', '90_历史版本_仅供追溯/机器学习_全蛋白Kla占比', '历史目标',
             '20260926 全 RNA 模型与 generate_extended_figures.R；目标已由导师改为四个组学内目标')
decisions.extend([
 ['四目标机器学习', 'outputs/20260929_assay_composition', '当前',
  '9443340 训练；ba51602 扩展图；training_importance_annotated.csv；comparison_metrics.csv',
  '保留 Symbol+Ensembl、材料实测参考菱形、三面板 UMAP、训练/测试颜色、四目标汇总和基线指标'],
 ['两目标机器学习', 'outputs/20260928_report_revision', '历史目标的有效修订版',
  'c8e5b25；workflow/report_revision/README.md',
  '旧外部 log2(TPM+1) 图不收录；历史结果不混入四目标当前结果'],
 ['全蛋白占比机器学习', 'outputs/20260926_full_rna_models', '历史目标',
  '1e12b82；generate_extended_figures.R；feature_importance_top25.csv',
  'Top25 保留最新版全量拟合模型重要性，明确其不同于当前四目标的训练子集重要性'],
 ['20261004 图片包', 'outputs/20261004_all_figures', '已被本包替代',
  '79187de；本次版本核对',
  '旧重绘丢失符号/参考菱形/分类 UMAP 面板并漏收四目标汇总；本包恢复这些内容'],
])

csv_write(dest / '图片索引与来源.csv', ['包内路径', '版本状态', '原始来源', 'SHA256', '选择依据'], manifest)
csv_write(dest / '版本选择依据.csv', ['图组', '权威来源', '状态', '版本依据', '选择或替换说明'], decisions)
csv_write(out / 'selected_figure_manifest.csv', ['包内路径', '版本状态', '原始来源', 'SHA256', '选择依据'], manifest)
csv_write(out / 'version_selection.csv', ['图组', '权威来源', '状态', '版本依据', '选择或替换说明'], decisions)
(dest / '先看这里_版本说明.md').write_text('''# 乳酸化图片：当前版本汇总（2026-10-06）

本包替代 2026-10-04 图片包。选图以发布清单、实际输入和最新生成脚本为依据，不按文件修改时间简单排序。

- 01：当前蛋白组 Figure 1–3、Figure S1–S4。每个面板只保留一个标准文件名；无框热图、ACTB/TUBB 对照和紧凑 S1b 放在“有效备选与对照”。
- 02：当前修复版 RNA 图。核心组图计 28 个独立 RNA 材料；蛋白组映射参考视图可显示 31 行，不能据此声称 31 份独立 RNA。
- 03：第一阶段乳酸代谢表达分析。
- 04：当前机器学习四个目标，分别是普通蛋白组 DDR/全部普通蛋白、Kla 组 DDR/全部 Kla、普通蛋白组 DNA repair/全部普通蛋白、Kla 组 DNA repair/全部 Kla。每项都在同一组学内计数。
- 90：历史模型和旧配套图，保留用于追溯；不能将它们的标签、指标或外部结果用于当前四目标结论。

## 本次修正

机器学习保留最新图的 Symbol + Ensembl 标签、材料参考菱形、类别/参考值/拟合值三面板 UMAP、训练与测试颜色，并补齐四目标内部留出、来源留出和外部比较概览。Top25 的回归重要性使用方差减少术语，不标为分类 Gini 重要性。

本次没有训练模型或修改预测值。占比、预测值、特征重要性、绝对误差及记录数等非负轴从 0 开始。分类轴保留标签。带正负号的残差另放 signed_diagnostic_appendix，保留负值；裁掉负值会隐藏低估误差。

UMAP 从同一矩阵按最新脚本的特征/主成分/邻居参数、seed=25 计算，计算后分别减去两轴最小值，使坐标从 0 开始。保存计算前和计算后的坐标；平移不改变该次嵌入的点间距离。因原图没有保存坐标且运行线程可能不同，不声称新嵌入与旧图片逐点坐标完全相同。UMAP 全数据拟合着色仅作描述，不代表验证性能。

外部新 RNA 使用已校正的 log2(TPM+0.5)；质谱与 RNA 是组织层级比较，并非患者配对。RNA 治疗前史未完整确认。复用 TCGA RNA 的图单独标明，不能当作独立外部 RNA 测试。

全部图提供 PNG/PDF。图片索引列出每个文件的版本状态、源路径和 SHA256；版本选择依据表记录当前/历史判定。PNG 适合查看与插图，PDF 适合排版与放大。
''', encoding='utf-8')
scripts = dest / '重绘脚本'
scripts.mkdir(exist_ok=True)
for name in ['render_latest_ml_zero.R', 'package_latest_figures.py', 'README.md']:
    shutil.copy2(root / 'workflow/figure_delivery' / name, scripts / name)
shutil.copy2(ml / 'axis_lower_bound_audit.csv', dest / '机器学习坐标起点检查.csv')
shutil.copy2(root / 'audit/20261006_figure_versions/REPORT.md', dest / '版本核对记录.md')
counts = {}
for relative, *_ in manifest:
    key = Path(relative).parts[0]
    counts[key] = counts.get(key, 0) + 1
summary = {'destination': str(dest), 'files': len(manifest), 'counts': counts}
(out / 'delivery_summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding='utf-8')
shutil.make_archive(str(dest), 'zip', root_dir=dest.parent, base_dir=dest.name)
print(json.dumps(summary, ensure_ascii=False))
