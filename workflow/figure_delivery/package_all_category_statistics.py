"""Deliver current individual figures and an explicit per-figure statistics inventory."""
from pathlib import Path
import csv
import hashlib
import shutil
import zipfile

root = Path(__file__).resolve().parents[2]
out = root / 'outputs/20261008_all_category_statistics'
paired = root / 'outputs/20261008_rna_figure2_tumor_reference'
audit = root / 'audit/20261008_all_category_statistics'
dest = Path('/Users/gzy2520/Desktop/乳酸化当前单图_统计补齐_20261008')
slim = Path('/Users/gzy2520/Desktop/统计补齐_15张单图_20261008')
dest.mkdir(exist_ok=True)
slim.mkdir(exist_ok=True)
mapping = {
 'Figure_1a_DDR_fraction_boxplot': out / 'figures/Figure_1a_DDR_fraction_boxplot',
 'Figure_1b_MKI67_over_H3C1_boxplot': out / 'figures/MKI67_over_H3C1_boxplot',
 'Figure_1_companion_MKI67_over_ACTB_boxplot': out / 'figures/MKI67_over_ACTB_boxplot',
 'Figure_1_companion_MKI67_over_TUBB_boxplot': out / 'figures/MKI67_over_TUBB_boxplot',
 'Figure_1a_RNA_DDR_expression_28materials_boxplot': out / 'figures/RNA_DDR_expression_28materials_boxplot',
 'Figure_1a_RNA_DDR_fraction_28materials_boxplot': out / 'figures/RNA_DDR_fraction_28materials_boxplot',
 'Figure_1a_RNA_DDR_and_lactylated_gene_expression_28materials_boxplot': paired / 'figures/Figure_2a_RNA_DDR_and_Kla_gene_expression',
 'Figure_1b_RNA_proliferation_hallmark_28materials_boxplot': paired / 'figures/Figure_2b_RNA_proliferation',
}
for letter, pathway in zip('cdefghi', ['BER', 'NER', 'MMR', 'FA', 'HR', 'AEJ', 'NHEJ']):
    mapping[f'Figure_2{letter}_DDR_pathway_summary_{pathway}_barplot'] = out / f'figures/Kla_{pathway}_barplot'
source_rows = list(csv.DictReader((root / 'outputs/20261006_teacher_figures/delivery_manifest.csv').open(encoding='utf-8-sig')))
manifest, coverage = [], []

def collect(source, relative, note, destination=dest):
    target = destination / relative
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, target)
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    assert hashlib.sha256(target.read_bytes()).hexdigest() == digest
    if destination == dest:
        manifest.append([str(relative), note, str(source), digest])

def disposition(relative):
    name = relative.stem
    if name in mapping:
        return ('已完成统计', '肿瘤组织对其余三类；来源聚类检验和BH校正；不可估计者NE')
    if relative.parts[0].startswith('04_'):
        return ('保留模型评价', '预测/误差/重要性等不作为新增独立生物学实测；保留MAE、R²等既有指标')
    if 'by_PXD' in name or 'detected_only' in name or 'by_tissue_' in name:
        return ('明细图未新增检验', '横轴为具体材料；四类别推断已放在对应汇总图，不把横线连到任意材料')
    if name.startswith('Focus_samples') or name == 'LDHA_LDHB_samples_28':
        return ('描述图未新增检验', '具体器官/材料和供者关系需另定比较设计；不能合并不同器官癌组织作为任意对照')
    if name == 'Kla_detection_context':
        return ('无可用重复检验', '材料汇总/并集检出计数，不是具有独立重复的类别实验')
    return ('不适用类别横线', '集合、热图、机制图或描述性映射/差值；保留图本身含义')

for row in source_rows:
    relative = Path(row['包内路径'])
    if relative.parts[0].startswith(('90_', '00_')) or relative.suffix not in ['.png', '.pdf']:
        continue
    status, reason = disposition(relative)
    source = mapping[relative.stem].with_suffix(relative.suffix) if relative.stem in mapping else Path(row['本次来源'])
    collect(source, relative, status)
    if relative.stem in mapping:
        collect(source, Path(source.name), status, slim)
    if relative.suffix == '.png':
        coverage.append([str(relative), status, reason, str(source)])
assert len(coverage) == 111
assert sum(row[1] == '已完成统计' for row in coverage) == 15
shutil.copy2(paired / 'tables/RNA_proliferation_pairwise_source_clustered.csv',
             out / 'statistics/RNA_proliferation_existing_comparisons.csv')
with (out / '逐图统计处理清单.csv').open('w', encoding='utf-8-sig', newline='') as f:
    w = csv.writer(f, lineterminator='\n')
    w.writerow(['图片路径', '处理状态', '依据与范围', '本版来源'])
    w.writerows(coverage)
extras = [(p, Path('统计结果') / p.name) for p in sorted((out / 'statistics').glob('*.csv'))]
extras += [(audit / 'README.md', Path('交付说明.md')),
           (out / '逐图统计处理清单.csv', Path('逐图统计处理清单.csv')),
           (out / 'qa/render_validation.csv', Path('核对结果.csv')),
           (out / 'sessionInfo.txt', Path('软件版本.txt'))]
for name in ['render_all_category_statistics.R', 'package_all_category_statistics.py']:
    extras.append((root / 'workflow/figure_delivery' / name, Path('复现脚本') / name))
for source, relative in extras:
    collect(source, relative, '统计与复现')
    collect(source, relative, '统计与复现', slim)
with (dest / '图片索引与SHA256.csv').open('w', encoding='utf-8-sig', newline='') as f:
    w = csv.writer(f, lineterminator='\n')
    w.writerow(['包内路径', '状态', '来源', 'SHA256']); w.writerows(manifest)
shutil.copy2(dest / '图片索引与SHA256.csv', out / 'delivery_manifest.csv')
for folder in [dest, slim]:
    archive = shutil.make_archive(str(folder), 'zip', root_dir=folder.parent, base_dir=folder.name)
    with zipfile.ZipFile(archive) as z:
        assert z.testzip() is None
        for p in folder.rglob('*'):
            if p.is_file():
                assert z.read(f'{folder.name}/{p.relative_to(folder)}') == p.read_bytes()
    print(folder)
