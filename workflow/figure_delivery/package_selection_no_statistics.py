"""Package current individual figures without biological significance annotations."""
from pathlib import Path
import csv
import hashlib
import re
import shutil
import subprocess
import zipfile

root=Path(__file__).resolve().parents[2]
out=root/'outputs/20261008_figure_selection_no_statistics'
dest=Path('/Users/gzy2520/Desktop/乳酸化选图包_无统计标注_20261008')
dest.mkdir(exist_ok=True)
plan=list(csv.DictReader((root/'audit/20261008_figure_selection_no_statistics/render_plan.csv').open()))
mapping={r['Stem']:out/'figures'/r['Stem'] for r in plan}
paired=root/'outputs/20261007_rna_figure2_refined/figures'
mapping.update({
 'Figure_1a_RNA_DDR_and_lactylated_gene_expression_28materials_boxplot':paired/'Figure_2a_RNA_DDR_and_Kla_gene_expression',
 'Figure_1b_RNA_proliferation_hallmark_28materials_boxplot':paired/'Figure_2b_RNA_proliferation'})
rows=list(csv.DictReader((root/'outputs/20261006_teacher_figures/delivery_manifest.csv').open(encoding='utf-8-sig')))
manifest=[]
for row in rows:
    relative=Path(row['包内路径'])
    if relative.parts[0].startswith(('00_','90_')) or relative.suffix not in ['.png','.pdf']:
        continue
    source=mapping[relative.stem].with_suffix(relative.suffix) if relative.stem in mapping else Path(row['本次来源'])
    target=dest/relative
    target.parent.mkdir(parents=True,exist_ok=True)
    shutil.copy2(source,target)
    digest=hashlib.sha256(source.read_bytes()).hexdigest()
    assert hashlib.sha256(target.read_bytes()).hexdigest()==digest
    if relative.suffix=='.pdf' and not relative.parts[0].startswith('04_'):
        text=subprocess.check_output(['pdftotext','-layout',str(target),'-'],text=True)
        assert not re.search(r'ANOVA|\bp\s*[=<>]|\bq\s*[=<>]|\bns\b|\bNE\b|CR2|BH q',text),relative
    manifest.append([str(relative),str(source),digest])
assert sum(Path(r[0]).suffix=='.png' for r in manifest)==111
assert sum(Path(r[0]).suffix=='.pdf' for r in manifest)==111
for p in [dest/'图片索引与来源.csv',out/'delivery_manifest.csv']:
    with p.open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.writer(f,lineterminator='\n');w.writerow(['图片路径','来源','SHA256']);w.writerows(manifest)
(dest/'选图说明.txt').write_text(
    '用于确认正文与补充材料最终选图。\n'
    '共111张当前独立单图，每张提供PNG和PDF；没有组合图和历史模型版本。\n'
    '已去除生物学组间比较横线、星号、ns/NE、P/q值及总体检验副标题。\n'
    '原始数据点、箱体、均值、中位数、误差线和样本量保留；所有图保持无背景网格。\n'
    '机器学习保留MAE、R²等既有模型评价指标。\n'
    '两张RNA主图采用已确认的28份材料及最新无横线样式。\n'
    '请按图片文件名或截图圈选最终需要的单图，再统一确定统计展示。\n',encoding='utf-8')
archive=shutil.make_archive(str(dest),'zip',root_dir=dest.parent,base_dir=dest.name)
with zipfile.ZipFile(archive) as z:
    assert z.testzip() is None
    for rel,_,digest in manifest:
        assert hashlib.sha256(z.read(f'{dest.name}/{rel}')).hexdigest()==digest
print(dest)
