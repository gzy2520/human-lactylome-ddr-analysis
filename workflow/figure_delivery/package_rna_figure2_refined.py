"""Package the two confirmed RNA panels separately, without a montage."""
from pathlib import Path
import csv
import hashlib
import shutil

root = Path(__file__).resolve().parents[2]
out = root / 'outputs/20261007_rna_figure2_refined'
dest = Path('/Users/gzy2520/Desktop/图二_RNA单图美化_20261007')
dest.mkdir(exist_ok=True)
files = [(p, Path(p.name)) for p in sorted((out / 'figures').iterdir()) if p.suffix in ['.png', '.pdf', '.svg']]
files += [(p, Path('数据与统计') / p.name) for p in sorted((out / 'tables').glob('*.csv')) if 'displayed_brackets' not in p.name]
files += [(root / 'audit/20261007_rna_figure2_refined/README.md', Path('图注与说明.md')),
          (root / 'workflow/figure_delivery/render_rna_figure2_refined.R', Path('重绘脚本/render_rna_figure2_refined.R'))]
records = []
for source, relative in files:
    target = dest / relative
    target.parent.mkdir(exist_ok=True, parents=True)
    shutil.copy2(source, target)
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    assert hashlib.sha256(target.read_bytes()).hexdigest() == digest
    records.append([str(relative), str(source), digest])
for path in [dest / '文件来源与SHA256.csv', out / 'delivery_manifest.csv']:
    with path.open('w', encoding='utf-8-sig', newline='') as handle:
        writer = csv.writer(handle, lineterminator='\n')
        writer.writerow(['包内路径', '原始来源', 'SHA256'])
        writer.writerows(records)
shutil.make_archive(str(dest), 'zip', root_dir=dest.parent, base_dir=dest.name)
print(dest)
