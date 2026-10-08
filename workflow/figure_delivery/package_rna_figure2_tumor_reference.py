"""Package separate RNA panels with tumor-reference comparison annotations."""
from pathlib import Path
import csv
import hashlib
import shutil
import zipfile

root = Path(__file__).resolve().parents[2]
out = root / 'outputs/20261008_rna_figure2_tumor_reference'
dest = Path('/Users/gzy2520/Desktop/图二_RNA_癌组织对照统计_20261008')
dest.mkdir(exist_ok=True)
files = [(p, Path(p.name)) for p in sorted((out / 'figures').iterdir())
         if p.suffix in ['.png', '.pdf', '.svg']]
files += [(p, Path('数据与统计') / p.name) for p in sorted((out / 'tables').glob('*.csv'))]
files += [(root / 'audit/20261008_rna_figure2_tumor_reference/README.md', Path('图注与说明.md')),
          (out / 'sessionInfo.txt', Path('软件版本.txt')),
          (root / 'workflow/figure_delivery/render_rna_figure2_tumor_reference.R',
           Path('重绘脚本/render_rna_figure2_tumor_reference.R'))]
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
archive = Path(shutil.make_archive(str(dest), 'zip', root_dir=dest.parent, base_dir=dest.name))
with zipfile.ZipFile(archive) as handle:
    assert handle.testzip() is None
    for relative, _, digest in records:
        assert hashlib.sha256(handle.read(f'{dest.name}/{relative}')).hexdigest() == digest
print(dest)
