#!/usr/bin/env python3
from pathlib import Path
import csv,hashlib,shutil,datetime,sys,subprocess
out=Path(sys.argv[1]);cfg=Path('config/lactate_metabolism')
def digest(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def writecsv(p,fields,rows):
 with p.open('w') as f:
  w=csv.writer(f,lineterminator='\n');w.writerow(fields);w.writerows(rows)
(out/'evidence').mkdir(exist_ok=True)
for p in cfg.iterdir():
 if p.is_file() and p.suffix in ('.csv','.tsv'):
  shutil.copy2(p,out/'evidence'/p.name)
shutil.copy2('docs/LACTATE_METABOLISM_METHODS_20260920.md',out/'METHODS.md')
figures=list(csv.DictReader((out/'tables/figure_inventory.csv').open()))
sources=[]
for r in figures:
 n=r['Figure']
 if n.startswith('GO_') or n.startswith('Mechanism_'):src='tables/material_gene_expression.csv; tables/membership_coverage.csv'
 elif n.startswith('Focus_samples') or n.startswith('LDHA_'):src='tables/sample_gene_expression.csv.gz; tables/materials_28.csv'
 elif n=='Focus_contrasts':src='tables/focus_expression_differences.csv'
 elif n=='Kla_detection_context':src='tables/focus_kla_detection_context.csv'
 else:src='evidence/literature_evidence.tsv; evidence/curated_mechanism_evidence.csv; workflow/lactate_metabolism/analyze.R'
 sources.append([n,src,'workflow/lactate_metabolism/analyze.R'])
writecsv(out/'tables/figure_sources.csv',['Figure','SourceTables','Renderer'],sources)
inputs=set(p for p in cfg.rglob('*') if p.is_file())
inputs.update(p for p in Path('workflow/lactate_metabolism').iterdir() if p.is_file())
inputs.update(Path(x) for x in ['config/rna_current_release.csv','data/publication_input/group_summary_31.csv','audit/20260916_full_31_rna_status/rna_31_group_status.csv','docs/LACTATE_METABOLISM_METHODS_20260920.md'])
roles={r['Role']:Path(r['Path']) for r in csv.DictReader(open('config/rna_current_release.csv'))}
for n in ['qsmooth_A_collapsed_log2tpm.tsv.gz','qsmooth_B_full_log2tpm.tsv.gz','unified_groupmedian_log2tpm.tsv.gz']:inputs.add(roles['qsmooth']/'matrices'/n)
inputs.add(roles['qsmooth']/'group_expansion_31.csv');inputs.add(roles['expression']/'group_sample_qc.csv')
inputs.update((roles['expression']/'matrices').glob('*_log2tpm.tsv.gz'))
writecsv(out/'input_sha256.csv',['Path','SHA256'],[[str(p),digest(p)] for p in sorted(inputs)])
writecsv(out/'annotation_snapshot_files.csv',['Path','FetchedFileMtimeUTC','SHA256'],[[str(p),datetime.datetime.fromtimestamp(p.stat().st_mtime,datetime.timezone.utc).isoformat(),digest(p)] for p in sorted((cfg/'raw').iterdir()) if p.is_file()])
(out/'README.md').write_text('# 乳酸代谢探索交付\n\n开始：2026-09-20；完成：2026-09-21。\n\n先读 [中文汇报](REPORT_老师汇报.md)，方法见 [METHODS.md](METHODS.md)。\n\n23套图在 figures/；逐图源数据位置见 tables/figure_sources.csv；全部表在 tables/；冻结证据在 evidence/。源数据库响应在仓库 config/lactate_metabolism/raw/，由 input_sha256.csv 追溯。\n\n本包是表达探索，不代表乳酸浓度、代谢通量、Kla修饰强度或因果机制。显著性方法未调整；未训练预测模型。\n\n数值复现入口：bash workflow/lactate_metabolism/run.sh <新空输出目录>。\n')
files=sorted(p for p in out.rglob('*') if p.is_file() and p.name!='release_sha256.csv')
writecsv(out/'release_sha256.csv',['Path','SHA256'],[[str(p.relative_to(out)),digest(p)] for p in files])
print('FROZEN',len(inputs),'inputs;',len(files),'release files')
