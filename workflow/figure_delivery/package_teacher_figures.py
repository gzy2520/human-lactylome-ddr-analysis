"""Reuse the audited version selection and substitute reviewed no-grid renders."""
from pathlib import Path
import csv, hashlib, shutil, json

root=Path(__file__).resolve().parents[2]
out=root/'outputs/20261006_teacher_figures'
dest=Path('/Users/gzy2520/Desktop/乳酸化全部图片_无网格与统计标注_20261006')
dest.mkdir(exist_ok=True)
previous=root/'outputs/20261006_figures_latest/selected_figure_manifest.csv'
records=list(csv.DictReader(previous.open(encoding='utf-8-sig')))
mapping={
 'Figure_1a_DDR_fraction_boxplot':'protein/fraction/Figure_1_DDR_fraction_candidate_category_boxplot_refined',
 'Figure_1b_MKI67_over_H3C1_boxplot':'protein/ratios/Figure_1_MKI67_over_H3C1_boxplot',
 'Figure_1_companion_MKI67_over_ACTB_boxplot':'protein/ratios/Figure_1_MKI67_over_ACTB_boxplot',
 'Figure_1_companion_MKI67_over_TUBB_boxplot':'protein/ratios/Figure_1_MKI67_over_TUBB_boxplot',
 'Figure_S1a_DDR_fraction_by_PXD':'protein/s1/Figure_S1a_DDR_fraction_by_PXD',
 'Figure_S1b_MKI67_over_H3C1_by_PXD':'protein/s1/Figure_S1b_MKI67_over_H3C1_by_PXD',
 'Figure_S1b_MKI67_over_H3C1_11_detected_only':'protein/s1/Figure_S1b_MKI67_over_H3C1_11_detected_only',
 'Figure_S2a_whole_proteome_UpSet':'protein/publication/Supplementary_Figure_S1a_whole_proteome_UpSet',
 'Figure_S2b_Kla_proteome_UpSet':'protein/publication/Supplementary_Figure_S1b_Kla_proteome_UpSet',
 'Supplementary_Figure_S2b_DDR_pathway_summary_cancer_cell_lines':'protein/publication/Supplementary_Figure_S2b_DDR_pathway_summary_cancer_cell_lines',
 'Supplementary_Figure_S2c_DDR_pathway_summary_normal_cell_lines':'protein/publication/Supplementary_Figure_S2c_DDR_pathway_summary_normal_cell_lines',
}
for letter,pa in zip('cdefghi',['BER','NER','MMR','FA','HR','AEJ','NHEJ']):
 mapping[f'Figure_2{letter}_DDR_pathway_summary_{pa}_barplot']=f'protein/pathways/Figure_2_DDR_pathway_summary_{pa}_barplot'
for stem in ['Figure_S3a_DDR_pathway_matrix_tumor_tissues','Figure_S3b_DDR_pathway_matrix_non_tumor_tissues',
             'Figure_S4a_DDR_pathway_matrix_cancer_cell_lines','Figure_S4b_DDR_pathway_matrix_normal_cell_lines']:
 folder='Supplementary_Figure_S3' if stem.startswith('Figure_S3') else 'Supplementary_Figure_S4'
 mapping[stem]=f'protein/matrix_root/results/final_figures_and_tables/{folder}/{stem}'
manifest=[]
def collect(source,relative,status,old=''):
 assert source.is_file(),source
 target=dest/relative;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(source,target)
 manifest.append([str(relative),status,str(source),hashlib.sha256(source.read_bytes()).hexdigest(),old])
for r in records:
 old=Path(r['原始来源']);relative=Path(r['包内路径'])
 if '20261006_figures_latest/ML/' in str(old):
  new=out/'ML'/old.relative_to(root/'outputs/20261006_figures_latest/ML')
 elif '20260919_repair_release' in str(old):
  new=out/'RNA'/old.relative_to(root/'outputs/20260919_repair_release')
 elif '20260920_lactate_metabolism/figures' in str(old):
  new=out/'lactate/figures'/old.name
 else:
  new=(out/mapping[old.stem]).with_suffix(old.suffix) if old.stem in mapping else out/'protein/publication'/old.name
 if relative.parts[0]=='01_当前蛋白组主图与附图':relative=Path('01_蛋白组单图与附图',*relative.parts[1:])
 collect(new,relative,r['版本状态'],str(old))
for source in sorted((out/'main_figures').iterdir()):
 if source.suffix in ['.png','.pdf','.docx','.csv']:collect(source,Path('00_正文主图与ML结果表')/source.name,'当前正文')
for source in sorted((out/'statistics').iterdir()):
 if source.suffix=='.csv':collect(source,Path('统计结果')/source.name,'当前统计')
for name in ['numerical_validation.csv','grid_audit.csv','final_export_validation.csv']:
 source=out/'qa'/name
 collect(source,Path('核对结果')/name,'交付核对')
for name in ['REPORT.md','figure_captions.md']:
 source=root/'audit/20261006_teacher_figures'/name
 collect(source,Path('说明与图注')/name,'当前说明')
for source in sorted((root/'workflow/figure_delivery').iterdir()):
 if source.name in ['render_teacher_nogrid.R','teacher_pairwise.R','assemble_teacher_main.R',
                    'render_latest_ml_zero.R','package_teacher_figures.py','build_ml_body_table.py',
                    'qa_teacher_release.R','README.md']:
  collect(source,Path('重绘脚本')/source.name,'复现脚本')
with (dest/'图片索引与来源.csv').open('w',encoding='utf-8-sig',newline='') as f:
 w=csv.writer(f,lineterminator='\n');w.writerow(['包内路径','版本状态','本次来源','SHA256','上版权威来源']);w.writerows(manifest)
with (out/'delivery_manifest.csv').open('w',encoding='utf-8-sig',newline='') as f:
 w=csv.writer(f,lineterminator='\n');w.writerow(['包内路径','版本状态','本次来源','SHA256','上版权威来源']);w.writerows(manifest)
summary={'Destination':str(dest),'PNG':sum(Path(r[0]).suffix=='.png' for r in manifest),
         'PDF':sum(Path(r[0]).suffix=='.pdf' for r in manifest),'DOCX':sum(Path(r[0]).suffix=='.docx' for r in manifest)}
(out/'delivery_summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2))
shutil.make_archive(str(dest),'zip',root_dir=dest.parent,base_dir=dest.name)
print(json.dumps(summary,ensure_ascii=False))
