#!/usr/bin/env python3
"""Collect publication figures and rerendered ML figures without modifying frozen outputs."""
from pathlib import Path
import shutil,csv,hashlib,json
root=Path(__file__).resolve().parents[2]
out=root/'outputs/20261004_all_figures'
dest=Path('/Users/gzy2520/Desktop/乳酸化全部图片_20261004')
dest.mkdir(exist_ok=True)
manifest=[]
def collect(source,label,exclude=()):
 for p in sorted(source.rglob('*')):
  if p.suffix.lower() not in ('.png','.pdf','.svg'):continue
  if any(x in p.relative_to(source).parts for x in exclude):continue
  rel=Path(label)/p.relative_to(source);target=dest/rel;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,target)
  manifest.append([str(rel),str(p),hashlib.sha256(p.read_bytes()).hexdigest()])
formal=root/'corrected_final_result_20260905_sample_only'
for folder in ['Figure_1','Figure_2','Figure_3','Supplementary_Figure_S1','Supplementary_Figure_S2','Supplementary_Figure_S3','Supplementary_Figure_S4']:
 collect(formal/folder,'01_正式蛋白组主图与附图/'+folder)
collect(root/'outputs/20260919_repair_release/core','02_RNA正式图')
collect(root/'outputs/20260919_repair_release/reference','03_RNA补充参考图')
collect(root/'outputs/20260920_lactate_metabolism/figures','04_乳酸代谢图')
collect(out/'ML','05_机器学习_零起点')
with (dest/'图片索引与来源.csv').open('w',encoding='utf-8-sig',newline='') as f:
 w=csv.writer(f);w.writerow(['包内路径','原始来源','SHA256']);w.writerows(manifest)
text='''# 乳酸化项目图片汇总

## 内容

01：正式蛋白组 Figure 1–3、Supplementary Figure S1–S4，包括既有配套图和无框版本。
02：校正后的正式 RNA 图（28 个独立材料 / 1,898 条 RNA）。
03：RNA 补充参考图，部分以历史 31 个材料组展示，勿与 28 个独立 RNA 材料混淆。
04：乳酸代谢、GO 与机制相关图。
05：三版主要机器学习结果：当前四个组学内比例目标、上一版 DDR/DNA repair 跨组学比例、此前全蛋白 Kla 占比模型。保留各自目标定义和评价范围。

PNG 适合直接查看和插入报告，PDF 适合放大、排版和矢量编辑。图片索引记录每个文件的来源与校验值。未收录试跑、预览、失败下载及已知错误转换的外部图片，避免与有效结果混淆。

## 机器学习坐标

占比、预测值、特征重要性、绝对误差及误差记录数等非负数值轴均从 0 开始，取消下端留白；材料、类别、基因等分类轴保留标签。

UMAP 两轴各减去该轴最小值，坐标从 0 开始，图内明确标为平移坐标。该平移不改变点间距离；UMAP 的全数据拟合着色只作描述性展示。

主包的误差图使用绝对误差，因此可从 0 开始。signed_diagnostic_appendix 中保留带正负号的残差，以便识别高估或低估；此类坐标不能裁掉负值。该附录是唯一保留负数范围的机器学习图。

本次只重绘图片，未重新训练或修改预测值。原始模型和旧图片保留在项目目录。外部当前四目标及上一版两目标使用已校正为 log2(TPM+0.5) 的 RNA，未采用旧的错误转换结果。
'''
(dest/'图片说明.md').write_text(text)
shutil.copytree(root/'workflow/figure_delivery',dest/'重绘脚本',dirs_exist_ok=True)
counts={}
for rel,_,_ in manifest:
 key=Path(rel).parts[0];counts[key]=counts.get(key,0)+1
(out/'delivery_summary.json').write_text(json.dumps({'destination':str(dest),'files':len(manifest),'counts':counts},ensure_ascii=False,indent=2))
shutil.make_archive(str(dest),'zip',root_dir=dest.parent,base_dir=dest.name)
print(json.dumps({'files':len(manifest),'counts':counts,'zip':str(dest)+'.zip'},ensure_ascii=False))
