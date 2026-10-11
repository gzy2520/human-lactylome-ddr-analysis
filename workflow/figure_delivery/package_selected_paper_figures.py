#!/usr/bin/env python3
"""Copy selected source figures into the project delivery folder; build vector previews."""
from pathlib import Path
import csv,hashlib,io,shutil,re
from reportlab.pdfgen import canvas
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.cidfonts import UnicodeCIDFont
from pypdf import PdfReader,PdfWriter,Transformation
ROOT=Path(__file__).resolve().parents[2]
PLAN=ROOT/'deliverables/20261011_paper_figure_plan'
OUT=Path('/Users/gzy2520/Desktop/Research/kla/results/交付文件/乳酸化论文选图方案_20261011')
rows=list(csv.DictReader((PLAN/'面板选图与更新清单.csv').open(encoding='utf-8-sig')))
pdfmetrics.registerFont(UnicodeCIDFont('STSong-Light'))
manifest=[]; pending=[]
for row in rows:
 panel=row['面板'];sources=[Path(x) for x in row['现有图或样式来源'].split(' | ') if x.endswith('.pdf')]
 if panel in ['Fig1A','Fig2C','S1A']:
  pending.append((panel,row['内容'],row['保留与处理']));continue
 if panel=='Fig3B':continue # A and B are already paired in the chosen joint figure.
 ident='Fig3AB' if panel=='Fig3A' else panel
 group=re.match(r'(Fig\d+|S\d+)',panel).group(1)
 dest=OUT/('01_主图' if panel.startswith('Fig') else '02_补充图')/group
 dest.mkdir(parents=True,exist_ok=True)
 assert len(sources)==1,(panel,sources)
 src=sources[0]
 status=('数据可沿用，统一排版' if panel=='Fig3A' else
         '数据可沿用，待拆分RNA面板' if panel.startswith('S5') else
         '数据可沿用，核对统计说明' if panel=='S3D' else
         '保留此图，更新数据或统计')
 label=re.sub(r'[/\\:*?"<>|]','_',row['内容'])
 target=dest/f'{ident}_{label}.pdf'
 shutil.copy2(src,target)
 png=src.with_suffix('.png')
 if png.exists():shutil.copy2(png,target.with_suffix('.png'))
 assert hashlib.sha256(src.read_bytes()).digest()==hashlib.sha256(target.read_bytes()).digest()
 manifest.append({'Panel':ident,'Title':row['内容'],'Status':status,'Source':str(src),'PDF':str(target),'PNG':str(target.with_suffix('.png')) if png.exists() else '', 'SHA256':hashlib.sha256(target.read_bytes()).hexdigest(),'Pages':len(PdfReader(target).pages)})
with (OUT/'图片清单.csv').open('w',encoding='utf-8-sig',newline='') as f:
 w=csv.DictWriter(f,fieldnames=manifest[0],lineterminator='\n');w.writeheader();w.writerows(manifest)
# The preview places original vector PDF pages below a separate, clearly labeled header.
def preview(items,path,title):
 writer=PdfWriter();num=0
 for item in items:
  reader=PdfReader(item['PDF'])
  for i,srcpage in enumerate(reader.pages):
   num+=1; W,H=1200,860; buf=io.BytesIO();c=canvas.Canvas(buf,pagesize=(W,H));c.setFont('STSong-Light',18)
   c.drawString(28,H-32,f"{item['Panel']}  {item['Title']}")
   c.setFont('STSong-Light',12);c.drawString(28,H-55,item['Status']+'；本页为选定源图预览，未在此重算数据。')
   c.setFont('STSong-Light',10);c.drawString(28,16,f"{title} | {num} | 源图第{i+1}/{len(reader.pages)}页")
   c.save();page=PdfReader(buf).pages[0]
   sw,sh=float(srcpage.mediabox.width),float(srcpage.mediabox.height)
   scale=min((W-48)/sw,(H-100)/sh);x=(W-sw*scale)/2;y=36+(H-100-sh*scale)/2
   page.merge_transformed_page(srcpage,Transformation().scale(scale).translate(x,y))
   writer.add_page(page)
 writer.write(path)
 return num
nmain=preview([x for x in manifest if x['Panel'].startswith('Fig')],OUT/'00_主图预览.pdf','论文主图选择')
nall=preview(manifest,OUT/'00_全部选定图预览.pdf','主图与补充图选择')
text='# 先看图片\n\n本文件夹已放入选定的现有PDF与PNG，按主图、补充图编号。\n\n'
text+='- [先看主图](<'+str(OUT/'00_主图预览.pdf')+'>)\n- [连续浏览全部选定图](<'+str(OUT/'00_全部选定图预览.pdf')+'>)\n\n'
text+='保留的旧图已经复制进来；预览页顶部标明后续更新状态。源图内容未修改。Fig3AB是10月11日全组审计后的联合图；其他需更新的图可先确认用途和版式。\n\n'
text+='## 尚无成图的面板\n\n'
for p,t,s in pending:text+=f'- {p} {t}：{s}。未用旧图冒充新面板。\n'
text+='\nFig2C的七通路原样本图已经放在S8A–G；新类别并集面板还需制作。S5暂放选中的完整源图，后续提取RNA重新排版。\n\n## 逐图打开\n\n'
for x in manifest:text+=f"- [{x['Panel']} {x['Title']}](<{x['PDF']}>) — {x['Status']}\n"
(OUT/'先看图片.md').write_text(text)
# Track the lightweight inventory and QA record; large exported images stay under results.
shutil.copy2(OUT/'图片清单.csv',PLAN/'实际图片交付清单.csv')
(PLAN/'图片交付核验.md').write_text(f'# 实际图片交付\n\n位置：{OUT}\n\n已复制{len(manifest)}份PDF、{sum(bool(x["PNG"]) for x in manifest)}份PNG。PDF逐文件与源SHA256一致。主图预览{nmain}页，全部预览{nall}页。\n\n待制作：'+', '.join(x[0] for x in pending)+'。Fig3A/B合并保存为Fig3AB；没有重复复制相同DDR面板。\n')
print(f'Copied {len(manifest)} PDFs and {sum(bool(x["PNG"]) for x in manifest)} PNGs; preview pages {nmain}/{nall}')
