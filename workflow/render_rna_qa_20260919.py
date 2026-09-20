from pathlib import Path
import subprocess, csv
from PIL import Image, ImageDraw
root=Path('outputs/20260919_repair_release')
qa=Path('audit/20260919_rna_repair/visual')
files=sorted(root.glob('*/*.pdf'))
assert len(files)==18,len(files)
rows=[]
for i,p in enumerate(files):
 dest=qa/f'{i+1:02d}'
 subprocess.run(['pdftoppm','-scale-to','1600','-png','-singlefile',str(p),str(dest)],check=True)
 rows.append({'Index':i+1,'File':str(p),'Preview':str(dest)+'.png'})
with (qa/'index.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
for batch in range(3):
 sheet=Image.new('RGB',(1800,1320),'white');draw=ImageDraw.Draw(sheet)
 for j,r in enumerate(rows[batch*6:(batch+1)*6]):
  im=Image.open(r['Preview']);im.thumbnail((595,620))
  x=(j%3)*600;y=(j//3)*660
  sheet.paste(im,(x+(600-im.width)//2,y+25))
  draw.text((x+10,y+5),f"{r['Index']:02d} {Path(r['File']).stem[:65]}",fill='black')
 sheet.save(qa/f'contact_{batch+1}.jpg',quality=95)
print('Rendered',len(rows),'PDFs')
