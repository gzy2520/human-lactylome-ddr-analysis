from pathlib import Path
import subprocess,csv,math,sys
from PIL import Image,ImageDraw
root=Path(sys.argv[1] if len(sys.argv)>1 else 'outputs/20260920_lactate_metabolism');qa=Path(sys.argv[2] if len(sys.argv)>2 else 'audit/20260920_lactate_metabolism/visual');qa.mkdir(parents=True,exist_ok=True)
rows=[]
for i,p in enumerate(sorted((root/'figures').glob('*.pdf')),1):
 dest=qa/f'{i:02d}'
 subprocess.run(['pdftoppm','-scale-to','1500','-png','-singlefile',str(p),str(dest)],check=True)
 rows.append({'Index':i,'File':str(p),'Preview':str(dest)+'.png'})
with (qa/'index.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=rows[0],lineterminator='\n');w.writeheader();w.writerows(rows)
for b in range(math.ceil(len(rows)/6)):
 sheet=Image.new('RGB',(1800,1320),'white');draw=ImageDraw.Draw(sheet)
 for j,r in enumerate(rows[b*6:(b+1)*6]):
  im=Image.open(r['Preview']);im.thumbnail((595,620));x=j%3*600;y=j//3*660
  sheet.paste(im,(x+(600-im.width)//2,y+25));draw.text((x+5,y+5),f"{r['Index']:02d} {Path(r['File']).stem}",fill='black')
 sheet.save(qa/f'contact_{b+1}.jpg',quality=95)
print(len(rows),'PDFs rendered')
