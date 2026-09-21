#!/usr/bin/env python3
import pathlib,csv,json,urllib.request,urllib.parse,hashlib,time
root=pathlib.Path(__file__).resolve().parents[2];cfg=root/'config/lactate_metabolism'
rows=list(csv.DictReader((cfg/'curated_id_mapping.csv').open()))
accs=sorted(set(r['UniProt'] for r in rows if r['UniProt']))
manifest=[]
for start in range(0,len(accs),30):
 query='organism_id:9606 AND reviewed:true AND ('+' OR '.join('accession:'+a for a in accs[start:start+30])+')'
 url='https://rest.uniprot.org/uniprotkb/search?'+urllib.parse.urlencode({'query':query,'format':'json','size':500})
 p=cfg/'raw'/('uniprot_'+hashlib.sha256(query.encode()).hexdigest()[:12]+'.json')
 if not p.exists():
  for attempt in range(4):
   try:
    with urllib.request.urlopen(url,timeout=90) as r:
     data=r.read();headers=dict(r.headers)
    json.loads(data);p.write_bytes(data);p.with_suffix('.headers.json').write_text(json.dumps(headers,indent=2));break
   except Exception:
    if attempt==3:raise
    time.sleep(2)
 manifest.append([p.name,url,hashlib.sha256(p.read_bytes()).hexdigest()])
 print(p.name,len(json.loads(p.read_text())['results']),flush=True)
with (cfg/'protein_download_manifest.csv').open('w') as f:
 w=csv.writer(f,lineterminator='\n');w.writerow(['File','URL','SHA256']);w.writerows(manifest)
