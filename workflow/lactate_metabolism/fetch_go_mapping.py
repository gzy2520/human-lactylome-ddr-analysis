#!/usr/bin/env python3
import pathlib,csv,json,urllib.request,urllib.parse,time,hashlib
cfg=pathlib.Path('config/lactate_metabolism');rows=list(csv.DictReader((cfg/'go_annotations.csv').open()));accs=sorted(set(r['UniProt'] for r in rows if ':' not in r['UniProt']));out=[];manifest=[]
for i in range(0,len(accs),60):
 q='organism_id:9606 AND ('+' OR '.join('accession:'+a for a in accs[i:i+60])+')'
 url='https://rest.uniprot.org/uniprotkb/search?'+urllib.parse.urlencode({'query':q,'format':'json','size':500,'fields':'accession,gene_primary,xref_ensembl,xref_geneid'})
 p=cfg/'raw'/('go_uniprot_'+hashlib.sha256(q.encode()).hexdigest()[:12]+'.json')
 if not p.exists():
  for attempt in range(4):
   try:
    data=urllib.request.urlopen(url,timeout=90).read();json.loads(data);p.write_bytes(data);break
   except Exception:
    if attempt==3:raise
    time.sleep(2)
 r=json.loads(p.read_text())
 for x in r['results']:
  acc=x['primaryAccession'];sym=';'.join(g.get('geneName',{}).get('value','') for g in x.get('genes',[]))
  ent=';'.join(y['id'] for y in x.get('uniProtKBCrossReferences',[]) if y['database']=='GeneID')
  for y in x.get('uniProtKBCrossReferences',[]):
   if y['database']=='Ensembl':
    for z in y.get('properties',[]):
     if z['key']=='GeneId':out.append((acc,z['value'].split('.')[0],ent,sym))
 manifest.append([p.name,url,hashlib.sha256(p.read_bytes()).hexdigest()]);print(i,len(r['results']),flush=True)
with (cfg/'go_current_uniprot_mapping.csv').open('w') as f:
 w=csv.writer(f,lineterminator='\n');w.writerow(['UniProt','Ensembl','Entrez','Symbol']);w.writerows(sorted(set(out)))
with (cfg/'go_mapping_download_manifest.csv').open('w') as f:
 w=csv.writer(f,lineterminator='\n');w.writerow(['File','URL','SHA256']);w.writerows(manifest)
