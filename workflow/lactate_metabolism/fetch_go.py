#!/usr/bin/env python3
"""Freeze complete human GO annotations; no expression analysis in Python."""
import urllib.request,urllib.parse,json,time,pathlib,hashlib,csv
root=pathlib.Path(__file__).resolve().parents[2]; out=root/'config/lactate_metabolism/raw';out.mkdir(parents=True,exist_ok=True)
terms=['GO:0019244','GO:0006089','GO:0006099','GO:0006094','GO:0006096','GO:0015727']
base='https://www.ebi.ac.uk/QuickGO/services/'
manifest=[]
def fetch(name,url):
 p=out/name
 if not p.exists():
  for attempt in range(5):
   try:
    req=urllib.request.Request(url,headers={'Accept':'application/json'})
    with urllib.request.urlopen(req,timeout=90) as r: data=r.read()
    json.loads(data);p.write_bytes(data);break
   except Exception:
    if attempt==4:raise
    time.sleep(2*(attempt+1))
 data=p.read_bytes();manifest.append([name,url,hashlib.sha256(data).hexdigest()]);return json.loads(data)
fetch('go_release.json',base+'ontology/go/about')
fetch('annotation_release.json',base+'annotation/about')
for term in terms:
 key=term.replace(':','_')
 fetch(key+'_term.json',base+'ontology/go/terms/'+urllib.parse.quote(term))
 fetch(key+'_descendants.json',base+'ontology/go/terms/'+urllib.parse.quote(term)+'/descendants?relations=is_a,part_of')
 page=1
 while True:
  params={'goId':term,'taxonId':'9606','goUsage':'descendants','goUsageRelationships':'is_a,part_of','limit':200,'page':page}
  r=fetch(key+'_annotations_'+str(page)+'.json',base+'annotation/search?'+urllib.parse.urlencode(params))
  assert r.get('numberOfHits',0)>=0
  if page>=r['pageInfo']['total'] or not r['results']:break
  page+=1
 print(term,r['numberOfHits'],'annotations',page,'pages',flush=True)
with (out.parent/'download_manifest.csv').open('w') as f:
 w=csv.writer(f,lineterminator='\n');w.writerow(['File','URL','SHA256']);w.writerows(manifest)
