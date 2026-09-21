#!/usr/bin/env python3
import pathlib,json,csv,datetime,hashlib
cfg=pathlib.Path('config/lactate_metabolism')
manifest=list(csv.DictReader((cfg/'protein_download_manifest.csv').open()))
entries={}
for row in manifest:
 for r in json.loads((cfg/'raw'/row['File']).read_text())['results']:entries[r['primaryAccession']]=r
rows=[];maps=[]
for acc,r in entries.items():
 gene=';'.join(g.get('geneName',{}).get('value','') for g in r.get('genes',[]))
 loc=[];functions=[];reactions=[];refs=[]
 for c in r.get('comments',[]):
  if c['commentType']=='FUNCTION':functions.extend(t['value'] for t in c.get('texts',[]))
  if c['commentType']=='SUBCELLULAR LOCATION':
   for l in c.get('subcellularLocations',[]):loc.append(l.get('location',{}).get('value',''))
  if c['commentType']=='CATALYTIC ACTIVITY':reactions.append(c.get('reaction',{}).get('name',''))
 for x in r.get('uniProtKBCrossReferences',[]):
  if x['database']=='Ensembl':
   for p in x.get('properties',[]):
    if p['key']=='GeneId':maps.append([acc,gene,p['value'].split('.')[0]])
 rows.append([acc,gene,'; '.join(sorted(set(loc))),' | '.join(functions),' | '.join(reactions),r['entryAudit']['lastAnnotationUpdateDate']])
with (cfg/'uniprot_function_evidence.csv').open('w') as f:
 w=csv.writer(f,lineterminator='\n');w.writerow(['UniProt','GeneName','Compartment','Function','CatalyticReaction','LastAnnotationUpdate']);w.writerows(rows)
with (cfg/'uniprot_current_ensembl.csv').open('w') as f:
 w=csv.writer(f,lineterminator='\n');w.writerow(['UniProt','GeneName','Ensembl']);w.writerows(sorted(set(tuple(r) for r in maps)))
print(len(rows),'reviewed protein records frozen')
