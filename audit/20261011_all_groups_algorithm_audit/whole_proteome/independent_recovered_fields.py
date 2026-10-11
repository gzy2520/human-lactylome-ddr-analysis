import csv,gzip,re,zipfile,io,json
from pathlib import Path
from collections import Counter,defaultdict
O=Path(__file__).parent; S=Path('/Users/gzy2520/Desktop/Research/kla/data')
def rows(p,sep=','):
 with (gzip.open(p,'rt') if str(p).endswith('.gz') else open(p)) as f:return list(csv.DictReader(f,delimiter=sep))
def save(n,v):
 with open(O/n,'w') as f:w=csv.DictWriter(f,fieldnames=list(v[0]));w.writeheader();w.writerows(v)
pat=re.compile(r'(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9](?:[A-Z0-9]{3}[0-9])?)')
def parse(v):
 a=set()
 for t in re.split('[;,:]',v):
  t=t.strip()
  if re.search('REV__|CON__|DECOY|CONTAM',t,re.I):continue
  if t.startswith(('sp|','tr|')):t=t[3:]
  t=t.split('|')[0].removeprefix('NX_');t=re.sub('-[0-9]+$','',t)
  if pat.fullmatch(t):a.add(t)
 return a
logs=[];bad=[]
def scan(label,rr,selected=None):
 rr=list(rr);cols=[j for j in rr[0] if re.search('protein|accession|razor|reverse|contaminant',j,re.I)]
 for j in cols:
  vals=[r[j] for r in rr];pref=[(i,r) for i,r in enumerate(rr) if re.search('REV__|CON__|DECOY|CONTAM',r[j],re.I)];logs.append(dict(Source=label,Field=j,Rows=len(rr),Colon=sum(':'in v for v in vals),Pipe=sum('|'in v for v in vals),PrefixRows=len(pref),TenCharIds=len({a for v in vals for a in parse(v) if len(a)==10})))
  for i,r in pref:
   if not any(r.get(k,'').lower() in ['+','1','true','yes'] for k in ['Reverse','Potential contaminant','Contaminant','Only identified by site']):bad.append(dict(Source=label,Row=i+1,Field=j,Value=r[j]))
 return rr
u=rows(O/'recovered_unit_values.csv');seen=set()
for x in u:
 p=Path(x['SourceFile'])
 if p in seen:continue
 seen.add(p)
 if 'PXD010154' in str(p):scan(str(p),rows(p,'\t'))
 elif 'PXD022005' in str(p):
  with zipfile.ZipFile(p) as z:scan(str(p)+'!proteinGroups.txt',csv.DictReader(io.TextIOWrapper(z.open('proteinGroups.txt')),delimiter='\t'))
 elif 'PXD002400' in str(p):
  with zipfile.ZipFile(p) as z:scan(str(p)+'!evidence.txt MCF10A rows',(r for r in csv.DictReader(io.TextIOWrapper(z.open('evidence.txt')),delimiter='\t') if '_10A_' in r['Raw file']))
 elif 'PXD030304' in str(p):
  with open(p) as f:h=next(csv.reader(f,delimiter='\t'))[1:]
  vals=[s.split(';')[0] for s in h];logs.append(dict(Source=str(p),Field='protein-column header accession before ;symbol',Rows=len(h),Colon=sum(':'in s for s in vals),Pipe=sum('|'in s for s in vals),PrefixRows=sum(bool(re.search('REV__|CON__',s)) for s in vals),TenCharIds=len({a for v in vals for a in parse(v) if len(a)==10})))
  expected={a for v in vals for a in parse(v)};actual={x['BaseAccession'] for x in rows(O/'recovered_protein_membership.csv.gz') if (x['PXD'],x['SampleGroup']) in [(y['PXD'],y['SampleGroup']) for y in u if y['ReferencePXD']=='PXD030304']};assert actual<=expected
save('recovered_all_identifier_fields.csv',logs)
if bad:save('recovered_prefix_without_flags.csv',bad)
reg=rows(O/'run_to_unit_registry.csv');mp={x['Automatic_MS_filename']:x for x in rows(S/'PXD030304/search_results/ProCan-DepMapSanger_mapping_file_replicates.txt','\t')};checks=[]
for uid in sorted({x['SampleID'] for x in reg if x['SampleID'].startswith('PXD030304_')}):
 rr=[x for x in reg if x['SampleID']==uid];lines={mp[r['RawFile']]['Cell_line'] for r in rr};toks={r['PreparationID'] for r in rr};checks.append(dict(SampleID=uid,RegistryRows=len(rr),UniqueRuns=len({r['RawFile']for r in rr}),CellLines=';'.join(sorted(lines)),PreparationTokens=';'.join(sorted(toks)),AllMapped=all(r['RawFile']in mp for r in rr)))
save('independent_procan_registry_checks.csv',checks)
p=S/'PXD030304/metadata/project.json';data=json.load(open(p));texts=[]
def rec(v,k=''):
 if isinstance(v,dict):
  for kk,vv in v.items():rec(vv,kk)
 elif isinstance(v,list):
  for vv in v:rec(vv,k)
 elif isinstance(v,str) and re.search('triplic|duplic|mass spectrom',v,re.I):texts.append(k+': '+v)
rec(data);(O/'procan_metadata_evidence.txt').write_text('\n'.join(texts))
print('fields',len(logs),'prefix_without_flags',len(bad),'registry units',len(checks))
