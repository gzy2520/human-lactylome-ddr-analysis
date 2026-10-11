from pathlib import Path
import pandas as pd,re,json
out=Path(__file__).resolve().parent;raw=Path('/Users/gzy2520/Desktop/Research/kla/data')
headers=json.load(open(out/'raw_headers.json'));paths=set(headers)
for f in ['pxd028488_source_sites.csv','pxd053474_source_sites.csv']:
 d=pd.read_csv(out/f,dtype=str).fillna('')
 for p in d.SourceFile.unique():paths.add(str(raw.parent/p)if not p.startswith('/')else p)
# PSM/protein mapping supports are separately assessed.
for px in ['PXD028488','PXD053474']:
 for p in (raw/px/'search_results').rglob('protein-peptides.csv'):
  if not any(x in str(p)for x in ['BV2','RAW-']):paths.add(str(p))
rows=[];summary=[]
cols=['Protein','Proteins','Leading proteins','Majority protein IDs','Protein IDs','Protein accession','Accession','Protein Accession','Protein.Group','PTM.ProteinId','PG.ProteinAccessions']
for path in sorted(paths):
 p=Path(path)
 if not p.exists():continue
 if p.suffix=='.xlsx':
  if p.name=='LactylSites.xlsx':sheet='Lactyl (K)Sites';skip=0
  elif p.name=='SA206LPLaB1_Annotation.xlsx':sheet='Annotation_Combine';skip=0
  elif p.name=='41419_2025_8113_MOESM2_ESM.xlsx':sheet=0;skip=1
  else:continue
  d=pd.read_excel(p,sheet_name=sheet,skiprows=skip,dtype=str).fillna('')
 else:
  d=pd.read_csv(p,sep='\t'if p.suffix in ['.txt','.tsv']else ',',dtype=str,low_memory=False).fillna('')
 anyflag=pd.Series(False,index=d.index)
 for c in ['Reverse','Potential contaminant','Contaminant','Only identified by site']:
  if c in d:anyflag|=d[c].str.lower().isin(['+','1','true','yes'])
 n=0;unflag=0
 for c in cols:
  if c not in d:continue
  hit=d[c].str.contains(r'REV__|CON__|DECOY|CONTAM',regex=True)
  n+=hit.sum();unflag+=(hit&~anyflag).sum()
  for ix,v in d.loc[hit].iterrows():rows.append(dict(SourceFile=str(p),DataRow=ix+2,Column=c,RawAccession=v[c],AnyExclusionFlag=bool(anyflag[ix]),Reverse=v.get('Reverse',''),PotentialContaminant=v.get('Potential contaminant',''),Protein=v.get('Protein',''),Proteins=v.get('Proteins',''),SiteId=v.get('id',''),IdentificationCtr1=v.get('Identification type ctr_1',''),IdentificationCtr2=v.get('Identification type ctr_2','')))
 summary.append(dict(SourceFile=str(p),Rows=len(d),AccessionColumns=';'.join(c for c in cols if c in d),PrefixFieldRecords=n,UnflaggedPrefixFieldRecords=unflag))
pd.DataFrame(rows).to_csv(out/'all_raw_prefixed_identifiers.csv.gz',index=False)
pd.DataFrame(summary).to_csv(out/'raw_prefix_scan_coverage.csv',index=False)
u=pd.DataFrame(rows);u=u[~u.AnyExclusionFlag];u.to_csv(out/'unflagged_raw_prefixed_identifiers.csv',index=False)
print(u[['SourceFile','DataRow','Column','RawAccession']].to_string(index=False));print('files scanned',len(summary))
