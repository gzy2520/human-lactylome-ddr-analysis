from pathlib import Path
import pandas as pd,numpy as np,re,json
root=Path(__file__).resolve().parents[3];out=Path(__file__).resolve().parent;raw=Path('/Users/gzy2520/Desktop/Research/kla/data')
release=root/'deliverables/20261011_lactylation_ddr_tall_verified/tables'
meta=pd.read_csv(release/'display_mapping_31.csv');meta=meta[meta.X>12]
cur=pd.read_csv(release/'ddr_source_unit_values.csv');cur=cur[(cur.Dataset=='Lactylome (Kla)')&(cur.X>12)]
go=pd.read_csv(root/'data/publication_input/human_ddr_go_annotations.tsv',sep='\t',dtype=str).fillna('')
pat=re.compile(r'^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9](?:[A-Z0-9]{3}[0-9])?)$')
def ids(v):
 ans=[]
 for t in re.split('[;,:]',str(v)):
  t=t.strip()
  if re.match(r'^(REV__|CON__|DECOY|CONTAM|CRAP)',t,re.I):continue
  if t.startswith(('sp|','tr|')):t=t.split('|')[1]
  else:t=t.split('|')[0]
  t=re.sub(r'-\d+$','',re.sub(r'^NX_','',t))
  if pat.fullmatch(t) and t not in ans:ans.append(t)
 return ans
keep=(go['GENE PRODUCT DB']=='UniProtKB')&go['TAXON ID'].isin(['9606',''])&~go.QUALIFIER.str.split('|').apply(lambda x:'NOT' in x)
ddr={a for v in go.loc[keep,'GENE PRODUCT ID'] for a in ids(v)}
def oldids(v):
 ans=[]
 for t in re.split('[;,]',str(v)):
  t=re.sub(r'^.*\|([^|]+)\|.*$',r'\1',t.strip());t=re.sub(r'^([^|;]+)\|.*$',r'\1',t);t=re.sub(r'-\d+$','',re.sub(r'^NX_','',t))
  if pat.fullmatch(t):ans.append(t)
 return sorted(set(ans))
def positive(s):return pd.to_numeric(s,errors='coerce').fillna(0)>0
def flagged(s):return s.astype(str).str.lower().isin(['+','true','1','yes','y'])
def clean(d):
 k=pd.Series(True,index=d.index)
 for c in ['Reverse','Potential contaminant','Contaminant','Only identified by site']:
  if c in d:k &= ~flagged(d[c])
 for c in ['Protein','Leading proteins']:
  if c in d:k &= ~d[c].str.match(r'^(?:REV__|CON__|DECOY|CONTAM|CRAP)',case=False)
 return k
rows=[];evidence=[];headers={};rules=[]
def add(d,k,px,g,s,col,path,rule,first=False,oldfirst=False):
 chosen=d.loc[k]
 acc=set();old=set()
 for ix,v in chosen[col].items():
  new=ids(v);prev=oldids(v)
  if first:new=new[:1]
  if first and px in ['PXD078736','PXD078013']:
   prev=ids(re.sub(r'(?:REV__|CON__)+','',str(v)))[:1]
  elif first:prev=prev[:1]
  acc.update(new);old.update(prev)
  evidence.append(dict(PXD=px,SampleGroup=g,SampleID=s,SourceFile=str(path),SourceRow=int(ix)+2,Column=col,RawAccession=v,Parsed=';'.join(new),OldParsed=';'.join(prev),Changed=set(new)!=set(prev)))
 for a in acc:rows.append(dict(PXD=px,SampleGroup=g,SampleID=s,BaseAccession=a,IsDdr=a in ddr,SourceFile=str(path),Rule=rule))
 rules.append(dict(PXD=px,SampleGroup=g,SampleID=s,SourceFile=str(path),RawRows=len(d),EligibleRows=len(chosen),ProteinColumn=col,Rule=rule,NewProteins=len(acc),OldIdentifierProteins=len(old),OnlyNew=';'.join(sorted(acc-old)),OnlyOld=';'.join(sorted(old-acc))))
 headers[str(path)]=list(d.columns)
def table(rel,sheet=None,skip=0):
 p=raw/rel
 if sheet is None:d=pd.read_csv(p,sep='\t' if p.suffix in ['.txt','.tsv']else ',',dtype=str,low_memory=False).fillna('')
 else:d=pd.read_excel(p,sheet_name=sheet,skiprows=skip,dtype=str).fillna('')
 return p,d
# MaxQuant site tables, retaining source-defined protein representatives where the original method does.
for px,rel,col,mode,groups in [
 ('PXD050147','PXD050147/search_results/Lactyl_K_Sites.txt','Proteins','loc',{'HepG2 WT and SIRT1 or SIRT3 KO':[f'{a}_Kla_rep{i}' for a in ['WT','SIRT1KO','SIRT3KO'] for i in [1,2,3]]}),
 ('PXD028737','PXD028737/search_results/extracted_pairing/combined/combined/txt/La(K)Sites.txt','Proteins','loc',{'HMC3':['H0','H24']}),
 ('PXD060185','PXD060185/search_results/RESULT/combined/txt/La (K)Sites.txt','Protein','ident',{'MCF7':['A'],'MDA-MB-468':['B'],'MCF10A':['C'],'T-47D':['D']}),
 ('PXD078736','PXD078736/search_results/txt/La(K)Sites.txt','Protein','ident',{'HK-2 control and mannitol':[f'{x}_{i}' for x in ['ctr','man'] for i in [1,2,3]]}),
 ('PXD058534','PXD058534/search_results/extracted_pairing/txt/La (K)Sites.txt','Proteins','single',{'pretreated HK-2':['HK2']})]:
 p,d=table(rel);base=clean(d)
 if 'id'in d:base &= d.id.ne('')
 if mode=='ident' and 'Amino acid'in d:base &= d['Amino acid'].eq('K')
 for g,ss in groups.items():
  for s in ss:
   c=('Localization prob '+s) if mode=='loc'else ('Identification type '+s) if mode=='ident'else 'Localization prob'
   k=base & (d[c].ne('') & d[c].ne('Not Found') if mode=='ident'else positive(d[c]))
   add(d,k,px,g,s.replace('_Kla',''),col,p,f'{c}: '+('detected text'if mode=='ident'else '>0')+'; filter contaminant/reverse flags; unique base UniProt',first=mode=='ident')
# Author intensity/localization workbooks.
p,d=table('PXD054919/supplementary/41419_2025_8113_MOESM2_ESM.xlsx',0,1)
for s in ['A549-1','A549-2','A549-3']:add(d,positive(d['Intensity '+s]),'PXD054919','A549',s,'Protein accession',p,'Intensity '+s+'>0')
p,d=table('PXD063266/search_results/LactylSites.xlsx','Lactyl (K)Sites')
base=clean(d)&d.id.ne('')
for s in ['1','2','3']:add(d,base&positive(d['Localization prob '+s]),'PXD063266','PC-3M',s,'Proteins',p,'Localization prob '+s+'>0; valid id; contaminant/reverse excluded')
p,d=table('PXD070007/search_results/SA206LPLaB1_Annotation.xlsx','Annotation_Combine')
for g,ss in [('neural stem cells',['ENSA','HMP1']),('glioblastoma stem cells',['G2907','G3028','G3264','GSC23','MES28','RKI'])]:
 for s in ss:add(d,positive(d['Localization probability'])&positive(d[s]),'PXD070007',g,s,'Protein accession',p,'Localization probability>0; named model intensity>0; model identity retained')
# Spectronaut site report and raw-file precursor matrix.
p,d=table('PXD066351/search_results/XB01472B1DPLa-MSstats_Input.csv')
base=d['PTM.ModificationTitle'].eq('Lac (K)')&d['PTM.SiteAA'].eq('K')&positive(d['PTM.SiteProbability'])&positive(d['PTM.Quantity'])
for s in ['NC116','R116']:add(d,base&d['R.Condition'].eq(s),'PXD066351','HCT116 control and Roseburia co-culture',s,'PTM.ProteinId',p,'K lactylation; site probability and quantity>0; one raw file per condition')
p,d=table('PXD073311/search_results/extracted_pairing/IPX0015307003_Database_search_result/Database_search_result/report.pr_matrix.tsv')
for s in [f'A{x}h_{i}'for x in [0,6]for i in [1,2,3]]:
 c=[x for x in d.columns if x.endswith('_'+s+'.raw')];assert len(c)==1
 add(d,d['Modified.Sequence'].str.contains('K(UniMod:2114)',regex=False)&positive(d[c[0]]),'PXD073311','HUVEC control and Pg infection',s,'Protein.Group',p,'K(UniMod:2114) precursor; per-raw positive quantity')
# RKO: independently expand positive-evidence site IDs and join validated protein-site-position pairs.
p,e=table('PXD078013/search_results/evidence.txt');pp,d=table('PXD078013/search_results/proteinGroups.txt')
e=e.loc[clean(e)&positive(e['La (K)'])&e['La (K) site IDs'].ne('')]
es={}
for _,v in e.iterrows():
 for t in v['La (K) site IDs'].split(';'):es.setdefault(t,set()).add(v.Experiment)
valid=clean(d);rko={s:set()for s in ['sg_1','sg_2','sg737_1','sg737_2']};mismatch=[]
for ix,v in d.loc[valid].iterrows():
 ts=v['La (K) site IDs'].split(';') if v['La (K) site IDs'] else [];ps=v['La (K) site positions'].split(';') if v['La (K) site positions'] else []
 if len(ts)!=len(ps):mismatch.append(dict(Row=ix+2,SiteIDs=v['La (K) site IDs'],Positions=v['La (K) site positions']));continue
 for t,pos in zip(ts,ps):
  if not re.fullmatch(r'\d+',pos):continue
  for s in es.get(t,set()):rko[s].add(ix)
for s,rr in rko.items():
 z=d.copy();z['Representative']=z['Majority protein IDs'].where(z['Majority protein IDs'].ne(''),z['Protein IDs'])
 add(z,z.index.isin(rr),'PXD078013','RKO WT and GSK3B KO',s,'Representative',pp,'Kla-positive evidence site IDs; matching site-position pairs; first source-order majority representative',first=True)
pd.DataFrame(mismatch).to_csv(out/'rko_site_position_mismatch.csv',index=False)
# PEAKS and subcellular studies: fresh original parser, real PSM/site/position evidence.
for px in ['PXD028488','PXD053474']:
 d=pd.read_csv(out/(px.lower()+'_source_sites.csv'),dtype=str).fillna('')
 if px=='PXD028488':d=d[d.PrimaryIncluded.eq('True')]
 for g,sub in d.groupby('CellOrTissueType'):
  if g=='TALL-104':sub=sub[sub.EnrichmentStatus.eq('Enriched')];sub['AuditSample']=sub.SampleName
  else:sub['AuditSample']=g+('_1'if px=='PXD028488'and g=='HCT116'else '_technical_merged')
  for s,x in sub.groupby('AuditSample'):
   acc=set(x.BaseAccession);acc={a for a in acc if pat.fullmatch(a)}
   for a in acc:rows.append(dict(PXD=px,SampleGroup=g,SampleID=s,BaseAccession=a,IsDdr=a in ddr,SourceFile=';'.join(sorted(set(x.SourceFile))),Rule='original source parser: site-position mapped, technical acquisitions merged'))
   rules.append(dict(PXD=px,SampleGroup=g,SampleID=s,SourceFile=';'.join(sorted(set(x.SourceFile))),RawRows=len(x),EligibleRows=len(x),ProteinColumn='BaseAccession (reparsed source)',Rule='PEAKS PSM+protein-peptide start mapping; S3 or DDA/DIA site agreement for PXD053474',NewProteins=len(acc),OldIdentifierProteins=len(acc),OnlyNew='',OnlyOld=''))
mem=pd.DataFrame(rows).drop_duplicates(['PXD','SampleGroup','SampleID','BaseAccession']);mem=mem.merge(meta[['PXD','SampleGroup','GroupID']],on=['PXD','SampleGroup'],validate='many_to_one');mem['Dataset']='Lactylome (Kla)'
mem.to_csv(out/'independent_membership.csv.gz',index=False)
pd.DataFrame(evidence).to_csv(out/'raw_row_identifier_evidence.csv.gz',index=False);pd.DataFrame(rules).to_csv(out/'raw_source_rules.csv',index=False)
with (out/'raw_headers.json').open('w')as f:json.dump(headers,f,indent=2)
st=mem.groupby(['GroupID','PXD','SampleGroup','SampleID']).agg(Proteins=('BaseAccession','nunique'),DDR=('IsDdr','sum')).reset_index();st['Percent']=100*st.DDR/st.Proteins
st=st.merge(cur[['GroupID','SampleID','ProteinCount','DdrProteinCount','Fraction','UnitType','Condition']],on=['GroupID','SampleID'],validate='one_to_one');st['DeltaProteins']=st.Proteins-st.ProteinCount;st['DeltaDDR']=st.DDR-st.DdrProteinCount;st['DeltaPercent']=st.Percent-st.Fraction
st.to_csv(out/'sample_count_comparison.csv',index=False)
# Compare exact independent sets with the actual source-recovery/base memberships that fed the release.
base=pd.read_csv(root/'audit/20261011_ddr_sample_recovery/base_rebuild/all_membership.csv.gz');base=base[base.Dataset.eq('Lactylome (Kla)')]
rec=pd.read_csv(root/'audit/20261011_ddr_sample_recovery/kla/recovered_kla_membership.csv.gz')
tall=pd.read_csv(root/'audit/20261011_tall_algorithm_history/arm_scope/arm_sample_membership.csv.gz');tall=tall[tall.EnrichmentStatus.eq('Enriched')].rename(columns={'SampleName':'SampleID'});tall['PXD']='PXD028488';tall['SampleGroup']='TALL-104'
current=pd.concat([base,rec,tall],ignore_index=True);sets=[]
for key,x in mem.groupby(['PXD','SampleGroup','SampleID']):
 old=current[(current.PXD==key[0])&(current.SampleGroup==key[1])&(current.SampleID.astype(str)==key[2])];a=set(x.BaseAccession);b=set(old.BaseAccession)
 sets.append(dict(PXD=key[0],SampleGroup=key[1],SampleID=key[2],OnlyIndependent=';'.join(sorted(a-b)),OnlyCurrent=';'.join(sorted(b-a)),ExactSet=a==b))
pd.DataFrame(sets).to_csv(out/'exact_set_comparison.csv',index=False)
summary=st.groupby(['GroupID','PXD','SampleGroup']).agg(Units=('SampleID','count'),MaxAbsProteinDelta=('DeltaProteins',lambda x:abs(x).max()),MaxAbsDdrDelta=('DeltaDDR',lambda x:abs(x).max()),MaxAbsPercentDelta=('DeltaPercent',lambda x:abs(x).max())).reset_index();summary['Dataset']='Lactylome (Kla)';summary['Status']=np.where((summary.MaxAbsProteinDelta==0)&(summary.MaxAbsDdrDelta==0),'source_reparse_counts_match','correction_required');summary.to_csv(out/'group_algorithm_audit.csv',index=False)
print(summary.to_string(index=False));print('membership',len(mem),'source evidence rows',len(evidence));print(st[st.DeltaProteins.ne(0)|st.DeltaDDR.ne(0)].to_string(index=False))
