from pathlib import Path
import sys,tempfile,csv
import pandas as pd
root=Path(__file__).resolve().parents[3]
out=Path(__file__).resolve().parent
raw=Path('/Users/gzy2520/Desktop/Research/kla/data')
sys.path.insert(0,str(root/'workflow/source_legacy'))
from python.utils.extractors import extract_pxd028488
from python.utils.common import base_accession
with tempfile.TemporaryDirectory(prefix='tall-audit-') as td:
 td=Path(td)
 for collection,folder in [('Enrichment-Search files','TALL-NALAC-Search files'),('Nonenrichment-Search files','TALL-Nonenrichment-Search files')]:
  dst=td/'PXD028488/search_results'/collection/folder
  dst.parent.mkdir(parents=True,exist_ok=True);dst.symlink_to(raw/'PXD028488/search_results'/collection/folder,target_is_directory=True)
 sites,logs,dirs=extract_pxd028488(td,root,'10.1038/s41592-022-01523-1')
 sites['SourceFile']=sites.SourceFile.str.replace(str(td),str(raw),regex=False)
go=pd.read_csv(root/'data/publication_input/human_ddr_go_annotations.tsv',sep='\t',dtype=str).fillna('')
ddr={base_accession(x) for x in go.loc[(go['GENE PRODUCT DB']=='UniProtKB')&go['TAXON ID'].isin(['9606',''])&~go.QUALIFIER.str.split('|').apply(lambda x:'NOT' in x),'GENE PRODUCT ID']}
mem=sites[['EnrichmentStatus','SampleName','BaseAccession','SourceFile']].drop_duplicates()
mem['IsDdr']=mem.BaseAccession.isin(ddr)
mem.to_csv(out/'arm_sample_membership.csv.gz',index=False)
rows=[]
for (arm,sample),d in mem.groupby(['EnrichmentStatus','SampleName']):
 acc=set(d.BaseAccession);rows.append(dict(Scope=arm,SampleID=sample,Proteins=len(acc),DDR=len(acc&ddr),Percent=100*len(acc&ddr)/len(acc)))
for arm in ['Enriched','Unenriched','Both_arms']:
 acc=set(mem.BaseAccession if arm=='Both_arms' else mem.loc[mem.EnrichmentStatus==arm,'BaseAccession'])
 rows.append(dict(Scope=arm,SampleID='GROUP_UNION',Proteins=len(acc),DDR=len(acc&ddr),Percent=100*len(acc&ddr)/len(acc)))
pd.DataFrame(rows).to_csv(out/'arm_sample_and_union.csv',index=False)
f=pd.read_csv(root/'data/publication_input/kla_protein_membership_31.csv')
frozen=set(f.loc[(f.PXD=='PXD028488')&(f.SampleGroup=='TALL-104'),'BaseAccession'])
acc=set(mem.BaseAccession)
pd.DataFrame([dict(OnlyFrozen=';'.join(sorted(frozen-acc)),OnlyRebuilt=';'.join(sorted(acc-frozen)),ExactSet=frozen==acc)]).to_csv(out/'frozen_set_check.csv',index=False)
print(pd.DataFrame(rows).to_string(index=False));print('Matches frozen',frozen==acc)
# Show which arm contributes the group-level proteins used in old DDR labels.
a=set(mem.loc[mem.EnrichmentStatus=='Enriched','BaseAccession']);b=set(mem.loc[mem.EnrichmentStatus=='Unenriched','BaseAccession'])
pd.DataFrame([dict(Set=n,Proteins=len(x),DDR=len(x&ddr)) for n,x in [('Enriched_only',a-b),('Unenriched_only',b-a),('Both',a&b)]]).to_csv(out/'arm_overlap.csv',index=False)
