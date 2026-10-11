from pathlib import Path
import sys
import pandas as pd
project=Path(__file__).resolve().parents[3]
source=Path('/Users/gzy2520/Desktop/Research/kla/data')
out=Path(__file__).resolve().parent
sys.path.insert(0,str(project/'workflow/source_legacy'))
from python.utils.extractors import extract_pxd028488,extract_pxd053474_dda,extract_pxd053474_dia,extract_pxd053474_supplementary,reconcile_pxd053474
s,l,_=extract_pxd028488(source,source.parent,'10.1038/s41592-022-01523-1')
s=s[s.PrimaryIncluded & s.CellOrTissueType.isin(['HCT116','HEK293T','T-ALL'])]
s['CellOrTissueType']=s.CellOrTissueType.replace({'T-ALL':'TALL-104'})
s.to_csv(out/'pxd028488_source_sites.csv',index=False)
l.to_csv(out/'pxd028488_extraction_log.csv',index=False)
frames=[];logs=[]
for f in [extract_pxd053474_dda,extract_pxd053474_dia,extract_pxd053474_supplementary]:
 a,b=f(source,source.parent,'10.1021/jasms.4c00366');frames.append(a);logs.append(b)
s,c=reconcile_pxd053474(*frames)
s[s.PrimaryIncluded].to_csv(out/'pxd053474_source_sites.csv',index=False)
c.to_csv(out/'pxd053474_site_reconciliation.csv',index=False)
pd.concat(logs).to_csv(out/'pxd053474_extraction_log.csv',index=False)
print(s[s.PrimaryIncluded].groupby('CellOrTissueType').BaseAccession.nunique())

from python.utils.extractors import extract_pxd078013,extract_maxquant_site_table
s,l=extract_pxd078013(source,source.parent,'',{})
s.to_csv(out/'pxd078013_legacy_sites.csv',index=False)
l.to_csv(out/'pxd078013_legacy_log.csv',index=False)
for px,path in [('PXD060185','PXD060185/search_results/RESULT/combined/txt/La (K)Sites.txt'),('PXD078736','PXD078736/search_results/txt/La(K)Sites.txt')]:
 s,l=extract_maxquant_site_table(px,'',source/path,source.parent,'cell',{},'legacy_maxquant')
 s.to_csv(out/(px.lower()+'_legacy_sites.csv'),index=False)
 l.to_csv(out/(px.lower()+'_legacy_log.csv'),index=False)
