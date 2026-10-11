from pathlib import Path
import pandas as pd,numpy as np,json,re
out=Path(__file__).resolve().parent;root=out.parents[2];raw=Path('/Users/gzy2520/Desktop/Research/kla/data')
m=pd.read_csv(out/'independent_membership.csv.gz',dtype={'SampleID':str});e=pd.read_csv(out/'raw_row_identifier_evidence.csv.gz',keep_default_na=False);st=pd.read_csv(out/'sample_count_comparison.csv',dtype={'SampleID':str});rules=pd.read_csv(out/'raw_source_rules.csv');audit=pd.read_csv(out/'group_algorithm_audit.csv')
# Raw-record ID syntax audit, across actual retained source observations.
syn=[]
for (px,g),x in e.groupby(['PXD','SampleGroup']):
 vals=x.RawAccession
 syn.append(dict(PXD=px,SampleGroup=g,SourceRecords=len(x),DistinctRawIDStrings=vals.nunique(),ColonRecords=vals.str.contains(':',regex=False).sum(),SemicolonRecords=vals.str.contains(';',regex=False).sum(),CommaRecords=vals.str.contains(',',regex=False).sum(),PipeRecords=vals.str.contains('|',regex=False).sum(),IsoformRecords=vals.str.contains(r'-\d+(?:[|;,:]|$)',regex=True).sum(),PrefixedReverseOrContaminantRecords=vals.str.contains(r'(?:REV__|CON__)',regex=True).sum(),IdentifierParserChangedRecords=x.Changed.sum()))
pd.DataFrame(syn).to_csv(out/'identifier_syntax_audit.csv',index=False)
e[e.Changed].to_csv(out/'identifier_changed_records.csv',index=False)
# Compare independent MaxQuant/RKO parsing against freshly executed original Python parsers.
lc=[]
for px in ['PXD060185','PXD078736','PXD078013']:
 leg=pd.read_csv(out/(px.lower()+'_legacy_sites.csv'),dtype=str).fillna('')
 for sample,x in m[m.PXD.eq(px)].groupby('SampleID'):
  a=set(x.BaseAccession);b=set(leg.loc[leg.SampleName.eq(sample),'BaseAccession'])
  lc.append(dict(PXD=px,SampleID=sample,IndependentCount=len(a),OriginalParserCount=len(b),OnlyIndependent=';'.join(sorted(a-b)),OnlyOriginalParser=';'.join(sorted(b-a)),ExactSet=a==b))
pd.DataFrame(lc).to_csv(out/'fresh_legacy_set_comparison.csv',index=False)
# Exact offending raw records, preserving the source flags and identity columns.
p=raw/'PXD078736/search_results/txt/La(K)Sites.txt';d=pd.read_csv(p,sep='\t',dtype=str).fillna('');z=d[d.Protein.eq('REV__Q14093')].copy();z.insert(0,'FileLine',z.index+2);z.insert(0,'SourceFile',str(p));z.to_csv(out/'confirmed_reverse_record.csv',index=False)
correct=m[m.PXD.eq('PXD078736')].copy();correct.to_csv(out/'corrected_kla31_28_membership.csv.gz',index=False)
st[st.PXD.eq('PXD078736')].to_csv(out/'corrected_kla31_28_sample_counts.csv',index=False)
impact=[]
for label,vals in [('current',st[st.PXD.eq('PXD078736')].Fraction),('corrected',st[st.PXD.eq('PXD078736')].Percent)]:impact.append(dict(Version=label,N=len(vals),Mean=vals.mean(),SD=vals.std(ddof=1),SEM=vals.std(ddof=1)/np.sqrt(len(vals))))
pd.DataFrame(impact).to_csv(out/'kla31_28_group_impact.csv',index=False)
# Detailed one-row-per-position deliverable; status reflects independent reconstruction, not membership-only arithmetic.
notes={
'KLA31_13':'Protein representative; sample A; exact original MaxQuant parser set',
'KLA31_14':'PEAKS PSM/protein-peptide position mapping; one enrichment run',
'KLA31_15':'Original S3 plus DDA/DIA site corroboration; compartments/fractions collapsed',
'KLA31_16':'R.Condition NC116/R116 each one raw file; no within-condition replicates',
'KLA31_17':'Full original PEAKS site parser; enriched arm 196/201/198 proteins, 21/25/22 DDR; unenriched arm separately audited, not combined into enriched sample IDs',
'KLA31_18':'Unflagged Protein REV__O14732 hidden by canonical Proteins ID; remove original row before ID parsing; SIRT3KO rep1/2 each denominator -1',
'KLA31_19':'Author Protein accession; per-sample Intensity>0',
'KLA31_20':'Protein representative; sample B; exact original MaxQuant parser set',
'KLA31_21':'Protein representative; sample D; exact original MaxQuant parser set',
'KLA31_22':'Proteins list expanded; localization>0 and valid site row; replicates 1/2/3',
'KLA31_23':'Six model columns; model variation, no within-model replicate identity',
'KLA31_24':'Original-order majority representative; evidence site IDs joined to matching numeric position pairs; original parser exact set',
'KLA31_25':'HCD collision energies and enriched/nonenriched evidence consolidated to one technical sourceunit',
'KLA31_26':'Unflagged Protein REV__Q14093; remove original row; H0/H24 each denominator -1; distinct conditions',
'KLA31_27':'Unflagged Protein REV__Q14093; remove original row; denominator -1; n=1',
'KLA31_28':'Unflagged REV__Q14093 rescued by prefix-stripping in current/legacy; exclude ctr_1/ctr_2 records. man condition biological meaning unresolved in dataset metadata',
'KLA31_29':'Protein representative; sample C; exact original MaxQuant parser set',
'KLA31_30':'Two model columns; no within-model replicate identity',
'KLA31_31':'Six raw-file precursor columns; K(UniMod:2114); 3 A0h and 3 A6h samples',
}
audit['Evidence']=audit.GroupID.map(notes);audit['IdentifierPolicy']='valid 6/10-character UniProt; preserve source-order representative when required; split ;,: and handle sp/tr pipes; collapse isoforms; exclude prefixed decoys'
audit['NumeratorDenominator']='DDR is subset of independently recovered per-unit protein set; same denominator'
audit['SourceFiles']=[';'.join(sorted(set(rules.loc[(rules.PXD==r.PXD)&(rules.SampleGroup==r.SampleGroup),'SourceFile'])))for _,r in audit.iterrows()]
audit.to_csv(out/'group_algorithm_audit.csv',index=False)
assert len(audit)==19 and len(st)==54
assert (m.groupby(['PXD','SampleGroup','SampleID','BaseAccession']).size()==1).all()
assert (st.DDR<=st.Proteins).all()
print(pd.DataFrame(lc).to_string(index=False));print(pd.DataFrame(impact).to_string(index=False));print('changed record:',z[['FileLine','id','Protein','Reverse','Potential contaminant']].to_string(index=False))

changed=st[st.DeltaProteins.ne(0)|st.DeltaDDR.ne(0)]
affected=sorted(changed.GroupID.unique())
m[m.GroupID.isin(affected)].to_csv(out/'corrected_affected_membership.csv.gz',index=False)
changed.to_csv(out/'corrected_affected_sample_counts.csv',index=False)
impact=[]
for gid,x in st[st.GroupID.isin(affected)].groupby('GroupID'):
 for label,c in [('current','Fraction'),('corrected','Percent')]:
  v=x[c];impact.append(dict(GroupID=gid,Version=label,N=len(v),Mean=v.mean(),SD=v.std(ddof=1),SEM=v.std(ddof=1)/np.sqrt(len(v))))
pd.DataFrame(impact).to_csv(out/'affected_group_mean_sem.csv',index=False)
confirmed=[];support=[]
for px,rel,gene in [('PXD028737','search_results/extracted_pairing/combined/combined/txt/La(K)Sites.txt','Q14093'),('PXD050147','search_results/Lactyl_K_Sites.txt','O14732'),('PXD058534','search_results/extracted_pairing/txt/La (K)Sites.txt','Q14093'),('PXD078736','search_results/txt/La(K)Sites.txt','Q14093')]:
 p=raw/px/rel;d=pd.read_csv(p,sep='\t',dtype=str).fillna('')
 rr=d[d.Protein.str.startswith('REV__') & d.Protein.str.contains(gene)].copy();rr.insert(0,'FileLine',rr.index+2);rr.insert(0,'SourceFile',str(p));confirmed.append(rr)
 for _,c in changed[changed.PXD.eq(px)].iterrows():
  x=e[(e.PXD==px)&(e.SampleID==c.SampleID)]
  forward=x[x.Parsed.str.split(';').apply(lambda ids:gene in ids)]
  support.append(dict(PXD=px,SampleID=c.SampleID,RemovedID=gene,ForwardSupportingEligibleRows=len(forward),ForwardSourceRows=';'.join(forward.SourceRow.astype(str))))
pd.concat(confirmed,ignore_index=True).to_csv(out/'confirmed_all_reverse_rows.csv',index=False)
pd.DataFrame(support).to_csv(out/'affected_id_forward_support.csv',index=False)
old=pd.read_csv(root/'data/candidate/sample_boxplot_values.csv',dtype={'SampleID':str})
old.loc[(old.PXD=='PXD028488')&(old.SampleGroup=='HCT116'),'SampleID']='HCT116_1'
old.loc[(old.PXD=='PXD053474')&(old.SampleGroup=='HCT116'),'SampleID']='HCT116_technical_merged'
old.loc[(old.PXD=='PXD028488')&(old.SampleGroup=='HEK293T'),'SampleID']='HEK293T_technical_merged'
for i in [1,2,3]:old.loc[(old.PXD=='PXD028488')&(old.SampleGroup=='TALL-104')&(old.SampleID==f'Sample {i}'),'SampleID']=f'TALL-NaLAC-enrichment-{i}.raw'
old=old.rename(columns={'KlaProteinCount':'OldR_ProteinCount','KlaDdrProteinCount':'OldR_DdrProteinCount','KlaDdrFractionPercentage':'OldR_Percent'})
st.merge(old[['PXD','SampleGroup','SampleID','OldR_ProteinCount','OldR_DdrProteinCount','OldR_Percent']],on=['PXD','SampleGroup','SampleID'],how='left',validate='one_to_one').to_csv(out/'old_R_current_independent_threeway.csv',index=False)
