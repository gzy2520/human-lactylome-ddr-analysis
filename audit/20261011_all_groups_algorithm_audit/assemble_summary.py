from pathlib import Path
import csv
root=Path(__file__).resolve().parent
release=root.parents[1]/'deliverables/20261011_lactylation_ddr_all_groups_audited'
def read(p):
 with Path(p).open(encoding='utf-8-sig') as f:return list(csv.DictReader(f))
base=read(root/'cross_group/historical_union_vs_current_mean_62.csv')
new={(x['GroupID'],x['Dataset']):x for x in read(release/'tables/ddr_mean_sem_62.csv')}
reviews={}
for area,file in [('kla_tissues','group_review.csv'),('kla_cells','group_algorithm_audit.csv'),('whole_proteome','all_31_whole_proteome_review.csv')]:
 for x in read(root/area/file):reviews[(x['GroupID'],x['Dataset'])]=(area,file,x)
assert len(reviews)==62 and len(new)==62
rows=[]
for b in sorted(base,key=lambda x:(int(x['X']),x['Dataset'])):
 key=(b['GroupID'],b['Dataset']);area,file,x=reviews[key];n=new[key]
 action='原表复算一致，保留当前集合与单位'
 if n['UnitMark']=='p':action='完整处理矩阵改为1个panel；保留来源列ID；不画样本SEM'
 elif key[0] in ['KLA31_18','KLA31_26','KLA31_27','KLA31_28'] and key[1]=='Lactylome (Kla)':action='剔除无forward支持的反向代表序列，重新计算样本分母'
 if key[0]=='KLA31_28':action+='；显示条件采用源标签ctr/man，生物释义待核'
 rows.append(dict(GroupID=key[0],PXD=b['PXD'],SampleGroup=b['SampleGroup'],Dataset=key[1],RawReviewStatus=x['Status'],HistoricalUnionProteins=b['HistoricalN'],HistoricalUnionDDR=b['HistoricalDdr'],HistoricalUnionPercent=b['HistoricalFraction'],PreAuditMean=b['Mean'],PostAuditMean=n['Mean'],MeanChange=float(n['Mean'])-float(b['Mean']),AnalysisUnits=n['N'],UnitTypes=n['UnitTypes'],PlottedSEM=n['PlotSEM'],FinalAction=action,IdentifierReview=x.get('IdentifierConclusion',x.get('IdentifierPolicy','独立检查原始ID全部字段：多分隔符、pipe、isoform、6/10字符、污染物/反向代表')),NumeratorDenominatorReview=x.get('NumeratorDenominatorConclusion',x.get('NumeratorDenominator','同一原表重建集合的DDR子集/全部canonical UniProt')),SampleUnitReview=x.get('ObservationUnitConclusion',x.get('Recommendation',x.get('Evidence',''))),SourceEvidence=x.get('SourceFiles',x.get('Evidence','')),DetailFile=str(Path(area)/file)))
with (root/'all_62_algorithm_comparison.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
assert sum(abs(r['MeanChange'])>1e-10 for r in rows)==4
print('All 62 positions reviewed; 4 mean changes; 4 panel conversions.')
