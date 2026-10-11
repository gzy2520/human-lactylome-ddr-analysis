import csv,gzip,statistics,hashlib
from pathlib import Path
from collections import defaultdict
O=Path(__file__).parent;W=O.parents[2]
def rd(p):
 with (gzip.open(p,'rt') if str(p).endswith('.gz') else open(p)) as f:yield from csv.DictReader(f)
def wr(n,r):
 r=list(r)
 with (gzip.open(O/n,'wt') if n.endswith('.gz') else open(O/n,'w')) as f:w=csv.DictWriter(f,fieldnames=list(r[0]));w.writeheader();w.writerows(r)
def key(x):return (x['PXD'],x['SampleGroup'],x['SampleID'])
a=defaultdict(set);ad=defaultdict(set)
for n in ['corrected_direct_source_membership.csv.gz','recovered_protein_membership.csv.gz']:
 for x in rd(O/n):
  if x.get('UseRole','primary')!='primary':continue
  a[key(x)].add(x['BaseAccession'])
  if x['IsDdr'].lower()=='true':ad[key(x)].add(x['BaseAccession'])
b=defaultdict(set);bd=defaultdict(set)
for x in rd(W/'deliverables/20261011_lactylation_ddr_tall_verified/inputs/ddr_protein_membership.csv.gz'):
 if x['Dataset']!='Whole proteome':continue
 b[key(x)].add(x['BaseAccession'])
 if x['IsDdr'].lower()=='true':bd[key(x)].add(x['BaseAccession'])
cur=[x for x in rd(W/'deliverables/20261011_lactylation_ddr_tall_verified/tables/ddr_source_unit_values.csv') if x['Dataset']=='Whole proteome']
old=[x for x in rd(W/'data/candidate/figure1_sample_boxplot_values.csv') if x['Dataset']=='Whole proteome']
frozen={(x['PXD'],x['SampleGroup']):x for x in rd(W/'data/publication_input/group_summary_31.csv')}
notes={
1:'全部有效鉴定行3984/154；六reporter通道正值并集3574/140；373行六通道均无正值。3 RCT列的positive mask相同，零SD来自同一检出集合。',
2:'原计数对象为ENSP及其他原始SourceProteinID，现为映射后去重UniProt；一对多映射逐UniProt direct GO分类；一个组织profile合并36采集。',
3:'PD Found in Sample HSP1–4判断检出，Accession稳定ID；未见隐藏前缀代表污染。',
4:'PD Found in Sample NSP1–4判断检出，与scar分别筛列。',
5:'作者完整处理矩阵6082行三列全正，最小值0.00262465560808778重复84/69/594次；原文列说明三生物重复relative intensity。按一个panel集合展示并省略样本SEM。',
6:'原计数对象为ENSP及其他原始SourceProteinID，现为映射后去重UniProt；一个组织profile合并36采集；旧IsDdr一对多传播会多计68 DDR。',
7:'原全部行10680/461，33个L/N定量列正值并集9492/385，1188 ID无任何正定量；保留L/N条件身份，跨条件分散度需明确。',
8:'NAT1–5正值；蛋白组展开有效稳定ID，独立于PCa列。',
9:'ANTs 4254行8列全正；Non-rec1全列1，为作者比值矩阵。按一个panel集合展示，省略样本SEM；同名列与CISs需用Sheet区分。',
10:'原workbook 2.a protein raw information，94个ESCC-xxxT列；正值ID直查冻结GO。',
11:'PCa1–5正值，独立于NAT列。',
12:'CISs 5198行8列全正；Non-rec1全列1，为作者比值矩阵。按一个panel集合展示，省略样本SEM；同名列与ANTs需用Sheet区分。',
13:'ProCan三个技术制备，每制备两个仪器采集合并；替代旧averaged feature matrix集合。',
14:'ProCan三个技术制备，每制备两个仪器采集合并；与位置15共享真实参考。',
15:'ProCan三个技术制备，每制备两个仪器采集合并；与位置14共享真实参考。',
16:'NC116与R116两列为control/co-culture两个条件，各一个run；跨条件SD不应标重复SEM。',
17:'PEAKS nonenrichment proteins.csv Area Sample1–3，Whole accession字段无colon解析差异；与Kla PSM列colon问题不同。',
18:'WT/SIRT1KO/SIRT3KO各3定量列；旧鉴定并集6772/450，当前正值并集6713/448；59 ID在9列均无正值；总体分散度含条件差异。',
19:'ProCan三个技术制备，每制备两个仪器采集合并。',
20:'ProCan三个技术制备，每制备两个仪器采集合并。',
21:'ProCan三个技术制备，每制备两个仪器采集合并。',
22:'PC3M heavy Intensity H 1st/2nd两生物重复，10/26采集先合并；Majority有效ID，union8049/262。',
23:'G2907/G3028/G3264/GSC23/MES28/RKI六个细胞模型；模型差异应标model dispersion。',
24:'ProCan三个技术制备；来源RKO cell line，不将GSK3B KO标签移植到参考测量。',
25:'ProCan Control_HEK293T_lys 401 QC采集合并为一个source pool；H002/H003分开作替代；不限制到旧6692 averaged universe。',
26:'H0与H24两个条件的LFQ正值，非两独立重复；LFQ缺失与原Intensity检出需区分。',
27:'HK2 amostra1/3/4正值三参考；与位置28真实参考复用。',
28:'HK2 amostra1/3/4正值三参考；与位置27真实参考复用。',
29:'MCF10A一个preparation，10fraction各2采集合并；Leading Razor代表，6/10字符4920/242；旧6字符4839/242。',
30:'ENSA/HMP1两个NSC模型，按model dispersion标注。',
31:'HUVEC A0h1–3 control列，不混入A6h感染列。'}
reviews=[];unit=[]
for k in sorted(set(a)|set(b)):
 unit.append(dict(PXD=k[0],SampleGroup=k[1],SampleID=k[2],RawRebuiltN=len(a[k]),CurrentN=len(b[k]),RawDDR=len(ad[k]),CurrentDDR=len(bd[k]),Added=';'.join(sorted(a[k]-b[k])),Removed=';'.join(sorted(b[k]-a[k])),MembershipEqual=a[k]==b[k],DdrEqual=ad[k]==bd[k]))
for g in sorted({x['GroupID']for x in cur}):
 cc=[x for x in cur if x['GroupID']==g];p,s=cc[0]['PXD'],cc[0]['SampleGroup'];n=int(g.split('_')[1]);ff=frozen[(p,s)];oo=[x for x in old if x['PXD']==p and x['SampleGroup']==s];kk=[key(x)for x in cc];un=set().union(*(a[k]for k in kk));du=set().union(*(ad[k]for k in kk));vs=[100*len(ad[k])/len(a[k])for k in kk];same=all(a[k]==b[k]and ad[k]==bd[k]for k in kk)
 reviews.append(dict(GroupID=g,PXD=p,SampleGroup=s,Dataset='Whole proteome',Status='matrix_panel_correction' if n in [5,9,12] else 'unit_interpretation_review' if n in [7,16,18,23,26,30] else 'source_reconstruction_pass',IdentifierConclusion='Stable canonical UniProt; full valid 6/10-character IDs; frozen human non-NOT GO direct matching',NumeratorDenominatorConclusion='Unique DDR UniProt / unique quantified or source-identified UniProt; historical source-union differences explained in Review',ObservationUnitConclusion=notes[n],HistoricalUnionN=ff['ReferenceProteinCount'],HistoricalUnionDdr=ff['ReferenceDdrProteinCount'],HistoricalFraction=100*float(ff['ReferenceDdrFraction']),OldSampleRows=len(oo),OldObservationTypes=';'.join(sorted({x['ObservationType']for x in oo})),OldSampleMean=statistics.mean(float(x['DdrFractionPercentage'])for x in oo)if oo else'',RebuiltUnits=len(kk),RebuiltUnionN=len(un),RebuiltUnionDdr=len(du),RebuiltMean=statistics.mean(vs),RebuiltSD=statistics.stdev(vs)if len(vs)>1 else'',CurrentMembershipAndDDREqual=same,RecommendedUnits=1 if n in [5,9,12]else len(kk),RecommendedSem='omit_processed_complete_matrix'if n in [5,9,12]else'omit_across_conditions'if n in [16,26]else'technical_preparation_SEM'if n in [13,14,15,19,20,21,24]else'condition_or_model_dispersion'if n in [7,18,23,30]else'no_SEM_single_source'if len(kk)==1 else'source_unit_SEM',SourceFiles=';'.join(sorted({x['SourceFile']for x in cc})),Review=notes[n]))
wr('all_31_whole_proteome_review.csv',reviews);wr('raw_rebuild_vs_current_unit_comparison.csv',unit)
wr('raw_rebuilt_whole_membership.csv.gz',(dict(PXD=k[0],SampleGroup=k[1],SampleID=k[2],BaseAccession=v,IsDdr=v in ad[k])for k in sorted(a)for v in sorted(a[k])))
print('Groups',len(reviews),'units',len(unit),'equal',sum(x['MembershipEqual']and x['DdrEqual']for x in unit));print([x for x in unit if not x['MembershipEqual']or not x['DdrEqual']])
