#!/usr/bin/env python3
"""Recover source-resolved whole-proteome detection units without changing releases."""
import csv, gzip, io, json, re, zipfile, hashlib, statistics
from pathlib import Path
from collections import defaultdict, Counter

ROOT = Path(__file__).resolve().parents[3]
SOURCE = Path('/Users/gzy2520/Desktop/Research/kla/data')
OUT = Path(__file__).resolve().parent
csv.field_size_limit(10000000)
def read(path, sep=','):
    return list(csv.DictReader(open(path, encoding='utf-8-sig'), delimiter=sep))
def write(name, rows):
    if not rows: return
    f = gzip.open(OUT/name, 'wt') if name.endswith('.gz') else open(OUT/name, 'w')
    with f:
        w=csv.DictWriter(f, fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def positive(x):
    try: return float(x)>0
    except (ValueError,TypeError): return False
def valid(row):
    return not any(row.get(k,'').lower() in ['+','1','true','yes'] for k in ['Reverse','Potential contaminant','Contaminant','Only identified by site'])
pattern = re.compile(r'^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9](?:[A-Z0-9]{3}[0-9])?)$')
def accession(token):
    token=token.strip()
    if token.startswith(('CON__','REV__')): return ''
    if '|' in token: token=token.split('|')[1] if token.startswith(('sp|','tr|')) else token.split('|')[0]
    token=re.sub(r'-[0-9]+$','',token.removeprefix('NX_'))
    return token if pattern.match(token) else ''
def split(value):
    return {a for x in re.split('[;,:]',value) if (a:=accession(x))}
groups=read(ROOT/'data/publication_input/group_summary_31.csv')
membership=read(ROOT/'data/publication_input/reference_protein_membership_31.csv')
go=read(ROOT/'data/publication_input/human_ddr_go_annotations.tsv','\t')
ddr=set().union(*(split(x['GENE PRODUCT ID']) for x in go if x['GENE PRODUCT DB']=='UniProtKB' and x['TAXON ID'] in ['','9606'] and 'NOT' not in x['QUALIFIER'].split('|')))
ensmap=defaultdict(set)
for x in membership:
    if x['IdentifierType']=='ENSEMBLPROT':ensmap[x['SourceProteinID']].update(split(x['MappedBaseAccessions']))
def protein_ids(value):
    out=split(value)
    for t in re.split('[;,:]',value):
        t=re.sub(r'\.[0-9]+$','',re.sub(r'_RNA$','',t.strip()))
        if re.fullmatch(r'ENSP[0-9]+',t):out.update(ensmap[t])
    return out
frozen=defaultdict(set)
for x in membership:frozen[(x['PXD'],x['SampleGroup'])].update(split(x['MappedBaseAccessions']))
units=[];members=[];runs=[];inventory=[];manifests=[];mapping_issues=[];method_comparisons=[]
def add_manifest(path):
    if any(x['Path']==str(path) for x in manifests):return
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    manifests.append(dict(Path=str(path),Bytes=path.stat().st_size,SHA256=h.hexdigest()))
def emit(g,uid,ids,kind,source,columns,nruns,nfractions,note,role='primary'):
    ids=set(ids);assert ids
    row=dict(PXD=g['PXD'],SampleGroup=g['SampleGroup'],ReferencePXD=g['ReferencePXD'],Dataset='Whole proteome',SampleID=uid,ObservationType=kind,ProteinCount=len(ids),DdrProteinCount=len(ids&ddr),DdrFractionPercentage=100*len(ids&ddr)/len(ids),AcquisitionCount=nruns,FractionCount=nfractions,SourceFile=str(source),SourceColumns=columns,UnitDefinition=note,UseRole=role)
    units.append(row)
    members.extend(dict(PXD=g['PXD'],SampleGroup=g['SampleGroup'],SampleID=uid,BaseAccession=a,IsDdr=a in ddr,ObservationType=kind,UseRole=role) for a in sorted(ids))
    return row
def ziprows(z,name):return csv.DictReader(io.TextIOWrapper(z.open(name)),delimiter='\t')

# PXD010154: one named tissue specimen profile; all fraction acquisitions already combined.
for material,token in [('normal human lung','P013163'),('normal pregnancy placenta','P013680')]:
    g=next(x for x in groups if x['SampleGroup']==material)
    path=SOURCE/g['ReferenceEvidenceFile'].removeprefix('data/')
    rows=read(path,'\t');col='Intensity '+token;ids=set()
    seen_ens=set()
    for x in rows:
        if valid(x):
            ids.update(protein_ids(x['Majority protein IDs']))
            for t in re.split('[;,]',x['Majority protein IDs']):
                t=re.sub(r'\.[0-9]+$','',re.sub(r'_RNA$','',t.strip()))
                if re.fullmatch(r'ENSP[0-9]+',t):seen_ens.add(t)
    for t in sorted(seen_ens):
        aa=ensmap[t]
        if len(aa)!=1:mapping_issues.append(dict(PXD=g['PXD'],SampleGroup=material,SourceIdentifier=t,MappedBaseAccessions=';'.join(sorted(aa)),MappingCount=len(aa),Issue='unmapped' if not aa else 'one_to_many'))
    zipname=next((SOURCE/'PXD010154/search_results').glob(token+'*zip'))
    z=zipfile.ZipFile(zipname);sr=list(ziprows(z,'txt/summary.txt'));raws=[x for x in sr if x.get('Experiment')==token]
    emit(g,token,ids,'single_source_profile',path,'valid proteinGroups; '+col+' identifies sole experiment',len(raws),'merged','One named source tissue profile; valid identified proteinGroups, contaminants/reverse/site-only excluded. Original experiment combines all LC-MS fractions. Mapped stable UniProt denominator; unmapped ENSP retained in audit. No dispersion estimate.')
    propagated={a for x in membership if x['PXD']==g['PXD'] and x['SampleGroup']==material and x['IsDdr']=='TRUE' for a in split(x['MappedBaseAccessions'])}
    method_comparisons.append(dict(PXD=g['PXD'],SampleGroup=material,Method='legacy_ENSP_IsDdr_propagation',ProteinCount=len(ids),DdrProteinCount=len(ids&propagated),DifferenceFromRecommended=len(ids&propagated)-len(ids&ddr),Explanation='Legacy any-DDR ENSP boolean propagated to all mapped accessions; recommended classifies each UniProt directly against frozen GO'))
    inventory.append(dict(PXD=g['PXD'],SampleGroup=material,SourceFile=str(path),ResolvedColumns=col,RecoveredUnits=1,ObservationType='single_source_profile',Evidence='proteinGroups has one Intensity experiment column; archive summary lists '+str(len(raws))+' acquisition rows under '+token))
    for x in raws:runs.append(dict(PXD=g['PXD'],SampleGroup=material,SampleID=token,RawFile=x['Raw file'],PreparationID=token,SourceColumn=col,UnitEvidence='summary Experiment '+token))
    add_manifest(path);add_manifest(zipname)

# ProCan: three technical preparation IDs, each measured twice. Keep QC controls separate.
base=SOURCE/'PXD030304/search_results'
mapfile=base/'ProCan-DepMapSanger_mapping_file_replicates.txt'
mapping=read(mapfile,'\t'); byrun={x['Automatic_MS_filename']:x for x in mapping}
target={'MCF7','HCT-116','A549','MDA-MB-468','T47D','RKO'}
selected={x['Automatic_MS_filename']:x for x in mapping if x['Cell_line'] in target or 'HEK293T' in x['Cell_line']}
detected={}; path=base/'ProCan-DepMapSanger_peptide_counts_per_protein_per_sample.txt'
with open(path) as f:
    rd=csv.reader(f,delimiter='\t');header=next(rd);acc=[accession(x.split(';')[0]) for x in header[1:]]
    for x in rd:
        if x[0] in selected:detected[x[0]]={a for a,v in zip(acc,x[1:]) if a and positive(v)}
assert set(detected)==set(selected)
alias={'MCF7':'MCF7','HCT116':'HCT-116','A549':'A549','MDA-MB-468':'MDA-MB-468','T-47D':'T47D','RKO WT and GSK3B KO':'RKO'}
for g in groups:
    if g['ReferencePXD']!='PXD030304':continue
    line=alias.get(g['SampleGroup'])
    sample_runs=defaultdict(list)
    if line:
        for run,x in selected.items():
            if x['Cell_line']!=line:continue
            parts=run.split('_')
            prep=parts[4] if '_e0022_' in run else parts[2]
            sample_runs[prep].append(run)
        assert len(sample_runs)==3 and all(len(rr)==2 for rr in sample_runs.values())
    else:
        for run,x in selected.items():
            if 'HEK293T' in x['Cell_line']:sample_runs[x['Cell_line']].append(run)
    for prep,rr in sorted(sample_runs.items()):
        ids=set().union(*(detected[r] for r in rr));uid='PXD030304_'+(line+'_'+prep if line else prep)
        kind='technical_preparation' if line else 'qc_source_pool'
        role='primary' if line or prep=='Control_HEK293T_lys' else 'qc_alternative'
        note='Union of duplicate acquisitions sharing preparation token; three technical preparations per cell line.' if line else 'Union across QC acquisitions for this explicitly labelled control source; lys is primary, H002/H003 standards retained as alternatives.'
        emit(g,uid,ids,kind,path,'positive peptide counts',len(rr),1,note,role)
        for run in rr:
            runs.append(dict(PXD=g['PXD'],SampleGroup=g['SampleGroup'],SampleID=uid,RawFile=run,PreparationID=prep,SourceColumn='Run='+run,UnitEvidence='replicate mapping and shared preparation token' if line else selected[run]['Cell_line']))
    inventory.append(dict(PXD=g['PXD'],SampleGroup=g['SampleGroup'],SourceFile=str(path),ResolvedColumns='Run with positive peptide counts; mapping Cell_line + Automatic_MS_filename',RecoveredUnits=3 if line else 1,ObservationType='technical_preparation' if line else 'qc_source_pool',Evidence='PXD030304 projectDescription: technical triplicate with duplicate instrument acquisition; HEK controls separated by lys/H002/H003 source label'))
add_manifest(path);add_manifest(mapfile)

# PC-3M: heavy channel; 1st and 2nd are biological repeats in published experiment.
g=next(x for x in groups if x['SampleGroup']=='PC-3M');path=SOURCE/g['ReferenceEvidenceFile'].removeprefix('data/');z=zipfile.ZipFile(path)
rows=list(ziprows(z,'proteinGroups.txt'));sr=list(ziprows(z,'summary.txt'))
for rep in ['1st','2nd']:
    col='Intensity H '+rep;ids=set()
    for x in rows:
        if valid(x) and positive(x[col]):ids.update(split(x['Majority protein IDs']))
    rr=[x for x in sr if x.get('Experiment')==rep]
    emit(g,'PXD022005_PC3M_H_'+rep,ids,'sample',path,col,len(rr),'merged','Heavy-labelled PC-3M biological repeat '+rep+'; proteinGroups experiment intensity integrates fractions. Publication: PMC10372916, results SILAC paragraph confirms two biological repeats.')
    for x in rr:runs.append(dict(PXD=g['PXD'],SampleGroup=g['SampleGroup'],SampleID='PXD022005_PC3M_H_'+rep,RawFile=x['Raw file'],PreparationID=rep,SourceColumn=col,UnitEvidence='summary Experiment '+rep+'; PC-3M heavy channel confirmed in PMC10372916'))
inventory.append(dict(PXD=g['PXD'],SampleGroup=g['SampleGroup'],SourceFile=str(path),ResolvedColumns='Intensity H 1st;Intensity H 2nd',RecoveredUnits=2,ObservationType='sample',Evidence='https://pmc.ncbi.nlm.nih.gov/articles/PMC10372916/ confirms PC-3M heavy and two biological replicates; archive summary links all fraction runs to 1st/2nd'))
add_manifest(path)

# MCF10A: one source preparation, ten IEF fractions, each measured twice.
g=next(x for x in groups if x['SampleGroup']=='MCF10A');path=SOURCE/g['ReferenceEvidenceFile'].removeprefix('data/');z=zipfile.ZipFile(path)
sr=list(ziprows(z,'summary.txt'));rr=[x['Raw file'] for x in sr if '_10A_' in x['Raw file']]
assert len(rr)==20
rrset=set(rr);ids=set();perrun=defaultdict(set)
for x in ziprows(z,'evidence.txt'):
    if x['Raw file'] in rrset and valid(x):
        aa=split(x['Leading Razor Protein']);ids.update(aa);perrun[x['Raw file']].update(aa)
emit(g,'PXD002400_MCF10A_20111109',ids,'single_source_profile',path,'evidence.txt: Raw file *_10A_*; Leading Razor Protein; valid identified rows',20,10,'Union of 10 IEF fractions and two technical LC-MS injections per fraction. Leading Razor Protein represents peptide assignment; accept both 6- and 10-character stable UniProt IDs. One preparation. Protein FDR parameter is 1; evidence-based detection profile.')
old_ids={a for a in ids if len(a)==6}
method_comparisons.append(dict(PXD=g['PXD'],SampleGroup=g['SampleGroup'],Method='legacy_6char_UniProt_regex',ProteinCount=len(old_ids),DdrProteinCount=len(old_ids&ddr),DifferenceFromRecommended=len(old_ids)-len(ids),Explanation='Legacy regex excluded valid 10-character UniProt IDs; recommended retains all valid leading-razor stable IDs'))
for run in rr:runs.append(dict(PXD=g['PXD'],SampleGroup=g['SampleGroup'],SampleID='PXD002400_MCF10A_20111109',RawFile=run,PreparationID='20111109_10A',SourceColumn='Raw file='+run,UnitEvidence='sampleProcessingProtocol duplicate acquisitions; F1-F10 fractions'))
inventory.append(dict(PXD=g['PXD'],SampleGroup=g['SampleGroup'],SourceFile=str(path),ResolvedColumns='evidence.txt Raw file/Leading Razor Protein; 20 *_10A_* runs',RecoveredUnits=1,ObservationType='single_source_profile',Evidence='project.json describes IEF and duplicate LC-MS; summary has 10 fractions x 2 runs; Protein FDR parameter 1'))
add_manifest(path)

write('recovered_unit_values.csv',units);write('recovered_protein_membership.csv.gz',members);write('run_to_unit_registry.csv',runs);write('source_inventory.csv',inventory)
write('identifier_mapping_issues.csv',mapping_issues);write('legacy_method_comparison.csv',method_comparisons)
comparisons=[]
for key in sorted({(x['PXD'],x['SampleGroup']) for x in units}):
    uu=[x for x in units if (x['PXD'],x['SampleGroup'])==key and x['UseRole']=='primary']
    sets={x['BaseAccession'] for x in members if (x['PXD'],x['SampleGroup'])==key and x['UseRole']=='primary'}
    fr=frozen[key];vv=[x['DdrFractionPercentage'] for x in uu];sd=statistics.stdev(vv) if len(vv)>1 else ''
    comparisons.append(dict(PXD=key[0],SampleGroup=key[1],N=len(uu),ObservationType=uu[0]['ObservationType'],Mean=sum(vv)/len(vv),SD=sd,SEM=sd/(len(vv)**.5) if sd!='' else '',RecoveredUnion=len(sets),FrozenUnion=len(fr),NewVsFrozen=len(sets-fr),MissingVsFrozen=len(fr-sets),RecoveredDDR=len(sets&ddr)))
write('recovery_summary.csv',comparisons)
for p in [ROOT/'data/publication_input/group_summary_31.csv',ROOT/'data/publication_input/human_ddr_go_annotations.tsv', ROOT/'data/publication_input/reference_protein_membership_31.csv', *[SOURCE/x/'metadata/project.json' for x in ['PXD002400','PXD022005','PXD030304']]]:add_manifest(p)
write('source_sha256.csv',manifests)
assert len(inventory)==12
assert len(members)==sum(x['ProteinCount'] for x in units)
assert all(x['DdrProteinCount']<=x['ProteinCount'] for x in units)
print(json.dumps(dict(units=len(units),primary_units=sum(x['UseRole']=='primary' for x in units),groups=len(inventory),membership_rows=len(members),DDR_annotation_accessions=len(ddr)),indent=2))
for x in comparisons:print(x)
