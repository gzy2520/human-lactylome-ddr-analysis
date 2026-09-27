"""Bounded metadata inventory; eligibility is NOT inferred from keyword matches."""
import csv,json,pathlib,hashlib
root=pathlib.Path('audit/20260927_ddr_fraction');raw=root/'raw'
old=set()
with open('outputs/20260920_lactate_metabolism/tables/rna_proteome_mapping_31.csv') as f:
 for r in csv.DictReader(f):
  for k,v in r.items():
   if 'PXD' in k:old.update(v.split(';'))
projects={}
for name in ['pride_search','pride_search_page1','pride_lactylome']:
 for p in json.loads((raw/(name+'.json')).read_text()):projects[p['accession']]=p
reasons={
'PXD083647':'Targeted in-vitro phosphorylation reaction; not global Kla / ordinary proteome.',
'PXD080849':'Lactate-induced PGK1 targeted study; no complete untreated denominator/RNA established.',
'PXD077071':'FLAG-IP interactome with HNRNPU overexpression; not target fraction assay.',
'PXD073998':'Lactyl adenosine/tRNA metabolite assay; not protein lactylome.',
'PXD076410':'AARS1 inhibitor/structural experiment; not established global target fraction.',
'PXD070427':'Mixed human/mouse ZCWPW2 perturbation study; human untreated matched three-omics not established.',
'PXD068838':'Mixed human/mouse ZCWPW2 perturbation study; human untreated matched three-omics not established.',
'PXD069158':'TIGAR perturbation/mixed models; untreated human global denominator and RNA not established.',
'PXD045404':'Tufm IP interactome; not global target fraction.',
'PXD062720':'Candidate: bladder cancer EPI experiment; need untreated controls plus ordinary proteome and independent RNA.'}
rows=[]
for acc,p in sorted(projects.items()):
 human='sapiens' in str(p.get('organisms',[]));used=acc in old
 reason='Already used in original proteome/reference mapping.' if used else reasons.get(acc,'Not individually eligibility-reviewed; keyword match only.')
 rows.append(dict(Accession=acc,Title=p['title'],PublicationDate=p.get('publicationDate',''),HumanMetadata=human,AlreadyUsedProteome=used,Assessment=reason,URL='https://www.ebi.ac.uk/pride/archive/projects/'+acc))
with (root/'pride_inventory.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
urls={
'ddr_descendants':'https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/GO:0006974/descendants?relations=is_a,part_of',
'repair_descendants':'https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/GO:0006281/descendants?relations=is_a,part_of',
'pride_search':'https://www.ebi.ac.uk/pride/ws/archive/v2/search/projects?keyword=lactylation&pageSize=100&page=0',
'pride_search_page1':'https://www.ebi.ac.uk/pride/ws/archive/v2/search/projects?keyword=lactylation&pageSize=100&page=1',
'pride_lactylome':'https://www.ebi.ac.uk/pride/ws/archive/v2/search/projects?keyword=lactylome&pageSize=100&page=0'}
with (root/'retrieval_manifest.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=['File','URL','AccessDate','ValidJSON','SHA256']);w.writeheader()
 for name,url in urls.items():
  p=raw/(name+'.json');json.loads(p.read_text());w.writerow(dict(File=str(p),URL=url,AccessDate='2026-09-27',ValidJSON=True,SHA256=hashlib.sha256(p.read_bytes()).hexdigest()))
print('Inventory:',len(rows),'unique projects;',sum(x['HumanMetadata'] for x in rows),'with human metadata. Eligibility not inferred.')
