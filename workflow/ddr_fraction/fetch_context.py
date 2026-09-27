import urllib.request,json,pathlib,concurrent.futures
root=pathlib.Path('audit/20260927_ddr_fraction/raw')
jobs={
 'repair_descendants':'https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/GO:0006281/descendants?relations=is_a,part_of',
 'ddr_descendants':'https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/GO:0006974/descendants?relations=is_a,part_of',
 'pride_search':'https://www.ebi.ac.uk/pride/ws/archive/v2/search/projects?keyword=lactylation&pageSize=100&page=0',
}
def fetch(item):
 name,url=item
 try:
  req=urllib.request.Request(url,headers={'Accept':'application/json','User-Agent':'Kla-research-audit/1.0'})
  data=urllib.request.urlopen(req,timeout=45).read();obj=json.loads(data)
  (root/(name+'.json')).write_bytes(data)
  return {'name':name,'url':url,'ok':True,'bytes':len(data),'keys':list(obj)[:8] if isinstance(obj,dict) else ['list']}
 except Exception as e:return {'name':name,'url':url,'ok':False,'error':str(e)}
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex: result=list(ex.map(fetch,jobs.items()))
(root/'fetch_log.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
