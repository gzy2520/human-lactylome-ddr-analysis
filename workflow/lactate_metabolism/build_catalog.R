#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table);library(jsonlite);library(AnnotationDbi);library(org.Hs.eg.db)})
set.seed(25)
cfg <- 'config/lactate_metabolism'; raw <- file.path(cfg,'raw')
files <- list.files(raw,pattern='_annotations_[0-9]+.json$',full.names=TRUE)
stopifnot(length(files)>=6)
anns <- rbindlist(lapply(files,function(f){
 x<-fromJSON(f,simplifyVector=FALSE);t<-sub('_annotations_.*','',basename(f));t<-sub('_',':',t)
 rbindlist(lapply(x$results,function(y)data.table(Term=t,AnnotationID=y$id,UniProt=sub('UniProtKB:','',y$geneProductId),Qualifier=y$qualifier,AnnotatedGO=y$goId,Evidence=y$goEvidence,ECO=y$evidenceCode,Reference=y$reference,Taxon=y$taxonId,AssignedBy=y$assignedBy,AnnotationDate=y$date,SourceSymbol=y$symbol)))
}))
stopifnot(all(anns$Taxon==9606))
anns[,ExcludedNOT:=grepl('NOT',Qualifier)];anns[,DirectAnnotation:=Term==AnnotatedGO]
fwrite(anns,file.path(cfg,'go_annotations.csv'))
uni<-unique(anns$UniProt);valid<-intersect(uni,keys(org.Hs.eg.db,keytype='UNIPROT'))
map<-as.data.table(select(org.Hs.eg.db,keys=valid,keytype='UNIPROT',columns=c('ENSEMBL','ENTREZID','SYMBOL')))
setnames(map,c('UNIPROT','ENSEMBL','ENTREZID','SYMBOL'),c('UniProt','Ensembl','Entrez','Symbol'))
map<-unique(map)
current_path<-file.path(cfg,'go_current_uniprot_mapping.csv')
if(file.exists(current_path)){
 current<-fread(current_path,colClasses='character')
 map[,Entrez:=as.character(Entrez)]
 map<-unique(rbind(current,map[!UniProt %in% current$UniProt],fill=TRUE))
}
fwrite(map,file.path(cfg,'go_id_mapping.csv'))
# Freeze curated candidates independently of GO annotations; symbols are resolved once here.
mods<-list(
 'L-lactate conversion'=c('LDHA','LDHB','LDHC'),
 'Lactate transport'=c('SLC16A1','SLC16A3','SLC16A7','SLC16A8','SLC5A8','SLC5A12','BSG','EMB'),
 'Pyruvate mitochondrial entry'=c('MPC1','MPC2','PDHA1','PDHB','DLAT','DLD','PDK1','PDK2','PDK3','PDK4','PDP1','PDP2'),
 'Gluconeogenic branch'=c('PC','PCK1','PCK2','FBP1','FBP2','G6PC1','G6PC3'),
 'Lactyl substrate formation'=c('ACSS2','SUCLG1','SUCLG2','AARS1','AARS2'),
 'D-lactate and glyoxalase'=c('GLO1','HAGH','LDHD'),
 'Glycogen upstream'=c('PYGL','PYGM','PYGB','AGL','PGM1'),
 'Fructose upstream'=c('KHK','ALDOB','TKFC'),
 'Glycerol upstream'=c('GK','GPD1','GPD2'),
 'Alanine pyruvate link'=c('GPT','GPT2'))
reactions<-c('Pyruvate + NADH <-> L-lactate + NAD+','Transmembrane monocarboxylate transport or transporter support','Pyruvate import; acetyl-CoA formation; PDH regulation','Pyruvate to glucose pathway: bypass enzymes; members are not equivalent','L-lactate activation / lactyl donor formation reported in specific models','Methylglyoxal detoxification; D-lactate generation or oxidation','Glycogen -> glucose-1-P -> glucose-6-P','Fructose -> triose intermediates','Glycerol -> glycerol-3-P <-> DHAP','Alanine <-> pyruvate')
refs<-c('UniProt curated function','UniProt curated function; PMID:11101640','UniProt curated function; Reactome R-HSA-70268','UniProt curated function; PMID:11788376','PMID:39561764; PMID:39642882; DOI:10.1038/s41586-024-07992-y','PMID:31767537; PMID:30931947','Reactome:R-HSA-70221','Reactome:R-HSA-70350; PMID:28425966','PMID:35562085; UniProt curated function','UniProt curated function')
cur<-rbindlist(lapply(seq_along(mods),function(i)data.table(CandidateSymbol=mods[[i]],Module=names(mods)[i],Reaction=reactions[i],EvidenceSource=refs[i])))
# Resolve current official labels using the installed frozen annotation snapshot.
cur[,LookupSymbol:=CandidateSymbol]
cm<-as.data.table(select(org.Hs.eg.db,keys=unique(cur$LookupSymbol),keytype='SYMBOL',columns=c('ENSEMBL','ENTREZID','UNIPROT')))
setnames(cm,c('SYMBOL','ENSEMBL','ENTREZID','UNIPROT'),c('LookupSymbol','Ensembl','Entrez','UniProt'))
cur<-merge(cur,cm,by='LookupSymbol',all.x=TRUE,allow.cartesian=TRUE)
cur[,Symbol:=CandidateSymbol]
cur[,Role:='enzyme'];cur[Symbol %in% c('BSG','EMB'),Role:='transporter auxiliary'];cur[grepl('^SLC|^MPC',Symbol),Role:='transporter'];cur[grepl('^PDK|^PDP',Symbol),Role:='PDH regulator']
cur[,Compartment:='see frozen UniProt entry; localization not measured by RNA']
cur[,DirectionCaveat:='Expression does not measure reaction rate, net direction or metabolite concentration']
cur[Module=='Lactyl substrate formation',DirectionCaveat:='Noncanonical activity depends on localization and biochemical context; total RNA is not activity']
cur[Module=='D-lactate and glyoxalase',DirectionCaveat:='D-lactate branch; do not combine with L-lactate or presume MS isomer resolution']
cur[Symbol=='G6PC3',DirectionCaveat:='Ubiquitous glucose-6-phosphatase-related enzyme; not evidence of systemic glucose export']
cur[,InclusionReason:=paste('Mechanistic candidate:',Module)]
fwrite(cur,file.path(cfg,'curated_id_mapping.csv'))
# Preserve ambiguous ID mappings in catalog; no arbitrary first-match selection.
terms<-rbindlist(lapply(list.files(raw,pattern='_term.json$',full.names=TRUE),function(f){x<-fromJSON(f)$results;data.table(Term=x$id,Module=x$name,Definition=x$definition$text,Obsolete=x$isObsolete)}),fill=TRUE)
fwrite(terms,file.path(cfg,'go_terms.csv'))
ga<-merge(anns[ExcludedNOT==FALSE],map,by='UniProt',all.x=TRUE,allow.cartesian=TRUE)
ga<-merge(ga,terms[,.(Term,Module)],by='Term')
gmem<-unique(ga[!is.na(Ensembl),.(Ensembl,Symbol,Module,Term,Origin='GO annotation')])
cmen<-unique(cur[!is.na(Ensembl),.(Ensembl,Symbol,Module,Term='',Origin='curated mechanism')])
mem<-unique(rbind(gmem,cmen));fwrite(mem,file.path(cfg,'gene_module_membership.csv'))
fwrite(as.data.table(metadata(org.Hs.eg.db)),file.path(cfg,'annotation_database_metadata.csv'))
writeLines(capture.output(sessionInfo()),file.path(cfg,'catalog_sessionInfo.txt'))
cat('GO annotations',nrow(anns),'GO terms',nrow(terms),'union Ensembl genes',uniqueN(mem$Ensembl),'\n')
