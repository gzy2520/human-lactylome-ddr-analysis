#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table);library(jsonlite)})
cfg<-'config/lactate_metabolism';args<-commandArgs(TRUE);out<-if(length(args))file.path(args[1],'tables') else 'outputs/20260920_lactate_metabolism/tables'
c<-fread(file.path(cfg,'curated_id_mapping.csv'));u<-fread(file.path(cfg,'uniprot_current_ensembl.csv'));f<-fread(file.path(cfg,'uniprot_function_evidence.csv'));mem<-fread(file.path(out,'membership_coverage.csv'))
measured<-unique(mem[Origin=='curated mechanism' & InCommonMatrix==TRUE,.(Ensembl,Symbol)])
x<-unique(c[,.(Ensembl,Symbol,Module,Role,Reaction,EvidenceSource,DirectionCaveat)])
u<-unique(u[,.(Ensembl,UniProt)]);f<-f[,.(UniProt,Compartment,Function,CatalyticReaction,LastAnnotationUpdate)]
x<-merge(x,merge(u,f,by='UniProt'),by='Ensembl',all.x=TRUE,allow.cartesian=TRUE)
x[,Measured:=Ensembl %in% measured$Ensembl]
x[,CurrentUniProtConfirmed:=!is.na(UniProt)]
stopifnot(all(x[Measured==TRUE,CurrentUniProtConfirmed]))
fwrite(x,file.path(cfg,'curated_mechanism_evidence.csv'))
a<-fread(file.path(cfg,'go_annotations.csv'));m<-fread(file.path(cfg,'go_id_mapping.csv'))
a<-merge(a,m,by='UniProt',all.x=TRUE,allow.cartesian=TRUE)
a[,Experimental:=Evidence %in% c('EXP','IDA','IPI','IMP','IGI','IEP','HTP','HDA','HMP','HGI','HEP')]
a[,CommonMatrix:=Ensembl %in% mem[InCommonMatrix==TRUE,Ensembl]]
fwrite(a,file.path(cfg,'go_annotations_mapped.csv'))
summary<-a[,.(AnnotationRecords=uniqueN(AnnotationID),ProteinProducts=uniqueN(UniProt),UnmappedProducts=uniqueN(UniProt[is.na(Ensembl)|Ensembl=='']),DirectRecords=uniqueN(AnnotationID[DirectAnnotation==TRUE]),ExperimentalRecords=uniqueN(AnnotationID[Experimental==TRUE]),ExperimentalMeasuredIDs=uniqueN(Ensembl[Experimental==TRUE & CommonMatrix==TRUE & ExcludedNOT==FALSE])),by=Term]
fwrite(summary,file.path(out,'go_evidence_coverage.csv'))
# Check download completion and annotation relation contract against cached term descendants.
checks<-rbindlist(lapply(summary$Term,function(term){
 key<-gsub(':','_',term);fs<-list.files(file.path(cfg,'raw'),pattern=paste0('^',key,'_annotations_'),full.names=TRUE)
 pages<-lapply(fs,fromJSON);total<-pages[[1]]$numberOfHits
 desc<-fromJSON(file.path(cfg,'raw',paste0(key,'_descendants.json')))$results
 allowed<-unique(c(term,unlist(desc$descendants)))
 data.table(Term=term,ExpectedAnnotations=total,Downloaded=sum(vapply(pages,function(z)nrow(z$results),integer(1))),Pages=length(fs),ExpectedPages=pages[[1]]$pageInfo$total,RelationScopeValid=all(a[Term==term,AnnotatedGO] %in% allowed))
}));stopifnot(all(checks$ExpectedAnnotations==checks$Downloaded),all(checks$Pages==checks$ExpectedPages),all(checks$RelationScopeValid));fwrite(checks,file.path(out,'go_download_validation.csv'))
cat('CATALOG_PASS: all GO pages complete; all measured curated genes confirmed by current reviewed UniProt cross-reference\n')
