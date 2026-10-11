suppressPackageStartupMessages({library(data.table);library(readxl)})
out<-'audit/20261011_all_groups_algorithm_audit/whole_proteome'; d<-as.data.table(read.csv(gzfile(file.path(out,'original_direct_source_membership.csv.gz')))); sources<-unique(d$SourceFile); logs<-list();hits<-list()
for(p in sources){
 if(grepl('xlsx$',p)){
 if(grepl('tables4',p)){ds<-list(Sheet1=read_excel(p,skip=5))} else if(grepl('patients_iTRAQ',p)){ds<-lapply(c('ANTs','CISs'),function(s)read_excel(p,sheet=s));names(ds)<-c('ANTs','CISs')} else if(grepl('Dataset1',p)){ds<-list(ESCC=read_excel(p,sheet='2.a protein raw information',skip=1))}else{ds<-list(Annotation_Combine=read_excel(p,sheet="Annotation_Combine"))}
 }else{ds<-list(table=fread(p,skip=if(grepl('HK-2_Spectronaut',p))2 else 0,integer64='double'))}
 for(sh in names(ds)){x<-as.data.table(ds[[sh]]);cols<-grep('protein|accession|majority|razor|reverse|contaminant|decoy',names(x),ignore.case=TRUE,value=TRUE)
 for(j in cols){v<-as.character(x[[j]]);h<-which(grepl('(^|[;,: ])(REV__|CON__|DECOY|CONTAM)',v,ignore.case=TRUE));logs[[length(logs)+1]]<-data.table(SourceFile=p,Sheet=sh,Field=j,Rows=nrow(x),PrefixedRows=length(h));if(length(h))for(i in h){hits[[length(hits)+1]]<-data.table(SourceFile=p,Sheet=sh,Row=i,Field=j,Value=v[i],Reverse=if('Reverse'%in%names(x))as.character(x$Reverse[i])else'',PotentialContaminant=if('Potential contaminant'%in%names(x))as.character(x[['Potential contaminant']][i])else'')}}
 }}
fwrite(rbindlist(logs),file.path(out,'all_direct_identifier_fields.csv'));if(length(hits))fwrite(rbindlist(hits),file.path(out,'all_direct_prefixed_rows.csv'))
