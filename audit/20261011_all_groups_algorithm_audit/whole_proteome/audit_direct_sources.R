suppressPackageStartupMessages({library(data.table);library(readxl)})
project_root<-getwd();source_root<-'/Users/gzy2520/Desktop/Research/kla/data'
out<-'audit/20261011_all_groups_algorithm_audit/whole_proteome'
txt<-readLines('R/candidate/prepare_sample_boxplot_inputs.R')
expr<-parse(text=txt)
# Load only pure helper definitions; never execute original top-level pipeline.
for(e in expr)if(is.call(e)&&as.character(e[[1]])=='<-'&&is.call(e[[3]])&&as.character(e[[3]][[1]])=='function')eval(e)
original_split<-split_accessions;original_records<-records
diagnostics<-list();parse_changes<-list();mode<-'original'
strict_split<-function(v){
 token<-trimws(unlist(strsplit(as.character(v),'[;,:]')));token<-token[!is.na(token)]
 token<-token[!grepl('^(REV__|CON__|DECOY|CONTAM)',token,ignore.case=TRUE)]
 token<-sub('^(sp|tr)\\|','',token);token<-sub('\\|.*$','',token)
 token<-sub('^NX_','',token);token<-sub('-[0-9]+$','',token)
 sort(unique(token[is_uniprot(token)]))
}
split_accessions<-function(values){
 old<-original_split(values);new<-strict_split(values)
 result<-if(mode=='original')old else new
 attr(result,'audit')<-list(InputValues=length(values),ColonValues=sum(grepl(':',values),na.rm=TRUE),PipeValues=sum(grepl('\\|',values),na.rm=TRUE),MultipleGroupValues=sum(grepl('[;,:]',values),na.rm=TRUE),ContaminantTokens=sum(grepl('REV__|CON__|DECOY|CONTAM',values,ignore.case=TRUE),na.rm=TRUE),OldCount=length(old),CorrectedCount=length(new),Added=paste(setdiff(new,old),collapse=';'),Removed=paste(setdiff(old,new),collapse=';'),TenCharCount=sum(nchar(new)==10))
 result
}
records<-function(pxd,group,sample_id,accessions,source_mode,source_file,sample_class,observation_type='sample'){
 d<-attr(accessions,'audit');if(!is.null(d))diagnostics[[length(diagnostics)+1]]<<-c(list(Mode=mode,PXD=pxd,SampleGroup=group,SampleID=sample_id,SourceFile=source_file),d)
 original_records(pxd,group,sample_id,as.character(accessions),source_mode,source_file,sample_class,observation_type)
}
reference_records_list<-list();add_reference_records<-function(value,reference_pxd){value[,ReferencePXD:=reference_pxd];reference_records_list[[length(reference_records_list)+1]]<<-value}
lo<-grep('^# Same-study, sample-resolved references',txt);hi<-grep('^direct_reference_keys <-',txt)-1L
sources<-parse(text=txt[lo:hi])
all<-list()
for(mm in c('original','corrected')){
 mode<-mm;reference_records_list<-list();eval(sources)
 all[[mm]]<-rbindlist(reference_records_list,fill=TRUE)
}
# The ESCC reference was added by a separate renderer; read its original workbook directly.
p<-file.path(source_root,'PXD065830/supplementary/Dataset1.xlsx')
x<-read_excel(p,sheet='2.a protein raw information',col_names=FALSE,.name_repair='minimal');headers<-trimws(as.character(unlist(x[2,])));ri<-3:nrow(x);cols<-which(grepl('^ESCC-[0-9]+T$',headers));stopifnot(length(cols)==94)
for(mm in c('original','corrected')){mode<-mm;v<-rbindlist(lapply(cols,function(j)records('PXD064038','MEC and NEC ESCC groups',headers[j],split_accessions(x[[1]][ri][positive_numeric(x[[j]][ri])]),'original_author_workbook',p,'tumor_sample')));all[[mm]]<-rbindlist(list(all[[mm]],v),fill=TRUE)}
go<-fread('data/publication_input/human_ddr_go_annotations.tsv');ddr<-strict_split(go[(is.na(get('TAXON ID'))|get('TAXON ID')==9606)&!grepl('(^|\\|)NOT($|\\|)',QUALIFIER)&get('GENE PRODUCT DB')=='UniProtKB'][['GENE PRODUCT ID']])
for(mm in names(all)){z<-unique(all[[mm]]);z[,`:=`(IsDdr=BaseAccession %in% ddr,Dataset='Whole proteome')];fwrite(z,file.path(out,paste0(mm,'_direct_source_membership.csv.gz')));s<-z[,.(ProteinCount=uniqueN(BaseAccession),DdrProteinCount=uniqueN(BaseAccession[IsDdr]),SourceFile=first(SourceFile),SampleClass=first(SampleClass)),by=.(PXD,SampleGroup,SampleID)];s[,Fraction:=100*DdrProteinCount/ProteinCount];fwrite(s,file.path(out,paste0(mm,'_direct_source_scores.csv')))}
fwrite(rbindlist(diagnostics,fill=TRUE),file.path(out,'direct_source_parser_evidence.csv'))
# Preserve cell masks from complete-channel sources, including source sheet and column identity.
masks<-list()
inspect<-function(d,path,sheet,cols,idcol){for(col in cols){v<-suppressWarnings(as.numeric(d[[col]]));masks[[length(masks)+1]]<<-data.table(SourceFile=path,Sheet=sheet,Column=col,Rows=nrow(d),Positive=sum(is.finite(v)&v>0),Zero=sum(is.finite(v)&v==0),Negative=sum(is.finite(v)&v<0),Nonfinite=sum(!is.finite(v)),UniqueNumeric=uniqueN(v[is.finite(v)]),Minimum=min(v,na.rm=TRUE),Maximum=max(v,na.rm=TRUE))}}
p<-file.path(source_root,'PXD033146/search_results/extracted_pairing/search_result-HA119TQ/proteinGroups.txt');d<-fread(p,integer64="double");d<-d[valid_maxquant_rows(d)];inspect(d,p,'proteinGroups valid rows',paste('Reporter intensity corrected',1:3),'Majority protein IDs')
p<-file.path(source_root,'PXD050470/supplementary/prca2331-sup-0006-tables4.xlsx');d<-read_excel(p,sheet='Sheet1',skip=5);inspect(d,p,'Sheet1;skip=5',paste0('Intensity_',c('H072','H081','H187')),'Protein accession')
p<-file.path(source_root,'PXD065775/search_results/20170330_01-24_patients_iTRAQ.xlsx');for(sh in c('ANTs','CISs')){d<-read_excel(p,sheet=sh);inspect(d,p,sh,c(paste0('Non-rec',1:4),paste0('Rec',1:4)),'Accession')}
fwrite(rbindlist(masks),file.path(out,'complete_channel_source_values_audit.csv'))
cat('Read original tables for',uniqueN(all$corrected[,.(PXD,SampleGroup)]),'groups; corrected direct units',nrow(unique(all$corrected[,.(PXD,SampleGroup,SampleID)])),'\n')
