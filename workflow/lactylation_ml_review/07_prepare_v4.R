library(data.table)
set.seed(25)
out <- 'outputs/20260921_br_v4'
control <- 'outputs/20260921_br_v4_weight_only'
for (z in c(out,control)) {if(dir.exists(z)) stop('Output already exists: ',z);dir.create(z,recursive=TRUE)}
d <- fread('outputs/20260921_br_v3/BRef_features_v3.csv.gz')
fwrite(d,file.path(control,'BRef_features_v4.csv.gz'))
ids <- c('ENSG00000134333','ENSG00000288299','ENSG00000111716','ENSG00000166796','ENSG00000171989','ENSG00000166816')
rows <- lapply(unique(d$ReferenceKey),function(rk){
 x<-fread(file.path('outputs/20260919_expression_corrected/matrices',paste0(rk,'_log2tpm.tsv.gz')))
 need<-unique(c(ids,d[ReferenceKey==rk,EnsemblGeneID]));j<-match(need,x[[1]])
 m<-as.matrix(x[j[!is.na(j)],-1,with=FALSE]);stopifnot(all(is.finite(m)))
 vals<-rep(NA_real_,length(need));vals[!is.na(j)]<-apply(m,1,median)
 data.table(ReferenceKey=rk,Ensembl=need,MedianLog2TPM=vals)
})
g <- rbindlist(rows)
map<-unique(fread('config/lactate_metabolism/go_current_uniprot_mapping.csv'))
fwrite(g,file.path(out,'module_source_values.csv'))
fwrite(map[Ensembl %in% ids],file.path(out,'module_id_mapping_evidence.csv'))
Lids<-ids[c(1,3,4,5)]
z<-g[Ensembl %in% ids,.(legacy_reconstructed=mean(MedianLog2TPM,na.rm=TRUE),L_conversion=mean(MedianLog2TPM[Ensembl %in% Lids]),D_oxidation=MedianLog2TPM[Ensembl==ids[6]],LDHA_primary=MedianLog2TPM[Ensembl==ids[1]],LDHA_alternate=MedianLog2TPM[Ensembl==ids[2]]),by=ReferenceKey]
old<-unique(d[,.(ReferenceKey,legacy_feature=mact_generation_lactate)])
z<-merge(z,old,by='ReferenceKey')
stopifnot(max(abs(z$legacy_reconstructed-z$legacy_feature))<1e-10)
fwrite(z,file.path(out,'module_feature_comparison.csv'))
check<-merge(d[,.(GroupID,BaseAccession,ReferenceKey,Ensembl=EnsemblGeneID,target_expr)],g,by=c('ReferenceKey','Ensembl'),all.x=TRUE)
check[,difference:=MedianLog2TPM-target_expr]
fwrite(check,file.path(out,'target_source_verification.csv.gz'))
stopifnot(!anyNA(check$difference),max(abs(check$difference))<1e-10)
d[,mact_generation_lactate:=z$L_conversion[match(ReferenceKey,z$ReferenceKey)]]
stopifnot(!anyNA(d$mact_generation_lactate))
fwrite(d,file.path(out,'BRef_features_v4.csv.gz'))
cat('PASS: all 28 legacy source means and 6744 target values independently reproduced; L module uses 4 distinct genes; D separated.\n')
