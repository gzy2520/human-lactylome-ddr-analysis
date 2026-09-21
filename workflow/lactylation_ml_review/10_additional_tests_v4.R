library(data.table)
O<-'outputs/20260921_br_v4'
res<-list();test<-function(id,ok,details){stopifnot(isTRUE(ok));res[[id]]<<-data.table(TestID=id,Status='PASS',Details=details)}
loss<-fread(file.path(O,'inner_connected_group_losses.csv'))
re<-loss[,.(loss=mean(LogLoss),groups=.N),by=.(OuterConnGroup,Model,candidate_lambda)]
tune<-fread(file.path(O,'BRv4_inner_tuning_log.csv'));setnames(tune,'ConnGroup','OuterConnGroup')
z<-merge(re,tune,by=c('OuterConnGroup','Model','candidate_lambda'))
test('B1_inner_group_mean',nrow(z)==19*6*4&&all(z$groups==18)&&max(abs(z$loss-z$inner_mean_logloss))<1e-12,'Each candidate risk reproduced from 18 actual connected-group losses.')
test('B2_inner_no_outer',all(loss$OuterConnGroup!=loss$ConnGroup),'Outer test groups absent from every inner risk table.')
g<-fread(file.path(O,'module_source_values.csv'))
ids<-c('ENSG00000134333','ENSG00000111716','ENSG00000166796','ENSG00000171989')
l<-g[Ensembl%in%ids,.(module=mean(MedianLog2TPM),genes=uniqueN(Ensembl)),by=ReferenceKey]
d<-fread(file.path(O,'BRef_features_v4.csv.gz'))
test('B3_four_gene_L_module',all(l$genes==4)&&max(abs(d$mact_generation_lactate-l$module[match(d$ReferenceKey,l$ReferenceKey)]))<1e-12,'Exact four-gene means; LDHD absent from L feature.')
test('B4_alternate_LDHA_missing',all(is.na(g[Ensembl=='ENSG00000288299',MedianLog2TPM])),'Alternate LDHA mapping missing in 28 sources, not counted as measured.')
check<-fread(file.path(O,'target_source_verification.csv.gz'))
test('B5_source_target',nrow(check)==6744&&!anyNA(check$difference)&&max(abs(check$difference))<1e-10,'All target features independently traced to source matrices.')
# Match each old/new comparison by connected-group ID, never row position.
old<-fread('audit/20260921_br_v4/v3_module_fold_deltas.csv')
test('B6_old_worse_count',nrow(old)==19&&sum(old$delta>1e-12)==9,'Old 15/19 statement reproduced as 9/19 from predictions.')
s<-fread(file.path(O,'S2_completion_checklist.csv'))
test('B7_S2_scope',nrow(s)==19&&all(s$models_evaluated=='MP;ME_P;ME2P')&&all(s$status=='COMPLETED'),'Exactly declared three models completed for all deletions.')
fwrite(rbindlist(res),file.path(O,'additional_acceptance.csv'))
print(rbindlist(res))
