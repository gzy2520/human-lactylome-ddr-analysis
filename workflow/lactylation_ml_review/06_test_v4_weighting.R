library(data.table)
source('workflow/lactylation_ml_review/05e_training_engine_v4.R')
# Regression fixture: row risk and equal connected-group risk must differ.
d<-data.table(g=c(rep(1,100),2),y=c(rep(0,100),1),p=.1)
d[,loss:=-y*log(p)-(1-y)*log1p(-p)]
stopifnot(abs(mean(d[,.(v=mean(loss)),by=g]$v)-mean(d$loss))>.5)
# Exercise actual engine, comparing stored group losses with tuning means.
x<-fread('outputs/20260921_br_v3/BRef_features_v3.csv.gz')
a<-train_predict_cv_engine_v4(x[ConnGroup!=1],x[ConnGroup==1],models='ME2P')
g<-a$inner_group_losses
z<-g[,.(computed=mean(LogLoss),n=.N),by=.(Model,candidate_lambda)]
z<-merge(z,a$inner_tuning,by=c('Model','candidate_lambda'))
stopifnot(all(z$n==18),max(abs(z$computed-z$inner_mean_logloss))<1e-12)
cat('PASS: actual inner selection equals 18 connected-group losses, independent of row and inner-fold sizes.\n')

p<-fread('outputs/20260921_br_v3/BRv3_oos_preds.csv.gz')
delta<-p[,.(delta=logloss(y,p_ME2P)-logloss(y,p_ME_P)),by=ConnGroup]
stopifnot(sum(delta$delta>1e-12)==9)
dir.create('audit/20260921_br_v4',recursive=TRUE,showWarnings=FALSE)
fwrite(delta,'audit/20260921_br_v4/v3_module_fold_deltas.csv')
