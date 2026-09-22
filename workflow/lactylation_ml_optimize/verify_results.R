library(data.table)
out <- 'outputs/20260922_br_v5'
p <- fread(file.path(out,'oos_predictions.csv.gz'))
s <- fread(file.path(out,'performance_summary.csv'))
tests <- fread('audit/20260922_br_v5/acceptance_tests.csv')
check <- function(n,ok){stopifnot(isTRUE(ok));tests <<- rbind(tests,data.table(Test=n,Pass=TRUE))}
for(m in s$Model){
 prob<-p[[paste0('p_',m)]];loss<- -p$y*log(prob)-(1-p$y)*log1p(-prob)
 independent<-mean(vapply(split(loss,p$ConnGroup),mean,numeric(1)))
 check(paste('base R grouped LogLoss',m),abs(independent-s[Model==m,LogLoss])<1e-12)
}
old<-fread('outputs/20260921_br_v4/BRv4_performance_summary.csv')[Summary=='equal_weight_19folds']
for(m in intersect(old$Model,s$Model))check(paste('unchanged v4 metric',m),max(abs(unlist(old[Model==m,.(LogLoss,AP,ROC,Brier)])-unlist(s[Model==m,.(LogLoss,AP,ROC,Brier)])))<1e-12)
z<-fread(file.path(out,'s2_per_group.csv'))
check('complete S2 nested retraining',nrow(z)==19*18*2 && all(z$Dropped!=z$Test) && all(is.finite(z$LogLoss)))
check('S2 unique drop-test-model combinations',!anyDuplicated(z[,.(Dropped,Test,Model)]))
r<-fread(file.path(out,'inner_group_losses.csv'));k<-fread(file.path(out,'inner_risks.csv'))
check('no outer group in inner validation',all(r$OuterGroup!=r$ConnGroup))
recalc<-r[,.(Recomputed=mean(LogLoss)),by=.(OuterGroup,Model,Candidate)]
j<-merge(k,recalc,by=c('OuterGroup','Model','Candidate'));check('all inner risks independently regrouped',nrow(j)==19*2*9 && max(abs(j$Risk-j$Recomputed))<1e-12)
fwrite(tests,'audit/20260922_br_v5/acceptance_tests.csv')
cat(nrow(tests),'checks PASS\n')
