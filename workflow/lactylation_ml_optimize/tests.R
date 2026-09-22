source('workflow/lactylation_ml_optimize/engine.R')
setDTthreads(1)
dir.create('audit/20260922_br_v5',recursive=TRUE,showWarnings=FALSE)
records <- list()
check <- function(name,ok){records[[name]] <<- data.table(Test=name,Pass=isTRUE(ok));stopifnot(isTRUE(ok));cat('PASS',name,'\n')}
d <- fread('outputs/20260921_br_v4/BRef_features_v4.csv.gz')
cfg <- v5_grid()[ID=='group_2']
groups <- sort(unique(d$ConnGroup));a <- d[!ConnGroup%in%groups[1:3]];v <- d[ConnGroup%in%groups[1:3]]
q <- v5_freq(a,v,cfg)$p
check('group evidence invariant to duplicate source rows',identical(q,v5_freq(rbind(a,a[ConnGroup==groups[4]]),v,cfg)$p))
check('unseen uses training equal group prevalence',abs(v5_freq(a,data.table(BaseAccession='UNSEEN'),cfg)$p-a[,.(r=mean(y)),by=ConnGroup][,mean(r)])<1e-12)
b <- copy(a);held <- v5_folds(a$ConnGroup,4)[[1]];b[ConnGroup%in%held,y:=1-y]
check('crossfit predictions exclude own held group labels',identical(v5_cf(a,cfg)[a$ConnGroup%in%held],v5_cf(b,cfg)[b$ConnGroup%in%held]))
check('calibration row weights balance source groups',max(abs(tapply(v5_weights(a),a$ConnGroup,sum)-1/uniqueN(a$ConnGroup)))<1e-12)
fit <- v5_cal(a,v5_cf(a,cfg));check('calibration stationary solution',fit$gradient<1e-5)
f <- v5_fit(a,v);changed <- copy(v);changed[,y:=1-y];f2 <- v5_fit(a,changed)
check('outer labels cannot change fitted choices or predictions',identical(f,f2))
check('inner groups exclude outer groups',!any(f$inner$ConnGroup%in%v$ConnGroup))
check('inner loss counts each training group once per candidate',all(f$risks$Groups==uniqueN(a$ConnGroup)))
if(file.exists('outputs/20260922_br_v5/oos_predictions.csv.gz')){
 p <- fread('outputs/20260922_br_v5/oos_predictions.csv.gz')
 check('complete outer coverage',nrow(p)==6744 && uniqueN(p$ConnGroup)==19 && sum(p$y)==2708)
 check('original MP exactly reproduced',max(abs(p$p_MP-p$p_MP_check))<1e-12)
 check('probabilities finite and bounded',all(is.finite(p$p_SHRINK_CAL)) && all(p$p_SHRINK_CAL>0 & p$p_SHRINK_CAL<1))
}
fwrite(rbindlist(records),'audit/20260922_br_v5/acceptance_tests.csv')
