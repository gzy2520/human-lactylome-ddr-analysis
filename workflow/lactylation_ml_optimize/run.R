source('workflow/lactylation_ml_optimize/engine.R')
setDTthreads(1)
out <- 'outputs/20260922_br_v5'
dir.create(out, recursive=TRUE, showWarnings=FALSE)
stopifnot(!file.exists(file.path(out,'oos_predictions.csv.gz')))
d <- fread('outputs/20260921_br_v4/BRef_features_v4.csv.gz')
old <- fread('outputs/20260921_br_v4/BRv4_oos_preds.csv.gz')
keys <- c('GroupID','ReferenceKey','ConnGroup','BaseAccession','EnsemblGeneID')
stopifnot(nrow(d)==6744,sum(d$y)==2708,uniqueN(d$ConnGroup)==19,!anyDuplicated(d[,..keys]))
preds <- choices <- inner <- risks <- list()
for(g in sort(unique(d$ConnGroup))){
  cat('Outer group',g,'\n');flush.console()
  f <- v5_fit(d[ConnGroup!=g],d[ConnGroup==g])
  preds[[as.character(g)]] <- f$pred
  for(n in c('selected','inner','risks')) f[[n]][,OuterGroup:=g]
  choices[[as.character(g)]] <- f$selected
  inner[[as.character(g)]] <- f$inner
  risks[[as.character(g)]] <- f$risks
}
p <- merge(old,rbindlist(preds),by=keys,all=TRUE)
stopifnot(nrow(p)==6744,all(is.finite(p$p_SHRINK)),max(abs(p$p_MP-p$p_MP_check))<1e-12)
fwrite(p,file.path(out,'oos_predictions.csv.gz'))
fwrite(rbindlist(choices),file.path(out,'selected_parameters.csv'))
fwrite(rbindlist(inner),file.path(out,'inner_group_losses.csv'))
fwrite(rbindlist(risks),file.path(out,'inner_risks.csv'))
fwrite(v5_grid(),file.path(out,'candidate_grid.csv'))
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
