#!/usr/bin/env Rscript
# 03b_BRef_train.R — 外层留连通分组（19折），折内标准化，M0-M3比较。
suppressPackageStartupMessages({library(data.table)})
set.seed(25)
W <- "/Users/gzy2520/Desktop/Research/kla/.worktrees/kla31-rnaseq-20260912"
setwd(W); out <- "outputs/20260921_lactylation_ml_review"
d <- fread(file.path(out,"BRef_features.csv.gz"))
d[, y := as.integer(y)]
feats <- list(M0=character(0), M1="target_expr", M2=c("target_expr","lac_mean"), M3=c("target_expr","lac_mean","reg_mean"))
ap <- function(y, p) { o <- order(-p); yy <- y[o]; P <- sum(yy); if (P==0||P==length(yy)) return(NA_real_); pr <- cumsum(yy)/seq_along(yy); rc <- cumsum(yy)/P; mean(pr[yy==1]) }
auc_mw <- function(y, p) { if (length(unique(y))<2) return(NA_real_); r <- rank(p); n1 <- sum(y==1); n0 <- sum(y==0); (sum(r[y==1]) - n1*(n1+1)/2)/(n1*n0) }
brier <- function(y, p) mean((p-y)^2)
conns <- sort(unique(d$ConnGroup))
preds <- list(); per <- list()
for (cg in conns) {
  te <- d[ConnGroup==cg]; tr <- d[ConnGroup!=cg]
  # 折内中心/尺度（训练估计）
  mus <- list(); sds <- list()
  for (f in unique(unlist(feats))) { mus[[f]] <- mean(tr[[f]]); sds[[f]] <- sd(tr[[f]]); if (!is.finite(sds[[f]])||sds[[f]]==0) sds[[f]] <- 1 }
  zs <- function(dd) { o <- dd; for (f in names(mus)) o[, (paste0("z_",f)) := (get(f)-mus[[f]])/sds[[f]]]; o }
  tr <- zs(copy(tr)); te <- zs(copy(te))
  p0 <- mean(tr$y)
  te[, p_M0 := p0]
  for (m in c("M1","M2","M3")) {
    fz <- paste0("z_", feats[[m]])
    fm <- as.formula(paste("y ~", paste(fz, collapse=" + ")))
    fit <- tryCatch(suppressWarnings(glm(fm, data=tr, family=binomial())), error=function(e) NULL)
    pv <- tryCatch({ if (is.null(fit)) rep(NA_real_, nrow(te)) else as.numeric(predict(fit, newdata=te, type="response")) }, error=function(e) rep(NA_real_, nrow(te)))
    if (any(!is.finite(pv))) pv[!is.finite(pv)] <- p0
    te[, (paste0("p_",m)) := pv]
    if (m!="M0" && fit$converged==FALSE) te[, note := "nonconverged"]  # 诊断列（若存在）
  }
  te[, ConnGroup_out := cg]
  preds[[as.character(cg)]] <- te[, .(GroupID, ReferenceKey, ConnGroup, ConnGroup_out, PXD, BaseAccession, EnsemblGeneID, y, p_M0, p_M1, p_M2, p_M3)]
}
P <- rbindlist(preds)
fwrite(P, file.path(out,"BRef_oos_preds.csv.gz"))
pooled <- rbindlist(lapply(c("M0","M1","M2","M3"), function(m) {
  v <- paste0("p_",m)
  data.table(Model=m, N=nrow(P), PosRate=round(mean(P$y),4),
    PRAUC=round(ap(P$y, P[[v]]),4), ROCAUC=round(auc_mw(P$y, P[[v]]),4), Brier=round(brier(P$y, P[[v]]),4))
}))
fwrite(pooled, file.path(out,"BRef_pooled_perf.csv"))
byfold <- rbindlist(lapply(conns, function(cg) rbindlist(lapply(c("M0","M1","M2","M3"), function(m) {
  s <- P[ConnGroup_out==cg]; v <- paste0("p_",m)
  data.table(ConnGroup=cg, Model=m, N=nrow(s), Pos=sum(s$y), PRAUC=round(ap(s$y,s[[v]]),4), ROCAUC=round(auc_mw(s$y,s[[v]]),4), Brier=round(brier(s$y,s[[v]]),4))
}))))
fwrite(byfold, file.path(out,"BRef_perfold_perf.csv"))
# 系数稳定性（M3每折系数）
cat("POOLED:\n"); print(pooled)
cat("PER-FOLD (head):\n"); print(head(byfold, 12))
