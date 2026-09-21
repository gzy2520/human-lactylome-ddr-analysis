# 09_sensitivity_v2.R — 两类敏感性 + 材料级校准 + 权重检查 + 系数稳定性
suppressPackageStartupMessages({library(data.table)})
source("workflow/lactylation_ml_review/05_RR_functions.R")
set.seed(25)
W <- "/Users/gzy2520/Desktop/Research/kla/.worktrees/kla31-rnaseq-20260912"
setwd(W); out <- "outputs/20260922_br_v2"
P <- fread(file.path(out, "BRv2_oos_preds.csv.gz"))
MODELS <- c("M0","MP","MP_only","ME","ME2","ME2P","XEN")
conns <- sort(unique(P$ConnGroup))

## ---- S1: 评价汇总对单个测试组的敏感性（删一测试折后重算等权均值；不证明训练不受影响）----
s1 <- rbindlist(lapply(conns, function(cg) {
  sub <- P[ConnGroup_out != cg]
  rbindlist(lapply(MODELS, function(m) {
    per <- sub[, metrics_row(get("y"), get(paste0("p_", m)), m, ConnGroup)]
    data.table(dropped_test = cg, Model = m,
      eqw_AP = round(mean(per$AP, na.rm=TRUE),4), eqw_LL = round(mean(per$LogLoss, na.rm=TRUE),4))
  }))
}))
s1_base <- s1[dropped_test == conns[length(conns)]][, .(Model, eqw_AP_ref18 = eqw_AP)]
s1 <- merge(s1, s1_base, by = "Model")
fwrite(s1, file.path(out, "S1_eval_sensitivity_drop_one_test_group.csv"))
cat("S1 MP eqw_AP range:", min(s1[Model=="MP", eqw_AP]), "-", max(s1[Model=="MP", eqw_AP]), "\n")

## ---- S2: 训练影响（删一连通组后重新执行分组验证，仅在剩余18组上与原结果对比）----
d <- fread(file.path(out, "BRef_features_v2.csv.gz")); d[, y := as.integer(y)]
mod_name <- "mact_generation_lactate"; LAMS <- c(0.05, 0.1, 0.3, 1.0)
spec_me2p <- c("z_target_expr", "z_mod", "z_mp")
s2 <- rbindlist(lapply(conns, function(drop) {
  dsub <- d[ConnGroup != drop]
  rbindlist(lapply(setdiff(conns, drop), function(cg) {
    te <- copy(dsub[ConnGroup == cg]); tr <- copy(dsub[ConnGroup != cg])
    refs_tr <- unique(tr$ReferenceKey)
    tv <- tapply(tr$target_expr, tr$ReferenceKey, median)[refs_tr]
    mu_t <- mean(tv, na.rm=TRUE); sd_t <- sd(tv, na.rm=TRUE); if (!is.finite(sd_t)||sd_t==0) sd_t <- 1
    imp_t <- median(tv, na.rm=TRUE); zt <- function(dd){x<-dd$target_expr;x[is.na(x)]<-imp_t;(x-mu_t)/sd_t}
    mv <- tapply(tr[[mod_name]], tr$ReferenceKey, median)[refs_tr]
    mu_m <- mean(mv, na.rm=TRUE); sd_m <- sd(mv, na.rm=TRUE); if (!is.finite(sd_m)||sd_m==0) sd_m <- 1
    imp_m <- median(mv, na.rm=TRUE); zm <- function(dd){x<-dd[[mod_name]];x[is.na(x)]<-imp_m;(x-mu_m)/sd_m}
    st <- tr[, .(n=.N, k=sum(y)), by=BaseAccession]; st[, freq_full := (k+1)/(n+2)]
    te[, freq_use := st[match(te$BaseAccession, st$BaseAccession), freq_full]]
    te[is.na(freq_use), freq_use := (sum(tr$y)+1)/(nrow(tr)+2)]
    frf <- st[match(tr$BaseAccession, st$BaseAccession), freq_full]
    frf[is.na(frf)] <- (sum(tr$y)+1)/(nrow(tr)+2)
    fv <- tapply(frf, tr$ReferenceKey, median)[refs_tr]
    mu_f <- mean(fv, na.rm=TRUE); sd_f <- sd(fv, na.rm=TRUE); if (!is.finite(sd_f)||sd_f==0) sd_f <- 1
    inner_conns <- sort(unique(tr$ConnGroup)); set.seed(25)
    ia <- split(sample(inner_conns), cut(seq_along(inner_conns), 5))
    tr[, inner_f := sapply(ConnGroup, function(g) which(vapply(ia, function(s) g %in% s, logical(1))))]
    tr[, freq_use := NA_real_]
    for (f in seq_len(5)) {
      sb <- tr[inner_f != f]; stt <- sb[, .(n=.N, k=sum(y)), by=BaseAccession]
      ov <- (sum(sb$y)+1)/(nrow(sb)+2)
      fr <- stt[match(tr[inner_f == f, BaseAccession], stt$BaseAccession), (k+1)/(n+2)]
      fr[is.na(fr)] <- ov
      tr[inner_f == f, freq_use := fr]
    }
    zf_tr <- function(dd) (dd$freq_use - mu_f)/sd_f
    zf_te <- function(dd) (dd$freq_use - mu_f)/sd_f
    build_X <- function(dd, spec, zf) {
      X <- matrix(1, nrow(dd), 1)
      for (f in spec) X <- cbind(X, switch(f, z_target_expr=zt(dd), z_mod=zm(dd), z_mp=zf(dd)))
      if (ncol(X)>1) colnames(X)[2:ncol(X)] <- spec
      X
    }
    te[, p_M0 := mean(tr$y)]
    te[, p_MP := freq_use]
    fitE <- ridge_glm_v2(build_X(tr, "z_target_expr", zf_tr), tr$y, LAMS[2])
    te[, p_ME := pred_ridge_v2(fitE, build_X(te, "z_target_expr", zf_te))]
    fit2 <- ridge_glm_v2(build_X(tr, spec_me2p, zf_tr), tr$y, LAMS[2])
    te[, p_ME2P := pred_ridge_v2(fit2, build_X(te, spec_me2p, zf_te))]
    rbindlist(lapply(c("M0","MP","ME","ME2P"), function(m) {
      te[, metrics_row(get("y"), get(paste0("p_", m)), m, cg)][, dropped_train := drop][]
    }))
  }))
}))
fwrite(s2, file.path(out, "S2_train_sensitivity_drop_one_conn.csv"))
s2_sum <- s2[, .(eqw_AP = round(mean(AP, na.rm=TRUE),4), eqw_LL = round(mean(LogLoss, na.rm=TRUE),4)), by=.(dropped_train, Model)]
fwrite(s2_sum, file.path(out, "S2_train_sensitivity_summary.csv"))
cat("S2 done rows:", nrow(s2), "\n")

## ---- 蛋白行数多的材料是否主导：材料行数 vs 正例率 ----
mat_w <- P[, .(n_rows=.N, pos_rate=round(mean(y),4)), by=ReferenceKey]
fwrite(mat_w[order(-n_rows)], file.path(out, "material_row_weights_check.csv"))

## ---- 材料级校准：预测均值 vs 实际检出率（材料级模块的主要评估口径）----
matcal <- rbindlist(lapply(c("M0","MP","ME","ME2","ME2P"), function(m) {
  P[, .(Model=m, ReferenceKey, mean_pred=mean(get(paste0("p_",m))), obs_rate=mean(y), n_rows=.N)]
}))
matcal[, calib_err := round(mean_pred - obs_rate, 4)]
fwrite(matcal[order(Model, ReferenceKey)], file.path(out, "BRv2_material_level_calibration.csv"))
cat("Material-level calib err (mean |err|):\n")
print(matcal[, .(mean_abs_err=round(mean(abs(calib_err)),4), max_abs_err=round(max(abs(calib_err)),4)), by=Model])

## ---- 系数稳定性（z_mod 跨折符号/大小）----
coefs <- fread(file.path(out, "BRv2_coefficients.csv"))
zs <- coefs[Feature %in% c("z_mod","z_target_expr","z_mp")]
fwrite(zs, file.path(out, "BRv2_module_coef_stability.csv"))
cat("z_mod coef sign (ME2):\n"); print(zs[Model=="ME2" & Feature=="z_mod", .(pos=sum(Coef>0), neg=sum(Coef<0), sd=round(sd(Coef),4))])
cat("S1/S2/calib/coef done\n")
