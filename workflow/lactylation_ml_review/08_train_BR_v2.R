# 08_train_BR_v2.R — B-R v2 外层留连通组（19折），冻结配置见 config/lactylation_ml_review/model_config_v2.md
# 主指标 logloss；同报 Brier/AP/ROC。材料级模块主要影响材料级概率（加法逻辑回归），主指标不设ROC。
suppressPackageStartupMessages({library(data.table)})
source("workflow/lactylation_ml_review/05_RR_functions.R")
set.seed(25)
W <- "/Users/gzy2520/Desktop/Research/kla/.worktrees/kla31-rnaseq-20260912"
setwd(W); out <- "outputs/20260922_br_v2"
d <- fread(file.path(out, "BRef_features_v2.csv.gz"))
d[, y := as.integer(y)]
mods <- fread("config/lactylation_ml_review/modules_frozen_v2.csv")
mod_names <- unique(mods$ModuleID)
mod_name <- "mact_generation_lactate"  # ME2/ME2P/XEN 使用（冻结配置）
conns <- sort(unique(d$ConnGroup))
LAMS <- c(0.05, 0.1, 0.3, 1.0)  # 预固定候选
MODELS <- c("M0","MP","MP_only","ME","ME2","ME2P","XEN")
featspec <- list(
  M0 = character(0), MP = character(0), MP_only = character(0),
  ME = "z_target_expr",
  ME2 = c("z_target_expr", "z_mod"),
  ME2P = c("z_target_expr", "z_mod", "z_mp"),
  XEN = c("z_target_expr", "z_mod", "z_ix")
)
build_X <- function(dd, spec, zt, zm, zf) {
  X <- matrix(1, nrow = nrow(dd), ncol = 1, dimnames = list(NULL, "(Intercept)"))
  for (f in spec) {
    v <- switch(f, z_target_expr = zt(dd), z_mod = zm(dd), z_mp = zf(dd),
      z_ix = zt(dd) * zm(dd), stop("unknown feat ", f))
    X <- cbind(X, v)
  }
  if (ncol(X) > 1) colnames(X)[2:ncol(X)] <- spec
  X
}
preds_list <- list(); diag_list <- list(); coef_list <- list(); scal_list <- list(); mp_cov_list <- list()
for (cg in conns) {
  te <- copy(d[ConnGroup == cg]); tr <- copy(d[ConnGroup != cg])
  ## 标准化统计：仅训练折不重复参考材料（材料权重1，非按行重复加权）
  refs_tr <- unique(tr$ReferenceKey)
  med_by_ref <- function(col) tapply(tr[[col]], tr$ReferenceKey, median)[refs_tr]
  tv <- med_by_ref("target_expr"); mu_t <- mean(tv, na.rm=TRUE); sd_t <- sd(tv, na.rm=TRUE)
  if (!is.finite(sd_t) || sd_t == 0) sd_t <- 1
  imp_t <- median(tv, na.rm=TRUE)
  zt <- function(dd) { x <- dd$target_expr; x[is.na(x)] <- imp_t; (x - mu_t)/sd_t }
  mv <- med_by_ref(mod_name); mu_m <- mean(mv, na.rm=TRUE); sd_m <- sd(mv, na.rm=TRUE)
  if (!is.finite(sd_m) || sd_m == 0) sd_m <- 1
  imp_m <- median(mv, na.rm=TRUE)
  zm <- function(dd) { x <- dd[[mod_name]]; x[is.na(x)] <- imp_m; (x - mu_m)/sd_m }
  ## MP：测试行用全部外层训练频率；训练行用内层5折连通组交叉拟合（行自身标签不进自身特征）
  te[, freq_use := 0]  # placeholder, replaced below
  st_full <- tr[, .(n = .N, k = sum(y)), by = BaseAccession]
  st_full[, freq_full := (k + 1)/(n + 2)]
  te[, freq_use := st_full[match(te$BaseAccession, st_full$BaseAccession), freq_full]]
  te[is.na(freq_use), freq_use := (sum(tr$y) + 1)/(nrow(tr) + 2)]  # 回退=平滑总体率，已记录
  te[, mp_source := fifelse(is.na(st_full[match(te$BaseAccession, st_full$BaseAccession), freq_full]), "fallback_overall", "full_train")]
  freq_row_full <- st_full[match(tr$BaseAccession, st_full$BaseAccession), freq_full]
  freq_row_full[is.na(freq_row_full)] <- (sum(tr$y) + 1)/(nrow(tr) + 2)
  fv_refs <- tapply(freq_row_full, tr$ReferenceKey, median)[refs_tr]
  mu_f <- mean(fv_refs, na.rm = TRUE); sd_f <- sd(fv_refs, na.rm = TRUE)
  if (!is.finite(sd_f) || sd_f == 0) sd_f <- 1
  zf_te <- function(dd) (dd$freq_use - mu_f)/sd_f
  mp_cov_list[[as.character(cg)]] <- data.table(ConnGroup = cg,
    fallback_rate_test = mean(te$mp_source == "fallback_overall"),
    n_train_prot = uniqueN(tr$BaseAccession), n_test_prot = uniqueN(te$BaseAccession))
  ## 内层：训练折连通组切5组（共享来源整组不拆）
  inner_conns <- sort(unique(tr$ConnGroup))
  set.seed(25)
  inner_assign <- split(sample(inner_conns), cut(seq_along(inner_conns), 5))
  tr[, inner_f := sapply(ConnGroup, function(g) which(vapply(inner_assign, function(s) g %in% s, logical(1))))]
  tr[, freq_use := NA_real_]
  for (f in seq_len(5)) {
    sub_tr <- tr[inner_f != f]
    st <- sub_tr[, .(n = .N, k = sum(y)), by = BaseAccession]
    ov <- (sum(sub_tr$y) + 1)/(nrow(sub_tr) + 2)
    fr <- st[match(tr[inner_f == f, BaseAccession], st$BaseAccession), (k + 1)/(n + 2)]
    fr[is.na(fr)] <- ov
    tr[inner_f == f, freq_use := fr]
  }
  zf_tr <- function(dd) (dd$freq_use - mu_f)/sd_f
  ## 内层调参（ME2P/XEN；ME/ME2 预先声明固定 lambda=LAMS[2]）
  pick_lambda <- function(spec) {
    if (length(spec) <= 2) return(LAMS[2])
    ics <- sapply(LAMS, function(l) {
      ll <- 0
      for (f in seq_len(5)) {
        tr_i <- tr[inner_f != f]; ev_i <- tr[inner_f == f]
        Xi <- build_X(tr_i, spec, zt, zm, zf_tr); Xi_e <- build_X(ev_i, spec, zt, zm, zf_tr)
        fit <- ridge_glm_v2(Xi, tr_i$y, l)
        ll <- ll + logloss(ev_i$y, pred_ridge_v2(fit, Xi_e))
      }
      ll
    })
    names(ics) <- LAMS
    as.numeric(names(which.min(ics)))
  }
  te[, p_M0 := mean(tr$y)]
  te[, p_MP := freq_use]
  te[, p_MP_only := freq_use]
  for (m in c("ME","ME2","ME2P","XEN")) {
    spec <- featspec[[m]]
    lam <- pick_lambda(spec)
    Xtr <- build_X(tr, spec, zt, zm, zf_tr)
    Xte <- build_X(te, spec, zt, zm, zf_te)
    fit <- ridge_glm_v2(Xtr, tr$y, lam)
    te[, (paste0("p_", m)) := pred_ridge_v2(fit, Xte)]
    diag_list[[paste0(cg, "_", m)]] <- data.table(ConnGroup = cg, Model = m, lambda = lam,
      converged = fit$converged, iterations = fit$iterations, separation_suspected = fit$separation_suspected,
      n_train = nrow(Xtr), n_feat = ncol(Xtr) - 1, fallback_rate_test = mean(te$mp_source == "fallback_overall"))
    coef_list[[paste0(cg, "_", m)]] <- data.table(ConnGroup = cg, Model = m, lambda = lam,
      Feature = c("(Intercept)", spec), Coef = round(fit$beta, 5))
  }
  te[, ConnGroup_out := cg]
  keep <- c("GroupID","ReferenceKey","ConnGroup","ConnGroup_out","PXD","BaseAccession","EnsemblGeneID","y", paste0("p_", MODELS))
  preds_list[[as.character(cg)]] <- te[, ..keep]
  scal_list[[as.character(cg)]] <- data.table(ConnGroup = cg, mu_target = mu_t, sd_target = sd_t, imp_target = imp_t,
    mu_mod = mu_m, sd_mod = sd_m, imp_mod = imp_m, mu_mp = mu_f, sd_mp = sd_f)
}
P <- rbindlist(preds_list)
fwrite(P, file.path(out, "BRv2_oos_preds.csv.gz"))
diag <- rbindlist(diag_list); fwrite(diag, file.path(out, "BRv2_fit_diagnostics.csv"))
coefs <- rbindlist(coef_list); fwrite(coefs, file.path(out, "BRv2_coefficients.csv"))
fwrite(rbindlist(scal_list), file.path(out, "BRv2_scaler_params.csv"))
fwrite(rbindlist(mp_cov_list), file.path(out, "BRv2_MP_coverage.csv"))
## 汇总：逐折 / 19折等权 / 行数加权 / pooled（pooled仅参考；禁跨折常数比较的pooled ROC）
perf_list <- list()
for (m in MODELS) {
  per <- P[, metrics_row(get("y"), get(paste0("p_", m)), m, ConnGroup)]
  eqw <- per[, lapply(.SD, function(x) if (all(is.na(x))) NA_real_ else round(mean(x, na.rm=TRUE), 4)), .SDcols=c("AP","PRAUC_trap","ROC","Brier","LogLoss")]
  wgt <- per[, lapply(.SD, function(x) if (all(is.na(x))) NA_real_ else round(weighted.mean(x, N, na.rm=TRUE), 4)), .SDcols=c("AP","PRAUC_trap","ROC","Brier","LogLoss")]
  pl <- P[, metrics_row(get("y"), get(paste0("p_", m)), m, NA)]
  perf_list[[m]] <- rbind(
    cbind(data.table(Model=m, Summary="equal_weight_19folds"), eqw),
    cbind(data.table(Model=m, Summary="n_weighted"), wgt),
    cbind(data.table(Model=m, Summary="pooled_reference_only"), pl[, .(AP, PRAUC_trap, ROC, Brier, LogLoss)])
  )
  fwrite(per, file.path(out, paste0("BRv2_perfold_", m, ".csv")))
}
perf <- rbindlist(perf_list)
fwrite(perf, file.path(out, "BRv2_performance_summary.csv"))
print(perf[Summary == "equal_weight_19folds"])
cat("V2 TRAIN DONE\n")

