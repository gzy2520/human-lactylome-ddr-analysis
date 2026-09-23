#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
out <- args[[1L]]
read <- function(name) fread(file.path(out, name))
base <- fread("outputs/20260923_rna_only_kla_fraction/material_labels_28.csv")
p <- read("material_oof_predictions.csv")
s <- read("selected_parameters.csv")
perf <- read("performance_summary.csv")
catalog <- read("module_catalog.csv")
hash <- read("input_md5.csv")
stopifnot(nrow(base) == 28L, nrow(p) == 28L, nrow(s) == 15L,
          nrow(catalog) == 16L, all(catalog$MeasuredGenes > 0L),
          !anyDuplicated(p$ReferenceKey), !anyDuplicated(s$OuterConnGroup),
          setequal(p$ReferenceKey, base$ReferenceKey),
          setequal(p$ConnGroup, s$OuterConnGroup),
          all(file.exists(hash$File)),
          identical(unname(tools::md5sum(hash$File)), hash$MD5),
          all(is.finite(as.matrix(p[, setdiff(names(p),
                                           c("ReferenceKey", "ConnGroup")),
                                    with = FALSE]))))
for (i in seq_len(nrow(p))) {
  tr <- base[ConnGroup != p$ConnGroup[[i]]]
  stopifnot(abs(p$MedianBaseline[[i]] - median(tr$ObservedPercent)) < 1e-10)
  w <- 1 / as.vector(table(tr$ConnGroup)[as.character(tr$ConnGroup)])
  ix <- order(tr$ObservedPercent)
  wm <- tr$ObservedPercent[ix][which(cumsum(w[ix]) >= sum(w) / 2)[1L]]
  stopifnot(abs(p$BalancedMedianBaseline[[i]] - wm) < 1e-10)
  chosen <- s[OuterConnGroup == p$ConnGroup[[i]], Model]
  stopifnot(abs(p$SelectedWithFallback[[i]] - p[[chosen]][[i]]) < 1e-10)
}
for (i in seq_len(nrow(perf))) {
  nm <- perf$Model[[i]]
  err <- p[[nm]] - p$ObservedPercent
  macro <- p[, .(MAE = mean(abs(get(nm) - ObservedPercent))),
             by = ConnGroup]
  stopifnot(abs(perf$MaterialMAE_pp[[i]] - mean(abs(err))) < 1e-10,
            abs(perf$MaterialRMSE_pp[[i]] - sqrt(mean(err^2))) < 1e-10,
            abs(perf$ConnectedGroupMeanMAE_pp[[i]] - mean(macro$MAE)) < 1e-10)
}
cat("PASS: frozen module coverage, 15 source-held-out folds,",
    "held-out medians, fold selections, metrics, and input fingerprints.\n")
