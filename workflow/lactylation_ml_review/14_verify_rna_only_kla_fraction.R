#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
out <- args[[1L]]
read <- function(name) fread(file.path(out, name))
t <- read("protein_group_labels_31.csv")
m <- read("material_labels_28.csv")
p <- read("material_oof_predictions_28.csv")
gp <- read("protein_group_oof_predictions_31.csv")
components <- read("source_components_31.csv")
selected <- read("selected_parameters.csv")
metrics <- read("performance_summary.csv")
hashes <- read("input_md5.csv")

stopifnot(nrow(t) == 31L, nrow(m) == 28L, nrow(p) == 28L,
          nrow(gp) == 31L, nrow(components) == 31L,
          uniqueN(components$ConnGroup) == 15L,
          !anyDuplicated(p$ReferenceKey),
          setequal(p$ReferenceKey, m$ReferenceKey),
          setequal(t$GroupID, gp$GroupID),
          all(file.exists(hashes$File)),
          identical(unname(tools::md5sum(hashes$File)), hashes$MD5),
          all(abs(t$DetectedKlaPercent -
                    100 * t$MatchedUniqueAccessions / t$ReferenceUniqueAccessions) < 1e-10),
          all(t$KlaOutsideReference == t$KlaUniqueAccessions -
                t$MatchedUniqueAccessions))
for (key in c("PXD", "ReferencePXD", "ReferenceKey"))
  stopifnot(components[, all(uniqueN(ConnGroup) == 1L), by = key][, all(V1)])
for (i in seq_len(nrow(m))) {
  z <- t[ReferenceKey == m$ReferenceKey[[i]]]
  stopifnot(nrow(z) == m$ProteomeGroups[[i]],
            abs(mean(z$DetectedKlaPercent) - m$ObservedPercent[[i]]) < 1e-10)
}
for (i in seq_len(nrow(p))) {
  train <- m[ConnGroup != p$ConnGroup[[i]]]
  stopifnot(abs(p$Baseline[[i]] - mean(train$ObservedPercent)) < 1e-10,
            abs(p$MedianBaseline[[i]] - median(train$ObservedPercent)) < 1e-10)
  choice <- selected[OuterConnGroup == p$ConnGroup[[i]], Model]
  stopifnot(length(choice) == 1L,
            abs(p$SelectedRNA[[i]] - p[[choice]][[i]]) < 1e-10)
}
for (j in seq_len(nrow(metrics))) {
  name <- metrics$Model[[j]]
  err <- p[[name]] - p$ObservedPercent
  ge <- gp[[name]] - gp$DetectedKlaPercent
  comp_mae <- p[, .(MAE = mean(abs(get(name) - ObservedPercent))),
                by = ConnGroup]
  stopifnot(abs(metrics$MaterialMAE_pp[[j]] - mean(abs(err))) < 1e-10,
            abs(metrics$MaterialRMSE_pp[[j]] - sqrt(mean(err^2))) < 1e-10,
            abs(metrics$ConnectedGroupMeanMAE_pp[[j]] - mean(comp_mae$MAE)) < 1e-10,
            abs(metrics$Original31GroupMAE_pp[[j]] - mean(abs(ge))) < 1e-10)
}
cat("PASS: 31 labels, 28 independent RNA materials, 15 source-connected folds,",
    "held-out baselines, selected predictions, metrics, and input fingerprints.\n")
