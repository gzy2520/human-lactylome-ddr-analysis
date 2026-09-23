#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
out <- args[[1L]]
if (dir.exists(out) && length(list.files(out, all.files = TRUE, no.. = TRUE)))
  stop("Refusing to overwrite nonempty output")
paths <- c(
  original = "outputs/20260923_rna_only_kla_fraction_optimization",
  GSE178411_replace = "outputs/20260924_rna_ref_expansion_replace",
  equal_source = "outputs/20260924_rna_ref_expansion_equal_source"
)
files <- as.vector(unlist(lapply(paths, function(path)
  file.path(path, c("performance_summary.csv", "material_oof_predictions.csv")))))
stopifnot(all(file.exists(files)))
perf <- rbindlist(lapply(names(paths), function(nm) {
  d <- fread(file.path(paths[[nm]], "performance_summary.csv"))
  d[, Scenario := nm]
  d
}))
preds <- lapply(paths, function(path)
  fread(file.path(path, "material_oof_predictions.csv")))
original <- preds[["original"]]
for (nm in names(preds)) {
  p <- preds[[nm]]
  setorder(p, ReferenceKey)
  stopifnot(nrow(p) == 28L,
            identical(p$ReferenceKey, original$ReferenceKey),
            identical(p$ConnGroup, original$ConnGroup),
            identical(p$ObservedPercent, original$ObservedPercent),
            identical(p$BalancedMedianBaseline,
                      original$BalancedMedianBaseline))
}
paired <- rbindlist(lapply(names(preds), function(nm) {
  p <- preds[[nm]]
  p[, .(Scenario = nm, ReferenceKey, ConnGroup, ObservedPercent,
        PredictedPercent = All16_ridge,
        AbsoluteError_pp = abs(All16_ridge - ObservedPercent),
        BalancedBaselineAbsoluteError_pp =
          abs(BalancedMedianBaseline - ObservedPercent))]
}))
scar <- paired[ReferenceKey == "GSE181540_HS_input"]
setorder(perf, Scenario, Model)
setorder(paired, Scenario, ReferenceKey)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
fwrite(perf, file.path(out, "performance_comparison.csv"))
fwrite(paired, file.path(out, "paired_material_errors_all16.csv"))
fwrite(scar, file.path(out, "scar_heldout_prediction.csv"))
fwrite(data.table(File = files,
                  MD5 = unname(tools::md5sum(files))),
       file.path(out, "input_md5.csv"))
print(perf[Model %in% c("All16_ridge", "BalancedMedianBaseline",
                        "SelectedWithFallback"),
           .(Scenario, Model, MaterialMAE_pp,
             ConnectedGroupMeanMAE_pp, MaterialRMSE_pp)])
