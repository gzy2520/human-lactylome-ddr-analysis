#!/usr/bin/env Rscript
# Source-held-out comparison: RNA reference updates never add protein labels.
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
out <- args[[1L]]
if (dir.exists(out) && length(list.files(out, all.files = TRUE, no.. = TRUE)))
  stop("Refusing to overwrite nonempty output")
paths <- c(
  original = "outputs/20260923_rna_only_kla_fraction_optimization",
  scar_replace = "outputs/20260924_rna_ref_expansion_replace",
  scar_equal = "outputs/20260924_rna_ref_expansion_equal_source",
  tendon_replace = "outputs/20260924_rna_ref_tendon_replace",
  tendon_equal = "outputs/20260924_rna_ref_tendon_equal",
  tendon_degenerative = "outputs/20260924_rna_ref_tendon_degenerate",
  scar_tendon_equal = "outputs/20260924_rna_ref_combined_equal"
)
files <- unlist(lapply(paths, function(path)
  file.path(path, c("performance_summary.csv", "material_oof_predictions.csv",
                    "input_md5.csv"))), use.names = FALSE)
stopifnot(all(file.exists(files)))
preds <- lapply(paths, function(path)
  fread(file.path(path, "material_oof_predictions.csv")))
base <- copy(preds[["original"]])
setorder(base, ReferenceKey)
for (name in names(preds)) {
  setorder(preds[[name]], ReferenceKey)
  p <- preds[[name]]
  stopifnot(nrow(p) == 28L,
            identical(p$ReferenceKey, base$ReferenceKey),
            identical(p$ConnGroup, base$ConnGroup),
            identical(p$ObservedPercent, base$ObservedPercent),
            identical(p$BalancedMedianBaseline,
                      base$BalancedMedianBaseline))
}
models <- c("All16_ridge", "Mechanism10_ridge",
            "SelectedWithFallback", "BalancedMedianBaseline")
performance <- rbindlist(lapply(names(paths), function(name) {
  p <- fread(file.path(paths[[name]], "performance_summary.csv"))
  p[Model %in% models, Scenario := name][!is.na(Scenario),
     .(Scenario, Model, MaterialMAE_pp, ConnectedGroupMeanMAE_pp,
       MaterialRMSE_pp)]
}))
heldout <- rbindlist(lapply(names(preds), function(name) {
  p <- preds[[name]]
  p[ReferenceKey %in% c("GSE181540_HS_input", "tendon_nondiabetic"),
    .(Scenario = name, ReferenceKey, ObservedPercent,
      All16PredictedPercent = All16_ridge,
      All16AbsoluteError_pp = abs(All16_ridge - ObservedPercent),
      BaselinePredictedPercent = BalancedMedianBaseline,
      BaselineAbsoluteError_pp =
        abs(BalancedMedianBaseline - ObservedPercent))]
}))
component <- rbindlist(lapply(names(preds), function(name) {
  preds[[name]][,
    .(All16MAE_pp = mean(abs(All16_ridge - ObservedPercent)),
      BaselineMAE_pp = mean(abs(BalancedMedianBaseline - ObservedPercent)),
      N = .N), by = ConnGroup][, Scenario := name]
}))
setorder(performance, Scenario, Model)
setorder(heldout, Scenario, ReferenceKey)
setorder(component, Scenario, ConnGroup)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
fwrite(performance, file.path(out, "performance_comparison.csv"))
fwrite(heldout, file.path(out, "updated_material_heldout.csv"))
fwrite(component, file.path(out, "component_mae.csv"))
fwrite(data.table(File = files,
                  MD5 = unname(tools::md5sum(files))),
       file.path(out, "input_md5.csv"))
print(performance[Model %in% c("All16_ridge", "BalancedMedianBaseline")])
