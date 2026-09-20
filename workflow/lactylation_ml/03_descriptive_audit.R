#!/usr/bin/env Rscript
# 03_descriptive_audit.R — leakage-free descriptive audit ONLY (no training, no p-hacking).
# Uses only local frozen aggregates; every number is a description of detection process,
# NOT a predictive performance claim. Seed 25 where randomness is involved (none here).
suppressPackageStartupMessages({library(data.table)})
set.seed(25)
root <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE))[1L]), "..", ".."))
setwd(root)
out_dir <- "outputs/20260920_lactylation_ml"
obs <- fread("outputs/20260920_lactylation_ml/observation_contract.csv")

# D1: Kla counts vs detection depth (confounding diagnostic)
obs[, log10Ref := log10(ReferenceProteinCount)]
ct <- cor.test(obs$KlaProteinCount, obs$log10Ref, method="spearman", exact=FALSE)
d1 <- data.table(Metric=c("n_materials", "spearman_Kla_vs_log10RefDepth", "p_value",
  "min_RefDepth", "max_RefDepth", "min_Kla", "max_Kla"),
  Value=c(nrow(obs), round(unname(ct$estimate), 4), signif(ct$p.value, 4),
    min(obs$ReferenceProteinCount), max(obs$ReferenceProteinCount),
    min(obs$KlaProteinCount), max(obs$KlaProteinCount)))
fwrite(d1, file.path(out_dir, "audit_detection_confounding.csv"))

# D2: Kla counts by material category (study==material confounding note)
d2 <- obs[, .(n=.N, median_Kla=median(KlaProteinCount), min_Kla=min(KlaProteinCount),
  max_Kla=max(KlaProteinCount), median_RefDepth=median(ReferenceProteinCount)), by=Category]
d2[, Note := "Each Category pools distinct PXDs; study and material effects inseparable"]
fwrite(d2, file.path(out_dir, "audit_by_category.csv"))

# D3: matched-modality coverage (which rows can even define a denominator-based label)
m <- fread("data/candidate/sample_resolved_matched_modalities.csv")
d3 <- data.table(Metric=c("rows_with_same_sample_Kla_and_Ref", "distinct_PXD_with_match", "distinct_groups_covered"),
  Value=c(nrow(m), length(unique(m$PXD)), length(unique(m$PublicationGroup))))
fwrite(d3, file.path(out_dir, "audit_matched_coverage.csv"))

# D4: regulator RNA spread across 28 refs (feature-variation diagnostic; descriptive only)
reg <- fread("outputs/20260919_repair_release/core/regulator_rna_qsmooth_28materials.csv")
d4 <- reg[, .(n_refs=uniqueN(ReferenceKey), min=min(QsmoothLog2TPM), max=max(QsmoothLog2TPM),
  median=median(QsmoothLog2TPM), sd=sd(QsmoothLog2TPM)), by=.(DisplayName, Role, EnsemblGeneID)]
setorder(d4, DisplayName)
fwrite(d4, file.path(out_dir, "audit_regulator_spread.csv"))

# D5: connectivity (independent-unit accounting)
conn <- fread("outputs/20260920_lactylation_ml/connectivity_groups.csv")
d5 <- conn[, .(ConnGroup, n_rows, members, refs, pxds)]
fwrite(d5, file.path(out_dir, "audit_connectivity.csv"))
cat(sprintf("AUDIT_OK: D1 spearman=%.3f p=%s; matched_rows=%d; conn_groups=%d\n",
  unname(ct$estimate), signif(ct$p.value, 3), nrow(m), nrow(d5)))

# D6: figure — Kla counts vs detection depth (per-figure source data saved alongside)
if (requireNamespace("ggplot2", quietly=TRUE)) {
  suppressPackageStartupMessages(library(ggplot2))
  p <- ggplot(obs, aes(log10Ref, KlaProteinCount, color=Category)) +
    geom_point(size=2.6) +
    geom_text(aes(label=GroupID), size=2.4, vjust=-0.7, show.legend=FALSE) +
    labs(x="log10 whole-proteome detected proteins (detection depth)",
      y="Kla detected protein count (detection feature; NOT concentration)",
      title="Kla detection counts covary with detection depth",
      subtitle="31 materials; cross-study counts NOT comparable as unified label; seed=25 (no randomness in plot)") +
    theme_minimal(base_size=10)
  ggsave(file.path(out_dir, "fig_audit_detection_confounding.png"), p, width=7.5, height=5, dpi=200)
  ggsave(file.path(out_dir, "fig_audit_detection_confounding.pdf"), p, width=7.5, height=5)
  fwrite(obs[, .(GroupID, PXD, SampleGroup, Category, KlaProteinCount, ReferenceProteinCount, log10Ref)],
    file.path(out_dir, "fig_audit_detection_confounding_sourcedata.csv"))
  cat("FIG_OK: fig_audit_detection_confounding.{png,pdf} + sourcedata\n")
}
