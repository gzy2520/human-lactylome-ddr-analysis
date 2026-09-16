#!/usr/bin/env Rscript
# Exploratory: what the transcriptome can say about the proteome's binary DDR readout.
#
# The proteome side records, per material group, only *which* proteins were captured - in the
# lactylome (Kla) and in the reference proteome. `data/publication_input/kla_protein_membership_31.csv`
# has exactly four columns (PXD, SampleGroup, BaseAccession, IsDdr) and no quantitative value at
# all, and the reference membership table is the same shape. So an absent protein is ambiguous:
# it may not be expressed, or it may be expressed and simply not captured as lactylated.
#
# The 31-group qsmooth reference profile supplies a continuous value per gene per material class
# and can separate those two cases. This script builds that comparison and reports the three
# quantities it supports:
#
#   1. how much of the proteome's absence is explained by absence of the transcript
#   2. per group, the lactylation excess over what the group's expression would predict
#   3. per pathway and per protein, where the expression/detection relationship breaks down
#
# Read as exploratory. The RNA reference is a material-class profile from different studies than
# the proteome, so this compares classes, not paired samples.
#
# Usage: explore_rna_assisted_ddr_20260917.R <project_root> [out_dir]
args <- commandArgs(TRUE)
stopifnot(length(args) >= 1L)
root <- normalizePath(args[[1L]], mustWork = TRUE)
out_dir <- if (length(args) >= 2L) args[[2L]] else file.path(root, "outputs", "20260917_rna_assisted_ddr")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages(library(data.table))
set.seed(25)

## ---- inputs --------------------------------------------------------------

expr_dir <- file.path(root, "outputs", "20260916_qsmooth_31group")
panel_dir <- file.path(root, "outputs", "20260916_ddr_panel_31group")

mat <- fread(file.path(expr_dir, "matrices", "qsmooth_A_collapsed_log2tpm.tsv.gz"))
expansion <- fread(file.path(expr_dir, "group_expansion_31.csv"))
ledger <- fread(file.path(root, "audit", "20260916_full_31_rna_status", "rna_31_group_status.csv"))

rna <- as.matrix(mat[, -1L]); rownames(rna) <- mat[[1L]]
rna31 <- rna[, match(expansion$ReferenceKey, colnames(rna)), drop = FALSE]
colnames(rna31) <- expansion$GroupID

kla <- fread(file.path(root, "data", "publication_input", "kla_protein_membership_31.csv"))
ref <- fread(file.path(root, "data", "publication_input", "reference_protein_membership_31.csv"))

mapping <- fread(file.path(panel_dir, "ddr_uniprot_to_ensembl.tsv"))
annotation <- fread(file.path(panel_dir, "kla_ddr_annotation.csv"))

# only the unambiguous accession -> gene pairs whose gene is actually in the profile
panel <- mapping[NGeneIDs == 1L & EnsemblGeneIDs != "" & EnsemblGeneIDs %in% rownames(rna31) &
                   grepl("KlaDdr", PanelMembership)]
cat(sprintf("Kla n DDR panel mapped 1:1 into the profile: %d proteins\n", nrow(panel)))

## ---- the pair table: group x panel protein -------------------------------

ledger[, GroupKey := paste(PXD, SampleGroup, sep = "__")]
kla[, GroupKey := paste(PXD, SampleGroup, sep = "__")]
ref[, GroupKey := paste(PXD, SampleGroup, sep = "__")]
ref[, BaseAccession := MappedBaseAccessions]

pairs <- rbindlist(lapply(seq_len(nrow(ledger)), function(i) {
  g <- ledger$GroupID[i]; key <- ledger$GroupKey[i]
  kset <- kla[GroupKey == key, BaseAccession]
  rset <- unique(ref[GroupKey == key, BaseAccession])
  data.table(
    GroupID = g, BaseAccession = panel$BaseAccession, EnsemblGeneID = panel$EnsemblGeneIDs,
    KlaDetected = panel$BaseAccession %in% kset,
    RefDetected = panel$BaseAccession %in% rset,
    RNA = rna31[panel$EnsemblGeneIDs, match(g, colnames(rna31))]
  )
}))
pairs[, ExprPct := frank(RNA, ties.method = "average") / .N, by = GroupID]

## ---- 1. absence of protein versus absence of transcript -------------------

not_kla <- pairs[KlaDetected == FALSE]
absence <- data.table(
  Metric = c("group_x_protein_pairs", "kla_detected", "reference_detected",
             "not_kla_pairs", "not_kla_below_group_median_expression",
             "not_kla_above_group_median_expression",
             "not_kla_above_median_and_reference_detected"),
  Value = c(nrow(pairs), sum(pairs$KlaDetected), sum(pairs$RefDetected),
            nrow(not_kla), sum(not_kla$ExprPct <= 0.5), sum(not_kla$ExprPct > 0.5),
            nrow(not_kla[ExprPct > 0.5 & RefDetected == TRUE]))
)
absence[, Percent := round(100 * Value / nrow(pairs), 2)]

# sensitivity: the headline share depends on where the expression cut is put, so report the sweep
sweep <- rbindlist(lapply(c(0, 0.5, 1, 2, 3, 4), function(th) {
  data.table(AbsoluteLog2TPMCut = th, ShareOfAbsentPairsAbove = round(100 * mean(not_kla$RNA > th), 1))
}))
quantile_sweep <- rbindlist(lapply(c(0.25, 0.5, 0.75), function(q) {
  data.table(WithinGroupPercentileCut = q, ShareOfAbsentPairsAbove = round(100 * mean(not_kla$ExprPct > q), 1))
}))

## ---- 2. per-group lactylation excess over expression ----------------------

# empirical detection probability as a function of expression decile, pooled over all pairs:
# deliberately a binned rate rather than a fitted model, so nothing is extrapolated
pairs[, RNAbin := cut(RNA, breaks = unique(quantile(RNA, probs = seq(0, 1, 0.1))),
                      include.lowest = TRUE, labels = FALSE)]
bin_rate <- pairs[, .(pDetect = mean(KlaDetected), nPairs = .N, RNAmid = median(RNA)), by = RNAbin]
pairs <- merge(pairs, bin_rate[, .(RNAbin, pDetect)], by = "RNAbin")

by_group <- pairs[, .(Observed = sum(KlaDetected), Expected = sum(pDetect), N = .N,
                      MeanRNA = mean(RNA)), by = GroupID]
by_group[, `:=`(Excess = Observed - Expected, ExcessRate = (Observed - Expected) / N)]
by_group <- merge(by_group, ledger[, .(GroupID, PXD, SampleGroup, Category)], by = "GroupID")
setorder(by_group, -ExcessRate)

cat("\nPer-group lactylation excess over expression expectation:\n")
print(by_group[, .(GroupID, SampleGroup = substr(SampleGroup, 1, 26), Category,
                   Obs = Observed, Exp = round(Expected, 1), ExcessRate = round(ExcessRate, 3))])
cat("\nBy category:\n")
print(by_group[, .(n = .N, MeanExcessRate = round(mean(ExcessRate), 3)), by = Category][order(-MeanExcessRate)])
cat(sprintf("\nSpearman(excess rate, mean RNA) = %.3f\n",
            cor(by_group$ExcessRate, by_group$MeanRNA, method = "spearman")))

## ---- 3. per pathway and per protein --------------------------------------

pathways <- c("BER", "NER", "MMR", "FA", "HR", "AEJ", "NHEJ")
long <- melt(annotation[, c("BaseAccession", pathways), with = FALSE], id.vars = "BaseAccession",
             variable.name = "Pathway", value.name = "State")
long <- long[State %in% c("1", "-1")]
by_pathway <- merge(long[, .(BaseAccession, Pathway)], pairs, by = "BaseAccession", allow.cartesian = TRUE)[
  , .(Pairs = .N, Proteins = uniqueN(BaseAccession),
      KlaDetectionRate = round(100 * mean(KlaDetected), 1),
      ReferenceDetectionRate = round(100 * mean(RefDetected), 1),
      HighExpressionRate = round(100 * mean(ExprPct > 0.5), 1)),
  by = Pathway][order(-KlaDetectionRate)]

by_protein <- pairs[, .(EnsemblGeneID = EnsemblGeneID[1L], MeanRNA = mean(RNA), MinRNA = min(RNA),
                        KlaRate = mean(KlaDetected), RefRate = mean(RefDetected), Groups = .N),
                    by = BaseAccession]
by_protein <- merge(by_protein, annotation[, c("BaseAccession", pathways), with = FALSE],
                    by = "BaseAccession", all.x = TRUE)
by_protein[, Pathway := apply(.SD, 1L, function(v) paste(pathways[!is.na(v) & v != ""], collapse = ";")),
           .SDcols = pathways]
setorder(by_protein, -MeanRNA)

## ---- write ---------------------------------------------------------------

fwrite(pairs, file.path(out_dir, "group_by_protein_pairs.tsv.gz"), sep = "\t")
fwrite(absence, file.path(out_dir, "absence_decomposition.csv"))
fwrite(sweep, file.path(out_dir, "expression_threshold_sweep.csv"))
fwrite(quantile_sweep, file.path(out_dir, "within_group_percentile_sweep.csv"))
fwrite(by_group, file.path(out_dir, "group_lactylation_excess.csv"))
fwrite(by_pathway, file.path(out_dir, "pathway_summary.csv"))
fwrite(by_protein, file.path(out_dir, "protein_summary.csv.gz"), sep = "\t")
writeLines(trimws(capture.output(sessionInfo()), which = "right"),
           file.path(out_dir, "sessionInfo.txt"))
writeLines(c(
  "# RNA-assisted view of the DDR proteome (exploratory, 2026-09-17)", "",
  "The proteome records only which proteins were captured, so an absent protein is ambiguous.",
  "The 31-group qsmooth profile supplies a continuous value per gene and separates",
  "'not expressed' from 'expressed but not captured as lactylated'.", "",
  "Exploratory: the RNA reference is a material-class profile drawn from different studies than",
  "the proteome, so groups are compared as classes, not as paired samples.",
  "`expression_threshold_sweep.csv` and `within_group_percentile_sweep.csv` exist because the",
  "headline share depends on where the expression cut is placed."
), file.path(out_dir, "README.md"))
cat(sprintf("\nEXPLORE_DONE -> %s\n", out_dir))
