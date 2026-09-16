#!/usr/bin/env Rscript
# Cross-material reference profile for the 31 group rows (2026-09-16).
#
# Takes the per-group log2(TPM + 0.5) matrices produced by the server-side extraction and
# reduces them to one comparable reference profile per material class using qsmooth
# (smooth quantile normalization, the method YARN wraps as normalizeTissueAware).
#
# Two inputs are smoothed so they can be cross-checked:
#   A  collapsed : one median profile per reference matrix          (17,340 x 28)
#   B  full      : every sample, with the reference matrix as group factor (17,340 x 1898)
# They answer different questions. A weights each material class equally, which is what a
# cross-tissue comparison wants; B keeps within-class biological variation but lets the
# classes with hundreds of samples (GTEx lung, TCGA PRAD) dominate the quantile estimates.
#
# This is a normalisation, not a test: no group is tested against another, and a low weight
# in one quantile region means only that the groups genuinely differ there.
#
# Usage: qsmooth_31group_20260916.R <expression_dir> <out_dir> <contract_csv>
args <- commandArgs(TRUE)
stopifnot(length(args) == 3L)
expr_dir <- normalizePath(args[[1L]], mustWork = TRUE)
out_dir <- args[[2L]]
contract_path <- normalizePath(args[[3L]], mustWork = TRUE)

suppressPackageStartupMessages({
  library(data.table)
  library(qsmooth)
})
set.seed(25)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
if (length(list.files(out_dir, all.files = TRUE, no.. = TRUE))) {
  stop("output directory already exists and is non-empty: ", out_dir, call. = FALSE)
}
dir.create(file.path(out_dir, "matrices"), showWarnings = FALSE)

## ---- assemble the common gene space --------------------------------------

index <- fread(file.path(expr_dir, "group_index.csv"))
refs <- unique(index$ReferenceKey)
stopifnot(length(refs) == 28L)

read_one <- function(key) {
  f <- file.path(expr_dir, "matrices", paste0(key, "_log2tpm.tsv.gz"))
  stopifnot(file.exists(f))
  d <- fread(f, sep = "\t", header = TRUE, check.names = FALSE)
  ids <- d[[1L]]
  m <- as.matrix(d[, -1L, with = FALSE])
  rownames(m) <- ids
  stopifnot(!anyDuplicated(ids), all(is.finite(m)))
  m
}

mats <- lapply(refs, read_one)
names(mats) <- refs
common <- Reduce(intersect, lapply(mats, rownames))
cat(sprintf("reference matrices: %d, common genes: %d\n", length(mats), length(common)))
stopifnot(length(common) > 10000L)

full <- do.call(cbind, lapply(mats, function(m) m[common, , drop = FALSE]))
group_factor <- factor(rep(refs, vapply(mats, ncol, integer(1L))))
stopifnot(ncol(full) == length(group_factor), all(table(group_factor) >= 1L))
cat(sprintf("full matrix: %d genes x %d samples across %d groups\n",
            nrow(full), ncol(full), nlevels(group_factor)))

collapsed <- sapply(refs, function(k) {
  sub <- mats[[k]][common, , drop = FALSE]
  if (ncol(sub) == 1L) sub[, 1L] else apply(sub, 1L, median)
})
stopifnot(nrow(collapsed) == length(common), ncol(collapsed) == length(refs))
cat(sprintf("collapsed matrix: %d genes x %d groups\n", nrow(collapsed), ncol(collapsed)))

write_matrix <- function(mat, path) {
  con <- gzfile(path, "wt"); on.exit(close(con))
  write.table(data.frame(stable_gene_id = rownames(mat), mat, check.names = FALSE), con,
              sep = "\t", quote = FALSE, row.names = FALSE)
}
write_matrix(full, file.path(out_dir, "matrices", "unified_full_log2tpm.tsv.gz"))
write_matrix(collapsed, file.path(out_dir, "matrices", "unified_groupmedian_log2tpm.tsv.gz"))

## ---- smooth ---------------------------------------------------------------

# qsmooth returns a "qsmooth" object holding the smoothed matrix and one weight per gene:
# the weight is 1 where a global quantile normalization is applied and 0 where the groups are
# left to their own quantiles, and it is smoothed across neighbouring genes.
run_qsmooth <- function(mat, groups, label) {
  cat(sprintf("qsmooth on %s (%d x %d) ...\n", label, nrow(mat), ncol(mat)))
  res <- qsmooth(object = mat, group_factor = groups)
  sm <- qsmoothData(res)
  w <- qsmoothWeights(res)
  stopifnot(identical(dim(sm), dim(mat)), !anyNA(sm), length(w) == nrow(mat))
  list(smoothed = sm, weights = w)
}

a <- run_qsmooth(collapsed, factor(refs), "A collapsed")
b <- run_qsmooth(full, group_factor, "B full")

write_matrix(a$smoothed, file.path(out_dir, "matrices", "qsmooth_A_collapsed_log2tpm.tsv.gz"))
write_matrix(b$smoothed, file.path(out_dir, "matrices", "qsmooth_B_full_log2tpm.tsv.gz"))
fwrite(data.table(stable_gene_id = rownames(a$smoothed), Weight = a$weights),
       file.path(out_dir, "qsmooth_weights_A.csv.gz"))
fwrite(data.table(stable_gene_id = rownames(b$smoothed), Weight = b$weights),
       file.path(out_dir, "qsmooth_weights_B.csv.gz"))

## ---- A versus B -----------------------------------------------------------

# collapse B the same way A was built, so the comparison is like for like
b_collapsed <- sapply(refs, function(k) {
  sub <- b$smoothed[, group_factor == k, drop = FALSE]
  if (ncol(sub) == 1L) sub[, 1L] else apply(sub, 1L, median)
})
delta <- a$smoothed - b_collapsed
cmp <- data.table(
  stable_gene_id = rownames(delta),
  MeanAbsDelta = rowMeans(abs(delta)),
  MaxAbsDelta = apply(abs(delta), 1L, max),
  Correlation = vapply(seq_len(nrow(delta)),
                       function(i) suppressWarnings(cor(a$smoothed[i, ], b_collapsed[i, ])), numeric(1L))
)
fwrite(cmp, file.path(out_dir, "AB_comparison.csv"))
fwrite(data.table(stable_gene_id = rownames(delta), delta),
       file.path(out_dir, "AB_delta_by_gene_and_group.csv.gz"))

## ---- expand to the 31 group rows -----------------------------------------

contract <- fread(contract_path, sep = ",", header = TRUE, check.names = FALSE)
ref2groups <- tapply(contract$GroupIDs, contract$ReferenceKey,
                     function(x) paste(unique(unlist(strsplit(x, ";", fixed = TRUE))), collapse = ";"))
by_ref <- data.table(ReferenceKey = refs,
                     Groups = unname(ref2groups[refs]),
                     Samples = vapply(mats, ncol, integer(1L)))
by_ref[, SharedReference := grepl(";", Groups)]
expanded <- rbindlist(lapply(seq_len(nrow(by_ref)), function(i) {
  data.table(GroupID = strsplit(by_ref$Groups[i], ";", fixed = TRUE)[[1L]],
             ReferenceKey = by_ref$ReferenceKey[i],
             Samples = by_ref$Samples[i],
             SharedReference = by_ref$SharedReference[i])
}))
setorder(expanded, GroupID)
stopifnot(nrow(expanded) == 31L, !anyDuplicated(expanded$GroupID))
fwrite(expanded, file.path(out_dir, "group_expansion_31.csv"))

## ---- summary --------------------------------------------------------------

gene_sd <- apply(a$smoothed, 1L, sd)
summary_lines <- c(
  sprintf("genes,%d", length(common)),
  sprintf("reference_matrices,%d", length(refs)),
  sprintf("group_rows,%d", nrow(expanded)),
  sprintf("samples_total,%d", ncol(full)),
  sprintf("weight_A_median,%.4f", median(a$weights)),
  sprintf("weight_A_min,%.4f", min(a$weights)),
  sprintf("weight_A_max,%.4f", max(a$weights)),
  sprintf("AB_mean_abs_delta,%.5f", mean(cmp$MeanAbsDelta)),
  sprintf("AB_median_abs_delta,%.5f", median(cmp$MeanAbsDelta)),
  sprintf("AB_p99_abs_delta,%.5f", quantile(cmp$MeanAbsDelta, 0.99)),
  sprintf("AB_genes_delta_gt_0.5,%d", sum(cmp$MeanAbsDelta > 0.5)),
  sprintf("AB_median_gene_correlation,%.4f", median(cmp$Correlation))
)
writeLines(summary_lines, file.path(out_dir, "summary.csv"))
cat(paste(summary_lines, collapse = "\n"), "\n")
writeLines(trimws(capture.output(sessionInfo()), which = "right"),
           file.path(out_dir, "sessionInfo.txt"))
cat("QSMOOTH_DONE\n")
