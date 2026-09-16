#!/usr/bin/env Rscript
# Verification of the 31-group expression matrices. Fails loudly rather than reporting a
# matrix as usable when it is not: every claim below is checked against the stored objects.
#
# Usage: verify_31_expression_matrices_20260916.R <out_dir>
args <- commandArgs(TRUE)
stopifnot(length(args) == 1L)
out_dir <- normalizePath(args[[1L]], mustWork = TRUE)

obj_files <- sort(list.files(file.path(out_dir, "objects"), pattern = "\\.rds$", full.names = TRUE))
stopifnot(length(obj_files) > 0L)

fail <- character(0)
note <- function(ok, msg) if (!isTRUE(ok)) fail <<- c(fail, msg)

rows <- list()
for (f in obj_files) {
  o <- readRDS(f)
  k <- o$ReferenceKey
  tpm <- o$tpm; l2 <- o$log2_tpm; cnt <- o$counts

  note(is.matrix(tpm), paste(k, "tpm is not a matrix"))
  note(!anyDuplicated(rownames(tpm)), paste(k, "duplicate gene identifiers"))
  note(!anyDuplicated(colnames(tpm)), paste(k, "duplicate sample identifiers"))
  note(all(grepl("^ENSG[0-9]+(_PAR_Y)?$", rownames(tpm))), paste(k, "row keys are not bare Ensembl gene IDs"))
  note(all(is.finite(tpm)) && all(tpm >= 0), paste(k, "tpm has non-finite or negative values"))
  note(identical(dim(tpm), dim(l2)), paste(k, "log2 TPM dimensions differ from TPM"))
  note(max(abs(log2(tpm + 0.5) - l2)) < 1e-8, paste(k, "log2 TPM is not log2(TPM + 0.5)"))
  totals <- colSums(tpm)
  note(all(abs(totals - 1e6) < 1), paste(k, "TPM columns do not each sum to 1e6"))

  dropped_no_length <- NA_integer_
  if (!is.null(cnt)) {
    # the TPM matrix holds only genes that have a length in the shared annotation, so it is
    # a subset of the archived counts rather than dimensionally identical to it
    note(identical(colnames(cnt), colnames(tpm)), paste(k, "counts and TPM sample order differ"))
    note(all(rownames(tpm) %in% rownames(cnt)), paste(k, "TPM has genes absent from counts"))
    note(all(is.finite(cnt)) && all(cnt >= 0), paste(k, "counts have non-finite or negative values"))
    note(all(cnt == floor(cnt)), paste(k, "counts are not integers"))
    dropped_no_length <- nrow(cnt) - nrow(tpm)
    # only meaningful where the source values really are integer counts: for a source whose
    # archived "counts" are rounded fractional values, small values legitimately round to 0
    if (isTRUE(o$Provenance$NativeValuesInteger)) {
      sub <- cnt[rownames(tpm), , drop = FALSE]
      nz_cnt <- sum(rowSums(sub) > 0); nz_tpm <- sum(rowSums(tpm) > 0)
      note(abs(nz_cnt - nz_tpm) <= max(5, 0.01 * nz_cnt),
           sprintf("%s: %d genes carry counts but %d carry TPM", k, nz_cnt, nz_tpm))
    }
  }

  q <- o$qc
  note(all(is.finite(q$TPMTotal)), paste(k, "QC table has non-finite values"))
  rows[[k]] <- data.frame(ReferenceKey = k, GroupIDs = o$GroupIDs, Samples = ncol(tpm),
                          Genes = nrow(tpm),
                          GenesWithoutLength = dropped_no_length,
                          MedianLibrarySize = if (is.null(cnt)) NA_real_ else median(colSums(cnt)),
                          # Inf here means the group has a single sample, so no pairwise
                          # correlation exists to report
                          MinPairwiseSpearman = suppressWarnings(min(q$PairwiseSpearmanMin, na.rm = TRUE)),
                          MedianZeroFraction = median(q$ZeroFraction),
                          ValueRoute = o$Provenance$NativeValueRoute,
                          stringsAsFactors = FALSE)
}

if (length(fail)) {
  cat("VERIFICATION_FAILED\n"); cat(paste0(" - ", fail, collapse = "\n"), "\n")
  quit(status = 1L)
}
summary_tab <- do.call(rbind, rows)
rownames(summary_tab) <- NULL
write.csv(summary_tab, file.path(out_dir, "verification_summary.csv"), row.names = FALSE)
cat(sprintf("VERIFICATION_PASS matrices=%d samples=%d genes_min=%d genes_max=%d\n",
            nrow(summary_tab), sum(summary_tab$Samples),
            min(summary_tab$Genes), max(summary_tab$Genes)))
