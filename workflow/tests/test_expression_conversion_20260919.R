#!/usr/bin/env Rscript
source("workflow/lib_kla31_expression_20260916.R")
set.seed(25)
ids <- sprintf("ENSG%011d", seq_len(1200))
x <- matrix(rep(seq_len(1200), 2), ncol = 2,
            dimnames = list(ids, c("sample_a", "sample_b")))
lengths <- setNames(rep(c(1000, 10000), 600), ids)
expected <- sweep(x, 2, colSums(x), "/") * 1e6
stopifnot(max(abs(fpkm_to_tpm(x, lengths) - expected)) < 1e-9,
          max(abs(fpkm_to_tpm(x, rev(lengths)) - expected)) < 1e-9)
bad <- sweep(x, 1, lengths / 1000, "/")
stopifnot(max(abs(rescale_to_tpm(bad) - expected)) > 100)
expect_error <- function(expr) stopifnot(inherits(tryCatch({force(expr); NULL}, error = identity), "error"))
expect_error(rescale_to_tpm(x * 0))
x_bad <- x; x_bad[1,1] <- NA_real_; expect_error(rescale_to_tpm(x_bad))
d <- tempfile(); dir.create(d); assert_empty_output(d)
writeLines("existing cache", file.path(d, "old.rds")); expect_error(assert_empty_output(d))
unlink(d, recursive = TRUE)
cat("PASS: FPKM formula, length-order invariance, invalid inputs, stale-output rejection\n")
# Existing Stage2 output must fail before annotation loading and leave its contents intact.
d <- tempfile(); dir.create(d); writeLines('preserved', file.path(d,'sentinel'))
log <- tempfile()
rc <- system2(file.path(R.home('bin'),'Rscript'),c('--vanilla',
  'workflow/build_31_expression_matrices_20260916.R','.',
  'config/rnaseq_expression_extraction_contract_20260916.csv',shQuote(d)),stdout=log,stderr=log)
stopifnot(rc!=0L,identical(readLines(file.path(d,'sentinel')),'preserved'),
          any(grepl('Refusing non-empty output',readLines(log))))
unlink(c(d,log),recursive=TRUE)
# Incomplete object sets must fail before final manifests are written.
d <- tempfile(); dir.create(file.path(d,'objects'),recursive=TRUE)
file.create(file.path(d,'objects','partial.rds')); log<-tempfile()
rc <- system2(file.path(R.home('bin'),'Rscript'),c('--vanilla',
  'workflow/finalize_31_expression_matrices_20260916.R',shQuote(d),
  'config/rnaseq_expression_extraction_contract_20260916.csv'),stdout=log,stderr=log)
stopifnot(rc!=0L,!file.exists(file.path(d,'group_manifest.csv')))
unlink(c(d,log),recursive=TRUE)
cat('PASS: existing output protected; incomplete finalize rejected before metadata writes\n')
