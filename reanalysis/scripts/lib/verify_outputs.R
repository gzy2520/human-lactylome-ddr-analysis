#!/usr/bin/env Rscript
# 用法: Rscript lib/verify_outputs.R <project_root> <baseline.csv> [排除正则...]
# 校验基线列出的文件是否仍与基线字节一致；用于全 R 重构期间的回归检查。
# 注意：本文件是重构期间工具，不属于合并目标，重构完成后保留在 lib/。
args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(args[[1]])
baseline <- read.csv(args[[2]], check.names = FALSE, stringsAsFactors = FALSE)
exclude <- if (length(args) >= 3) args[[3]] else NULL
sha1 <- function(p) {
  digest::digest(file = p, algo = "sha256", serialize = FALSE)
}
changed <- list()
for (i in seq_len(nrow(baseline))) {
  rel <- baseline$RelativePath[i]
  p <- file.path(root, rel)
  if (!is.null(exclude) && grepl(exclude, rel, perl = TRUE)) next
  if (!file.exists(p)) {
    changed[[length(changed) + 1]] <- data.frame(
      RelativePath = rel, Status = "MISSING", stringsAsFactors = FALSE
    )
    next
  }
  if (baseline$SizeBytes[i] != file.info(p)$size ||
      baseline$SHA256[i] != sha1(p)) {
    changed[[length(changed) + 1]] <- data.frame(
      RelativePath = rel, Status = "CHANGED", stringsAsFactors = FALSE
    )
  }
}
if (length(changed)) {
  out <- do.call(rbind, changed)
  write.csv(out, file.path(dirname(args[[2]]), "diff_report.csv"), row.names = FALSE)
  cat("DIFF FOUND:", nrow(out), "files\n")
  print(out)
  quit(status = 1)
}
cat("OK: all", nrow(baseline), "baseline files unchanged\n")
