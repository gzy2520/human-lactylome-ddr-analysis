#!/usr/bin/env Rscript
# Task 11: R port of build_final_manifest.py
# Produces reanalysis/reports/final_file_manifest_sha256.csv with columns
# RelativePath,SizeBytes,SHA256 — byte-identical to the Python implementation.
#
# Usage: Rscript build_final_manifest.R <project_root>
suppressMessages(library(digest))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript build_final_manifest.R <project_root>", call. = FALSE)
}
ROOT <- normalizePath(args[1], winslash = "/", mustWork = TRUE)
OUTPUT <- file.path(ROOT, "reanalysis", "reports", "final_file_manifest_sha256.csv")

sha256_file <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

# Emulate Python pathlib suffix == ".pyc": a leading-dot name with no other
# dot (e.g. ".pyc") has an empty suffix in pathlib, so it is NOT excluded.
has_pyc_suffix <- function(name) {
  dots <- gregexpr("\\.", name)[[1]]
  if (dots[1] < 0) return(FALSE)
  last <- dots[length(dots)]
  if (last <= 1) return(FALSE)
  substr(name, last + 1, nchar(name)) == "pyc"
}

# Recursive file listing equivalent to Python's Path.rglob("*") filtered to
# files: all.files = TRUE includes dotfiles (rglob does too).
collect_tree <- function(dir) {
  if (!dir.exists(dir)) return(character(0))
  list.files(dir, recursive = TRUE, all.files = TRUE, full.names = TRUE)
}

paths <- character(0)
# Root-level documents (only if present)
for (f in c("PROJECT_INDEX.md", "NEW_CHAT_PROJECT_PROMPT.md")) {
  p <- file.path(ROOT, f)
  if (file.exists(p)) paths <- c(paths, p)
}
# Full trees
paths <- c(paths, collect_tree(file.path(ROOT, "reanalysis")))
paths <- c(paths, collect_tree(file.path(ROOT, "previous_umap")))
# data/PXD*/metadata trees, PXD dirs visited in sorted order
data_dir <- file.path(ROOT, "data")
if (dir.exists(data_dir)) {
  pxd_dirs <- list.files(data_dir, pattern = "^PXD", full.names = TRUE)
  pxd_dirs <- pxd_dirs[dir.exists(pxd_dirs)]
  for (pxd in sort(pxd_dirs, method = "radix")) {
    paths <- c(paths, collect_tree(file.path(pxd, "metadata")))
  }
}

paths <- unique(paths)
# Keep regular files only (Python: path.is_file())
keep <- file.exists(paths) & !dir.exists(paths)
paths <- paths[keep]
# Exclude .DS_Store, *.pyc and the manifest itself
paths <- paths[basename(paths) != ".DS_Store"]
paths <- paths[!vapply(basename(paths), has_pyc_suffix, logical(1))]
paths <- paths[paths != OUTPUT]

# Sort key: Python sorts by absolute posix path string; since every path
# shares ROOT as a prefix, sorting by the relative path string is identical.
# radix sort is byte-order based (== Unicode code point order for UTF-8),
# matching Python's str ordering.
prefix <- paste0(ROOT, "/")
rel <- ifelse(startsWith(paths, prefix),
              substring(paths, nchar(prefix) + 1L),
              paths)
ord <- order(rel, method = "radix")
paths <- paths[ord]
rel <- rel[ord]

info <- file.info(paths)
sizes <- as.character(info$size) # character: safe for sizes > 2^31, no rounding
hashes <- unname(vapply(paths, sha256_file, character(1)))

df <- data.frame(
  RelativePath = rel,
  SizeBytes = sizes,
  SHA256 = hashes,
  stringsAsFactors = FALSE
)
dir.create(dirname(OUTPUT), recursive = TRUE, showWarnings = FALSE)
# eol = "\r\n" matches Python csv module's default line terminator;
# quote = FALSE matches QUOTE_MINIMAL output (no field here needs quoting).
write.csv(df, OUTPUT, row.names = FALSE, quote = FALSE, eol = "\r\n",
          fileEncoding = "UTF-8")

cat(sprintf("Final files hashed: %d\n", length(paths)))
cat(sprintf("Manifest: %s\n", substring(OUTPUT, nchar(ROOT) + 2L)))
