#!/usr/bin/env Rscript
# Build the publication-deliverable SHA256 manifest.
#
# Usage: Rscript build_final_manifest.R <project_root>
suppressMessages(library(digest))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript build_final_manifest.R <project_root>", call. = FALSE)
}
ROOT <- normalizePath(args[1], winslash = "/", mustWork = TRUE)
OUTPUT <- file.path(
  ROOT, "results", "reports", "publication_manifest_sha256.csv"
)

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
for (f in c(
  "README.md",
  "PROJECT_INDEX.md",
  "NEW_CHAT_PROJECT_PROMPT.md",
  "DESCRIPTION",
  "requirements.txt"
)) {
  p <- file.path(ROOT, f)
  if (file.exists(p)) paths <- c(paths, p)
}

# Version-controlled publication source.
for (dir in c(
  "R", "python", "config", "workflow", "tests", "manuscript", "docs"
)) {
  paths <- c(paths, collect_tree(file.path(ROOT, dir)))
}

# Small versioned scientific inputs; raw PXD files remain outside the manifest.
for (f in c(
  "data/annotations/GO-repair+damage(human).tsv",
  "data/identifier/260810乳酸化DDR基因评分表.xlsx",
  "data/identifier/乳酸化调控因子_Writer-Eraser-Reader.xlsx"
)) {
  p <- file.path(ROOT, f)
  if (file.exists(p)) paths <- c(paths, p)
}

# Publication results only. Historical exploratory results are intentionally
# excluded even if they remain available in a local archive.
result_trees <- c(
  "results/figures/bp_semantic_umap",
  "results/figures/five_set_embeddings",
  "results/figures/five_set_pathway_matrix",
  "results/figures/four_class_venn",
  "results/figures/pathway_specific_umap",
  "results/tables/bp_semantic_umap",
  "results/tables/five_set_embeddings",
  "results/tables/five_set_pathway_matrix",
  "results/tables/four_class_venn",
  "results/tables/pathway_specific_umap",
  "results/tables/protein_function_inputs",
  "results/provenance"
)
for (dir in result_trees) {
  paths <- c(paths, collect_tree(file.path(ROOT, dir)))
}
result_files <- c(
  "results/figures/cell_type_kla_vs_reference_ddr_fraction_en.pdf",
  "results/figures/cell_type_kla_vs_reference_ddr_fraction_en.png",
  "results/figures/cell_type_kla_vs_reference_ddr_fraction_zh.pdf",
  "results/figures/cell_type_kla_vs_reference_ddr_fraction_zh.png",
  "results/figures/kla_regulator_cross_study_relative_intensity_heatmap_en.pdf",
  "results/figures/kla_regulator_cross_study_relative_intensity_heatmap_en.png",
  "results/figures/kla_regulator_cross_study_relative_intensity_heatmap_zh.pdf",
  "results/figures/kla_regulator_cross_study_relative_intensity_heatmap_zh.png",
  "results/figures/kla_regulator_whole_proteome_relative_intensity_heatmap_en.pdf",
  "results/figures/kla_regulator_whole_proteome_relative_intensity_heatmap_en.png",
  "results/figures/kla_regulator_whole_proteome_relative_intensity_heatmap_zh.pdf",
  "results/figures/kla_regulator_whole_proteome_relative_intensity_heatmap_zh.png",
  "results/tables/cell_type_kla_vs_reference_ddr_statistics_accession_only.csv",
  "results/tables/cell_type_kla_vs_reference_ddr_statistics_accession_only_paired_33.csv",
  "results/tables/cell_type_kla_vs_reference_ddr_statistics_accession_only_paired_33_zh.csv",
  "results/tables/kla_and_reference_teacher_review_zh.csv",
  "results/tables/kla_regulator_heatmap_axis_order.csv",
  "results/tables/kla_regulator_intensity_availability_audit.csv",
  "results/tables/kla_regulator_whole_proteome_heatmap_rows.csv",
  "results/tables/strict_reference_material_identity_audit.csv",
  "results/tables/strict_reference_material_identity_audit_zh.csv",
  "results/reports/FIVE_SET_UMAP_TSNE_PCA_PATHWAY_GRIDS_33GROUP_V3.md",
  "results/reports/UMAP_PATHWAY_PIE_33GROUP_V4_BP_SEMANTIC.md",
  "results/reports/UMAP_PATHWAY_SPECIFIC_33GROUP_V1.md"
)
for (f in result_files) {
  p <- file.path(ROOT, f)
  if (file.exists(p)) paths <- c(paths, p)
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
# Use repository-native LF endings so git whitespace checks remain clean.
# quote = FALSE matches QUOTE_MINIMAL output (no field here needs quoting).
write.csv(df, OUTPUT, row.names = FALSE, quote = FALSE, eol = "\n",
          fileEncoding = "UTF-8")

cat(sprintf("Final files hashed: %d\n", length(paths)))
cat(sprintf("Manifest: %s\n", substring(OUTPUT, nchar(ROOT) + 2L)))
