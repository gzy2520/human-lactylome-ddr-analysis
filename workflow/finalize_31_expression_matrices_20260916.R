#!/usr/bin/env Rscript
# Stage 3 of the 31-group expression extraction: rebuild the whole-group manifest, QC and
# README from the per-group objects, so the summaries always describe every group on disk
# rather than only the groups touched by the most recent run.
#
# Usage: finalize_31_expression_matrices_20260916.R <out_dir> <contract_csv>
args <- commandArgs(TRUE)
stopifnot(length(args) == 2L)
out_dir <- normalizePath(args[[1L]], mustWork = TRUE)
contract_path <- normalizePath(args[[2L]], mustWork = TRUE)

# A reference matrix can back several group rows (the three HCT116 rows and the two HK-2
# rows share one matrix each). Such a row is computed once and then read from cache, so the
# group list is taken from the contract rather than from whichever row happened to build it.
contract <- read.csv(contract_path, stringsAsFactors = FALSE, check.names = FALSE)
ref2groups <- tapply(contract$GroupIDs, contract$ReferenceKey,
                     function(x) paste(unique(unlist(strsplit(x, ";", fixed = TRUE))), collapse = ";"))

obj_files <- sort(list.files(file.path(out_dir, "objects"), pattern = "\\.rds$", full.names = TRUE))
stopifnot(length(obj_files) > 0L)

manifest_rows <- list(); qc_rows <- list(); gene_rows <- list()
for (f in obj_files) {
  o <- readRDS(f)
  prov <- o$Provenance
  if (!is.null(ref2groups[[o$ReferenceKey]])) o$GroupIDs <- ref2groups[[o$ReferenceKey]]
  mat_file <- file.path("matrices", paste0(o$ReferenceKey, "_log2tpm.tsv.gz"))
  cnt_file <- file.path("matrices", paste0(o$ReferenceKey, "_counts.tsv.gz"))
  manifest_rows[[o$ReferenceKey]] <- data.frame(
    ReferenceKey = o$ReferenceKey, GroupIDs = o$GroupIDs, GroupLabel = o$GroupLabel,
    SourceKind = prov$SourceKind, SourcePath = prov$SourcePath, SelectorValue = prov$SelectorValue,
    Samples = ncol(o$tpm), Genes = nrow(o$tpm),
    GenesWithCounts = if (is.null(o$counts)) NA_integer_ else nrow(o$counts),
    ValueRoute = prov$NativeValueRoute, IDRule = prov$IDRule, IDNote = prov$IDNote,
    NativeValuesInteger = prov$NativeValuesInteger, SourceFileMD5 = prov$SourceFileMD5,
    SourceNote = if (is.null(prov$Note)) "" else prov$Note,
    Log2TPMFile = if (file.exists(file.path(out_dir, mat_file))) mat_file else NA_character_,
    CountsFile = if (file.exists(file.path(out_dir, cnt_file))) cnt_file else NA_character_,
    stringsAsFactors = FALSE)

  q <- o$qc; q$ReferenceKey <- o$ReferenceKey; q$GroupIDs <- o$GroupIDs
  qc_rows[[o$ReferenceKey]] <- q

  detected <- rowSums(o$tpm >= 1)
  gene_rows[[o$ReferenceKey]] <- data.frame(
    ReferenceKey = o$ReferenceKey, GroupIDs = o$GroupIDs, StableGeneID = rownames(o$tpm),
    MeanLog2TPM = rowMeans(o$log2_tpm), SamplesWithTPM1 = detected,
    ExpressedInAtLeastHalf = detected >= ceiling(ncol(o$tpm) / 2),
    stringsAsFactors = FALSE)
}

manifest <- do.call(rbind, manifest_rows)
qc <- do.call(rbind, qc_rows)
genes <- do.call(rbind, gene_rows)
rownames(manifest) <- NULL; rownames(qc) <- NULL; rownames(genes) <- NULL

write.csv(manifest, file.path(out_dir, "group_manifest.csv"), row.names = FALSE)
write.csv(qc, file.path(out_dir, "group_sample_qc.csv"), row.names = FALSE)
write.csv(genes, file.path(out_dir, "group_gene_summary.csv"), row.names = FALSE)

# one row per group row (31), so a shared matrix is listed once per group it backs
group_rows <- do.call(rbind, lapply(seq_len(nrow(manifest)), function(i) {
  data.frame(GroupID = strsplit(manifest$GroupIDs[i], ";", fixed = TRUE)[[1L]],
             ReferenceKey = manifest$ReferenceKey[i],
             GroupLabel = manifest$GroupLabel[i],
             Samples = manifest$Samples[i], Genes = manifest$Genes[i],
             Log2TPMFile = manifest$Log2TPMFile[i],
             CountsFile = manifest$CountsFile[i],
             stringsAsFactors = FALSE)
}))
group_rows <- group_rows[order(group_rows$GroupID), ]
write.csv(group_rows, file.path(out_dir, "group_index.csv"), row.names = FALSE)
stopifnot(nrow(group_rows) == nrow(contract), !anyDuplicated(group_rows$GroupID))

fmt <- function(x, d = 3L) formatC(x, format = "f", digits = d)
lines <- c(
  "# 31-group RNA expression matrices (server-side extraction, 2026-09-16)", "",
  sprintf("Groups: %d unique reference matrices covering all %d proteome/Kla group rows (three rows share the HCT116 reference, two share the HK-2 reference).",
          nrow(manifest), length(unique(unlist(strsplit(manifest$GroupIDs, ";", fixed = TRUE))))),
  sprintf("Samples: %d in total. Genes per matrix: %d to %d.", nrow(qc), min(manifest$Genes), max(manifest$Genes)),
  "",
  "## What this is",
  "Per-group expression profiles, not a differential-expression comparison. Each of the 31 group",
  "rows was reduced to the expression profile of its own selected, unperturbed samples. No group",
  "was compared against another group here; the cross-tissue comparison (qsmooth) is a separate",
  "downstream step that runs on these matrices.",
  "",
  "## Conventions",
  "* Primary metric: log2(TPM + 0.5). Cleaned integer counts are archived alongside wherever the",
  "  source actually provides counts.",
  "* TPM is computed from counts against one shared merged-exon length table",
  "  (Ensembl release 111, metadata/annotation/human_gene_lengths_ensembl111.tsv), so every",
  "  count-based group is normalised the same way. Matrices built from a source that ships only",
  "  FPKM/RPKM are converted with the same table and are marked as converted in the manifest.",
  "* Row keys are stable Ensembl gene identifiers with the version suffix removed (ENSG...,",
  "  ENSG..._PAR_Y for pseudoautosomal duplicates). Gene symbols are never a join key: the two",
  "  author matrices whose rows are symbols (GSE171750, GSE283812) go through the official NCBI",
  "  Symbol -> GeneID -> Ensembl route, and symbols with more than one official GeneID are dropped.",
  "  A symbol that official nomenclature has since renamed is retried through the HGNC rename",
  "  history (metadata/annotation/hgnc_complete_set.tsv.gz, prev_symbol), and an alias-derived row",
  "  is dropped rather than merged when the Ensembl gene it lands on is already claimed by another",
  "  row of the same matrix. Without that fallback every renamed gene was silently discarded, and",
  "  because the combined gene space is an intersection, one matrix losing a gene removed it from",
  "  every cross-material comparison downstream.",
  "* Columns are named by the source's own sample identifier (GSM, GTEx SAMPID, TCGA barcode).",
  "* Gene sets differ between matrices because the sources were quantified against different",
  "  annotations; intersect on the stable gene ID before combining them.",
  "",
  "## Files",
  "* `matrices/<ReferenceKey>_log2tpm.tsv.gz` - genes x samples, the primary metric.",
  "* `matrices/<ReferenceKey>_counts.tsv.gz` - archived integer counts, where available.",
  "* `group_index.csv` - group row -> reference key -> sample and gene counts.",
  "* `group_manifest.csv` - full provenance per matrix: source file, selector, value route, ID rule,",
  "  mapping yield and source MD5.",
  "* `group_sample_qc.csv` - per sample: library size, TPM total, detected genes, zero fraction and",
  "  within-group pairwise Spearman range.",
  "* `group_gene_summary.csv` - per gene per group: mean log2(TPM+0.5) and detection.",
  "* `objects/<ReferenceKey>.rds` - per-group object with counts, TPM, log2 TPM, QC and provenance.",
  "",
  "## Selection rule applied",
  "Selected samples carry no experimental drug, knockdown, overexpression or infection. DMSO/vehicle",
  "arms are marked as vehicle and are never merged with untreated samples. One group row accepted an",
  "n = 1 descriptive reference (KLA31_17 TALL-104) where no second public unperturbed library exists;",
  "no within-group variance is estimable for that row.",
  "Known limitations that were not resolvable from public metadata are recorded per group in",
  "`group_manifest.csv` (`SourceNote` and `IDNote`) and in the project ledger.",
  "",
  "## Session",
  ""
)
writeLines(lines, file.path(out_dir, "README.md"))
writeLines(trimws(capture.output(sessionInfo()), which = "right"), file.path(out_dir, "sessionInfo.txt"))
cat("FINALIZE_DONE matrices=", nrow(manifest), " samples=", nrow(qc), "\n", sep = "")
