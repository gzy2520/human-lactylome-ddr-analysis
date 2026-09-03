#!/usr/bin/env Rscript

# Update Figure 3a and Figure 3b regulator percentile heatmaps for the 31-group
# candidate scope by calculating and appending ESCC (PXD064038 and PXD065830)
# percentiles into candidate publication inputs, updating the input manifest,
# and re-rendering high-resolution publication heatmaps.
#
# All biological joins use UniProt BaseAccessions; Gene Symbols are display annotations.

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
  library(dplyr)
  library(ggplot2)
})

set.seed(25)

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1L) {
  normalizePath(args[[1L]], mustWork = TRUE)
} else {
  tryCatch(
    rprojroot::find_root(rprojroot::is_git_root),
    error = function(e) normalizePath(".", mustWork = TRUE)
  )
}
if (!dir.exists(file.path(project_root, "data"))) {
  project_root <- normalizePath(".", mustWork = TRUE)
}

scope_tag <- "escc_inclusion_20260903_pxd065830_tumor_reference"
scope_dir <- file.path(project_root, "data", "candidate", scope_tag)
pub_input_dir <- file.path(scope_dir, "publication_input")
formal_output_dir <- file.path(project_root, "results", scope_tag, "formal_figures")
candidate_heatmap_dir <- file.path(project_root, "results", "candidate", scope_tag, "heatmaps")
dir.create(candidate_heatmap_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(formal_output_dir, recursive = TRUE, showWarnings = FALSE)

stop_if <- function(cond, msg) if (!isTRUE(cond)) stop(msg, call. = FALSE)
base_accession <- function(x) sub("-[0-9]+$", "", sub("^[a-z]+[|]([A-Za-z0-9_]+)[|].*$", "\\1", x))
safe_numeric <- function(v) suppressWarnings(as.numeric(gsub(",", "", as.character(v), fixed = TRUE)))

raw_data_root <- if (dir.exists(file.path(project_root, "data", "PXD065830"))) {
  project_root
} else {
  "/Users/gzy2520/Desktop/Research/kla"
}

# -------------------------------------------------------------------------
# 1. Reference Whole-Proteome Percentiles for PXD065830 (ESCC-T, 94 samples)
# -------------------------------------------------------------------------
message(">>> Computing PXD065830 reference whole-proteome percentiles...")

ref_file <- file.path(raw_data_root, "data", "PXD065830", "supplementary", "Dataset1.xlsx")
stop_if(file.exists(ref_file), paste("Missing reference file:", ref_file))

raw_ref <- as.data.frame(read_excel(ref_file, sheet = "2.a protein raw information", col_names = FALSE))
ref_headers <- trimws(as.character(unlist(raw_ref[2, ], use.names = FALSE)))
t_cols <- which(grepl("^ESCC-[0-9]+T$", ref_headers))
stop_if(length(t_cols) == 94L, paste("Expected 94 ESCC tumor samples; found", length(t_cols)))

row_indices <- seq.int(3L, nrow(raw_ref))
raw_accs <- trimws(as.character(raw_ref[row_indices, 1]))
base_accs <- base_accession(raw_accs)

# Compute within-sample percentile for each of 94 tumor samples
ref_sample_pcts <- list()
for (col_idx in t_cols) {
  s_name <- ref_headers[col_idx]
  vals <- safe_numeric(raw_ref[row_indices, col_idx])
  valid <- !is.na(vals) & vals > 0
  sub_df <- data.frame(
    BaseAccession = base_accs[valid],
    Signal = vals[valid],
    stringsAsFactors = FALSE
  )
  sub_dt <- as.data.table(sub_df)[, .(Signal = sum(Signal)), by = BaseAccession]
  sub_dt[, Log2Signal := log2(Signal + 1)]
  n_prot <- nrow(sub_dt)
  sub_dt[, WithinSamplePercentile := 100 * (rank(Log2Signal, ties.method = "average") - 1) / (n_prot - 1)]
  ref_sample_pcts[[s_name]] <- sub_dt
}

# Template for 48 regulators in reference percentiles table
ref_csv_path <- file.path(pub_input_dir, "regulator_reference_percentiles_30.csv")
dt_ref_orig <- fread(ref_csv_path)

# Extract unique 48 regulators template from an existing group
ref_reg_template <- unique(dt_ref_orig[PXD == "PXD066054" & SampleGroup == "BPH", .(RegulatorBaseAccession, GeneSymbol)])
stop_if(nrow(ref_reg_template) == 48L, "Expected 48 unique regulators in reference table.")

# Calculate summary across 94 samples for each regulator
ref_escc_rows <- list()
for (i in seq_len(nrow(ref_reg_template))) {
  acc <- ref_reg_template$RegulatorBaseAccession[[i]]
  sym <- ref_reg_template$GeneSymbol[[i]]
  
  pcts <- numeric(length(t_cols))
  log2s <- numeric(length(t_cols))
  detected <- logical(length(t_cols))
  
  for (j in seq_along(t_cols)) {
    s_name <- ref_headers[t_cols[j]]
    hit <- ref_sample_pcts[[s_name]][BaseAccession == acc]
    if (nrow(hit) > 0) {
      pcts[j] <- max(hit$WithinSamplePercentile)
      log2s[j] <- max(hit$Log2Signal)
      detected[j] <- TRUE
    } else {
      pcts[j] <- 0
      log2s[j] <- NA_real_
      detected[j] <- FALSE
    }
  }
  
  n_det <- sum(detected)
  med_pct <- median(pcts)
  med_log2 <- if (n_det > 0) median(log2s[detected], na.rm = TRUE) else NA_real_
  
  ref_escc_rows[[i]] <- data.table(
    PXD = "PXD064038",
    SampleGroup = "MEC and NEC ESCC groups",
    RegulatorBaseAccession = acc,
    GeneSymbol = sym,
    WholeProteomeRelativePercentile = med_pct,
    DetectedSampleCount = as.integer(n_det),
    QuantSampleCount = 94L,
    DetectedSampleFraction = n_det / 94,
    MedianLog2SignalDetected = med_log2,
    Measurement = "ordinary iProX T-sample protein raw intensity",
    SourceFile = "data/PXD065830/supplementary/Dataset1.xlsx"
  )
}
ref_escc_dt <- rbindlist(ref_escc_rows)
message("Computed ", nrow(ref_escc_dt), " reference regulator rows for ESCC.")

# -------------------------------------------------------------------------
# 2. Lactylome (Kla) Percentiles for PXD064038 (6 samples: MEC_1-3, NEC_1-3)
# -------------------------------------------------------------------------
message(">>> Computing PXD064038 lactylome (Kla) percentiles...")

kla_file <- file.path(raw_data_root, "data", "PXD064038", "search_results", "extracted_pairing", "txt", "txt", "La (K)Sites.txt")
stop_if(file.exists(kla_file), paste("Missing Kla file:", kla_file))

dt_kla_raw <- fread(kla_file, sep = "\t", fill = TRUE)
if ("Reverse" %in% names(dt_kla_raw)) dt_kla_raw <- dt_kla_raw[is.na(Reverse) | Reverse != "+"]
if ("Potential contaminant" %in% names(dt_kla_raw)) dt_kla_raw <- dt_kla_raw[is.na(`Potential contaminant`) | `Potential contaminant` != "+"]
if ("id" %in% names(dt_kla_raw)) dt_kla_raw <- dt_kla_raw[!is.na(id)]

sample_tokens <- c("MEC_1", "MEC_2", "MEC_3", "NEC_1", "NEC_2", "NEC_3")
stop_if(all(paste("Intensity", sample_tokens) %in% names(dt_kla_raw)), "Missing intensity columns in Kla file.")

kla_csv_path <- file.path(pub_input_dir, "regulator_kla_percentiles_30.csv")
dt_kla_orig <- fread(kla_csv_path)

# Extract 49 unique regulator rows template from an existing group
kla_reg_template <- unique(dt_kla_orig[PXD == "PXD066054" & SampleGroup == "BPH", .(
  Role, GeneSymbol, RegulatorBaseAccession, RegulatorDisplayName, RoleEntryOrder
)])
stop_if(nrow(kla_reg_template) == 49L, "Expected 49 regulator rows in Kla template.")
target_kla_accs <- unique(kla_reg_template$RegulatorBaseAccession)

kla_sample_features <- list()
for (s in sample_tokens) {
  int_col <- paste("Intensity", s)
  prob_col <- paste("Localization prob", s)
  sub_dt <- copy(dt_kla_raw)
  if (prob_col %in% names(sub_dt)) {
    sub_dt <- sub_dt[safe_numeric(get(prob_col)) > 0]
  }
  sub_dt[, Signal := safe_numeric(get(int_col))]
  sub_dt <- sub_dt[!is.na(Signal) & Signal > 0]
  
  sub_dt[, AccList := lapply(strsplit(as.character(Proteins), "[;, ]+"), function(toks) {
    unique(base_accession(toks[nzchar(toks)]))
  })]
  
  sub_dt[, TargetAcc := vapply(AccList, function(toks) {
    hits <- toks[toks %in% target_kla_accs]
    if (length(hits)) hits[[1]] else NA_character_
  }, character(1))]
  
  sub_dt[, FeatureID := vapply(AccList, function(toks) {
    if (length(toks)) paste0("ACC:", toks[[1]]) else NA_character_
  }, character(1))]
  
  sub_dt[, CanonicalFeature := ifelse(!is.na(TargetAcc) & nzchar(TargetAcc), paste0("ACC:", TargetAcc), FeatureID)]
  sub_dt <- sub_dt[!is.na(CanonicalFeature) & nzchar(CanonicalFeature)]
  
  collapsed <- sub_dt[, .(Signal = sum(Signal)), by = .(CanonicalFeature, TargetAcc)]
  collapsed[, QuantSample := s]
  collapsed[, Log2Signal := log2(Signal + 1)]
  n_feat <- nrow(collapsed)
  collapsed[, WithinSamplePercentile := 100 * (rank(Log2Signal, ties.method = "average") - 1) / (n_feat - 1)]
  kla_sample_features[[s]] <- collapsed
}

all_kla_samples_dt <- rbindlist(kla_sample_features)

kla_escc_rows <- list()
for (i in seq_len(nrow(kla_reg_template))) {
  role <- kla_reg_template$Role[[i]]
  sym <- kla_reg_template$GeneSymbol[[i]]
  acc <- kla_reg_template$RegulatorBaseAccession[[i]]
  disp <- kla_reg_template$RegulatorDisplayName[[i]]
  ord <- kla_reg_template$RoleEntryOrder[[i]]
  
  pcts <- numeric(length(sample_tokens))
  log2s <- numeric(length(sample_tokens))
  detected <- logical(length(sample_tokens))
  
  for (j in seq_along(sample_tokens)) {
    s <- sample_tokens[j]
    hit <- all_kla_samples_dt[QuantSample == s & TargetAcc == acc]
    if (nrow(hit) > 0) {
      pcts[j] <- max(hit$WithinSamplePercentile)
      log2s[j] <- max(hit$Log2Signal)
      detected[j] <- TRUE
    } else {
      pcts[j] <- 0
      log2s[j] <- NA_real_
      detected[j] <- FALSE
    }
  }
  
  n_det <- sum(detected)
  med_pct <- median(pcts)
  med_log2 <- if (n_det > 0) median(log2s[detected], na.rm = TRUE) else NA_real_
  
  kla_escc_rows[[i]] <- data.table(
    Role = role,
    GeneSymbol = sym,
    PXD = "PXD064038",
    SampleGroup = "MEC and NEC ESCC groups",
    SampleGroupID = "PXD064038__MEC and NEC ESCC groups",
    RowLabel = "食管鳞癌 · PXD064038",
    BiologicalMaterial = "ESCC tumor tissue",
    GeneLevelAuditStatus = "逐蛋白可审计",
    IdentificationDetected = n_det > 0,
    RegulatorBaseAccession = acc,
    RegulatorDisplayName = disp,
    RoleEntryOrder = as.integer(ord),
    Category = "cancer_tissue",
    CategoryZh = "肿瘤组织",
    CategoryEn = "tumor tissues",
    ComparisonRowOrder = 10L,
    WholeProteomeDisplayRowOrder = 10L,
    RowLabelEn = "ESCC MEC/NEC groups · PXD064038",
    RelativeKlaPercentile = med_pct,
    DetectedSampleCount = as.integer(n_det),
    QuantSampleCount = 6L,
    DetectedSampleFraction = n_det / 6,
    MedianLog2SignalDetected = med_log2,
    MeasurementLevel = "Kla_site",
    QuantField = "MaxQuant site Intensity",
    SourceFile = "data/PXD064038/search_results/extracted_pairing/txt/txt/La (K)Sites.txt",
    QuantificationAvailable = TRUE,
    QuantState = if (n_det > 0) "有定量信号" else "定量样本中未检出"
  )
}
kla_escc_dt <- rbindlist(kla_escc_rows)
message("Computed ", nrow(kla_escc_dt), " Kla regulator rows for ESCC.")

# -------------------------------------------------------------------------
# 3. Update candidate publication_input tables and manifest
# -------------------------------------------------------------------------
message(">>> Updating publication input tables in candidate scope...")

# Remove any existing PXD064038 rows to make this idempotent
clean_ref <- dt_ref_orig[!(PXD == "PXD064038" & SampleGroup == "MEC and NEC ESCC groups")]
updated_ref <- rbind(clean_ref, ref_escc_dt)
fwrite(updated_ref, ref_csv_path, na = "")

clean_kla <- dt_kla_orig[!(PXD == "PXD064038" & SampleGroup == "MEC and NEC ESCC groups")]
updated_kla <- rbind(clean_kla, kla_escc_dt)
fwrite(updated_kla, kla_csv_path, na = "")

# Rebuild INPUT_MANIFEST.csv
manifest_files <- sort(list.files(pub_input_dir, full.names = FALSE, no.. = TRUE))
manifest_files <- setdiff(manifest_files, "INPUT_MANIFEST.csv")
manifest <- data.table(
  File = manifest_files,
  Bytes = vapply(manifest_files, function(fn) file.info(file.path(pub_input_dir, fn))$size, numeric(1)),
  MD5 = vapply(manifest_files, function(fn) digest::digest(file = file.path(pub_input_dir, fn), algo = "md5", serialize = FALSE), character(1))
)
fwrite(manifest, file.path(pub_input_dir, "INPUT_MANIFEST.csv"), na = "")
message("Updated INPUT_MANIFEST.csv with new checksums.")

# -------------------------------------------------------------------------
# 4. Re-render publication figures using build_publication_outputs.R
# -------------------------------------------------------------------------
message(">>> Re-rendering formal figures via build_publication_outputs.R...")

cmd <- sprintf(
  "KLA_PUBLICATION_INPUT=%s KLA_PUBLICATION_OUTPUT=%s KLA_PUBLICATION_EXPECTED_GROUPS=31 KLA_PUBLICATION_CATEGORY_COUNTS='normal_tissue=9;cancer_tissue=3;cancer_cells=12;normal_cells=7' Rscript %s %s",
  shQuote(pub_input_dir),
  shQuote(formal_output_dir),
  shQuote(file.path(project_root, "R", "publication", "build_publication_outputs.R")),
  shQuote(project_root)
)
res <- system(cmd)
stop_if(res == 0, "Failed to run build_publication_outputs.R")

# Copy the updated heatmaps to candidate heatmaps directory
file.copy(
  file.path(formal_output_dir, "Figure_3a_reference_regulator_percentiles.png"),
  file.path(candidate_heatmap_dir, "Figure_3a_reference_regulator_percentiles.png"),
  overwrite = TRUE
)
file.copy(
  file.path(formal_output_dir, "Figure_3a_reference_regulator_percentiles.pdf"),
  file.path(candidate_heatmap_dir, "Figure_3a_reference_regulator_percentiles.pdf"),
  overwrite = TRUE
)
file.copy(
  file.path(formal_output_dir, "Figure_3b_Kla_regulator_percentiles.png"),
  file.path(candidate_heatmap_dir, "Figure_3b_Kla_regulator_percentiles.png"),
  overwrite = TRUE
)
file.copy(
  file.path(formal_output_dir, "Figure_3b_Kla_regulator_percentiles.pdf"),
  file.path(candidate_heatmap_dir, "Figure_3b_Kla_regulator_percentiles.pdf"),
  overwrite = TRUE
)

message("SUCCESS: Figure 3a and Figure 3b updated successfully with 31 datasets.")
