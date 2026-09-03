#!/usr/bin/env Rscript

# Build an isolated 2026-09-03 scope that adds the representative ESCC
# lactylome group PXD064038. The frozen 30-group publication input is never
# edited. The six MEC/NEC samples are represented by one publication group,
# while their source-defined sample observations are retained for candidate
# plots. The reference is the ordinary, non-lactylated ESCC tumor arm of
# iProX/PXD065830; the 24 non-tumor samples in that project are excluded.

suppressPackageStartupMessages({
  library(data.table)
  library(openxlsx2)
  library(readxl)
})

set.seed(25)

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1L) {
  normalizePath(args[[1L]], mustWork = TRUE)
} else {
  normalizePath(".", mustWork = TRUE)
}
source_root <- if (length(args) >= 2L) {
  normalizePath(args[[2L]], mustWork = TRUE)
} else {
  stop(
    "Usage: Rscript R/candidate/prepare_escc_inclusion_inputs.R <project-root> <source-data-root>",
    call. = FALSE
  )
}

stop_if <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

require_file <- function(path) {
  stop_if(file.exists(path), paste0("Required source file is missing: ", path))
  path
}

base_accession <- function(values) {
  values <- trimws(as.character(values))
  values[is.na(values)] <- ""
  values <- sub("^.*\\|([^|]+)\\|.*$", "\\1", values)
  values <- sub("^([^|;]+)\\|.*$", "\\1", values)
  values <- sub("^NX_", "", values)
  sub("-[0-9]+$", "", values)
}

is_uniprot <- function(values) {
  grepl(
    "^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9](?:[A-Z0-9]{3}[0-9])?)$",
    values
  )
}

split_accessions <- function(values) {
  tokens <- unlist(strsplit(as.character(values), "[;,]"), use.names = FALSE)
  accessions <- base_accession(tokens)
  sort(unique(accessions[is_uniprot(accessions)]))
}

is_true_flag <- function(values) {
  values <- tolower(trimws(as.character(values)))
  values[is.na(values)] <- ""
  values %in% c("+", "1", "true", "yes", "y")
}

positive_numeric <- function(values) {
  values <- suppressWarnings(as.numeric(values))
  is.finite(values) & values > 0
}

valid_maxquant_rows <- function(data) {
  keep <- rep(TRUE, nrow(data))
  for (column in c("Reverse", "Potential contaminant", "Contaminant", "Only identified by site")) {
    if (column %in% names(data)) keep <- keep & !is_true_flag(data[[column]])
  }
  keep[is.na(keep)] <- FALSE
  keep
}

first_nonempty <- function(values) {
  values <- trimws(unlist(strsplit(as.character(values), ";"), use.names = FALSE))
  values <- values[!is.na(values) & nzchar(values)]
  if (length(values)) values[[1L]] else NA_character_
}

extract_pxd065830_tumor_reference <- function(path) {
  data <- read_excel(
    path,
    sheet = "2.a protein raw information",
    col_names = FALSE
  )
  stop_if(nrow(data) >= 3L && ncol(data) >= 14L,
    "PXD065830 Dataset1 protein sheet is unexpectedly small.")
  headers <- trimws(as.character(unlist(data[2, ], use.names = FALSE)))
  tumor_columns <- which(grepl("^ESCC-[0-9]+T$", headers))
  stop_if(length(tumor_columns) == 94L,
    paste0("PXD065830 must contribute exactly 94 ESCC tumor T columns; found ", length(tumor_columns), "."))
  row_indices <- seq.int(3L, nrow(data))
  present_matrix <- do.call(cbind, lapply(tumor_columns, function(column) {
    positive_numeric(data[[column]][row_indices])
  }))
  keep <- rowSums(present_matrix, na.rm = TRUE) > 0L
  accessions <- split_accessions(data[[1L]][row_indices][keep])
  stop_if(length(accessions) == 8083L,
    paste0("PXD065830 ESCC tumor union must contain 8083 BaseAccessions; found ", length(accessions), "."))
  stop_if(!anyDuplicated(accessions), "PXD065830 ESCC tumor reference contains duplicate BaseAccessions.")

  display <- rbindlist(lapply(row_indices[keep], function(index) {
    ids <- split_accessions(data[[1L]][[index]])
    if (!length(ids)) return(NULL)
    data.table(
      BaseAccession = ids,
      GeneSymbolAudit = first_nonempty(data[[3L]][[index]]),
      ProteinNameAudit = first_nonempty(data[[2L]][[index]])
    )
  }), fill = TRUE)
  display <- display[, .(
    GeneSymbolAudit = first_nonempty(GeneSymbolAudit),
    ProteinNameAudit = first_nonempty(ProteinNameAudit)
  ), by = BaseAccession]

  # Preserve each deposited ESCC-T column as an individual source observation
  # for the Figure 1 sample-point layer.  Membership itself remains keyed by
  # stable UniProt BaseAccessions; no Gene Symbol is used for set operations.
  sample_accessions <- setNames(lapply(seq_along(tumor_columns), function(index) {
    column <- tumor_columns[[index]]
    sample_keep <- positive_numeric(data[[column]][row_indices])
    split_accessions(data[[1L]][row_indices][sample_keep])
  }), headers[tumor_columns])
  stop_if(all(vapply(sample_accessions, length, integer(1)) > 0L),
    "At least one PXD065830 ESCC-T source column has no positive BaseAccessions.")

  list(
    accessions = accessions,
    sample_ids = headers[tumor_columns],
    sample_accessions = sample_accessions,
    protein_table_rows = sum(keep),
    display = display
  )
}

scope_tag <- Sys.getenv(
  "KLA_ESCC_SCOPE_TAG",
  unset = "escc_inclusion_20260903_pxd065830_tumor_reference"
)
scope_dir <- file.path(project_root, "data", "candidate", scope_tag)
publication_input_dir <- file.path(scope_dir, "publication_input")
candidate_input_dir <- file.path(scope_dir, "candidate_input")
dir.create(publication_input_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(candidate_input_dir, recursive = TRUE, showWarnings = FALSE)

current_publication_dir <- file.path(project_root, "data", "publication_input")
current_candidate_dir <- file.path(project_root, "data", "candidate")
current_groups <- fread(file.path(current_publication_dir, "group_summary_30.csv"), check.names = FALSE)
current_kla_membership <- fread(file.path(current_publication_dir, "kla_protein_membership_30.csv"), check.names = FALSE)
current_reference_membership <- fread(file.path(current_publication_dir, "reference_protein_membership_30.csv"), check.names = FALSE)
current_all_venn <- fread(file.path(current_publication_dir, "venn_all_kla.csv"), check.names = FALSE)
current_kla_ddr_venn <- fread(file.path(current_publication_dir, "venn_kla_ddr.csv"), check.names = FALSE)
current_reference_venn <- fread(file.path(current_publication_dir, "venn_reference.csv"), check.names = FALSE)
current_s4_path <- require_file(file.path(
  current_publication_dir,
  "Supplementary_Table_S4_Pathway_Protein_Ranking.xlsx"
))

source_file_rel <- file.path(
  "data", "PXD064038", "search_results", "extracted_pairing", "txt", "txt", "La (K)Sites.txt"
)
source_file <- require_file(file.path(source_root, source_file_rel))
clinical_file_rel <- file.path("data", "PXD064038", "metadata", "Clinical information of samples.docx")
clinical_file <- file.path(source_root, clinical_file_rel)
reference_pxd <- "PXD065830"
reference_file_rel <- file.path("data", reference_pxd, "supplementary", "Dataset1.xlsx")
reference_file <- require_file(file.path(source_root, reference_file_rel))

raw <- fread(source_file, check.names = FALSE, showProgress = FALSE)
required_raw_columns <- c("Proteins", "Reverse", "Potential contaminant", "id")
stop_if(all(required_raw_columns %in% names(raw)), "PXD064038 La (K)Sites.txt is missing required MaxQuant columns.")
sample_ids <- c("MEC_1", "MEC_2", "MEC_3", "NEC_1", "NEC_2", "NEC_3")
localization_columns <- intersect(paste("Localization prob", sample_ids), names(raw))
stop_if(length(localization_columns) == length(sample_ids), "PXD064038 does not contain all six MEC/NEC localization columns.")

base_keep <- valid_maxquant_rows(raw) & !is.na(raw$id)
sample_records <- rbindlist(lapply(sample_ids, function(sample_id) {
  keep <- base_keep & positive_numeric(raw[[paste("Localization prob", sample_id)]])
  data.table(
    SampleID = sample_id,
    BaseAccession = split_accessions(raw$Proteins[keep])
  )
}), fill = TRUE)
sample_records <- unique(sample_records)
stop_if(all(sample_records[, .N, by = SampleID]$N > 0L), "At least one PXD064038 sample has no localized Kla proteins.")

escc_accessions <- sort(unique(sample_records$BaseAccession))
stop_if(length(escc_accessions) > 0L, "PXD064038 produced no valid UniProt BaseAccessions.")

go_annotations <- fread(
  file.path(current_publication_dir, "human_ddr_go_annotations.tsv"),
  sep = "\t",
  quote = "",
  encoding = "UTF-8",
  check.names = FALSE
)
go_keep <- (
  is.na(go_annotations[["TAXON ID"]]) | go_annotations[["TAXON ID"]] == 9606
) &
  !grepl("(^|\\|)NOT($|\\|)", go_annotations$QUALIFIER) &
  go_annotations[["GENE PRODUCT DB"]] == "UniProtKB"
ddr_accessions <- sort(unique(base_accession(go_annotations[["GENE PRODUCT ID"]][go_keep])))
ddr_accessions <- ddr_accessions[is_uniprot(ddr_accessions)]
sample_records[, IsDdr := BaseAccession %in% ddr_accessions]
escc_ddr_accessions <- sort(intersect(escc_accessions, ddr_accessions))

reference_result <- extract_pxd065830_tumor_reference(reference_file)
reference_accessions <- reference_result$accessions
reference_ddr_accessions <- sort(intersect(reference_accessions, ddr_accessions))
reference_fraction <- length(reference_ddr_accessions) / length(reference_accessions)
reference_sample_counts <- rbindlist(lapply(reference_result$sample_ids, function(sample_id) {
  sample_accessions <- reference_result$sample_accessions[[sample_id]]
  ddr_sample_accessions <- intersect(sample_accessions, ddr_accessions)
  data.table(
    SampleID = sample_id,
    ProteinCount = length(sample_accessions),
    DdrProteinCount = length(ddr_sample_accessions),
    DdrFraction = length(ddr_sample_accessions) / length(sample_accessions)
  )
}))
stop_if(nrow(reference_sample_counts) == 94L &&
          all(is.finite(reference_sample_counts$DdrFraction)) &&
          all(reference_sample_counts$DdrFraction >= 0 & reference_sample_counts$DdrFraction <= 1),
  "PXD065830 ESCC-T per-sample DDR fractions are invalid.")

sample_counts <- sample_records[, .(
  KlaProteinCount = uniqueN(BaseAccession),
  KlaDdrProteinCount = uniqueN(BaseAccession[IsDdr]),
  KlaDdrFraction = uniqueN(BaseAccession[IsDdr]) / uniqueN(BaseAccession)
), by = SampleID]
sample_counts[, SampleOrder := match(SampleID, sample_ids)]
setorder(sample_counts, SampleOrder)
sample_counts[, SampleOrder := NULL]

group_pxd <- "PXD064038"
group_name <- "MEC and NEC ESCC groups"
group_category <- "cancer_tissue"
group_row_order <- 10L
group_label_en <- "ESCC MEC/NEC groups"
group_label_zh <- "ESCC MEC/NEC组"
group_reference_note <- paste(
  "External ESCC tumor ordinary-proteome reference from iProX/PXD065830 Dataset1.xlsx, sheet 2.a protein raw information:",
  "only the 94 ESCC T samples are used; the 24 N non-tumor samples are excluded.",
  paste0("The T union contains ", length(reference_accessions), " BaseAccessions and ", length(reference_ddr_accessions), " DDR BaseAccessions."),
  "This is an independent ESCC tumor cohort, not a same-specimen paired measurement; comparison is by BaseAccession.",
  "PXD053809 is not selected because its deposited processed XLS is encrypted and is not currently auditable; PXD010154 remains healthy-esophagus background only."
)

new_group <- data.table(
  PXD = group_pxd,
  SampleGroup = group_name,
  DisplayLabel = paste0(group_label_en, " · ", group_pxd),
  BiologicalMaterial = "ESCC tumor tissue",
  KlaEvidenceFile = source_file_rel,
  KlaProteinCount = length(escc_accessions),
  KlaDdrProteinCount = length(escc_ddr_accessions),
  KlaDdrFraction = length(escc_ddr_accessions) / length(escc_accessions),
  PairedAnalysisIncluded = TRUE,
  ReferencePXD = reference_pxd,
  ReferenceEvidenceFile = reference_file_rel,
  ReferenceProteinCount = length(reference_accessions),
  ReferenceDdrProteinCount = length(reference_ddr_accessions),
  ReferenceDdrFraction = reference_fraction,
  DdrFractionPercentagePointDifference = length(escc_ddr_accessions) / length(escc_accessions) - reference_fraction,
  ReferenceMatchNote = group_reference_note,
  MatchMode = "BaseAccession_external_ESCC_tumor",
  SymbolFallbackCount = 0L,
  RowOrder = group_row_order,
  KlaLabelEn = group_label_en,
  KlaLabelZh = group_label_zh,
  ReferenceLabelEn = "Independent ESCC tumor whole proteome",
  ReferenceLabelZh = "独立ESCC肿瘤非乳酸化全蛋白组",
  NamingBasis = "TCR paper defines MEC/NEC as two ESCC tumor subgroups; combined as one representative ESCC tissue group",
  Category = group_category,
  CategoryZh = "肿瘤组织",
  CategoryEn = "tumor tissues",
  EnglishDisplayLabel = paste0(group_label_en, " · ", group_pxd)
)

stop_if(!anyDuplicated(current_groups[, .(PXD, SampleGroup)]), "Current publication groups are not unique.")
stop_if(!anyDuplicated(new_group[, .(PXD, SampleGroup)]), "The ESCC group key is duplicated.")
expanded_groups <- rbindlist(list(current_groups, new_group), fill = TRUE)
setorder(expanded_groups, RowOrder)
stop_if(nrow(expanded_groups) == nrow(current_groups) + 1L, "Expanded group summary did not add exactly one ESCC group.")

new_kla_membership <- data.table(
  PXD = group_pxd,
  SampleGroup = group_name,
  BaseAccession = escc_accessions,
  IsDdr = escc_accessions %in% ddr_accessions
)
expanded_kla_membership <- unique(rbindlist(list(current_kla_membership, new_kla_membership), fill = TRUE))
stop_if(!anyDuplicated(expanded_kla_membership[, .(PXD, SampleGroup, BaseAccession)]), "Expanded Kla membership has duplicate group/accession rows.")

new_reference_membership <- data.table(
  PXD = group_pxd,
  SampleGroup = group_name,
  ReferencePXD = reference_pxd,
  SourceProteinID = reference_accessions,
  IdentifierType = "UniProtKB",
  MappedBaseAccessions = reference_accessions,
  IsDdr = reference_accessions %in% ddr_accessions
)
expanded_reference_membership <- unique(rbindlist(list(current_reference_membership, new_reference_membership), fill = TRUE))
stop_if(!anyDuplicated(expanded_reference_membership[, .(PXD, SampleGroup, SourceProteinID)]),
  "Expanded reference membership has duplicate group/accession rows.")

display_from_raw <- rbindlist(lapply(seq_len(nrow(raw)), function(index) {
  accessions <- split_accessions(raw$Proteins[[index]])
  if (!length(accessions)) return(NULL)
  data.table(
    BaseAccession = accessions,
    GeneSymbolAudit = first_nonempty(raw[["Gene names"]][[index]]),
    ProteinNameAudit = first_nonempty(raw[["Protein names"]][[index]])
  )
}), fill = TRUE)
display_from_raw <- display_from_raw[, .(
  GeneSymbolAudit = first_nonempty(GeneSymbolAudit),
  ProteinNameAudit = first_nonempty(ProteinNameAudit)
), by = BaseAccession]
new_display <- display_from_raw[BaseAccession %in% escc_accessions]
new_display[, `:=`(
  ReviewedStatus = "source_display_audit_only",
  GeneSymbolAuditSource = source_file_rel,
  AnnotationMappingSource = ifelse(
    BaseAccession %in% escc_ddr_accessions,
    paste(source_file_rel, "human_ddr_go_annotations.tsv", sep = "; "),
    source_file_rel
  )
)]
display_columns <- c(
  "BaseAccession", "GeneSymbolAudit", "ProteinNameAudit", "ReviewedStatus",
  "GeneSymbolAuditSource", "AnnotationMappingSource"
)
new_reference_display <- copy(reference_result$display)
new_reference_display[, `:=`(
  ReviewedStatus = "source_display_audit_only",
  GeneSymbolAuditSource = reference_file_rel,
  AnnotationMappingSource = ifelse(
    BaseAccession %in% reference_ddr_accessions,
    paste(reference_file_rel, "human_ddr_go_annotations.tsv", sep = "; "),
    reference_file_rel
  )
)]
current_display <- unique(rbindlist(list(
  current_all_venn[, ..display_columns],
  current_reference_venn[, ..display_columns]
), fill = TRUE), by = "BaseAccession")
display <- unique(rbindlist(list(
  current_display,
  new_display[, ..display_columns],
  new_reference_display[, ..display_columns]
), fill = TRUE), by = "BaseAccession")

venn_categories <- c("normal_tissue", "cancer_tissue", "cancer_cells", "normal_cells")
region_names <- vapply(seq_len(15L), function(mask) {
  members <- venn_categories[bitwAnd(mask, bitwShiftL(1L, seq_along(venn_categories) - 1L)) != 0L]
  if (length(members) == length(venn_categories)) "all_four" else paste0(paste(members, collapse = "_and_"), "_only")
}, character(1))

expand_reference_membership_for_venn <- function(membership) {
  membership <- membership[!is.na(MappedBaseAccessions) & nzchar(trimws(MappedBaseAccessions))]
  membership[, .(
    BaseAccession = trimws(unlist(strsplit(MappedBaseAccessions, ";", fixed = TRUE)))
  ), by = .(PXD, SampleGroup, IsDdr)]
}

build_venn <- function(membership, ddr_only = FALSE) {
  if (ddr_only) membership <- membership[is_true_flag(IsDdr)]
  category_lookup <- expanded_groups[, .(PXD, SampleGroup, Category)]
  membership <- unique(merge(membership, category_lookup, by = c("PXD", "SampleGroup"), all = FALSE)[, .(BaseAccession, Category)])
  flags <- dcast(
    membership[, Present := TRUE],
    BaseAccession ~ Category,
    value.var = "Present",
    fill = FALSE
  )
  for (category in venn_categories) {
    if (!category %in% names(flags)) flags[, (category) := FALSE]
  }
  flags <- flags[, c("BaseAccession", venn_categories), with = FALSE]
  setnames(flags, venn_categories, paste0("In_", venn_categories))
  flag_matrix <- as.matrix(flags[, paste0("In_", venn_categories), with = FALSE])
  masks <- as.integer(flag_matrix %*% (2^(seq_along(venn_categories) - 1L)))
  flags[, Region := region_names[masks]]
  output <- merge(display, flags, by = "BaseAccession", all.y = TRUE, sort = FALSE)
  setcolorder(output, c(display_columns, paste0("In_", venn_categories), "Region"))
  setorder(output, BaseAccession)
  output
}

expanded_reference_venn_membership <- expand_reference_membership_for_venn(expanded_reference_membership)
expanded_all_venn <- build_venn(expanded_kla_membership, FALSE)
expanded_kla_ddr_venn <- build_venn(expanded_kla_membership, TRUE)
expanded_reference_venn <- build_venn(expanded_reference_venn_membership, FALSE)
expanded_reference_ddr_venn <- build_venn(expanded_reference_venn_membership, TRUE)
stop_if(nrow(expanded_all_venn) == uniqueN(expanded_kla_membership$BaseAccession), "Expanded all-Kla union has an unexpected size.")
stop_if(nrow(expanded_kla_ddr_venn) == uniqueN(expanded_kla_membership[is_true_flag(IsDdr), BaseAccession]), "Expanded Kla-DDR union has an unexpected size.")
stop_if(nrow(expanded_reference_venn) == uniqueN(expanded_reference_venn_membership$BaseAccession), "Expanded reference union has an unexpected size.")
stop_if(nrow(expanded_reference_ddr_venn) == uniqueN(expanded_reference_venn_membership[is_true_flag(IsDdr), BaseAccession]), "Expanded reference-DDR union has an unexpected size.")

write_data_workbook <- function(sheets, output_path) {
  workbook <- wb_workbook(creator = "klaPublicationAnalysis-escc-inclusion")
  for (sheet_name in names(sheets)) {
    data <- as.data.frame(sheets[[sheet_name]])
    table_name <- paste0("S4", gsub("[^A-Za-z0-9]", "", sheet_name))
    workbook <- wb_add_worksheet(workbook, sheet = sheet_name, grid_lines = FALSE)
    workbook <- wb_add_data_table(
      workbook,
      sheet = sheet_name,
      x = data,
      table_name = table_name,
      table_style = "TableStyleLight1",
      banded_rows = FALSE,
      with_filter = TRUE,
      inline_strings = FALSE
    )
    workbook <- wb_freeze_pane(workbook, sheet = sheet_name, first_row = TRUE)
  }
  wb_save(workbook, output_path, overwrite = TRUE)
}

pathway_order <- c("BER", "NER", "MMR", "FA", "HR", "AEJ", "NHEJ")
s4_sheets <- lapply(excel_sheets(current_s4_path), function(sheet_name) {
  as.data.table(read_excel(current_s4_path, sheet = sheet_name))
})
names(s4_sheets) <- excel_sheets(current_s4_path)
tumor_sheet <- s4_sheets[["TumorTissues"]]
pathway_columns <- c("BER", "NER", "MMR", "FA", "HR", "AEJ", "NHEJ")
s4_columns <- c("BaseAccession", "GeneSymbol", "ProteinName", pathway_columns, "SignedScore", "Note")
all_s4_rows <- rbindlist(lapply(s4_sheets, function(sheet) sheet[, ..s4_columns]), fill = TRUE)
all_s4_rows <- all_s4_rows[!duplicated(BaseAccession)]
new_tumor_ids <- setdiff(
  expanded_kla_ddr_venn[In_cancer_tissue == TRUE, BaseAccession],
  tumor_sheet$BaseAccession
)
new_s4_rows <- all_s4_rows[BaseAccession %in% new_tumor_ids, ..s4_columns]
missing_curated_rows <- setdiff(new_tumor_ids, new_s4_rows$BaseAccession)
if (length(missing_curated_rows)) {
  fallback_display <- display[match(missing_curated_rows, BaseAccession)]
  stop_if(all(!is.na(fallback_display$GeneSymbolAudit) & nzchar(fallback_display$GeneSymbolAudit)),
    "Missing display audit names for new tumor S4 rows.")
  fallback_rows <- data.table(
    BaseAccession = missing_curated_rows,
    GeneSymbol = fallback_display$GeneSymbolAudit,
    ProteinName = fallback_display$ProteinNameAudit,
    BER = 0, NER = 0, MMR = 0, FA = 0, HR = 0, AEJ = 0, NHEJ = 0,
    SignedScore = 0,
    Note = "Broad DNA-repair/DNA-damage-response GO annotation; no specific seven-pathway seed hit"
  )
  new_s4_rows <- rbindlist(list(new_s4_rows, fallback_rows), fill = TRUE)
}
new_s4_rows <- new_s4_rows[match(new_tumor_ids, BaseAccession), ..s4_columns]
s4_sheets[["TumorTissues"]] <- rbind(tumor_sheet, new_s4_rows, fill = TRUE)
setorder(s4_sheets[["TumorTissues"]], SignedScore, BaseAccession)
stop_if(nrow(s4_sheets[["TumorTissues"]]) == nrow(tumor_sheet) + length(new_tumor_ids),
  "Expanded tumor S4 did not add all new cancer-tissue Kla-DDR proteins.")

current_publication_files <- list.files(current_publication_dir, full.names = FALSE, no.. = TRUE)
for (filename in setdiff(current_publication_files, "INPUT_MANIFEST.csv")) {
  stop_if(file.copy(file.path(current_publication_dir, filename), file.path(publication_input_dir, filename), overwrite = TRUE),
    paste0("Could not copy publication input: ", filename))
}
fwrite(expanded_groups, file.path(publication_input_dir, "group_summary_30.csv"), na = "")
fwrite(expanded_kla_membership, file.path(publication_input_dir, "kla_protein_membership_30.csv"), na = "")
fwrite(expanded_reference_membership, file.path(publication_input_dir, "reference_protein_membership_30.csv"), na = "")
fwrite(expanded_all_venn, file.path(publication_input_dir, "venn_all_kla.csv"), na = "")
fwrite(expanded_kla_ddr_venn, file.path(publication_input_dir, "venn_kla_ddr.csv"), na = "")
fwrite(expanded_reference_venn, file.path(publication_input_dir, "venn_reference.csv"), na = "")
fwrite(expanded_reference_ddr_venn, file.path(publication_input_dir, "venn_reference_ddr.csv"), na = "")
write_data_workbook(s4_sheets, file.path(publication_input_dir, "Supplementary_Table_S4_Pathway_Protein_Ranking.xlsx"))

manifest_files <- sort(list.files(publication_input_dir, full.names = FALSE, no.. = TRUE))
manifest_files <- setdiff(manifest_files, "INPUT_MANIFEST.csv")
manifest <- data.table(
  File = manifest_files,
  Bytes = vapply(manifest_files, function(filename) file.info(file.path(publication_input_dir, filename))$size, numeric(1)),
  MD5 = vapply(manifest_files, function(filename) digest::digest(file = file.path(publication_input_dir, filename), algo = "md5", serialize = FALSE), character(1))
)
fwrite(manifest, file.path(publication_input_dir, "INPUT_MANIFEST.csv"), na = "")

candidate_files <- c(
  "figure1_sample_boxplot_values.csv",
  "figure1_sample_boxplot_source_registry.csv",
  "figure1_pathway_summary_sample_boxplot_values.csv",
  "sample_boxplot_values.csv",
  "sample_boxplot_reconciliation.csv",
  "sample_boxplot_source_registry.csv",
  "biological_sample_count_record.csv",
  "sample_design_30.csv"
)
for (filename in candidate_files) {
  stop_if(file.copy(file.path(current_candidate_dir, filename), file.path(candidate_input_dir, filename), overwrite = TRUE),
    paste0("Could not copy candidate input: ", filename))
}

sample_common <- merge(
  sample_records,
  sample_counts,
  by = "SampleID",
  all.x = TRUE,
  sort = FALSE
)
sample_common[, `:=`(
  RowOrder = group_row_order,
  PXD = group_pxd,
  SampleGroup = group_name,
  Category = group_category,
  DisplayGroup = group_label_en,
  ConditionLabel = SampleID,
  SampleClass = fifelse(grepl("^MEC", SampleID), "ESCC with lymph-node metastasis", "ESCC without lymph-node metastasis"),
  ObservationType = "sample",
  SourceMode = "deposited_sample_table",
  SourceFile = source_file_rel,
  ReferenceFraction = reference_fraction * 100,
  ReferenceDdr = length(reference_ddr_accessions),
  ReferenceTotal = length(reference_accessions),
  KlaDdrFractionPercentage = KlaDdrFraction * 100
)]
setcolorder(sample_common, c(
  "RowOrder", "PXD", "SampleGroup", "Category", "DisplayGroup", "SampleID", "ConditionLabel",
  "SampleClass", "ObservationType", "SourceMode", "SourceFile", "KlaProteinCount",
  "KlaDdrProteinCount", "KlaDdrFraction", "KlaDdrFractionPercentage", "ReferenceFraction",
  "ReferenceDdr", "ReferenceTotal", "BaseAccession", "IsDdr"
))

sample_meta <- unique(sample_common[, !c("BaseAccession", "IsDdr"), with = FALSE])
stop_if(nrow(sample_meta) == length(sample_ids), "ESCC sample metadata did not collapse to six source observations.")

new_sample_boxplot_values <- sample_meta[, .(
  RowOrder, PXD, SampleGroup, Category, DisplayGroup, SampleID, ConditionLabel,
  SampleClass, ObservationType, SourceMode, SourceFile, KlaProteinCount,
  KlaDdrProteinCount, KlaDdrFraction, KlaDdrFractionPercentage, ReferenceFraction,
  ReferenceDdr, ReferenceTotal
)]
existing_sample_boxplot_values <- fread(file.path(candidate_input_dir, "sample_boxplot_values.csv"), check.names = FALSE)
updated_sample_boxplot_values <- rbindlist(list(existing_sample_boxplot_values, new_sample_boxplot_values), fill = TRUE)
setorder(updated_sample_boxplot_values, RowOrder, PXD, SampleGroup, SampleID)
fwrite(updated_sample_boxplot_values, file.path(candidate_input_dir, "sample_boxplot_values.csv"), na = "")

new_figure1_values <- sample_meta[, .(
  RowOrder, PXD, SampleGroup, Category, DisplayGroup,
  Dataset = "Lactylome (Kla)", SampleID, ConditionLabel, SampleClass, ObservationType,
  SourceMode, SourceFile, ReferencePXD = reference_pxd, ProteinCount = KlaProteinCount,
  DdrProteinCount = KlaDdrProteinCount, DdrFraction = KlaDdrFraction,
  DdrFractionPercentage = KlaDdrFractionPercentage,
  FrozenKlaProteinCount = length(escc_accessions),
  FrozenKlaDdrProteinCount = length(escc_ddr_accessions),
  FrozenKlaDdrFraction = length(escc_ddr_accessions) / length(escc_accessions) * 100,
  FrozenReferenceProteinCount = length(reference_accessions),
  FrozenReferenceDdrProteinCount = length(reference_ddr_accessions),
  FrozenReferenceDdrFraction = reference_fraction * 100
)]
new_figure1_reference_values <- reference_sample_counts[, .(
  RowOrder = group_row_order,
  PXD = group_pxd,
  SampleGroup = group_name,
  Category = group_category,
  DisplayGroup = group_label_en,
  Dataset = "Whole proteome",
  SampleID,
  ConditionLabel = SampleID,
  SampleClass = "ESCC tumor tissue",
  ObservationType = "sample",
  SourceMode = "external_tumor_reference_sample",
  SourceFile = reference_file_rel,
  ReferencePXD = reference_pxd,
  ProteinCount,
  DdrProteinCount,
  DdrFraction,
  DdrFractionPercentage = DdrFraction * 100,
  FrozenKlaProteinCount = length(escc_accessions),
  FrozenKlaDdrProteinCount = length(escc_ddr_accessions),
  FrozenKlaDdrFraction = length(escc_ddr_accessions) / length(escc_accessions) * 100,
  FrozenReferenceProteinCount = length(reference_accessions),
  FrozenReferenceDdrProteinCount = length(reference_ddr_accessions),
  FrozenReferenceDdrFraction = reference_fraction * 100
)]
existing_figure1_values <- fread(file.path(candidate_input_dir, "figure1_sample_boxplot_values.csv"), check.names = FALSE)
updated_figure1_values <- rbindlist(list(existing_figure1_values, new_figure1_values, new_figure1_reference_values), fill = TRUE)
setorder(updated_figure1_values, RowOrder, PXD, SampleGroup, Dataset, SampleID)
fwrite(updated_figure1_values, file.path(candidate_input_dir, "figure1_sample_boxplot_values.csv"), na = "")

new_registry <- data.table(
  PXD = group_pxd,
  SampleGroup = group_name,
  SampleID = sample_ids,
  SampleClass = sample_common$SampleClass[match(sample_ids, sample_common$SampleID)],
  Parser = "MaxQuant La (K)Sites / sample localization",
  SourceFile = source_file_rel,
  Dataset = "Lactylome (Kla)"
)
new_reference_registry <- data.table(
  PXD = group_pxd,
  SampleGroup = group_name,
  SampleID = reference_result$sample_ids,
  SampleClass = "ESCC tumor tissue",
  Parser = "whole-proteome / external PXD065830 individual ESCC-T columns",
  SourceFile = reference_file_rel,
  Dataset = "Whole proteome"
)
for (filename in "sample_boxplot_source_registry.csv") {
  registry <- fread(file.path(candidate_input_dir, filename), check.names = FALSE)
  registry <- rbindlist(list(registry, new_registry), fill = TRUE)
  setorder(registry, Dataset, PXD, SampleGroup, SampleID)
  fwrite(registry, file.path(candidate_input_dir, filename), na = "")
}
registry <- fread(file.path(candidate_input_dir, "figure1_sample_boxplot_source_registry.csv"), check.names = FALSE)
registry <- rbindlist(list(registry, new_registry, new_reference_registry), fill = TRUE)
setorder(registry, Dataset, PXD, SampleGroup, SampleID)
fwrite(registry, file.path(candidate_input_dir, "figure1_sample_boxplot_source_registry.csv"), na = "")

reconciliation <- fread(file.path(candidate_input_dir, "sample_boxplot_reconciliation.csv"), check.names = FALSE)
new_reconciliation <- data.table(
  RowOrder = group_row_order,
  PXD = group_pxd,
  SampleGroup = group_name,
  FrozenKlaProteinCount = length(escc_accessions),
  FrozenKlaDdrProteinCount = length(escc_ddr_accessions),
  FrozenKlaDdrFraction = length(escc_ddr_accessions) / length(escc_accessions) * 100,
  ObservedKlaProteinCount = length(escc_accessions),
  ObservedKlaDdrProteinCount = length(escc_ddr_accessions),
  ObservedSampleCount = length(sample_ids),
  SourceModes = "deposited_sample_table",
  KlaProteinCountMatchesFrozen = TRUE,
  KlaDdrProteinCountMatchesFrozen = TRUE,
  GroupUnionStatus = "PASS"
)
reconciliation <- rbindlist(list(reconciliation, new_reconciliation), fill = TRUE)
setorder(reconciliation, RowOrder, PXD, SampleGroup)
fwrite(reconciliation, file.path(candidate_input_dir, "sample_boxplot_reconciliation.csv"), na = "")

sample_count_record <- fread(file.path(candidate_input_dir, "biological_sample_count_record.csv"), check.names = FALSE)
new_sample_count_record <- data.table(
  RowOrder = group_row_order,
  PXD = group_pxd,
  SampleGroup = group_name,
  KlaSampleCount = length(sample_ids),
  KlaSampleIDs = paste(sample_ids, collapse = ";"),
  ReferencePXD = reference_pxd,
  ReferenceSampleCount = length(reference_result$sample_ids),
  ReferenceSampleIDs = paste(reference_result$sample_ids, collapse = ";"),
  CountingBasis = "six deposited Kla observations plus 94 external ESCC-T whole-proteome observations",
  Notes = paste0("All 94 PXD065830 Dataset1 ESCC-T source columns are plotted individually; the 24 N non-tumor columns are excluded. They are independent ESCC tumor observations, not matched to the six PXD064038 Kla samples.")
)
sample_count_record <- rbindlist(list(sample_count_record, new_sample_count_record), fill = TRUE)
setorder(sample_count_record, RowOrder, PXD, SampleGroup)
fwrite(sample_count_record, file.path(candidate_input_dir, "biological_sample_count_record.csv"), na = "")

sample_design <- fread(file.path(candidate_input_dir, "sample_design_30.csv"), check.names = FALSE, na.strings = c("", "NA"))
new_sample_design <- data.table(
  RowOrder = group_row_order,
  PXD = group_pxd,
  SampleGroup = group_name,
  Category = group_category,
  KlaN = length(sample_ids),
  ReferenceN = length(reference_result$sample_ids),
  KlaSampleDesign = paste(sample_ids, collapse = ";"),
  ReferenceSampleDesign = paste(reference_result$sample_ids, collapse = ";"),
  Aggregation = "six Kla samples versus 94 external ESCC-T source observations",
  MatchClass = "external_disease"
)
sample_design <- rbindlist(list(sample_design, new_sample_design), fill = TRUE)
setorder(sample_design, RowOrder, PXD, SampleGroup)
fwrite(sample_design, file.path(candidate_input_dir, "sample_design_30.csv"), na = "")

pathway_scores <- rbindlist(lapply(s4_sheets, function(sheet) {
  sheet[, c("BaseAccession", pathway_order), with = FALSE]
}), fill = TRUE)
pathway_scores[, BaseAccession := base_accession(BaseAccession)]
pathway_scores <- pathway_scores[, lapply(.SD, function(values) values[[1L]]),
  by = BaseAccession, .SDcols = pathway_order]

new_pathway_values <- rbindlist(lapply(sample_ids, function(sample_id) {
  ddr_ids <- sample_records[SampleID == sample_id & IsDdr == TRUE, BaseAccession]
  scored <- pathway_scores[BaseAccession %in% ddr_ids]
  stop_if(setequal(ddr_ids, scored$BaseAccession), paste0("Missing expanded S4 scores for ", sample_id))
  rbindlist(lapply(pathway_order, function(pathway) {
    states <- scored[[pathway]]
    positive_count <- sum(states == 1, na.rm = TRUE)
    negative_count <- sum(states == -1, na.rm = TRUE)
    data.table(
      RowOrder = group_row_order,
      PXD = group_pxd,
      SampleGroup = group_name,
      Category = group_category,
      DisplayGroup = group_label_en,
      Dataset = "Lactylome (Kla)",
      SampleID = sample_id,
      ConditionLabel = sample_id,
      SampleClass = sample_common$SampleClass[match(sample_id, sample_common$SampleID)],
      ObservationType = "sample",
      SourceMode = "deposited_sample_table",
      SourceFile = source_file_rel,
      Pathway = pathway,
      PositiveProteinCount = positive_count,
      NegativeProteinCount = negative_count,
      AnyPathwayProteinCount = sum(states != 0, na.rm = TRUE),
      KlaDdrProteinCount = length(ddr_ids),
      PositiveFraction = positive_count / length(ddr_ids),
      NegativeFraction = negative_count / length(ddr_ids),
      SignedFraction = (positive_count - negative_count) / length(ddr_ids)
    )
  }))
}), fill = TRUE)
existing_pathway_values <- fread(file.path(candidate_input_dir, "figure1_pathway_summary_sample_boxplot_values.csv"), check.names = FALSE)
updated_pathway_values <- rbindlist(list(existing_pathway_values, new_pathway_values), fill = TRUE)
setorder(updated_pathway_values, RowOrder, PXD, SampleGroup, SampleID)
updated_pathway_values[, PathwayOrder := match(Pathway, pathway_order)]
setorder(updated_pathway_values, RowOrder, PXD, SampleGroup, SampleID, PathwayOrder)
updated_pathway_values[, PathwayOrder := NULL]
fwrite(updated_pathway_values, file.path(candidate_input_dir, "figure1_pathway_summary_sample_boxplot_values.csv"), na = "")

dataset_boxplot_helper <- file.path(
  project_root, "R", "candidate", "build_dataset_level_boxplot_inputs.R"
)
stop_if(file.exists(dataset_boxplot_helper),
  paste0("Missing dataset-level boxplot helper: ", dataset_boxplot_helper))
source(dataset_boxplot_helper, local = TRUE)
dataset_level_figure1_values <- build_dataset_level_figure1_values(expanded_groups)
dataset_level_pathway_values <- build_dataset_level_pathway_summary(
  expanded_groups,
  expanded_kla_membership,
  expanded_reference_membership,
  pathway_scores
)
fwrite(dataset_level_figure1_values,
  file.path(candidate_input_dir, "figure1_dataset_boxplot_values.csv"), na = "")
fwrite(dataset_level_pathway_values,
  file.path(candidate_input_dir, "pathway_summary_dataset_boxplot_values.csv"), na = "")

source_summary <- sample_counts[, .(
  PXD = group_pxd,
  SampleGroup = group_name,
  SampleID,
  KlaProteinCount,
  KlaDdrProteinCount,
  KlaDdrFraction,
  SourceFile = source_file_rel,
  ClinicalMetadataAvailable = file.exists(clinical_file)
)]
fwrite(source_summary, file.path(scope_dir, "escc_sample_summary.csv"), na = "")

audit <- data.table(
  Item = c(
    "selected_group", "excluded_escc_related_dataset_1", "excluded_escc_related_dataset_2",
    "escc_source_pxd", "escc_sample_count", "escc_kla_union_count", "escc_kla_ddr_count",
    "new_all_kla_union_count", "new_kla_ddr_union_count", "previous_catalog_kla_count",
    "strict_reference_status", "reference_pxd", "reference_tumor_sample_count",
    "reference_protein_union_count", "reference_ddr_count", "reference_ddr_fraction",
    "PXD053809_control_status", "PXD010154_status",
    "expanded_group_count", "expanded_tumor_S4_count", "clinical_docx_available"
  ),
  Value = c(
    group_name,
    "PXD048995 (KYSE30 histone-focused; excluded)",
    "PXD063945 (ESCC neoadjuvant; excluded)",
    group_pxd,
    length(sample_ids),
    length(escc_accessions),
    length(escc_ddr_accessions),
    length(setdiff(escc_accessions, current_kla_membership$BaseAccession)),
    length(setdiff(escc_ddr_accessions, current_kla_ddr_venn$BaseAccession)),
    "1246 (legacy catalog note; updated count uses current parser contract)",
    "available; independent ESCC tumor ordinary-proteome union by BaseAccession",
    reference_pxd,
    length(reference_result$sample_ids),
    length(reference_accessions),
    length(reference_ddr_accessions),
    sprintf("%.10f", reference_fraction),
    "not selected; PXD053809 processed XLS is encrypted and not auditable in the current source cache",
    "healthy esophagus background only; not used as matched ordinary-proteome reference",
    nrow(expanded_groups),
    nrow(s4_sheets[["TumorTissues"]]),
    as.character(file.exists(clinical_file))
  )
)
fwrite(audit, file.path(scope_dir, "escc_inclusion_audit.csv"), na = "")

report_lines <- c(
  "# ESCC inclusion audit (2026-09-03)",
  "",
  paste0("- Selected representative group: `", group_pxd, " / ", group_name, "`.",
         " MEC and NEC are retained as six source observations under one cancer-tissue group."),
  "- Excluded from this update: PXD048995 (KYSE30 histone-focused) and PXD063945 (neoadjuvant ESCC).",
  paste0("- Parsed PXD064038: ", length(escc_accessions), " Kla BaseAccessions; ", length(escc_ddr_accessions), " Kla-DDR BaseAccessions."),
  paste0("- Non-lactylated control: iProX/PXD065830 Dataset1, sheet 2.a, using only ", length(reference_result$sample_ids),
         " ESCC tumor T samples (the 24 N non-tumor samples are excluded). The union contains ",
         length(reference_accessions), " BaseAccessions and ", length(reference_ddr_accessions),
         " DDR BaseAccessions (", sprintf("%.3f", reference_fraction * 100), "%)."),
  "- PXD053809 is not selected for this reproducible run because its deposited processed XLS is encrypted; PXD010154 remains a healthy-esophagus background and is not the control.",
  paste0("- Added ", length(setdiff(escc_accessions, current_kla_membership$BaseAccession)), " new all-Kla union proteins and ", length(setdiff(escc_ddr_accessions, current_kla_ddr_venn$BaseAccession)), " new Kla-DDR union proteins."),
  paste0("- The expanded S4 tumor panel adds ", length(new_tumor_ids),
         " cancer-tissue Kla-DDR proteins: existing seven-pathway states are retained where frozen S4 already curated the accession; otherwise broad GO-only accessions receive zero pathway states."),
  "- Dataset-level boxplot inputs are regenerated from the expanded group unions: one point per PXD/sample-group and seven pathways for each of Kla and whole proteome.",
  "",
  "The dated input and candidate output directories are isolated from the frozen publication input and approved result directories."
)
writeLines(report_lines, file.path(scope_dir, "README.md"), useBytes = TRUE)

message(
  "PASS: prepared ", scope_tag, " with ", length(escc_accessions),
  " ESCC Kla proteins, ", length(escc_ddr_accessions), " Kla-DDR proteins, and ",
  nrow(expanded_groups), " publication groups."
)
