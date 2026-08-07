#!/usr/bin/env Rscript
# =============================================================================
# analyze_regulators.R — merged regulator analyses (Task 4 of the R refactor)
#
# Merges three previously separate scripts into run_stage blocks:
#   stage "landscape"      : plot_kla_regulator_landscape.R
#   stage "heatmap"        : analyze_kla_regulator_intensity.R
#   stage "whole_proteome" : analyze_kla_regulator_whole_proteome_intensity.R
#
# Shared helpers come from reanalysis/scripts/lib/ (sourced below); definitions
# duplicated by the lib (base_accession, match_target_accession, accession_feature,
# safe_numeric, relative_path, read_delimited) are removed from the stages and
# replaced by lib calls. match_target_accession(x) -> match_target_accession(x,
# target_accessions) and relative_path(x) -> relative_path(x, project_root) were
# adapted because the lib versions take explicit parameters instead of closures.
# Stage-private functions that differ in signature/behaviour from lib keep their
# private definitions here (see notes per stage).
#
# CLI: Rscript analyze_regulators.R <project_root> [--stage landscape|heatmap|whole_proteome]
#   stage defaults to "all" (run every stage in order: landscape -> heatmap -> whole_proteome).
#   The positional form `<project_root> <stage>` is also accepted.
# =============================================================================
suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(ggplot2)
  library(readxl); library(stringr); library(tidyr)
  library(patchwork)   # landscape stage (plot_kla_regulator_landscape.R) requires patchwork
})
args <- commandArgs(trailingOnly = TRUE)
project_root <- normalizePath(if (length(args)) args[[1]] else ".")
stage <- if (length(args) >= 2) args[[2]] else "all"
# accept the documented "--stage <name>" form as well as the positional form
if (stage == "--stage" && length(args) >= 3) stage <- args[[3]]
lib <- file.path(project_root, "reanalysis", "scripts", "lib")
for (f in c("accession_utils.R", "io_utils.R", "extractors.R"))
  source(file.path(lib, f))
run_stage <- function(name, code) {
  if (stage == "all" || stage == name) { message("[stage] ", name); code }
}
table_dir <- file.path(project_root, "reanalysis", "results", "tables")
figure_dir <- file.path(project_root, "reanalysis", "results", "figures")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

# NOTE on scoping: run_stage evaluates its code argument in the calling (global)
# environment, so stage-local variables are in fact script-level. Each stage block
# re-defines its complete variable set (regulators, sample_catalog, target_accessions,
# quant_parts, audit, plot_data, ...) at the top of the block before use, so the
# same names in different stages never leak into one another; stages communicate
# only through the CSV outputs in table_dir. The landscape stage's private
# extract_maxquant_sites() (signature: path, pxd, sample_group, localization_pattern,
# sheet) is defined inside the landscape block and therefore masks the lib version
# (path, sample_tokens, sheet) from that point on — expected and harmless, because
# no other stage uses the lib version.

# --- Stage: landscape (from plot_kla_regulator_landscape.R) ------------------
run_stage("landscape", {
regulator_path <- file.path(
  project_root, "data", "identifier",
  "乳酸化调控因子_Writer-Eraser-Reader.xlsx"
)
mapping_path <- file.path(
  project_root, "reanalysis", "config",
  "lactylation_regulator_uniprot_mapping.csv"
)
pairing_path <- file.path(
  table_dir,
  "lactylome_and_reference_proteome_pairing_zh.csv"
)
primary_path <- file.path(
  project_root, "reanalysis", "intermediate", "kla_by_dataset",
  "all_primary_sample_level_kla_sites.csv"
)

required_files <- c(regulator_path, mapping_path, pairing_path, primary_path)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) {
  stop("Missing required input files: ", paste(missing_files, collapse = ", "))
}

regulators <- read_excel(regulator_path) |>
  transmute(
    Role = trimws(as.character(分组)),
    GeneSymbol = toupper(trimws(as.character(基因))),
    Reference = as.character(参考文献),
    RoleEntryOrder = row_number()
  ) |>
  filter(
    Role %in% c("Writer", "Eraser", "Reader"),
    !is.na(GeneSymbol),
    nzchar(GeneSymbol)
  )
target_genes <- unique(regulators$GeneSymbol)

accession_map <- read.csv(
  mapping_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
) |>
  mutate(
    GeneSymbol = toupper(GeneSymbol),
    BaseAccession = sub("-[0-9]+$", "", BaseAccession)
  )
accession_to_gene <- setNames(accession_map$GeneSymbol, accession_map$BaseAccession)

pairing <- read.csv(
  pairing_path,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  na.strings = c("", "NA")
)
sample_catalog <- pairing |>
  filter(当前已具备成对计数条件) |>
  transmute(
    PXD = 乳酸化PXD,
    SampleGroup = 样本组,
    BiologicalMaterial = 材料类型,
    LactylomeEvidenceFile = 乳酸化证据文件,
    LactylomeProteinCount = 乳酸化蛋白数,
    RowOrder = row_number(),
    SampleGroupID = paste(PXD, SampleGroup, sep = "__")
  )
if (nrow(sample_catalog) != 40) {
  stop("Expected 40 pair-ready sample groups, found ", nrow(sample_catalog))
}

display_names <- c(
  "human hippocampus" = "人海马组织",
  "pathological rotator cuff tendon" = "病理性肩袖肌腱",
  "normal human lung" = "正常人肺组织",
  "normal liver" = "正常肝组织",
  "nonmetastatic HCC" = "非转移性HCC",
  "lung-metastatic HCC" = "肺转移HCC",
  "hypertrophic scar" = "增生性瘢痕",
  "adjacent skin" = "邻近皮肤",
  "HepG2 WT and SIRT1 or SIRT3 KO" = "HepG2 WT/SIRT1-KO/SIRT3-KO",
  "human fibroblasts mock and HCMV or HSV-1" = "人成纤维细胞 mock/HCMV/HSV-1",
  "human fibroblasts mock and HCMV" = "人成纤维细胞 mock/HCMV",
  "pretreated HK-2" = "预处理HK-2",
  "bladder cancer cells treated with EPI" = "EPI处理膀胱癌细胞",
  "normal pregnancy placenta" = "正常妊娠胎盘",
  "severe preeclampsia placenta" = "重度子痫前期胎盘",
  "MEC and NEC ESCC groups" = "ESCC MEC/NEC组",
  "human sperm" = "人精子",
  "BPH" = "良性前列腺增生",
  "prostate cancer" = "前列腺癌",
  "HCT116 control and Roseburia co-culture" = "HCT116对照/Roseburia共培养",
  "glioblastoma stem cells" = "胶质母细胞瘤干细胞",
  "neural stem cells" = "神经干细胞",
  "HUVEC control and Pg infection" = "HUVEC对照/Pg感染",
  "AC16 control and hypoxia" = "AC16对照/低氧",
  "HCC" = "肝细胞癌",
  "adjacent liver" = "邻近肝组织",
  "RKO WT and GSK3B KO" = "RKO WT/GSK3B-KO",
  "HK-2 control and mannitol" = "HK-2对照/甘露醇"
)
sample_catalog <- sample_catalog |>
  mutate(
    DisplayName = ifelse(
      SampleGroup %in% names(display_names),
      unname(display_names[SampleGroup]),
      SampleGroup
    ),
    RowLabel = paste0(DisplayName, " · ", PXD)
  )

split_target_genes <- function(values) {
  values <- as.character(values)
  rows <- lapply(seq_along(values), function(index) {
    genes <- unlist(strsplit(values[[index]], "[;,/ ]+"))
    genes <- toupper(trimws(genes))
    genes <- unique(genes[genes %in% target_genes])
    if (!length(genes)) {
      return(NULL)
    }
    data.frame(SourceRow = index, GeneSymbol = genes)
  })
  bind_rows(rows)
}


empty_evidence <- function() {
  data.frame(
    PXD = character(),
    SampleGroup = character(),
    GeneSymbol = character(),
    KlaSite = character(),
    EvidenceResolution = character(),
    SourceFile = character(),
    stringsAsFactors = FALSE
  )
}

extract_maxquant_sites <- function(
  path,
  pxd,
  sample_group,
  localization_pattern = NULL,
  sheet = NULL
) {
  if (!file.exists(path)) {
    return(empty_evidence())
  }
  data <- if (grepl("\\.xlsx$", path, ignore.case = TRUE)) {
    read_excel(path, sheet = sheet)
  } else {
    read.delim(
      path,
      check.names = FALSE,
      stringsAsFactors = FALSE,
      quote = "",
      comment.char = ""
    )
  }
  if (!"Gene names" %in% names(data)) {
    return(empty_evidence())
  }
  keep <- rep(TRUE, nrow(data))
  if ("Reverse" %in% names(data)) {
    keep <- keep & (is.na(data$Reverse) | data$Reverse != "+")
  }
  if ("Potential contaminant" %in% names(data)) {
    keep <- keep &
      (is.na(data$`Potential contaminant`) | data$`Potential contaminant` != "+")
  }
  if ("id" %in% names(data)) {
    keep <- keep & !is.na(data$id)
  }
  localization_columns <- if (is.null(localization_pattern)) {
    if ("Localization prob" %in% names(data)) "Localization prob" else character()
  } else {
    grep(
      paste0("^Localization prob .*", localization_pattern),
      names(data),
      value = TRUE,
      ignore.case = TRUE
    )
  }
  if (length(localization_columns)) {
    localized <- rowSums(
      sapply(
        data[localization_columns],
        function(values) suppressWarnings(as.numeric(values)) > 0
      ),
      na.rm = TRUE
    ) > 0
    keep <- keep & localized
  }
  data <- data[keep, , drop = FALSE]
  if (!nrow(data)) {
    return(empty_evidence())
  }
  gene_rows <- split_target_genes(data$`Gene names`)
  if (!nrow(gene_rows)) {
    return(empty_evidence())
  }
  position_column <- intersect(
    c("Position", "Positions within proteins", "Positions"),
    names(data)
  )
  positions <- if (length(position_column)) {
    as.character(data[[position_column[[1]]]])
  } else {
    rep(NA_character_, nrow(data))
  }
  gene_rows |>
    transmute(
      PXD = pxd,
      SampleGroup = sample_group,
      GeneSymbol,
      KlaSite = ifelse(
        is.na(positions[SourceRow]) | !nzchar(positions[SourceRow]),
        NA_character_,
        paste0("K", sub(";.*$", "", positions[SourceRow]))
      ),
      EvidenceResolution = "site",
      SourceFile = sub(paste0("^", project_root, "/"), "", path)
    )
}

primary <- read.csv(
  primary_path,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  na.strings = c("", "NA")
)
core_alias <- data.frame(
  PXD = c(
    "PXD014870", "PXD028488", "PXD028488", "PXD028488",
    "PXD050470", "PXD053474", "PXD060185", "PXD060185",
    "PXD060185", "PXD060185", "PXD078013", "PXD078736"
  ),
  CellOrTissueType = c(
    "MCF7", "HEK293T", "HCT116", "T-ALL",
    "Human hippocampus", "HCT116", "MCF10A", "MCF7",
    "MDA-MB-468", "T-47D", "RKO", "HK-2"
  ),
  SampleGroup = c(
    "MCF7", "HEK293T", "HCT116", "TALL-104",
    "human hippocampus", "HCT116", "MCF10A", "MCF7",
    "MDA-MB-468", "T-47D", "RKO WT and GSK3B KO",
    "HK-2 control and mannitol"
  ),
  stringsAsFactors = FALSE
)
core_evidence <- primary |>
  filter(
    PrimaryIncluded %in% c(TRUE, "TRUE", "True", 1, "1"),
    !is.na(GeneSymbol),
    toupper(GeneSymbol) %in% target_genes
  ) |>
  inner_join(core_alias, by = c("PXD", "CellOrTissueType")) |>
  transmute(
    PXD,
    SampleGroup,
    GeneSymbol = toupper(GeneSymbol),
    KlaSite,
    EvidenceResolution = "site",
    SourceFile
  )

evidence_parts <- list(core_evidence)
audit <- sample_catalog |>
  transmute(
    PXD,
    SampleGroup,
    SampleGroupID,
    GeneLevelAuditStatus = "逐蛋白可审计",
    Parser = "pending",
    AuditReason = "",
    SourceFile = LactylomeEvidenceFile
  )

add_evidence <- function(data) {
  evidence_parts[[length(evidence_parts) + 1]] <<- data
}
set_audit <- function(pxd, group, status, parser, reason = "") {
  rows <- audit$PXD == pxd & audit$SampleGroup == group
  audit$GeneLevelAuditStatus[rows] <<- status
  audit$Parser[rows] <<- parser
  audit$AuditReason[rows] <<- reason
}

core_rows <- unique(core_alias[, c("PXD", "SampleGroup")])
for (i in seq_len(nrow(core_rows))) {
  set_audit(
    core_rows$PXD[[i]],
    core_rows$SampleGroup[[i]],
    "逐蛋白可审计",
    "core_sample_level_long_table"
  )
}

maxquant_jobs <- list(
  list("PXD033146", "pathological rotator cuff tendon",
       "data/PXD033146/search_results/extracted_pairing/search_result-HA119TPLa/La (K)Sites.txt", NULL, NULL),
  list("PXD036307", "normal human lung",
       "data/PXD036307/search_results/extracted/txt/La (K)Sites.txt", NULL, NULL),
  list("PXD050147", "HepG2 WT and SIRT1 or SIRT3 KO",
       "data/PXD050147/search_results/Lactyl_K_Sites.txt", NULL, NULL),
  list("PXD058534", "pretreated HK-2",
       "data/PXD058534/search_results/extracted_pairing/txt/La (K)Sites.txt", NULL, NULL),
  list("PXD062720", "bladder cancer cells treated with EPI",
       "data/PXD062720/search_results/extracted_pairing/txt/La (K)Sites.txt", NULL, NULL),
  list("PXD063047", "normal pregnancy placenta",
       "data/PXD063047/search_results/extracted/combined/txt/La (K)Sites.txt", "Con_", NULL),
  list("PXD063047", "severe preeclampsia placenta",
       "data/PXD063047/search_results/extracted/combined/txt/La (K)Sites.txt", "PE_", NULL),
  list("PXD063266", "PC-3M",
       "data/PXD063266/search_results/LactylSites.xlsx", NULL, "Lactyl (K)Sites"),
  list("PXD064038", "MEC and NEC ESCC groups",
       "data/PXD064038/search_results/extracted_pairing/txt/txt/La (K)Sites.txt", NULL, NULL)
)
for (job in maxquant_jobs) {
  path <- file.path(project_root, job[[3]])
  add_evidence(extract_maxquant_sites(
    path, job[[1]], job[[2]], job[[4]], job[[5]]
  ))
  set_audit(job[[1]], job[[2]], "逐蛋白可审计", "maxquant_site_table")
}

extract_pxd046800 <- function() {
  path <- file.path(
    project_root,
    "data/PXD046800/search_results/HFX2_LFQ_QB001_Lacty_PeptideGroups.txt"
  )
  data <- read.delim(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    quote = "\"",
    comment.char = ""
  )
  data <- data[
    grepl("Lacty|Lactyl|La \\(K\\)", data$Modifications, ignore.case = TRUE),
    ,
    drop = FALSE
  ]
  accessions <- base_accession(data$`Master Protein Accessions`)
  genes <- unname(accession_to_gene[accessions])
  modifications <- as.character(data$Modifications)
  out <- list()
  group_specs <- list(
    "hypertrophic scar" = "HSP",
    "adjacent skin" = "NSP"
  )
  for (group in names(group_specs)) {
    columns <- grep(
      paste0("Found in Sample:.*", group_specs[[group]]),
      names(data),
      value = TRUE
    )
    keep <- genes %in% target_genes &
      rowSums(sapply(data[columns], function(values) {
        !is.na(values) & nzchar(values) & values != "Not Found"
      })) > 0
    indexes <- which(keep)
    if (!length(indexes)) next
    out[[group]] <- bind_rows(lapply(indexes, function(index) {
      sites <- unique(unlist(str_extract_all(modifications[[index]], "K[0-9]+")))
      if (!length(sites)) sites <- NA_character_
      data.frame(
        PXD = "PXD046800",
        SampleGroup = group,
        GeneSymbol = genes[[index]],
        KlaSite = sites,
        EvidenceResolution = ifelse(is.na(sites), "modified_peptide", "site"),
        SourceFile = "data/PXD046800/search_results/HFX2_LFQ_QB001_Lacty_PeptideGroups.txt"
      )
    }))
  }
  bind_rows(out)
}
add_evidence(extract_pxd046800())
set_audit("PXD046800", "hypertrophic scar", "逐蛋白可审计", "proteome_discoverer_lactylated_peptide_table")
set_audit("PXD046800", "adjacent skin", "逐蛋白可审计", "proteome_discoverer_lactylated_peptide_table")

extract_a549 <- function() {
  path <- file.path(
    project_root,
    "data/PXD054919/supplementary/41419_2025_8113_MOESM2_ESM.xlsx"
  )
  data <- read_excel(path, skip = 1)
  data |>
    transmute(
      PXD = "PXD054919",
      SampleGroup = "A549",
      GeneSymbol = toupper(`Gene name`),
      KlaSite = paste0("K", Position),
      EvidenceResolution = "site",
      SourceFile = "data/PXD054919/supplementary/41419_2025_8113_MOESM2_ESM.xlsx"
    ) |>
    filter(GeneSymbol %in% target_genes)
}
add_evidence(extract_a549())
set_audit("PXD054919", "A549", "逐蛋白可审计", "author_supplementary_site_table")

extract_spectronaut <- function(path, pxd, sample_group, require_site_confidence) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("data.table is required for Spectronaut report parsing")
  }
  columns <- c(
    "Protein.Group", "Genes", "Modified.Sequence",
    "Q.Value", "PG.Q.Value"
  )
  header <- names(data.table::fread(path, nrows = 0, data.table = FALSE))
  columns <- intersect(columns, header)
  if (require_site_confidence && "PTM.Site.Confidence" %in% header) {
    columns <- c(columns, "PTM.Site.Confidence")
  }
  data <- data.table::fread(
    path,
    select = columns,
    showProgress = FALSE,
    data.table = FALSE
  )
  keep <- grepl("K\\(UniMod:378\\)", data$Modified.Sequence) &
    suppressWarnings(as.numeric(data$Q.Value)) <= 0.01 &
    suppressWarnings(as.numeric(data$PG.Q.Value)) <= 0.01
  if ("PTM.Site.Confidence" %in% names(data)) {
    keep <- keep &
      suppressWarnings(as.numeric(data$PTM.Site.Confidence)) > 0
  }
  data <- data[keep, , drop = FALSE]
  gene_rows <- split_target_genes(data$Genes)
  if (!nrow(gene_rows)) return(empty_evidence())
  gene_rows |>
    transmute(
      PXD = pxd,
      SampleGroup = sample_group,
      GeneSymbol,
      KlaSite = NA_character_,
      EvidenceResolution = "modified_peptide",
      SourceFile = sub(paste0("^", project_root, "/"), "", path)
    )
}
add_evidence(extract_spectronaut(
  file.path(project_root, "data/PXD055230/search_results/LaIP_HSV1_DIA.tsv"),
  "PXD055230", "human fibroblasts mock and HCMV or HSV-1", TRUE
))
add_evidence(extract_spectronaut(
  file.path(project_root, "data/PXD057709/search_results/LaIP_report.tsv"),
  "PXD057709", "human fibroblasts mock and HCMV", FALSE
))
set_audit(
  "PXD055230", "human fibroblasts mock and HCMV or HSV-1",
  "逐蛋白可审计", "spectronaut_modified_peptide_report"
)
set_audit(
  "PXD057709", "human fibroblasts mock and HCMV",
  "逐蛋白可审计", "spectronaut_modified_peptide_report"
)

extract_pxd064912 <- function() {
  path <- file.path(
    project_root,
    "data/PXD064912/supplementary/europepmc/mmc1.xlsx"
  )
  data <- read_excel(path, skip = 1)
  prob_columns <- grep("^PTM.SiteProbability", names(data), value = TRUE)
  keep <- tolower(data$PTM.ModificationTitle) == "lactylation" &
    data$PTM.SiteAA == "K" &
    rowSums(sapply(data[prob_columns], function(values) {
      suppressWarnings(as.numeric(values)) > 0
    }), na.rm = TRUE) > 0
  data <- data[keep, , drop = FALSE]
  gene_rows <- split_target_genes(data$PG.Genes)
  if (!nrow(gene_rows)) return(empty_evidence())
  gene_rows |>
    transmute(
      PXD = "PXD064912",
      SampleGroup = "human sperm",
      GeneSymbol,
      KlaSite = paste0("K", data$PTM.SiteLocation[SourceRow]),
      EvidenceResolution = "site",
      SourceFile = "data/PXD064912/supplementary/europepmc/mmc1.xlsx"
    )
}
add_evidence(extract_pxd064912())
set_audit("PXD064912", "human sperm", "逐蛋白可审计", "author_supplementary_site_table")

extract_pxd066054 <- function(group) {
  path <- file.path(
    project_root,
    "data/PXD066054/search_results/extracted/PLa/XB08700B1DPLa_-PTMSiteReport.tsv"
  )
  data <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  prefix <- if (group == "BPH") "^NAT" else "^PCa"
  data <- data[
    data$PTM.ModificationTitle == "L-Lac(K)" &
      data$PTM.SiteAA == "K" &
      grepl(prefix, data$R.Condition) &
      suppressWarnings(as.numeric(data$PTM.SiteProbability)) > 0,
  ]
  gene_rows <- split_target_genes(data$PG.Genes)
  if (!nrow(gene_rows)) return(empty_evidence())
  gene_rows |>
    transmute(
      PXD = "PXD066054",
      SampleGroup = group,
      GeneSymbol,
      KlaSite = paste0("K", data$PTM.SiteLocation[SourceRow]),
      EvidenceResolution = "site",
      SourceFile = "data/PXD066054/search_results/extracted/PLa/XB08700B1DPLa_-PTMSiteReport.tsv"
    )
}
for (group in c("BPH", "prostate cancer")) {
  add_evidence(extract_pxd066054(group))
  set_audit("PXD066054", group, "逐蛋白可审计", "spectronaut_ptm_site_report")
}

extract_pxd066351 <- function() {
  path <- file.path(
    project_root,
    "data/PXD066351/search_results/XB01472B1DPLa-MSstats_Input.csv"
  )
  data <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  data <- data[
    data$PTM.ModificationTitle == "Lac (K)" &
      data$PTM.SiteAA == "K" &
      suppressWarnings(as.numeric(data$PTM.SiteProbability)) > 0 &
      suppressWarnings(as.numeric(data$PTM.Quantity)) > 0,
    ,
    drop = FALSE
  ]
  accessions <- unique(base_accession(data$PTM.ProteinId))
  genes <- unname(accession_to_gene[accessions])
  genes <- unique(genes[!is.na(genes) & genes %in% target_genes])
  data.frame(
    PXD = "PXD066351",
    SampleGroup = "HCT116 control and Roseburia co-culture",
    GeneSymbol = genes,
    KlaSite = NA_character_,
    EvidenceResolution = "site",
    SourceFile = "data/PXD066351/search_results/XB01472B1DPLa-MSstats_Input.csv"
  )
}
add_evidence(extract_pxd066351())
set_audit(
  "PXD066351", "HCT116 control and Roseburia co-culture",
  "逐蛋白可审计", "msstats_lactylation_ptm_table"
)

extract_pxd070007 <- function(group) {
  path <- file.path(
    project_root,
    "data/PXD070007/search_results/SA206LPLaB1_Annotation.xlsx"
  )
  data <- read_excel(path, sheet = "Annotation_Combine")
  sample_columns <- if (group == "glioblastoma stem cells") {
    c("G2907", "G3028", "G3264", "GSC23", "MES28", "RKI")
  } else {
    c("ENSA", "HMP1")
  }
  detected <- rowSums(sapply(data[sample_columns], function(values) {
    suppressWarnings(as.numeric(values)) > 0
  }), na.rm = TRUE) > 0
  keep <- detected &
    suppressWarnings(as.numeric(data$`Localization probability`)) > 0
  data <- data[keep, , drop = FALSE]
  data.frame(
    PXD = "PXD070007",
    SampleGroup = group,
    GeneSymbol = toupper(data$`Gene name`),
    KlaSite = paste0("K", data$Position),
    EvidenceResolution = "site",
    SourceFile = "data/PXD070007/search_results/SA206LPLaB1_Annotation.xlsx"
  ) |>
    filter(GeneSymbol %in% target_genes)
}
for (group in c("glioblastoma stem cells", "neural stem cells")) {
  add_evidence(extract_pxd070007(group))
  set_audit("PXD070007", group, "逐蛋白可审计", "author_annotation_site_table")
}

extract_pxd028737 <- function() {
  path <- file.path(
    project_root,
    "data/PXD028737/search_results/extracted_pairing/combined/combined/txt/La(K)Sites.txt"
  )
  extract_maxquant_sites(
    path,
    "PXD028737",
    "HMC3"
  )
}
add_evidence(extract_pxd028737())
set_audit(
  "PXD028737", "HMC3",
  "逐蛋白可审计", "maxquant_site_table",
  ""
)

extract_pxd073311 <- function() {
  path <- file.path(
    project_root,
    paste0(
      "data/PXD073311/search_results/extracted_pairing/",
      "IPX0015307003_Database_search_result/Database_search_result/",
      "report.pr_matrix.tsv"
    )
  )
  data <- read.delim(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = ""
  )
  data <- data[
    grepl("K\\(UniMod:2114\\)", data$Modified.Sequence),
    ,
    drop = FALSE
  ]
  gene_rows <- split_target_genes(data$Genes)
  if (!nrow(gene_rows)) return(empty_evidence())
  gene_rows |>
    transmute(
      PXD = "PXD073311",
      SampleGroup = "HUVEC control and Pg infection",
      GeneSymbol,
      KlaSite = NA_character_,
      EvidenceResolution = "modified_precursor",
      SourceFile = sub(paste0("^", project_root, "/"), "", path)
    )
}
add_evidence(extract_pxd073311())
set_audit(
  "PXD073311", "HUVEC control and Pg infection",
  "逐蛋白可审计", "spectronaut_lactyl_precursor_matrix",
  ""
)

extract_pxd075014 <- function() {
  path <- file.path(
    project_root,
    "data/PXD075014/supplementary/Table2.XLSX"
  )
  data <- read_excel(path, sheet = "Kla peptides")
  data <- data[
    grepl("lactylation", data$Modifications, ignore.case = TRUE),
    ,
    drop = FALSE
  ]
  genes <- toupper(str_match(
    data$`Master Protein Descriptions`,
    "GN=([^ ]+)"
  )[, 2])
  data.frame(
    PXD = "PXD075014",
    SampleGroup = "AC16 control and hypoxia",
    GeneSymbol = genes,
    KlaSite = NA_character_,
    EvidenceResolution = "modified_peptide",
    SourceFile = "data/PXD075014/supplementary/Table2.XLSX"
  ) |>
    filter(GeneSymbol %in% target_genes)
}
add_evidence(extract_pxd075014())
set_audit(
  "PXD075014", "AC16 control and hypoxia",
  "逐蛋白可审计", "author_supplementary_modified_peptide_table",
  ""
)

extract_pxd075377 <- function(group) {
  path <- file.path(
    project_root,
    "data/PXD075377/search_results/extracted/2-Basic_analysis/MS_identified_information.txt"
  )
  data <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  intensity_column <- if (group == "HCC") "Intensity HCC" else "Intensity Control"
  keep <- suppressWarnings(as.numeric(data[[intensity_column]])) > 0 &
    suppressWarnings(as.numeric(data$`Localization probability`)) > 0
  data <- data[keep, , drop = FALSE]
  data.frame(
    PXD = "PXD075377",
    SampleGroup = group,
    GeneSymbol = toupper(data$`Gene name`),
    KlaSite = paste0("K", data$Position),
    EvidenceResolution = "site",
    SourceFile = "data/PXD075377/search_results/extracted/2-Basic_analysis/MS_identified_information.txt"
  ) |>
    filter(GeneSymbol %in% target_genes)
}
for (group in c("HCC", "adjacent liver")) {
  add_evidence(extract_pxd075377(group))
  set_audit("PXD075377", group, "逐蛋白可审计", "author_site_intensity_table")
}

for (group in sample_catalog$SampleGroup[sample_catalog$PXD == "PXD037371"]) {
  set_audit(
    "PXD037371", group, "组别无法拆分", "dataset_union_only",
    "本地位点表有蛋白身份，但当前文件未提供三类临床组与通道的可靠映射。"
  )
}
if (any(audit$Parser == "pending")) {
  stop(
    "Unassigned parser for: ",
    paste(
      audit$SampleGroupID[audit$Parser == "pending"],
      collapse = ", "
    )
  )
}

all_evidence <- bind_rows(evidence_parts) |>
  filter(GeneSymbol %in% target_genes) |>
  distinct()

evidence_summary <- all_evidence |>
  group_by(PXD, SampleGroup, GeneSymbol) |>
  summarise(
    Detected = TRUE,
    ExactSiteCount = if (any(EvidenceResolution == "site" & !is.na(KlaSite))) {
      n_distinct(KlaSite[EvidenceResolution == "site" & !is.na(KlaSite)])
    } else {
      NA_integer_
    },
    EvidenceResolution = paste(
      sort(unique(EvidenceResolution)),
      collapse = ";"
    ),
    SourceFile = paste(sort(unique(SourceFile)), collapse = ";"),
    .groups = "drop"
  )

plot_data <- crossing(
  regulators |>
    select(Role, GeneSymbol, Reference, RoleEntryOrder),
  sample_catalog |>
    select(
      PXD, SampleGroup, SampleGroupID, RowOrder, RowLabel,
      BiologicalMaterial, LactylomeProteinCount
    )
) |>
  left_join(
    audit |>
      select(
        PXD, SampleGroup, GeneLevelAuditStatus,
        Parser, AuditReason
      ),
    by = c("PXD", "SampleGroup")
  ) |>
  left_join(
    evidence_summary,
    by = c("PXD", "SampleGroup", "GeneSymbol")
  ) |>
  mutate(
    Detected = replace_na(Detected, FALSE),
    CellState = case_when(
      GeneLevelAuditStatus != "逐蛋白可审计" ~ "无法逐蛋白判断",
      !Detected ~ "未检出",
      is.na(ExactSiteCount) ~ "检出，位点数不可得",
      ExactSiteCount == 1 ~ "1个位点",
      ExactSiteCount <= 3 ~ "2-3个位点",
      TRUE ~ "≥4个位点"
    ),
    TileLabel = case_when(
      GeneLevelAuditStatus != "逐蛋白可审计" ~ "?",
      !Detected ~ "",
      is.na(ExactSiteCount) ~ "P",
      TRUE ~ as.character(ExactSiteCount)
    )
  ) |>
  arrange(RoleEntryOrder, RowOrder)

write.csv(
  plot_data |>
    select(
      Role, GeneSymbol, PXD, SampleGroup, SampleGroupID,
      RowLabel, BiologicalMaterial, LactylomeProteinCount,
      GeneLevelAuditStatus, Parser, AuditReason,
      Detected, ExactSiteCount, EvidenceResolution,
      CellState, SourceFile, Reference
    ),
  file.path(table_dir, "kla_regulator_cell_type_long.csv"),
  row.names = FALSE
)
write.csv(
  audit,
  file.path(table_dir, "kla_regulator_dataset_audit.csv"),
  row.names = FALSE
)

matrix_table <- plot_data |>
  select(Role, RoleEntryOrder, GeneSymbol, RowLabel, CellState) |>
  pivot_wider(names_from = RowLabel, values_from = CellState) |>
  arrange(RoleEntryOrder) |>
  select(-RoleEntryOrder)
write.csv(
  matrix_table,
  file.path(table_dir, "kla_regulator_cell_type_site_matrix.csv"),
  row.names = FALSE
)

gene_summary <- plot_data |>
  filter(!duplicated(paste(PXD, SampleGroup, GeneSymbol))) |>
  group_by(GeneSymbol) |>
  summarise(
    Roles = paste(unique(regulators$Role[regulators$GeneSymbol == first(GeneSymbol)]), collapse = "/"),
    AuditableGroupCount = sum(GeneLevelAuditStatus == "逐蛋白可审计"),
    DetectedGroupCount = sum(Detected & GeneLevelAuditStatus == "逐蛋白可审计"),
    GroupsWithExactSites = sum(!is.na(ExactSiteCount)),
    MaximumExactSiteCount = if (all(is.na(ExactSiteCount))) 0 else max(ExactSiteCount, na.rm = TRUE),
    DetectedGroups = paste(RowLabel[Detected], collapse = ";"),
    .groups = "drop"
  ) |>
  arrange(desc(DetectedGroupCount), GeneSymbol)
write.csv(
  gene_summary,
  file.path(table_dir, "kla_regulator_detection_summary.csv"),
  row.names = FALSE
)

role_colors <- c(
  Writer = "#1976B9",
  Eraser = "#2A9D4B",
  Reader = "#E85D0F"
)
state_colors <- c(
  "无法逐蛋白判断" = "#9CA3AF",
  "未检出" = "#F3F4F6",
  "检出，位点数不可得" = "#E9B949",
  "1个位点" = "#CFE8F3",
  "2-3个位点" = "#5FA8C6",
  "≥4个位点" = "#155B7A"
)
state_levels <- names(state_colors)
font_family <- "Arial Unicode MS"
auditable_row_labels <- sample_catalog |>
  inner_join(
    audit |>
      filter(GeneLevelAuditStatus == "逐蛋白可审计") |>
      select(PXD, SampleGroup),
    by = c("PXD", "SampleGroup")
  ) |>
  pull(RowLabel)

make_role_panel <- function(role, show_y_axis) {
  role_data <- plot_data |>
    filter(
      Role == role,
      GeneLevelAuditStatus == "逐蛋白可审计"
    ) |>
    mutate(
      GeneEntry = factor(
        paste0(GeneSymbol, "__", RoleEntryOrder),
        levels = paste0(
          regulators$GeneSymbol[regulators$Role == role],
          "__",
          regulators$RoleEntryOrder[regulators$Role == role]
        )
      ),
      RowLabel = factor(
        RowLabel,
        levels = rev(auditable_row_labels)
      ),
      CellState = factor(CellState, levels = state_levels),
      TextColor = case_when(
        CellState %in% c("≥4个位点", "2-3个位点") ~ "white",
        CellState == "无法逐蛋白判断" ~ "white",
        TRUE ~ "#17324D"
      )
    )
  gene_labels <- setNames(
    role_data$GeneSymbol[!duplicated(role_data$GeneEntry)],
    as.character(role_data$GeneEntry[!duplicated(role_data$GeneEntry)])
  )

  ggplot(role_data, aes(x = GeneEntry, y = RowLabel)) +
    geom_tile(
      aes(fill = CellState),
      color = "#D8DEE5",
      linewidth = 0.18,
      width = 0.97,
      height = 0.97
    ) +
    geom_text(
      aes(label = TileLabel, color = TextColor),
      size = 1.9,
      family = font_family,
      show.legend = FALSE
    ) +
    scale_color_identity() +
    scale_fill_manual(
      values = state_colors,
      limits = state_levels,
      drop = FALSE,
      name = "Kla证据状态",
      labels = c(
        "无法逐蛋白判断",
        "未检出",
        "检出（位点数不可得）",
        "1个位点",
        "2-3个位点",
        "≥4个位点"
      )
    ) +
    scale_x_discrete(labels = gene_labels, expand = expansion(mult = 0)) +
    scale_y_discrete(drop = FALSE, expand = expansion(mult = 0)) +
    labs(title = role, x = NULL, y = NULL) +
    theme_minimal(base_size = 8, base_family = font_family) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(
        angle = 90, hjust = 1, vjust = 0.5,
        color = "#30343B", size = 7
      ),
      axis.text.y = if (show_y_axis) {
        element_text(color = "#30343B", size = 6.7)
      } else {
        element_blank()
      },
      axis.ticks = element_blank(),
      plot.title = element_text(
        color = unname(role_colors[[role]]),
        face = "bold", hjust = 0.5, size = 13,
        margin = margin(b = 5)
      ),
      plot.margin = margin(2, 2, 2, 2),
      legend.position = "right",
      legend.title = element_text(face = "bold"),
      legend.text = element_text(size = 8),
      legend.key.height = grid::unit(0.46, "cm")
    )
}

combined_plot <- (
  make_role_panel("Writer", TRUE) +
    make_role_panel("Eraser", FALSE) +
    make_role_panel("Reader", FALSE) +
    plot_layout(
      widths = c(
        sum(regulators$Role == "Writer"),
        sum(regulators$Role == "Eraser"),
        sum(regulators$Role == "Reader")
      ),
      guides = "collect"
    )
) +
  plot_annotation(
    title = "乳酸化调控因子在37个可审计细胞系与组织样本组中的Kla证据图谱",
    subtitle = paste0(
      "数字为可定位的唯一Kla位点数；P表示蛋白/修饰肽层面检出但无法得到精确位点数；",
      "无法可靠拆分的样本组不进入本图。颜色不代表表达量或Log2 FC。"
    ),
    caption = paste0(
      "数据范围：40个候选样本组中37个具有逐蛋白证据；PXD037371的3组因TMT通道映射不明而排除。",
      "主分析定位阈值为>0；多重角色蛋白在相应Writer/Eraser/Reader分区重复展示。"
    ),
    theme = theme(
      text = element_text(family = font_family, color = "#20252B"),
      plot.title = element_text(face = "bold", size = 16, hjust = 0),
      plot.subtitle = element_text(size = 9.5, hjust = 0, margin = margin(b = 8)),
      plot.caption = element_text(size = 8, hjust = 0, color = "#5A626B")
    )
  ) &
  theme(legend.position = "right")

png_path <- file.path(
  figure_dir,
  "kla_regulator_cell_type_detection_landscape.png"
)
pdf_path <- file.path(
  figure_dir,
  "kla_regulator_cell_type_detection_landscape.pdf"
)
ggsave(
  png_path,
  combined_plot,
  width = 17,
  height = 16.5,
  units = "in",
  dpi = 300,
  bg = "white",
  limitsize = FALSE
)
ggsave(
  pdf_path,
  combined_plot,
  width = 17,
  height = 16.5,
  units = "in",
  device = cairo_pdf,
  bg = "white",
  limitsize = FALSE
)

message("Built expanded Kla regulator landscape:")
message("- ", png_path)
message("- ", pdf_path)
message("- ", nrow(sample_catalog), " sample groups")
message("- ", sum(audit$GeneLevelAuditStatus == "逐蛋白可审计"), " auditable groups")
})

# --- Stage: heatmap (from analyze_kla_regulator_intensity.R) ------------------
run_stage("heatmap", {
regulator_path <- file.path(
  project_root, "data", "identifier",
  "乳酸化调控因子_Writer-Eraser-Reader.xlsx"
)
detection_path <- file.path(table_dir, "kla_regulator_cell_type_long.csv")
mapping_path <- file.path(
  project_root, "reanalysis", "config",
  "lactylation_regulator_uniprot_mapping.csv"
)

required_files <- c(regulator_path, detection_path, mapping_path)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) {
  stop("Missing required files: ", paste(missing_files, collapse = ", "))
}

regulators <- read_excel(regulator_path) |>
  transmute(
    Role = trimws(as.character(分组)),
    GeneSymbol = toupper(trimws(as.character(基因))),
    RoleEntryOrder = row_number()
  ) |>
  filter(
    Role %in% c("Writer", "Eraser", "Reader"),
    !is.na(GeneSymbol),
    nzchar(GeneSymbol)
  )
target_genes <- unique(regulators$GeneSymbol)
plot_font <- "Arial Unicode MS"

accession_map <- read.csv(
  mapping_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
) |>
  transmute(
    BaseAccession = sub("-[0-9]+$", "", BaseAccession),
    GeneSymbol = toupper(GeneSymbol)
  ) |>
  filter(GeneSymbol %in% target_genes)
accession_to_gene <- setNames(accession_map$GeneSymbol, accession_map$BaseAccession)
regulators <- regulators |>
  left_join(accession_map, by = "GeneSymbol")
if (any(is.na(regulators$BaseAccession))) {
  stop(
    "Missing UniProt mapping for regulator genes: ",
    paste(regulators$GeneSymbol[is.na(regulators$BaseAccession)], collapse = ", ")
  )
}
target_accessions <- unique(regulators$BaseAccession)

detection <- read.csv(
  detection_path,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  na.strings = c("", "NA")
)
sample_catalog <- detection |>
  distinct(
    PXD, SampleGroup, SampleGroupID, RowLabel,
    BiologicalMaterial, LactylomeProteinCount,
    GeneLevelAuditStatus, Parser, AuditReason
  ) |>
  mutate(RowOrder = row_number())
if (nrow(sample_catalog) != 40) {
  stop("Expected 40 sample groups, found ", nrow(sample_catalog))
}


quant_parts <- list()

add_quant <- function(
  data,
  pxd,
  sample_group,
  feature_values,
  target_values,
  sample_columns,
  sample_labels = sample_columns,
  measurement_level,
  quant_field,
  source_file
) {
  sample_columns <- intersect(sample_columns, names(data))
  if (!length(sample_columns) || !nrow(data)) return(invisible(NULL))
  if (length(sample_labels) != length(sample_columns)) {
    stop("sample_labels and sample_columns length mismatch for ", source_file)
  }
  base <- data.frame(
    FeatureID = as.character(feature_values),
    TargetAccession = as.character(target_values),
    stringsAsFactors = FALSE
  ) |>
    mutate(
      TargetGene = unname(accession_to_gene[TargetAccession])
    )
  values <- as.data.frame(lapply(data[sample_columns], safe_numeric))
  names(values) <- sample_labels
  long <- bind_cols(base, values) |>
    pivot_longer(
      cols = all_of(sample_labels),
      names_to = "QuantSample",
      values_to = "Signal"
    ) |>
    filter(
      !is.na(FeatureID),
      nzchar(FeatureID),
      is.finite(Signal),
      Signal > 0
    ) |>
    mutate(
      PXD = pxd,
      SampleGroup = sample_group,
      MeasurementLevel = measurement_level,
      QuantField = quant_field,
      SourceFile = relative_path(source_file, project_root),
      CanonicalFeature = ifelse(
        !is.na(TargetAccession) & nzchar(TargetAccession),
        paste0("ACC:", TargetAccession),
        FeatureID
      )
    ) |>
    group_by(
      PXD, SampleGroup, QuantSample, CanonicalFeature,
      TargetAccession, TargetGene,
      MeasurementLevel, QuantField, SourceFile
    ) |>
    summarise(Signal = sum(Signal, na.rm = TRUE), .groups = "drop")
  if (nrow(long)) quant_parts[[length(quant_parts) + 1]] <<- long
  invisible(NULL)
}


extract_maxquant <- function(
  path,
  pxd,
  sample_group,
  sample_tokens = NULL,
  sample_labels = sample_tokens,
  sheet = NULL
) {
  if (!file.exists(path)) return(invisible(NULL))
  data <- if (grepl("\\.xlsx$", path, ignore.case = TRUE)) {
    read_excel(path, sheet = sheet)
  } else {
    read_delimited(path)
  }
  if (!"Gene names" %in% names(data)) return(invisible(NULL))
  keep <- rep(TRUE, nrow(data))
  if ("Reverse" %in% names(data)) {
    keep <- keep & (is.na(data$Reverse) | data$Reverse != "+")
  }
  if ("Potential contaminant" %in% names(data)) {
    keep <- keep &
      (is.na(data$`Potential contaminant`) | data$`Potential contaminant` != "+")
  }
  if ("id" %in% names(data)) {
    keep <- keep & !is.na(data$id)
  }
  data <- data[keep, , drop = FALSE]
  if (!nrow(data)) return(invisible(NULL))

  if (is.null(sample_tokens)) {
    sample_columns <- "Intensity"
    sample_labels <- basename(dirname(dirname(path)))
    localization_columns <- "Localization prob"
  } else {
    sample_columns <- paste("Intensity", sample_tokens)
    localization_columns <- paste("Localization prob", sample_tokens)
  }
  valid <- sample_columns %in% names(data)
  sample_columns <- sample_columns[valid]
  sample_labels <- sample_labels[valid]
  localization_columns <- localization_columns[valid]
  for (index in seq_along(sample_columns)) {
    local_keep <- rep(TRUE, nrow(data))
    if (localization_columns[[index]] %in% names(data)) {
      local_keep <- safe_numeric(data[[localization_columns[[index]]]]) > 0
      local_keep[is.na(local_keep)] <- FALSE
    }
    subset <- data[local_keep, , drop = FALSE]
    add_quant(
      subset, pxd, sample_group,
      feature_values = accession_feature(subset$Proteins),
      target_values = match_target_accession(subset$Proteins, target_accessions),
      sample_columns = sample_columns[[index]],
      sample_labels = sample_labels[[index]],
      measurement_level = "Kla_site",
      quant_field = "MaxQuant site Intensity",
      source_file = path
    )
  }
}

# MaxQuant site tables.
pxd014_files <- list.files(
  file.path(project_root, "data", "PXD014870", "search_results"),
  pattern = "^Lactyl \\(K\\)Sites\\.txt$",
  recursive = TRUE,
  full.names = TRUE
)
for (path in pxd014_files) {
  extract_maxquant(path, "PXD014870", "MCF7")
}

extract_maxquant(
  file.path(
    project_root,
    "data/PXD033146/search_results/extracted_pairing/search_result-HA119TPLa/La (K)Sites.txt"
  ),
  "PXD033146", "pathological rotator cuff tendon"
)
extract_maxquant(
  file.path(project_root, "data/PXD036307/search_results/extracted/txt/La (K)Sites.txt"),
  "PXD036307", "normal human lung",
  c("PTB340", "PTB342", "PTB344", "PTB346", "PTB364", "PTB372")
)
extract_maxquant(
  file.path(project_root, "data/PXD050147/search_results/Lactyl_K_Sites.txt"),
  "PXD050147", "HepG2 WT and SIRT1 or SIRT3 KO",
  c(
    "SIRT1KO_Kla_rep1", "SIRT1KO_Kla_rep2", "SIRT1KO_Kla_rep3",
    "SIRT3KO_Kla_rep1", "SIRT3KO_Kla_rep2", "SIRT3KO_Kla_rep3",
    "WT_Kla_rep1", "WT_Kla_rep2", "WT_Kla_rep3"
  )
)
extract_maxquant(
  file.path(project_root, "data/PXD058534/search_results/extracted_pairing/txt/La (K)Sites.txt"),
  "PXD058534", "pretreated HK-2", "HK2"
)
for (spec in list(
  list("A", "MCF7"),
  list("B", "MDA-MB-468"),
  list("C", "MCF10A"),
  list("D", "T-47D")
)) {
  extract_maxquant(
    file.path(project_root, "data/PXD060185/search_results/RESULT/combined/txt/La (K)Sites.txt"),
    "PXD060185", spec[[2]], spec[[1]], spec[[2]]
  )
}
extract_maxquant(
  file.path(project_root, "data/PXD062720/search_results/extracted_pairing/txt/La (K)Sites.txt"),
  "PXD062720", "bladder cancer cells treated with EPI",
  c("A_1", "A_3")
)
extract_maxquant(
  file.path(project_root, "data/PXD063047/search_results/extracted/combined/txt/La (K)Sites.txt"),
  "PXD063047", "normal pregnancy placenta",
  c("Con_1", "Con_2", "Con_3")
)
extract_maxquant(
  file.path(project_root, "data/PXD063047/search_results/extracted/combined/txt/La (K)Sites.txt"),
  "PXD063047", "severe preeclampsia placenta",
  c("PE_1", "PE_2", "PE_3")
)
extract_maxquant(
  file.path(project_root, "data/PXD063266/search_results/LactylSites.xlsx"),
  "PXD063266", "PC-3M", c("1", "2", "3"),
  sheet = "Lactyl (K)Sites"
)
extract_maxquant(
  file.path(project_root, "data/PXD064038/search_results/extracted_pairing/txt/txt/La (K)Sites.txt"),
  "PXD064038", "MEC and NEC ESCC groups",
  c("MEC_1", "MEC_2", "MEC_3", "NEC_1", "NEC_2", "NEC_3")
)
extract_maxquant(
  file.path(project_root, "data/PXD078736/search_results/txt/La(K)Sites.txt"),
  "PXD078736", "HK-2 control and mannitol",
  c("ctr_1", "ctr_2", "ctr_3", "man_1", "man_2", "man_3")
)
extract_maxquant(
  file.path(
    project_root,
    "data/PXD028737/search_results/extracted_pairing/combined/combined/txt/La(K)Sites.txt"
  ),
  "PXD028737", "HMC3",
  c("H0", "H24"),
  c("HMC3_normoxia", "HMC3_hypoxia")
)

# PEAKS PXD028488: sum lactylated peptide areas to protein accession.
peaks_specs <- list(
  list(
    "HEK293T",
    "data/PXD028488/search_results/Enrichment-Search files/HEK293T-Enrichment-all HCD-Search files/protein-peptides.csv"
  ),
  list(
    "HCT116",
    "data/PXD028488/search_results/Enrichment-Search files/HCT116-Enrichment-Search files/protein-peptides.csv"
  ),
  list(
    "TALL-104",
    "data/PXD028488/search_results/Enrichment-Search files/TALL-NALAC-Search files/protein-peptides.csv"
  )
)
for (spec in peaks_specs) {
  path <- file.path(project_root, spec[[2]])
  data <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  keep <- grepl("Lactyl|\\+72\\.02", data$PTM, ignore.case = TRUE) |
    grepl("K\\(\\+72\\.02\\)|K\\(Lactyl\\)", data$Peptide, ignore.case = TRUE)
  data <- data[keep, , drop = FALSE]
  area_columns <- grep("^Area", names(data), value = TRUE)
  add_quant(
    data, "PXD028488", spec[[1]],
    feature_values = accession_feature(data$`Protein Accession`),
    target_values = match_target_accession(data$`Protein Accession`, target_accessions),
    sample_columns = area_columns,
    sample_labels = paste(spec[[1]], sub("^Area ?", "", area_columns), sep = "__"),
    measurement_level = "modified_peptide",
    quant_field = "PEAKS peptide Area",
    source_file = path
  )
}

# Proteome Discoverer normalized abundance from explicitly lactylated peptides.
path <- file.path(
  project_root,
  "data/PXD046800/search_results/HFX2_LFQ_QB001_Lacty_PeptideGroups.txt"
)
data <- read.delim(
  path, check.names = FALSE, stringsAsFactors = FALSE,
  quote = "\"", comment.char = ""
)
data <- data[
  grepl("Lacty|Lactyl|La \\(K\\)", data$Modifications, ignore.case = TRUE),
  ,
  drop = FALSE
]
for (spec in list(
  list("hypertrophic scar", "HSP"),
  list("adjacent skin", "NSP")
)) {
  columns <- grep(
    paste0('^"Abundances \\(Normalized\\):.*', spec[[2]]),
    names(data), value = TRUE
  )
  if (!length(columns)) {
    columns <- grep(
      paste0("^Abundances \\(Normalized\\):.*", spec[[2]]),
      names(data), value = TRUE
    )
  }
  add_quant(
    data, "PXD046800", spec[[1]],
    feature_values = accession_feature(data$`Master Protein Accessions`),
    target_values = match_target_accession(data$`Master Protein Accessions`, target_accessions),
    sample_columns = columns,
    sample_labels = sub("^.*Sample, |\"$", "", columns),
    measurement_level = "modified_peptide",
    quant_field = "Proteome Discoverer lactylated peptide normalized abundance",
    source_file = path
  )
}
rm(data)

# PXD053474: use enriched DDA peptide areas and enriched DIA PTM quantities.
dda_files <- list.files(
  file.path(
    project_root, "data", "PXD053474", "search_results", "extracted",
    "Subcellular-Whole-cell-lysates-Enriched-DDA-Search-files"
  ),
  pattern = "^protein-peptides\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)
for (path in dda_files) {
  data <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  keep <- grepl("Lactyl|\\+72\\.02", data$PTM, ignore.case = TRUE) |
    grepl("K\\(\\+72\\.02\\)|K\\(Lactyl\\)", data$Peptide, ignore.case = TRUE)
  data <- data[keep, , drop = FALSE]
  columns <- grep("^Area", names(data), value = TRUE)
  labels <- paste("DDA", basename(dirname(path)), sub("^Area ?", "", columns), sep = "__")
  add_quant(
    data, "PXD053474", "HCT116",
    feature_values = accession_feature(data$`Protein Accession`),
    target_values = match_target_accession(data$`Protein Accession`, target_accessions),
    sample_columns = columns,
    sample_labels = labels,
    measurement_level = "modified_peptide",
    quant_field = "PEAKS peptide Area",
    source_file = path
  )
}

dia_files <- list.files(
  file.path(
    project_root, "data", "PXD053474", "search_results", "extracted",
    "Subcellular-Whole-cell-lysates-Enriched-DIA-Search-files"
  ),
  pattern = "ptm-site_+Report\\.tsv$",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)
for (path in dia_files) {
  data <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  keep <- grepl("Lac|Lactyl", data$PTM.ModificationTitle, ignore.case = TRUE) &
    data$PTM.SiteAA == "K"
  probability_columns <- grep("PTM\\.SiteProbability$", names(data), value = TRUE)
  if (length(probability_columns)) {
    keep <- keep & rowSums(
      sapply(data[probability_columns], function(x) safe_numeric(x) > 0),
      na.rm = TRUE
    ) > 0
  }
  data <- data[keep, , drop = FALSE]
  columns <- grep("PTM\\.Quantity$", names(data), value = TRUE)
  add_quant(
    data, "PXD053474", "HCT116",
    feature_values = accession_feature(data$PTM.ProteinId),
    target_values = match_target_accession(data$PTM.ProteinId, target_accessions),
    sample_columns = columns,
    sample_labels = paste(
      "DIA",
      basename(dirname(path)),
      seq_along(columns),
      sep = "__"
    ),
    measurement_level = "Kla_site",
    quant_field = "Spectronaut PTM.Quantity",
    source_file = path
  )
}

# Author supplementary tables with explicit site quantities.
path <- file.path(
  project_root,
  "data/PXD050470/supplementary/prca2331-sup-0005-tables3.xlsx"
)
data <- read_excel(path, sheet = 1, skip = 12)
keep <- safe_numeric(data$`Localization probability`) > 0
data <- data[keep, , drop = FALSE]
columns <- grep("^Intensity_H", names(data), value = TRUE)
add_quant(
  data, "PXD050470", "human hippocampus",
  feature_values = accession_feature(data$`Proteins accession`),
  target_values = match_target_accession(data$`Proteins accession`, target_accessions),
  sample_columns = columns,
  sample_labels = sub("^Intensity_", "", columns),
  measurement_level = "Kla_site",
  quant_field = "author supplementary site Intensity",
  source_file = path
)
rm(data)

path <- file.path(
  project_root,
  "data/PXD054919/supplementary/41419_2025_8113_MOESM2_ESM.xlsx"
)
data <- read_excel(path, skip = 1)
columns <- grep("^Intensity A549-", names(data), value = TRUE)
add_quant(
  data, "PXD054919", "A549",
  feature_values = accession_feature(data$`Protein accession`),
  target_values = match_target_accession(data$`Protein accession`, target_accessions),
  sample_columns = columns,
  sample_labels = sub("^Intensity ", "", columns),
  measurement_level = "Kla_site",
  quant_field = "author site Intensity",
  source_file = path
)
rm(data)

path <- file.path(
  project_root,
  "data/PXD064912/supplementary/europepmc/mmc1.xlsx"
)
data <- read_excel(path, skip = 1)
keep <- tolower(data$PTM.ModificationTitle) == "lactylation" &
  data$PTM.SiteAA == "K"
data <- data[keep, , drop = FALSE]
columns <- grep("^PTM\\.Quantity", names(data), value = TRUE)
add_quant(
  data, "PXD064912", "human sperm",
  feature_values = accession_feature(data$PTM.ProteinId),
  target_values = match_target_accession(data$PTM.ProteinId, target_accessions),
  sample_columns = columns,
  sample_labels = sub("^PTM\\.Quantity\\[|\\]$", "", columns),
  measurement_level = "Kla_site",
  quant_field = "author PTM.Quantity",
  source_file = path
)
rm(data)

path <- file.path(
  project_root,
  "data/PXD070007/search_results/SA206LPLaB1_Annotation.xlsx"
)
data <- read_excel(path, sheet = "Annotation_Combine")
for (spec in list(
  list("glioblastoma stem cells", c("G2907", "G3028", "G3264", "GSC23", "MES28", "RKI")),
  list("neural stem cells", c("ENSA", "HMP1"))
)) {
  keep <- safe_numeric(data$`Localization probability`) > 0
  subset <- data[keep, , drop = FALSE]
  add_quant(
    subset, "PXD070007", spec[[1]],
    feature_values = accession_feature(subset$`Protein accession`),
    target_values = match_target_accession(subset$`Protein accession`, target_accessions),
    sample_columns = spec[[2]],
    measurement_level = "Kla_site",
    quant_field = "author annotation sample intensity",
    source_file = path
  )
}
rm(data)

path <- file.path(
  project_root,
  "data/PXD075377/search_results/extracted/2-Basic_analysis/MS_identified_information.txt"
)
data <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
for (spec in list(
  list("HCC", "Intensity HCC", "HCC"),
  list("adjacent liver", "Intensity Control", "Control")
)) {
  keep <- safe_numeric(data$`Localization probability`) > 0
  subset <- data[keep, , drop = FALSE]
  add_quant(
    subset, "PXD075377", spec[[1]],
    feature_values = accession_feature(subset$`Protein accession`),
    target_values = match_target_accession(subset$`Protein accession`, target_accessions),
    sample_columns = spec[[2]],
    sample_labels = spec[[3]],
    measurement_level = "Kla_site",
    quant_field = "author site Intensity",
    source_file = path
  )
}
rm(data)

# Spectronaut modified-peptide reports. Read only required columns from large files.
extract_spectronaut_precursor <- function(path, pxd, sample_group, require_site_confidence) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("data.table is required for large Spectronaut reports")
  }
  header <- names(data.table::fread(path, nrows = 0, data.table = FALSE))
  columns <- intersect(
    c(
      "Protein.Group", "Genes", "Modified.Sequence", "Q.Value", "PG.Q.Value",
      "PTM.Site.Confidence", "Run", "Precursor.Quantity"
    ),
    header
  )
  data <- data.table::fread(
    path, select = columns, showProgress = FALSE, data.table = FALSE
  )
  keep <- grepl("K\\(UniMod:378\\)", data$Modified.Sequence) &
    safe_numeric(data$Q.Value) <= 0.01 &
    safe_numeric(data$PG.Q.Value) <= 0.01
  if (require_site_confidence && "PTM.Site.Confidence" %in% names(data)) {
    keep <- keep & safe_numeric(data$PTM.Site.Confidence) > 0
  }
  data <- data[keep, , drop = FALSE]
  runs <- unique(data$Run)
  for (run in runs) {
    subset <- data[data$Run == run, , drop = FALSE]
    add_quant(
      subset, pxd, sample_group,
      feature_values = accession_feature(subset$Protein.Group),
      target_values = match_target_accession(subset$Protein.Group, target_accessions),
      sample_columns = "Precursor.Quantity",
      sample_labels = run,
      measurement_level = "modified_precursor",
      quant_field = "Spectronaut Precursor.Quantity",
      source_file = path
    )
  }
  rm(data)
  gc(verbose = FALSE)
}

extract_spectronaut_precursor(
  file.path(project_root, "data/PXD055230/search_results/LaIP_HSV1_DIA.tsv"),
  "PXD055230", "human fibroblasts mock and HCMV or HSV-1", TRUE
)
extract_spectronaut_precursor(
  file.path(project_root, "data/PXD057709/search_results/LaIP_report.tsv"),
  "PXD057709", "human fibroblasts mock and HCMV", FALSE
)

# Spectronaut PTM site report with condition labels.
path <- file.path(
  project_root,
  "data/PXD066054/search_results/extracted/PLa/XB08700B1DPLa_-PTMSiteReport.tsv"
)
data <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
keep <- data$PTM.ModificationTitle == "L-Lac(K)" &
  data$PTM.SiteAA == "K" &
  safe_numeric(data$PTM.SiteProbability) > 0
data <- data[keep, , drop = FALSE]
for (spec in list(
  list("BPH", "^NAT"),
  list("prostate cancer", "^PCa")
)) {
  subset <- data[grepl(spec[[2]], data$R.Condition), , drop = FALSE]
  for (condition in unique(subset$R.Condition)) {
    condition_data <- subset[subset$R.Condition == condition, , drop = FALSE]
    add_quant(
      condition_data, "PXD066054", spec[[1]],
      feature_values = accession_feature(condition_data$PTM.ProteinId),
      target_values = match_target_accession(condition_data$PTM.ProteinId, target_accessions),
      sample_columns = "PTM.Quantity",
      sample_labels = condition,
      measurement_level = "Kla_site",
      quant_field = "Spectronaut PTM.Quantity",
      source_file = path
    )
  }
}
rm(data)

# Explicit Lac(K) PTM quantities with accession-based target mapping.
path <- file.path(
  project_root,
  "data/PXD066351/search_results/XB01472B1DPLa-MSstats_Input.csv"
)
data <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
keep <- data$PTM.ModificationTitle == "Lac (K)" &
  data$PTM.SiteAA == "K" &
  safe_numeric(data$PTM.SiteProbability) > 0 &
  safe_numeric(data$PTM.Quantity) > 0
data <- data[keep, , drop = FALSE]
for (condition in unique(data$R.Condition)) {
  condition_data <- data[data$R.Condition == condition, , drop = FALSE]
  add_quant(
    condition_data, "PXD066351",
    "HCT116 control and Roseburia co-culture",
    feature_values = accession_feature(condition_data$PTM.ProteinId),
    target_values = match_target_accession(condition_data$PTM.ProteinId, target_accessions),
    sample_columns = "PTM.Quantity",
    sample_labels = condition,
    measurement_level = "Kla_site",
    quant_field = "MSstats Lac(K) PTM.Quantity",
    source_file = path
  )
}
rm(data)

# HUVEC Kla precursor matrix. UniMod:2114 is L-lactyllysine.
path <- file.path(
  project_root,
  paste0(
    "data/PXD073311/search_results/extracted_pairing/",
    "IPX0015307003_Database_search_result/Database_search_result/",
    "report.pr_matrix.tsv"
  )
)
data <- read.delim(
  path,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = ""
)
data <- data[
  grepl("K\\(UniMod:2114\\)", data$Modified.Sequence),
  ,
  drop = FALSE
]
columns <- names(data)[grepl("\\.raw$", names(data), ignore.case = TRUE)]
labels <- basename(gsub("\\\\", "/", columns))
labels <- sub("^FWMS20253547_02_", "", labels)
labels <- sub("\\.raw$", "", labels, ignore.case = TRUE)
add_quant(
  data, "PXD073311", "HUVEC control and Pg infection",
  feature_values = accession_feature(data$Protein.Group),
  target_values = match_target_accession(data$Protein.Group, target_accessions),
  sample_columns = columns,
  sample_labels = labels,
  measurement_level = "modified_precursor",
  quant_field = "Spectronaut precursor matrix quantity",
  source_file = path
)
rm(data)

# AC16 TMT Kla peptide supplementary table.
path <- file.path(
  project_root,
  "data/PXD075014/supplementary/Table2.XLSX"
)
data <- read_excel(path, sheet = "Kla peptides")
data <- data[
  grepl("lactylation", data$Modifications, ignore.case = TRUE),
  ,
  drop = FALSE
]
columns <- grep("^Abundance: F1:", names(data), value = TRUE)
labels <- c(
  "hypoxia_126", "hypoxia_127", "hypoxia_128",
  "control_129", "control_130", "control_131"
)
add_quant(
  data, "PXD075014", "AC16 control and hypoxia",
  feature_values = accession_feature(data$`Master Protein Accessions`),
  target_values = match_target_accession(data$`Master Protein Accessions`, target_accessions),
  sample_columns = columns,
  sample_labels = labels,
  measurement_level = "modified_peptide",
  quant_field = "author TMT peptide abundance",
  source_file = path
)
rm(data)

# PXD078013 has peptide evidence intensity and explicit Kla site IDs.
path <- file.path(project_root, "data/PXD078013/search_results/evidence.txt")
data <- read.delim(
  path, check.names = FALSE, stringsAsFactors = FALSE,
  quote = "", comment.char = ""
)
keep <- safe_numeric(data$`La (K)`) > 0 &
  !is.na(data$`La (K) site IDs`) &
  nzchar(as.character(data$`La (K) site IDs`)) &
  (is.na(data$Reverse) | data$Reverse != "+") &
  (is.na(data$`Potential contaminant`) | data$`Potential contaminant` != "+")
data <- data[keep, , drop = FALSE]
for (raw_file in unique(data$`Raw file`)) {
  subset <- data[data$`Raw file` == raw_file, , drop = FALSE]
  add_quant(
    subset, "PXD078013", "RKO WT and GSK3B KO",
    feature_values = accession_feature(subset$Proteins),
    target_values = match_target_accession(subset$Proteins, target_accessions),
    sample_columns = "Intensity",
    sample_labels = raw_file,
    measurement_level = "modified_peptide",
    quant_field = "MaxQuant evidence Intensity",
    source_file = path
  )
}
rm(data)

if (!length(quant_parts)) stop("No quantitative data were extracted")
quant_features <- bind_rows(quant_parts) |>
  group_by(
    PXD, SampleGroup, QuantSample, CanonicalFeature,
    TargetAccession, TargetGene,
    MeasurementLevel, QuantField
  ) |>
  summarise(
    Signal = sum(Signal, na.rm = TRUE),
    SourceFile = paste(sort(unique(SourceFile)), collapse = ";"),
    .groups = "drop"
  ) |>
  group_by(PXD, SampleGroup, QuantSample) |>
  mutate(
    Log2Signal = log2(Signal + 1),
    SampleMedianLog2 = median(Log2Signal, na.rm = TRUE),
    SampleNormalizedLog2 = Log2Signal - SampleMedianLog2,
    WithinSamplePercentile = {
      values <- Log2Signal
      if (length(values) == 1) {
        100
      } else {
        100 * (rank(values, ties.method = "average") - 1) /
          (length(values) - 1)
      }
    }
  ) |>
  ungroup()

sample_registry <- quant_features |>
  group_by(PXD, SampleGroup, QuantSample) |>
  summarise(
    SampleMedianLog2 = first(SampleMedianLog2),
    .groups = "drop"
  )

target_sample_grid <- sample_registry |>
  crossing(
    regulators |>
      distinct(GeneSymbol, BaseAccession) |>
      rename(RegulatorBaseAccession = BaseAccession)
  ) |>
  left_join(
    quant_features |>
      filter(!is.na(TargetAccession), TargetAccession %in% target_accessions) |>
      group_by(
        PXD, SampleGroup, QuantSample, TargetAccession, TargetGene
      ) |>
      summarise(
        Signal = sum(Signal),
        Log2Signal = log2(Signal + 1),
        SampleNormalizedLog2 = median(SampleNormalizedLog2),
        WithinSamplePercentile = max(WithinSamplePercentile),
        MeasurementLevel = paste(sort(unique(MeasurementLevel)), collapse = ";"),
        QuantField = paste(sort(unique(QuantField)), collapse = ";"),
        SourceFile = paste(sort(unique(SourceFile)), collapse = ";"),
        .groups = "drop"
      ),
    by = c(
      "PXD", "SampleGroup", "QuantSample",
      "RegulatorBaseAccession" = "TargetAccession"
    )
  ) |>
  mutate(
    Detected = !is.na(Signal) & Signal > 0,
    Signal = replace_na(Signal, 0),
    Log2Signal = replace_na(Log2Signal, 0),
    SampleNormalizedLog2 = ifelse(
      is.na(SampleNormalizedLog2),
      -SampleMedianLog2,
      SampleNormalizedLog2
    ),
    WithinSamplePercentile = replace_na(WithinSamplePercentile, 0)
  ) |>
  mutate(
    IdentityMatchMode = ifelse(
      Detected,
      "BaseAccession",
      "not_detected"
    )
  )

group_quant <- target_sample_grid |>
  group_by(PXD, SampleGroup, RegulatorBaseAccession, GeneSymbol) |>
  summarise(
    RelativeKlaPercentile = median(WithinSamplePercentile),
    DetectedSampleCount = sum(Detected),
    QuantSampleCount = n(),
    DetectedSampleFraction = mean(Detected),
    MedianLog2SignalDetected = ifelse(
      any(Detected),
      median(Log2Signal[Detected]),
      NA_real_
    ),
    MeasurementLevel = paste(sort(unique(na.omit(MeasurementLevel))), collapse = ";"),
    QuantField = paste(sort(unique(na.omit(QuantField))), collapse = ";"),
    SourceFile = paste(sort(unique(na.omit(SourceFile))), collapse = ";"),
    .groups = "drop"
  )

heatmap_data <- detection |>
  select(
    Role, GeneSymbol, PXD, SampleGroup, SampleGroupID, RowLabel,
    BiologicalMaterial, GeneLevelAuditStatus, Detected
  ) |>
  rename(IdentificationDetected = Detected) |>
  left_join(
    regulators |>
      select(Role, GeneSymbol, BaseAccession) |>
      rename(RegulatorBaseAccession = BaseAccession),
    by = c("Role", "GeneSymbol")
  ) |>
  left_join(
    group_quant,
    by = c("PXD", "SampleGroup", "GeneSymbol", "RegulatorBaseAccession")
  ) |>
  mutate(
    QuantificationAvailable = !is.na(QuantSampleCount),
    QuantState = case_when(
      !QuantificationAvailable ~ "无可比定量值",
      DetectedSampleCount == 0 ~ "定量样本中未检出",
      TRUE ~ "有定量信号"
    ),
    RelativeKlaPercentile = ifelse(
      QuantificationAvailable,
      replace_na(RelativeKlaPercentile, 0),
      NA_real_
    )
  )

write.csv(
  target_sample_grid,
  file.path(table_dir, "kla_regulator_intensity_sample_level_long.csv"),
  row.names = FALSE,
  na = ""
)
write.csv(
  heatmap_data,
  file.path(table_dir, "kla_regulator_normalized_intensity_long.csv"),
  row.names = FALSE,
  na = ""
)
id_mapping_audit <- target_sample_grid |>
  group_by(GeneSymbol, RegulatorBaseAccession) |>
  summarise(
    QuantitativeSampleRows = n(),
    DetectedSampleRows = sum(Detected),
    IdentityMatchMode = "UniProt_BaseAccession_only",
    GeneSymbolFallbackCount = 0L,
    MappingSource = "UniProtKB reviewed Homo sapiens",
    .groups = "drop"
  ) |>
  arrange(GeneSymbol)
write.csv(
  id_mapping_audit,
  file.path(table_dir, "kla_regulator_intensity_id_mapping_audit.csv"),
  row.names = FALSE,
  na = ""
)

quant_summary <- quant_features |>
  group_by(PXD, SampleGroup) |>
  summarise(
    QuantSampleCount = n_distinct(QuantSample),
    QuantifiedFeatureCount = n_distinct(CanonicalFeature),
    QuantifiedTargetGeneCount = n_distinct(TargetGene[!is.na(TargetGene)]),
    MeasurementLevel = paste(sort(unique(MeasurementLevel)), collapse = ";"),
    QuantField = paste(sort(unique(QuantField)), collapse = ";"),
    SourceFile = paste(sort(unique(SourceFile)), collapse = ";"),
    .groups = "drop"
  )

unavailable_reason <- function(pxd, status) {
  if (pxd == "PXD037371") {
    return("位点表含TMT reporter intensity，但当前缺少三类临床组与通道的可靠映射，不能拆分到对应组织组。")
  }
  if (status == "无逐蛋白明细") {
    return("当前文件不能把定量值可靠关联到Kla蛋白身份。")
  }
  "未从当前证据文件中提取到可解释的Kla定量字段，需进一步人工核对。"
}

audit <- sample_catalog |>
  left_join(quant_summary, by = c("PXD", "SampleGroup")) |>
  mutate(
    定量可用 = !is.na(QuantSampleCount),
    定量样本数 = replace_na(QuantSampleCount, 0L),
    定量层级 = ifelse(定量可用, MeasurementLevel, ""),
    使用的定量字段 = ifelse(定量可用, QuantField, ""),
    定量Kla特征数 = replace_na(QuantifiedFeatureCount, 0L),
    定量调控蛋白数 = replace_na(QuantifiedTargetGeneCount, 0L),
    调控蛋白身份匹配 = "UniProt BaseAccession only",
    GeneSymbol回退数 = 0L,
    跨研究热图数值 = ifelse(
      定量可用,
      "每个定量样本内log2(signal+1)后按全部Kla特征计算百分位，再对本样本组取中位数",
      ""
    ),
    可做PXD内Z分数 = 定量可用 & 定量样本数 >= 2,
    限制或不可用原因 = ifelse(
      定量可用,
      case_when(
        grepl(";", 定量层级, fixed = TRUE) ~
          "同一PXD包含多个测量层级；保留来源并仅用组内百分位汇总，不直接比较原始强度。",
        TRUE ~ "不同PXD的原始数值不可直接比较；跨研究只解释为组内相对Kla信号。"
      ),
      mapply(unavailable_reason, PXD, GeneLevelAuditStatus)
    ),
    定量来源文件 = ifelse(定量可用, SourceFile, "")
  ) |>
  transmute(
    PXD,
    样本组 = SampleGroup,
    行标签 = RowLabel,
    材料类型 = BiologicalMaterial,
    逐蛋白身份审计状态 = GeneLevelAuditStatus,
    定量可用,
    定量样本数,
    定量层级,
    使用的定量字段,
    定量Kla特征数,
    定量调控蛋白数,
    调控蛋白身份匹配,
    GeneSymbol回退数,
    跨研究热图数值,
    可做PXD内Z分数,
    限制或不可用原因,
    定量来源文件
  )
write.csv(
  audit,
  file.path(table_dir, "kla_regulator_intensity_availability_audit.csv"),
  row.names = FALSE,
  na = ""
)
write.csv(
  audit |>
    filter(!定量可用) |>
    select(PXD, 样本组, 限制或不可用原因, 定量来源文件),
  file.path(table_dir, "kla_regulator_intensity_plot_exclusions.csv"),
  row.names = FALSE,
  na = ""
)

# Main cross-study heatmap.
role_levels <- c("Writer", "Eraser", "Reader")
gene_levels <- regulators |>
  arrange(factor(Role, levels = role_levels), RoleEntryOrder) |>
  pull(GeneSymbol) |>
  unique()
row_levels <- rev(sample_catalog$RowLabel)
plot_data <- heatmap_data |>
  filter(QuantificationAvailable) |>
  mutate(
    Role = factor(Role, levels = role_levels),
    GeneSymbol = factor(GeneSymbol, levels = gene_levels),
    RowLabel = factor(RowLabel, levels = row_levels)
  )

main_plot <- ggplot(plot_data, aes(GeneSymbol, RowLabel, fill = RelativeKlaPercentile)) +
  geom_tile(color = "white", linewidth = 0.22) +
  geom_text(
    data = plot_data |> filter(!QuantificationAvailable),
    aes(label = "?"),
    color = "#4D4D4D",
    size = 2.4,
    fontface = "bold",
    family = plot_font
  ) +
  facet_grid(. ~ Role, scales = "free_x", space = "free_x") +
  scale_fill_gradientn(
    colours = c("#FFFFFF", "#FFF3E0", "#FDBB84", "#FC8D59", "#B2182B"),
    values = scales::rescale(c(0, 20, 50, 80, 100)),
    limits = c(0, 100),
    na.value = "#D9D9D9",
    name = "组内Kla信号\n百分位"
  ) +
  labs(
    title = "乳酸化调控蛋白在细胞系与组织中的相对Kla信号",
    subtitle = paste(
      "每个定量样本内按全部Kla特征计算信号百分位，再在样本组内取中位数；",
      "信号由白色向暖色递增；0表示定量样本中未检出；无可靠定量值的样本组不进入本图。",
      "该数值不是表达量或Log FC。"
    ),
    x = NULL,
    y = NULL,
    caption = paste0(
      "调控蛋白身份仅按人源UniProt reviewed BaseAccession匹配；",
      "GeneSymbol只作为图中显示标签，不参与蛋白命中或强度汇总。"
    )
  ) +
  theme_minimal(base_size = 8.5, base_family = plot_font) +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold", size = 9),
    strip.background = element_rect(fill = "#F2F2F2", color = NA),
    axis.text.x = element_text(angle = 55, hjust = 1, vjust = 1, size = 7),
    axis.text.y = element_text(size = 7.2),
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 8.2, lineheight = 1.12),
    plot.caption = element_text(size = 7.5, hjust = 0, color = "#5C626A"),
    legend.position = "bottom",
    legend.key.width = grid::unit(35, "mm"),
    plot.margin = margin(8, 10, 8, 8)
  )

ggsave(
  file.path(figure_dir, "kla_regulator_cross_study_relative_intensity_heatmap.png"),
  main_plot,
  width = 15.5,
  height = 11.5,
  dpi = 320,
  bg = "white"
)
ggsave(
  file.path(figure_dir, "kla_regulator_cross_study_relative_intensity_heatmap.pdf"),
  main_plot,
  width = 15.5,
  height = 11.5,
  device = cairo_pdf,
  bg = "white"
)

# Within-PXD gene-wise z-scores are valid only across samples from the same PXD.
zscore_data <- target_sample_grid |>
  group_by(PXD, GeneSymbol) |>
  mutate(
    WithinPXDZ = {
      values <- SampleNormalizedLog2
      value_sd <- sd(values, na.rm = TRUE)
      if (
        n_distinct(QuantSample) >= 2 &&
          sum(Detected, na.rm = TRUE) >= 2 &&
          is.finite(value_sd) &&
          value_sd > 0
      ) {
        as.numeric(scale(values))
      } else {
        rep(NA_real_, length(values))
      }
    }
  ) |>
  ungroup() |>
  filter(!is.na(WithinPXDZ))
write.csv(
  zscore_data,
  file.path(table_dir, "kla_regulator_within_pxd_zscore_long.csv"),
  row.names = FALSE,
  na = ""
)

z_pxd <- zscore_data |>
  group_by(PXD) |>
  summarise(
    Samples = n_distinct(QuantSample),
    Genes = n_distinct(GeneSymbol),
    .groups = "drop"
  ) |>
  filter(Samples >= 2, Genes >= 2) |>
  pull(PXD)

pdf_path <- file.path(figure_dir, "kla_regulator_within_pxd_zscore_heatmaps.pdf")
grDevices::cairo_pdf(pdf_path, width = 12, height = 8.5, onefile = TRUE)
for (pxd in z_pxd) {
  data <- zscore_data |>
    filter(PXD == pxd) |>
    mutate(
      GeneSymbol = factor(GeneSymbol, levels = gene_levels),
      QuantSample = factor(QuantSample, levels = rev(unique(QuantSample)))
    )
  p <- ggplot(data, aes(GeneSymbol, QuantSample, fill = pmax(-2.5, pmin(2.5, WithinPXDZ)))) +
    geom_tile(color = "white", linewidth = 0.25) +
    scale_fill_gradient2(
      low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
      midpoint = 0, limits = c(-2.5, 2.5),
      name = "PXD内Z分数"
    ) +
    labs(
      title = paste0(pxd, "：乳酸化调控蛋白的PXD内相对变化"),
      subtitle = "log2(signal+1)后先做样本中位数中心化，再对每个蛋白跨本PXD样本计算Z分数；不跨PXD比较。",
      x = NULL,
      y = NULL
    ) +
    theme_minimal(base_size = 9, base_family = plot_font) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 55, hjust = 1, size = 7),
      axis.text.y = element_text(size = 7),
      plot.title = element_text(face = "bold"),
      legend.position = "bottom"
    )
  print(p)
}
grDevices::dev.off()

message(
  "Quantitative sample groups: ",
  sum(audit$定量可用),
  "/",
  nrow(audit),
  "; within-PXD heatmaps: ",
  length(z_pxd)
)
})

# --- Stage: whole_proteome (from analyze_kla_regulator_whole_proteome_intensity.R) --
run_stage("whole_proteome", {
# Stage-private base_accession: the lib version (accession_utils.R) lacks the
# left-token rule sub("^([^|;]+)\\|.*$", "\\1") that the original
# analyze_kla_regulator_whole_proteome_intensity.R private version had. That rule
# is REQUIRED by this project's data: PEAKS proteins.csv Accession values are
# "P78527|PRKDC_HUMAN" (accession|gene-name, no sp| prefix) and the original
# script extracted the left token "P78527" to match regulator accessions; without
# it the merged run mapped 0 regulators for PXD028488/PXD053474 (original:
# 16/13/16/32) and its tables no longer matched the baseline. This private
# definition masks the lib one from this point on (whole_proteome is the last
# stage, so nothing downstream is affected).
base_accession <- function(values) {
  values <- as.character(values)
  values <- sub("^.*\\|([^|]+)\\|.*$", "\\1", values)
  values <- sub("^([^|;]+)\\|.*$", "\\1", values)
  sub("-[0-9]+$", "", values)
}
regulator_path <- file.path(
  project_root, "data", "identifier",
  "乳酸化调控因子_Writer-Eraser-Reader.xlsx"
)
mapping_path <- file.path(
  project_root, "reanalysis", "config",
  "lactylation_regulator_uniprot_mapping.csv"
)
detection_path <- file.path(
  table_dir, "kla_regulator_cell_type_long.csv"
)
reviewed_uniprot_path <- file.path(
  project_root, "reanalysis", "config",
  "uniprot_human_reviewed_2026-08-05.tsv"
)

required_files <- c(
  regulator_path,
  mapping_path,
  detection_path,
  reviewed_uniprot_path
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) {
  stop("Missing required files: ", paste(missing_files, collapse = ", "))
}

regulators <- read_excel(regulator_path) |>
  transmute(
    Role = trimws(as.character(分组)),
    GeneSymbol = toupper(trimws(as.character(基因))),
    RoleEntryOrder = row_number()
  ) |>
  filter(
    Role %in% c("Writer", "Eraser", "Reader"),
    !is.na(GeneSymbol),
    nzchar(GeneSymbol)
  )

accession_map <- read.csv(
  mapping_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
) |>
  transmute(
    BaseAccession = sub("-[0-9]+$", "", BaseAccession),
    GeneSymbol = toupper(GeneSymbol)
  )
accession_to_gene <- setNames(accession_map$GeneSymbol, accession_map$BaseAccession)
regulators <- regulators |>
  left_join(accession_map, by = "GeneSymbol")
if (any(is.na(regulators$BaseAccession))) {
  stop(
    "Missing UniProt mapping for: ",
    paste(regulators$GeneSymbol[is.na(regulators$BaseAccession)], collapse = ", ")
  )
}
target_accessions <- unique(regulators$BaseAccession)

detection <- read.csv(
  detection_path,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  na.strings = c("", "NA")
)
sample_catalog <- detection |>
  distinct(
    PXD, SampleGroup, SampleGroupID, RowLabel,
    BiologicalMaterial, GeneLevelAuditStatus
  ) |>
  mutate(RowOrder = row_number())


reviewed_uniprot <- data.table::fread(
  reviewed_uniprot_path,
  sep = "\t",
  quote = "",
  data.table = FALSE,
  check.names = FALSE,
  showProgress = FALSE
)
reviewed_primary_summary <- reviewed_uniprot |>
  transmute(
    BaseAccession = base_accession(Entry),
    Symbols = `Gene Names (primary)`
  ) |>
  separate_rows(Symbols, sep = "[;[:space:]]+") |>
  transmute(
    SourceSymbol = toupper(trimws(Symbols)),
    BaseAccession
  ) |>
  filter(nzchar(SourceSymbol)) |>
  distinct() |>
  group_by(SourceSymbol) |>
  summarise(
    CandidateCount = n_distinct(BaseAccession),
    CandidateAccessions = paste(
      sort(unique(BaseAccession)),
      collapse = ";"
    ),
    .groups = "drop"
  )
reviewed_alias_summary <- reviewed_uniprot |>
  transmute(
    BaseAccession = base_accession(Entry),
    Symbols = `Gene Names`
  ) |>
  separate_rows(Symbols, sep = "[;[:space:]]+") |>
  transmute(
    SourceSymbol = toupper(trimws(Symbols)),
    BaseAccession
  ) |>
  filter(nzchar(SourceSymbol)) |>
  distinct() |>
  group_by(SourceSymbol) |>
  summarise(
    CandidateCount = n_distinct(BaseAccession),
    CandidateAccessions = paste(
      sort(unique(BaseAccession)),
      collapse = ";"
    ),
    .groups = "drop"
  )

map_reviewed_symbol_features <- function(source_features) {
  data.frame(
    SourceFeature = as.character(source_features),
    stringsAsFactors = FALSE
  ) |>
    mutate(SourceSymbol = strsplit(SourceFeature, "_", fixed = TRUE)) |>
    unnest(SourceSymbol) |>
    mutate(SourceSymbol = toupper(trimws(SourceSymbol))) |>
    filter(nzchar(SourceSymbol)) |>
    distinct() |>
    left_join(
      reviewed_primary_summary |>
        rename(
          PrimaryCandidateCount = CandidateCount,
          PrimaryCandidates = CandidateAccessions
        ),
      by = "SourceSymbol"
    ) |>
    left_join(
      reviewed_alias_summary |>
        rename(
          AliasCandidateCount = CandidateCount,
          AliasCandidates = CandidateAccessions
        ),
      by = "SourceSymbol"
    ) |>
    mutate(
      MappingMode = case_when(
        !is.na(PrimaryCandidateCount) & PrimaryCandidateCount == 1 ~
          "reviewed_primary_symbol",
        is.na(PrimaryCandidateCount) & AliasCandidateCount == 1 ~
          "reviewed_unique_alias",
        !is.na(PrimaryCandidateCount) & PrimaryCandidateCount > 1 ~
          "ambiguous_reviewed_primary",
        AliasCandidateCount > 1 ~ "ambiguous_reviewed_alias",
        TRUE ~ "unmapped"
      ),
      BaseAccession = case_when(
        MappingMode == "reviewed_primary_symbol" ~ PrimaryCandidates,
        MappingMode == "reviewed_unique_alias" ~ AliasCandidates,
        TRUE ~ ""
      ),
      TargetAccession = ifelse(
        BaseAccession %in% target_accessions,
        BaseAccession,
        NA_character_
      ),
      MappingSource =
        "UniProtKB reviewed Homo sapiens snapshot 2026-08-05"
    )
}

quant_parts <- list()
quant_audit <- sample_catalog |>
  transmute(
    PXD,
    SampleGroup,
    SampleGroupID,
    WholeProteomeQuantAvailable = FALSE,
    WholeProteomeSource = "",
    WholeProteomeMeasurement = "",
    WholeProteomeStatus = "unavailable",
    WholeProteomeReason =
      "未配置同研究/同样本全蛋白强度文件；不使用Kla信号替代"
  )

update_audit <- function(
  pxd,
  sample_group,
  source_file,
  measurement,
  reason = "同研究或同样本普通全蛋白定量可用"
) {
  rows <- quant_audit$PXD == pxd & quant_audit$SampleGroup == sample_group
  if (!any(rows)) return(invisible(NULL))
  quant_audit$WholeProteomeQuantAvailable[rows] <<- TRUE
  old_sources <- quant_audit$WholeProteomeSource[rows]
  old_measurements <- quant_audit$WholeProteomeMeasurement[rows]
  quant_audit$WholeProteomeSource[rows] <<- vapply(
    old_sources,
    function(old) paste(
      sort(unique(c(
        strsplit(old, ";", fixed = TRUE)[[1]][nzchar(old)],
        source_file
      ))),
      collapse = ";"
    ),
    character(1)
  )
  quant_audit$WholeProteomeMeasurement[rows] <<- vapply(
    old_measurements,
    function(old) paste(
      sort(unique(c(
        strsplit(old, ";", fixed = TRUE)[[1]][nzchar(old)],
        measurement
      ))),
      collapse = ";"
    ),
    character(1)
  )
  quant_audit$WholeProteomeStatus[rows] <<- "available"
  quant_audit$WholeProteomeReason[rows] <<- reason
  invisible(NULL)
}

add_total_quant <- function(
  data,
  pxd,
  sample_group,
  feature_values,
  target_values,
  sample_columns,
  sample_labels = sample_columns,
  measurement,
  source_file
) {
  sample_columns <- intersect(sample_columns, names(data))
  if (!length(sample_columns) || !nrow(data)) return(invisible(NULL))
  if (length(sample_labels) != length(sample_columns)) {
    stop("sample_labels and sample_columns length mismatch for ", source_file)
  }
  base <- data.frame(
    FeatureID = as.character(feature_values),
    TargetAccession = as.character(target_values),
    stringsAsFactors = FALSE
  ) |>
    mutate(
      TargetGene = unname(accession_to_gene[TargetAccession])
    )
  values <- as.data.frame(lapply(data[sample_columns], safe_numeric))
  names(values) <- sample_labels
  long <- bind_cols(base, values) |>
    pivot_longer(
      cols = all_of(sample_labels),
      names_to = "QuantSample",
      values_to = "Signal"
    ) |>
    filter(
      !is.na(FeatureID),
      nzchar(FeatureID),
      is.finite(Signal),
      Signal > 0
    ) |>
    mutate(
      PXD = pxd,
      SampleGroup = sample_group,
      Measurement = measurement,
      SourceFile = relative_path(source_file, project_root),
      CanonicalFeature = ifelse(
        !is.na(TargetAccession) & nzchar(TargetAccession),
        paste0("ACC:", TargetAccession),
        FeatureID
      )
    ) |>
    group_by(
      PXD, SampleGroup, QuantSample, CanonicalFeature,
      TargetAccession, TargetGene, Measurement, SourceFile
    ) |>
    summarise(Signal = sum(Signal, na.rm = TRUE), .groups = "drop")
  if (nrow(long)) {
    quant_parts[[length(quant_parts) + 1]] <<- long
    update_audit(
      pxd,
      sample_group,
      relative_path(source_file, project_root),
      measurement
    )
  }
  invisible(NULL)
}

add_maxquant_proteome <- function(
  path,
  pxd,
  sample_group,
  token_map,
  measurement = "MaxQuant proteinGroups"
) {
  if (!file.exists(path)) return(invisible(NULL))
  data <- read.delim(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = ""
  )
  keep <- rep(TRUE, nrow(data))
  if ("Reverse" %in% names(data)) keep <- keep & data$Reverse != "+"
  if ("Potential contaminant" %in% names(data)) {
    keep <- keep & data$`Potential contaminant` != "+"
  }
  if ("Only identified by site" %in% names(data)) {
    keep <- keep & data$`Only identified by site` != "+"
  }
  if ("id" %in% names(data)) keep <- keep & !is.na(data$id)
  data <- data[keep, , drop = FALSE]
  accession_col <- intersect(
    c("Majority protein IDs", "Protein IDs"),
    names(data)
  )
  if (!length(accession_col)) return(invisible(NULL))
  accession_col <- accession_col[[1]]
  for (spec in token_map) {
    token <- spec[[1]]
    label <- spec[[2]]
    candidate_columns <- c(
      token,
      paste("LFQ intensity", token),
      paste("Intensity", token),
      paste("iBAQ", token)
    )
    column <- candidate_columns[candidate_columns %in% names(data)]
    if (!length(column)) next
    add_total_quant(
      data,
      pxd,
      sample_group,
      accession_feature(data[[accession_col]]),
      match_target_accession(data[[accession_col]], target_accessions),
      column[[1]],
      label,
      measurement,
      path
    )
  }
}

add_pd_proteome <- function(path, pxd, group_map) {
  if (!file.exists(path)) return(invisible(NULL))
  data <- read.delim(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    quote = "\"",
    comment.char = ""
  )
  for (spec in group_map) {
    group <- spec[[1]]
    token <- spec[[2]]
    columns <- grep(
      paste0("^Abundances \\(Normalized\\):.*", token),
      names(data),
      value = TRUE
    )
    add_total_quant(
      data,
      pxd,
      group,
      accession_feature(data$Accession),
      match_target_accession(data$Accession, target_accessions),
      columns,
      sub("^.*Sample, ", "", columns),
      "Proteome Discoverer normalized total-protein abundance",
      path
    )
  }
}

add_spectronaut_report <- function(path, pxd, sample_group, measurement) {
  if (!file.exists(path)) return(invisible(NULL))
  data <- data.table::fread(path, data.table = FALSE, showProgress = FALSE)
  if (!all(c("Protein.Group", "PG.Quantity") %in% names(data))) {
    return(invisible(NULL))
  }
  if ("PG.Q.Value" %in% names(data)) {
    data <- data[safe_numeric(data$PG.Q.Value) <= 0.01, , drop = FALSE]
  }
  run_column <- if ("Run" %in% names(data)) "Run" else "File.Name"
  for (run in unique(data[[run_column]])) {
    subset <- data[data[[run_column]] == run, , drop = FALSE] |>
      group_by(Protein.Group) |>
      summarise(PG.Quantity = max(safe_numeric(PG.Quantity)), .groups = "drop")
    add_total_quant(
      subset,
      pxd,
      sample_group,
      accession_feature(subset$Protein.Group),
      match_target_accession(subset$Protein.Group, target_accessions),
      "PG.Quantity",
      as.character(run),
      measurement,
      path
    )
  }
}

add_spectronaut_standard_report <- function(
  path,
  pxd,
  sample_group,
  measurement
) {
  if (!file.exists(path)) return(invisible(NULL))
  data <- data.table::fread(path, data.table = FALSE, showProgress = FALSE)
  accession_column <- intersect(
    c("PG.ProteinAccessions", "Protein.Group"),
    names(data)
  )
  quantity_columns <- grep("PG\\.Quantity$", names(data), value = TRUE)
  if (!length(accession_column) || !length(quantity_columns)) {
    return(invisible(NULL))
  }
  accession_column <- accession_column[[1]]
  for (column in quantity_columns) {
    sample_label <- sub("^\\[[0-9]+\\] ", "", column)
    sample_label <- sub("\\.PG\\.Quantity$", "", sample_label)
    subset <- data |>
      transmute(
        ProteinAccessions = .data[[accession_column]],
        Quantity = safe_numeric(.data[[column]])
      ) |>
      group_by(ProteinAccessions) |>
      summarise(Quantity = max(Quantity, na.rm = TRUE), .groups = "drop")
    add_total_quant(
      subset,
      pxd,
      sample_group,
      accession_feature(subset$ProteinAccessions),
      match_target_accession(subset$ProteinAccessions, target_accessions),
      "Quantity",
      sample_label,
      measurement,
      path
    )
  }
}

add_spectronaut_matrix <- function(
  path,
  pxd,
  sample_group,
  measurement,
  column_pattern = NULL
) {
  if (!file.exists(path)) return(invisible(NULL))
  data <- read.delim(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = ""
  )
  columns <- grep("PG\\.Quantity$", names(data), value = TRUE)
  if (!is.null(column_pattern)) {
    columns <- columns[grepl(column_pattern, columns)]
  }
  if (!length(columns)) return(invisible(NULL))
  add_total_quant(
    data,
    pxd,
    sample_group,
    accession_feature(data$PG.ProteinGroups),
    match_target_accession(data$PG.ProteinGroups, target_accessions),
    columns,
    sub("^.*\\] |\\.PG\\.Quantity$", "", columns),
    measurement,
    path
  )
}

add_peaks_proteins <- function(paths, pxd, sample_group) {
  for (path in paths) {
    data <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
    area_columns <- grep("^Area", names(data), value = TRUE)
    if (!length(area_columns) || !"Accession" %in% names(data)) next
    add_total_quant(
      data,
      pxd,
      sample_group,
      accession_feature(data$Accession),
      match_target_accession(data$Accession, target_accessions),
      area_columns,
      paste0(sample_group, "__", basename(dirname(path)), "__", sub("^Area ", "", area_columns)),
      "PEAKS non-enriched protein Area",
      path
    )
  }
}

add_pxd043880_hippocampus <- function(path) {
  if (!file.exists(path)) return(invisible(NULL))
  data <- read_excel(
    path,
    sheet = "Source Data Proteins",
    col_names = FALSE,
    col_types = "text",
    .name_repair = "minimal"
  )
  if (nrow(data) < 3 || ncol(data) < 9) {
    stop("Unexpected PXD043880 protein matrix layout: ", path)
  }

  donor_rows <- 3:nrow(data)
  donor_ids <- trimws(as.character(data[[1]][donor_rows]))
  donor_keep <- !is.na(donor_ids) &
    nzchar(donor_ids) &
    !grepl("\\((removed|outlier)\\)", donor_ids, ignore.case = TRUE)
  donor_rows <- donor_rows[donor_keep]
  donor_ids <- donor_ids[donor_keep]
  if (length(donor_ids) != 74) {
    stop("Expected 74 retained PXD043880 donors, found ", length(donor_ids))
  }

  source_features <- trimws(as.character(
    unlist(data[2, 9:ncol(data)], use.names = FALSE)
  ))
  feature_keep <- !is.na(source_features) & nzchar(source_features)
  source_features <- source_features[feature_keep]
  intensity <- as.matrix(
    data[donor_rows, 9:ncol(data), drop = FALSE]
  )[, feature_keep, drop = FALSE]
  quant_data <- as.data.frame(
    t(intensity),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  sample_labels <- make.unique(paste0("CA1_", donor_ids))
  names(quant_data) <- sample_labels

  mapping <- map_reviewed_symbol_features(source_features)
  target_by_feature <- mapping |>
    group_by(SourceFeature) |>
    summarise(
      TargetHitCount = n_distinct(TargetAccession[!is.na(TargetAccession)]),
      TargetAccession = paste(
        sort(unique(TargetAccession[!is.na(TargetAccession)])),
        collapse = ";"
      ),
      .groups = "drop"
    )
  if (any(target_by_feature$TargetHitCount > 1)) {
    stop(
      "PXD043880 source feature maps to multiple regulator accessions: ",
      paste(
        target_by_feature$SourceFeature[
          target_by_feature$TargetHitCount > 1
        ],
        collapse = ", "
      )
    )
  }
  feature_targets <- target_by_feature$TargetAccession[
    match(source_features, target_by_feature$SourceFeature)
  ]
  feature_targets[!nzchar(feature_targets)] <- NA_character_

  mapping |>
    mutate(
      SourceFile = relative_path(path, project_root),
      RetainedDonorCount = length(donor_ids),
      TargetRegulatorMatch = !is.na(TargetAccession)
    ) |>
    write.csv(
      file.path(
        table_dir,
        "kla_regulator_whole_proteome_hippocampus_id_mapping_audit.csv"
      ),
      row.names = FALSE,
      na = ""
    )

  add_total_quant(
    quant_data,
    "PXD050470",
    "human hippocampus",
    paste0("PXD043880_SYMBOL_FEATURE:", source_features),
    feature_targets,
    sample_labels,
    sample_labels,
    paste0(
      "PXD043880 CA1 LFQ intensity after reviewed UniProt ",
      "symbol conversion"
    ),
    path
  )
  update_audit(
    "PXD050470",
    "human hippocampus",
    relative_path(path, project_root),
    paste0(
      "PXD043880 CA1 LFQ intensity after reviewed UniProt ",
      "symbol conversion"
    ),
    paste0(
      "独立正常CA1海马全蛋白参照；74名供体；",
      "symbol先转换为人源reviewed UniProt BaseAccession"
    )
  )
}

# Same-study total-proteome protein quantities. Kla-IP-only files are excluded.
add_pxd043880_hippocampus(file.path(
  project_root,
  "data/PXD043880/supplementary/13024_2023_650_MOESM1_ESM.xlsx"
))

add_peaks_proteins(
  list.files(
    file.path(project_root, "data/PXD028488/search_results/Nonenrichment-Search files"),
    pattern = "^proteins\\.csv$",
    recursive = TRUE,
    full.names = TRUE
  )[grepl("HEK293T", list.files(
    file.path(project_root, "data/PXD028488/search_results/Nonenrichment-Search files"),
    pattern = "^proteins\\.csv$",
    recursive = TRUE,
    full.names = TRUE
  ))],
  "PXD028488",
  "HEK293T"
)
add_peaks_proteins(
  list.files(
    file.path(project_root, "data/PXD028488/search_results/Nonenrichment-Search files"),
    pattern = "^proteins\\.csv$",
    recursive = TRUE,
    full.names = TRUE
  )[grepl("HCT116", list.files(
    file.path(project_root, "data/PXD028488/search_results/Nonenrichment-Search files"),
    pattern = "^proteins\\.csv$",
    recursive = TRUE,
    full.names = TRUE
  ))],
  "PXD028488",
  "HCT116"
)
add_peaks_proteins(
  list.files(
    file.path(project_root, "data/PXD028488/search_results/Nonenrichment-Search files"),
    pattern = "^proteins\\.csv$",
    recursive = TRUE,
    full.names = TRUE
  )[grepl("TALL", list.files(
    file.path(project_root, "data/PXD028488/search_results/Nonenrichment-Search files"),
    pattern = "^proteins\\.csv$",
    recursive = TRUE,
    full.names = TRUE
  ))],
  "PXD028488",
  "TALL-104"
)

maxquant_jobs <- list(
  list(
    "PXD033146", "pathological rotator cuff tendon",
    "data/PXD033146/search_results/extracted_pairing/search_result-HA119TQ/proteinGroups.txt",
    lapply(1:6, function(i) list(
      paste("Reporter intensity corrected", i),
      paste0("channel_", i)
    ))
  ),
  list(
    "PXD050147", "HepG2 WT and SIRT1 or SIRT3 KO",
    "data/PXD050147/search_results/SIRT_proteinGroups.txt",
    lapply(
      c(
        "SIRT1KO_pro_rep1", "SIRT1KO_pro_rep2", "SIRT1KO_pro_rep3",
        "SIRT3KO_pro_rep1", "SIRT3KO_pro_rep2", "SIRT3KO_pro_rep3",
        "WT_pro_rep1", "WT_pro_rep2", "WT_pro_rep3"
      ),
      function(x) list(x, x)
    )
  ),
  list(
    "PXD062720", "bladder cancer cells treated with EPI",
    "data/PXD062720/search_results/extracted_pairing/txt/proteinGroups.txt",
    list(list("A_1", "A_1"), list("A_3", "A_3"))
  )
)
for (job in maxquant_jobs) {
  add_maxquant_proteome(
    file.path(project_root, job[[3]]),
    job[[1]],
    job[[2]],
    job[[4]]
  )
}

add_pd_proteome(
  file.path(project_root, "data/PXD046800/search_results/HFX2_LFQ_QB002_Proteins.txt"),
  "PXD046800",
  list(
    list("hypertrophic scar", "HSP"),
    list("adjacent skin", "NSP")
  )
)

add_spectronaut_report(
  file.path(project_root, "data/PXD055230/search_results/WP_HSV1_DIA.tsv"),
  "PXD055230",
  "human fibroblasts mock and HCMV or HSV-1",
  "Spectronaut whole-proteome PG.Quantity"
)
add_spectronaut_report(
  file.path(project_root, "data/PXD057709/search_results/WP_report.tsv"),
  "PXD057709",
  "human fibroblasts mock and HCMV",
  "Spectronaut whole-proteome PG.Quantity"
)
add_spectronaut_matrix(
  file.path(project_root, "data/PXD066054/search_results/extracted/DA/Protein_Quant.tsv"),
  "PXD066054",
  "BPH",
  "Spectronaut whole-proteome PG.Quantity",
  "NAT"
)
add_spectronaut_matrix(
  file.path(project_root, "data/PXD066054/search_results/extracted/DA/Protein_Quant.tsv"),
  "PXD066054",
  "prostate cancer",
  "Spectronaut whole-proteome PG.Quantity",
  "PCa"
)
add_spectronaut_matrix(
  file.path(project_root, "data/PXD066351/search_results/XB01472B1DA-Protein_Quant.tsv"),
  "PXD066351",
  "HCT116 control and Roseburia co-culture",
  "Spectronaut whole-proteome PG.Quantity"
)

add_peaks_proteins(
  list.files(
    file.path(
      project_root,
      "data/PXD053474/search_results/extracted",
      "Subcellular-Whole-cell-lysates-Unenriched-DDA-Search-files"
    ),
    pattern = "^proteins\\.csv$",
    recursive = TRUE,
    full.names = TRUE
  )[grepl("Whole-cell-lysates-DDA-input", list.files(
    file.path(
      project_root,
      "data/PXD053474/search_results/extracted",
      "Subcellular-Whole-cell-lysates-Unenriched-DDA-Search-files"
    ),
    pattern = "^proteins\\.csv$",
    recursive = TRUE,
    full.names = TRUE
  ))],
  "PXD053474",
  "HCT116"
)
add_spectronaut_standard_report(
  file.path(
    project_root,
    "data/PXD053474/search_results/extracted",
    "Subcellular-Whole-cell-lysates-Unenriched-DIA-Search-files",
    "Subcellular-Unenriched-DIA-Search files",
    "Whole-cell-lysates-DIA-input",
    "Whole-cell-DIA-input-standard-Report.tsv"
  ),
  "PXD053474",
  "HCT116",
  "Spectronaut whole-proteome PG.Quantity"
)

pxd073311_total_path <- file.path(
  project_root,
  "data/PXD073311/search_results/extracted_pairing",
  "IPX0015307001_Database_search_result/Database_search_result",
  "report.pg_matrix.tsv"
)
pxd073311_total <- read.delim(
  pxd073311_total_path,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = ""
)
pxd073311_total_columns <- grep(
  "\\.raw$",
  names(pxd073311_total),
  value = TRUE,
  ignore.case = TRUE
)
add_total_quant(
  pxd073311_total,
  "PXD073311",
  "HUVEC control and Pg infection",
  accession_feature(pxd073311_total$Protein.Group),
  match_target_accession(pxd073311_total$Protein.Group, target_accessions),
  pxd073311_total_columns,
  c("A0h_1", "A0h_2", "A0h_3", "A6h_1", "A6h_2", "A6h_3"),
  "Spectronaut whole-proteome PG matrix",
  pxd073311_total_path
)

if (!length(quant_parts)) stop("No whole-proteome quantitative data extracted")

quant_features <- bind_rows(quant_parts) |>
  group_by(
    PXD, SampleGroup, QuantSample, CanonicalFeature,
    TargetAccession, TargetGene, Measurement
  ) |>
  summarise(
    Signal = sum(Signal, na.rm = TRUE),
    SourceFile = paste(sort(unique(SourceFile)), collapse = ";"),
    .groups = "drop"
  ) |>
  group_by(PXD, SampleGroup, QuantSample) |>
  mutate(
    Log2Signal = log2(Signal + 1),
    WholeProteomePercentile = ifelse(
      n() <= 1,
      100,
      100 * (rank(Log2Signal, ties.method = "average") - 1) / (n() - 1)
    )
  ) |>
  ungroup()

target_quant <- quant_features |>
  filter(
    !is.na(TargetAccession),
    TargetAccession %in% target_accessions
  ) |>
  group_by(
    PXD, SampleGroup, QuantSample,
    TargetAccession, TargetGene
  ) |>
  summarise(
    Signal = sum(Signal),
    Log2Signal = log2(Signal + 1),
    WholeProteomePercentile = max(WholeProteomePercentile),
    Measurement = paste(sort(unique(Measurement)), collapse = ";"),
    SourceFile = paste(sort(unique(SourceFile)), collapse = ";"),
    .groups = "drop"
  )

sample_registry <- quant_features |>
  distinct(PXD, SampleGroup, QuantSample)
target_sample_grid <- sample_registry |>
  crossing(
    regulators |>
      distinct(GeneSymbol, BaseAccession) |>
      rename(RegulatorBaseAccession = BaseAccession)
  ) |>
  left_join(
    target_quant,
    by = c(
      "PXD", "SampleGroup", "QuantSample",
      "RegulatorBaseAccession" = "TargetAccession"
    )
  ) |>
  mutate(
    GeneSymbol = unname(accession_to_gene[RegulatorBaseAccession]),
    IdentityMatchMode = "UniProt_BaseAccession_only",
    Detected = !is.na(Signal) & Signal > 0,
    Signal = replace_na(Signal, 0),
    Log2Signal = replace_na(Log2Signal, 0),
    WholeProteomePercentile = replace_na(WholeProteomePercentile, 0)
  )

group_summary <- target_sample_grid |>
  group_by(PXD, SampleGroup, RegulatorBaseAccession, GeneSymbol) |>
  summarise(
    WholeProteomeRelativePercentile = median(WholeProteomePercentile),
    DetectedSampleCount = sum(Detected),
    QuantSampleCount = n(),
    DetectedSampleFraction = mean(Detected),
    MedianLog2SignalDetected = ifelse(
      any(Detected),
      median(Log2Signal[Detected]),
      NA_real_
    ),
    Measurement = paste(sort(unique(na.omit(Measurement))), collapse = ";"),
    SourceFile = paste(sort(unique(na.omit(SourceFile))), collapse = ";"),
    .groups = "drop"
  )

quant_summary <- target_sample_grid |>
  group_by(PXD, SampleGroup) |>
  summarise(
    WholeProteomeSampleCount = n_distinct(QuantSample),
    WholeProteomeMappedRegulatorCount = n_distinct(
      RegulatorBaseAccession[Detected]
    ),
    .groups = "drop"
  )
whole_feature_summary <- quant_features |>
  group_by(PXD, SampleGroup) |>
  summarise(
    WholeProteomeFeatureCount = n_distinct(CanonicalFeature),
    .groups = "drop"
  )
quant_summary <- quant_summary |>
  left_join(whole_feature_summary, by = c("PXD", "SampleGroup"))

quant_audit <- quant_audit |>
  left_join(quant_summary, by = c("PXD", "SampleGroup")) |>
  mutate(
    WholeProteomeSampleCount = replace_na(WholeProteomeSampleCount, 0L),
    WholeProteomeFeatureCount = replace_na(
      WholeProteomeFeatureCount, 0L
    ),
    WholeProteomeMappedRegulatorCount = replace_na(
      WholeProteomeMappedRegulatorCount, 0L
    ),
    WholeProteomeQuantAvailable =
      WholeProteomeSampleCount > 0,
    WholeProteomeStatus = ifelse(
      WholeProteomeQuantAvailable, "available", "unavailable"
    )
  )

write.csv(
  quant_audit,
  file.path(
    table_dir,
    "kla_regulator_whole_proteome_intensity_availability_audit.csv"
  ),
  row.names = FALSE,
  na = ""
)
write.csv(
  quant_features,
  file.path(
    table_dir,
    "kla_regulator_whole_proteome_sample_level_long.csv"
  ),
  row.names = FALSE,
  na = ""
)
write.csv(
  target_sample_grid,
  file.path(
    table_dir,
    "kla_regulator_whole_proteome_regulator_sample_level_long.csv"
  ),
  row.names = FALSE,
  na = ""
)
write.csv(
  group_summary,
  file.path(
    table_dir,
    "kla_regulator_whole_proteome_normalized_long.csv"
  ),
  row.names = FALSE,
  na = ""
)

available_sample_groups <- paste(
  quant_audit$PXD[quant_audit$WholeProteomeQuantAvailable],
  quant_audit$SampleGroup[quant_audit$WholeProteomeQuantAvailable],
  sep = "__"
)
plot_data <- sample_catalog |>
  select(PXD, SampleGroup, SampleGroupID, RowLabel, RowOrder) |>
  crossing(
    regulators |>
      select(Role, GeneSymbol, BaseAccession, RoleEntryOrder) |>
      rename(RegulatorBaseAccession = BaseAccession)
  ) |>
  left_join(
    group_summary,
    by = c("PXD", "SampleGroup", "GeneSymbol", "RegulatorBaseAccession")
  ) |>
  mutate(
    QuantificationAvailable = paste(PXD, SampleGroup, sep = "__") %in%
      available_sample_groups,
    WholeProteomeRelativePercentile = ifelse(
      QuantificationAvailable,
      replace_na(WholeProteomeRelativePercentile, 0),
      NA_real_
    ),
    Role = factor(Role, levels = c("Writer", "Eraser", "Reader")),
    GeneSymbol = factor(
      GeneSymbol,
      levels = unique(regulators$GeneSymbol[order(regulators$RoleEntryOrder)])
    ),
    RowLabel = factor(RowLabel, levels = rev(sample_catalog$RowLabel))
  )

plot_font <- "Arial Unicode MS"
main_plot <- ggplot(
  plot_data,
  aes(GeneSymbol, RowLabel, fill = WholeProteomeRelativePercentile)
) +
  geom_tile(color = "white", linewidth = 0.22) +
  geom_text(
    data = plot_data |> filter(!QuantificationAvailable),
    aes(label = "?"),
    color = "#4D4D4D",
    size = 2.3,
    fontface = "bold",
    family = plot_font
  ) +
  facet_grid(. ~ Role, scales = "free_x", space = "free_x") +
  scale_fill_gradientn(
    colours = c("#FFFFFF", "#FFF3E0", "#FDBB84", "#FC8D59", "#B2182B"),
    values = scales::rescale(c(0, 20, 50, 80, 100)),
    limits = c(0, 100),
    na.value = "#D9D9D9",
    name = "全蛋白组信号\n百分位"
  ) +
  labs(
    title = "乳酸化调控蛋白在全蛋白组中的相对信号",
    subtitle = paste(
      "信号来自对应的普通全蛋白定量文件（同研究或独立正常参照），不使用Kla富集信号；",
      "颜色由白色向暖色递增；?表示没有可对应的全蛋白强度。"
    ),
    x = NULL,
    y = NULL,
    caption = paste0(
      "身份按人源UniProt reviewed BaseAccession匹配；",
      "全蛋白组百分位只在各自样本内计算，不是Log FC。"
    )
  ) +
  theme_minimal(base_size = 8.5, base_family = plot_font) +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold", size = 9),
    strip.background = element_rect(fill = "#F2F2F2", color = NA),
    axis.text.x = element_text(angle = 55, hjust = 1, vjust = 1, size = 7),
    axis.text.y = element_text(size = 7.2),
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 8.2, lineheight = 1.12),
    plot.caption = element_text(size = 7.5, hjust = 0, color = "#5C626A"),
    legend.position = "bottom",
    legend.key.width = grid::unit(35, "mm"),
    plot.margin = margin(8, 10, 8, 8)
  )

ggsave(
  file.path(
    figure_dir,
    "kla_regulator_whole_proteome_relative_intensity_heatmap.png"
  ),
  main_plot,
  width = 15.5,
  height = 11.5,
  dpi = 320,
  bg = "white"
)
ggsave(
  file.path(
    figure_dir,
    "kla_regulator_whole_proteome_relative_intensity_heatmap.pdf"
  ),
  main_plot,
  width = 15.5,
  height = 11.5,
  device = cairo_pdf,
  bg = "white"
)

message(
  "Whole-proteome regulator heatmap: ",
  sum(quant_audit$WholeProteomeQuantAvailable),
  "/",
  nrow(quant_audit),
  " sample groups with usable total-proteome intensity."
)
})
