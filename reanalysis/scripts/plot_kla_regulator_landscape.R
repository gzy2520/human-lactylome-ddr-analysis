#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(readxl)
  library(stringr)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) normalizePath(args[[1]]) else normalizePath(".")
table_dir <- file.path(project_root, "reanalysis", "results", "tables")
figure_dir <- file.path(project_root, "reanalysis", "results", "figures")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

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
scope_exclusion_path <- file.path(
  project_root, "reanalysis", "config",
  "final_sample_group_exclusions.csv"
)

required_files <- c(
  regulator_path, mapping_path, pairing_path, primary_path,
  scope_exclusion_path
)
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
    Role %in% c("Writer", "Eraser", "Writer-Eraser", "Reader"),
    !is.na(GeneSymbol),
    nzchar(GeneSymbol)
  ) |>
  distinct(Role, GeneSymbol, .keep_all = TRUE)
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
scope_exclusions <- read.csv(
  scope_exclusion_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
) |>
  select(PXD, SampleGroup) |>
  distinct()
sample_catalog <- pairing |>
  filter(当前Kla证据可分析) |>
  transmute(
    PXD = 乳酸化PXD,
    SampleGroup = 样本组,
    BiologicalMaterial = 材料类型,
    LactylomeEvidenceFile = 乳酸化证据文件,
    LactylomeProteinCount = 乳酸化蛋白数,
    RowOrder = row_number(),
    SampleGroupID = paste(PXD, SampleGroup, sep = "__")
  ) |>
  anti_join(scope_exclusions, by = c("PXD", "SampleGroup"))
if (nrow(sample_catalog) != 33) {
  stop("Expected 33 final sample groups, found ", nrow(sample_catalog))
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

base_accession <- function(values) {
  values <- sub("^.*\\|([^|]+)\\|.*$", "\\1", as.character(values))
  sub("-[0-9]+$", "", values)
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
  semi_join(sample_catalog, by = c("PXD", "SampleGroup")) |>
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
  `Writer-Eraser` = "#C83E4D",
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
    make_role_panel("Writer-Eraser", FALSE) +
    make_role_panel("Reader", FALSE) +
    plot_layout(
      widths = c(
        sum(regulators$Role == "Writer"),
        sum(regulators$Role == "Eraser"),
        sum(regulators$Role == "Writer-Eraser"),
        sum(regulators$Role == "Reader")
      ),
      guides = "collect"
    )
) +
  plot_annotation(
    title = "乳酸化调控因子在33个正式细胞系与组织样本组中的Kla证据图谱",
    subtitle = paste0(
      "数字为可定位的唯一Kla位点数；P表示蛋白/修饰肽层面检出但无法得到精确位点数；",
      "无法可靠拆分的样本组不进入本图。颜色不代表表达量或Log2 FC。"
    ),
    caption = paste0(
      "数据范围：最终正式分析固定为33个样本组；已删除的7个样本组不进入活动数据、图表或后续分析。",
      "主分析定位阈值为>0；保留老师表中的多重角色记录，Reader分区置于最后。"
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
