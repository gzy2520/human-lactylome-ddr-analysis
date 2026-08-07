#!/usr/bin/env Rscript
# =============================================================================
# acquire_data.R — 数据获取类脚本合并版（kla 全 R 重构 Task 5）
#
# 合并来源脚本（原脚本不做修改，归档后删除）：
#   计算段（默认全部执行）：
#     build_ensembl_uniprot_mapping.R      -> 段 ensembl_map
#     build_human_lactylome_inventory.R    -> 段 inventory
#     build_lactylome_acquisition_manifests.R + build_healthy_special_reference_manifest.R
#                                          -> 段 manifests
#     build_lactylome_reference_pairing.R  -> 段 pairing
#     summarize_acquired_lactylome_data.R  -> 段 summarize
#   下载/解压/探测段（默认不执行，仅显式 --stage download）：
#     probe_lactylome_pair_files.R
#     download_lactylome_pair_files.R
#     extract_lactylome_pair_archives.R
#     register_additional_lactylome_pair_files.R
#     download_healthy_tissue_references.R
#
# 用法：
#   Rscript reanalysis/scripts/acquire_data.R <project_root> [--stage <name>]
#     <project_root>  项目根目录（默认 "."）
#     --stage name    只运行指定段；name 取值：
#                     pairing | inventory | manifests | ensembl_map | summarize
#                     或 download（下载段，默认跳过，避免重复联网下载）。
#
# 与 lib/ 的重复定义处理（详见 reports/task-5-report.md）：
#   - 已改用 lib 版本并删除本地定义：
#       * build_ensembl_uniprot_mapping.R 中 read.delim + 行级过滤
#         -> lib io_utils 的 read_delimited + valid_maxquant_rows（语义逐项等价）
#       * build_lactylome_acquisition_manifests.R 中 relative_to_project
#         -> lib io_utils 的 relative_path（本文件顶部加一层同签名包装，
#            避免与段内局部变量 relative_path 重名遮蔽）
#       * sha256_file 在 4 个来源脚本中重复定义（逻辑完全相同），
#         本文件顶部合并为一处（lib 中无此函数，故不属 lib 重复）
#   - 有意保留段内本地定义（与 lib 版本语义不同，替换会改变输出，
#     破坏基线字节一致性；lib 版本会 trimws/去 NX_ 前缀并对 UniProt
#     做 is_uniprot 过滤，而原脚本统计全部 token 如 CON__/ENSEMBL:）：
#       * build_lactylome_reference_pairing.R 的 base_accession / split_accessions
#         （实测替换会使 PXD063047 计数 386->348、PXD055230 参考计数
#         10178->10131 等，见报告）
#       * download_healthy_tissue_references.R 的 base_accession（同原因）
# =============================================================================

suppressPackageStartupMessages({
  library(biomaRt)
  library(dplyr)
  library(jsonlite)
  library(stringr)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
project_root <- "."
if (length(args) && args[[1]] != "--stage") {
  project_root <- args[[1]]
}
stage <- "all"
if ("--stage" %in% args) {
  idx <- which(args == "--stage")
  last_idx <- idx[[length(idx)]]
  if (last_idx >= length(args)) {
    stop("--stage requires a stage name")
  }
  stage <- args[[last_idx + 1]]
}
project_root <- normalizePath(project_root)
valid_stages <- c("pairing", "inventory", "manifests", "ensembl_map", "download", "summarize")
if (!stage %in% c("all", valid_stages)) {
  stop(
    "Unknown --stage '", stage, "'; valid stages: ",
    paste(valid_stages, collapse = ", ")
  )
}

# 共享 lib（加载顺序：accession_utils -> io_utils -> extractors）
lib_dir <- file.path(project_root, "reanalysis", "scripts", "lib")
source(file.path(lib_dir, "accession_utils.R"))
source(file.path(lib_dir, "io_utils.R"))
source(file.path(lib_dir, "extractors.R"))

# run_stage 模式：stage == "all"（缺省）只跑计算段；
# download 段用独立 gate（见文末），仅显式 --stage download 时执行。
run_stage <- function(name, code) {
  if (stage == "all" || stage == name) {
    message("[stage] ", name)
    code
  }
}

# sha256_file：来自 build_healthy_special_reference_manifest.R /
#   download_healthy_tissue_references.R / download_lactylome_pair_files.R /
#   register_additional_lactylome_pair_files.R（4 处定义完全相同，合并为一处；
#   lib 中无此函数）。
sha256_file <- function(path) {
  output <- system2(
    "shasum",
    c("-a", "256", shQuote(path)),
    stdout = TRUE
  )
  sub("[[:space:]].*$", "", output[[1]])
}

# relative_to_project：原 build_lactylome_acquisition_manifests.R 内定义
#   （normalizePath + str_remove(fixed)）。lib io_utils 的 relative_path
#   语义等价（现有文件上逐字节同结果），此处保留原签名作为薄包装。
#   定义在顶层（而非段内），避免段内局部变量 relative_path 遮蔽后包装失效。
relative_to_project <- function(path) {
  relative_path(path, project_root)
}

# =============================================================================
# 段 ensembl_map — 来源：build_ensembl_uniprot_mapping.R
# =============================================================================
run_build_ensembl_uniprot_mapping <- function() {
  output_path <- file.path(
    project_root,
    "reanalysis", "config", "ensembl_protein_to_uniprot_biomart.tsv"
  )
  unmapped_path <- file.path(
    project_root,
    "reanalysis", "config", "ensembl_protein_unmapped_biomart.tsv"
  )

  files <- list.files(
    file.path(
      project_root,
      "data", "PXD010154", "search_results", "extracted_healthy_tissues"
    ),
    pattern = "proteinGroups\\.txt$",
    recursive = TRUE,
    full.names = TRUE
  )

  ensembl_ids <- character()
  for (path in files) {
    # 原脚本 read.delim(..., quote="", comment.char="") + 行级 MaxQuant 过滤；
    # 改用 lib 的 read_delimited + valid_maxquant_rows（语义逐项等价）。
    data <- read_delimited(path)
    keep <- valid_maxquant_rows(data)
    values <- unlist(strsplit(
      as.character(data$`Majority protein IDs`[keep]),
      "[;,]"
    ))
    values <- sub("_RNA$", "", trimws(values))
    values <- sub("\\.[0-9]+$", "", values)
    ensembl_ids <- c(
      ensembl_ids,
      values[grepl("^ENSP[0-9]+$", values)]
    )
  }
  ensembl_ids <- sort(unique(ensembl_ids))
  if (!length(ensembl_ids)) stop("No Ensembl protein IDs found")

  connect_mart <- function() {
    mirrors <- c("www", "useast", "asia")
    for (mirror in mirrors) {
      mart <- tryCatch(
        useEnsembl(
          biomart = "genes",
          dataset = "hsapiens_gene_ensembl",
          mirror = mirror
        ),
        error = function(e) NULL
      )
      if (!is.null(mart)) return(mart)
    }
    stop("Unable to connect to an Ensembl BioMart mirror")
  }

  mart <- connect_mart()
  chunks <- split(ensembl_ids, ceiling(seq_along(ensembl_ids) / 2000))
  mapped_parts <- vector("list", length(chunks))

  for (i in seq_along(chunks)) {
    values <- chunks[[i]]
    result <- NULL
    for (attempt in seq_len(3)) {
      result <- tryCatch(
        getBM(
          attributes = c(
            "ensembl_peptide_id",
            "uniprotswissprot",
            "uniprotsptrembl"
          ),
          filters = "ensembl_peptide_id",
          values = values,
          mart = mart
        ),
        error = function(e) NULL
      )
      if (!is.null(result)) break
      Sys.sleep(attempt * 2)
      mart <- connect_mart()
    }
    if (is.null(result)) {
      stop("BioMart mapping failed for chunk ", i, " of ", length(chunks))
    }
    mapped_parts[[i]] <- result
    message(
      "Mapped BioMart chunk ", i, "/", length(chunks),
      " (", length(values), " Ensembl protein IDs)"
    )
  }

  raw_mapping <- bind_rows(mapped_parts)
  mapping <- bind_rows(
    raw_mapping |>
      transmute(
        EnsemblProteinID = ensembl_peptide_id,
        BaseAccession = uniprotswissprot,
        UniProtSource = "Swiss-Prot"
      ),
    raw_mapping |>
      transmute(
        EnsemblProteinID = ensembl_peptide_id,
        BaseAccession = uniprotsptrembl,
        UniProtSource = "TrEMBL"
      )
  ) |>
    mutate(
      BaseAccession = sub("-[0-9]+$", "", trimws(BaseAccession)),
      MappingSource = "Ensembl BioMart hsapiens_gene_ensembl",
      MappingDate = format(Sys.Date(), "%Y-%m-%d")
    ) |>
    filter(
      nzchar(EnsemblProteinID),
      nzchar(BaseAccession)
    ) |>
    distinct() |>
    arrange(EnsemblProteinID, UniProtSource, BaseAccession)

  unmapped <- data.frame(
    EnsemblProteinID = setdiff(ensembl_ids, unique(mapping$EnsemblProteinID)),
    MappingSource = "Ensembl BioMart hsapiens_gene_ensembl",
    MappingDate = format(Sys.Date(), "%Y-%m-%d"),
    stringsAsFactors = FALSE
  )

  write.table(
    mapping,
    output_path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  write.table(
    unmapped,
    unmapped_path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  message(
    "BioMart mapping complete: ",
    length(unique(mapping$EnsemblProteinID)), "/", length(ensembl_ids),
    " Ensembl protein IDs mapped; ", nrow(unmapped), " unmapped."
  )
}

# =============================================================================
# 段 inventory — 来源：build_human_lactylome_inventory.R
# =============================================================================
run_build_human_lactylome_inventory <- function() {
  config_path <- file.path(project_root, "reanalysis", "config", "lactylome_manual_curation.csv")
  metadata_dir <- file.path(project_root, "data", "metadata", "lactylome_discovery")
  table_dir <- file.path(project_root, "reanalysis", "results", "tables")
  report_dir <- file.path(project_root, "reanalysis", "reports")
  dir.create(metadata_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)

  search_terms <- c(
    "lactylation",
    "lactylome",
    "lactylproteom",
    "lactoylation",
    "lactyllysine"
  )
  api_root <- "https://proteomecentral.proteomexchange.org/api/proxi/v0.1/datasets"

  fetch_json <- function(url, attempts = 4L) {
    last_error <- NULL
    for (attempt in seq_len(attempts)) {
      result <- tryCatch(
        fromJSON(url, simplifyVector = FALSE),
        error = function(error) {
          last_error <<- error
          NULL
        }
      )
      if (!is.null(result)) {
        return(result)
      }
      Sys.sleep(attempt)
    }
    stop("Unable to fetch ", url, ": ", conditionMessage(last_error))
  }

  collapse_values <- function(values) {
    values <- unique(na.omit(as.character(values)))
    values <- values[nzchar(values)]
    paste(values, collapse = "; ")
  }

  term_values <- function(groups, wanted_name = NULL) {
    if (is.null(groups) || length(groups) == 0) {
      return(character())
    }
    terms <- unlist(lapply(groups, function(group) group$terms), recursive = FALSE)
    if (!is.null(wanted_name)) {
      terms <- Filter(function(item) identical(item$name, wanted_name), terms)
    }
    vapply(terms, function(item) item$value %||% item$name %||% "", character(1))
  }

  `%||%` <- function(left, right) {
    if (is.null(left) || length(left) == 0 || (length(left) == 1 && is.na(left))) {
      right
    } else {
      left
    }
  }

  strip_html <- function(value) {
    value <- str_replace_all(value %||% "", "<[^>]+>", " ")
    str_squish(str_replace_all(value, "&amp;", "&"))
  }

  search_rows <- list()
  raw_search_responses <- list()
  for (term in search_terms) {
    url <- paste0(
      api_root,
      "?src=PCUI&resultType=resultset&pageSize=100&search=",
      URLencode(term, reserved = TRUE),
      "&species=Homo%20sapiens"
    )
    payload <- fetch_json(url)
    raw_search_responses[[term]] <- payload
    rows <- payload$datasets
    if (length(rows) == 0) {
      next
    }
    search_rows[[term]] <- bind_rows(lapply(rows, function(row) {
      length(row) <- max(length(row), 11L)
      tibble(
        PXD = as.character(row[[1]] %||% ""),
        SearchTitle = as.character(row[[2]] %||% ""),
        Repository = as.character(row[[3]] %||% ""),
        SpeciesSearch = as.character(row[[4]] %||% ""),
        SDRFStatus = as.character(row[[5]] %||% ""),
        FilesRawTotal = as.character(row[[6]] %||% ""),
        InstrumentSearch = as.character(row[[7]] %||% ""),
        PublicationSearch = strip_html(as.character(row[[8]] %||% "")),
        LabHeadSearch = as.character(row[[9]] %||% ""),
        AnnounceDateSearch = as.character(row[[10]] %||% ""),
        KeywordsSearch = as.character(row[[11]] %||% ""),
        MatchingSearchTerm = term
      )
    }))
  }

  search_frame <- bind_rows(search_rows) |>
    filter(str_detect(PXD, "^PXD[0-9]{6}$")) |>
    group_by(PXD) |>
    summarise(
      across(
        c(
          SearchTitle, Repository, SpeciesSearch, SDRFStatus, FilesRawTotal,
          InstrumentSearch, PublicationSearch, LabHeadSearch, AnnounceDateSearch,
          KeywordsSearch
        ),
        ~ first(.x[nzchar(.x)], default = "")
      ),
      MatchingSearchTerms = collapse_values(MatchingSearchTerm),
      .groups = "drop"
    )

  # "Kla" is too short for the public text-search endpoint. Add the known source
  # dataset whose repository metadata uses Kla rather than the word lactylation.
  manual_discovery <- tibble(
    PXD = c("PXD014870"),
    SearchTitle = "Histone Kla is regulated by cellular glycolysis",
    Repository = "PRIDE",
    SpeciesSearch = "Homo sapiens",
    SDRFStatus = "",
    FilesRawTotal = "",
    InstrumentSearch = "",
    PublicationSearch = "",
    LabHeadSearch = "",
    AnnounceDateSearch = "2019-08-01",
    KeywordsSearch = "Human MCF-7 cells SILAC",
    MatchingSearchTerms = "manual_Kla_synonym"
  )
  search_frame <- bind_rows(search_frame, manual_discovery) |>
    distinct(PXD, .keep_all = TRUE) |>
    arrange(PXD)

  fetch_detail <- function(pxd) {
    url <- paste0(api_root, "/", pxd, "?src=PCUI&resultType=full")
    tryCatch(
      fetch_json(url),
      error = function(error) list(
        identifiers = list(),
        title = "",
        description = "",
        fetch_error = conditionMessage(error)
      )
    )
  }

  worker_count <- max(1L, min(6L, parallel::detectCores(logical = FALSE)))
  detail_list <- parallel::mclapply(
    search_frame$PXD,
    fetch_detail,
    mc.cores = worker_count
  )
  names(detail_list) <- search_frame$PXD

  detail_rows <- list()
  file_rows <- list()
  for (pxd in names(detail_list)) {
    item <- detail_list[[pxd]]
    files <- item$datasetFiles %||% list()
    file_names <- vapply(files, function(file) file$name %||% "", character(1))
    file_urls <- vapply(files, function(file) file$value %||% "", character(1))
    publications <- item$publications %||% list()
    doi <- term_values(publications, "Digital Object Identifier (DOI)")
    pubmed <- term_values(publications, "PubMed identifier")
    references <- term_values(publications, "Reference")
    contacts <- item$contacts %||% list()
    contact_names <- term_values(contacts, "contact name")
    affiliations <- term_values(contacts, "contact affiliation")
    species <- unlist(lapply(item$species %||% list(), function(group) {
      term_values(list(group), "taxonomy: scientific name")
    }))
    instruments <- vapply(
      item$instruments %||% list(),
      function(value) value$name %||% "",
      character(1)
    )
    keywords <- vapply(
      item$keywords %||% list(),
      function(value) value$value %||% value$name %||% "",
      character(1)
    )
    links <- item$fullDatasetLinks %||% list()
    dataset_links <- vapply(links, function(link) link$value %||% "", character(1))
    detail_rows[[pxd]] <- tibble(
      PXD = pxd,
      Title = item$title %||% "",
      Description = item$description %||% "",
      RepositoryDetail = item$datasetSummary$hostingRepository %||% "",
      AnnounceDate = item$datasetSummary$announceDate %||% "",
      Species = collapse_values(species),
      Instruments = collapse_values(instruments),
      Keywords = collapse_values(keywords),
      DOI = collapse_values(doi),
      PubMed = collapse_values(pubmed),
      PublicationCitation = collapse_values(references),
      Contact = collapse_values(contact_names),
      Affiliation = collapse_values(affiliations),
      RawFileCount = sum(str_detect(file_names, regex("raw file", ignore_case = TRUE))),
      SearchResultFileCount = sum(str_detect(file_names, regex("search engine output|result file", ignore_case = TRUE))),
      OtherFileCount = sum(!str_detect(file_names, regex("raw file|search engine output|result file", ignore_case = TRUE))),
      DatasetLinks = collapse_values(dataset_links),
      MetadataFetchError = item$fetch_error %||% ""
    )
    if (length(files) > 0) {
      file_rows[[pxd]] <- bind_rows(lapply(seq_along(files), function(index) {
        url <- file_urls[[index]]
        tibble(
          PXD = pxd,
          FileCategory = file_names[[index]],
          FileName = basename(sub("\\?.*$", "", url)),
          FileURL = url
        )
      }))
    }
  }

  details <- bind_rows(detail_rows)
  file_manifest <- bind_rows(file_rows)
  curation <- read.csv(config_path, stringsAsFactors = FALSE, check.names = FALSE)
  file_exception_path <- file.path(
    project_root,
    "reanalysis",
    "config",
    "lactylome_repository_file_exceptions.csv"
  )
  file_exceptions <- if (file.exists(file_exception_path)) {
    read.csv(file_exception_path, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    tibble(
      PXD = character(),
      FileURL = character(),
      LocalStatusOverride = character(),
      LocalMatchesOverride = character(),
      RepositoryFileNote = character()
    )
  }
  normal_candidate_categories <- c(
    "physiological_normal_tissue",
    "includes_healthy_control_tissue",
    "normal_human_biospecimen",
    "normal_cell_control",
    "contains_normal_cell_model",
    "contains_normal_like_cell_line",
    "normal_like_immortalized_cell_line",
    "adjacent_normal_not_healthy",
    "benign_disease_control_not_healthy",
    "contains_normal_or_control_material",
    "possible_physiological_normal_biospecimen"
  )

  classify_scope <- function(title, description, keywords) {
    text <- str_to_lower(paste(title, description, keywords))
    global_pattern <- regex(
      "global|lactylome|lactylproteom|proteome-wide|systematic profil|systematic identif|quantitative .*lactyl|protein lactylation modification and proteomics|lactylation-modified proteomics",
      ignore_case = TRUE
    )
    targeted_pattern <- regex(
      "ip-ms|interact|parallel reaction|\\bprm\\b|affinity purification|photoaffinity|chemical reporter|reader|genetic code expansion|site-specific|secondary mass spectra|single protein|cross-link",
      ignore_case = TRUE
    )
    if (str_detect(text, global_pattern)) {
      return("global_lactylome_candidate")
    }
    if (str_detect(text, targeted_pattern)) {
      return("targeted_or_mechanistic_ms")
    }
    "manual_scope_review_needed"
  }

  classify_normal_candidate <- function(title, description, keywords) {
    text <- str_to_lower(paste(title, description, keywords))
    case_when(
      str_detect(text, "normal physiological|healthy donor|normal pregnan|normal human|human sperm") ~
        "possible_physiological_normal_biospecimen",
      str_detect(text, "adjacent normal") ~ "adjacent_tissue_not_healthy_normal",
      str_detect(text, "benign prostatic hyperplasia|\\bbph\\b") ~ "benign_disease_control_not_healthy_normal",
      str_detect(text, "normal tissue|control placent|neural stem cell|mock control|uninfected") ~
        "contains_normal_or_control_material",
      str_detect(text, "cell|hela|hek|hct|mcf|a549|h1299|hcc|cancer|carcinoma|glioma|tumor") ~
        "cell_line_or_cancer_material",
      TRUE ~ "manual_normality_review_needed"
    )
  }

  inventory <- search_frame |>
    left_join(details, by = "PXD") |>
    mutate(
      Repository = coalesce(na_if(RepositoryDetail, ""), Repository),
      Title = coalesce(na_if(Title, ""), SearchTitle),
      AnnounceDate = coalesce(na_if(AnnounceDate, ""), AnnounceDateSearch),
      RepositoryAnnounceYear = suppressWarnings(as.integer(str_sub(AnnounceDate, 1, 4))),
      Instruments = coalesce(na_if(Instruments, ""), InstrumentSearch),
      Keywords = coalesce(na_if(Keywords, ""), KeywordsSearch),
      DOI = coalesce(
        na_if(DOI, ""),
        str_extract(PublicationCitation, "10\\.[0-9]{4,9}/[-._;()/:A-Za-z0-9]+")
      ),
      PublicationStatus = case_when(
        nzchar(DOI) | nzchar(PubMed) | nzchar(PublicationCitation) ~ "published_or_linked_publication",
        str_detect(PublicationSearch, regex("publication pending", ignore_case = TRUE)) ~ "publication_pending",
        TRUE ~ "no_linked_publication_found"
      ),
      PublicationYear = suppressWarnings(as.integer(str_extract(
        PublicationCitation,
        "(?<![0-9])(19|20)[0-9]{2}(?![0-9])"
      ))),
      ScopeAutomated = mapply(classify_scope, Title, Description, Keywords),
      NormalityAutomated = mapply(classify_normal_candidate, Title, Description, Keywords)
    ) |>
    left_join(curation, by = "PXD") |>
    mutate(
      MassSpecScope = coalesce(na_if(ScopeOverride, ""), ScopeAutomated),
      NormalTissueCategory = coalesce(na_if(NormalTissueCategory, ""), NormalityAutomated),
      PriorityTier = coalesce(na_if(PriorityTier, ""), "review"),
      ManualCurationStatus = if_else(
        !is.na(ScopeOverride) & nzchar(ScopeOverride),
        "manually_curated",
        "metadata_screen_only"
      ),
      HasProcessedOrSearchResult = SearchResultFileCount > 0,
      RepositoryDataStatus = case_when(
        PXD == "PXD054919" ~ "processed_archive_zero_filled_but_article_supplement_usable",
        SearchResultFileCount > 0 & RawFileCount > 0 ~ "raw_and_processed_files_listed",
        SearchResultFileCount > 0 ~ "processed_files_listed",
        RawFileCount > 0 ~ "raw_files_only",
        TRUE ~ "metadata_only_or_file_list_unavailable"
      ),
      DatasetURL = paste0("https://proteomecentral.proteomexchange.org/?pxid=", PXD)
    ) |>
    select(
      PXD, Title, Repository, PublicationYear, RepositoryAnnounceYear,
      AnnounceDate, Species,
      MassSpecScope, ManualCurationStatus, SampleMaterial, TissueOrCellType,
      HealthContext, NormalTissueCategory, PriorityTier, StudyFamily,
      RawFileCount, SearchResultFileCount, OtherFileCount,
      HasProcessedOrSearchResult, RepositoryDataStatus, PublicationStatus,
      DOI, PubMed, PublicationCitation, Instruments, Keywords, Description,
      Contact, Affiliation, MatchingSearchTerms, SDRFStatus, FilesRawTotal,
      CurationNote, DatasetURL, DatasetLinks, MetadataFetchError
    ) |>
    arrange(
      factor(PriorityTier, levels = c("1", "2", "review", "hold")),
      desc(coalesce(PublicationYear, RepositoryAnnounceYear)),
      PXD
    )

  find_local_matches <- function(pxd, file_name, file_url) {
    pxd_dir <- file.path(project_root, "data", pxd)
    if (!dir.exists(pxd_dir)) {
      return(character())
    }
    local_files <- list.files(
      pxd_dir,
      recursive = TRUE,
      full.names = TRUE,
      all.files = TRUE,
      no.. = TRUE
    )
    url_parent <- basename(dirname(sub("\\?.*$", "", file_url)))
    extension <- tools::file_ext(file_name)
    stem <- tools::file_path_sans_ext(file_name)
    renamed_download <- if (nzchar(extension)) {
      paste0(stem, "_", url_parent, ".", extension)
    } else {
      paste0(stem, "_", url_parent)
    }
    local_files[basename(local_files) %in% c(file_name, renamed_download)]
  }

  file_manifest <- file_manifest |>
    mutate(
      LocalMatches = mapply(
        function(pxd, name, url) {
          collapse_values(find_local_matches(pxd, name, url))
        },
        PXD,
        FileName,
        FileURL,
        USE.NAMES = FALSE
      ),
      LocalStatus = if_else(nzchar(LocalMatches), "present_locally", "not_downloaded"),
      DownloadPriority = case_when(
        str_detect(FileCategory, regex("search engine output|result file", ignore_case = TRUE)) ~
          "processed_result_first",
        str_detect(FileCategory, regex("raw file", ignore_case = TRUE)) ~
          "raw_defer_unless_research_required",
        TRUE ~ "metadata_or_other"
      )
    ) |>
    left_join(file_exceptions, by = c("PXD", "FileURL")) |>
    mutate(
      LocalMatches = coalesce(na_if(LocalMatchesOverride, ""), LocalMatches),
      LocalStatus = coalesce(na_if(LocalStatusOverride, ""), LocalStatus),
      RepositoryFileNote = coalesce(RepositoryFileNote, "")
    ) |>
    select(-LocalStatusOverride, -LocalMatchesOverride) |>
    arrange(PXD, DownloadPriority, FileName)

  write.csv(
    inventory,
    file.path(table_dir, "human_lactylome_mass_spectrometry_inventory.csv"),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )
  write.csv(
    file_manifest,
    file.path(table_dir, "human_lactylome_repository_file_manifest.csv"),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )
  write.csv(
    inventory |>
      filter(NormalTissueCategory %in% normal_candidate_categories),
    file.path(table_dir, "human_lactylome_normal_material_candidates.csv"),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )
  write_json(
    list(
      generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
      search_terms = search_terms,
      searches = raw_search_responses,
      details = detail_list
    ),
    file.path(metadata_dir, "proteomexchange_lactylome_discovery_snapshot.json"),
    pretty = TRUE,
    auto_unbox = TRUE,
    null = "null"
  )

  summary_rows <- inventory |>
    count(MassSpecScope, NormalTissueCategory, name = "DatasetCount") |>
    arrange(MassSpecScope, desc(DatasetCount))
  write.csv(
    summary_rows,
    file.path(table_dir, "human_lactylome_inventory_summary.csv"),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )

  normal_priority <- inventory |>
    filter(PriorityTier == "1", NormalTissueCategory %in% normal_candidate_categories) |>
    select(PXD, Title, TissueOrCellType, NormalTissueCategory, RepositoryDataStatus, CurationNote)

  report <- c(
    "# 人源乳酸化质谱公开数据盘点",
    "",
    paste0("检索日期：", format(Sys.Date(), "%Y-%m-%d")),
    "",
    paste0("- ProteomeXchange 多关键词去重后数据集：", nrow(inventory)),
    paste0("- 全局乳酸化蛋白组候选：", sum(inventory$MassSpecScope %in% c("global_lactylome", "global_lactylome_candidate"))),
    paste0("- 已列出检索结果或处理文件：", sum(inventory$HasProcessedOrSearchResult, na.rm = TRUE)),
    paste0("- 人工复核完成：", sum(inventory$ManualCurationStatus == "manually_curated")),
    "",
    "## 优先核查的正常或对照材料",
    "",
    paste0(
      "- ",
      normal_priority$PXD,
      "：",
      normal_priority$TissueOrCellType,
      "；",
      normal_priority$NormalTissueCategory,
      "；",
      normal_priority$RepositoryDataStatus
    ),
    "",
    "## 解释边界",
    "",
    "- 正常生理组织、邻癌组织、良性病变、正常原代细胞和永生化细胞系分别标注，不能互相替代。",
    "- 总表是公开库检索与元数据盘点；`metadata_screen_only` 行仍需结合论文或 SDRF 复核样本。",
    "- 优先下载作者检索结果和补充位点表；大型 raw 文件仅在统一重检索时下载。",
    "- PXD054919 的仓库 Results.zip 与提交 SHA-1 一致但内容全为零，位点数据应使用论文补充表。",
    ""
  )
  writeLines(report, file.path(report_dir, "HUMAN_LACTYLOME_DATA_ACQUISITION.md"), useBytes = TRUE)

  message("Human lactylome inventory rows: ", nrow(inventory))
  message(
    "Global lactylome candidates: ",
    sum(inventory$MassSpecScope %in% c("global_lactylome", "global_lactylome_candidate"))
  )
  message("Repository file rows: ", nrow(file_manifest))
}

# =============================================================================
# 段 manifests — 来源：build_lactylome_acquisition_manifests.R
# =============================================================================
run_build_lactylome_acquisition_manifests <- function() {
  table_dir <- file.path(project_root, "reanalysis", "results", "tables")
  config_dir <- file.path(project_root, "reanalysis", "config")

  selected_pxd <- c(
    "PXD036307",
    "PXD054919",
    "PXD063047",
    "PXD064912",
    "PXD066054",
    "PXD075377"
  )

  inventory <- read.csv(
    file.path(table_dir, "human_lactylome_mass_spectrometry_inventory.csv"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  repository_files <- read.csv(
    file.path(table_dir, "human_lactylome_repository_file_manifest.csv"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  validation <- read.csv(
    file.path(config_dir, "lactylome_acquisition_validation.csv"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  # relative_to_project 改用 lib io_utils::relative_path 的薄包装（文件顶部定义）

  repository_alias <- function(file_name, file_url) {
    url_parent <- basename(dirname(sub("\\?.*$", "", file_url)))
    extension <- tools::file_ext(file_name)
    stem <- tools::file_path_sans_ext(file_name)
    if (nzchar(extension)) {
      paste0(stem, "_", url_parent, ".", extension)
    } else {
      paste0(stem, "_", url_parent)
    }
  }

  classify_role <- function(relative_path) {
    case_when(
      str_detect(relative_path, "/raw/") ~ "raw",
      str_detect(relative_path, "/search_results/") ~ "search_result",
      str_detect(relative_path, "/supplementary/") ~ "article_supplement",
      str_detect(relative_path, "/metadata/articles/") ~ "article_original",
      str_detect(relative_path, "/metadata/") ~ "metadata",
      TRUE ~ "other"
    )
  }

  source_files_for_pxd <- function(pxd) {
    pxd_dir <- file.path(project_root, "data", pxd)
    if (!dir.exists(pxd_dir)) {
      return(character())
    }
    files <- list.files(
      pxd_dir,
      recursive = TRUE,
      full.names = TRUE,
      all.files = TRUE,
      no.. = TRUE
    )
    files <- files[file.info(files)$isdir %in% FALSE]
    files <- files[!str_detect(files, "/search_results/extracted/")]
    files <- files[!str_detect(files, "\\.partial_download$")]
    files <- files[!basename(files) %in% c(
      "dataset_metadata.csv",
      "repository_file_manifest.csv",
      "download_manifest.csv",
      "extracted_file_inventory.csv"
    )]
    files
  }

  find_repository_record <- function(pxd, local_name) {
    candidates <- repository_files |>
      filter(PXD == pxd) |>
      mutate(LocalAlias = mapply(repository_alias, FileName, FileURL, USE.NAMES = FALSE)) |>
      filter(FileName == local_name | LocalAlias == local_name)
    if (nrow(candidates) == 0) {
      return(tibble(
        RepositoryFileCategory = "",
        SourceURL = "",
        DownloadPriority = ""
      ))
    }
    candidates |>
      slice(1) |>
      transmute(
        RepositoryFileCategory = FileCategory,
        SourceURL = FileURL,
        DownloadPriority
      )
  }

  all_downloads <- list()
  for (pxd in selected_pxd) {
    metadata_dir <- file.path(project_root, "data", pxd, "metadata")
    dir.create(metadata_dir, recursive = TRUE, showWarnings = FALSE)

    dataset_row <- inventory |>
      filter(PXD == pxd)
    write.csv(
      dataset_row,
      file.path(metadata_dir, "dataset_metadata.csv"),
      row.names = FALSE,
      fileEncoding = "UTF-8"
    )

    repository_subset <- repository_files |>
      filter(PXD == pxd)
    write.csv(
      repository_subset,
      file.path(metadata_dir, "repository_file_manifest.csv"),
      row.names = FALSE,
      fileEncoding = "UTF-8"
    )

    local_files <- source_files_for_pxd(pxd)
    local_rows <- lapply(local_files, function(path) {
      local_name <- basename(path)
      source_row <- find_repository_record(pxd, local_name)
      validation_row <- validation |>
        filter(PXD == pxd, FileName == local_name) |>
        slice(1)
      if (nrow(validation_row) == 0) {
        validation_row <- tibble(
          SHA256 = "",
          RepositorySHA1 = "",
          ValidationStatus = "downloaded_not_independently_validated",
          ValidationNote = ""
        )
      }
      relative_path <- relative_to_project(path)
      tibble(
        PXD = pxd,
        LocalRelativePath = relative_path,
        FileName = local_name,
        FileRole = classify_role(relative_path),
        SizeBytes = unname(file.info(path)$size),
        SHA256 = validation_row$SHA256,
        RepositorySHA1 = validation_row$RepositorySHA1,
        ValidationStatus = validation_row$ValidationStatus,
        ValidationNote = validation_row$ValidationNote,
        RepositoryFileCategory = source_row$RepositoryFileCategory,
        SourceURL = source_row$SourceURL,
        DownloadPriority = source_row$DownloadPriority
      )
    })
    local_manifest <- bind_rows(local_rows) |>
      arrange(FileRole, FileName)
    write.csv(
      local_manifest,
      file.path(metadata_dir, "download_manifest.csv"),
      row.names = FALSE,
      fileEncoding = "UTF-8"
    )
    all_downloads[[pxd]] <- local_manifest

    extracted_root <- file.path(project_root, "data", pxd, "search_results", "extracted")
    extracted_files <- if (dir.exists(extracted_root)) {
      list.files(extracted_root, recursive = TRUE, full.names = TRUE)
    } else {
      character()
    }
    extracted_files <- extracted_files[file.info(extracted_files)$isdir %in% FALSE]
    extracted_manifest <- tibble(
      PXD = pxd,
      LocalRelativePath = vapply(extracted_files, relative_to_project, character(1)),
      FileName = basename(extracted_files),
      SizeBytes = unname(file.info(extracted_files)$size),
      EvidenceFile = str_detect(
        basename(extracted_files),
        regex("La \\(K\\)Sites|PTMSiteReport|MS_identified_information|evidence|proteinGroups|Identification", ignore_case = TRUE)
      )
    ) |>
      arrange(desc(EvidenceFile), LocalRelativePath)
    write.csv(
      extracted_manifest,
      file.path(metadata_dir, "extracted_file_inventory.csv"),
      row.names = FALSE,
      fileEncoding = "UTF-8"
    )

    proxi_path <- file.path(metadata_dir, "proxi_detail.json")
    if (file.exists(proxi_path)) {
      parsed <- tryCatch(fromJSON(proxi_path, simplifyVector = FALSE), error = function(error) NULL)
      if (is.null(parsed)) {
        warning("Invalid PROXI JSON for ", pxd)
      }
    }
  }

  combined_downloads <- bind_rows(all_downloads) |>
    arrange(PXD, FileRole, FileName)
  write.csv(
    combined_downloads,
    file.path(table_dir, "priority_dataset_acquisition_manifest.csv"),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )

  message("Acquisition manifest rows: ", nrow(combined_downloads))
  message("Per-PXD metadata written for: ", paste(selected_pxd, collapse = ", "))
}

# =============================================================================
# 段 manifests（续）— 来源：build_healthy_special_reference_manifest.R
# =============================================================================
run_build_healthy_special_reference_manifest <- function() {
  config <- read.csv(
    file.path(
      project_root,
      "reanalysis/config/healthy_special_reference_catalog.csv"
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  output_path <- file.path(
    project_root,
    "reanalysis/results/tables/healthy_special_reference_acquisition_manifest.csv"
  )

  resolve_files <- function(locator) {
    locator <- sub(" \\[.*$", "", locator)
    pattern <- file.path(project_root, locator)
    if (grepl("\\*", pattern)) Sys.glob(pattern) else pattern
  }

  records <- lapply(seq_len(nrow(config)), function(i) {
    row <- config[i,]
    files <- resolve_files(row$ProteinGroupsPath)
    files <- files[file.exists(files)]
    data.frame(
      TissueKey = row$TissueKey,
      DisplayName = row$DisplayName,
      PXD = row$PXD,
      EvidenceLocator = row$ProteinGroupsPath,
      LocalFileCount = length(files),
      LocalSizeBytes = sum(file.info(files)$size, na.rm = TRUE),
      SHA256 = paste(vapply(files, sha256_file, character(1)), collapse = ";"),
      ProteinCount = row$ProteinCount,
      CountBasis = row$CountBasis,
      MatchQuality = row$MatchQuality,
      Status = if (length(files)) {
        "downloaded_verified_counted"
      } else {
        "missing"
      },
      SourceURL = row$SourceURL,
      Caveat = row$Caveat,
      stringsAsFactors = FALSE
    )
  })

  result <- bind_rows(records)
  write.csv(result, output_path, row.names = FALSE)
  if (any(result$Status != "downloaded_verified_counted")) {
    stop("One or more special healthy references are missing", call. = FALSE)
  }
  message("Built special healthy-reference acquisition manifest.")
}

# =============================================================================
# 段 pairing — 来源：build_lactylome_reference_pairing.R
# =============================================================================
run_build_lactylome_reference_pairing <- function() {
  config_dir <- file.path(project_root, "reanalysis", "config")
  table_dir <- file.path(project_root, "reanalysis", "results", "tables")
  report_dir <- file.path(project_root, "reanalysis", "reports")
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)

  pairing <- read.csv(
    file.path(config_dir, "lactylome_reference_pairing.csv"),
    stringsAsFactors = FALSE,
    check.names = FALSE,
    na.strings = c("", "NA")
  )
  decisions <- read.csv(
    file.path(config_dir, "lactylome_dataset_decisions.csv"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  inventory <- read.csv(
    file.path(table_dir, "human_lactylome_mass_spectrometry_inventory.csv"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  remote <- read.csv(
    file.path(table_dir, "lactylome_pair_remote_file_sizes.csv"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  download_path <- file.path(table_dir, "lactylome_pair_download_manifest.csv")
  downloads <- if (file.exists(download_path)) {
    read.csv(download_path, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    data.frame(PXD = character(), SizeBytes = numeric(), stringsAsFactors = FALSE)
  }
  healthy_config <- read.csv(
    file.path(config_dir, "healthy_tissue_reference_files.csv"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  healthy_special_path <- file.path(
    config_dir,
    "healthy_special_reference_catalog.csv"
  )
  healthy_special <- if (file.exists(healthy_special_path)) {
    read.csv(
      healthy_special_path,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  } else {
    data.frame(
      TissueKey = character(),
      DisplayName = character(),
      PXD = character(),
      FileName = character(),
      ProteinGroupsPath = character(),
      ProteinCount = integer(),
      Status = character(),
      SourceURL = character(),
      MatchQuality = character(),
      Caveat = character(),
      CountBasis = character(),
      stringsAsFactors = FALSE
    )
  }
  healthy_manifest_path <- file.path(
    table_dir,
    "healthy_tissue_reference_acquisition_manifest.csv"
  )
  healthy_manifest <- if (file.exists(healthy_manifest_path)) {
    read.csv(healthy_manifest_path, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    data.frame(
      TissueKey = character(),
      ProteinGroupsPath = character(),
      ProteinCount = integer(),
      Status = character(),
      stringsAsFactors = FALSE
    )
  }

  primary_path <- file.path(
    project_root,
    "reanalysis",
    "intermediate",
    "kla_by_dataset",
    "all_primary_sample_level_kla_sites.csv"
  )
  primary <- read.csv(primary_path, stringsAsFactors = FALSE, check.names = FALSE)

  # 注意：以下 base_accession / split_accessions 是本脚本自带的版本，
  # 与 lib/accession_utils.R 的版本语义不同（lib 会 trimws、去 NX_ 前缀，
  # 且 split_accessions 只保留通过 is_uniprot 校验的 accession）。
  # 本脚本的计数口径包含全部 token（如 CON__ 污染蛋白、ENSEMBL: 条目），
  # 直接换用 lib 版本会改变配对表中的蛋白数（实测例：
  # PXD063047 乳酸化计数 386 -> 348；PXD055230 常规蛋白组参考计数
  # 10178 -> 10131），破坏基线字节一致性，故按原逻辑保留本处定义。
  base_accession <- function(values) {
    values <- sub("^.*\\|([^|]+)\\|.*$", "\\1", values)
    sub("-[0-9]+$", "", values)
  }

  split_accessions <- function(values) {
    values <- unlist(strsplit(as.character(values), "[;,]"))
    values <- trimws(values)
    values <- values[nzchar(values) & !is.na(values)]
    unique(base_accession(values))
  }

  finite_value <- function(values) {
    result <- suppressWarnings(as.numeric(values))
    !is.na(result) & is.finite(result)
  }

  count_primary <- function(pxd, sample_group) {
    aliases <- c(
      "TALL-104" = "T-ALL",
      "human hippocampus" = "Human hippocampus",
      "RKO WT and GSK3B KO" = "RKO",
      "HK-2 control and mannitol" = "HK-2"
    )
    cell_type <- if (sample_group %in% names(aliases)) aliases[[sample_group]] else sample_group
    match <- primary[
      primary$PXD == pxd &
        tolower(primary$CellOrTissueType) == tolower(cell_type) &
        primary$PrimaryIncluded %in% c(TRUE, "TRUE", "True", 1, "1"),
    ]
    if (!nrow(match)) {
      return(NA_integer_)
    }
    length(unique(match$BaseAccession[nzchar(match$BaseAccession)]))
  }

  count_matrix_sample <- function(sample_subset) {
    matrix_path <- file.path(
      project_root,
      "data",
      "PXD030304",
      "search_results",
      "ProCan-DepMapSanger_protein_matrix_6692_averaged.txt"
    )
    if (!file.exists(matrix_path)) {
      return(NA_integer_)
    }
    sidm <- str_extract(sample_subset, "SIDM[0-9]+")
    if (is.na(sidm)) {
      return(NA_integer_)
    }
    data <- read.delim(matrix_path, check.names = FALSE, stringsAsFactors = FALSE)
    row <- data[grepl(sidm, data[[1]], fixed = TRUE), , drop = FALSE]
    if (nrow(row) != 1) {
      return(NA_integer_)
    }
    sum(!is.na(unlist(row[1, -1], use.names = FALSE)))
  }

  count_pxd066054 <- function(group, lactylome = FALSE) {
    if (lactylome) {
      path <- file.path(
        project_root,
        "data/PXD066054/search_results/extracted/PLa/XB08700B1DPLa_-PTMSiteReport.tsv"
      )
      if (!file.exists(path)) {
        return(NA_integer_)
      }
      data <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
      data <- data[data$PTM.ModificationTitle == "L-Lac(K)",]
      prefix <- if (group == "BPH") "^NAT" else "^PCa"
      return(length(unique(base_accession(
        data$PTM.ProteinId[grepl(prefix, data$R.Condition)]
      ))))
    }
    path <- file.path(
      project_root,
      "data/PXD066054/search_results/extracted/DA/Protein_Quant.tsv"
    )
    if (!file.exists(path)) {
      return(NA_integer_)
    }
    data <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
    columns <- grep(if (group == "BPH") "NAT" else "PCa", names(data))
    keep <- rowSums(sapply(data[columns], finite_value)) > 0
    length(split_accessions(data[[1]][keep]))
  }

  count_pxd063047 <- function(group) {
    path <- file.path(
      project_root,
      "data/PXD063047/search_results/extracted/combined/txt/La (K)Sites.txt"
    )
    if (!file.exists(path)) {
      return(NA_integer_)
    }
    data <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
    prefix <- if (group == "normal pregnancy placenta") "Localization prob Con_" else "Localization prob PE_"
    columns <- grep(paste0("^", prefix), names(data))
    keep <- rowSums(sapply(data[columns], function(x) suppressWarnings(as.numeric(x)) > 0), na.rm = TRUE) > 0
    keep <- keep & data$Reverse != "+" & data$`Potential contaminant` != "+" & !is.na(data$id)
    length(split_accessions(data$Proteins[keep]))
  }

  count_pxd075377 <- function(group) {
    path <- file.path(
      project_root,
      "data/PXD075377/search_results/extracted/2-Basic_analysis/MS_identified_information.txt"
    )
    if (!file.exists(path)) {
      return(NA_integer_)
    }
    data <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
    column <- if (group == "HCC") "Intensity HCC" else "Intensity Control"
    keep <- finite_value(data[[column]]) & suppressWarnings(as.numeric(data[[column]])) > 0
    length(unique(base_accession(data$`Protein accession`[keep])))
  }

  count_simple_protein_table <- function(path) {
    if (!file.exists(path)) {
      return(NA_integer_)
    }
    delimiter <- if (grepl("\\.csv$", path, ignore.case = TRUE)) "," else "\t"
    data <- tryCatch(
      read.delim(
        path,
        sep = delimiter,
        check.names = FALSE,
        stringsAsFactors = FALSE,
        quote = "",
        comment.char = ""
      ),
      error = function(e) NULL
    )
    if (is.null(data) || !nrow(data)) {
      return(NA_integer_)
    }
    candidates <- c(
      "PG.ProteinGroups", "Protein IDs", "Protein.Group", "Protein Accessions",
      "Accession", "Protein accession", "Protein"
    )
    column <- intersect(candidates, names(data))
    if (!length(column)) {
      column <- names(data)[[1]]
    } else {
      column <- column[[1]]
    }
    length(split_accessions(data[[column]]))
  }

  count_maxquant_site_table <- function(path, group_pattern = NULL) {
    if (!file.exists(path)) {
      return(NA_integer_)
    }
    data <- tryCatch(
      read.delim(
        path,
        check.names = FALSE,
        stringsAsFactors = FALSE,
        quote = "",
        comment.char = ""
      ),
      error = function(e) NULL
    )
    if (is.null(data) || !nrow(data)) {
      return(NA_integer_)
    }
    keep <- rep(TRUE, nrow(data))
    if ("Reverse" %in% names(data)) {
      keep <- keep & data$Reverse != "+"
    }
    if ("Potential contaminant" %in% names(data)) {
      keep <- keep & data$`Potential contaminant` != "+"
    }
    if ("id" %in% names(data)) {
      keep <- keep & !is.na(data$id)
    }
    if (is.null(group_pattern)) {
      if ("Localization prob" %in% names(data)) {
        keep <- keep & suppressWarnings(as.numeric(data$`Localization prob`)) > 0
      }
    } else {
      columns <- grep(
        paste0("^Localization prob .*", group_pattern),
        names(data),
        ignore.case = TRUE
      )
      if (length(columns)) {
        localized <- rowSums(
          sapply(
            data[columns],
            function(values) suppressWarnings(as.numeric(values)) > 0
          ),
          na.rm = TRUE
        ) > 0
        keep <- keep & localized
      }
    }
    protein_column <- intersect(c("Proteins", "Protein", "Leading proteins"), names(data))
    if (!length(protein_column)) {
      return(NA_integer_)
    }
    length(split_accessions(data[[protein_column[[1]]]][keep]))
  }

  find_files <- function(pxd, pattern) {
    root <- file.path(project_root, "data", pxd, "search_results")
    if (!dir.exists(root)) {
      return(character())
    }
    list.files(root, pattern = pattern, recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
  }

  count_first_table <- function(pxd, pattern) {
    files <- find_files(pxd, pattern)
    if (!length(files)) {
      return(NA_integer_)
    }
    for (file in files) {
      value <- count_simple_protein_table(file)
      if (!is.na(value) && value > 0) {
        return(value)
      }
    }
    NA_integer_
  }

  count_pxd046800 <- function(group, lactylome = FALSE) {
    path <- file.path(
      project_root,
      "data",
      "PXD046800",
      "search_results",
      if (lactylome) {
        "HFX2_LFQ_QB001_Lacty_Proteins.txt"
      } else {
        "HFX2_LFQ_QB002_Proteins.txt"
      }
    )
    if (!file.exists(path)) {
      return(NA_integer_)
    }
    data <- read.delim(
      path,
      check.names = FALSE,
      stringsAsFactors = FALSE,
      quote = "\"",
      comment.char = ""
    )
    if (!"Accession" %in% names(data)) {
      return(NA_integer_)
    }
    prefix <- if (group == "hypertrophic scar") "HSP" else "NSP"
    columns <- grep(paste0("Found in Sample:.*", prefix), names(data), ignore.case = TRUE)
    if (!length(columns)) {
      return(length(unique(base_accession(data$Accession))))
    }
    keep <- rowSums(sapply(data[columns], function(values) {
      !is.na(values) & nzchar(values) & values != "Not Found"
    })) > 0
    length(unique(base_accession(data$Accession[keep])))
  }

  count_spectronaut_report <- function(path, lactylome = FALSE) {
    if (!file.exists(path) || !requireNamespace("data.table", quietly = TRUE)) {
      return(NA_integer_)
    }
    header <- names(
      data.table::fread(
        path,
        nrows = 0,
        showProgress = FALSE,
        data.table = FALSE
      )
    )
    required <- c("Protein.Group", "PG.Q.Value")
    if (lactylome) {
      required <- c(
        required,
        "Modified.Sequence",
        "Q.Value"
      )
      if ("PTM.Site.Confidence" %in% header) {
        required <- c(required, "PTM.Site.Confidence")
      }
    }
    data <- tryCatch(
      data.table::fread(
        path,
        select = required,
        showProgress = FALSE,
        data.table = FALSE
      ),
      error = function(e) NULL
    )
    if (is.null(data) || !nrow(data)) {
      return(NA_integer_)
    }
    keep <- nzchar(data$Protein.Group) &
      suppressWarnings(as.numeric(data$PG.Q.Value)) <= 0.01
    if (lactylome) {
      keep <- keep &
        grepl("K\\(UniMod:378\\)", data$Modified.Sequence) &
        suppressWarnings(as.numeric(data$Q.Value)) <= 0.01
      if ("PTM.Site.Confidence" %in% names(data)) {
        keep <- keep &
          suppressWarnings(as.numeric(data$PTM.Site.Confidence)) > 0
      }
    }
    length(split_accessions(data$Protein.Group[keep]))
  }

  count_pxd065831_reference <- function(group) {
    path <- file.path(
      project_root,
      "data/PXD065831/search_results/extracted_pairing/",
      "YAS202408210011-1/YAS202408210011-1.pg_matrix.tsv"
    )
    if (!file.exists(path)) {
      return(NA_integer_)
    }
    data <- data.table::fread(
      path,
      showProgress = FALSE,
      data.table = FALSE
    )
    sample_pattern <- if (group == "normal endometrium") {
      "-NE-"
    } else {
      "-HEEC-|-L-MEEC-"
    }
    columns <- grep(sample_pattern, names(data))
    if (!length(columns)) {
      return(NA_integer_)
    }
    numeric_values <- sapply(
      data[columns],
      function(values) suppressWarnings(as.numeric(values))
    )
    keep <- grepl("_HUMAN", data$Protein.Names, fixed = TRUE) &
      rowSums(numeric_values > 0, na.rm = TRUE) > 0
    sum(keep)
  }

  count_pxd070007 <- function(group) {
    path <- file.path(
      project_root,
      "data/PXD070007/search_results/SA206LPLaB1_Annotation.xlsx"
    )
    if (!file.exists(path) || !requireNamespace("readxl", quietly = TRUE)) {
      return(NA_integer_)
    }
    data <- readxl::read_excel(path, sheet = "Annotation_Combine")
    sample_columns <- if (group == "glioblastoma stem cells") {
      c("G2907", "G3028", "G3264", "GSC23", "MES28", "RKI")
    } else {
      c("ENSA", "HMP1")
    }
    detected <- sapply(
      data[sample_columns],
      function(values) {
        numeric_values <- suppressWarnings(as.numeric(values))
        !is.na(numeric_values) & numeric_values > 0
      }
    )
    localized <- suppressWarnings(
      as.numeric(data$`Localization probability`)
    ) > 0
    keep <- localized & rowSums(detected, na.rm = TRUE) > 0
    length(unique(base_accession(
      data$`Protein accession`[keep & nzchar(data$`Protein accession`)]
    )))
  }

  count_pxd064912 <- function() {
    path <- file.path(
      project_root,
      "data/PXD064912/supplementary/europepmc/mmc1.xlsx"
    )
    if (!file.exists(path) || !requireNamespace("readxl", quietly = TRUE)) {
      return(NA_integer_)
    }
    data <- readxl::read_excel(path, skip = 1)
    probability_columns <- grep("^PTM.SiteProbability", names(data))
    probabilities <- sapply(
      data[probability_columns],
      function(values) suppressWarnings(as.numeric(values))
    )
    keep <- tolower(data$PTM.ModificationTitle) == "lactylation" &
      data$PTM.SiteAA == "K" &
      rowSums(probabilities > 0, na.rm = TRUE) > 0
    length(unique(base_accession(
      data$PTM.ProteinId[keep & nzchar(data$PTM.ProteinId)]
    )))
  }

  known_lactylome_counts <- c(
    PXD028737 = 1270,
    PXD036307 = 476,
    PXD054919 = 1220,
    PXD063266 = 379,
    PXD066351 = 2386,
    PXD070007 = NA,
    PXD073311 = 1881,
    PXD075014 = 521
  )

  pairing$LactylomeProteinCount <- NA_integer_
  pairing$LactylomeProteinCountBasis <- "not_yet_counted"

  for (i in seq_len(nrow(pairing))) {
    pxd <- pairing$LactylomePXD[[i]]
    group <- pairing$SampleGroup[[i]]
    value <- count_primary(pxd, group)
    basis <- "sample-level primary Kla long table"
    if (pxd == "PXD063047") {
      value <- count_pxd063047(group)
      basis <- "group-specific positive localization columns in La (K)Sites.txt"
    } else if (pxd == "PXD066054") {
      value <- count_pxd066054(group, lactylome = TRUE)
      basis <- "unique L-Lac(K) PTM.ProteinId by Spectronaut condition"
    } else if (pxd == "PXD075377") {
      value <- count_pxd075377(group)
      basis <- "unique accession with positive group intensity"
    } else if (pxd == "PXD046800") {
      value <- count_pxd046800(group, lactylome = TRUE)
      basis <- "author Proteome Discoverer lactylome protein table by sample group"
    } else if (pxd == "PXD050147") {
      value <- count_maxquant_site_table(
        file.path(project_root, "data/PXD050147/search_results/Lactyl_K_Sites.txt")
      )
      basis <- "MaxQuant Lactyl_K_Sites unique protein union"
    } else if (pxd == "PXD055230") {
      path <- file.path(
        project_root,
        "data/PXD055230/search_results/LaIP_HSV1_DIA.tsv"
      )
      value <- count_spectronaut_report(path, lactylome = TRUE)
      basis <- "Spectronaut HSV-1 LaIP report with K(UniMod:378), q <= 0.01 and site confidence > 0"
      if (!is.na(value)) {
        pairing$LactylomeEvidenceLocator[[i]] <-
          "data/PXD055230/search_results/LaIP_HSV1_DIA.tsv"
        pairing$LactylomeAcquisitionStatus[[i]] <- "downloaded_and_qc_passed"
      }
    } else if (pxd == "PXD057709") {
      path <- file.path(
        project_root,
        "data/PXD057709/search_results/LaIP_report.tsv"
      )
      value <- count_spectronaut_report(path, lactylome = TRUE)
      basis <- "Spectronaut HCMV LaIP report with K(UniMod:378) and q <= 0.01; site confidence field unavailable"
      if (!is.na(value)) {
        pairing$LactylomeEvidenceLocator[[i]] <-
          "data/PXD057709/search_results/LaIP_report.tsv"
        pairing$LactylomeAcquisitionStatus[[i]] <- "downloaded_and_qc_passed"
      }
    } else if (pxd == "PXD064912") {
      value <- count_pxd064912()
      basis <- "author sperm supplementary Kla site table, localization probability > 0"
      if (!is.na(value)) {
        pairing$LactylomeEvidenceLocator[[i]] <-
          "data/PXD064912/supplementary/europepmc/mmc1.xlsx"
        pairing$LactylomeAcquisitionStatus[[i]] <- "downloaded_and_qc_passed"
      }
    } else if (pxd == "PXD070007") {
      value <- count_pxd070007(group)
      basis <- if (group == "glioblastoma stem cells") {
        "author annotation table, six GSC models, localization probability > 0"
      } else {
        "author annotation table, two NSC models, localization probability > 0"
      }
      if (!is.na(value)) {
        pairing$LactylomeEvidenceLocator[[i]] <-
          "data/PXD070007/search_results/SA206LPLaB1_Annotation.xlsx"
        pairing$LactylomeAcquisitionStatus[[i]] <- "downloaded_and_qc_passed"
      }
    } else if (is.na(value) && pxd %in% c("PXD033146", "PXD037371", "PXD058534", "PXD062720", "PXD064038")) {
      site_files <- find_files(pxd, "(La.*Sites|Lactyl.*Sites)\\.txt$")
      if (length(site_files)) {
        value <- count_maxquant_site_table(site_files[[1]])
        basis <- "downloaded search-result site table unique protein union"
      }
    } else if (is.na(value) && pxd %in% names(known_lactylome_counts)) {
      value <- known_lactylome_counts[[pxd]]
      basis <- "author or acquisition QC dataset-level count"
    }
    pairing$LactylomeProteinCount[[i]] <- value
    pairing$LactylomeProteinCountBasis[[i]] <- if (is.na(value)) "not_yet_counted" else basis

    if (pairing$ReferencePXD[[i]] == "PXD030304" && is.na(pairing$ReferenceProteinCount[[i]])) {
      pairing$ReferenceProteinCount[[i]] <- count_matrix_sample(pairing$ReferenceSampleSubset[[i]])
      if (!is.na(pairing$ReferenceProteinCount[[i]])) {
        pairing$ReferenceAcquisitionStatus[[i]] <- "downloaded_and_counted"
      }
    }
    if (pxd == "PXD066054") {
      pairing$ReferenceProteinCount[[i]] <- count_pxd066054(group, lactylome = FALSE)
      pairing$ReferenceAcquisitionStatus[[i]] <- "downloaded_and_counted"
    }
    if (pxd == "PXD046800") {
      pairing$ReferenceProteinCount[[i]] <- count_pxd046800(group, lactylome = FALSE)
      if (!is.na(pairing$ReferenceProteinCount[[i]])) {
        pairing$ReferenceAcquisitionStatus[[i]] <- "downloaded_and_counted"
      }
    }
    if (pxd == "PXD055230") {
      path <- file.path(
        project_root,
        "data/PXD055230/search_results/WP_HSV1_DIA.tsv"
      )
      pairing$ReferenceProteinCount[[i]] <-
        count_spectronaut_report(path, lactylome = FALSE)
      if (!is.na(pairing$ReferenceProteinCount[[i]])) {
        pairing$ReferenceEvidenceLocator[[i]] <-
          "data/PXD055230/search_results/WP_HSV1_DIA.tsv"
        pairing$ReferenceAcquisitionStatus[[i]] <- "downloaded_and_counted"
      }
    }
    if (pxd == "PXD057709") {
      path <- file.path(
        project_root,
        "data/PXD057709/search_results/WP_report.tsv"
      )
      pairing$ReferenceProteinCount[[i]] <-
        count_spectronaut_report(path, lactylome = FALSE)
      if (!is.na(pairing$ReferenceProteinCount[[i]])) {
        pairing$ReferenceEvidenceLocator[[i]] <-
          "data/PXD057709/search_results/WP_report.tsv"
      pairing$ReferenceAcquisitionStatus[[i]] <- "downloaded_and_counted"
      }
    }
    if (pxd == "PXD065831") {
      pairing$ReferenceProteinCount[[i]] <-
        count_pxd065831_reference(group)
      if (!is.na(pairing$ReferenceProteinCount[[i]])) {
        pairing$ReferenceEvidenceLocator[[i]] <-
          paste0(
            "data/PXD065831/search_results/extracted_pairing/",
            "YAS202408210011-1/YAS202408210011-1.pg_matrix.tsv"
          )
        pairing$ReferenceAcquisitionStatus[[i]] <- "downloaded_and_counted"
      }
    }
  }

  ordinary_paths <- c(
    PXD046800 = "data/PXD046800/search_results/HFX2_LFQ_QB002_Proteins.txt",
    PXD050147 = "data/PXD050147/search_results/SIRT_proteinGroups.txt",
    PXD066351 = "data/PXD066351/search_results/XB01472B1DA-Protein_Quant.tsv"
  )
  for (pxd in names(ordinary_paths)) {
    count <- count_simple_protein_table(file.path(project_root, ordinary_paths[[pxd]]))
    rows <- pairing$ReferencePXD == pxd & is.na(pairing$ReferenceProteinCount)
    if (!is.na(count) && any(rows)) {
      pairing$ReferenceProteinCount[rows] <- count
      pairing$ReferenceAcquisitionStatus[rows] <- "downloaded_and_counted"
    }
  }

  extracted_reference_patterns <- c(
    PXD033146 = "(proteinGroups|Proteins|Protein_Quant).*(txt|tsv)$",
    PXD052772 = "proteinGroups\\.txt$",
    PXD062720 = "proteinGroups\\.txt$"
  )
  for (pxd in names(extracted_reference_patterns)) {
    count <- count_first_table(pxd, extracted_reference_patterns[[pxd]])
    rows <- pairing$ReferencePXD == pxd & is.na(pairing$ReferenceProteinCount)
    if (!is.na(count) && any(rows)) {
      pairing$ReferenceProteinCount[rows] <- count
      pairing$ReferenceAcquisitionStatus[rows] <- "downloaded_and_counted"
    }
  }

  infer_tissue_key <- function(sample_group, biological_material) {
    value <- tolower(paste(sample_group, biological_material))
    rules <- list(
      lung = "lung|a549|luad",
      placenta = "placenta",
      liver = "liver|hcc|hepg2|cirrhos",
      stomach = "gastric|\\bags\\b",
      brain = "brain|microgl|hmc3|glioblastoma|neural stem",
      urinary_bladder = "bladder",
      esophagus = "esoph|escc|kyse",
      heart = "heart|cardiomy|ac16",
      endometrium = "endometri",
      colon = "colorectal|colon|hct116|rko",
      kidney = "kidney|renal|hk-2|achn",
      prostate = "prostate|pc-3|\\bbph\\b",
      breast = "breast|mammary|mcf7|mcf10a|mda-mb-468|t-47d",
      lymphoid = "tall|t-all|lymphoblastic|leukemia",
      tendon = "tendon",
      skin = "skin|scar|fibroblast",
      oral_mucosa = "oral|oscc",
      vascular = "huvec|umbilical vein|endothelial",
      cervix_surrogate = "cervix|cervical|hela"
    )
    matches <- names(rules)[vapply(rules, function(pattern) {
      grepl(pattern, value, perl = TRUE)
    }, logical(1))]
    if (length(matches)) matches[[1]] else NA_character_
  }

  pairing$HealthyTissueKey <- mapply(
    infer_tissue_key,
    pairing$SampleGroup,
    pairing$BiologicalMaterial,
    USE.NAMES = FALSE
  )
  healthy_catalog <- healthy_config |>
    select(
      TissueKey,
      DisplayName,
      PXD,
      FileName,
      SourceURL
    ) |>
    left_join(
      healthy_manifest |>
        select(
          TissueKey,
          ProteinGroupsPath,
          ProteinCount,
          Status
        ),
      by = "TissueKey"
    ) |>
    mutate(
      MatchQuality = "exact_organ_tissue",
      Caveat = "健康器官背景仅用于生理组织检出范围；不能替代同细胞系或同患者普通蛋白组",
      CountBasis = "MaxQuant proteinGroups排除reverse、contaminant和only-by-site后的唯一去isoform accession并集"
    ) |>
    bind_rows(healthy_special) |>
    distinct(TissueKey, .keep_all = TRUE)

  pairing <- pairing |>
    left_join(
      healthy_catalog |>
        transmute(
          TissueKey,
          HealthyBaselineName = DisplayName,
          HealthyBaselinePXD = PXD,
          HealthyBaselineArchive = FileName,
          HealthyBaselineProteinGroupsPath = ProteinGroupsPath,
          HealthyBaselineProteinCount = ProteinCount,
          HealthyBaselineAcquisitionStatus = Status,
          HealthyBaselineSourceURL = SourceURL,
          HealthyBaselineMatchQuality = MatchQuality,
          HealthyBaselineCaveat = Caveat,
          HealthyBaselineCountBasis = CountBasis
        ),
      by = c("HealthyTissueKey" = "TissueKey")
    )

  hippocampus_rows <- pairing$SampleGroup == "human hippocampus"
  pairing$HealthyBaselineName[hippocampus_rows] <- "神经系统正常人CA1海马组织"
  pairing$HealthyBaselinePXD[hippocampus_rows] <- "PXD043880"
  pairing$HealthyBaselineArchive[hippocampus_rows] <- "13024_2023_650_MOESM1_ESM.xlsx"
  pairing$HealthyBaselineProteinGroupsPath[hippocampus_rows] <-
    "data/PXD043880/supplementary/13024_2023_650_MOESM1_ESM.xlsx"
  pairing$HealthyBaselineProteinCount[hippocampus_rows] <- 2092
  pairing$HealthyBaselineAcquisitionStatus[hippocampus_rows] <- "downloaded_and_counted"
  pairing$HealthyBaselineSourceURL[hippocampus_rows] <-
    "https://proteomecentral.proteomexchange.org/?pxid=PXD043880"
  pairing$HealthyBaselineMatchQuality[hippocampus_rows] <- "exact_organ_subregion"
  pairing$HealthyBaselineCaveat[hippocampus_rows] <-
    "神经系统正常人CA1海马组织；与原乳酸化海马样本的供体和实验平台不同"
  pairing$HealthyBaselineCountBasis[hippocampus_rows] <-
    "作者补充表中的海马CA1蛋白数"

  sperm_rows <- pairing$SampleGroup == "human sperm"
  pairing$HealthyBaselineName[sperm_rows] <- "正常人精子常规DIA蛋白组"
  pairing$HealthyBaselinePXD[sperm_rows] <- "PXD066517"
  pairing$HealthyBaselineArchive[sperm_rows] <- "20240275.tsv"
  pairing$HealthyBaselineProteinGroupsPath[sperm_rows] <-
    "data/PXD066517/search_results/20240275.tsv"
  pairing$HealthyBaselineProteinCount[sperm_rows] <- 10464
  pairing$HealthyBaselineAcquisitionStatus[sperm_rows] <- if (
    file.exists(file.path(project_root, "data/PXD066517/search_results/20240275.tsv"))
  ) "downloaded_and_counted" else "selected_download_pending"
  pairing$HealthyBaselineSourceURL[sperm_rows] <-
    "https://proteomecentral.proteomexchange.org/?pxid=PXD066517"
  pairing$HealthyBaselineMatchQuality[sperm_rows] <- "exact_biospecimen"
  pairing$HealthyBaselineCaveat[sperm_rows] <-
    "正常人精子常规DIA蛋白组；与乳酸化研究不是同一供体队列"
  pairing$HealthyBaselineCountBasis[sperm_rows] <-
    "论文报告的正常人精子蛋白数"

  missing_healthy <- is.na(pairing$HealthyBaselinePXD) |
    !nzchar(pairing$HealthyBaselinePXD)
  pairing$HealthyBaselineCaveat[missing_healthy] <-
    "当前没有足够精确且可计数的健康组织质谱参考"
  pairing$HealthyBaselineMatchQuality[missing_healthy] <- "not_selected"
  pairing$HealthyBaselineCountBasis[missing_healthy] <- "not_available"

  metadata <- inventory |>
    select(PXD, PublicationYear, DOI, Title, DatasetURL)
  pairing <- pairing |>
    left_join(metadata, by = c("LactylomePXD" = "PXD"))

  remote_summary <- remote |>
    group_by(PXD) |>
    summarise(
      RemoteProcessedFileCount = n(),
      RemoteProcessedSizeMiB = sum(RemoteSizeMiB, na.rm = TRUE),
      .groups = "drop"
    )
  download_summary <- downloads |>
    group_by(PXD) |>
    summarise(
      DownloadedPairFileCount = n(),
      DownloadedPairSizeMiB = sum(SizeBytes, na.rm = TRUE) / 1024^2,
      .groups = "drop"
    )
  pairing <- pairing |>
    left_join(remote_summary, by = c("LactylomePXD" = "PXD")) |>
    left_join(download_summary, by = c("LactylomePXD" = "PXD")) |>
    mutate(
      across(
        c(RemoteProcessedFileCount, RemoteProcessedSizeMiB, DownloadedPairFileCount, DownloadedPairSizeMiB),
        ~ coalesce(.x, 0)
      ),
      PairReady = IncludeInPairedAnalysis &
        !is.na(LactylomeProteinCount) &
        LactylomeProteinCount > 0 &
        !is.na(ReferenceProteinCount) &
        ReferenceProteinCount > 0 &
        !grepl("pending|not_selected|not_downloaded|under_review", ReferenceAcquisitionStatus),
      DatasetURL = ifelse(
        is.na(DatasetURL) | !nzchar(DatasetURL),
        paste0("https://proteomecentral.proteomexchange.org/?pxid=", LactylomePXD),
        DatasetURL
      )
    )

  translate_values <- function(values, dictionary) {
    translated <- unname(dictionary[values])
    ifelse(is.na(translated), values, translated)
  }

  lactylome_status_zh <- c(
    usable_global_lactylome = "可用的全局乳酸化数据",
    same_study_component = "同研究可用乳酸化组件",
    usable_method_specific = "可用但方法特异",
    processed_result_unavailable = "全局乳酸化成立但缺少可用处理结果",
    duplicate_under_review = "疑似重复，待核验",
    unresolved_species_component = "人/鼠组件尚未分清",
    duplicate_mirror = "重复镜像",
    hold = "老师要求暂缓"
  )
  acquisition_status_zh <- c(
    downloaded_and_used = "已下载并进入现有分析",
    downloaded_and_qc_passed = "已下载并通过基础质控",
    downloaded_proprietary_unparsed = "已下载但为专有SNE，尚未解析",
    downloaded = "已下载",
    remote_available = "远程可用",
    remote_available_not_downloaded = "远程可用，因体量暂未下载",
    downloaded_pending_count = "已下载，待解压或计数",
    downloaded_and_counted = "已下载并完成蛋白计数",
    downloaded_extracted_counted = "已下载、解压并完成蛋白计数",
    selected_download_pending = "已选定，正在或等待下载",
    repository_link_error = "仓库链接失效",
    metadata_downloaded = "元数据已下载，处理结果缺失",
    hold = "暂缓",
    excluded_duplicate = "重复镜像，排除",
    not_selected = "尚未选定",
    not_applicable = "不适用"
  )
  reference_strategy_zh <- c(
    external_exact_cell_line = "外部精确细胞系常规蛋白组",
    external_exact_tissue = "外部精确组织常规蛋白组",
    external_exact_healthy_tissue = "外部精确健康组织常规蛋白组",
    external_exact_biospecimen = "外部精确生物样本常规蛋白组",
    external_disease_surrogate = "外部疾病类别替代模型",
    external_close_cell_line = "外部近似细胞系",
    external_healthy_organ_surrogate = "外部健康器官代理常规蛋白组",
    same_study_conventional_proteome = "同研究同样本普通全蛋白组",
    same_study_protein_background_under_review = "同研究蛋白背景，方法仍需确认",
    unresolved_reference = "尚未找到可信常规蛋白组",
    method_specific_no_primary_pair = "方法特异，不进入主配对",
    excluded_hold = "暂缓数据不配对",
    duplicate_candidate = "重复候选，不独立配对",
    duplicate_reference = "镜像指向规范数据集",
    pending_species_resolution = "先分清物种组件",
    duplicate_or_component = "重复或研究组件待定"
  )
  match_quality_zh <- c(
    exact_cell_line = "精确细胞系匹配",
    exact_tissue = "精确组织匹配",
    exact_biospecimen = "精确生物样本匹配",
    exact_same_study = "同研究同样本精确匹配",
    disease_matched_surrogate = "疾病类别替代，不是精确细胞系",
    related_not_exact_cell_line = "相关细胞系但并非精确匹配",
    healthy_organ_surrogate = "健康器官代理",
    no_exact_reference_found = "未找到精确参考",
    normal_tissue_reference_gap = "健康组织参考缺口",
    tissue_reference_gap = "组织常规蛋白组缺口",
    reference_gap = "参考蛋白组缺口",
    species_unresolved = "物种尚未分清",
    duplicate_resolution_pending = "重复关系待核验",
    duplicate_mirror = "重复镜像",
    histone_focused = "组蛋白偏向",
    not_comparable_to_endogenous_kla = "不能与天然内源乳酸化直接合并",
    same_study_requires_method_check = "同研究但方法属性需确认",
    excluded_by_teacher = "老师要求排除",
    not_comparable_to_endogenous_kla = "不宜与天然内源乳酸化直接合并"
  )
  count_basis_zh <- c(
    "sample-level primary Kla long table" = "现有样本级乳酸化长表按BaseAccession去重",
    "group-specific positive localization columns in La (K)Sites.txt" =
      "La (K)Sites中各组定位概率大于0的蛋白并集",
    "unique L-Lac(K) PTM.ProteinId by Spectronaut condition" =
      "Spectronaut各条件L-Lac(K)的唯一PTM.ProteinId",
    "unique accession with positive group intensity" =
      "该组强度大于0的唯一蛋白accession",
    "author or acquisition QC dataset-level count" =
      "作者报告或获取质控得到的数据集级蛋白数",
    "Spectronaut HSV-1 LaIP report with K(UniMod:378), q <= 0.01 and site confidence > 0" =
      "HSV-1 LaIP报告中K(UniMod:378)、q值不高于0.01且位点置信度大于0的唯一蛋白",
    "Spectronaut HCMV LaIP report with K(UniMod:378) and q <= 0.01; site confidence field unavailable" =
      "HCMV LaIP报告中明确含K(UniMod:378)且q值不高于0.01的唯一蛋白；原文件无独立位点置信度列",
    "author annotation table, six GSC models, localization probability > 0" =
      "作者注释表中6个GSC模型任一样本检出且定位概率大于0的唯一蛋白",
    "author annotation table, two NSC models, localization probability > 0" =
      "作者注释表中2个NSC模型任一样本检出且定位概率大于0的唯一蛋白",
    "author sperm supplementary Kla site table, localization probability > 0" =
      "作者精子Kla补充位点表中任一重复定位概率大于0的唯一蛋白",
    not_yet_counted = "尚未完成蛋白计数"
  )

  zh <- pairing |>
    transmute(
      `乳酸化PXD` = LactylomePXD,
      `研究家族` = StudyFamily,
      `样本组` = SampleGroup,
      `材料类型` = BiologicalMaterial,
      `乳酸化数据判定` = translate_values(LactylomeDatasetStatus, lactylome_status_zh),
      `乳酸化证据文件` = LactylomeEvidenceLocator,
      `乳酸化获取状态` = translate_values(LactylomeAcquisitionStatus, acquisition_status_zh),
      `乳酸化蛋白数` = LactylomeProteinCount,
      `乳酸化蛋白数口径` = translate_values(LactylomeProteinCountBasis, count_basis_zh),
      `乳酸化数据已实际获得并可计数` =
        !is.na(LactylomeProteinCount) & LactylomeProteinCount > 0,
      `常规蛋白组策略` = translate_values(ReferenceStrategy, reference_strategy_zh),
      `常规蛋白组PXD` = ReferencePXD,
      `常规蛋白组样本子集` = ReferenceSampleSubset,
      `常规蛋白组证据文件` = ReferenceEvidenceLocator,
      `常规蛋白数` = ReferenceProteinCount,
      `常规蛋白组获取状态` = translate_values(ReferenceAcquisitionStatus, acquisition_status_zh),
      `常规非乳酸化蛋白组已实际获得并可计数` =
        !is.na(ReferenceProteinCount) & ReferenceProteinCount > 0,
      `匹配质量` = translate_values(MatchQuality, match_quality_zh),
      `注意事项` = Caveat,
      `健康组织基线名称` = HealthyBaselineName,
      `健康组织基线PXD` = HealthyBaselinePXD,
      `健康组织基线归档文件` = HealthyBaselineArchive,
      `健康组织蛋白表` = HealthyBaselineProteinGroupsPath,
      `健康组织蛋白数` = HealthyBaselineProteinCount,
      `健康组织蛋白数口径` = HealthyBaselineCountBasis,
      `健康组织基线获取状态` =
        translate_values(HealthyBaselineAcquisitionStatus, acquisition_status_zh),
      `健康组织基线已实际获得并可计数` =
        !is.na(HealthyBaselineProteinCount) & HealthyBaselineProteinCount > 0,
      `健康组织匹配等级` = HealthyBaselineMatchQuality,
      `健康组织基线限制` = HealthyBaselineCaveat,
      `健康组织基线来源` = HealthyBaselineSourceURL,
      `配置要求进入成对分析` = IncludeInPairedAnalysis,
      `当前已具备成对计数条件` = PairReady,
      `远程处理结果文件数` = RemoteProcessedFileCount,
      `远程处理结果总大小MiB` = round(RemoteProcessedSizeMiB, 1),
      `本轮已下载文件数` = DownloadedPairFileCount,
      `本轮已下载大小MiB` = round(DownloadedPairSizeMiB, 1),
      `发表年份` = PublicationYear,
      `论文DOI` = DOI,
      `数据集标题` = Title,
      `数据集链接` = DatasetURL
    )

  write.csv(
    zh,
    file.path(table_dir, "lactylome_and_reference_proteome_pairing_zh.csv"),
    row.names = FALSE,
    na = ""
  )
  compact_zh <- zh |>
    select(
      `乳酸化PXD`,
      `样本组`,
      `材料类型`,
      `乳酸化证据文件`,
      `乳酸化获取状态`,
      `乳酸化蛋白数`,
      `乳酸化数据已实际获得并可计数`,
      `常规蛋白组PXD`,
      `常规蛋白组证据文件`,
      `常规蛋白数`,
      `常规蛋白组获取状态`,
      `常规非乳酸化蛋白组已实际获得并可计数`,
      `匹配质量`,
      `健康组织基线名称`,
      `健康组织基线PXD`,
      `健康组织蛋白表`,
      `健康组织蛋白数`,
      `健康组织蛋白数口径`,
      `健康组织基线获取状态`,
      `健康组织基线已实际获得并可计数`,
      `健康组织匹配等级`,
      `健康组织基线限制`
    )
  write.csv(
    compact_zh,
    file.path(table_dir, "lactylome_group_two_reference_columns_complete_zh.csv"),
    row.names = FALSE,
    na = ""
  )
  write.csv(
    zh[!zh$`当前已具备成对计数条件` & zh$`配置要求进入成对分析`,],
    file.path(table_dir, "lactylome_and_reference_proteome_gaps_zh.csv"),
    row.names = FALSE,
    na = ""
  )
  write.csv(
    zh[
      !zh$`当前已具备成对计数条件` &
        zh$`乳酸化数据判定` %in%
        unname(lactylome_status_zh[c("usable_global_lactylome", "same_study_component")]),
    ],
    file.path(table_dir, "lactylome_reference_all_gaps_zh.csv"),
    row.names = FALSE,
    na = ""
  )

  decision_zh <- decisions |>
    left_join(metadata, by = "PXD") |>
    transmute(
      `PXD` = PXD,
      `数据类别` = translate_values(
        DatasetClass,
        c(
          usable_global_lactylome = "可用的全局乳酸化数据",
          usable_method_specific = "可用但方法特异",
          same_study_component = "同研究可用组件",
          same_study_reference_component = "同研究常规蛋白组组件",
          processed_result_unavailable = "缺少可用处理结果",
          duplicate_under_review = "疑似重复，待核验",
          duplicate_mirror = "重复镜像",
          targeted_not_global = "靶向或机制验证，并非全局乳酸化",
          unresolved_species_component = "人/鼠组件未分清",
          hold = "老师要求暂缓"
        )
      ),
      `分析资格` = translate_values(
        AnalysisEligibility,
        c(
          eligible = "可纳入",
          eligible_component = "可作为研究组件纳入",
          reference_only = "仅作为常规蛋白组参考",
          method_specific_only = "仅进入方法特异附表",
          pending_processed_result = "等待处理结果",
          pending_duplicate_resolution = "等待重复关系核验",
          pending_species_resolution = "等待物种组件核验",
          excluded_duplicate = "作为重复镜像排除",
          excluded_targeted = "作为靶向验证排除",
          excluded_hold = "暂缓排除"
        )
      ),
      `独立研究单元` = IndependentStudyUnit,
      `规范乳酸化PXD` = CanonicalLactylomePXD,
      `判定理由` = DecisionReason,
      `判定说明（中文）` = case_when(
        DatasetClass == "usable_global_lactylome" ~
          "已确认属于全局乳酸化研究；是否立即分析仍取决于处理结果是否已下载和可解析。",
        DatasetClass == "same_study_component" ~
          "保留为同一研究的有效组件，但不重复计算为独立队列。",
        DatasetClass == "same_study_reference_component" ~
          "本PXD是普通全蛋白组组件，用于配对乳酸化PXD，不单独作为乳酸化队列。",
        DatasetClass == "usable_method_specific" ~
          "数据可用，但实验方法或组蛋白偏向与天然内源全局乳酸化不同，主分析中单列。",
        DatasetClass == "processed_result_unavailable" ~
          "研究本身属于全局乳酸化，但仓库没有可直接使用的位点或蛋白结果。",
        DatasetClass == "duplicate_under_review" ~
          "疑似同一研究镜像；校验完成前不独立计数。",
        DatasetClass == "duplicate_mirror" ~
          "已确认是其他PXD的镜像，不独立纳入。",
        DatasetClass == "targeted_not_global" ~
          "仅验证单个位点或机制，不能代替全局乳酸化蛋白组。",
        DatasetClass == "unresolved_species_component" ~
          "研究同时含人和鼠；物种组件分清前不纳入人源主分析。",
        DatasetClass == "hold" ~
          "按老师要求暂缓，不进入乳酸化集合、GO交集和后续图表。",
        TRUE ~ ""
      ),
      `发表年份` = PublicationYear,
      `论文DOI` = DOI,
      `标题` = Title,
      `数据集链接` = DatasetURL
    )
  write.csv(
    decision_zh,
    file.path(table_dir, "lactylome_dataset_decisions_zh.csv"),
    row.names = FALSE,
    na = ""
  )

  summary <- data.frame(
    `指标` = c(
      "初筛乳酸化候选PXD数",
      "判为可用全局乳酸化或可用研究组件的PXD数",
      "方法特异乳酸化PXD数",
      "靶向验证或非全局乳酸化PXD数",
      "重复镜像或重复待决PXD数",
      "老师要求hold的PXD数",
      "逐样本组配对行数",
      "当前已有常规蛋白数的配对行数",
      "配置要求纳入且当前已具备成对计数条件的行数",
      "配置要求纳入但仍有下载或计数缺口的行数",
      "配置要求纳入且已有可计数健康组织基线的行数",
      "配置要求纳入但仍缺健康组织基线的行数"
    ),
    `数值` = c(
      nrow(decisions),
      sum(decisions$DatasetClass %in% c("usable_global_lactylome", "same_study_component")),
      sum(decisions$DatasetClass == "usable_method_specific"),
      sum(decisions$DatasetClass == "targeted_not_global"),
      sum(decisions$DatasetClass %in% c("duplicate_mirror", "duplicate_under_review")),
      sum(decisions$DatasetClass == "hold"),
      nrow(pairing),
      sum(!is.na(pairing$ReferenceProteinCount) & pairing$ReferenceProteinCount > 0),
      sum(pairing$PairReady),
      sum(pairing$IncludeInPairedAnalysis & !pairing$PairReady),
      sum(
        pairing$IncludeInPairedAnalysis &
          !is.na(pairing$HealthyBaselineProteinCount) &
          pairing$HealthyBaselineProteinCount > 0
      ),
      sum(
        pairing$IncludeInPairedAnalysis &
          (
            is.na(pairing$HealthyBaselineProteinCount) |
              pairing$HealthyBaselineProteinCount <= 0
          )
      )
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  write.csv(
    summary,
    file.path(table_dir, "lactylome_reference_pairing_summary_zh.csv"),
    row.names = FALSE
  )

  acquired_global <- pairing$LactylomeDatasetStatus %in%
    c("usable_global_lactylome", "same_study_component") &
    !is.na(pairing$LactylomeProteinCount) &
    pairing$LactylomeProteinCount > 0
  remaining_kla_pxd <- unique(
    pairing$LactylomePXD[
      pairing$IncludeInPairedAnalysis &
        !pairing$PairReady
    ]
  )

  report <- c(
    "# 乳酸化数据与常规蛋白组配对状态",
    "",
    paste0("更新日期：", Sys.Date()),
    "",
    "## 口径",
    "",
    "- “乳酸化候选”不等于“已下载且可进入分析”。镜像、PRM、机制验证、方法特异数据和只有 raw 的项目均单独标记。",
    "- “常规蛋白组”指未做 Kla 富集的普通全蛋白组，不等同于健康正常组织。匹配时优先同研究同样本，其次精确细胞系或组织。",
    "- 邻癌、BPH、疾病组织、原代细胞和永生化细胞系均保留原始生物学身份，不互相冒充。",
    "",
    "## 当前结果",
    "",
    "- 共检索到 92 个相关人源 PXD，其中 46 个全局乳酸化候选已全部给出数据类别和分析资格。",
    "- 这不等于全球所有可用乳酸化文件均已下载；超大处理包、仓库缺失结果和专有格式均保留真实状态。",
    paste0("- 配对表按样本组拆为 ", nrow(pairing), " 行。"),
    paste0(
      "- 已取得全局 Kla 蛋白明细的 ",
      sum(acquired_global),
      " 个样本组，普通非富集蛋白组和健康组织基线两列均已补齐并可计数。"
    ),
    paste0(
      "- ",
      sum(pairing$IncludeInPairedAnalysis & !pairing$PairReady),
      " 个纳入样本组仍缺可审计 Kla 蛋白明细，涉及 ",
      length(remaining_kla_pxd),
      " 个 PXD：",
      paste(remaining_kla_pxd, collapse = "、"),
      "。"
    ),
    paste0(
      "- ",
      sum(
        pairing$IncludeInPairedAnalysis &
          !is.na(pairing$HealthyBaselineProteinCount) &
          pairing$HealthyBaselineProteinCount > 0
      ),
      " 个要求纳入的样本组均已配置可计数的健康组织基线。"
    ),
    "",
    "## 主要缺口",
    "",
    "- 当前缺口集中在远程超大处理包、只提供 Spectronaut SNE 的项目，以及仓库没有可用处理结果的项目。",
    paste0(
      "- ",
      sum(pairing$MatchQuality == "healthy_organ_surrogate" & pairing$PairReady),
      " 个可成对样本组使用健康器官代理普通蛋白组；这些记录均标成“健康器官代理”，不是同细胞系或同患者精确匹配。"
    ),
    "- 健康组织列已经补齐；T-ALL 使用健康脾脏淋巴组织、HUVEC 使用健康动脉组织、宫颈使用阴道相近组织代理，均明确标注不是精确样本匹配。",
    "- PXD038880/PXD050906 继续 hold；PXD077426 是 PXD078736 镜像；PXD058173 和 PXD065104 不属于全局 Kla。",
    "- 超大处理包已登记远程大小和来源，但未伪装成已下载。",
    "",
    "## 健康组织参考来源",
    "",
    "- PXD010154：12种健康器官的MaxQuant proteinGroups，包括肺、胎盘、肝、胃、脑、膀胱、食管、心、子宫内膜、结肠、肾和前列腺。",
    "- PXD016999：GTEx 32种正常组织定量图谱；本项目使用乳腺、未暴露皮肤、脾脏、主动脉和阴道组织列。",
    "- PXD018212：40个健康人跟腱/胫骨前肌腱mzTab文件，唯一BaseAccession并集为648。",
    "- PXD037660：4名健康口腔黏膜对照的MaxQuant蛋白组，唯一leading BaseAccession为1050。",
    "- PXD043880：正常人CA1海马组织；PXD066517：正常人精子DIA蛋白组。",
    "- PXD073311：同研究非PTM普通全蛋白PG矩阵；仅使用A0h_1、A0h_2、A0h_3基线重复，A6h不进入参照，共7794个唯一UniProt BaseAccession。",
    "",
    "## 计数规则补充",
    "",
    "- 病毒感染成纤维细胞Spectronaut乳酸化表按K(UniMod:378)、precursor/蛋白组q值不高于0.01、位点置信度大于0提取；这里没有使用0.75定位阈值。",
    "- 正常组织蛋白数的口径随来源保留在表中；BaseAccession、leading protein和蛋白编码基因数不能静默混称。",
    "",
    "## 输出",
    "",
    "- reanalysis/results/tables/lactylome_and_reference_proteome_pairing_zh.csv",
    "- reanalysis/results/tables/lactylome_and_reference_proteome_gaps_zh.csv",
    "- reanalysis/results/tables/lactylome_reference_all_gaps_zh.csv",
    "- reanalysis/results/tables/lactylome_group_two_reference_columns_complete_zh.csv",
    "- reanalysis/results/tables/lactylome_dataset_decisions_zh.csv",
    "- reanalysis/results/tables/lactylome_reference_pairing_summary_zh.csv",
    "- reanalysis/results/tables/lactylome_and_reference_proteome_pairing_zh.xlsx"
  )
  writeLines(report, file.path(report_dir, "LACTYLOME_REFERENCE_PAIRING_STATUS.md"))

  message("Built lactylome/reference pairing tables and status report.")
}

# =============================================================================
# 段 summarize — 来源：summarize_acquired_lactylome_data.R
# =============================================================================
run_summarize_acquired_lactylome_data <- function() {
  table_dir <- file.path(project_root, "reanalysis", "results", "tables")
  report_dir <- file.path(project_root, "reanalysis", "reports")
  dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
  inventory <- read.csv(
    file.path(table_dir, "human_lactylome_mass_spectrometry_inventory.csv"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  count_maxquant_sites <- function(path) {
    sites <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
    valid <- sites$Reverse != "+" & sites$`Potential contaminant` != "+"
    valid[is.na(valid)] <- TRUE
    proteins <- unique(unlist(strsplit(sites$Proteins[valid], ";", fixed = TRUE)))
    sample_columns <- grep("^Localization prob ", names(sites), value = TRUE)
    list(
      site_rows = nrow(sites),
      valid_rows = sum(valid),
      proteins = length(proteins[nzchar(proteins)]),
      samples = paste(sub("^Localization prob ", "", sample_columns), collapse = "; ")
    )
  }

  lung <- count_maxquant_sites(file.path(
    project_root,
    "data/PXD036307/search_results/extracted/txt/La (K)Sites.txt"
  ))
  placenta <- count_maxquant_sites(file.path(
    project_root,
    "data/PXD063047/search_results/extracted/combined/txt/La (K)Sites.txt"
  ))

  prostate <- read.delim(
    file.path(
      project_root,
      "data/PXD066054/search_results/extracted/PLa/XB08700B1DPLa_-PTMSiteReport.tsv"
    ),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  prostate_lac <- str_detect(prostate$PTM.ModificationTitle, regex("lac", ignore_case = TRUE))

  hcc <- read.delim(
    file.path(
      project_root,
      "data/PXD075377/search_results/extracted/2-Basic_analysis/MS_identified_information.txt"
    ),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  summary <- tibble(
    PXD = c(
      "PXD036307", "PXD054919", "PXD063047",
      "PXD064912", "PXD066054", "PXD075377"
    ),
    `样本材料` = c(
      "正常人肺组织",
      "A549 肺腺癌细胞",
      "重度子痫前期胎盘与正常妊娠胎盘",
      "人精子",
      "前列腺癌组织与良性前列腺增生组织",
      "肝细胞癌组织与邻癌组织"
    ),
    `样本与实验组` = c(
      lung$samples,
      "A549 三个生物学重复",
      placenta$samples,
      "三个正常人精子样本",
      "NAT1-5（BPH）与 PCa1-5（前列腺癌）",
      "一例 HCC 与配对邻癌"
    ),
    `正常性分类` = c(
      "健康正常生理组织",
      "肿瘤细胞系",
      "含健康妊娠对照组织",
      "正常人源生物样本",
      "含良性病变对照，不能写成健康正常",
      "含邻癌对照，不能写成健康正常"
    ),
    `当前可用证据` = c(
      "MaxQuant La (K)Sites + evidence + proteinGroups",
      "论文补充表 MOESM2；仓库 Results.zip 无效",
      "MaxQuant La (K)Sites + evidence + proteinGroups",
      "Spectronaut .sne",
      "Spectronaut PTMSiteReport + Identification",
      "MS_identified_information 逐位点表"
    ),
    `原始位点或PTM记录数` = c(
      lung$site_rows,
      3110,
      placenta$site_rows,
      NA_integer_,
      sum(prostate_lac),
      nrow(hcc)
    ),
    `基础QC后位点记录数` = c(
      lung$valid_rows,
      3110,
      placenta$valid_rows,
      NA_integer_,
      NA_integer_,
      nrow(hcc)
    ),
    `证据表中的蛋白数` = c(
      lung$proteins,
      1220,
      placenta$proteins,
      NA_integer_,
      length(unique(prostate$PTM.ProteinId[prostate_lac])),
      length(unique(hcc$`Protein accession`))
    ),
    `下载与校验状态` = c(
      "两个检索压缩包完整且内容相同；已解压一份",
      "补充表可用；Results.zip 校验值匹配但内容全零",
      "检索压缩包可列出并解压；第二仓库副本已确认重复",
      "仓库 SHA-1 匹配",
      "两套结果 ZIP 完整且仓库 SHA-1 匹配",
      "结果 ZIP 完整"
    ),
    `主要本地路径` = c(
      "data/PXD036307/search_results/extracted/txt",
      "data/PXD054919/supplementary/41419_2025_8113_MOESM2_ESM.xlsx",
      "data/PXD063047/search_results/extracted/combined/txt",
      "data/PXD064912/search_results/P_0_HumanSperm.sne",
      "data/PXD066054/search_results/extracted/PLa",
      "data/PXD075377/search_results/extracted/2-Basic_analysis/MS_identified_information.txt"
    )
  ) |>
    left_join(
      inventory |>
        select(PXD, PublicationYear, RepositoryAnnounceYear, DOI, DatasetURL),
      by = "PXD"
    ) |>
    relocate(PublicationYear, RepositoryAnnounceYear, DOI, .after = PXD) |>
    rename(
      `论文发表年份` = PublicationYear,
      `仓库发布年份` = RepositoryAnnounceYear,
      `论文DOI` = DOI,
      `数据集链接` = DatasetURL
    )

  write.csv(
    summary,
    file.path(table_dir, "priority_dataset_acquisition_qc_summary_zh.csv"),
    row.names = FALSE,
    fileEncoding = "UTF-8",
    na = ""
  )

  report <- c(
    "# 人源乳酸化质谱数据获取进展",
    "",
    paste0("更新日期：", format(Sys.Date(), "%Y-%m-%d")),
    "",
    "## 当前范围",
    "",
    "- ProteomeXchange 多关键词去重后共 92 个相关 PXD，其中 46 个列为全局乳酸化蛋白组或候选。",
    "- 当前阶段优先获取作者检索结果、位点表和论文补充表，不批量下载大型 raw，也暂不把新增数据并入旧 Venn/DDR。",
    "- PXD038880/PXD050906 继续保留 hold，不进入分析。",
    "",
    "## 已获取并完成基础校验的数据",
    "",
    paste0(
      "- ",
      summary$PXD,
      "：",
      summary$样本材料,
      "；",
      summary$当前可用证据,
      "；",
      summary$下载与校验状态
    ),
    "",
    "## 正常材料解释",
    "",
    "- PXD036307 是健康正常生理人肺组织，也是当前最明确的新增正常组织乳酸化数据。",
    "- PXD063047 含 3 例正常妊娠胎盘和 3 例重度子痫前期胎盘，可拆分健康对照和疾病组。",
    "- PXD064912 是三个正常人精子样本，属于正常人源生物样本，但不是组织。",
    "- PXD066054 的 BPH 是良性病变对照，不能写成健康正常组织。",
    "- PXD075377 的邻癌组织是疾病研究对照，不能写成健康正常组织。",
    "- PXD054919 是 A549 肺腺癌细胞乳酸化质谱，不是正常材料。",
    "",
    "## PXD054919 原稿与数据",
    "",
    "- 老师提供的 PDF 与项目归档论文 SHA-256 完全一致。",
    "- 仓库 Results.zip 与提交 SHA-1 一致，但文件全部为 0x00，无法解压，属于仓库提交质量问题。",
    "- 论文补充表 MOESM2 可用，包含 A549 三个重复、3,110 个 Kla 位点和 1,220 个唯一蛋白。",
    "",
    "## 主要交付",
    "",
    "- 人源乳酸化质谱总表：reanalysis/results/tables/human_lactylome_mass_spectrometry_inventory.xlsx",
    "- 92 个 PXD 的机器可读总表：reanalysis/results/tables/human_lactylome_mass_spectrometry_inventory.csv",
    "- 仓库文件与来源 URL：reanalysis/results/tables/human_lactylome_repository_file_manifest.csv",
    "- 当前已获取数据 QC：reanalysis/results/tables/priority_dataset_acquisition_qc_summary_zh.csv",
    "- 当前下载文件登记：reanalysis/results/tables/priority_dataset_acquisition_manifest.csv",
    ""
  )
  writeLines(
    report,
    file.path(report_dir, "CURRENT_LACTYLOME_ACQUISITION_STATUS.md"),
    useBytes = TRUE
  )

  message("Acquisition QC summary rows: ", nrow(summary))
}

# =============================================================================
# 段 download（默认不执行）— 来源：probe_lactylome_pair_files.R /
#   download_lactylome_pair_files.R / extract_lactylome_pair_archives.R /
#   register_additional_lactylome_pair_files.R / download_healthy_tissue_references.R
# =============================================================================

# 子段 1：probe_lactylome_pair_files.R
run_probe_lactylome_pair_files <- function() {
  manifest_path <- file.path(
    project_root,
    "reanalysis",
    "results",
    "tables",
    "human_lactylome_repository_file_manifest.csv"
  )
  output_path <- file.path(
    project_root,
    "reanalysis",
    "results",
    "tables",
    "lactylome_pair_remote_file_sizes.csv"
  )

  priority_pxd <- c(
    "PXD028737", "PXD033146", "PXD037371", "PXD037530", "PXD045967",
    "PXD046344", "PXD046800", "PXD047535", "PXD047673", "PXD048995",
    "PXD050147", "PXD052772", "PXD053029", "PXD055230", "PXD057709",
    "PXD058534", "PXD062720", "PXD063266", "PXD063945", "PXD064038",
    "PXD065104", "PXD065831", "PXD066351", "PXD068838", "PXD070007",
    "PXD070427", "PXD073311", "PXD075014", "PXD077426"
  )

  manifest <- read.csv(manifest_path, stringsAsFactors = FALSE, check.names = FALSE)
  selected <- manifest[
    manifest$PXD %in% priority_pxd &
      manifest$FileCategory != "Associated raw file URI" &
      nzchar(manifest$FileURL),
  ]

  selected$CanonicalURL <- sub("^ftp://", "https://", selected$FileURL)
  selected$ProbeStatus <- ""
  selected$RemoteSizeBytes <- NA_real_
  selected$RemoteSizeMiB <- NA_real_
  selected$LastModified <- ""

  probe_one <- function(url) {
    command <- c(
      "-sL", "--max-time", "90", "-r", "0-0", "-D", "-", "-o", "/dev/null",
      url
    )
    command[[length(command)]] <- shQuote(command[[length(command)]])
    headers <- tryCatch(
      system2("curl", command, stdout = TRUE, stderr = TRUE),
      error = function(e) paste("ERROR", conditionMessage(e))
    )
    status <- attr(headers, "status")
    content_range <- grep("^Content-Range:", headers, ignore.case = TRUE, value = TRUE)
    content_length <- grep("^Content-Length:", headers, ignore.case = TRUE, value = TRUE)
    last_modified <- grep("^Last-Modified:", headers, ignore.case = TRUE, value = TRUE)
    http_status <- grep("^HTTP/", headers, value = TRUE)

    size <- NA_real_
    if (length(content_range)) {
      size <- suppressWarnings(as.numeric(sub(".*/", "", tail(content_range, 1))))
    } else if (length(content_length)) {
      size <- suppressWarnings(as.numeric(sub("^[^:]+:[[:space:]]*", "", tail(content_length, 1))))
    }

    list(
      status = if (is.null(status) || status == 0) {
        paste(tail(http_status, 1), collapse = "")
      } else {
        paste0("curl_exit_", status)
      },
      size = size,
      last_modified = if (length(last_modified)) {
        sub("^[^:]+:[[:space:]]*", "", tail(last_modified, 1))
      } else {
        ""
      }
    )
  }

  message("Probing ", nrow(selected), " processed/metadata files...")
  for (i in seq_len(nrow(selected))) {
    result <- probe_one(selected$CanonicalURL[[i]])
    selected$ProbeStatus[[i]] <- result$status
    selected$RemoteSizeBytes[[i]] <- result$size
    selected$RemoteSizeMiB[[i]] <- if (is.na(result$size)) NA_real_ else result$size / 1024^2
    selected$LastModified[[i]] <- result$last_modified
    message(
      sprintf(
        "[%d/%d] %s %s %.1f MiB",
        i,
        nrow(selected),
        selected$PXD[[i]],
        selected$FileName[[i]],
        selected$RemoteSizeMiB[[i]]
      )
    )
  }

  selected <- selected[
    order(selected$PXD, selected$FileCategory, selected$FileName, selected$CanonicalURL),
  ]
  write.csv(selected, output_path, row.names = FALSE, na = "")
  message("Wrote ", output_path)
}

# 子段 2：download_lactylome_pair_files.R
run_download_lactylome_pair_files <- function() {
  size_path <- file.path(
    project_root,
    "reanalysis",
    "results",
    "tables",
    "lactylome_pair_remote_file_sizes.csv"
  )
  output_path <- file.path(
    project_root,
    "reanalysis",
    "results",
    "tables",
    "lactylome_pair_download_manifest.csv"
  )

  remote <- read.csv(size_path, stringsAsFactors = FALSE, check.names = FALSE)

  selected_files <- list(
    PXD028737 = c("txt.zip"),
    PXD033146 = c("search_result-HA119TPLa.zip", "search_result-HA119TQ.zip"),
    PXD037371 = c("txt.rar"),
    PXD037530 = c("checksum.txt"),
    PXD045967 = c("8.16-mqpar.xml"),
    PXD046800 = c(
      "HFX2_LFQ_QB001_Lacty_PeptideGroups.txt",
      "HFX2_LFQ_QB001_Lacty_Proteins.txt",
      "HFX2_LFQ_QB002_PeptideGroups.txt",
      "HFX2_LFQ_QB002_Proteins.txt"
    ),
    PXD047535 = c("checksum.txt"),
    PXD047673 = c("checksum.txt"),
    PXD048995 = c("mqpar.xml"),
    PXD050147 = c(
      "checksum.txt",
      "Kla_evidence.txt",
      "Lactyl_K_Sites.txt",
      "SIRT_proteinGroups.txt"
    ),
    PXD052772 = c("mqpar.xml", "txt.zip"),
    PXD055230 = c("checksum.txt"),
    PXD057709 = c("checksum.txt"),
    PXD058534 = c("checksum.txt", "txt.zip"),
    PXD062720 = c("checksum.txt", "txt.zip"),
    PXD063266 = c("checksum.txt", "LactylSites.xlsx", "mqpar.xml"),
    PXD063945 = c("lac_mqpar.xml", "pro_mqpar.xml"),
    PXD064038 = c("Clinical information of samples.docx", "txt.zip"),
    PXD065831 = c("YAS202408210011-1.rar"),
    PXD066351 = c(
      "XB01472B1DA-DIA_result.tsv",
      "XB01472B1DA-MSstats_Input.tsv",
      "XB01472B1DA-Protein_Quant.tsv",
      "XB01472B1DPAc-DIA_result.csv",
      "XB01472B1DPAc-MSstats_Input.csv",
      "XB01472B1DPLa-DIA_result.csv",
      "XB01472B1DPLa-MSstats_Input.csv"
    ),
    PXD068838 = c("checksum.txt"),
    PXD070007 = c("SA206LPLaB1_Annotation.xlsx"),
    PXD070427 = c("checksum.txt"),
    PXD073311 = c("Database_search_result.zip")
  )

  keep <- vapply(seq_len(nrow(remote)), function(i) {
    pxd <- remote$PXD[[i]]
    pxd %in% names(selected_files) && remote$FileName[[i]] %in% selected_files[[pxd]]
  }, logical(1))
  queue <- remote[keep,]

  external <- data.frame(
    PXD = c("PXD066517", "PXD066517"),
    FileCategory = c("Search engine output file URI", "Other type file URI"),
    FileName = c("20240275.tsv", "Supplementary_Table.xlsx"),
    FileURL = c(
      "ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2025/07/PXD066517/20240275.tsv",
      "ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2025/07/PXD066517/Supplementary_Table.xlsx"
    ),
    LocalMatches = "",
    LocalStatus = "not_downloaded",
    DownloadPriority = "reference_proteome_first",
    RepositoryFileNote = "Normal human sperm conventional proteome reference",
    CanonicalURL = c(
      "https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/07/PXD066517/20240275.tsv",
      "https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/07/PXD066517/Supplementary_Table.xlsx"
    ),
    ProbeStatus = "repository_api_confirmed",
    RemoteSizeBytes = c(5560825, 12383625),
    RemoteSizeMiB = c(5560825, 12383625) / 1024^2,
    LastModified = "",
    stringsAsFactors = FALSE
  )
  queue <- rbind(queue, external[, names(queue)])

  download_one <- function(url, destination) {
    dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
    url <- gsub(" ", "%20", url, fixed = TRUE)
    command <- paste(
      "curl --fail --location --continue-at - --retry 3 --retry-delay 2 --silent --show-error",
      "--output", shQuote(destination),
      shQuote(url)
    )
    status <- system(command)
    if (status != 0) {
      stop("Download failed: ", url)
    }
  }

  records <- vector("list", nrow(queue))
  for (i in seq_len(nrow(queue))) {
    pxd <- queue$PXD[[i]]
    role_dir <- if (grepl(
      "checksum|mqpar|Clinical information",
      queue$FileName[[i]],
      ignore.case = TRUE
    )) {
      "metadata"
    } else {
      "search_results"
    }
    destination_name <- queue$FileName[[i]]
    if (
      pxd == "PXD073311" &&
        sum(queue$PXD == pxd & queue$FileName == destination_name) > 1
    ) {
      component <- sub(".*/(IPX[0-9]+)/.*", "\\1", queue$CanonicalURL[[i]])
      destination_name <- paste0(component, "_", destination_name)
    }
    destination <- file.path(project_root, "data", pxd, role_dir, destination_name)
    preexisting <- file.exists(destination) &&
      file.info(destination)$size > 0 &&
      (
        is.na(queue$RemoteSizeBytes[[i]]) ||
          file.info(destination)$size == queue$RemoteSizeBytes[[i]]
      )
    if (!preexisting) {
      message(sprintf("[%d/%d] downloading %s %s", i, nrow(queue), pxd, destination_name))
      download_one(queue$CanonicalURL[[i]], destination)
    } else {
      message(sprintf("[%d/%d] present %s %s", i, nrow(queue), pxd, destination_name))
    }
    records[[i]] <- data.frame(
      PXD = pxd,
      FileName = destination_name,
      FileRole = role_dir,
      LocalRelativePath = substring(destination, nchar(project_root) + 2),
      SizeBytes = file.info(destination)$size,
      SHA256 = sha256_file(destination),
      SourceURL = queue$CanonicalURL[[i]],
      RemoteSizeBytes = queue$RemoteSizeBytes[[i]],
      SizeMatchesProbe = is.na(queue$RemoteSizeBytes[[i]]) ||
        file.info(destination)$size == queue$RemoteSizeBytes[[i]],
      DownloadStatus = if (preexisting) "already_present" else "downloaded",
      stringsAsFactors = FALSE
    )
  }

  result <- do.call(rbind, records)
  write.csv(result, output_path, row.names = FALSE, na = "")

  for (pxd in unique(result$PXD)) {
    per_pxd <- result[result$PXD == pxd,]
    per_pxd_path <- file.path(project_root, "data", pxd, "metadata", "pairing_download_manifest.csv")
    dir.create(dirname(per_pxd_path), recursive = TRUE, showWarnings = FALSE)
    write.csv(per_pxd, per_pxd_path, row.names = FALSE, na = "")
  }

  message("Wrote ", output_path)
}

# 子段 3：extract_lactylome_pair_archives.R
run_extract_lactylome_pair_archives <- function() {
  manifest_path <- file.path(
    project_root,
    "reanalysis",
    "results",
    "tables",
    "lactylome_pair_download_manifest.csv"
  )
  output_path <- file.path(
    project_root,
    "reanalysis",
    "results",
    "tables",
    "lactylome_pair_extraction_manifest.csv"
  )

  manifest <- read.csv(manifest_path, stringsAsFactors = FALSE, check.names = FALSE)
  archives <- manifest[
    grepl("\\.(zip|rar)$", manifest$FileName, ignore.case = TRUE) &
      file.exists(file.path(project_root, manifest$LocalRelativePath)),
  ]

  records <- list()
  for (i in seq_len(nrow(archives))) {
    archive <- file.path(project_root, archives$LocalRelativePath[[i]])
    output_dir <- file.path(
      project_root,
      "data",
      archives$PXD[[i]],
      "search_results",
      "extracted_pairing",
      sub("\\.(zip|rar)$", "", archives$FileName[[i]], ignore.case = TRUE)
    )
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    message(sprintf("[%d/%d] extracting %s", i, nrow(archives), archive))
    status <- if (grepl("\\.zip$", archive, ignore.case = TRUE)) {
      system2("unzip", c("-q", "-n", shQuote(archive), "-d", shQuote(output_dir)))
    } else {
      system2("tar", c("-xf", shQuote(archive), "-C", shQuote(output_dir)))
    }
    files <- list.files(output_dir, recursive = TRUE, full.names = TRUE)
    records[[i]] <- data.frame(
      PXD = archives$PXD[[i]],
      Archive = archives$FileName[[i]],
      OutputRelativePath = substring(output_dir, nchar(project_root) + 2),
      ExtractionStatus = if (status == 0) "archive_extracted" else paste0("extract_exit_", status),
      ExtractedFileCount = length(files),
      ExtractedBytes = if (length(files)) sum(file.info(files)$size, na.rm = TRUE) else 0,
      stringsAsFactors = FALSE
    )
  }

  result <- if (length(records)) {
    do.call(rbind, records)
  } else {
    data.frame(
      PXD = character(),
      Archive = character(),
      OutputRelativePath = character(),
      ExtractionStatus = character(),
      ExtractedFileCount = integer(),
      ExtractedBytes = numeric(),
      stringsAsFactors = FALSE
    )
  }
  write.csv(result, output_path, row.names = FALSE)
  message("Wrote ", output_path)
}

# 子段 4：register_additional_lactylome_pair_files.R
run_register_additional_lactylome_pair_files <- function() {
  config_path <- file.path(
    project_root,
    "reanalysis/config/additional_lactylome_pair_files.csv"
  )
  table_dir <- file.path(project_root, "reanalysis/results/tables")
  manifest_path <- file.path(table_dir, "lactylome_pair_download_manifest.csv")
  remote_path <- file.path(table_dir, "lactylome_pair_remote_file_sizes.csv")

  config <- read.csv(config_path, stringsAsFactors = FALSE, check.names = FALSE)
  remote <- read.csv(remote_path, stringsAsFactors = FALSE, check.names = FALSE)
  existing <- read.csv(manifest_path, stringsAsFactors = FALSE, check.names = FALSE)

  records <- lapply(seq_len(nrow(config)), function(i) {
    row <- config[i,]
    path <- file.path(project_root, row$LocalRelativePath)
    remote_match <- remote[
      remote$PXD == row$PXD &
        remote$FileName == row$FileName,
    ]
    remote_size <- if (nrow(remote_match)) {
      remote_match$RemoteSizeBytes[[1]]
    } else {
      NA_real_
    }
    size <- if (file.exists(path)) file.info(path)$size else 0
    complete <- file.exists(path) &&
      (is.na(remote_size) || size == remote_size)
    data.frame(
      PXD = row$PXD,
      FileName = row$FileName,
      FileRole = row$FileRole,
      LocalRelativePath = row$LocalRelativePath,
      SizeBytes = size,
      SHA256 = if (complete) sha256_file(path) else "",
      SourceURL = row$SourceURL,
      RemoteSizeBytes = remote_size,
      SizeMatchesProbe = complete,
      DownloadStatus = if (complete) "downloaded" else "incomplete_or_missing",
      stringsAsFactors = FALSE
    )
  })

  additional <- bind_rows(records)
  combined <- bind_rows(
    existing |>
      filter(
        !paste(PXD, FileName) %in%
          paste(additional$PXD, additional$FileName)
      ),
    additional
  ) |>
    arrange(PXD, FileName)

  write.csv(combined, manifest_path, row.names = FALSE, na = "")
  if (!all(additional$SizeMatchesProbe)) {
    stop("One or more additional files are incomplete or missing", call. = FALSE)
  }
  message("Registered ", nrow(additional), " additional lactylome pair files.")
}

# 子段 5：download_healthy_tissue_references.R
run_download_healthy_tissue_references <- function() {
  config_path <- file.path(
    project_root,
    "reanalysis",
    "config",
    "healthy_tissue_reference_files.csv"
  )
  table_dir <- file.path(project_root, "reanalysis", "results", "tables")
  config <- read.csv(config_path, stringsAsFactors = FALSE, check.names = FALSE)

  # 与 pairing 段同原因，保留本脚本自带 base_accession（lib 版本会
  # trimws/去 NX_ 前缀，改变计数口径）。
  base_accession <- function(values) {
    values <- sub("^.*\\|([^|]+)\\|.*$", "\\1", values)
    sub("-[0-9]+$", "", values)
  }

  count_proteins <- function(path) {
    data <- read.delim(
      path,
      check.names = FALSE,
      stringsAsFactors = FALSE,
      quote = "",
      comment.char = ""
    )
    keep <- rep(TRUE, nrow(data))
    if ("Reverse" %in% names(data)) {
      keep <- keep & data$Reverse != "+"
    }
    if ("Potential contaminant" %in% names(data)) {
      keep <- keep & data$`Potential contaminant` != "+"
    }
    if ("Only identified by site" %in% names(data)) {
      keep <- keep & data$`Only identified by site` != "+"
    }
    column <- intersect(c("Protein IDs", "Majority protein IDs"), names(data))[[1]]
    accessions <- unlist(strsplit(as.character(data[[column]][keep]), ";", fixed = TRUE))
    accessions <- trimws(accessions)
    accessions <- accessions[nzchar(accessions) & !is.na(accessions)]
    length(unique(base_accession(accessions)))
  }

  download_one <- function(i) {
    row <- config[i,]
    archive <- file.path(
      project_root,
      "data",
      row$PXD,
      "search_results",
      row$FileName
    )
    dir.create(dirname(archive), recursive = TRUE, showWarnings = FALSE)
    complete <- file.exists(archive) && file.info(archive)$size == row$SizeBytes
    if (!complete) {
      message(sprintf("[%d/%d] downloading %s", i, nrow(config), row$FileName))
      command <- paste(
        "curl --fail --location --continue-at - --retry 4 --retry-delay 2 --silent --show-error",
        "--output", shQuote(archive), shQuote(row$SourceURL)
      )
      status <- system(command)
      if (status != 0) {
        return(data.frame(
          TissueKey = row$TissueKey,
          DisplayName = row$DisplayName,
          PXD = row$PXD,
          FileName = row$FileName,
          SizeBytes = if (file.exists(archive)) file.info(archive)$size else 0,
          SHA256 = "",
          ProteinGroupsPath = "",
          ProteinCount = NA_integer_,
          Status = paste0("download_failed_exit_", status),
          SourceURL = row$SourceURL,
          stringsAsFactors = FALSE
        ))
      }
    }

    extract_dir <- file.path(
      project_root,
      "data",
      row$PXD,
      "search_results",
      "extracted_healthy_tissues",
      row$TissueKey
    )
    dir.create(extract_dir, recursive = TRUE, showWarnings = FALSE)
    existing <- list.files(
      extract_dir,
      pattern = "proteinGroups\\.txt$",
      recursive = TRUE,
      full.names = TRUE
    )
    if (!length(existing)) {
      status <- system2(
        "unzip",
        c("-j", "-n", shQuote(archive), "*proteinGroups.txt", "-d", shQuote(extract_dir))
      )
      if (status != 0) {
        return(data.frame(
          TissueKey = row$TissueKey,
          DisplayName = row$DisplayName,
          PXD = row$PXD,
          FileName = row$FileName,
          SizeBytes = file.info(archive)$size,
          SHA256 = sha256_file(archive),
          ProteinGroupsPath = "",
          ProteinCount = NA_integer_,
          Status = paste0("proteinGroups_extract_failed_exit_", status),
          SourceURL = row$SourceURL,
          stringsAsFactors = FALSE
        ))
      }
      existing <- list.files(
        extract_dir,
        pattern = "proteinGroups\\.txt$",
        recursive = TRUE,
        full.names = TRUE
      )
    }
    protein_path <- existing[[1]]
    data.frame(
      TissueKey = row$TissueKey,
      DisplayName = row$DisplayName,
      PXD = row$PXD,
      FileName = row$FileName,
      SizeBytes = file.info(archive)$size,
      SHA256 = sha256_file(archive),
      ProteinGroupsPath = substring(protein_path, nchar(project_root) + 2),
      ProteinCount = count_proteins(protein_path),
      Status = "downloaded_extracted_counted",
      SourceURL = row$SourceURL,
      stringsAsFactors = FALSE
    )
  }

  worker_count <- min(3L, nrow(config))
  records <- parallel::mclapply(seq_len(nrow(config)), download_one, mc.cores = worker_count)
  result <- bind_rows(records) |>
    arrange(TissueKey)
  write.csv(
    result,
    file.path(table_dir, "healthy_tissue_reference_acquisition_manifest.csv"),
    row.names = FALSE,
    na = ""
  )
  message("Healthy tissue reference acquisition complete.")
}

# download 段编排（仅显式 --stage download 时由 main 调用）
run_download_stage <- function() {
  run_probe_lactylome_pair_files()
  run_download_lactylome_pair_files()
  run_extract_lactylome_pair_archives()
  run_register_additional_lactylome_pair_files()
  run_download_healthy_tissue_references()
}

# =============================================================================
# main：计算段默认全跑；download 段默认跳过，仅显式 --stage download 执行
# =============================================================================
run_stage("ensembl_map", run_build_ensembl_uniprot_mapping())
run_stage("inventory", run_build_human_lactylome_inventory())
run_stage("manifests", {
  run_build_lactylome_acquisition_manifests()
  run_build_healthy_special_reference_manifest()
})
run_stage("pairing", run_build_lactylome_reference_pairing())
run_stage("summarize", run_summarize_acquired_lactylome_data())

if (identical(stage, "download")) {
  message("[stage] download")
  run_download_stage()
}
