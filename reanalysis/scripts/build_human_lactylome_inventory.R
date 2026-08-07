#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(jsonlite)
  library(stringr)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) normalizePath(args[[1]]) else normalizePath(".")
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
