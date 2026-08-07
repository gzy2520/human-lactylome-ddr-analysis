#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) normalizePath(args[[1]]) else normalizePath(".")
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
