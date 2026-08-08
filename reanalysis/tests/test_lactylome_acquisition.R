#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) normalizePath(args[[1]]) else normalizePath(".")
table_dir <- file.path(project_root, "reanalysis", "results", "tables")

assert <- function(condition, message) {
  if (!isTRUE(condition)) {
    stop(message, call. = FALSE)
  }
}

inventory <- read.csv(
  file.path(table_dir, "human_lactylome_mass_spectrometry_inventory.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
assert(nrow(inventory) == 92, "Expected 92 unique PXD rows")
assert(!anyDuplicated(inventory$PXD), "Inventory contains duplicate PXD values")
assert(
  sum(inventory$MassSpecScope %in% c("global_lactylome", "global_lactylome_candidate")) == 46,
  "Expected 46 global lactylome datasets/candidates"
)
assert(
  all(inventory$PriorityTier[inventory$PXD %in% c("PXD038880", "PXD050906")] == "hold"),
  "PXD038880/PXD050906 must remain hold"
)

expected_evidence <- c(
  "data/PXD036307/search_results/extracted/txt/La (K)Sites.txt",
  "data/PXD054919/supplementary/41419_2025_8113_MOESM2_ESM.xlsx",
  "data/PXD063047/search_results/extracted/combined/txt/La (K)Sites.txt",
  "data/PXD064912/search_results/P_0_HumanSperm.sne",
  "data/PXD066054/search_results/extracted/PLa/XB08700B1DPLa_-PTMSiteReport.tsv",
  "data/PXD075377/search_results/extracted/2-Basic_analysis/MS_identified_information.txt"
)
evidence_paths <- file.path(project_root, expected_evidence)
assert(all(file.exists(evidence_paths)), "One or more key evidence files are missing")
assert(all(file.info(evidence_paths)$size > 0), "One or more key evidence files are empty")

selected_pxd <- c(
  "PXD036307", "PXD054919", "PXD063047",
  "PXD064912", "PXD066054", "PXD075377"
)
manifest_paths <- file.path(project_root, "data", selected_pxd, "metadata", "download_manifest.csv")
assert(all(file.exists(manifest_paths)), "One or more per-PXD download manifests are missing")

repository_manifest <- read.csv(
  file.path(table_dir, "human_lactylome_repository_file_manifest.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
duplicate_row <- repository_manifest[
  repository_manifest$PXD == "PXD063047" &
    grepl("IPX0011731002/combined.rar$", repository_manifest$FileURL),
]
assert(nrow(duplicate_row) == 1, "PXD063047 duplicate repository row is missing")
assert(
  duplicate_row$LocalStatus == "verified_repository_duplicate_not_redownloaded",
  "PXD063047 duplicate repository row is not documented correctly"
)

zero_archive <- file.path(project_root, "data/PXD054919/search_results/Results.zip")
bytes <- readBin(zero_archive, what = "raw", n = file.info(zero_archive)$size)
assert(length(bytes) > 0 && all(bytes == as.raw(0)), "PXD054919 Results.zip is no longer all-zero")

workbook <- file.path(table_dir, "human_lactylome_mass_spectrometry_inventory.xlsx")
assert(file.exists(workbook) && file.info(workbook)$size > 0, "Inventory workbook is missing")

decisions <- read.csv(
  file.path(project_root, "reanalysis/config/lactylome_dataset_decisions.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
candidate_pxd <- inventory$PXD[
  inventory$MassSpecScope %in% c("global_lactylome", "global_lactylome_candidate")
]
assert(nrow(decisions) == 46, "Expected 46 dataset decisions")
assert(!anyDuplicated(decisions$PXD), "Dataset decision table contains duplicate PXD rows")
assert(
  setequal(decisions$PXD, candidate_pxd),
  "Dataset decision table does not cover every lactylome candidate exactly once"
)
assert(
  all(decisions$DatasetClass[decisions$PXD %in% c("PXD038880", "PXD050906")] == "hold"),
  "PXD038880/PXD050906 must remain hold in the final decision table"
)
assert(
  decisions$DatasetClass[decisions$PXD == "PXD077426"] == "duplicate_mirror",
  "PXD077426 must remain a mirror of PXD078736"
)
assert(
  all(
    decisions$DatasetClass[decisions$PXD %in% c("PXD058173", "PXD065104")] ==
      "targeted_not_global"
  ),
  "Targeted or mechanism-only datasets must not be treated as global lactylomes"
)

pairing <- read.csv(
  file.path(table_dir, "lactylome_and_reference_proteome_pairing_zh.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
assert(nrow(pairing) >= 50, "Expected sample-group-level lactylome/reference pairing rows")
assert(
  all(nzchar(pairing$乳酸化数据判定)),
  "Every pairing row must have a lactylome status"
)
assert(
  all(nzchar(pairing$常规蛋白组获取状态)),
  "Every pairing row must have an explicit conventional-proteome status"
)
assert(
  all(
    pairing$当前已具备成对计数条件[
      pairing$乳酸化PXD %in% c("PXD038880", "PXD050906", "PXD077426")
    ] == FALSE
  ),
  "Hold and duplicate datasets must not be pair-ready"
)
assert(
  all(
    pairing$常规蛋白数[
      pairing$当前已具备成对计数条件
    ] > 0
  ),
  "Every pair-ready row must have a positive conventional-proteome count"
)
assert(
  all(
    pairing$乳酸化蛋白数[
      pairing$当前已具备成对计数条件
    ] > 0
  ),
  "Every pair-ready row must have a positive lactylome count"
)
assert(
  all(
    pairing$健康组织基线PXD[
      pairing$样本组 %in% c("normal human lung", "normal pregnancy placenta")
    ] == "PXD010154"
  ),
  "Normal lung and placenta must use the selected healthy tissue atlas baseline"
)
assert(
  all(
    nzchar(pairing$健康组织基线PXD[pairing$配置要求进入成对分析])
  ),
  "Every analysis-required lactylome group must have a healthy-tissue baseline"
)
assert(
  all(
    pairing$健康组织蛋白数[pairing$配置要求进入成对分析] > 0
  ),
  "Every analysis-required lactylome group must have a positive healthy-tissue count"
)
assert(
  all(
    pairing$健康组织基线PXD[
      pairing$样本组 %in% c("MCF7", "MCF10A", "MDA-MB-468", "T-47D")
    ] == "PXD016999"
  ),
  "Breast-derived cell lines must use the selected healthy breast tissue baseline"
)
assert(
  pairing$常规蛋白组PXD[pairing$乳酸化PXD == "PXD073311"] == "PXD073311" &&
    pairing$常规蛋白数[pairing$乳酸化PXD == "PXD073311"] == 7794 &&
    pairing$常规蛋白组样本子集[pairing$乳酸化PXD == "PXD073311"] ==
      "A0h_1;A0h_2;A0h_3" &&
    grepl(
      "report.pg_matrix.tsv$",
      pairing$常规蛋白组证据文件[pairing$乳酸化PXD == "PXD073311"]
    ),
  paste(
    "PXD073311 must use its same-study A0h ordinary whole-proteome",
    "matrix and exclude A6h"
  )
)
assert(
  pairing$乳酸化蛋白数[pairing$乳酸化PXD == "PXD064912"] == 231 &&
    pairing$乳酸化证据文件[pairing$乳酸化PXD == "PXD064912"] ==
      "data/PXD064912/supplementary/europepmc/mmc1.xlsx",
  "PXD064912 must use the open-access author supplementary Kla site table"
)
assert(
  all(
    pairing$乳酸化获取状态[pairing$乳酸化PXD == "PXD065831"] ==
      "远程可用，因体量暂未下载"
  ) &&
    all(
      pairing$常规蛋白组获取状态[pairing$乳酸化PXD == "PXD065831"] ==
        "已下载并完成蛋白计数"
    ),
  "PXD065831 lactylome and ordinary proteome components must not be swapped"
)
assert(
  pairing$乳酸化蛋白数[
    pairing$乳酸化PXD == "PXD070007" &
      pairing$样本组 == "glioblastoma stem cells"
  ] == 2564 &&
    pairing$乳酸化蛋白数[
      pairing$乳酸化PXD == "PXD070007" &
        pairing$样本组 == "neural stem cells"
    ] == 1527,
  "PXD070007 GSC/NSC lactylome counts are incorrect"
)

acquired_global <- pairing[
  pairing$乳酸化数据判定 %in%
    c("可用的全局乳酸化数据", "同研究可用乳酸化组件") &
    !is.na(pairing$乳酸化蛋白数),
]
strict_acquired <- acquired_global[
  acquired_global$配置要求进入严格参照分析,
]
strict_excluded <- acquired_global[
  !acquired_global$配置要求进入严格参照分析,
]
assert(
  nrow(acquired_global) > 0 &&
    all(strict_acquired$常规蛋白数 > 0) &&
    all(
      is.na(strict_excluded$常规蛋白数) |
        strict_excluded$常规蛋白数 <= 0
    ) &&
    all(acquired_global$健康组织蛋白数 > 0),
  paste(
    "Strictly paired groups must have countable ordinary proteomes;",
    "excluded groups must not retain a surrogate ordinary proteome;"
    ,"all groups must retain a countable healthy-tissue baseline"
  )
)
strict_excluded_groups <- c(
  "normal liver",
  "nonmetastatic HCC",
  "lung-metastatic HCC",
  "bladder cancer cells treated with EPI",
  "severe preeclampsia placenta",
  "MEC and NEC ESCC groups",
  "AC16 control and hypoxia"
)
assert(
  all(
    pairing$匹配质量[pairing$样本组 %in% strict_excluded_groups] ==
      "未找到精确参考"
  ),
  "Strictly excluded groups must remain explicitly labelled as lacking an exact reference"
)
assert(
  all(
    pairing$匹配质量[
      pairing$样本组 %in%
        c("HMC3", "glioblastoma stem cells", "neural stem cells", "HCC", "adjacent liver")
    ] %in% c(
      "同研究同样本精确匹配",
      "同一生物样本精确匹配",
      "疾病组织精确匹配",
      "邻近组织精确匹配"
    )
  ),
  "Resolved exact references must not remain labelled as healthy-organ surrogates"
)

compact_pairing <- file.path(
  table_dir,
  "lactylome_group_two_reference_columns_complete_zh.csv"
)
assert(
  file.exists(compact_pairing) && file.info(compact_pairing)$size > 0,
  "Compact two-reference-column table is missing"
)

healthy_config <- read.csv(
  file.path(project_root, "reanalysis/config/healthy_tissue_reference_files.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
assert(nrow(healthy_config) == 12, "Expected 12 healthy-tissue reference archives")
assert(!anyDuplicated(healthy_config$TissueKey), "Healthy-tissue keys must be unique")
assert(all(healthy_config$SizeBytes > 0), "Healthy-tissue archive sizes must be recorded")

healthy_special <- read.csv(
  file.path(
    project_root,
    "reanalysis/config/healthy_special_reference_catalog.csv"
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
assert(nrow(healthy_special) == 7, "Expected seven supplemental healthy references")
assert(
  all(healthy_special$ProteinCount > 0),
  "Supplemental healthy references must have positive protein counts"
)
healthy_special_manifest <- read.csv(
  file.path(
    table_dir,
    "healthy_special_reference_acquisition_manifest.csv"
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
assert(
  all(healthy_special_manifest$Status == "downloaded_verified_counted"),
  "Every supplemental healthy reference must be downloaded and verified"
)

cat("All lactylome acquisition tests passed.\n")
