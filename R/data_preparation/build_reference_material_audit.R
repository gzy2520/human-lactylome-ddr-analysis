#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args)) normalizePath(args[[1]]) else normalizePath(".")
config_dir <- file.path(project_root, "config")
table_dir <- file.path(project_root, "results", "tables")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

scope <- read.csv(
  file.path(
    table_dir,
    "kla_regulator_intensity_availability_audit.csv"
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE,
  na.strings = c("", "NA")
) |>
  filter(`定量可用` %in% c(TRUE, "TRUE", "True", 1, "1")) |>
  transmute(
    KlaPXD = PXD,
    SampleGroup = `样本组`,
    ScopeOrder = row_number()
  )
if (nrow(scope) != 37) {
  stop("Strict material audit requires exactly 37 quantifiable Kla groups")
}

pairing <- read.csv(
  file.path(config_dir, "lactylome_reference_pairing.csv"),
  check.names = FALSE,
  stringsAsFactors = FALSE,
  na.strings = c("", "NA")
) |>
  transmute(
    KlaPXD = LactylomePXD,
    SampleGroup,
    ReferencePXD,
    ReferenceFile = ReferenceEvidenceLocator,
    ReferenceSubset = ReferenceSampleSubset,
    ConfigStrictInclude = IncludeInStrictReferenceAnalysis,
    ReferenceMatchQuality = MatchQuality,
    PairingCaveat = Caveat
  )

review <- read.csv(
  file.path(config_dir, "strict_reference_material_identity_review.csv"),
  check.names = FALSE,
  stringsAsFactors = FALSE,
  na.strings = c("", "NA")
) |>
  rename(KlaPXD = PXD)

audit <- scope |>
  left_join(pairing, by = c("KlaPXD", "SampleGroup")) |>
  left_join(review, by = c("KlaPXD", "SampleGroup"))

if (nrow(audit) != 37 || anyDuplicated(audit[c("KlaPXD", "SampleGroup")])) {
  stop("Material-review table must map one-to-one to the 37 Kla groups")
}
review_columns <- c(
  "KlaMaterialDetail", "KlaGranularity",
  "MaterialIdentityMatch", "ExperimentalStateMatch", "Evidence"
)
if (any(vapply(audit[review_columns], function(values) any(is.na(values)), logical(1)))) {
  stop("Material-review table has missing required values")
}

configured <- audit$ConfigStrictInclude %in% c(TRUE, "TRUE", "True", 1, "1")
material_match <- audit$MaterialIdentityMatch %in%
  c(TRUE, "TRUE", "True", 1, "1")
if (!identical(configured, material_match)) {
  bad <- paste(
    audit$KlaPXD[configured != material_match],
    audit$SampleGroup[configured != material_match],
    sep = "__"
  )
  stop(
    "Strict-reference inclusion disagrees with material identity for: ",
    paste(bad, collapse = ", ")
  )
}

active <- audit |>
  filter(configured) |>
  mutate(
    SharedReferenceFile = duplicated(ReferenceFile) |
      duplicated(ReferenceFile, fromLast = TRUE)
  )
shared_lookup <- active |>
  filter(SharedReferenceFile) |>
  group_by(ReferenceFile) |>
  summarise(
    SharedReferenceRows = paste(
      paste(KlaPXD, SampleGroup, ReferenceSubset, sep = " / "),
      collapse = " | "
    ),
    .groups = "drop"
  )

audit <- audit |>
  left_join(
    active |>
      select(KlaPXD, SampleGroup, SharedReferenceFile),
    by = c("KlaPXD", "SampleGroup")
  ) |>
  left_join(shared_lookup, by = "ReferenceFile") |>
  mutate(
    SharedReferenceFile = coalesce(SharedReferenceFile, FALSE),
    SharedReferenceRows = coalesce(SharedReferenceRows, ""),
    ReuseClass = case_when(
      !SharedReferenceFile ~ "not_shared",
      ReferencePXD == "PXD072220" ~
        "same_exact_cell_line_baseline_reused_for_two_Kla_studies",
      ReferencePXD == "PXD030304" ~
        "cell_line_atlas_shared_file_distinct_named_rows",
      TRUE ~ "shared_file_distinct_material_subsets"
    ),
    DistinctSubsetVerified = case_when(
      !SharedReferenceFile ~ "not_applicable",
      ReferencePXD == "PXD072220" ~
        "same_HK2_baseline_intentionally_reused",
      ReferencePXD == "PXD046800" &
        SampleGroup == "hypertrophic scar" ~ "HSP1_HSP4_only",
      ReferencePXD == "PXD046800" &
        SampleGroup == "adjacent skin" ~ "NSP1_NSP4_only",
      ReferencePXD == "PXD065775" &
        SampleGroup == "HCC" ~ "CISs_sheet_only",
      ReferencePXD == "PXD065775" &
        SampleGroup == "adjacent liver" ~ "ANTs_sheet_only",
      ReferencePXD == "PXD066054" &
        SampleGroup == "BPH" ~ "NAT1_NAT5_only",
      ReferencePXD == "PXD066054" &
        SampleGroup == "prostate cancer" ~ "PCa1_PCa5_only",
      ReferencePXD == "PXD069969" &
        SampleGroup == "glioblastoma stem cells" ~
          "G2907_G3028_G3264_GSC23_MES28_RKI_only",
      ReferencePXD == "PXD069969" &
        SampleGroup == "neural stem cells" ~ "ENSA_HMP1_only",
      TRUE ~ paste0("named_subset:", ReferenceSubset)
    ),
    AnalysisDecision = ifelse(
      configured & material_match,
      "included_exact_material_identity",
      "excluded_no_exact_material_reference"
    ),
    DecisionBasis = ifelse(
      configured & material_match,
      "Biological material identity/granularity is exact; experimental-state limitations are reported separately",
      PairingCaveat
    )
  ) |>
  arrange(ScopeOrder) |>
  select(
    KlaPXD,
    SampleGroup,
    KlaMaterial = KlaMaterialDetail,
    KlaGranularity,
    ReferencePXD,
    ReferenceFile,
    ReferenceSubset,
    ReferenceMaterial = ReferenceMaterialDetail,
    ReferenceGranularity,
    SharedReferenceFile,
    ReuseClass,
    DistinctSubsetVerified,
    MaterialIdentityMatch,
    ExperimentalStateMatch,
    AnalysisDecision,
    ReferenceMatchQuality,
    Evidence,
    DecisionBasis,
    PairingCaveat,
    SharedReferenceRows
  )

write.csv(
  audit,
  file.path(table_dir, "strict_reference_material_identity_audit.csv"),
  row.names = FALSE,
  na = ""
)

audit_zh <- audit
names(audit_zh) <- c(
  "乳酸化PXD", "样本组", "Kla材料", "Kla取材粒度",
  "普通全蛋白PXD", "普通全蛋白文件", "实际读取子集",
  "参照材料", "参照取材粒度", "是否共用参照文件",
  "共用类型", "子集分离验证", "材料身份是否严格匹配",
  "实验状态匹配情况", "分析决定", "原匹配等级",
  "核对证据", "决定依据", "原注意事项", "共用文件内全部读取行"
)
write.csv(
  audit_zh,
  file.path(table_dir, "strict_reference_material_identity_audit_zh.csv"),
  row.names = FALSE,
  na = ""
)

message(
  "Wrote strict material-identity audit: ",
  sum(material_match), " included and ",
  sum(!material_match), " excluded."
)
