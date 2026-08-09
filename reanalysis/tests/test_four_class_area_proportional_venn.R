#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) normalizePath(args[[1]]) else normalizePath(".")
figure_root <- file.path(
  project_root,
  "reanalysis",
  "results",
  "figures",
  "four_class_venn"
)
table_root <- file.path(
  project_root,
  "reanalysis",
  "results",
  "tables",
  "four_class_venn"
)

assert <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

analyses <- c(
  "all_kla_four_class_venn",
  "kla_ddr_four_class_venn",
  "reference_proteome_four_class_venn",
  "reference_proteome_ddr_four_class_venn"
)
categories <- c(
  "normal_tissue",
  "cancer_tissue",
  "normal_cells",
  "cancer_cells"
)
flags <- paste0("In_", categories)

memberships <- list()
for (analysis in analyses) {
  output_dir <- file.path(table_root, analysis)
  required_tables <- file.path(
    output_dir,
    c(
      "set_counts.csv",
      "region_counts.csv",
      "membership.csv",
      "euler_fit_regions.csv",
      "euler_fit_summary.csv"
    )
  )
  assert(all(file.exists(required_tables)), paste("Missing tables for", analysis))

  set_counts <- read.csv(
    file.path(output_dir, "set_counts.csv"),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  region_counts <- read.csv(
    file.path(output_dir, "region_counts.csv"),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  membership <- read.csv(
    file.path(output_dir, "membership.csv"),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  annotation <- read.csv(
    file.path(output_dir, "id_annotation_mapping.csv"),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  fit <- read.csv(
    file.path(output_dir, "euler_fit_summary.csv"),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  assert(
    identical(set_counts$Category, categories),
    paste("Invalid category order for", analysis)
  )
  assert(
    identical(
      names(membership),
      c(
        "BaseAccession",
        "GeneSymbolAudit",
        "ProteinNameAudit",
        "ReviewedStatus",
        "GeneSymbolAuditSource",
        "AnnotationMappingSource",
        flags,
        "Region"
      )
    ),
    paste("Membership table must contain ID, audit mapping, and flags for", analysis)
  )
  assert(
    identical(
      names(annotation),
      c(
        "BaseAccession",
        "GeneSymbolAudit",
        "ProteinNameAudit",
        "ReviewedStatus",
        "GeneSymbolAuditSource",
        "AnnotationMappingSource"
      )
    ),
    paste("ID annotation table is malformed for", analysis)
  )
  assert(
    setequal(annotation$BaseAccession, membership$BaseAccession),
    paste("ID annotation table does not cover membership IDs for", analysis)
  )
  assert(
    !anyDuplicated(membership$BaseAccession),
    paste("Duplicate BaseAccession in", analysis)
  )
  assert(
    all(!is.na(membership$GeneSymbolAudit) & nzchar(membership$GeneSymbolAudit)),
    paste("Audit symbols must not be blank for", analysis)
  )
  assert(
    all(!is.na(membership$ProteinNameAudit) & nzchar(membership$ProteinNameAudit)),
    paste("Audit protein names must not be blank for", analysis)
  )
  assert(
    all(
      membership$ReviewedStatus %in%
        c("reviewed", "non-reviewed", "source-only", "unmapped")
    ),
    paste("Unexpected review status for", analysis)
  )
  assert(
    identical(membership$BaseAccession, sort(membership$BaseAccession)),
    paste("Membership IDs must remain deterministically sorted for", analysis)
  )
  assert(
    sum(region_counts$ProteinCount) == nrow(membership),
    paste("Region counts do not reconstruct membership for", analysis)
  )
  for (i in seq_along(categories)) {
    assert(
      sum(membership[[flags[[i]]]]) == set_counts$ProteinCount[[i]],
      paste("Set count mismatch for", analysis, categories[[i]])
    )
  }
  assert(
    is.finite(fit$Stress[[1]]) && fit$Stress[[1]] < 0.02,
    paste("Euler fit stress is too high for", analysis)
  )
  memberships[[analysis]] <- membership

  for (language in c("zh", "en")) {
    for (extension in c("png", "pdf")) {
      path <- file.path(
        figure_root,
        paste0(analysis, "_", language, ".", extension)
      )
      assert(file.exists(path), paste("Missing figure", path))
      assert(file.info(path)$size > 10000, paste("Figure is unexpectedly small", path))
    }
  }
}

all_kla_ids <- memberships[["all_kla_four_class_venn"]]$BaseAccession
kla_ddr_ids <- memberships[["kla_ddr_four_class_venn"]]$BaseAccession
reference_ids <- memberships[["reference_proteome_four_class_venn"]]$BaseAccession
reference_ddr_ids <-
  memberships[["reference_proteome_ddr_four_class_venn"]]$BaseAccession
assert(
  all(kla_ddr_ids %in% all_kla_ids),
  "Kla DDR IDs must be a subset of all Kla IDs"
)
assert(
  all(reference_ddr_ids %in% reference_ids),
  "Reference DDR IDs must be a subset of reference proteome IDs"
)

scope <- read.csv(
  file.path(table_root, "venn_sample_group_scope.csv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
assert(nrow(scope) == 37, "Venn scope audit must retain the 37-group source audit")
assert(sum(scope$KlaIncludedInVenn) == 33, "Kla Venn must use 33 strict-reference groups")
assert(sum(scope$ReferenceIncludedInVenn) == 33, "Reference Venn must use 33 groups")
excluded_keys <- c(
  "PXD062720__bladder cancer cells treated with EPI",
  "PXD063047__severe preeclampsia placenta",
  "PXD064038__MEC and NEC ESCC groups",
  "PXD075014__AC16 control and hypoxia"
)
assert(
  setequal(scope$SampleGroupKey[!scope$KlaIncludedInVenn], excluded_keys),
  "Venn scope must exclude exactly the four no-reference Kla groups"
)
assert(
  !any(scope$KlaIncludedInVenn & scope$SampleGroupKey %in% excluded_keys),
  "Excluded Kla groups entered the Venn scope"
)
paired_stats <- read.csv(
  file.path(
    project_root,
    "reanalysis",
    "results",
    "tables",
    "cell_type_kla_vs_reference_ddr_statistics_accession_only_paired_33.csv"
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
paired_stats_zh <- read.csv(
  file.path(
    project_root,
    "reanalysis",
    "results",
    "tables",
    "cell_type_kla_vs_reference_ddr_statistics_accession_only_paired_33_zh.csv"
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
assert(nrow(paired_stats) == 33, "Strict paired English table must contain 33 groups")
assert(nrow(paired_stats_zh) == 33, "Strict paired Chinese table must contain 33 groups")
assert(
  !any(paste(paired_stats$PXD, paired_stats$SampleGroup, sep = "__") %in% excluded_keys),
  "Strict paired table contains an excluded Kla group"
)

archive_root <- file.path(
  project_root,
  "archive",
  "reanalysis_2026-08-07_pre_37group_reference_fix",
  "results",
  "tables",
  "four_class_venn"
)
for (analysis in c(
  "all_kla_four_class_venn",
  "kla_ddr_four_class_venn"
)) {
  archived_path <- file.path(archive_root, analysis, "membership.csv")
  assert(file.exists(archived_path), paste("Missing archived baseline for", analysis))
  archived <- read.csv(
    archived_path,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  current <- memberships[[analysis]]
  assert(
    !length(setdiff(current$BaseAccession, archived$BaseAccession)),
    paste("Current Kla Venn added proteins outside the archived baseline for", analysis)
  )
}

stats <- read.csv(
  file.path(
    project_root,
    "reanalysis",
    "results",
    "tables",
    "cell_type_kla_vs_reference_ddr_statistics_accession_only.csv"
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
assert(nrow(stats) == 37, "Four-class Venn source audit must contain 37 groups")
assert(
  identical(
    as.integer(table(factor(stats$Category, levels = categories))),
    c(10L, 3L, 10L, 14L)
  ),
  "Four-class Venn input category counts must be 10/3/10/14"
)
assert(
  sum(stats$PairedAnalysisIncluded) == 33,
  "Reference-proteome Venn input must contain 33 exact-reference groups"
)
reference_by_group <- read.csv(
  file.path(
    project_root,
    "reanalysis",
    "intermediate",
    "expanded_ddr_by_accession",
    "reference_proteins_by_sample_group.csv"
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
assert(
  nrow(unique(reference_by_group[c("PXD", "SampleGroup")])) == 33 &&
    !any(reference_by_group$PXD == "PXD062720"),
  "Reference Venn membership must use only the 33 exact-reference groups"
)

message("Four-class area-proportional Venn tests passed.")
