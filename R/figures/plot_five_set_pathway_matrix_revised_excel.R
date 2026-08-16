#!/usr/bin/env Rscript

# Generate an exploratory 4+1 signed pathway-matrix preview from the supervisor's
# revised 2026-08-16 workbook. The direct-GO-term publication figures remain
# unchanged and are written to a different output directory.

full_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", full_args, value = TRUE)
script_path <- if (length(file_arg) > 0L) {
  normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE)
} else {
  normalizePath(
    "R/figures/plot_five_set_pathway_matrix_revised_excel.R",
    mustWork = TRUE
  )
}

Sys.setenv(
  KLA_SCORE_ANALYSIS_NAME =
    "five_set_pathway_matrix_revised_excel_20260816",
  KLA_SCORE_WORKBOOK_PATH =
    "data/identifier/乳酸化DDR基因评分表_Revised_20260816.xlsx",
  KLA_SCORE_PATHWAY_DISPLAY =
    "config/seven_pathway_display.csv"
)

source(
  file.path(dirname(script_path), "plot_five_set_pathway_matrix.R"),
  chdir = FALSE
)

previous_workbook_path <- file.path(
  project_root,
  "data/identifier/260810乳酸化DDR基因评分表.xlsx"
)
assert(
  file.exists(previous_workbook_path),
  "The previous score workbook is unavailable for comparison."
)

previous_raw <- as.data.table(
  read_excel(previous_workbook_path, sheet = "评分表")
)
previous_scores <- previous_raw[
  !is.na(BaseAccession) & nzchar(trimws(BaseAccession))
]
previous_scores[, BaseAccession := trimws(BaseAccession)]
assert(
  uniqueN(previous_scores$BaseAccession) == nrow(previous_scores) &&
    all(score_columns %in% names(previous_scores)) &&
    all(as.matrix(
      previous_scores[, score_columns, with = FALSE]
    ) %in% c(-1, 0, 1)),
  "The previous workbook does not contain a valid comparison matrix."
)
previous_scores[, PreviousSignedScore := as.numeric(
  as.matrix(previous_scores[, score_columns, with = FALSE]) %*% weights
)]
assert(
  !anyNA(previous_scores$score) &&
    max(abs(
      as.numeric(previous_scores$score) -
        previous_scores$PreviousSignedScore
    )) < 1e-12,
  "The previous workbook score does not match the weighted formula."
)

comparison <- merge(
  scores[
    ,
    c(
      "BaseAccession",
      "GeneSymbol",
      "ProteinName",
      "SignedScore",
      score_columns
    ),
    with = FALSE
  ],
  previous_scores[
    BaseAccession %in% membership$BaseAccession,
    c("BaseAccession", "PreviousSignedScore", score_columns),
    with = FALSE
  ],
  by = "BaseAccession",
  suffixes = c("_Revised", "_Previous"),
  all = TRUE,
  sort = FALSE
)
assert(
  nrow(comparison) == 399L && !anyNA(comparison$BaseAccession),
  "The revised and previous workbooks do not cover the same current 399 proteins."
)

changed_matrix <- sapply(
  score_columns,
  function(pathway) {
    comparison[[paste0(pathway, "_Revised")]] !=
      comparison[[paste0(pathway, "_Previous")]]
  }
)
comparison[, ChangedPathwayCount := rowSums(changed_matrix)]
comparison[, ChangedPathways := apply(
  changed_matrix,
  1L,
  function(changed) paste(score_columns[changed], collapse = ";")
)]
comparison[, `:=`(
  ScoreDelta = SignedScore - PreviousSignedScore,
  ScoreChanged = SignedScore != PreviousSignedScore
)]
setorder(comparison, BaseAccession)

comparison_summary <- rbindlist(lapply(score_columns, function(pathway) {
  revised <- comparison[[paste0(pathway, "_Revised")]]
  previous <- comparison[[paste0(pathway, "_Previous")]]
  data.table(
    Pathway = pathway,
    ChangedProteinCount = sum(revised != previous),
    PreviousSuppressingCount = sum(previous == -1L),
    PreviousUnassignedCount = sum(previous == 0L),
    PreviousPromotingCount = sum(previous == 1L),
    RevisedSuppressingCount = sum(revised == -1L),
    RevisedUnassignedCount = sum(revised == 0L),
    RevisedPromotingCount = sum(revised == 1L)
  )
}))
comparison_overview <- data.table(
  CurrentProteinCount = nrow(comparison),
  ProteinsWithAnyStateChange = sum(comparison$ChangedPathwayCount > 0L),
  ChangedProteinPathwayCells = sum(comparison$ChangedPathwayCount),
  ProteinsWithScoreChange = sum(comparison$ScoreChanged)
)
assert(
  comparison_overview$ProteinsWithAnyStateChange == 12L &&
    comparison_overview$ChangedProteinPathwayCells == 19L &&
    comparison_overview$ProteinsWithScoreChange == 12L,
  "The revised-versus-previous workbook comparison changed unexpectedly."
)

fwrite(
  comparison[ChangedPathwayCount > 0L],
  file.path(
    table_dir,
    "revised_vs_previous_score_changes_current399.csv"
  )
)
fwrite(
  comparison_summary,
  file.path(
    table_dir,
    "revised_vs_previous_pathway_summary_current399.csv"
  )
)
fwrite(
  comparison_overview,
  file.path(
    table_dir,
    "revised_vs_previous_comparison_overview.csv"
  )
)

input_audit_path <- file.path(table_dir, "input_file_audit.csv")
input_audit <- fread(input_audit_path)
input_audit <- rbind(
  input_audit,
  data.table(
    Input = "Previous DDR score workbook for comparison only",
    Path = sub(
      paste0("^", project_root, "/?"),
      "",
      previous_workbook_path
    ),
    SHA256 = digest::digest(
      file = previous_workbook_path,
      algo = "sha256",
      serialize = FALSE
    ),
    Role = "comparison baseline; not used to draw revised figures"
  ),
  fill = TRUE
)
fwrite(input_audit, input_audit_path)

message(
  "Compared revised and previous workbooks: 12 current proteins and ",
  "19 protein-pathway cells changed."
)
