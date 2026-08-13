#!/usr/bin/env Rscript

# Publication workflow for the configured human Kla proteomics analysis.
#
# Usage:
#   Rscript workflow/run_pipeline.R [all|core|embeddings|figures|selected_figures|validate]
#
# The default is "all". Set KLA_DATA_ROOT only when the raw data directory is
# stored outside the repository. All stochastic steps use seed 25 internally.

args <- commandArgs(trailingOnly = TRUE)
target <- if (length(args)) args[[1L]] else "selected_figures"
valid_targets <- c(
  "all", "core", "embeddings", "figures", "selected_figures", "validate"
)
if (!target %in% valid_targets) {
  stop("Target must be one of: ", paste(valid_targets, collapse = ", "))
}

cmd <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd, value = TRUE)
script_path <- if (length(file_arg)) {
  normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE)
} else {
  normalizePath("workflow/run_pipeline.R", mustWork = TRUE)
}
project_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)

data_root <- Sys.getenv("KLA_DATA_ROOT", unset = file.path(project_root, "data"))
if (!dir.exists(data_root)) stop("Raw-data directory does not exist: ", data_root)
if (normalizePath(data_root) != normalizePath(file.path(project_root, "data"))) {
  stop(
    "External KLA_DATA_ROOT is recorded but automatic path remapping is not yet ",
    "supported by dataset-specific extractors. Link it to <project>/data first."
  )
}

run_r_script <- function(relative_path, label) {
  path <- file.path(project_root, relative_path)
  message("\n[", label, "] ", relative_path)
  status <- system2("Rscript", c(path, project_root))
  if (!identical(status, 0L)) {
    stop("Workflow step failed: ", relative_path, call. = FALSE)
  }
}

run_python_script <- function(relative_path, label) {
  path <- file.path(project_root, relative_path)
  message("\n[", label, "] ", relative_path)
  status <- system2(
    "python3",
    c(path, "--project-root", project_root)
  )
  if (!identical(status, 0L)) {
    stop("Workflow step failed: ", relative_path, call. = FALSE)
  }
}

core_steps <- list(
  c("R/data_preparation/build_kla_regulator_landscape.R", "02 expanded group catalog"),
  c("R/analysis/analyze_regulator_reference_intensity.R", "03 reference intensity"),
  c("R/analysis/analyze_regulator_kla_intensity.R", "04 Kla intensity"),
  c("R/data_preparation/build_reference_material_audit.R", "05 reference audit"),
  c("R/analysis/analyze_ddr_fraction.R", "06 DDR fraction"),
  c("R/figures/plot_four_class_venn.R", "07 four-class Venn"),
  c("R/data_preparation/build_teacher_review_table.R", "08 review table"),
  c("R/data_preparation/prepare_protein_function_inputs.R", "09 protein function inputs"),
  c("R/figures/plot_bp_semantic_umap.R", "10 BP semantic UMAP")
)

embedding_steps <- list(
  c("R/analysis/tune_five_set_embeddings.R", "11 tune five-set embeddings"),
  c("R/figures/plot_five_set_embeddings.R", "12 five-set embeddings")
)

figure_steps <- list(
  c("R/figures/plot_pathway_specific_umap.R", "13 pathway-specific UMAP"),
  c("R/figures/plot_five_set_pathway_matrix.R", "14 five-set pathway matrix"),
  c("R/analysis/summarize_four_class_venn_counts.R", "15 four-Venn set-count summary")
)

selected_figure_steps <- list(
  c("R/figures/plot_five_set_pathway_matrix.R", "11 selected 4+1 pathway matrices"),
  c("R/analysis/summarize_four_class_venn_counts.R", "12 four-Venn set-count summary")
)

if (target %in% c("all", "core")) {
  run_python_script(
    "python/data_preparation/build_core_kla_inputs.py",
    "01 core Kla evidence"
  )
  for (step in core_steps) run_r_script(step[[1L]], step[[2L]])
}
if (target %in% c("all", "embeddings")) {
  for (step in embedding_steps) run_r_script(step[[1L]], step[[2L]])
}
if (target %in% c("all", "figures")) {
  for (step in figure_steps) run_r_script(step[[1L]], step[[2L]])
}
if (target == "selected_figures") {
  run_python_script(
    "python/data_preparation/build_core_kla_inputs.py",
    "01 core Kla evidence"
  )
  run_python_script(
    "python/data_preparation/build_reference_proteome_membership.py",
    "02 reference-proteome membership"
  )
  for (step in core_steps[seq_len(8L)]) {
    run_r_script(step[[1L]], step[[2L]])
  }
  for (step in selected_figure_steps) {
    run_r_script(step[[1L]], step[[2L]])
  }
}
if (target %in% c("all", "selected_figures", "validate")) {
  run_r_script("tests/validate_publication_contract.R", "16 publication contract")
}

message("\nWorkflow target completed: ", target)
