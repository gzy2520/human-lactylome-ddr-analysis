#!/usr/bin/env Rscript

# Build an isolated correction package from the latest 31-group candidate input.
# The plotting code remains the reviewed candidate renderer; only the input
# contract is changed here: plotted observations are restricted to
# ObservationType == "sample". The original final package is copied intact and
# is never overwritten.

suppressPackageStartupMessages({
  library(data.table)
})

set.seed(25)

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1L) normalizePath(args[[1L]], mustWork = TRUE) else normalizePath(".", mustWork = TRUE)

stop_if <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

source_root <- normalizePath(Sys.getenv(
  "KLA_SOURCE_ROOT", unset = "/Users/gzy2520/Desktop/Research/kla"
), mustWork = TRUE)
original_final <- file.path(source_root, "results", "final_figures_and_tables")
expanded_input <- file.path(
  source_root, "data", "candidate",
  "escc_inclusion_20260903_pxd065830_tumor_reference", "candidate_input"
)
output_dir <- Sys.getenv(
  "KLA_CORRECTED_OUTPUT",
  unset = file.path(project_root, "corrected_final_result_20260905")
)

stop_if(dir.exists(original_final), paste0("Missing original final package: ", original_final))
stop_if(dir.exists(expanded_input), paste0("Missing latest expanded candidate input: ", expanded_input))
stop_if(!dir.exists(output_dir), paste0("Output already exists; choose a new KLA_CORRECTED_OUTPUT: ", output_dir))

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
log_dir <- file.path(output_dir, "render_logs")
filtered_dir <- file.path(output_dir, "Corrected_Data", "sample_only_inputs")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(filtered_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------------
# Copy the frozen package first. Only explicitly listed corrected figures and
# statistical sidecars are replaced below.
# -------------------------------------------------------------------------
copy_tree <- function(source, destination) {
  source_files <- list.files(source, recursive = TRUE, full.names = TRUE, all.files = FALSE)
  source_files <- source_files[file.info(source_files)$isdir %in% FALSE]
  for (source_file in source_files) {
    relative <- substring(source_file, nchar(source) + 2L)
    destination_file <- file.path(destination, relative)
    dir.create(dirname(destination_file), recursive = TRUE, showWarnings = FALSE)
    stop_if(file.copy(source_file, destination_file, overwrite = FALSE),
      paste0("Unable to copy frozen package file: ", source_file))
  }
}
copy_tree(original_final, output_dir)

copy_exact <- function(source, destination) {
  stop_if(file.exists(source), paste0("Missing renderer output: ", source))
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
  stop_if(file.copy(source, destination, overwrite = TRUE),
    paste0("Unable to install corrected file: ", destination))
}

# -------------------------------------------------------------------------
# Filter the three plot input tables. Keep the unfiltered source rows in the
# provenance package through the source paths and the group/reuse sidecars.
# -------------------------------------------------------------------------
input_files <- c(
  figure1_fraction = "figure1_sample_boxplot_values.csv",
  figure1_ratio = "figure1_mki67_ratio_sample_values.csv",
  pathway = "figure1_pathway_summary_sample_boxplot_values.csv"
)
filter_manifest <- rbindlist(lapply(names(input_files), function(input_name) {
  file_name <- unname(input_files[[input_name]])
  source_file <- file.path(expanded_input, file_name)
  stop_if(file.exists(source_file), paste0("Missing expanded input table: ", source_file))
  source_values <- fread(source_file, check.names = FALSE)
  stop_if("ObservationType" %in% names(source_values),
    paste0("Input table lacks ObservationType: ", source_file))
  sample_values <- source_values[ObservationType == "sample"]
  destination_file <- file.path(filtered_dir, file_name)
  fwrite(sample_values, destination_file, na = "")
  data.table(
    Input = input_name,
    File = file_name,
    SourcePath = normalizePath(source_file, mustWork = TRUE),
    SourceRows = nrow(source_values),
    SampleRows = nrow(sample_values),
    RemovedNonSampleRows = nrow(source_values) - nrow(sample_values),
    Filter = 'ObservationType == "sample"'
  )
}), fill = TRUE)
fwrite(filter_manifest, file.path(output_dir, "Corrected_Data", "input_filter_manifest.csv"), na = "")

# Copy the source registries that explain the stable PXD/sample joins and the
# protein-ratio key used by the candidate input. These are provenance copies,
# not alternative plotting inputs.
provenance_dir <- file.path(output_dir, "Corrected_Data", "source_provenance")
dir.create(provenance_dir, recursive = TRUE, showWarnings = FALSE)
provenance_sources <- c(
  figure1_sample_boxplot_source_registry = file.path(expanded_input, "figure1_sample_boxplot_source_registry.csv"),
  sample_boxplot_source_registry = file.path(expanded_input, "sample_boxplot_source_registry.csv"),
  biological_sample_count_record = file.path(expanded_input, "biological_sample_count_record.csv"),
  expanded_candidate_README = file.path(dirname(expanded_input), "README.md"),
  expanded_escc_inclusion_audit = file.path(dirname(expanded_input), "escc_inclusion_audit.csv"),
  mki67_ratio_protein_key = file.path(source_root, "data", "candidate", "figure1_mki67_ratio_protein_key.csv")
)
provenance_manifest <- rbindlist(lapply(names(provenance_sources), function(provenance_name) {
  source_file <- unname(provenance_sources[[provenance_name]])
  if (!file.exists(source_file)) {
    return(data.table(Name = provenance_name, SourcePath = source_file, Copied = FALSE))
  }
  destination_file <- file.path(provenance_dir, basename(source_file))
  stop_if(file.copy(source_file, destination_file, overwrite = TRUE),
    paste0("Unable to copy provenance file: ", source_file))
  data.table(Name = provenance_name, SourcePath = normalizePath(source_file, mustWork = TRUE), Copied = TRUE)
}), fill = TRUE)
fwrite(provenance_manifest, file.path(provenance_dir, "provenance_manifest.csv"), na = "")

fraction_full <- fread(file.path(expanded_input, input_files[["figure1_fraction"]]), check.names = FALSE)
fraction_sample <- fraction_full[ObservationType == "sample"]
ratio_full <- fread(file.path(expanded_input, input_files[["figure1_ratio"]]), check.names = FALSE)
ratio_sample <- ratio_full[ObservationType == "sample"]
pathway_full <- fread(file.path(expanded_input, input_files[["pathway"]]), check.names = FALSE)
pathway_sample <- pathway_full[ObservationType == "sample"]

# All 31 registry entries are retained even when a specific panel has no
# source-level sample observation. This lets S1 show an aligned empty/ND slot
# rather than silently changing the dataset universe.
group_registry <- fraction_full[, .(
  Category = first(Category),
  RowOrder = min(RowOrder),
  FullInputRows = .N,
  SampleInputRows = sum(ObservationType == "sample")
), by = .(PXD, SampleGroup)]
fwrite(group_registry, file.path(output_dir, "Corrected_Data", "31_group_registry.csv"), na = "")

coverage <- fraction_full[, .(
  FullRows = .N,
  SampleRows = sum(ObservationType == "sample"),
  NonSampleRows = sum(ObservationType != "sample"),
  SampleIDs = uniqueN(SampleID[ObservationType == "sample"])
), by = .(PXD, SampleGroup, Category, Dataset)]
setorder(coverage, Category, PXD, SampleGroup, Dataset)
fwrite(coverage, file.path(output_dir, "Corrected_Data", "figure1_group_modality_coverage.csv"), na = "")

# Reference reuse is deliberately not silently deduplicated: these rows are
# source-level observations used by the existing per-dataset comparison. The
# sidecar makes the non-independence visible for interpretation/statistics.
wp_sample <- fraction_sample[Dataset == "Whole proteome"]
wp_sample[, SourceKey := paste(ReferencePXD, SampleID, SourceFile, sep = "|")]
reuse <- wp_sample[duplicated(SourceKey) | duplicated(SourceKey, fromLast = TRUE)]
if (nrow(reuse)) {
  reuse[, ReusedAcrossSampleGroups := uniqueN(SampleGroup) > 1L, by = SourceKey]
  setorder(reuse, SourceKey, PXD, SampleGroup, SampleID)
} else {
  reuse[, ReusedAcrossSampleGroups := logical()]
}
fwrite(reuse, file.path(output_dir, "Corrected_Data", "whole_proteome_reference_reuse_audit.csv"), na = "")

# -------------------------------------------------------------------------
# Run the exact reviewed renderers against the filtered tables.
# -------------------------------------------------------------------------
renderer_dir <- file.path(project_root, "workflow", "corrected_renderers")
stop_if(dir.exists(renderer_dir), paste0("Missing corrected renderer directory: ", renderer_dir))
generated_dir <- file.path(output_dir, "generated_renderers")
dir.create(generated_dir, recursive = TRUE, showWarnings = FALSE)

run_renderer <- function(label, script, environment) {
  stop_if(file.exists(script), paste0("Missing renderer script: ", script))
  rendered <- system2(
    "Rscript", c("--vanilla", script, project_root), env = environment,
    stdout = TRUE, stderr = TRUE
  )
  writeLines(as.character(rendered), file.path(log_dir, paste0(label, ".log")))
  exit_status <- attr(rendered, "status")
  if (is.null(exit_status)) exit_status <- 0L
  stop_if(identical(as.integer(exit_status), 0L),
    paste0("Renderer failed (", label, "). See ", file.path(log_dir, paste0(label, ".log"))))
}

fraction_render_dir <- file.path(generated_dir, "figure1_fraction")
ratio_render_dir <- file.path(generated_dir, "figure1_ratio")
pathway_render_dir <- file.path(generated_dir, "figure2_pathways")
s1_render_dir <- file.path(generated_dir, "supplementary_s1")
dir.create(fraction_render_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(ratio_render_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(pathway_render_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(s1_render_dir, recursive = TRUE, showWarnings = FALSE)

run_renderer(
  "figure1_fraction",
  file.path(renderer_dir, "build_figure1_category_boxplot.R"),
  c(
    paste0("KLA_CANDIDATE_INPUT=", filtered_dir),
    paste0("KLA_CANDIDATE_OUTPUT=", fraction_render_dir)
  )
)
run_renderer(
  "figure1_ratio",
  file.path(renderer_dir, "build_figure1_mki67_ratio_boxplots.R"),
  c(
    paste0("KLA_MKI67_INPUT=", filtered_dir),
    paste0("KLA_MKI67_OUTPUT=", ratio_render_dir)
  )
)
run_renderer(
  "figure2_pathways",
  file.path(renderer_dir, "build_ddr_pathway_summary_boxplots.R"),
  c(
    paste0("KLA_CANDIDATE_INPUT=", filtered_dir),
    paste0("KLA_CANDIDATE_OUTPUT=", pathway_render_dir)
  )
)
run_renderer(
  "supplementary_s1",
  file.path(renderer_dir, "build_figure_s1_dataset_plots.R"),
  c(
    paste0("KLA_CANDIDATE_INPUT=", filtered_dir),
    paste0("KLA_S1_OUTPUT=", s1_render_dir),
    "KLA_SAMPLE_ONLY=TRUE",
    paste0("KLA_FULL_FIGURE1_INPUT=", file.path(expanded_input, input_files[["figure1_fraction"]])),
    paste0("KLA_S1A_EXPECTED_ROWS=", nrow(fraction_sample)),
    paste0("KLA_S1B_EXPECTED_ROWS=", ratio_sample[Denominator == "H3C1", .N])
  )
)

# -------------------------------------------------------------------------
# Install corrected figure files. Existing aliases are kept, but all aliases
# now point to the same corrected render so the package cannot mix versions.
# -------------------------------------------------------------------------
figure1_category_source <- file.path(fraction_render_dir, "Figure_1_DDR_fraction_candidate_category_boxplot_refined")
for (extension in c("png", "pdf")) {
  copy_exact(
    paste0(figure1_category_source, ".", extension),
    file.path(output_dir, "Figure_1", paste0("Figure_1a_DDR_fraction_boxplot.", extension))
  )
}

for (denominator in c("ACTB", "TUBB", "H3C1")) {
  ratio_source <- file.path(ratio_render_dir, paste0("Figure_1_MKI67_over_", denominator, "_boxplot"))
  destinations <- if (denominator == "H3C1") {
    file.path(output_dir, "Figure_1", paste0("Figure_1b_MKI67_over_H3C1_boxplot."))
  } else {
    c(
      file.path(output_dir, "Figure_1", "companion_controls", paste0("Figure_1_MKI67_over_", denominator, "_boxplot.")),
      file.path(output_dir, "Figure_1", "companion_controls", paste0("Figure_1_companion_MKI67_over_", denominator, "_boxplot."))
    )
  }
  for (destination_stem in destinations) {
    for (extension in c("png", "pdf")) {
      copy_exact(paste0(ratio_source, ".", extension), paste0(destination_stem, extension))
    }
  }
}

pathway_letters <- c(BER = "c", NER = "d", MMR = "e", FA = "f", HR = "g", AEJ = "h", NHEJ = "i")
for (pathway in names(pathway_letters)) {
  source_stem <- file.path(pathway_render_dir, paste0("Figure_2_DDR_pathway_summary_", pathway, "_barplot"))
  for (extension in c("png", "pdf")) {
    copy_exact(
      paste0(source_stem, ".", extension),
      file.path(output_dir, "Figure_2", paste0("Figure_2_DDR_pathway_summary_", pathway, "_barplot.", extension))
    )
    copy_exact(
      paste0(source_stem, ".", extension),
      file.path(output_dir, "Figure_2", paste0("Figure_2", pathway_letters[[pathway]], "_DDR_pathway_summary_", pathway, "_barplot.", extension))
    )
  }
}

for (stem in c(
  "Figure_S1a_DDR_fraction_by_PXD",
  "Figure_S1b_MKI67_over_H3C1_by_PXD",
  "Figure_S1b_MKI67_over_H3C1_11_detected_only"
)) {
  for (extension in c("png", "pdf")) {
    copy_exact(
      file.path(s1_render_dir, paste0(stem, ".", extension)),
      file.path(output_dir, "Supplementary_Figure_S1", paste0(stem, ".", extension))
    )
  }
}

# The frozen S3/S4 PDFs contain the complete protein-count titles, whereas the
# old PNG raster exports dropped the digits. Re-rasterize those four unchanged
# data panels from their PDFs so the two delivery formats have identical text.
pdftoppm_path <- Sys.which("pdftoppm")
stop_if(nzchar(pdftoppm_path), "pdftoppm is required to refresh the S3/S4 PNG title rasterization.")
matrix_rel_stems <- c(
  "Supplementary_Figure_S3/Figure_S3a_DDR_pathway_matrix_tumor_tissues",
  "Supplementary_Figure_S3/Figure_S3b_DDR_pathway_matrix_non_tumor_tissues",
  "Supplementary_Figure_S4/Figure_S4a_DDR_pathway_matrix_cancer_cell_lines",
  "Supplementary_Figure_S4/Figure_S4b_DDR_pathway_matrix_normal_cell_lines"
)
matrix_raster_dir <- file.path(generated_dir, "matrix_pdf_raster")
dir.create(matrix_raster_dir, recursive = TRUE, showWarnings = FALSE)
for (relative_stem in matrix_rel_stems) {
  pdf_file <- file.path(output_dir, paste0(relative_stem, ".pdf"))
  raster_stem <- file.path(matrix_raster_dir, basename(relative_stem))
  rasterized <- system2(
    pdftoppm_path,
    c("-png", "-r", "300", "-scale-to-x", "3300", "-scale-to-y", "1440", "-singlefile", pdf_file, raster_stem),
    stdout = TRUE, stderr = TRUE
  )
  raster_status <- attr(rasterized, "status")
  if (is.null(raster_status)) raster_status <- 0L
  stop_if(identical(as.integer(raster_status), 0L),
    paste0("Failed to rasterize matrix PDF: ", pdf_file))
  copy_exact(
    paste0(raster_stem, ".png"),
    file.path(output_dir, paste0(relative_stem, ".png"))
  )
}

# -------------------------------------------------------------------------
# Replace only statistical sidecars corresponding to corrected plots.
# -------------------------------------------------------------------------
copy_exact(
  file.path(fraction_render_dir, "figure1_category_omnibus_anova.csv"),
  file.path(output_dir, "Tables", "statistical_tests", "figure1_category_omnibus_anova.csv")
)
copy_exact(
  file.path(fraction_render_dir, "figure1_category_one_way_anova.csv"),
  file.path(output_dir, "Tables", "statistical_tests", "figure1_category_one_way_anova.csv")
)
copy_exact(
  file.path(fraction_render_dir, "figure1_category_boxplot_mean_median.csv"),
  file.path(output_dir, "Tables", "statistical_tests", "figure1_category_boxplot_mean_median.csv")
)
copy_exact(
  file.path(pathway_render_dir, "pathway_summary_two_way_anova.csv"),
  file.path(output_dir, "Tables", "statistical_tests", "pathway_summary_two_way_anova.csv")
)
copy_exact(
  file.path(filtered_dir, "figure1_mki67_ratio_significance.csv"),
  file.path(output_dir, "Tables", "statistical_tests", "figure1_mki67_ratio_significance.csv")
)

# Keep a renderer manifest in the corrected package. The complete package file
# list is added after this manifest is written so it also covers all frozen
# figures/tables copied from the original package.
renderer_manifest <- rbindlist(list(
  data.table(Figure = "Figure_1a", Renderer = "build_figure1_category_boxplot.R", Input = input_files[["figure1_fraction"]], ObservationType = "sample", InputRows = nrow(fraction_sample)),
  data.table(Figure = c("Figure_1b_ACTB", "Figure_1b_TUBB", "Figure_1b_H3C1"), Renderer = "build_figure1_mki67_ratio_boxplots.R", Input = input_files[["figure1_ratio"]], ObservationType = "sample", InputRows = c(sum(ratio_sample$Denominator == "ACTB"), sum(ratio_sample$Denominator == "TUBB"), sum(ratio_sample$Denominator == "H3C1"))),
  data.table(Figure = paste0("Figure_2", pathway_letters, "_", names(pathway_letters)), Renderer = "build_ddr_pathway_summary_boxplots.R", Input = input_files[["pathway"]], ObservationType = "sample", InputRows = pathway_sample[, .N, by = Pathway]$N),
  data.table(Figure = c("Figure_S1a", "Figure_S1b", "Figure_S1b_detected_only"), Renderer = "build_figure_s1_dataset_plots.R", Input = c(input_files[["figure1_fraction"]], input_files[["figure1_ratio"]], input_files[["figure1_ratio"]]), ObservationType = "sample", InputRows = c(nrow(fraction_sample), ratio_sample[Denominator == "H3C1", .N], ratio_sample[Denominator == "H3C1", .N]))
), fill = TRUE)
renderer_manifest[, Filter := 'ObservationType == "sample"']
fwrite(renderer_manifest, file.path(output_dir, "Corrected_Data", "corrected_renderer_manifest.csv"), na = "")

validation <- rbindlist(list(
  data.table(Check = "31-group registry", Observed = nrow(group_registry), Expected = 31L),
  data.table(Check = "Figure 1 DDR fraction sample rows", Observed = nrow(fraction_sample), Expected = 272L),
  data.table(Check = "Figure 1 MKI67 ratio sample rows", Observed = nrow(ratio_sample), Expected = 264L),
  data.table(Check = "Figure 2 pathway sample rows", Observed = nrow(pathway_sample), Expected = 504L),
  data.table(Check = "Figure 2 rows per pathway", Observed = min(pathway_sample[, .N, by = Pathway]$N), Expected = 72L),
  data.table(Check = "All corrected input rows are sample", Observed = uniqueN(c(fraction_sample$ObservationType, ratio_sample$ObservationType, pathway_sample$ObservationType)), Expected = 1L),
  data.table(Check = "Reused whole-proteome source keys", Observed = uniqueN(reuse$SourceKey), Expected = 11L)
), fill = TRUE)
validation[, Pass := Observed == Expected]
fwrite(validation, file.path(output_dir, "Corrected_Data", "correction_validation.csv"), na = "")
stop_if(all(validation$Pass), "Corrected package validation failed.")

message("CORRECTED_PACKAGE_COMPLETE ", normalizePath(output_dir, mustWork = TRUE))
message("Figure 1 fraction sample rows: ", nrow(fraction_sample))
message("Figure 1 ratio sample rows: ", nrow(ratio_sample))
message("Figure 2 pathway sample rows per pathway: ", paste(pathway_sample[, .N, by = Pathway]$N, collapse = ","))
message("31-group registry rows: ", nrow(group_registry))
