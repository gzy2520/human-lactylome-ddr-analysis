#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

project_root <- normalizePath(".", mustWork = TRUE)
analysis_name <- "kla_ddr_four_category_three_embedding_pathway_grids_33groups_v1"
table_dir <- file.path(project_root, "reanalysis/results/tables", analysis_name)
figure_dir <- file.path(project_root, "reanalysis/results/figures", analysis_name)
v4_table_dir <- file.path(
  project_root,
  "reanalysis/results/tables/kla_ddr_pathway_pie_umap_33groups_v4_bp_semantic"
)

required_tables <- c(
  "embedding_coordinates_507.csv",
  "embedding_parameters.csv",
  "category_membership_507_long.csv",
  "category_pathway_plot_data.csv",
  "category_pathway_summary.csv",
  "figure_manifest_12_grids.csv",
  "input_file_audit.csv",
  "session_info.txt"
)
required_paths <- file.path(table_dir, required_tables)
if (any(!file.exists(required_paths))) {
  stop(
    "Missing required output(s): ",
    paste(required_paths[!file.exists(required_paths)], collapse = "; ")
  )
}

assert <- function(condition, message) {
  if (!isTRUE(condition)) {
    stop(message, call. = FALSE)
  }
}

coordinates <- fread(file.path(table_dir, "embedding_coordinates_507.csv"))
parameters <- fread(file.path(table_dir, "embedding_parameters.csv"))
membership <- fread(file.path(table_dir, "category_membership_507_long.csv"))
plot_data <- fread(file.path(table_dir, "category_pathway_plot_data.csv"))
summary_table <- fread(file.path(table_dir, "category_pathway_summary.csv"))
manifest <- fread(file.path(table_dir, "figure_manifest_12_grids.csv"))
raw_umap <- fread(
  file.path(v4_table_dir, "umap_raw_coordinates_v4_bp_semantic.csv")
)

assert(
  nrow(coordinates) == 507L &&
    uniqueN(coordinates$BaseAccession) == 507L &&
    all(
      is.finite(
        as.matrix(
          coordinates[
            ,
            .(UMAP_1, UMAP_2, TSNE_1, TSNE_2, PCA_1, PCA_2)
          ]
        )
      )
    ),
  "Expected 507 finite coordinates for UMAP, t-SNE, and PCA."
)

umap_check <- merge(
  coordinates[, .(BaseAccession, UMAP_1, UMAP_2)],
  raw_umap,
  by = "BaseAccession",
  suffixes = c("_new", "_saved")
)
assert(
  nrow(umap_check) == 507L &&
    isTRUE(
      all.equal(umap_check$UMAP_1_new, umap_check$UMAP_1_saved, tolerance = 0)
    ) &&
    isTRUE(
      all.equal(umap_check$UMAP_2_new, umap_check$UMAP_2_saved, tolerance = 0)
    ),
  "UMAP coordinates are not the exact audited V4 raw coordinates."
)

category_counts <- membership[
  DetectedInCategory == TRUE,
  .N,
  by = .(Category, CategoryOrder)
][order(CategoryOrder)]
assert(
  identical(category_counts$Category, c(
    "normal_tissue",
    "normal_cells",
    "cancer_tissue",
    "cancer_cells"
  )) &&
    identical(category_counts$N, c(183L, 471L, 178L, 383L)),
  "Four-category counts or requested category order are incorrect."
)

assert(
  nrow(plot_data) == 507L * 4L * 9L &&
    uniqueN(plot_data$BaseAccession) == 507L &&
    uniqueN(plot_data$Category) == 4L &&
    uniqueN(plot_data$DisplayLabel) == 9L,
  "Expected a complete 507 x four-category x nine-pathway plotting table."
)
assert(
  plot_data[
    ,
    .N,
    by = .(Category, DisplayLabel)
  ][, all(N == 507L)] &&
    summary_table[
      ,
      all(
        PromotingCount +
          SuppressingCount +
          DetectedUnassignedCount ==
          CategoryProteinCount
      )
    ],
  "Each category-pathway panel must contain all 507 context proteins and valid category totals."
)
assert(
  plot_data[
    DetectedInCategory == FALSE,
    all(Status == "Outside category")
  ] &&
    plot_data[
      DetectedInCategory & Score == 0L,
      all(Status == "Detected, not assigned (0)")
    ] &&
    plot_data[
      DetectedInCategory & Score == 1L,
      all(Status == "Promoting (+1)")
    ] &&
    plot_data[
      DetectedInCategory & Score == -1L,
      all(Status == "Suppressing (-1)")
    ],
  "Category detection and signed pathway status encoding are inconsistent."
)

assert(
  parameters[Method == "Shared" & Parameter == "RandomSeedRecorded", Value] == "25" &&
    parameters[Method == "UMAP" & Parameter == "RandomSeed", Value] == "25" &&
    parameters[Method == "t-SNE" & Parameter == "RandomSeed", Value] == "25" &&
    parameters[
      Method == "PCA" & Parameter == "RandomSeedRecorded",
      Value
    ] == "25",
  "Seed 25 is not consistently recorded for all three embeddings."
)

assert(
  nrow(manifest) == 12L &&
    uniqueN(manifest$Category) == 4L &&
    uniqueN(manifest$Method) == 3L,
  "Expected exactly four categories x three methods = 12 figure entries."
)
figure_paths <- unlist(
  lapply(
    c("PNG", "PDF", "SVG"),
    function(column) file.path(project_root, manifest[[column]])
  ),
  use.names = FALSE
)
assert(
  length(figure_paths) == 36L &&
    all(file.exists(figure_paths)) &&
    all(file.info(figure_paths)$size > 0),
  "Expected 36 nonempty files for 12 figures in PNG/PDF/SVG formats."
)

cat(
  paste0(
    "PASS: 12 category-by-embedding pathway grids reuse one coordinate set per ",
    "method, preserve the four fixed category counts, and record seed 25.\n"
  )
)
