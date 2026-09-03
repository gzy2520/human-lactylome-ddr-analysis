suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(dplyr)
})

# ------------------------------------------------------------------------------
# User Rules & Style Specification
# 1. Higher priority of R programming language.
# 2. Strict ID standard: Data matched and indexed by stable UniProt BaseAccession;
#    Gene Symbol used strictly for display labels and titles.
# 3. Random seed: 25.
# 4. Canvas: Clean publication standard, pastel category palette, charcoal accents.
# ------------------------------------------------------------------------------
set.seed(25)

repo_root <- tryCatch(
  rprojroot::find_root(rprojroot::is_git_root),
  error = function(e) getwd()
)
if (!dir.exists(file.path(repo_root, "data"))) {
  repo_root <- getwd()
}

publication_font <- "Arial Unicode MS"
charcoal <- "#2F3437"
muted_text <- "#65717D"
grid_colour <- "#E5E8EC"
panel_border_colour <- "#C8CED6"

category_order <- c("normal_tissue", "cancer_tissue", "normal_cells", "cancer_cells")
category_labels <- c(
  normal_tissue = "non-tumor\ntissues",
  cancer_tissue = "tumor\ntissues",
  normal_cells = "normal\ncell lines",
  cancer_cells = "cancer\ncell lines"
)
category_fills <- c(
  normal_tissue = "#DCE9E2",
  cancer_tissue = "#F0DEDE",
  normal_cells = "#E7E1EE",
  cancer_cells = "#EEE4D2"
)

# 10 target regulators specified by user:
# AARS1, ACAT2, KRT18, SIRT2, PARK7, HDAC1, HDAC2, BRD4, SMARCA4, TRIM33
target_regulators <- data.table(
  GeneSymbol = c("AARS1", "ACAT2", "KRT18", "SIRT2", "PARK7", "HDAC1", "HDAC2", "BRD4", "SMARCA4", "TRIM33"),
  BaseAccession = c("P49588", "Q9BWD1", "P05783", "Q8IXJ6", "Q99497", "Q13547", "Q92769", "O60885", "P51532", "Q9UPN9"),
  Role = c("Writer", "Writer", "Writer", "Eraser", "Writer-Eraser", "Writer-Eraser", "Writer-Eraser", "Reader", "Reader", "Reader")
)

# Load sample percentiles input
input_path <- file.path(
  repo_root,
  "data/candidate/escc_inclusion_20260903_pxd065830_tumor_reference/target_10_regulators_sample_percentiles.csv"
)

if (!file.exists(input_path)) {
  stop("Input file not found: ", input_path)
}

sample_data <- fread(input_path)

# Verify matching standard: BaseAccession must be present
if (!all(target_regulators$BaseAccession %in% unique(sample_data$BaseAccession))) {
  missing_acc <- setdiff(target_regulators$BaseAccession, unique(sample_data$BaseAccession))
  stop("Missing required BaseAccessions in data: ", paste(missing_acc, collapse = ", "))
}

# Output directories
out_dir_candidate <- file.path(
  repo_root,
  "results/candidate/escc_inclusion_20260903_pxd065830_tumor_reference/regulator_sample_percentiles"
)
out_dir_baseline <- file.path(
  out_dir_candidate,
  "baseline_30_datasets"
)
dir.create(out_dir_candidate, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir_baseline, recursive = TRUE, showWarnings = FALSE)

render_scatter_plot <- function(accession_id, data_subset, output_file_stem, scope_subtitle) {
  # Index strictly by BaseAccession
  reg_info <- target_regulators[BaseAccession == accession_id]
  gene_display <- reg_info$GeneSymbol
  role_display <- reg_info$Role
  
  df <- copy(data_subset[BaseAccession == accession_id])
  df[, Category := factor(Category, levels = category_order, labels = category_labels[category_order])]
  
  # Category sample count annotations
  sample_counts <- df[, .(
    N = .N,
    DetN = sum(Detected),
    Y = 103.5
  ), by = Category]
  
  p <- ggplot(df, aes(x = Category, y = WholeProteomePercentile)) +
    # Reference grid lines
    geom_hline(yintercept = seq(0, 100, 25), colour = grid_colour, linewidth = 0.5) +
    # Jittered individual sample points (strictly seed = 25)
    geom_jitter(
      aes(fill = Category),
      shape = 21,
      size = 2.8,
      stroke = 0.55,
      colour = charcoal,
      alpha = 0.75,
      position = position_jitter(width = 0.22, height = 0, seed = 25)
    ) +
    # Median horizontal indicator bar
    stat_summary(
      fun = median,
      geom = "crossbar",
      colour = charcoal,
      linewidth = 0.7,
      width = 0.45
    ) +
    # Top sample size label
    geom_text(
      data = sample_counts,
      aes(x = Category, y = Y, label = paste0("n = ", N)),
      inherit.aes = FALSE,
      size = 4.2,
      colour = muted_text,
      family = publication_font
    ) +
    scale_fill_manual(values = unname(category_fills[category_order]), guide = "none") +
    scale_y_continuous(
      limits = c(-2, 108),
      breaks = seq(0, 100, 25),
      labels = function(y) ifelse(y <= 100, paste0(y, "%"), ""),
      expand = expansion(mult = c(0, 0))
    ) +
    labs(
      title = paste0(gene_display, " (", accession_id, ")"),
      subtitle = paste0("UniProt: ", accession_id, " · Lactylation ", role_display, " · ", scope_subtitle),
      x = NULL,
      y = "Whole-proteome expression percentile (%)",
      caption = paste(
        "Each point represents an individual whole-proteome reference sample; horizontal black bar indicates category median.",
        "Percentile = 100 × (rank - 1) / (N_features - 1) within each sample's ordinary proteome.",
        "Undetected samples sit at 0%. Random jitter seed = 25.",
        sep = "\n"
      )
    ) +
    theme_minimal(base_size = 13, base_family = publication_font) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(colour = grid_colour, linewidth = 0.5),
      panel.border = element_rect(colour = panel_border_colour, fill = NA, linewidth = 0.6),
      axis.text.x = element_text(size = 13, colour = charcoal, face = "bold", margin = margin(t = 6)),
      axis.text.y = element_text(size = 12.5, colour = charcoal),
      axis.title.y = element_text(size = 13.5, face = "bold", colour = charcoal, margin = margin(r = 10)),
      plot.title = element_text(size = 18, face = "bold", colour = charcoal, hjust = 0.5, margin = margin(b = 2)),
      plot.subtitle = element_text(size = 11.5, colour = muted_text, hjust = 0.5, margin = margin(b = 10)),
      plot.caption = element_text(size = 9.5, colour = muted_text, hjust = 0.5, margin = margin(t = 10)),
      plot.margin = margin(12, 16, 12, 12),
      plot.background = element_rect(fill = "white", colour = NA)
    )
  
  png_file <- paste0(output_file_stem, ".png")
  pdf_file <- paste0(output_file_stem, ".pdf")
  ggsave(png_file, p, width = 8.5, height = 7, dpi = 300, bg = "white")
  ggsave(pdf_file, p, width = 8.5, height = 7, bg = "white", device = cairo_pdf)
  message("Generated: ", basename(png_file), " and ", basename(pdf_file))
}

# 1. Render primary candidate scope (31 datasets, including ESCC PXD065830 tumor reference)
message("Rendering candidate scope scatter plots (31 datasets)...")
for (i in seq_len(nrow(target_regulators))) {
  acc <- target_regulators$BaseAccession[[i]]
  gene <- target_regulators$GeneSymbol[[i]]
  file_stem <- file.path(out_dir_candidate, paste0("Figure_3_regulator_sample_percentile_", gene, "_", acc))
  render_scatter_plot(
    accession_id = acc,
    data_subset = sample_data,
    output_file_stem = file_stem,
    scope_subtitle = "Whole-proteome expression percentile across reference samples"
  )
}

# 2. Render baseline scope (30 frozen reference datasets)
message("Rendering baseline scope scatter plots (30 datasets)...")
baseline_subset <- sample_data[Scope == "baseline_30_datasets"]
for (i in seq_len(nrow(target_regulators))) {
  acc <- target_regulators$BaseAccession[[i]]
  gene <- target_regulators$GeneSymbol[[i]]
  file_stem <- file.path(out_dir_baseline, paste0("Figure_3_regulator_sample_percentile_", gene, "_", acc, "_30datasets"))
  render_scatter_plot(
    accession_id = acc,
    data_subset = baseline_subset,
    output_file_stem = file_stem,
    scope_subtitle = "30 frozen reference datasets"
  )
}

message("All 10 regulator sample-level scatter plots successfully generated.")
