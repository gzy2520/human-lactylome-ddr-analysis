#!/usr/bin/env Rscript

# Render the upright Figure 1 candidate layout: four biological-category
# panels (columns), whole proteome and Kla as two vertical boxplots in each panel,
# and each dot as an individual deposited source-sample observation.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

set.seed(25)

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1L) normalizePath(args[[1L]], mustWork = TRUE) else normalizePath(".", mustWork = TRUE)
stop_if <- function(condition, message) if (!isTRUE(condition)) stop(message, call. = FALSE)

publication_font <- "Arial Unicode MS"
whole_proteome_colour <- "#4E79A7"
kla_colour <- "#F28E2B"
mean_colour <- "#C0392B"
charcoal <- "#2F3437"
muted_text <- "#65717D"
grid_colour <- "#D9DDE3"
panel_border_colour <- "#C8CED6"

category_order <- c("normal_tissue", "cancer_tissue", "normal_cells", "cancer_cells")
category_labels <- c(
  normal_tissue = "non-tumor tissues",
  cancer_tissue = "tumor tissues",
  normal_cells = "normal cell lines",
  cancer_cells = "cancer cell lines"
)
category_fills <- c(
  normal_tissue = "#DCE9E2",
  cancer_tissue = "#F0DEDE",
  normal_cells = "#E7E1EE",
  cancer_cells = "#EEE4D2"
)
dataset_order <- c("Whole proteome", "Lactylome (Kla)")

candidate_dir <- normalizePath(Sys.getenv(
  "KLA_CANDIDATE_INPUT", unset = file.path(project_root, "data", "candidate")
), mustWork = TRUE)
output_dir <- Sys.getenv("KLA_CANDIDATE_OUTPUT", unset = file.path(project_root, "results", "candidate"))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
values_path <- file.path(candidate_dir, "figure1_sample_boxplot_values.csv")
stop_if(file.exists(values_path), paste0("Missing Figure 1 sample input: ", values_path))
values <- fread(values_path, check.names = FALSE)
stop_if("ObservationType" %in% names(values), "Figure 1 input lacks ObservationType.")
values <- values[ObservationType == "sample"]
stop_if(nrow(values) > 0L, "Figure 1 has no sample-level observations after filtering.")

required_values <- c("PXD", "SampleGroup", "Category", "Dataset", "SampleID", "DdrFractionPercentage")
stop_if(all(required_values %in% names(values)), "Figure 1 sample input schema is incomplete.")
stop_if(all(values$Dataset %in% dataset_order), "Figure 1 contains an unknown modality.")
stop_if(all(values$Category %in% category_order), "Figure 1 contains an unknown category.")
stop_if(all(is.finite(values$DdrFractionPercentage) & values$DdrFractionPercentage >= 0 & values$DdrFractionPercentage <= 100),
  "Figure 1 fractions must be finite percentages from 0 to 100.")
stop_if(!anyDuplicated(values[, .(PXD, SampleGroup, Dataset, SampleID)]),
  "Figure 1 source-sample points must be unique within modality.")

source(file.path(project_root, "R", "candidate", "boxplot_significance.R"), local = TRUE)

# Retain category-specific ANOVA sidecar for contract compatibility
anova <- compute_figure1_sample_one_way_anova(values, category_order, dataset_order)
stop_if(nrow(anova) == length(category_order) && all(is.finite(anova$PValue)),
  "Figure 1 one-way ANOVA did not produce four finite category tests.")
fwrite(anova, file.path(output_dir, "figure1_category_one_way_anova.csv"), na = "")

# Compute four-category omnibus ANOVA
omnibus_anova <- compute_figure1_category_omnibus_anova(values, category_order, dataset_order)
fwrite(omnibus_anova, file.path(output_dir, "figure1_category_omnibus_anova.csv"), na = "")

cat_row <- omnibus_anova[Test == "two-way ANOVA" & Term == "CategoryFactor"]
kla_row <- omnibus_anova[Dataset == "Lactylome (Kla)"]
wp_row <- omnibus_anova[Dataset == "Whole proteome"]

format_anova_p <- function(p) {
  if (!is.finite(p)) return("NA")
  if (p < 2.2e-16) return("p < 2.2e-16")
  paste0("p = ", formatC(p, format = "e", digits = 2))
}
format_f_stat <- function(f) {
  if (!is.finite(f)) return("NA")
  sprintf("%.2f", f)
}

subtitle_text <- paste0(
  "Four-category ANOVA ", format_anova_p(cat_row$PValue[[1L]]),
  " (F = ", format_f_stat(cat_row$FStatistic[[1L]]), ", ", cat_row$Significance[[1L]], ")"
)

values[, CategoryLabel := factor(Category, levels = category_order, labels = unname(category_labels[category_order]))]
values[, Dataset := factor(Dataset, levels = dataset_order)]
panel_counts <- values[, .(N = .N), by = .(CategoryLabel, Dataset)]
summary_stats <- values[, .(
  Mean = mean(DdrFractionPercentage),
  Median = median(DdrFractionPercentage),
  N = .N
), by = .(Category, CategoryLabel, Dataset)]
fwrite(summary_stats, file.path(output_dir, "figure1_category_boxplot_mean_median.csv"), na = "")

max_val <- max(values$DdrFractionPercentage, na.rm = TRUE)
y_limit <- max(20, ceiling((max_val + 3) / 5) * 5)

dodge_width <- 0.72
summary_stats[, CatIdx := as.numeric(CategoryLabel)]
summary_stats[, XPos := fifelse(Dataset == "Whole proteome", CatIdx - dodge_width / 4, CatIdx + dodge_width / 4)]

figure_plot <- ggplot(values, aes(x = CategoryLabel, y = DdrFractionPercentage, fill = Dataset)) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5), colour = "#E5E7EB", linetype = "dashed", linewidth = 0.5) +
  geom_boxplot(
    position = position_dodge(width = dodge_width),
    width = 0.58, outlier.shape = NA, colour = charcoal,
    linewidth = 0.55, median.linewidth = 0.75, alpha = 0.82, na.rm = TRUE
  ) +
  geom_segment(
    data = summary_stats,
    aes(x = XPos - 0.16, xend = XPos + 0.16, y = Mean, yend = Mean),
    inherit.aes = FALSE, colour = mean_colour, linewidth = 0.80
  ) +
  geom_point(
    aes(fill = Dataset),
    position = position_jitterdodge(jitter.width = 0.14, dodge.width = dodge_width, seed = 25),
    shape = 21, size = 2.6, stroke = 0.50, colour = "white", alpha = 0.88, na.rm = TRUE
  ) +
  geom_text(
    data = summary_stats,
    aes(x = XPos, y = y_limit * 0.94, label = paste0("n=", N)),
    inherit.aes = FALSE, size = 3.8, family = publication_font, colour = muted_text, fontface = "bold"
  ) +
  scale_fill_manual(values = c("Whole proteome" = whole_proteome_colour, "Lactylome (Kla)" = kla_colour), breaks = dataset_order) +
  scale_y_continuous(limits = c(0, y_limit), breaks = scales::pretty_breaks(n = 5), labels = function(y) paste0(y, "%"), expand = expansion(mult = c(0, 0))) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE, keyheight = grid::unit(0.55, "cm"), keywidth = grid::unit(0.85, "cm"))) +
  labs(
    title = "DDR annotated protein fraction across four biological categories",
    subtitle = subtitle_text,
    x = NULL, y = "GO-DDR annotated protein fraction (%)", fill = NULL,
    caption = NULL
  ) +
  theme_minimal(base_size = 14, base_family = publication_font) +
  theme(
    panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(colour = grid_colour, linewidth = 0.50),
    panel.border = element_blank(),
    axis.line.x = element_line(colour = "#8C939E", linewidth = 0.60),
    axis.line.y = element_line(colour = "#8C939E", linewidth = 0.60),
    axis.text.x = element_text(size = 13.0, colour = charcoal, face = "bold", margin = margin(t = 6)),
    axis.text.y = element_text(size = 13.0, colour = charcoal),
    axis.title.y = element_text(size = 15, face = "bold", colour = charcoal, margin = margin(r = 10)),
    plot.title = element_text(size = 16.5, face = "bold", colour = charcoal, hjust = 0.5, margin = margin(b = 4)),
    plot.subtitle = element_text(size = 12.0, colour = muted_text, hjust = 0.5, margin = margin(b = 10)),
    legend.position = "top", legend.direction = "horizontal",
    legend.text = element_text(size = 13.0, colour = charcoal),
    legend.key.spacing.x = grid::unit(0.35, "cm"), legend.background = element_rect(fill = "white", colour = NA),
    legend.margin = margin(1, 0, 4, 0),
    plot.margin = margin(10, 16, 12, 12), plot.background = element_rect(fill = "white", colour = NA)
  )

stem <- file.path(output_dir, "Figure_1_DDR_fraction_candidate_category_boxplot_refined")
ggsave(paste0(stem, ".png"), figure_plot, width = 8.5, height = 7.0, dpi = 300, bg = "white", device = ragg::agg_png)
ggsave(paste0(stem, ".pdf"), figure_plot, width = 8.5, height = 7.0, bg = "white", device = cairo_pdf)
message("Wrote upright Figure 1 source-sample boxplot: ", stem, ".png/.pdf")
