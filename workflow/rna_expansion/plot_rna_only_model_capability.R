#!/usr/bin/env Rscript
# Honest, source-held-out visualization of the frozen RNA-only Kla fraction task.
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})
set.seed(25)
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
out <- args[[1L]]
if (dir.exists(out) && length(list.files(out, all.files = TRUE, no.. = TRUE)))
  stop("Refusing to overwrite a nonempty output directory")

pred_file <- "outputs/20260924_rna_ref_combined_equal/material_oof_predictions.csv"
perf_file <- "outputs/20260924_rna_ref_combined_equal/performance_summary.csv"
labels_file <- "outputs/20260923_rna_only_kla_fraction/material_labels_28.csv"
names_file <- "config/lactate_metabolism/material_labels.csv"
input_files <- c(pred_file, perf_file, labels_file, names_file)
stopifnot(all(file.exists(input_files)))
pred <- fread(pred_file)
perf <- fread(perf_file)
labels <- fread(labels_file)
display <- fread(names_file)
stopifnot(nrow(pred) == 28L, nrow(labels) == 28L,
          uniqueN(pred$ConnGroup) == 15L,
          !anyDuplicated(pred$ReferenceKey),
          !anyDuplicated(labels$ReferenceKey),
          !anyDuplicated(display$ReferenceKey),
          setequal(pred$ReferenceKey, labels$ReferenceKey))
d <- merge(pred[, .(ReferenceKey, ConnGroup, ObservedPercent,
                     Baseline = BalancedMedianBaseline, RNA = All16_ridge)],
           labels[, .(ReferenceKey, GroupIDs, Category,
                      FrozenObserved = ObservedPercent,
                      FrozenConnGroup = ConnGroup)], by = "ReferenceKey")
d <- merge(d, display, by = "ReferenceKey")
stopifnot(nrow(d) == 28L, !anyNA(d),
          identical(d$ObservedPercent, d$FrozenObserved),
          identical(d$ConnGroup, d$FrozenConnGroup))
d[, c("FrozenObserved", "FrozenConnGroup") := NULL]
setorder(d, ReferenceKey)

error <- function(which) {
  x <- d[[which]] - d$ObservedPercent
  by_source <- d[, .(MAE = mean(abs(get(which) - ObservedPercent))),
                 by = ConnGroup]
  data.table(Model = which, MaterialMAE_pp = mean(abs(x)),
             SourceMacroMAE_pp = mean(by_source$MAE),
             RMSE_pp = sqrt(mean(x^2)),
             Within5pp = sum(abs(x) <= 5),
             PredictionMin = min(d[[which]]),
             PredictionMax = max(d[[which]]))
}
metrics <- rbind(error("Baseline"), error("RNA"))
for (i in seq_len(nrow(metrics))) {
  key <- if (metrics$Model[[i]] == "Baseline")
    "BalancedMedianBaseline" else "All16_ridge"
  recorded <- perf[Model == key]
  stopifnot(nrow(recorded) == 1L,
            abs(metrics$MaterialMAE_pp[[i]] - recorded$MaterialMAE_pp) < 1e-10,
            abs(metrics$SourceMacroMAE_pp[[i]] -
                  recorded$ConnectedGroupMeanMAE_pp) < 1e-10,
            abs(metrics$RMSE_pp[[i]] - recorded$MaterialRMSE_pp) < 1e-10)
}

baseline_name <- sprintf("No-RNA baseline\nMAE %.2f pp",
                         metrics[Model == "Baseline", MaterialMAE_pp])
rna_name <- sprintf("RNA-only, 16 modules\nMAE %.2f pp",
                    metrics[Model == "RNA", MaterialMAE_pp])
long <- melt(d[, .(ReferenceKey, RNALabel, Category, ObservedPercent,
                   Baseline, RNA)],
             id.vars = c("ReferenceKey", "RNALabel", "Category",
                         "ObservedPercent"),
             measure.vars = c("Baseline", "RNA"),
             variable.name = "Model", value.name = "PredictedPercent")
long[, Panel := factor(fifelse(Model == "Baseline", baseline_name, rna_name),
                       levels = c(baseline_name, rna_name))]
long[, Category := factor(Category,
  levels = c("normal_tissue", "cancer_tissue", "normal_cells", "cancer_cells"),
  labels = c("Normal tissue", "Cancer tissue", "Normal cells", "Cancer cells"))]
palette <- c("Normal tissue" = "#4778A8", "Cancer tissue" = "#E78637",
             "Normal cells" = "#54A06B", "Cancer cells" = "#A26BA6")
band <- data.table(ObservedPercent = seq(0, 40, length.out = 200))
band[, `:=`(Low = pmax(0, ObservedPercent - 5),
            High = pmin(40, ObservedPercent + 5))]
theme_pub <- theme_minimal(base_family = "sans", base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "#E8ECEF", linewidth = .3),
        plot.title = element_text(face = "bold", size = 15, colour = "#17344A"),
        plot.subtitle = element_text(size = 10, colour = "#465760"),
        strip.text = element_text(face = "bold", size = 10, colour = "#17344A"),
        strip.background = element_rect(fill = "#F1F5F7", colour = NA),
        axis.title = element_text(colour = "#253746"),
        axis.text = element_text(colour = "#465760"),
        legend.title = element_blank(),
        plot.margin = margin(9, 12, 9, 9))

scatter <- ggplot(long, aes(x = ObservedPercent, y = PredictedPercent)) +
  geom_ribbon(data = band, aes(x = ObservedPercent, ymin = Low, ymax = High),
              inherit.aes = FALSE, fill = "#C7D1D8", alpha = .28) +
  geom_abline(slope = 1, intercept = 0, linewidth = .65,
              linetype = "22", colour = "#56646C") +
  geom_point(aes(fill = Category), shape = 21, colour = "#FFFFFF",
             stroke = .42, size = 2.9, alpha = .95) +
  facet_wrap(~Panel, nrow = 1) +
  scale_fill_manual(values = palette, drop = FALSE) +
  scale_x_continuous(limits = c(0, 40), breaks = seq(0, 40, 10),
                     expand = expansion(mult = .01)) +
  scale_y_continuous(limits = c(0, 40), breaks = seq(0, 40, 10),
                     expand = expansion(mult = .01)) +
  coord_fixed() +
  labs(x = "Observed detected Kla-protein fraction (%)",
       y = "Held-out prediction (%)",
       subtitle = "Dashed line: perfect prediction; shading: ±5 percentage points (visual guide)") +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE,
                             override.aes = list(size = 3.2))) +
  theme_pub +
  theme(legend.position = "bottom",
        strip.text = element_text(lineheight = 1.12),
        legend.text = element_text(size = 9),
        panel.spacing.x = unit(1.3, "lines"))

source_names <- c(
  "Torn tendon", "Lung / placenta", "Scar / adjacent skin",
  "Hippocampus", "Sperm", "BPH / prostate cancer",
  "Adjacent liver / HCC", "ESCC", "Shared cell lines (9)",
  "HepG2", "PC-3M", "NSC / MES28", "HMC3", "HK-2", "HUVEC"
)
stopifnot(length(source_names) == 15L, setequal(d$ConnGroup, seq_len(15L)))
sources <- d[, .(BaselineMAE = mean(abs(Baseline - ObservedPercent)),
                 RNAMAE = mean(abs(RNA - ObservedPercent)), N = .N),
             by = ConnGroup]
sources[, Source := sprintf("%02d  %s", ConnGroup,
                            source_names[ConnGroup])]
setorder(sources, ConnGroup)
sources[, Source := factor(Source, levels = rev(Source))]
source_long <- melt(sources[, .(Source, BaselineMAE, RNAMAE)],
                    id.vars = "Source", variable.name = "Model",
                    value.name = "MAE")
source_long[, Model := factor(Model, levels = c("BaselineMAE", "RNAMAE"),
                             labels = c("No-RNA baseline", "RNA-only"))]
source_plot <- ggplot(sources, aes(y = Source)) +
  geom_segment(aes(x = BaselineMAE, xend = RNAMAE, yend = Source),
               colour = "#B8C3CA", linewidth = .75) +
  geom_point(data = source_long, aes(x = MAE, colour = Model),
             inherit.aes = FALSE, y = source_long$Source, size = 2.5) +
  scale_colour_manual(values = c("No-RNA baseline" = "#788792",
                                 "RNA-only" = "#087E8B")) +
  scale_x_continuous(limits = c(0, 21), breaks = seq(0, 20, 5),
                     expand = expansion(mult = c(.01, .02))) +
  labs(title = "Error by held-out source",
       subtitle = sprintf("Each source counts once: mean MAE %.2f pp without RNA vs %.2f pp with RNA",
                          metrics[Model == "Baseline", SourceMacroMAE_pp],
                          metrics[Model == "RNA", SourceMacroMAE_pp]),
       x = "Mean absolute error within source (percentage points)", y = NULL) +
  theme_pub +
  theme(legend.position = "top", legend.justification = "left",
        panel.grid.major.y = element_blank(), axis.text.y = element_text(size = 9),
        plot.title = element_text(size = 13))

combined <- scatter / source_plot + plot_layout(heights = c(1.15, 1)) +
  plot_annotation(
    title = "RNA-only predictions remain compressed toward the middle",
    subtitle = "28 independent materials; 15 source-connected held-out folds; expanded scar + tendon RNA references",
    caption = paste(
      "The outcome is a detected Kla-bearing protein fraction, not site occupancy.",
      "RNA reference expansion adds no new protein-labeled material."),
    theme = theme(plot.title = element_text(face = "bold", size = 18,
                                            colour = "#17344A"),
                  plot.subtitle = element_text(size = 10.5,
                                               colour = "#465760"),
                  plot.caption = element_text(size = 8.5, hjust = 0,
                                              colour = "#56646C"),
                  plot.margin = margin(12, 16, 12, 12)))

# A material-level companion lets readers inspect every prediction explicitly.
material_long <- melt(d[, .(ReferenceKey, RNALabel, ObservedPercent,
                            Baseline, RNA)],
                      id.vars = c("ReferenceKey", "RNALabel"),
                      measure.vars = c("ObservedPercent", "Baseline", "RNA"),
                      variable.name = "Series", value.name = "Percent")
material_long[, Series := factor(Series,
  levels = c("ObservedPercent", "Baseline", "RNA"),
  labels = c("Observed", "No-RNA baseline", "RNA-only"))]
ranked <- d[order(ObservedPercent), RNALabel]
material_long[, RNALabel := factor(RNALabel, levels = ranked)]
material_plot <- ggplot(material_long,
                        aes(x = Percent, y = RNALabel, colour = Series)) +
  geom_point(aes(shape = Series), size = 2.6, stroke = .7,
             position = position_dodge(width = .53)) +
  scale_colour_manual(values = c("Observed" = "#202B35",
                                 "No-RNA baseline" = "#89969F",
                                 "RNA-only" = "#087E8B")) +
  scale_shape_manual(values = c("Observed" = 16,
                                "No-RNA baseline" = 1,
                                "RNA-only" = 17)) +
  scale_x_continuous(limits = c(0, 40), breaks = seq(0, 40, 10),
                     expand = expansion(mult = c(.01, .02))) +
  labs(title = "Held-out predictions for every material",
       subtitle = "Materials sorted by observed fraction; the RNA model underestimates high values and overestimates low values",
       x = "Detected Kla-protein fraction (%)", y = NULL,
       caption = "One point per independent RNA material; the underlying 31 protein groups are not treated as independent here.") +
  theme_pub +
  theme(legend.position = "top", legend.justification = "left",
        panel.grid.major.y = element_blank(), axis.text.y = element_text(size = 9),
        plot.caption = element_text(size = 8, hjust = 0))

dir.create(out, recursive = TRUE, showWarnings = FALSE)
save_both <- function(plot, stem, width, height) {
  ggsave(file.path(out, paste0(stem, ".png")), plot,
         device = ragg::agg_png, width = width, height = height,
         units = "in", dpi = 300, bg = "white")
  ggsave(file.path(out, paste0(stem, ".pdf")), plot,
         device = grDevices::cairo_pdf, width = width, height = height,
         units = "in", bg = "white")
}
save_both(combined, "Figure_RNAonly_model_capability", 11.8, 10.2)
save_both(material_plot, "Figure_RNAonly_material_predictions", 9.2, 10.7)
fwrite(d, file.path(out, "figure_material_data.csv"))
fwrite(sources, file.path(out, "figure_source_error_data.csv"))
fwrite(metrics, file.path(out, "figure_metrics.csv"))
fwrite(data.table(File = input_files,
                  MD5 = unname(tools::md5sum(input_files))),
       file.path(out, "input_md5.csv"))
writeLines(trimws(capture.output(sessionInfo()), which = "right"),
           file.path(out, "sessionInfo.txt"))
print(metrics)
