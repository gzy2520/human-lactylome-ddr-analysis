# 11_figures_v4.R
# Publication-quality automated figure generation and source data export for B-R v4.
# Plots per-fold performance, material calibration, module incremental attribution, and S1/S2 sensitivity.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

pdf(file = tempfile(fileext = ".pdf")) # prevent an implicit Rplots.pdf in the repository
set.seed(25)
W <- "/Users/gzy2520/Desktop/Research/kla/.worktrees/kla31-rnaseq-20260912"
setwd(W)

out_dir <- Sys.getenv("BR_V4_OUT", "outputs/20260921_br_v4")

# --- Figure 1: Per-fold AP and LogLoss Comparison across Models ---
cat("Generating Figure 1: Per-fold Performance...\n")
models_f1 <- c("M0", "MP", "MP_cal", "ME", "ME2", "ME_P", "ME2P", "XEN")
pf_list <- list()
for (m in models_f1) {
  f <- file.path(out_dir, paste0("BRv4_perfold_", m, ".csv"))
  dt <- fread(f)
  pf_list[[m]] <- dt[, .(Model = m, ConnGroup = factor(ConnGroup, levels = 1:19), AP, LogLoss, Brier, ROC)]
}
pf_all <- rbindlist(pf_list)
pf_all[, Model := factor(Model, levels = models_f1)]

fwrite(pf_all, file.path(out_dir, "fig_v4_perfold_performance_sourcedata.csv"))

# Panel A: LogLoss (Primary Metric)
p1a <- ggplot(pf_all[Model %in% c("M0", "MP", "MP_cal", "ME_P", "ME2P")],
              aes(x = ConnGroup, y = LogLoss, fill = Model)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_brewer(palette = "Set2") +
  labs(x = "Left-out Connected Component (ConnGroup)", y = "LogLoss (lower is better)",
       title = "A. Primary Metric: Out-of-fold LogLoss across 19 Leave-Connected-Group-Out Folds",
       subtitle = "Comparing prevalence baseline (M0), raw MP, calibrated MP, MP+Target (ME_P), and MP+Target+Module (ME2P)") +
  theme_bw(base_size = 11) +
  theme(plot.margin = margin(10, 12, 10, 20), legend.position = "top", panel.grid.minor = element_blank())

# Panel B: Average Precision (AP)
p1b <- ggplot(pf_all[Model %in% c("M0", "MP", "MP_cal", "ME_P", "ME2P")],
              aes(x = ConnGroup, y = AP, fill = Model)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_hline(yintercept = 0.4015, linetype = "dashed", color = "gray40") +
  scale_fill_brewer(palette = "Set2") +
  labs(x = "Left-out Connected Component (ConnGroup)", y = "Average precision (higher is better)",
       title = "B. Secondary Metric: Threshold-Aggregated Average Precision",
       subtitle = "Dashed line = overall dataset prevalence (0.4015)") +
  theme_bw(base_size = 11) +
  theme(plot.margin = margin(10, 12, 10, 20), legend.position = "top", panel.grid.minor = element_blank())

# Save combined plot
p1_combined <- gridExtra::grid.arrange(p1a, p1b, ncol = 1)
ggsave(file.path(out_dir, "fig_v4_perfold_performance.png"), p1_combined, width = 11, height = 8.5, dpi = 300)
ggsave(file.path(out_dir, "fig_v4_perfold_performance.pdf"), p1_combined, width = 11, height = 8.5)

# --- Figure 2: Material-Level Calibration across 28 ReferenceKeys ---
cat("Generating Figure 2: Material-Level Calibration...\n")
cal_ref <- fread(file.path(out_dir, "BRv4_material_calibration_referencekey.csv"))
cal_sub <- cal_ref[Model %in% c("M0", "MP", "ME_P", "ME2P")]
cal_sub[, Model := factor(Model, levels = c("M0", "MP", "ME_P", "ME2P"))]

fwrite(cal_sub, file.path(out_dir, "fig_v4_material_calibration_sourcedata.csv"))

p2 <- ggplot(cal_sub, aes(x = mean_pred, y = obs_rate, color = Model, shape = Model)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
  geom_point(size = 2.8, alpha = 0.85) +
  scale_color_manual(values = c("M0" = "#666666", "MP" = "#2b83ba", "ME_P" = "#fdae61", "ME2P" = "#d7191c")) +
  facet_wrap(~Model, nrow = 2) +
  labs(x = "Mean Predicted Probability P(Kla detected | Ref)",
       y = "Observed Detection Rate",
       title = "B-R v4 Material-Level Calibration across 28 RNA Reference Materials",
       subtitle = "Each point represents one unique ReferenceKey; dashed diagonal indicates ideal calibration") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "none")

ggsave(file.path(out_dir, "fig_v4_material_calibration.png"), p2, width = 7.5, height = 7, dpi = 300)
ggsave(file.path(out_dir, "fig_v4_material_calibration.pdf"), p2, width = 7.5, height = 7)

# --- Figure 3: Module Incremental Effect Attribution ---
cat("Generating Figure 3: Module Incremental Effect...\n")
pf_mep <- fread(file.path(out_dir, "BRv4_perfold_ME_P.csv"))
pf_me2p <- fread(file.path(out_dir, "BRv4_perfold_ME2P.csv"))
pf_me <- fread(file.path(out_dir, "BRv4_perfold_ME.csv"))
pf_me2 <- fread(file.path(out_dir, "BRv4_perfold_ME2.csv"))

pair_mp <- data.table(
  Comparison = "ME2P vs ME_P (with MP background)",
  ConnGroup = pf_mep$ConnGroup,
  delta_LogLoss = pf_me2p$LogLoss - pf_mep$LogLoss,
  delta_AP = pf_me2p$AP - pf_mep$AP,
  delta_Brier = pf_me2p$Brier - pf_mep$Brier
)

pair_no_mp <- data.table(
  Comparison = "ME2 vs ME (without MP background)",
  ConnGroup = pf_me$ConnGroup,
  delta_LogLoss = pf_me2$LogLoss - pf_me$LogLoss,
  delta_AP = pf_me2$AP - pf_me$AP,
  delta_Brier = pf_me2$Brier - pf_me$Brier
)

inc_dt <- rbind(pair_mp, pair_no_mp)
fwrite(inc_dt, file.path(out_dir, "fig_v4_module_increment_sourcedata.csv"))

p3a <- ggplot(inc_dt, aes(x = factor(ConnGroup, levels = 1:19), y = delta_LogLoss, fill = Comparison)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.5) +
  scale_fill_manual(values = c("ME2P vs ME_P (with MP background)" = "#e41a1c",
                               "ME2 vs ME (without MP background)" = "#377eb8")) +
  labs(x = "Connected Component (ConnGroup)",
       y = expression(Delta * " LogLoss (Module Added - Baseline)"),
       title = "A. Module Incremental LogLoss Effect Across Folds",
       subtitle = "Values > 0 indicate that adding the L-lactate conversion module worsens prediction loss") +
  theme_bw(base_size = 11) +
  theme(plot.margin = margin(10, 12, 10, 20), legend.position = "top", panel.grid.minor = element_blank())

p3b <- ggplot(inc_dt, aes(x = factor(ConnGroup, levels = 1:19), y = delta_AP, fill = Comparison)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.5) +
  scale_fill_manual(values = c("ME2P vs ME_P (with MP background)" = "#e41a1c",
                               "ME2 vs ME (without MP background)" = "#377eb8")) +
  labs(x = "Connected Component (ConnGroup)",
       y = expression(Delta * " AP (Module Added - Baseline)"),
       title = "B. Module Incremental Average Precision Effect Across Folds",
       subtitle = "Values > 0 indicate improvement in ranking precision; 0 indicates neutral") +
  theme_bw(base_size = 11) +
  theme(plot.margin = margin(10, 12, 10, 20), legend.position = "top", panel.grid.minor = element_blank())

p3_combined <- gridExtra::grid.arrange(p3a, p3b, ncol = 1)
ggsave(file.path(out_dir, "fig_v4_module_incremental_effect.png"), p3_combined, width = 10, height = 8, dpi = 300)
ggsave(file.path(out_dir, "fig_v4_module_incremental_effect.pdf"), p3_combined, width = 10, height = 8)

# --- Figure 4: S1 and S2 Sensitivity Distributions ---
cat("Generating Figure 4: Sensitivity Distributions...\n")
s1_dt <- fread(file.path(out_dir, "S1_eval_sensitivity_drop_one_test_group.csv"))
s2_sum <- fread(file.path(out_dir, "S2_train_sensitivity_summary.csv"))

s1_sub <- s1_dt[Model %in% c("MP", "ME_P", "ME2P"), .(Type = "S1 (Eval Drop)", Model, dropped_group = dropped_test_group, LogLoss = eqw_LogLoss, AP = eqw_AP)]
s2_sub <- s2_sum[Model %in% c("MP", "ME_P", "ME2P"), .(Type = "S2 (Train Drop)", Model, dropped_group = dropped_train_group, LogLoss = eqw_LL_new, AP = eqw_AP_new)]

sens_all <- rbind(s1_sub, s2_sub)
sens_all[, Model := factor(Model, levels = c("MP", "ME_P", "ME2P"))]
fwrite(sens_all, file.path(out_dir, "fig_v4_sensitivity_sourcedata.csv"))

p4a <- ggplot(sens_all, aes(x = Model, y = LogLoss, color = Type)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, position = position_dodge(width = 0.6)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.6), size = 1.8, alpha = 0.7) +
  scale_color_manual(values = c("S1 (Eval Drop)" = "#2b83ba", "S2 (Train Drop)" = "#d7191c")) +
  labs(x = "Model", y = "Equal-Weight LogLoss",
       title = "A. Sensitivity of Equal-Weight LogLoss under Single-Group Deletion",
       subtitle = "S1 = evaluation aggregation drop; S2 = true training dataset group deletion") +
  theme_bw(base_size = 11) +
  theme(plot.margin = margin(10, 12, 10, 20), legend.position = "top", panel.grid.minor = element_blank())

p4b <- ggplot(sens_all, aes(x = Model, y = AP, color = Type)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, position = position_dodge(width = 0.6)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.6), size = 1.8, alpha = 0.7) +
  scale_color_manual(values = c("S1 (Eval Drop)" = "#2b83ba", "S2 (Train Drop)" = "#d7191c")) +
  labs(x = "Model", y = "Equal-Weight Average Precision (AP)",
       title = "B. Sensitivity of Equal-Weight AP under Single-Group Deletion",
       subtitle = "Evaluating stability across 19 leave-one-group-out perturbations") +
  theme_bw(base_size = 11) +
  theme(plot.margin = margin(10, 12, 10, 20), legend.position = "top", panel.grid.minor = element_blank())

p4_combined <- gridExtra::grid.arrange(p4a, p4b, ncol = 1)
ggsave(file.path(out_dir, "fig_v4_sensitivity_s1_s2.png"), p4_combined, width = 8.5, height = 7.5, dpi = 300)
ggsave(file.path(out_dir, "fig_v4_sensitivity_s1_s2.pdf"), p4_combined, width = 8.5, height = 7.5)

# Render verification
cat("Verifying rendered figure files...\n")
fig_files <- c(
  "fig_v4_perfold_performance.png", "fig_v4_perfold_performance.pdf",
  "fig_v4_material_calibration.png", "fig_v4_material_calibration.pdf",
  "fig_v4_module_incremental_effect.png", "fig_v4_module_incremental_effect.pdf",
  "fig_v4_sensitivity_s1_s2.png", "fig_v4_sensitivity_s1_s2.pdf"
)

for (ff in fig_files) {
  fpath <- file.path(out_dir, ff)
  stopifnot(file.exists(fpath))
  sz <- file.info(fpath)$size
  cat(sprintf("   -> %s: %d bytes (OK)\n", ff, sz))
  if (grepl("\\.png$", ff)) {
    stopifnot(sz > 50000)
  } else {
    stopifnot(sz > 3000)
  }
}

cat("ALL FIGURES AND SOURCE DATA WRITTEN AND VERIFIED SUCCESSFULLY.\n")
