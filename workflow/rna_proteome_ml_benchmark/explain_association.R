#!/usr/bin/env Rscript
# Descriptive two-step association: RNA -> ordinary-protein detection -> Kla detection.
# Independent RNA references are linked by material, never treated as matched patients.
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ragg)
})
set.seed(25)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: explain_association.R NEW_OUTPUT_DIR")
out <- args[[1L]]
if (dir.exists(out) && length(list.files(out, all.files = TRUE, no.. = TRUE))) {
  stop("Refusing nonempty output directory")
}
dir.create(file.path(out, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "figures"), recursive = TRUE, showWarnings = FALSE)
put <- function(x, stem) fwrite(x, file.path(out, "tables", paste0(stem, ".csv")))
inputs <- c(
  script = "workflow/rna_proteome_ml_benchmark/explain_association.R",
  labels = "outputs/20260923_rna_only_kla_fraction/protein_group_labels_31.csv",
  mapping = "config/lactate_metabolism/curated_mechanism_evidence.csv",
  rna = "outputs/20260920_lactate_metabolism/tables/material_gene_expression.csv",
  reference = "data/publication_input/reference_protein_membership_31.csv",
  kla = "data/publication_input/kla_protein_membership_31.csv"
)
stopifnot(all(file.exists(inputs)))
put(data.table(Role = names(inputs), Path = unname(inputs),
               SHA256 = vapply(inputs, digest::digest, character(1),
                               algo = "sha256", file = TRUE, USE.NAMES = FALSE)), "input_manifest")

groups <- fread(inputs[["labels"]])
map <- unique(fread(inputs[["mapping"]])[Measured == TRUE,
                                             .(Ensembl, UniProt, Symbol)])
expr <- fread(inputs[["rna"]])[, .(ReferenceKey, Ensembl,
                                    RNA = UnsmoothedLog2TPM)]
stopifnot(nrow(groups) == 31L, uniqueN(groups$ReferenceKey) == 28L,
          uniqueN(groups$ConnGroup) == 15L, nrow(map) == 51L,
          !anyDuplicated(map$Ensembl), !anyDuplicated(map$UniProt),
          all(map$Ensembl %in% expr$Ensembl))

# Stable Ensembl and UniProt identifiers are analysis keys. Symbols are labels.
groups[, Join := 1L]
map[, Join := 1L]
z <- merge(groups[, .(GroupID, PXD, SampleGroup, ReferenceKey, ConnGroup,
                      ReferenceUniqueAccessions, Join)], map,
           by = "Join", allow.cartesian = TRUE)
z[, Join := NULL]
z <- merge(z, expr, by = c("ReferenceKey", "Ensembl"))
stopifnot(nrow(z) == 31L * 51L, !anyDuplicated(z, by = c("GroupID", "Ensembl")),
          all(is.finite(z$RNA)))

ref <- fread(inputs[["reference"]])
ref <- unique(ref[nzchar(MappedBaseAccessions),
                  .(UniProt = trimws(unlist(strsplit(MappedBaseAccessions, ";", fixed = TRUE)))),
                  by = .(PXD, SampleGroup)])
ref[, ProteinDetected := 1L]
kla <- unique(fread(inputs[["kla"]])[, .(PXD, SampleGroup,
                                          UniProt = BaseAccession)])
kla[, KlaDetected := 1L]
z <- merge(z, ref, by = c("PXD", "SampleGroup", "UniProt"), all.x = TRUE)
z <- merge(z, kla, by = c("PXD", "SampleGroup", "UniProt"), all.x = TRUE)
z[is.na(ProteinDetected), ProteinDetected := 0L]
z[is.na(KlaDetected), KlaDetected := 0L]
z[, LogDepth := log1p(ReferenceUniqueAccessions)]
setorder(z, GroupID, Ensembl)
stopifnot(nrow(z) == 1581L, sum(z$ProteinDetected) == 928L,
          sum(z$KlaDetected) == 175L,
          sum(z$ProteinDetected == 0L & z$KlaDetected == 1L) == 8L)
put(z[, .(GroupID, PXD, SampleGroup, ReferenceKey, ConnGroup, Ensembl,
          UniProt, Symbol, RNA, ReferenceUniqueAccessions,
          ProteinDetected, KlaDetected)], "gene_material_evidence")

auc <- function(y, p) {
  positives <- sum(y == 1L)
  negatives <- sum(y == 0L)
  if (positives == 0L || negatives == 0L) return(NA_real_)
  (sum(rank(p, ties.method = "average")[y == 1L]) -
     positives * (positives + 1L) / 2) / (positives * negatives)
}
average_precision <- function(y, p) {
  ord <- order(-p)
  pp <- p[ord]
  yy <- y[ord]
  if (sum(yy) == 0L) return(NA_real_)
  ends <- which(c(diff(pp) != 0, TRUE))
  cumulative_positive <- cumsum(yy)[ends]
  sum(diff(c(0, cumulative_positive)) / sum(yy) *
        cumulative_positive / ends)
}
spec <- data.table(
  Model = c("Protein_depth", "Protein_RNA", "Protein_RNA_depth", "Protein_gene_depth",
            "Protein_gene_RNA_depth", "Kla_protein_depth", "Kla_RNA_protein_depth",
            "Kla_gene_protein_depth", "Kla_gene_RNA_protein_depth"),
  Outcome = c(rep("ProteinDetected", 5L), rep("KlaDetected", 4L)),
  Formula = c(
    "ProteinDetected ~ LogDepth",
    "ProteinDetected ~ RNA",
    "ProteinDetected ~ RNA + LogDepth",
    "ProteinDetected ~ factor(Ensembl) + LogDepth",
    "ProteinDetected ~ factor(Ensembl) + RNA + LogDepth",
    "KlaDetected ~ ProteinDetected + LogDepth",
    "KlaDetected ~ ProteinDetected + RNA + LogDepth",
    "KlaDetected ~ factor(Ensembl) + ProteinDetected + LogDepth",
    "KlaDetected ~ factor(Ensembl) + ProteinDetected + RNA + LogDepth"
  )
)
put(spec, "model_contract")

fitted_rows <- list()
metric_rows <- list()
fold_rows <- list()
for (i in seq_len(nrow(spec))) {
  model <- spec$Model[i]
  outcome <- spec$Outcome[i]
  response <- z[[outcome]]
  formula <- as.formula(spec$Formula[i])
  fit <- glm(formula, family = binomial(), data = z,
             control = glm.control(maxit = 100L))
  stopifnot(fit$converged)
  in_sample <- as.vector(predict(fit, type = "response"))
  heldout <- numeric(nrow(z))
  for (held in sort(unique(z$ConnGroup))) {
    tr <- z$ConnGroup != held
    te <- !tr
    local <- glm(formula, family = binomial(), data = z[tr],
                 control = glm.control(maxit = 100L))
    stopifnot(local$converged)
    heldout[te] <- as.vector(predict(local, newdata = z[te], type = "response"))
    fold_rows[[length(fold_rows) + 1L]] <- data.table(
      Model = model, ConnGroup = held, N = sum(te), Positives = sum(response[te]),
      AUC = auc(response[te], heldout[te]),
      AveragePrecision = average_precision(response[te], heldout[te]),
      Brier = mean((response[te] - heldout[te])^2))
  }
  stopifnot(all(is.finite(in_sample)), all(is.finite(heldout)),
            all(in_sample >= 0 & in_sample <= 1),
            all(heldout >= 0 & heldout <= 1))
  fitted_rows[[length(fitted_rows) + 1L]] <- data.table(
    Model = model, GroupID = z$GroupID, Ensembl = z$Ensembl,
    ConnGroup = z$ConnGroup, Observed = response,
    InSample = in_sample, Heldout = heldout)
  metric_rows[[length(metric_rows) + 1L]] <- data.table(
    Model = model, Outcome = outcome, N = length(response),
    Positives = sum(response), Prevalence = mean(response),
    InSampleAUC = auc(response, in_sample), HeldoutAUC = auc(response, heldout),
    InSampleAP = average_precision(response, in_sample),
    HeldoutAP = average_precision(response, heldout),
    InSampleBrier = mean((response - in_sample)^2),
    HeldoutBrier = mean((response - heldout)^2),
    RNACoefficient = if ("RNA" %in% names(coef(fit))) unname(coef(fit)["RNA"]) else NA_real_,
    ProteinCoefficient = if ("ProteinDetected" %in% names(coef(fit)))
      unname(coef(fit)["ProteinDetected"]) else NA_real_)
}
fitted <- rbindlist(fitted_rows)
metrics <- rbindlist(metric_rows)
folds <- rbindlist(fold_rows)
put(fitted, "fitted_probabilities")
put(metrics, "model_metrics")
put(folds, "source_fold_metrics")

# Training-data visualization is deliberately labelled as descriptive fitting.
main <- metrics[Model %in% c("Protein_RNA_depth", "Kla_RNA_protein_depth")]
stage <- data.table(
  Model = c("Protein_RNA_depth", "Kla_RNA_protein_depth"),
  Label = c("RNA to ordinary-protein detection",
            "RNA + ordinary protein to Kla detection")
)
main <- merge(main, stage, by = "Model")
main[, PlotLabel := sprintf("%s\nfit AUC %.2f | held-out AUC %.2f",
                            Label, InSampleAUC, HeldoutAUC)]
figure_data <- merge(fitted[Model %in% main$Model],
                     main[, .(Model, PlotLabel)], by = "Model")
figure_data[, PlotLabel := factor(PlotLabel,
                                 levels = main$PlotLabel[match(stage$Model, main$Model)])]
roc_points <- rbindlist(lapply(unique(figure_data$Model), function(model) {
  x <- figure_data[Model == model]
  yy <- x$Observed
  pp <- x$Heldout
  ord <- order(-pp, seq_along(pp))
  ends <- which(c(diff(pp[ord]) != 0, TRUE))
  data.table(Model = model, PlotLabel = x$PlotLabel[1L],
             FPR = c(0, cumsum(yy[ord] == 0L)[ends] / sum(yy == 0L)),
             TPR = c(0, cumsum(yy[ord] == 1L)[ends] / sum(yy == 1L)))
}))
p <- ggplot(roc_points, aes(FPR, TPR)) +
  geom_abline(intercept = 0, slope = 1, colour = "grey72") +
  geom_line(linewidth = 1.1, colour = "#2166AC") +
  facet_wrap(~PlotLabel, nrow = 1) +
  coord_equal() +
  labs(title = "A two-step link across RNA and proteomics",
       subtitle = "51 prespecified genes × 31 materials; models describe detection in these experiments",
       x = "False-positive rate", y = "True-positive rate",
       caption = "ROC curves use source-held-out predictions. The fit AUC is in-sample only.\nRNA references are not from the proteomics patients; detection is not modification occupancy.") +
  theme_minimal(base_size = 10, base_family = "Arial") +
  theme(plot.title = element_text(face = "bold"), strip.text = element_text(face = "bold"),
        plot.caption = element_text(hjust = 0), panel.grid.minor = element_blank())
ggsave(file.path(out, "figures", "two_step_association_ROC.pdf"), p,
       width = 10, height = 5.4, device = cairo_pdf)
ggsave(file.path(out, "figures", "two_step_association_ROC.png"), p,
       width = 10, height = 5.4, device = agg_png, dpi = 180)

# Binned observed-vs-fitted relationships make the model legible to readers.
z[, RNABin := cut(rank(RNA, ties.method = "first"), breaks = 10L,
                 labels = sprintf("%02d", 1:10))]
ordinary <- merge(z[, .(GroupID, Ensembl, RNA, RNABin, ProteinDetected)],
                  fitted[Model == "Protein_RNA_depth", .(GroupID, Ensembl, Fitted = InSample)],
                  by = c("GroupID", "Ensembl"))
obs1 <- ordinary[, .(MeanRNA = mean(RNA), Observed = mean(ProteinDetected),
                     Fitted = mean(Fitted), N = .N), by = RNABin]
obs1[, Stage := "RNA → ordinary protein"]
modified <- merge(z[ProteinDetected == 1L,
                    .(GroupID, Ensembl, RNA, RNABin, KlaDetected)],
                  fitted[Model == "Kla_RNA_protein_depth",
                         .(GroupID, Ensembl, Fitted = InSample)],
                  by = c("GroupID", "Ensembl"))
modified[, RNABin := cut(rank(RNA, ties.method = "first"), breaks = 10L,
                         labels = sprintf("%02d", 1:10))]
obs2 <- modified[, .(MeanRNA = mean(RNA), Observed = mean(KlaDetected),
                     Fitted = mean(Fitted), N = .N), by = RNABin]
obs2[, Stage := "RNA + ordinary protein → Kla"]
bins <- rbindlist(list(obs1, obs2))
put(bins, "binned_observed_and_fit")
long <- melt(bins, id.vars = c("RNABin", "MeanRNA", "N", "Stage"),
             measure.vars = c("Observed", "Fitted"),
             variable.name = "Series", value.name = "DetectionRate")
q <- ggplot(long, aes(MeanRNA, DetectionRate, colour = Series, group = Series)) +
  geom_line(linewidth = 0.85) + geom_point(size = 2) +
  facet_wrap(~Stage, nrow = 1, scales = "free_y") +
  scale_colour_manual(values = c(Observed = "#2166AC", Fitted = "#B2182B")) +
  scale_y_continuous(labels = function(x) paste0(round(100 * x), "%")) +
  labs(title = "Higher RNA expression tracks protein detection",
       subtitle = "Descriptive training-data fit across 51 stable-ID genes and 31 material groups",
       x = "Mean unsmoothed log2(TPM + 0.5) within each stage's RNA decile",
       y = "Fraction detected", colour = NULL,
       caption = "Left: ordinary-proteome detection. Right: Kla detection among ordinary-proteome-positive gene/material pairs.\nExternal RNA is linked by material, not patient; fitted values describe this assembled dataset only.") +
  theme_minimal(base_size = 10, base_family = "Arial") +
  theme(plot.title = element_text(face = "bold"), legend.position = "bottom",
        panel.grid.minor = element_blank(), plot.caption = element_text(hjust = 0))
ggsave(file.path(out, "figures", "RNA_protein_Kla_association.pdf"), q,
       width = 11, height = 5.5, device = cairo_pdf)
ggsave(file.path(out, "figures", "RNA_protein_Kla_association.png"), q,
       width = 11, height = 5.5, device = agg_png, dpi = 180)

files <- list.files(out, recursive = TRUE, full.names = TRUE)
files <- files[file.info(files)$isdir == FALSE]
fwrite(data.table(File = substring(files, nchar(out) + 2L),
                  SHA256 = vapply(files, digest::digest, character(1),
                                  algo = "sha256", file = TRUE, USE.NAMES = FALSE)),
       file.path(out, "release_sha256.csv"))
cat("PASS: 51 genes × 31 materials, 15 source blocks; two-stage detection association.\n")
