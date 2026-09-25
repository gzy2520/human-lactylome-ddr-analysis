#!/usr/bin/env Rscript
# Fixed comparison of RNA, whole-proteome detection, and their combination.
# The label is a cross-assay detection percentage, not Kla occupancy.
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(ggplot2)
  library(ragg)
})
set.seed(25)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: run.R NEW_OUTPUT_DIR")
out <- args[[1L]]
if (dir.exists(out) && length(list.files(out, all.files = TRUE, no.. = TRUE))) {
  stop("Refusing nonempty output directory")
}
dir.create(file.path(out, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "figures"), recursive = TRUE, showWarnings = FALSE)
put <- function(x, stem) fwrite(x, file.path(out, "tables", paste0(stem, ".csv")))
inputs <- c(
  script = "workflow/rna_proteome_ml_benchmark/run.R",
  labels = "outputs/20260923_rna_only_kla_fraction/protein_group_labels_31.csv",
  reference = "data/publication_input/reference_protein_membership_31.csv",
  rna = "outputs/20260920_lactate_metabolism/tables/material_gene_expression.csv",
  modules = "config/lactate_metabolism/gene_module_membership.csv",
  catalog = "outputs/20260923_rna_only_kla_fraction_optimization/module_catalog.csv"
)
if (!all(file.exists(inputs))) stop("Missing input: ", paste(inputs[!file.exists(inputs)], collapse = ", "))
put(data.table(Role = names(inputs), Path = unname(inputs),
               SHA256 = vapply(inputs, digest::digest, character(1),
                               algo = "sha256", file = TRUE, USE.NAMES = FALSE)), "input_manifest")

d <- fread(inputs[["labels"]])
setorder(d, GroupID)
stopifnot(nrow(d) == 31L, !anyDuplicated(d$GroupID), uniqueN(d$ReferenceKey) == 28L,
          uniqueN(d$ConnGroup) == 15L, all(is.finite(d$DetectedKlaPercent)))
y <- d$DetectedKlaPercent
block <- d$ConnGroup

# Use every stable, mapped accession in the ordinary-proteome denominator.
# Multiple accessions from an ambiguous protein group are split exactly as in
# the frozen label definition. Blank, unmapped ENSP rows contribute no feature.
ref <- fread(inputs[["reference"]])
ref <- merge(ref[, .(PXD, SampleGroup, MappedBaseAccessions)],
             d[, .(PXD, SampleGroup, GroupID)], by = c("PXD", "SampleGroup"),
             all.x = TRUE, sort = FALSE, allow.cartesian = TRUE)
stopifnot(nrow(ref) == 188218L, !anyNA(ref$GroupID))
parts <- strsplit(ref$MappedBaseAccessions, ";", fixed = TRUE)
acc <- unique(data.table(GroupID = rep(ref$GroupID, lengths(parts)),
                         Accession = trimws(unlist(parts, use.names = FALSE))))
stopifnot(all(nzchar(acc$Accession)), !anyDuplicated(acc, by = c("GroupID", "Accession")))
proteins <- sort(unique(acc$Accession))
P <- sparseMatrix(i = match(acc$GroupID, d$GroupID),
                  j = match(acc$Accession, proteins), x = 1,
                  dims = c(nrow(d), length(proteins)))
depth <- as.numeric(rowSums(P))
stopifnot(identical(as.integer(depth), d$ReferenceUniqueAccessions))
raw_gram <- as.matrix(tcrossprod(P))
log_depth <- log1p(depth)

# The 16 RNA modules were fixed in the preceding RNA-only optimization.
rna <- fread(inputs[["rna"]])
stopifnot(uniqueN(rna$Ensembl) == 128L,
          uniqueN(rna$ReferenceKey) == 28L,
          !anyDuplicated(rna, by = c("ReferenceKey", "Ensembl")))
genes <- sort(unique(rna$Ensembl))
wide <- dcast(rna[, .(ReferenceKey, Ensembl, Value = UnsmoothedLog2TPM)],
              ReferenceKey ~ Ensembl, value.var = "Value")
stopifnot(!anyNA(wide), identical(names(wide)[-1L], genes))
X <- as.matrix(wide[match(d$ReferenceKey, wide$ReferenceKey), ..genes])
storage.mode(X) <- "double"
stopifnot(all(is.finite(X)))
mem <- unique(fread(inputs[["modules"]])[Ensembl %chin% genes,
                                          .(Origin, Module, Ensembl)])
catalog <- fread(inputs[["catalog"]])
mem <- merge(mem, catalog[, .(Origin, Module, FeatureID)],
             by = c("Origin", "Module"))
modules <- split(mem$Ensembl, mem$FeatureID)
stopifnot(nrow(catalog) == 16L, setequal(names(modules), catalog$FeatureID),
          all(lengths(modules) > 0L),
          identical(as.integer(lengths(modules[catalog$FeatureID])), catalog$MeasuredGenes))
modules <- lapply(modules[catalog$FeatureID], function(ids) match(ids, genes))

# All transformations below are trained on the current training sources only.
rna_kernel <- function(tr, te) {
  center <- colMeans(X[tr, , drop = FALSE])
  spread <- apply(X[tr, , drop = FALSE], 2L, sd)
  spread[!is.finite(spread) | spread == 0] <- 1
  ztr <- sweep(sweep(X[tr, , drop = FALSE], 2L, center, "-"), 2L, spread, "/")
  zte <- sweep(sweep(X[te, , drop = FALSE], 2L, center, "-"), 2L, spread, "/")
  score <- function(z) matrix(vapply(modules, function(ids)
    rowMeans(z[, ids, drop = FALSE]), numeric(nrow(z))),
    nrow = nrow(z), ncol = length(modules),
    dimnames = list(NULL, names(modules)))
  a <- score(ztr)
  b <- score(zte)
  mc <- colMeans(a)
  ms <- apply(a, 2L, sd)
  ms[!is.finite(ms) | ms == 0] <- 1
  a <- sweep(sweep(a, 2L, mc, "-"), 2L, ms, "/") / sqrt(ncol(a))
  b <- sweep(sweep(b, 2L, mc, "-"), 2L, ms, "/") / sqrt(ncol(b))
  list(train = tcrossprod(a), test = b %*% t(a))
}
profile_kernel <- function(tr, te) {
  a <- raw_gram[tr, tr, drop = FALSE]
  m <- rowMeans(a)
  grand <- mean(a)
  kt <- sweep(sweep(a, 1L, m, "-"), 2L, m, "-") + grand
  b <- raw_gram[te, tr, drop = FALSE]
  kb <- sweep(sweep(b, 1L, rowMeans(b), "-"), 2L, m, "-") + grand
  scale <- mean(diag(kt))
  stopifnot(is.finite(scale), scale > 0)
  list(train = kt / scale, test = kb / scale)
}
depth_kernel <- function(tr, te) {
  mu <- mean(log_depth[tr])
  s <- sd(log_depth[tr])
  if (!is.finite(s) || s == 0) s <- 1
  a <- (log_depth[tr] - mu) / s
  b <- (log_depth[te] - mu) / s
  list(train = tcrossprod(matrix(a, ncol = 1L)),
       test = tcrossprod(matrix(b, ncol = 1L), matrix(a, ncol = 1L)))
}
kernels <- function(tr, te) {
  R <- rna_kernel(tr, te)
  H <- profile_kernel(tr, te)
  D <- depth_kernel(tr, te)
  combine <- function(a, b) list(train = (a$train + b$train) / 2,
                                  test = (a$test + b$test) / 2)
  prot <- combine(H, D)
  list(RNA16 = R, DepthOnly = D, WholeProteome = prot,
       RNAplusProteome = combine(R, prot))
}
predict_krr <- function(tr, kernel, lambda) {
  mu <- mean(y[tr])
  alpha <- solve(kernel$train + diag(lambda, length(tr)), y[tr] - mu)
  pmin(100, pmax(0, as.vector(mu + kernel$test %*% alpha)))
}
weighted_median <- function(values, sources) {
  w <- 1 / as.vector(table(sources)[as.character(sources)])
  ix <- order(values)
  values[ix][which(cumsum(w[ix]) >= sum(w) / 2)[1L]]
}

# Prespecified simple linear ridge models. Source-blocked inner selection of
# penalty; complete connected sources are held out in the outer evaluation.
models <- c("RNA16", "DepthOnly", "WholeProteome", "RNAplusProteome")
lambdas <- c(100, 30, 10, 3, 1, 0.3, 0.1)
preds <- list()
selected <- list()
inner_log <- list()
for (outer in sort(unique(block))) {
  tr <- which(block != outer)
  te <- which(block == outer)
  held_sources <- sort(unique(block[tr]))
  loss <- CJ(Model = models, Lambda = lambdas)
  loss[, LossSum := 0.0]
  for (inner in held_sources) {
    tr2 <- tr[block[tr] != inner]
    te2 <- tr[block[tr] == inner]
    ks <- kernels(tr2, te2)
    for (row in seq_len(nrow(loss))) {
      p <- predict_krr(tr2, ks[[loss$Model[row]]], loss$Lambda[row])
      loss$LossSum[row] <- loss$LossSum[row] + mean(abs(p - y[te2]))
    }
  }
  loss[, `:=`(InnerSourceMAE = LossSum / length(held_sources),
              OuterConnGroup = outer, LossSum = NULL)]
  inner_log[[length(inner_log) + 1L]] <- copy(loss)
  best <- loss[order(Model, InnerSourceMAE, -Lambda), .SD[1L], by = Model]
  selected[[length(selected) + 1L]] <- best
  ks <- kernels(tr, te)
  row <- d[te, .(GroupID, PXD, SampleGroup, ReferenceKey, ConnGroup,
                 ObservedPercent = DetectedKlaPercent, ReferenceUniqueAccessions)]
  row[, BalancedMedian := weighted_median(y[tr], block[tr])]
  for (i in seq_len(nrow(best))) {
    model <- best$Model[i]
    row[, (model) := predict_krr(tr, ks[[model]], best$Lambda[i])]
  }
  preds[[length(preds) + 1L]] <- row
}
pred <- rbindlist(preds)
selection <- rbindlist(selected)
inner_log <- rbindlist(inner_log)
setorder(pred, GroupID)
stopifnot(nrow(pred) == 31L, !anyDuplicated(pred$GroupID),
          identical(pred$GroupID, d$GroupID),
          !anyNA(pred[, c("ObservedPercent", "BalancedMedian", models), with = FALSE]),
          nrow(selection) == 15L * length(models))
put(pred, "heldout_predictions_31")
put(selection[, .(OuterConnGroup, Model, Lambda, InnerSourceMAE)], "selected_penalties")
put(inner_log[, .(OuterConnGroup, Model, Lambda, InnerSourceMAE)], "inner_source_losses")
put(data.table(GroupID = d$GroupID, ReferenceKey = d$ReferenceKey,
               ConnGroup = block, ReferenceUniqueAccessions = depth), "feature_inventory")

all_models <- c("BalancedMedian", models)
performance <- rbindlist(lapply(all_models, function(model) {
  errors <- pred[[model]] - y
  by_source <- pred[, .(MAE = mean(abs(get(model) - ObservedPercent))), by = ConnGroup]
  data.table(Model = model, GroupMAE = mean(abs(errors)),
             SourceMeanMAE = mean(by_source$MAE), GroupRMSE = sqrt(mean(errors^2)),
             BetterThanBaselineGroups = sum(abs(errors) < abs(pred$BalancedMedian - y)),
             BetterThanBaselineSources = sum(by_source$MAE < pred[, .(
               MAE = mean(abs(BalancedMedian - ObservedPercent))), by = ConnGroup]$MAE))
}))
put(performance, "performance")
source_error <- pred[, lapply(.SD, function(v) mean(abs(v - ObservedPercent))),
                     by = ConnGroup, .SDcols = all_models]
setorder(source_error, ConnGroup)
source_error[, JointMinusProteome := RNAplusProteome - WholeProteome]
put(source_error, "source_errors")

plot_pred <- melt(pred[, c("GroupID", "ConnGroup", "ObservedPercent", all_models), with = FALSE],
                  id.vars = c("GroupID", "ConnGroup", "ObservedPercent"),
                  variable.name = "Model", value.name = "PredictedPercent")
plot_pred[, Model := factor(Model, levels = all_models,
                            labels = c("Median baseline", "RNA: 16 modules", "Proteome depth",
                                       "Whole-proteome detection", "RNA + whole proteome"))]
p <- ggplot(plot_pred, aes(ObservedPercent, PredictedPercent)) +
  geom_abline(intercept = 0, slope = 1, colour = "grey65") +
  geom_point(colour = "#2166AC", size = 1.7, alpha = 0.8) +
  facet_wrap(~Model, ncol = 3) +
  coord_equal(xlim = c(0, 40), ylim = c(0, 40)) +
  labs(title = "Held-out prediction of the detected Kla-protein percentage",
       subtitle = "31 proteomic groups | 15 connected source blocks | outer leave-one-source-out",
       x = "Observed detection percentage", y = "Predicted detection percentage",
       caption = "Percentage = 100 × detected Kla proteins also in the ordinary proteome / ordinary-proteome proteins.\nExternal RNA is tissue-matched, not donor-matched; this is not site occupancy.") +
  theme_minimal(base_size = 10, base_family = "Arial") +
  theme(plot.title = element_text(face = "bold"), panel.grid.minor = element_blank(),
        plot.caption = element_text(hjust = 0))
ggsave(file.path(out, "figures", "heldout_predictions.pdf"), p,
       width = 10, height = 7.8, device = cairo_pdf)
ggsave(file.path(out, "figures", "heldout_predictions.png"), p,
       width = 10, height = 7.8, device = agg_png, dpi = 180)

se <- copy(source_error)
se[, ConnGroup := factor(ConnGroup, levels = as.character(ConnGroup[order(JointMinusProteome)]))]
q <- ggplot(se, aes(ConnGroup, JointMinusProteome, fill = JointMinusProteome < 0)) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_col() +
  scale_fill_manual(values = c(`FALSE` = "#B2182B", `TRUE` = "#2166AC"), guide = "none") +
  labs(title = "Does RNA improve on the whole-proteome model?",
       subtitle = "Held-out mean absolute error within each independent source block",
       x = "Source block", y = "Joint minus proteome-only MAE (percentage points)",
       caption = "Below zero favours adding RNA. Blocks differ in size; no patient-matched RNA is available.") +
  theme_minimal(base_size = 10, base_family = "Arial") +
  theme(plot.title = element_text(face = "bold"), panel.grid.minor = element_blank(),
        plot.caption = element_text(hjust = 0))
ggsave(file.path(out, "figures", "added_RNA_by_source.pdf"), q,
       width = 9, height = 5.6, device = cairo_pdf)
ggsave(file.path(out, "figures", "added_RNA_by_source.png"), q,
       width = 9, height = 5.6, device = agg_png, dpi = 180)

files <- list.files(out, recursive = TRUE, full.names = TRUE)
files <- files[file.info(files)$isdir == FALSE]
fwrite(data.table(File = substring(files, nchar(out) + 2L),
                  SHA256 = vapply(files, digest::digest, character(1),
                                  algo = "sha256", file = TRUE, USE.NAMES = FALSE)),
       file.path(out, "release_sha256.csv"))
cat("PASS: 31 groups, 28 RNA references, 15 independent source blocks, ",
    length(proteins), " whole-proteome accession features.\n", sep = "")
