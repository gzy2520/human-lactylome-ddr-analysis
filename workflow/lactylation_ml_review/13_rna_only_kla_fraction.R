#!/usr/bin/env Rscript
# Exploratory prediction of the material-level detected Kla-protein percentage.
# Predictors are frozen RNA measurements only. Protein data are used only to
# construct the training/evaluation label; neither proteome depth nor group ID
# enters the prediction function. All splits and tuning use connected sources.
suppressPackageStartupMessages(library(data.table))
set.seed(25)

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
out <- args[[1L]]
if (dir.exists(out) && length(list.files(out, all.files = TRUE, no.. = TRUE)))
  stop("Refusing to overwrite a nonempty output directory")

inputs <- c(
  groups = "data/publication_input/group_summary_31.csv",
  kla = "data/publication_input/kla_protein_membership_31.csv",
  reference = "data/publication_input/reference_protein_membership_31.csv",
  mapping = "outputs/20260920_lactate_metabolism/tables/rna_proteome_mapping_31.csv",
  rna = "outputs/20260920_lactate_metabolism/tables/material_gene_expression.csv"
)
stopifnot(all(file.exists(inputs)))
g <- fread(inputs[["groups"]])
k <- fread(inputs[["kla"]])
r <- fread(inputs[["reference"]])
m <- fread(inputs[["mapping"]])
e <- fread(inputs[["rna"]])

# One observation per unique accession, including every accession in a mapped
# protein group. This keeps numerator and denominator on the same stable-ID
# universe and prevents multi-accession groups from inflating the denominator.
stopifnot(nrow(g) == 31L, nrow(m) == 31L,
          !anyDuplicated(m$GroupID),
          !anyDuplicated(g[, .(PXD, SampleGroup)]),
          !anyDuplicated(k[, .(PXD, SampleGroup, BaseAccession)]),
          nrow(k) == sum(g$KlaProteinCount),
          nrow(r) == sum(g$ReferenceProteinCount),
          !anyNA(k$BaseAccession), !anyNA(r$MappedBaseAccessions))
r[, RowID := .I]
ra <- r[, .(BaseAccession = trimws(unlist(strsplit(MappedBaseAccessions, ";", fixed = TRUE)))),
        by = .(RowID, PXD, SampleGroup)]
stopifnot(!anyNA(ra$BaseAccession), all(nzchar(ra$BaseAccession)))
ra <- unique(ra[, .(PXD, SampleGroup, BaseAccession)])
ka <- unique(k[, .(PXD, SampleGroup, BaseAccession)])
common <- merge(ra, ka, by = c("PXD", "SampleGroup", "BaseAccession"))
nr <- ra[, .(ReferenceUniqueAccessions = uniqueN(BaseAccession)),
         by = .(PXD, SampleGroup)]
nk <- ka[, .(KlaUniqueAccessions = uniqueN(BaseAccession)),
         by = .(PXD, SampleGroup)]
ni <- common[, .(MatchedUniqueAccessions = uniqueN(BaseAccession)),
             by = .(PXD, SampleGroup)]
target <- merge(merge(nr, nk, by = c("PXD", "SampleGroup")),
                ni, by = c("PXD", "SampleGroup"), all.x = TRUE)
target[is.na(MatchedUniqueAccessions), MatchedUniqueAccessions := 0L]
target[, DetectedKlaPercent := 100 * MatchedUniqueAccessions / ReferenceUniqueAccessions]
target[, KlaOutsideReference := KlaUniqueAccessions - MatchedUniqueAccessions]
target <- merge(m[, .(GroupID, ReferenceKey, PXD, SampleGroup, Category,
                      KlaLabelEn, RNALabel)], target, by = c("PXD", "SampleGroup"))
target <- merge(target,
                g[, .(PXD, SampleGroup, ReferencePXD, KlaProteinCount, ReferenceProteinCount,
                      ReferenceMatchNote)], by = c("PXD", "SampleGroup"))
setorder(target, GroupID)

# Study-level leakage can pass through either the Kla project, the reference
# proteome project, or a shared RNA profile. Build the transitive closure of
# all three links. The older DDR-only 19-block split did not merge projects
# sharing a reference-proteome study, so it is not reused here.
n <- nrow(target)
parent <- seq_len(n)
root <- function(i) {
  while (parent[[i]] != i) i <- parent[[i]]
  i
}
for (i in seq_len(n - 1L)) for (j in (i + 1L):n) {
  studies_i <- c(target$PXD[[i]], target$ReferencePXD[[i]])
  studies_j <- c(target$PXD[[j]], target$ReferencePXD[[j]])
  if (target$ReferenceKey[[i]] == target$ReferenceKey[[j]] ||
      length(intersect(studies_i, studies_j)) > 0L) {
    parent[[root(j)]] <- root(i)
  }
}
roots <- vapply(seq_len(n), root, integer(1))
target[, ConnGroup := match(roots, unique(roots))]
stopifnot(nrow(target) == 31L, uniqueN(target$ReferenceKey) == 28L,
          uniqueN(target$ConnGroup) == 15L,
          all(target$MatchedUniqueAccessions <= target$KlaUniqueAccessions),
          all(target$MatchedUniqueAccessions <= target$ReferenceUniqueAccessions),
          all(target$DetectedKlaPercent >= 0 & target$DetectedKlaPercent <= 100),
          target[, all(uniqueN(ConnGroup) == 1L), by = ReferenceKey][, all(V1)])
for (col in c("PXD", "ReferenceKey"))
  stopifnot(target[, all(uniqueN(ConnGroup) == 1L), by = col][, all(V1)])
stopifnot(target[, all(uniqueN(ConnGroup) == 1L), by = ReferencePXD][, all(V1)])

# The 31 proteomic rows contain repeated RNA materials. Aggregate repeated
# protein labels within a material before fitting; report their spread and
# retain each original group label for a secondary group-level error check.
material <- target[, .(
  ObservedPercent = mean(DetectedKlaPercent),
  MinimumGroupPercent = min(DetectedKlaPercent),
  MaximumGroupPercent = max(DetectedKlaPercent),
  ProteomeGroups = .N,
  GroupIDs = paste(sort(GroupID), collapse = ";"),
  ConnGroup = unique(ConnGroup),
  Category = unique(Category)
), by = ReferenceKey]
stopifnot(nrow(material) == 28L, all(is.finite(material$ObservedPercent)))
setorder(material, ReferenceKey)

# The earlier GO/metabolism panel was fixed before this prediction task.
# Unsmoothed log2(TPM + 0.5) avoids across-material smoothing leakage.
stopifnot(!anyDuplicated(e[, .(ReferenceKey, Ensembl)]),
          uniqueN(e$ReferenceKey) == 28L, uniqueN(e$Ensembl) == 128L,
          nrow(e) == 28L * 128L,
          all(is.finite(e$UnsmoothedLog2TPM)),
          setequal(e$ReferenceKey, material$ReferenceKey))
genes <- sort(unique(e$Ensembl))
features <- dcast(e[, .(ReferenceKey, Ensembl, Value = UnsmoothedLog2TPM)],
                  ReferenceKey ~ Ensembl, value.var = "Value")
setorder(features, ReferenceKey)
stopifnot(identical(features$ReferenceKey, material$ReferenceKey),
          identical(names(features)[-1L], genes), !anyNA(features))
x <- as.matrix(features[, ..genes])
storage.mode(x) <- "double"
y <- material$ObservedPercent
components <- material$ConnGroup

# Three RNA-only candidate sets, declared independently of the labels. Gene
# symbols appear only here as annotations; all actual keys are Ensembl IDs.
sets <- list(
  LDHA_LDHB = c("ENSG00000134333", "ENSG00000111716"),
  L_conversion_4 = c("ENSG00000134333", "ENSG00000111716",
                     "ENSG00000166796", "ENSG00000171989"),
  GO_metabolism_128 = genes
)
stopifnot(all(vapply(sets, function(ids) all(ids %in% genes), logical(1))))
lambda_grid <- c(100, 30, 10, 3, 1, 0.3, 0.1)
candidate <- rbindlist(lapply(names(sets), function(nm)
  data.table(Model = nm, Lambda = lambda_grid)))

predict_ridge <- function(train, test, feature_ids, lambda) {
  xt <- x[train, feature_ids, drop = FALSE]
  xv <- x[test, feature_ids, drop = FALSE]
  center <- colMeans(xt)
  scale <- apply(xt, 2L, sd)
  scale[!is.finite(scale) | scale == 0] <- 1
  xt <- sweep(sweep(xt, 2L, center, "-"), 2L, scale, "/") / sqrt(ncol(xt))
  xv <- sweep(sweep(xv, 2L, center, "-"), 2L, scale, "/") / sqrt(ncol(xv))
  ym <- mean(y[train])
  alpha <- solve(tcrossprod(xt) + diag(lambda, nrow(xt)), y[train] - ym)
  pmin(100, pmax(0, as.vector(ym + xv %*% crossprod(xt, alpha))))
}

preds <- list()
tuning <- list()
selected <- list()
for (outer in sort(unique(components))) {
  train <- which(components != outer)
  test <- which(components == outer)
  inner_groups <- sort(unique(components[train]))
  losses <- candidate[, {
    ids <- sets[[Model]]
    inner_mae <- vapply(inner_groups, function(inner) {
      tr <- train[components[train] != inner]
      va <- train[components[train] == inner]
      mean(abs(predict_ridge(tr, va, ids, Lambda) - y[va]))
    }, numeric(1))
    .(InnerMeanGroupMAE = mean(inner_mae))
  }, by = .(Model, Lambda)]
  losses[, OuterConnGroup := outer]
  tuning[[length(tuning) + 1L]] <- losses
  # Candidate and lambda order break exact ties toward simpler models and
  # stronger regularization. Outer labels never enter the selection.
  setorder(losses, InnerMeanGroupMAE)
  best <- losses[1L]
  selected[[length(selected) + 1L]] <- best[, .(OuterConnGroup, Model,
                                               Lambda, InnerMeanGroupMAE)]
  row <- material[test, .(ReferenceKey, ConnGroup, ObservedPercent,
                           ProteomeGroups)]
  row[, Baseline := mean(y[train])]
  row[, MedianBaseline := median(y[train])]
  for (nm in names(sets)) {
    local <- losses[Model == nm][order(InnerMeanGroupMAE)][1L]
    row[, (nm) := predict_ridge(train, test, sets[[nm]], local$Lambda)]
  }
  row[, SelectedRNA := predict_ridge(train, test, sets[[best$Model]],
                                     best$Lambda)]
  preds[[length(preds) + 1L]] <- row
}
pred <- rbindlist(preds)
setorder(pred, ReferenceKey)
stopifnot(nrow(pred) == 28L, !anyNA(pred),
          identical(pred$ReferenceKey, material$ReferenceKey),
          all(pred$SelectedRNA >= 0 & pred$SelectedRNA <= 100))
group_pred <- merge(target, pred[, c("ReferenceKey", "Baseline", names(sets),
                                     "MedianBaseline", "SelectedRNA"), with = FALSE], by = "ReferenceKey")
setorder(group_pred, GroupID)

model_cols <- c("Baseline", "MedianBaseline", names(sets), "SelectedRNA")
metrics <- rbindlist(lapply(model_cols, function(nm) {
  err <- pred[[nm]] - pred$ObservedPercent
  ge <- group_pred[[nm]] - group_pred$DetectedKlaPercent
  per_comp <- pred[, .(MAE = mean(abs(get(nm) - ObservedPercent))),
                   by = ConnGroup]
  data.table(Model = nm,
             MaterialMAE_pp = mean(abs(err)),
             MaterialRMSE_pp = sqrt(mean(err^2)),
             ConnectedGroupMeanMAE_pp = mean(per_comp$MAE),
             Original31GroupMAE_pp = mean(abs(ge)),
             MaterialMeanBias_pp = mean(err))
}))
by_comp <- pred[, .(
  Materials = .N,
  BaselineMAE_pp = mean(abs(Baseline - ObservedPercent)),
  MedianBaselineMAE_pp = mean(abs(MedianBaseline - ObservedPercent)),
  SelectedRNAMAE_pp = mean(abs(SelectedRNA - ObservedPercent)),
  SelectedMinusBaseline_pp = mean(abs(SelectedRNA - ObservedPercent)) -
    mean(abs(Baseline - ObservedPercent)),
  SelectedMinusMedian_pp = mean(abs(SelectedRNA - ObservedPercent)) -
    mean(abs(MedianBaseline - ObservedPercent))
), by = ConnGroup]
setorder(by_comp, ConnGroup)

dir.create(out, recursive = TRUE, showWarnings = FALSE)
fwrite(target, file.path(out, "protein_group_labels_31.csv"))
fwrite(target[, .(GroupID, PXD, ReferencePXD, ReferenceKey, ConnGroup)],
       file.path(out, "source_components_31.csv"))
fwrite(material, file.path(out, "material_labels_28.csv"))
fwrite(pred, file.path(out, "material_oof_predictions_28.csv"))
fwrite(group_pred, file.path(out, "protein_group_oof_predictions_31.csv"))
fwrite(metrics, file.path(out, "performance_summary.csv"))
fwrite(by_comp, file.path(out, "per_connected_group.csv"))
fwrite(rbindlist(tuning), file.path(out, "inner_tuning.csv"))
fwrite(rbindlist(selected), file.path(out, "selected_parameters.csv"))
fwrite(data.table(File = unname(inputs), MD5 = unname(tools::md5sum(inputs))),
       file.path(out, "input_md5.csv"))
writeLines(trimws(capture.output(sessionInfo()), which = "right"),
           file.path(out, "sessionInfo.txt"))
print(metrics)
