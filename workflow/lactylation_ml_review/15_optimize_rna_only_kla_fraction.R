#!/usr/bin/env Rscript
# Bounded, source-held-out RNA-only optimization. No protein-derived predictor.
suppressPackageStartupMessages(library(data.table))
set.seed(25)
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
out <- args[[1L]]
if (dir.exists(out) && length(list.files(out, all.files = TRUE, no.. = TRUE)))
  stop("Refusing to overwrite a nonempty output directory")

base <- "outputs/20260923_rna_only_kla_fraction"
labels_file <- file.path(base, "material_labels_28.csv")
rna_file <- "outputs/20260920_lactate_metabolism/tables/material_gene_expression.csv"
module_file <- "config/lactate_metabolism/gene_module_membership.csv"
stopifnot(all(file.exists(c(labels_file, rna_file, module_file))))
d <- fread(labels_file)
e <- fread(rna_file)
mem <- fread(module_file)
setorder(d, ReferenceKey)
stopifnot(nrow(d) == 28L, uniqueN(d$ReferenceKey) == 28L,
          uniqueN(d$ConnGroup) == 15L, uniqueN(e$Ensembl) == 128L,
          nrow(e) == 28L * 128L,
          !anyDuplicated(e[, .(ReferenceKey, Ensembl)]),
          setequal(d$ReferenceKey, e$ReferenceKey))
genes <- sort(unique(e$Ensembl))
wide <- dcast(e[, .(ReferenceKey, Ensembl, Value = UnsmoothedLog2TPM)],
              ReferenceKey ~ Ensembl, value.var = "Value")
setorder(wide, ReferenceKey)
stopifnot(identical(wide$ReferenceKey, d$ReferenceKey),
          identical(names(wide)[-1L], genes), !anyNA(wide))
x <- as.matrix(wide[, ..genes])
storage.mode(x) <- "double"
y <- d$ObservedPercent
block <- d$ConnGroup

# Membership was assembled for the prior lactate-metabolism analysis. Symbols
# are annotations only; every feature is constructed from stable Ensembl IDs.
mem <- unique(mem[!is.na(Ensembl) & Ensembl %in% genes,
                  .(Origin, Module, Ensembl)])
catalog <- unique(mem[, .(Origin, Module)])
setorder(catalog, Origin, Module)
catalog[, FeatureID := sprintf("M%02d", .I)]
mem <- merge(mem, catalog, by = c("Origin", "Module"))
modules <- split(mem$Ensembl, mem$FeatureID)
stopifnot(nrow(catalog) == 16L, all(lengths(modules) > 0L))
go <- catalog[Origin == "GO annotation", FeatureID]
mechanism <- catalog[Origin == "curated mechanism", FeatureID]
core_names <- c("glycolysis", "pyruvate fermentation to lactate",
                "lactate transport", "tricarboxylic acid cycle")
core <- catalog[Origin == "GO annotation" & Module %in% core_names, FeatureID]
stopifnot(length(go) == 6L, length(mechanism) == 10L,
          length(core) == 4L)
catalog[, MeasuredGenes := lengths(modules[FeatureID])]

design <- function(train, test, set_name) {
  center <- colMeans(x[train, , drop = FALSE])
  spread <- apply(x[train, , drop = FALSE], 2L, sd)
  spread[!is.finite(spread) | spread == 0] <- 1
  ztr <- sweep(sweep(x[train, , drop = FALSE], 2L, center, "-"),
               2L, spread, "/")
  zte <- sweep(sweep(x[test, , drop = FALSE], 2L, center, "-"),
               2L, spread, "/")
  if (set_name == "RNA128") {
    a <- ztr
    b <- zte
  } else {
    ids <- switch(set_name, Core4 = core, GO6 = go,
                  Mechanism10 = mechanism, All16 = catalog$FeatureID)
    stopifnot(length(ids) > 0L)
    score <- function(z) vapply(ids, function(id)
      rowMeans(z[, modules[[id]], drop = FALSE]), numeric(nrow(z)))
    a <- matrix(score(ztr), nrow = length(train), ncol = length(ids))
    b <- matrix(score(zte), nrow = length(test), ncol = length(ids))
    mc <- colMeans(a)
    ms <- apply(a, 2L, sd)
    ms[!is.finite(ms) | ms == 0] <- 1
    a <- sweep(sweep(a, 2L, mc, "-"), 2L, ms, "/")
    b <- sweep(sweep(b, 2L, mc, "-"), 2L, ms, "/")
  }
  list(train = a / sqrt(ncol(a)), test = b / sqrt(ncol(b)))
}

weighted_median <- function(v, groups) {
  w <- 1 / as.vector(table(groups)[as.character(groups)])
  ix <- order(v)
  v[ix][which(cumsum(w[ix]) >= sum(w) / 2)[1L]]
}
predict_model <- function(train, test, model, tuning) {
  if (model == "MedianBaseline")
    return(rep(median(y[train]), length(test)))
  if (model == "BalancedMedianBaseline")
    return(rep(weighted_median(y[train], block[train]), length(test)))
  set_name <- sub("_(ridge|huber|knn)$", "", model)
  kind <- sub("^.*_", "", model)
  xz <- design(train, test, set_name)
  a <- xz$train
  b <- xz$test
  if (kind == "ridge") {
    ym <- mean(y[train])
    alpha <- solve(tcrossprod(a) + diag(tuning, nrow(a)), y[train] - ym)
    pred <- as.vector(ym + b %*% crossprod(a, alpha))
  } else if (kind == "huber") {
    delta <- 3  # smooth absolute-error loss, in percentage points
    n <- length(train)
    fn <- function(theta) {
      residual <- y[train] - theta[[1L]] - as.vector(a %*% theta[-1L])
      mean(sqrt(residual^2 + delta^2) - delta) +
        tuning * sum(theta[-1L]^2) / 2
    }
    gr <- function(theta) {
      residual <- y[train] - theta[[1L]] - as.vector(a %*% theta[-1L])
      w <- residual / sqrt(residual^2 + delta^2)
      c(-mean(w), -as.vector(crossprod(a, w)) / n +
          tuning * theta[-1L])
    }
    fit <- optim(c(median(y[train]), rep(0, ncol(a))), fn, gr,
                 method = "BFGS", control = list(maxit = 150L, reltol = 1e-8))
    if (fit$convergence != 0L) stop("Robust fit did not converge")
    pred <- as.vector(fit$par[[1L]] + b %*% fit$par[-1L])
  } else if (kind == "knn") {
    pred <- vapply(seq_len(nrow(b)), function(i) {
      distance <- rowSums((sweep(a, 2L, b[i, ], "-"))^2)
      nearest <- order(distance, train)[seq_len(as.integer(tuning))]
      median(y[train][nearest])
    }, numeric(1))
  } else stop("Unknown model")
  pmin(100, pmax(0, pred))
}

grid <- rbindlist(list(
  data.table(Model = c("BalancedMedianBaseline", "MedianBaseline"), Tuning = 0),
  CJ(Model = paste0(c("Core4", "GO6", "Mechanism10", "All16"), "_ridge"),
     Tuning = c(100, 30, 10, 3, 1, 0.3, 0.1)),
  CJ(Model = paste0(c("Core4", "GO6"), "_huber"),
     Tuning = c(1, 0.3, 0.1, 0.03, 0.01)),
  CJ(Model = paste0(c("GO6", "RNA128"), "_knn"),
     Tuning = c(3, 5, 7))
))
model_names <- setdiff(unique(grid$Model),
                       c("BalancedMedianBaseline", "MedianBaseline"))
predictions <- list()
selection <- list()
inner_losses <- list()
for (outer in sort(unique(block))) {
  tr <- which(block != outer)
  te <- which(block == outer)
  inner <- sort(unique(block[tr]))
  loss <- grid[, .(InnerMeanGroupMAE = mean(vapply(inner, function(held) {
    tr2 <- tr[block[tr] != held]
    te2 <- tr[block[tr] == held]
    mean(abs(predict_model(tr2, te2, Model, Tuning) - y[te2]))
  }, numeric(1)))), by = .(Model, Tuning)]
  loss[, OuterConnGroup := outer]
  inner_losses[[length(inner_losses) + 1L]] <- loss
  setorder(loss, InnerMeanGroupMAE, Model, -Tuning)
  best <- loss[1L]
  selection[[length(selection) + 1L]] <- best
  row <- d[te, .(ReferenceKey, ConnGroup, ObservedPercent)]
  row[, MedianBaseline := median(y[tr])]
  row[, BalancedMedianBaseline := weighted_median(y[tr], block[tr])]
  for (nm in model_names) {
    opt <- loss[Model == nm][order(InnerMeanGroupMAE, -Tuning)][1L]
    row[, (nm) := predict_model(tr, te, nm, opt$Tuning)]
  }
  row[, SelectedWithFallback := predict_model(tr, te, best$Model, best$Tuning)]
  predictions[[length(predictions) + 1L]] <- row
}
pred <- rbindlist(predictions)
setorder(pred, ReferenceKey)
stopifnot(nrow(pred) == 28L, !anyNA(pred),
          identical(pred$ReferenceKey, d$ReferenceKey))
cols <- c("MedianBaseline", "BalancedMedianBaseline", model_names,
          "SelectedWithFallback")
perf <- rbindlist(lapply(cols, function(nm) {
  err <- pred[[nm]] - pred$ObservedPercent
  mae <- pred[, .(MAE = mean(abs(get(nm) - ObservedPercent))),
              by = ConnGroup]
  data.table(Model = nm, MaterialMAE_pp = mean(abs(err)),
             MaterialRMSE_pp = sqrt(mean(err^2)),
             ConnectedGroupMeanMAE_pp = mean(mae$MAE))
}))
setorder(perf, MaterialMAE_pp)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
fwrite(catalog, file.path(out, "module_catalog.csv"))
fwrite(pred, file.path(out, "material_oof_predictions.csv"))
fwrite(perf, file.path(out, "performance_summary.csv"))
fwrite(rbindlist(selection), file.path(out, "selected_parameters.csv"))
fwrite(rbindlist(inner_losses), file.path(out, "inner_tuning.csv"))
fwrite(data.table(File = c(labels_file, rna_file, module_file),
                  MD5 = unname(tools::md5sum(c(labels_file, rna_file,
                                               module_file)))),
       file.path(out, "input_md5.csv"))
writeLines(trimws(capture.output(sessionInfo()), which = "right"),
           file.path(out, "sessionInfo.txt"))
print(perf)
