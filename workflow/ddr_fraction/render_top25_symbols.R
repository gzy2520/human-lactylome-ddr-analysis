#!/usr/bin/env Rscript
# render_top25_symbols.R
# Re-renders top25 feature importance plots with readable Gene Symbols (Symbol + Ensembl ID).

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ragg)
  library(org.Hs.eg.db)
})

models_dir <- "outputs/20260927_ddr_fraction_models"

for (target in c("DDR", "DNA_repair")) {
  target_dir <- file.path(models_dir, target)
  imp_file <- file.path(target_dir, "tables/training_importance.csv")
  stopifnot(file.exists(imp_file))

  imp <- fread(imp_file)[1:25]

  # Map Ensembl to Symbol
  syms <- suppressMessages(
    mapIds(org.Hs.eg.db, keys = imp$Ensembl, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
  )
  syms[is.na(syms)] <- "Unannotated"
  imp[, Symbol := syms]
  imp[, GeneLabel := paste0(Symbol, " (", Ensembl, ")")]
  imp[, GeneLabel := factor(GeneLabel, levels = rev(GeneLabel))]

  # Save annotated table
  fwrite(imp, file.path(target_dir, "tables/training_importance_annotated.csv"))

  target_title <- if (target == "DDR") "Kla-DDR / DDR" else "Kla-DNA repair / DNA repair"

  g <- ggplot(imp, aes(x = Importance, y = GeneLabel)) +
    geom_col(fill = "#2166AC", alpha = 0.88, width = 0.72) +
    geom_text(aes(label = sprintf("%.0f", Importance)), hjust = -0.15, size = 3.2, colour = "grey25") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
    labs(
      title = paste0(target_title, " model: Top 25 training-split features"),
      subtitle = "Ranked across all 18,332 input genes; regression variance-reduction importance (Gini impurity)",
      x = "Impurity importance (Gini decrease)",
      y = "Gene: Symbol (Ensembl ID)",
      caption = "Mapped using org.Hs.eg.db. Labels preserve stable Ensembl IDs alongside official symbols.\nModel trained on n=1,524 split across 18,332 common RNA genes (500 trees, mtry=135)."
    ) +
    theme_minimal(base_size = 11, base_family = "Arial") +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      axis.text.y = element_text(size = 9, face = "bold", colour = "grey20"),
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 9.5, colour = "grey35"),
      plot.caption = element_text(size = 8, colour = "grey40")
    )

  out_png <- file.path(target_dir, "figures/top25.png")
  out_pdf <- file.path(target_dir, "figures/top25.pdf")

  ggsave(out_png, g, width = 10, height = 8.5, device = agg_png, dpi = 200)
  ggsave(out_pdf, g, width = 10, height = 8.5, device = cairo_pdf)
  cat("Successfully re-rendered top25 figures for", target, "\n")
}

# Update release_sha256.csv for DDR and DNA_repair
for (target in c("DDR", "DNA_repair")) {
  target_dir <- file.path(models_dir, target)
  fs <- sort(list.files(target_dir, recursive = TRUE, full.names = TRUE))
  fs <- fs[!grepl("release_sha256", fs)]
  dt <- data.table(
    Path = substring(fs, nchar(target_dir) + 2),
    SHA256 = vapply(fs, digest::digest, character(1), algo = "sha256", file = TRUE)
  )
  fwrite(dt, file.path(target_dir, "release_sha256.csv"))
}

cat("ALL TOP 25 FIGURES RE-RENDERED WITH SYMBOLS AND MANIFESTS UPDATED!\n")
