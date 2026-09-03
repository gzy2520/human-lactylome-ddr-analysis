# Shared significance calculations for the isolated boxplot review figures.
#
# The four-category and pathway-summary candidates use dataset-level points:
# one PXD/sample-group union per point.  These ANOVA results are descriptive
# review statistics and do not replace the frozen publication analysis or pool
# protein-level observations.

significance_label <- function(q_value) {
  if (is.na(q_value)) return("NA")
  if (q_value < 0.0001) return("****")
  if (q_value < 0.001) return("***")
  if (q_value < 0.01) return("**")
  if (q_value < 0.05) return("*")
  "ns"
}

safe_rank_sum_p <- function(x, y) {
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  if (!length(x) || !length(y)) return(NA_real_)
  stats::wilcox.test(x, y, alternative = "two.sided", exact = FALSE)$p.value
}

safe_signed_rank_p <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]
  if (!length(x)) return(NA_real_)
  difference <- x - y
  if (all(difference == 0)) return(1)
  stats::wilcox.test(x, y, paired = TRUE, alternative = "two.sided", exact = FALSE)$p.value
}

safe_aov_table <- function(formula, data) {
  tryCatch(
    summary(stats::aov(formula, data = data))[[1L]],
    error = function(e) NULL
  )
}

aov_row_index <- function(anova_table, term) {
  if (is.null(anova_table)) return(integer())
  matches <- which(trimws(rownames(anova_table)) == term)
  if (length(matches)) matches[[1L]] else integer()
}

compute_category_one_way_anova <- function(values, category_order, dataset_order = c("Whole proteome", "Lactylome (Kla)"), value_column = "DdrFractionPercentage") {
  stopifnot(value_column %in% names(values))
  result <- rbindlist(lapply(dataset_order, function(dataset) {
    panel <- values[as.character(Dataset) == dataset & as.character(Category) %in% category_order]
    panel <- panel[is.finite(get(value_column))]
    panel[, CategoryFactor := factor(as.character(Category), levels = category_order)]
    panel <- panel[!is.na(CategoryFactor)]
    table <- if (nrow(panel) && uniqueN(panel$CategoryFactor) >= 2L) {
      safe_aov_table(stats::as.formula(paste(value_column, "~ CategoryFactor")), panel)
    } else {
      NULL
    }
    category_row <- aov_row_index(table, "CategoryFactor")
    residual_row <- aov_row_index(table, "Residuals")
    p_value <- if (length(category_row)) {
      as.numeric(table[category_row, "Pr(>F)"])
    } else {
      NA_real_
    }
    f_value <- if (length(category_row)) {
      as.numeric(table[category_row, "F value"])
    } else {
      NA_real_
    }
    df_value <- if (length(category_row)) {
      as.numeric(table[category_row, "Df"])
    } else {
      NA_real_
    }
    residual_df <- if (length(residual_row)) {
      as.numeric(table[residual_row, "Df"])
    } else {
      NA_real_
    }
    data.table(
      Plot = "Figure 1 dataset-level four-category boxplot",
      Dataset = dataset,
      Category = NA_character_,
      Pathway = NA_character_,
      Term = "Category",
      Comparison = "Four-category omnibus test",
      Test = "one-way ANOVA",
      Pairing = "unpaired dataset-level PXD/sample-group points",
      N = nrow(panel),
      NCategory = uniqueN(panel$CategoryFactor),
      Df = df_value,
      ResidualDf = residual_df,
      FStatistic = f_value,
      PValue = p_value
    )
  }), fill = TRUE)
  result[, QValueBH := p.adjust(PValue, method = "BH")]
  result[, Significance := vapply(QValueBH, significance_label, character(1))]
  result[, AdjustmentFamily := "Figure 1 dataset-level boxplot: two modality-specific one-way ANOVA tests"]
  result[]
}

compute_pathway_two_way_anova <- function(values, category_order, pathway_order, dataset_order = c("Lactylome (Kla)", "Whole proteome")) {
  required <- c("Dataset", "Category", "Pathway", "PositiveFraction", "NegativeFraction")
  stopifnot(all(required %in% names(values)))
  direction_order <- c("Up/positive", "Down/negative")
  result <- rbindlist(lapply(dataset_order, function(dataset) {
    rbindlist(lapply(pathway_order, function(pathway) {
      panel <- values[as.character(Dataset) == dataset & as.character(Pathway) == pathway & as.character(Category) %in% category_order]
      point_ids <- if ("DatasetPointID" %in% names(panel)) panel$DatasetPointID else seq_len(nrow(panel))
      up <- data.table(
        DatasetPointID = point_ids,
        Category = as.character(panel$Category),
        Direction = direction_order[[1L]],
        Value = panel$PositiveFraction * 100
      )
      down <- data.table(
        DatasetPointID = point_ids,
        Category = as.character(panel$Category),
        Direction = direction_order[[2L]],
        Value = panel$NegativeFraction * 100
      )
      long <- rbindlist(list(up, down), use.names = TRUE, fill = TRUE)
      long[, CategoryFactor := factor(Category, levels = category_order)]
      long[, DirectionFactor := factor(Direction, levels = direction_order)]
      long <- long[is.finite(Value) & !is.na(CategoryFactor) & !is.na(DirectionFactor)]
      table <- if (nrow(long) && uniqueN(long$CategoryFactor) >= 2L && uniqueN(long$DirectionFactor) >= 2L) {
        safe_aov_table(Value ~ CategoryFactor * DirectionFactor, long)
      } else {
        NULL
      }
      terms <- c("CategoryFactor", "DirectionFactor", "CategoryFactor:DirectionFactor")
      rbindlist(lapply(terms, function(term) {
        term_row <- aov_row_index(table, term)
        residual_row <- aov_row_index(table, "Residuals")
        p_value <- if (length(term_row)) as.numeric(table[term_row, "Pr(>F)"]) else NA_real_
        f_value <- if (length(term_row)) as.numeric(table[term_row, "F value"]) else NA_real_
        df_value <- if (length(term_row)) as.numeric(table[term_row, "Df"]) else NA_real_
        residual_df <- if (length(residual_row)) as.numeric(table[residual_row, "Df"]) else NA_real_
        data.table(
          Plot = paste0(dataset, " / ", pathway),
          Dataset = dataset,
          Category = NA_character_,
          Pathway = pathway,
          Term = term,
          Comparison = "Category × Direction",
          Test = "two-way ANOVA",
          Pairing = "one Up and one Down value per dataset-level PXD/sample-group point",
          N = nrow(long),
          NPoint = uniqueN(long$DatasetPointID),
          NCategory = uniqueN(long$CategoryFactor),
          NDirection = uniqueN(long$DirectionFactor),
          Df = df_value,
          ResidualDf = residual_df,
          FStatistic = f_value,
          PValue = p_value
        )
      }), fill = TRUE)
    }), fill = TRUE)
  }), fill = TRUE)
  result[, QValueBH := p.adjust(PValue, method = "BH")]
  result[, Significance := vapply(QValueBH, significance_label, character(1))]
  result[, AdjustmentFamily := "DDR pathway-summary dataset-level boxplots: 14 plots x 3 two-way ANOVA terms"]
  result[]
}

compute_figure1_significance <- function(values, category_order) {
  result <- rbindlist(lapply(category_order, function(category) {
    panel <- values[Category == category]
    whole <- panel[Dataset == "Whole proteome", DdrFractionPercentage]
    kla <- panel[Dataset == "Lactylome (Kla)", DdrFractionPercentage]
    data.table(
      Plot = "Figure 1",
      Category = category,
      Pathway = NA_character_,
      Comparison = "Whole proteome vs Lactylome (Kla)",
      Test = "Wilcoxon rank-sum, two-sided",
      Pairing = "unpaired",
      N1 = length(whole),
      N2 = length(kla),
      PValue = safe_rank_sum_p(whole, kla)
    )
  }), fill = TRUE)
  result[, QValueBH := p.adjust(PValue, method = "BH")]
  result[, Significance := vapply(QValueBH, significance_label, character(1))]
  result[, AdjustmentFamily := "Figure 1: four category comparisons"]
  result[]
}

compute_pathway_significance <- function(values, category_order, pathway_order) {
  result <- rbindlist(lapply(category_order, function(category) {
    rbindlist(lapply(pathway_order, function(pathway) {
      panel <- values[Category == category & Pathway == pathway]
      data.table(
        Plot = "DDR pathway-summary",
        Category = category,
        Pathway = pathway,
        Comparison = "Positive score vs negative score",
        Test = "Wilcoxon signed-rank, two-sided",
        Pairing = "same source observation",
        N1 = sum(is.finite(panel$PositiveFraction) & is.finite(panel$NegativeFraction)),
        N2 = NA_integer_,
        PValue = safe_signed_rank_p(panel$PositiveFraction, panel$NegativeFraction)
      )
    }), fill = TRUE)
  }), fill = TRUE)
  result[, QValueBH := p.adjust(PValue, method = "BH")]
  result[, Significance := vapply(QValueBH, significance_label, character(1))]
  result[, AdjustmentFamily := "DDR pathway-summary: 28 category-pathway comparisons"]
  result[]
}

compute_ratio_significance <- function(values, category_order, denominator_order) {
  comparisons <- combn(category_order, 2L, simplify = FALSE)
  result <- rbindlist(lapply(denominator_order, function(denominator) {
    rbindlist(lapply(comparisons, function(pair) {
      first <- values[Denominator == denominator & Category == pair[[1L]], Ratio]
      second <- values[Denominator == denominator & Category == pair[[2L]], Ratio]
      data.table(
        Plot = paste0("Figure 1 MKI67/", denominator),
        Category = NA_character_,
        Pathway = NA_character_,
        Denominator = denominator,
        Comparison = paste(pair[[1L]], "vs", pair[[2L]]),
        Group1 = pair[[1L]],
        Group2 = pair[[2L]],
        Test = "Wilcoxon rank-sum, two-sided",
        Pairing = "unpaired source observations",
        N1 = sum(is.finite(first)),
        N2 = sum(is.finite(second)),
        PValue = safe_rank_sum_p(first, second)
      )
    }), fill = TRUE)
  }), fill = TRUE)
  result[, QValueBH := p.adjust(PValue, method = "BH"), by = Denominator]
  result[, Significance := vapply(QValueBH, significance_label, character(1))]
  result[, AdjustmentFamily := "MKI67 ratio: six category comparisons per denominator"]
  result[]
}

compute_ratio_global_significance <- function(values, category_order, denominator_order) {
  result <- rbindlist(lapply(denominator_order, function(denominator) {
    panel <- values[Denominator == denominator]
    panel <- panel[is.finite(Ratio) & Category %in% category_order]
    panel[, CategoryFactor := factor(Category, levels = category_order)]
    anova_table <- if (nrow(panel) && uniqueN(panel$CategoryFactor) >= 2L) {
      safe_aov_table(Ratio ~ CategoryFactor, panel)
    } else {
      NULL
    }
    category_row <- aov_row_index(anova_table, "CategoryFactor")
    p_value <- if (length(category_row)) {
      as.numeric(anova_table[category_row, "Pr(>F)"])
    } else {
      NA_real_
    }
    f_value <- if (length(category_row)) {
      as.numeric(anova_table[category_row, "F value"])
    } else {
      NA_real_
    }
    data.table(
      Plot = paste0("Figure 1 MKI67/", denominator),
      Category = NA_character_,
      Pathway = NA_character_,
      Denominator = denominator,
      Comparison = "Four-category omnibus test",
      Group1 = NA_character_,
      Group2 = NA_character_,
      Test = "one-way ANOVA",
      Pairing = "unpaired source observations",
      N1 = nrow(panel),
      N2 = NA_integer_,
      PValue = p_value,
      FStatistic = f_value
    )
  }), fill = TRUE)
  result[, QValueBH := p.adjust(PValue, method = "BH")]
  result[, Significance := vapply(QValueBH, significance_label, character(1))]
  result[, AdjustmentFamily := "MKI67 ratio: three omnibus denominator tests"]
  result[]
}
