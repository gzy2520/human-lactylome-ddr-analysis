# Shared ANOVA calculations for the isolated, source-resolved review figures.
#
# These functions intentionally operate on deposited sample observations. All
# protein membership upstream is defined with stable UniProt BaseAccessions.

significance_label <- function(q_value) {
  if (!is.finite(q_value)) return("NA")
  if (q_value < 0.0001) return("****")
  if (q_value < 0.001) return("***")
  if (q_value < 0.01) return("**")
  if (q_value < 0.05) return("*")
  "ns"
}

safe_aov_table <- function(formula, data) {
  tryCatch(
    summary(stats::aov(formula, data = data))[[1L]],
    error = function(error) NULL
  )
}

aov_value <- function(anova_table, term, column) {
  if (is.null(anova_table) || !term %in% trimws(rownames(anova_table)) ||
      !column %in% colnames(anova_table)) return(NA_real_)
  as.numeric(anova_table[match(term, trimws(rownames(anova_table))), column])
}

compute_figure1_sample_one_way_anova <- function(
    values,
    category_order,
    dataset_order = c("Whole proteome", "Lactylome (Kla)")) {
  required <- c("Category", "Dataset", "DdrFractionPercentage")
  stopifnot(all(required %in% names(values)))
  result <- data.table::rbindlist(lapply(category_order, function(category) {
    panel <- data.table::copy(values[as.character(Category) == category &
                                      as.character(Dataset) %in% dataset_order])
    panel <- panel[is.finite(DdrFractionPercentage)]
    panel[, DatasetFactor := factor(as.character(Dataset), levels = dataset_order)]
    table <- if (data.table::uniqueN(panel$DatasetFactor) == length(dataset_order)) {
      safe_aov_table(DdrFractionPercentage ~ DatasetFactor, panel)
    } else {
      NULL
    }
    data.table::data.table(
      Plot = "Figure 1 source-sample category panel",
      Category = category,
      Dataset = NA_character_,
      Pathway = NA_character_,
      Term = "DatasetFactor",
      Comparison = "Whole proteome versus Lactylome (Kla)",
      Test = "one-way ANOVA",
      AnalysisScale = "DDR fraction percentage",
      N = nrow(panel),
      NWholeProteome = sum(panel$DatasetFactor == "Whole proteome"),
      NKla = sum(panel$DatasetFactor == "Lactylome (Kla)"),
      Df = aov_value(table, "DatasetFactor", "Df"),
      ResidualDf = aov_value(table, "Residuals", "Df"),
      FStatistic = aov_value(table, "DatasetFactor", "F value"),
      PValue = aov_value(table, "DatasetFactor", "Pr(>F)")
    )
  }), fill = TRUE)
  result[, QValueBH := stats::p.adjust(PValue, method = "BH")]
  result[, Significance := vapply(QValueBH, significance_label, character(1))]
  result[, AdjustmentFamily := "Four Figure 1 category-specific one-way ANOVA tests"]
  result[]
}

compute_figure1_original_dataset_one_way_anova <- function(
    values,
    dataset_order = c("Whole proteome", "Lactylome (Kla)")) {
  required <- c("RowOrder", "PXD", "SampleGroup", "Category", "Dataset", "DdrFractionPercentage")
  stopifnot(all(required %in% names(values)))
  group_columns <- c("RowOrder", "PXD", "SampleGroup", "Category")
  groups <- unique(data.table::copy(values[, ..group_columns]))
  data.table::setorder(groups, RowOrder, PXD, SampleGroup)
  result <- data.table::rbindlist(lapply(seq_len(nrow(groups)), function(index) {
    group <- groups[index]
    panel <- data.table::copy(values[
      RowOrder == group$RowOrder & PXD == group$PXD &
        SampleGroup == group$SampleGroup & as.character(Dataset) %in% dataset_order
    ])
    panel <- panel[is.finite(DdrFractionPercentage)]
    panel[, DatasetFactor := factor(as.character(Dataset), levels = dataset_order)]
    n_whole <- sum(panel$DatasetFactor == dataset_order[[1L]], na.rm = TRUE)
    n_kla <- sum(panel$DatasetFactor == dataset_order[[2L]], na.rm = TRUE)
    table <- if (n_whole > 0L && n_kla > 0L && nrow(panel) > 2L) {
      safe_aov_table(DdrFractionPercentage ~ DatasetFactor, panel)
    } else {
      NULL
    }
    data.table::data.table(
      Plot = "Figure 1 PXD-axis source-sample row",
      RowOrder = group$RowOrder,
      PXD = group$PXD,
      SampleGroup = group$SampleGroup,
      Category = group$Category,
      Dataset = NA_character_,
      Pathway = NA_character_,
      Term = "DatasetFactor",
      Comparison = "Whole proteome versus Lactylome (Kla) within one PXD/sample-group row",
      Test = "one-way ANOVA",
      AnalysisScale = "DDR fraction percentage",
      N = nrow(panel),
      NWholeProteome = n_whole,
      NKla = n_kla,
      Df = aov_value(table, "DatasetFactor", "Df"),
      ResidualDf = aov_value(table, "Residuals", "Df"),
      FStatistic = aov_value(table, "DatasetFactor", "F value"),
      PValue = aov_value(table, "DatasetFactor", "Pr(>F)")
    )
  }), fill = TRUE)
  result[, QValueBH := stats::p.adjust(PValue, method = "BH")]
  result[, Significance := vapply(QValueBH, significance_label, character(1))]
  result[, AdjustmentFamily := paste(
    nrow(groups),
    "Figure 1 PXD-axis row-specific one-way ANOVA tests"
  )]
  result[]
}

compute_pathway_sample_two_way_anova <- function(values, category_order, pathway_order,
                                                 direction_order = c("Pro", "Inh")) {
  required <- c("Category", "Pathway", "PositiveFraction", "NegativeFraction")
  stopifnot(all(required %in% names(values)))
  result <- data.table::rbindlist(lapply(pathway_order, function(pathway) {
    panel <- data.table::copy(values[as.character(Pathway) == pathway &
                                      as.character(Category) %in% category_order])
    id <- if (all(c("PXD", "SampleGroup", "SampleID") %in% names(panel))) {
      paste(panel$PXD, panel$SampleGroup, panel$SampleID, sep = "|")
    } else if ("SampleID" %in% names(panel)) {
      as.character(panel$SampleID)
    } else {
      as.character(seq_len(nrow(panel)))
    }
    long <- data.table::rbindlist(list(
      data.table::data.table(SourceSampleID = id, Category = panel$Category,
        Direction = direction_order[[1L]], ValuePercent = panel$PositiveFraction * 100),
      data.table::data.table(SourceSampleID = id, Category = panel$Category,
        Direction = direction_order[[2L]], ValuePercent = panel$NegativeFraction * 100)
    ))
    long[, CategoryFactor := factor(as.character(Category), levels = category_order)]
    long[, DirectionFactor := factor(Direction, levels = direction_order)]
    long <- long[is.finite(ValuePercent) & !is.na(CategoryFactor) & !is.na(DirectionFactor)]
    table <- if (data.table::uniqueN(long$CategoryFactor) == length(category_order) &&
                 data.table::uniqueN(long$DirectionFactor) == length(direction_order)) {
      safe_aov_table(ValuePercent ~ CategoryFactor * DirectionFactor, long)
    } else {
      NULL
    }
    terms <- c("CategoryFactor", "DirectionFactor", "CategoryFactor:DirectionFactor")
    data.table::rbindlist(lapply(terms, function(term) {
      data.table::data.table(
        Plot = paste0("DDR pathway source-sample panel: ", pathway),
        Category = NA_character_,
        Dataset = "Lactylome (Kla)",
        Pathway = pathway,
        Term = term,
        Comparison = "Four categories by Pro versus Inh pathway state",
        Test = "two-way ANOVA",
        AnalysisScale = "Kla-DDR pathway fraction percentage",
        N = nrow(long),
        NPoint = data.table::uniqueN(long$SourceSampleID),
        Df = aov_value(table, term, "Df"),
        ResidualDf = aov_value(table, "Residuals", "Df"),
        FStatistic = aov_value(table, term, "F value"),
        PValue = aov_value(table, term, "Pr(>F)")
      )
    }))
  }), fill = TRUE)
  result[, QValueBH := stats::p.adjust(PValue, method = "BH")]
  result[, Significance := vapply(QValueBH, significance_label, character(1))]
  result[, AdjustmentFamily := "Seven pathways times three two-way ANOVA terms"]
  result[]
}

compute_ratio_global_significance <- function(values, denominator_order, category_order) {
  required <- c("Denominator", "Category", "Ratio")
  stopifnot(all(required %in% names(values)))
  result <- data.table::rbindlist(lapply(denominator_order, function(denominator) {
    panel <- data.table::copy(values[as.character(Denominator) == denominator])
    panel <- panel[is.finite(Ratio) & Ratio > 0]
    panel[, CategoryFactor := factor(as.character(Category), levels = category_order)]
    table <- if (data.table::uniqueN(panel$CategoryFactor) == length(category_order)) {
      safe_aov_table(log10(Ratio) ~ CategoryFactor, panel)
    } else {
      NULL
    }
    data.table::data.table(
      Plot = paste0("Figure 1 MKI67/", denominator),
      Denominator = denominator,
      Term = "CategoryFactor",
      Comparison = "Four-category omnibus test",
      Test = "one-way ANOVA",
      AnalysisScale = "log10(MKI67 / denominator)",
      N = nrow(panel),
      NCategory = data.table::uniqueN(panel$CategoryFactor),
      Df = aov_value(table, "CategoryFactor", "Df"),
      ResidualDf = aov_value(table, "Residuals", "Df"),
      FStatistic = aov_value(table, "CategoryFactor", "F value"),
      PValue = aov_value(table, "CategoryFactor", "Pr(>F)")
    )
  }), fill = TRUE)
  result[, QValueBH := stats::p.adjust(PValue, method = "BH")]
  result[, Significance := vapply(QValueBH, significance_label, character(1))]
  result[, AdjustmentFamily := "MKI67 ratio: three omnibus one-way ANOVA tests"]
  result[]
}
