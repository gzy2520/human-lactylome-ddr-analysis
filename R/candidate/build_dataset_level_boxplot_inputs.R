# Helpers for the dataset-level boxplot review figures.
#
# The frozen publication tables remain the source of truth.  These functions
# make a review-only, one-row-per-PXD/sample-group representation so that a
# point in a four-category plot is not inadvertently treated as a biological
# replicate when a source contains several sample observations.

dataset_boxplot_pathway_order <- c("BER", "NER", "MMR", "FA", "HR", "AEJ", "NHEJ")

dataset_boxplot_stop_if <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

dataset_boxplot_base_accession <- function(values) {
  values <- trimws(as.character(values))
  values[is.na(values)] <- ""
  values <- sub("^.*\\|([^|]+)\\|.*$", "\\1", values)
  values <- sub("^([^|;]+)\\|.*$", "\\1", values)
  values <- sub("^NX_", "", values)
  sub("-[0-9]+$", "", values)
}

dataset_boxplot_is_true <- function(values) {
  values <- tolower(trimws(as.character(values)))
  values[is.na(values)] <- ""
  values %in% c("+", "1", "true", "yes", "y")
}

dataset_boxplot_split_accessions <- function(values) {
  tokens <- unlist(strsplit(as.character(values), "[;,]"), use.names = FALSE)
  tokens <- dataset_boxplot_base_accession(tokens)
  sort(unique(tokens[nzchar(tokens)]))
}

dataset_boxplot_read_pathway_scores <- function(path) {
  dataset_boxplot_stop_if(file.exists(path), paste0("Missing pathway-score workbook: ", path))
  dataset_boxplot_stop_if(requireNamespace("readxl", quietly = TRUE),
    "The readxl package is required to read the S4 pathway workbook.")
  sheets <- readxl::excel_sheets(path)
  score_tables <- lapply(sheets, function(sheet) {
    data <- data.table::as.data.table(readxl::read_excel(path, sheet = sheet))
    required <- c("BaseAccession", dataset_boxplot_pathway_order)
    dataset_boxplot_stop_if(all(required %in% names(data)),
      paste0("Missing pathway score columns in S4 sheet ", sheet, "."))
    data[, c("BaseAccession", dataset_boxplot_pathway_order), with = FALSE]
  })
  scores <- data.table::rbindlist(score_tables, fill = TRUE)
  scores[, BaseAccession := dataset_boxplot_base_accession(BaseAccession)]
  scores <- scores[nzchar(BaseAccession)]
  for (pathway in dataset_boxplot_pathway_order) {
    conflicts <- scores[, .(DistinctScoreCount = uniqueN(get(pathway))), by = BaseAccession][DistinctScoreCount > 1L]
    dataset_boxplot_stop_if(!nrow(conflicts),
      paste0("Frozen pathway scores disagree for ", pathway, "."))
  }
  score_matrix <- as.matrix(scores[, ..dataset_boxplot_pathway_order])
  storage.mode(score_matrix) <- "numeric"
  dataset_boxplot_stop_if(all(score_matrix %in% c(-1, 0, 1)),
    "Frozen pathway scores must be -1, 0 or +1.")
  scores <- scores[, lapply(.SD, function(values) values[[1L]]),
    by = BaseAccession, .SDcols = dataset_boxplot_pathway_order]
  scores[]
}

build_dataset_level_figure1_values <- function(groups) {
  required <- c(
    "RowOrder", "PXD", "SampleGroup", "Category", "KlaLabelEn", "BiologicalMaterial",
    "KlaEvidenceFile", "KlaProteinCount", "KlaDdrProteinCount", "KlaDdrFraction",
    "ReferencePXD", "ReferenceEvidenceFile", "ReferenceProteinCount",
    "ReferenceDdrProteinCount", "ReferenceDdrFraction"
  )
  dataset_boxplot_stop_if(all(required %in% names(groups)),
    "Group summary is missing columns needed for dataset-level Figure 1 values.")
  dataset_boxplot_stop_if(!anyDuplicated(groups[, .(PXD, SampleGroup)]),
    "Group summary contains duplicate PXD/sample-group keys.")

  meta <- groups[, .(
    RowOrder, PXD, SampleGroup, Category, DisplayGroup = KlaLabelEn,
    BiologicalMaterial, ReferencePXD,
    KlaEvidenceFile, ReferenceEvidenceFile,
    KlaProteinCount, KlaDdrProteinCount, KlaDdrFraction,
    ReferenceProteinCount, ReferenceDdrProteinCount, ReferenceDdrFraction
  )]

  kla <- meta[, .(
    RowOrder, PXD, SampleGroup, Category, DisplayGroup, BiologicalMaterial,
    Dataset = "Lactylome (Kla)",
    DatasetPXD = PXD,
    SampleID = paste0(PXD, "_", "Kla_dataset_union"),
    ConditionLabel = DisplayGroup,
    SampleClass = BiologicalMaterial,
    ObservationType = "dataset_union",
    SourceMode = "validated_publication_group_union",
    SourceFile = KlaEvidenceFile,
    ReferencePXD,
    ProteinCount = KlaProteinCount,
    DdrProteinCount = KlaDdrProteinCount,
    DdrFraction = KlaDdrFraction,
    DdrFractionPercentage = KlaDdrFraction * 100,
    FrozenKlaProteinCount = KlaProteinCount,
    FrozenKlaDdrProteinCount = KlaDdrProteinCount,
    FrozenKlaDdrFraction = KlaDdrFraction * 100,
    FrozenReferenceProteinCount = ReferenceProteinCount,
    FrozenReferenceDdrProteinCount = ReferenceDdrProteinCount,
    FrozenReferenceDdrFraction = ReferenceDdrFraction * 100
  )]
  reference <- meta[, .(
    RowOrder, PXD, SampleGroup, Category, DisplayGroup, BiologicalMaterial,
    Dataset = "Whole proteome",
    DatasetPXD = ReferencePXD,
    SampleID = paste0(ReferencePXD, "_", "reference_dataset_union"),
    ConditionLabel = paste0(ReferencePXD, " tumor/reference union"),
    SampleClass = BiologicalMaterial,
    ObservationType = "dataset_union",
    SourceMode = "validated_reference_group_union",
    SourceFile = ReferenceEvidenceFile,
    ReferencePXD,
    ProteinCount = ReferenceProteinCount,
    DdrProteinCount = ReferenceDdrProteinCount,
    DdrFraction = ReferenceDdrFraction,
    DdrFractionPercentage = ReferenceDdrFraction * 100,
    FrozenKlaProteinCount = KlaProteinCount,
    FrozenKlaDdrProteinCount = KlaDdrProteinCount,
    FrozenKlaDdrFraction = KlaDdrFraction * 100,
    FrozenReferenceProteinCount = ReferenceProteinCount,
    FrozenReferenceDdrProteinCount = ReferenceDdrProteinCount,
    FrozenReferenceDdrFraction = ReferenceDdrFraction * 100
  )]
  output <- data.table::rbindlist(list(kla, reference), use.names = TRUE, fill = TRUE)
  output[, DatasetPointID := paste(PXD, SampleGroup, Dataset, sep = "__")]
  output[, DatasetPointLabel := paste(PXD, SampleGroup, sep = " / ")]
  data.table::setorder(output, RowOrder, Dataset, PXD, SampleGroup)
  dataset_boxplot_stop_if(!anyDuplicated(output$DatasetPointID),
    "Dataset-level Figure 1 values contain duplicate point keys.")
  dataset_boxplot_stop_if(all(is.finite(output$DdrFractionPercentage)),
    "Dataset-level Figure 1 values contain non-finite fractions.")
  dataset_boxplot_stop_if(all(output$DdrFractionPercentage >= 0 & output$DdrFractionPercentage <= 100),
    "Dataset-level Figure 1 percentages must be between 0 and 100.")
  output[]
}

build_dataset_level_pathway_summary <- function(groups, kla_membership, reference_membership, pathway_scores) {
  required_groups <- c(
    "RowOrder", "PXD", "SampleGroup", "Category", "KlaLabelEn", "BiologicalMaterial",
    "KlaEvidenceFile", "ReferencePXD", "ReferenceEvidenceFile",
    "KlaDdrProteinCount", "ReferenceDdrProteinCount"
  )
  dataset_boxplot_stop_if(all(required_groups %in% names(groups)),
    "Group summary is missing pathway-summary metadata.")
  required_kla <- c("PXD", "SampleGroup", "BaseAccession", "IsDdr")
  required_reference <- c("PXD", "SampleGroup", "ReferencePXD", "MappedBaseAccessions", "IsDdr")
  dataset_boxplot_stop_if(all(required_kla %in% names(kla_membership)),
    "Kla membership is missing pathway-summary columns.")
  dataset_boxplot_stop_if(all(required_reference %in% names(reference_membership)),
    "Reference membership is missing pathway-summary columns.")
  dataset_boxplot_stop_if(all(c("BaseAccession", dataset_boxplot_pathway_order) %in% names(pathway_scores)),
    "Pathway scores are missing the seven pathway states.")

  kla_ids <- unique(kla_membership[
    dataset_boxplot_is_true(IsDdr),
    .(PXD, SampleGroup, BaseAccession = dataset_boxplot_base_accession(BaseAccession))
  ])
  kla_ids <- kla_ids[nzchar(BaseAccession)]

  reference_ids <- reference_membership[dataset_boxplot_is_true(IsDdr), .(
    PXD, SampleGroup, ReferencePXD,
    BaseAccession = unlist(lapply(MappedBaseAccessions, dataset_boxplot_split_accessions), use.names = FALSE)
  ), by = .(PXD, SampleGroup, ReferencePXD, IsDdr)]
  reference_ids <- unique(reference_ids[nzchar(BaseAccession), .(PXD, SampleGroup, ReferencePXD, BaseAccession)])

  make_rows <- function(meta_row, dataset, ids, denominator, source_file, score_basis) {
    ids <- sort(unique(ids[nzchar(ids)]))
    mapped <- pathway_scores[BaseAccession %in% ids]
    mapped_ids <- unique(mapped$BaseAccession)
    dataset_boxplot_stop_if(denominator > 0L,
      paste0("No DDR denominator for ", meta_row$PXD, " / ", meta_row$SampleGroup, " / ", dataset, "."))
    dataset_boxplot_stop_if(dataset != "Lactylome (Kla)" || length(setdiff(ids, mapped_ids)) == 0L,
      paste0("A Kla-DDR group accession is absent from the S4 pathway score table for ", meta_row$PXD, "."))
    data.table::rbindlist(lapply(dataset_boxplot_pathway_order, function(pathway) {
      states <- if (nrow(mapped)) suppressWarnings(as.numeric(mapped[[pathway]])) else numeric()
      positive_count <- sum(states == 1, na.rm = TRUE)
      negative_count <- sum(states == -1, na.rm = TRUE)
      any_count <- sum(states != 0, na.rm = TRUE)
      data.table::data.table(
        RowOrder = meta_row$RowOrder,
        PXD = meta_row$PXD,
        SampleGroup = meta_row$SampleGroup,
        Category = meta_row$Category,
        DisplayGroup = meta_row$DisplayGroup,
        BiologicalMaterial = meta_row$BiologicalMaterial,
        Dataset = dataset,
        DatasetPXD = if (dataset == "Lactylome (Kla)") meta_row$PXD else meta_row$ReferencePXD,
        DatasetPointID = paste(meta_row$PXD, meta_row$SampleGroup, dataset, sep = "__"),
        DatasetPointLabel = paste(meta_row$PXD, meta_row$SampleGroup, sep = " / "),
        SampleID = if (dataset == "Lactylome (Kla)") paste0(meta_row$PXD, "_Kla_dataset_union") else paste0(meta_row$ReferencePXD, "_reference_dataset_union"),
        ConditionLabel = meta_row$DisplayGroup,
        SampleClass = meta_row$BiologicalMaterial,
        ObservationType = "dataset_union",
        SourceMode = if (dataset == "Lactylome (Kla)") "validated_publication_group_union" else "validated_reference_group_union",
        SourceFile = source_file,
        ReferencePXD = meta_row$ReferencePXD,
        Pathway = pathway,
        PositiveProteinCount = positive_count,
        NegativeProteinCount = negative_count,
        AnyPathwayProteinCount = any_count,
        PathwayScoreMappedProteinCount = length(mapped_ids),
        DdrProteinCount = denominator,
        KlaDdrProteinCount = meta_row$KlaDdrProteinCount,
        ReferenceDdrProteinCount = meta_row$ReferenceDdrProteinCount,
        PathwayScoreCoverage = length(mapped_ids) / denominator,
        PathwayScoreBasis = score_basis,
        PositiveFraction = positive_count / denominator,
        NegativeFraction = negative_count / denominator,
        SignedFraction = (positive_count - negative_count) / denominator
      )
    }))
  }

  meta <- groups[, .(
    RowOrder, PXD, SampleGroup, Category, DisplayGroup = KlaLabelEn, BiologicalMaterial,
    ReferencePXD, KlaEvidenceFile, ReferenceEvidenceFile,
    KlaDdrProteinCount, ReferenceDdrProteinCount
  )]
  rows <- lapply(seq_len(nrow(meta)), function(index) {
    item <- meta[index]
    kla <- kla_ids[PXD == item$PXD & SampleGroup == item$SampleGroup, BaseAccession]
    reference <- reference_ids[PXD == item$PXD & SampleGroup == item$SampleGroup, BaseAccession]
    data.table::rbindlist(list(
      make_rows(item, "Lactylome (Kla)", kla, item$KlaDdrProteinCount, item$KlaEvidenceFile,
        "Supplementary Table S4 seven-pathway signed states; denominator is all Kla-DDR proteins"),
      make_rows(item, "Whole proteome", reference, item$ReferenceDdrProteinCount, item$ReferenceEvidenceFile,
        "Supplementary Table S4 seven-pathway signed states mapped by BaseAccession; denominator is all reference-DDR proteins")
    ), fill = TRUE)
  })
  output <- data.table::rbindlist(rows, fill = TRUE)
  output[, PathwayOrder := match(Pathway, dataset_boxplot_pathway_order)]
  data.table::setorder(output, RowOrder, Dataset, DatasetPointID, PathwayOrder)
  output[, PathwayOrder := NULL]
  dataset_boxplot_stop_if(nrow(output) == nrow(groups) * 2L * length(dataset_boxplot_pathway_order),
    "Dataset-level pathway summary is incomplete.")
  dataset_boxplot_stop_if(!anyDuplicated(output[, .(DatasetPointID, Pathway)]),
    "Dataset-level pathway summary contains duplicate point/pathway rows.")
  dataset_boxplot_stop_if(all(is.finite(output$PositiveFraction) & is.finite(output$NegativeFraction)),
    "Dataset-level pathway summary contains non-finite fractions.")
  dataset_boxplot_stop_if(all(output$PositiveFraction >= 0 & output$PositiveFraction <= 1),
    "A positive dataset-level pathway fraction is outside 0-1.")
  dataset_boxplot_stop_if(all(output$NegativeFraction >= 0 & output$NegativeFraction <= 1),
    "A negative dataset-level pathway fraction is outside 0-1.")
  output[]
}
