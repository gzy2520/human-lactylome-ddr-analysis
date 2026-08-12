#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(eulerr)
  library(readr)
  library(stringr)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) normalizePath(args[[1]]) else normalizePath(".")
figure_root <- file.path(
  project_root,
  "results", "figures",
  "four_class_venn"
)
table_root <- file.path(
  project_root,
  "results", "tables",
  "four_class_venn"
)
dir.create(figure_root, recursive = TRUE, showWarnings = FALSE)
dir.create(table_root, recursive = TRUE, showWarnings = FALSE)

stats_path <- file.path(
  project_root,
  "results", "tables",
  "cell_type_kla_vs_reference_ddr_statistics_accession_only.csv"
)
kla_path <- file.path(
  project_root,
  "work", "intermediate",
  "expanded_ddr_by_accession",
  "kla_proteins_by_sample_group.csv"
)
reference_path <- file.path(
  project_root,
  "work", "intermediate",
  "expanded_ddr_by_accession",
  "reference_proteins_by_sample_group.csv"
)
all_uniprot_path <- file.path(
  project_root,
  "config",
  "uniprot_human_all_2026-08-06.tsv"
)
go_path <- file.path(
  project_root,
  "data",
  "annotations",
  "GO-repair+damage(human).tsv"
)

required <- c(
  stats_path,
  kla_path,
  reference_path,
  all_uniprot_path,
  go_path
)
missing <- required[!file.exists(required)]
if (length(missing)) {
  stop("Missing required input files: ", paste(missing, collapse = ", "))
}

category_order <- c(
  "normal_tissue",
  "cancer_tissue",
  "normal_cells",
  "cancer_cells"
)
category_labels_zh <- c(
  normal_tissue = "正常/非肿瘤组织",
  cancer_tissue = "癌症组织",
  normal_cells = "正常/非肿瘤细胞",
  cancer_cells = "癌症细胞"
)
category_labels_en <- c(
  normal_tissue = "Normal/non-tumor tissues",
  cancer_tissue = "Cancer tissues",
  normal_cells = "Normal/non-tumor cells",
  cancer_cells = "Cancer cells"
)

base_accession <- function(values) {
  values <- trimws(as.character(values))
  values <- sub("^.*\\|([^|]+)\\|.*$", "\\1", values)
  values <- sub("^([^|;]+)\\|.*$", "\\1", values)
  values <- sub("-[0-9]+$", "", values)
  values
}

clean_id_vector <- function(values) {
  values <- unlist(strsplit(as.character(values), "[;,]"))
  values <- base_accession(values)
  sort(unique(values[!is.na(values) & nzchar(values)]))
}

stats <- read_csv(stats_path, show_col_types = FALSE) |>
  mutate(
    SampleGroupKey = paste(PXD, SampleGroup, sep = "__")
  )
if (nrow(stats) != 37) {
  stop("Expected 37 Kla sample groups, found ", nrow(stats))
}
paired_stats <- stats |>
  filter(PairedAnalysisIncluded %in% c(TRUE, "TRUE", "True", 1, "1"))
if (nrow(paired_stats) != 33) {
  stop(
    "Expected 33 exact quantitative reference groups, found ",
    nrow(paired_stats)
  )
}
if (!identical(
  unique(stats$Category),
  category_order
)) {
  stop("The four-category order in the statistics table is invalid")
}
category_by_group <- setNames(stats$Category, stats$SampleGroupKey)

go <- read.delim(
  go_path,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = ""
)
go_taxon <- suppressWarnings(as.numeric(go$`TAXON ID`))
go_not <- grepl("(^|\\|)NOT($|\\|)", go$QUALIFIER)
go_keep <- (is.na(go_taxon) | go_taxon == 9606) &
  !go_not &
  go$`GENE PRODUCT DB` == "UniProtKB"
go_ddr_ids <- sort(unique(
  clean_id_vector(go$`GENE PRODUCT ID`[go_keep])
))
if (!length(go_ddr_ids)) stop("No human GO-DDR UniProt IDs were parsed")

all_uniprot <- read_tsv(
  all_uniprot_path,
  show_col_types = FALSE
)
source_symbol_parts <- list()
source_kla_path <- file.path(
  project_root,
  "work", "intermediate",
  "kla_by_dataset",
  "all_primary_sample_level_kla_sites.csv"
)
source_reference_path <- file.path(
  project_root,
  "results", "tables",
  "reference_proteome_all_proteins.csv"
)
if (file.exists(source_kla_path)) {
  source_symbol_parts[[length(source_symbol_parts) + 1]] <-
    read_csv(source_kla_path, show_col_types = FALSE) |>
    transmute(
      BaseAccession = base_accession(BaseAccession),
      SourceGeneSymbol = as.character(GeneSymbol),
      SourceProteinName = as.character(ProteinName),
      SourceSymbolSource = "Original Kla sample-level table; audit only"
    )
}
if (file.exists(source_reference_path)) {
  source_symbol_parts[[length(source_symbol_parts) + 1]] <-
    read_csv(source_reference_path, show_col_types = FALSE) |>
    transmute(
      BaseAccession = base_accession(BaseAccession),
      SourceGeneSymbol = as.character(GeneSymbol),
      SourceProteinName = as.character(ProteinName),
      SourceSymbolSource = "Original reference-proteome table; audit only"
    )
}
source_annotation <- bind_rows(source_symbol_parts) |>
  filter(
    !is.na(BaseAccession),
    nzchar(BaseAccession),
    (
      (!is.na(SourceGeneSymbol) & nzchar(SourceGeneSymbol)) |
        (!is.na(SourceProteinName) & nzchar(SourceProteinName))
    )
  ) |>
  group_by(BaseAccession) |>
  summarise(
    SourceGeneSymbol = paste(
      sort(unique(SourceGeneSymbol[
        !is.na(SourceGeneSymbol) & nzchar(SourceGeneSymbol)
      ])),
      collapse = ";"
    ),
    SourceProteinName = paste(
      sort(unique(SourceProteinName[
        !is.na(SourceProteinName) & nzchar(SourceProteinName)
      ])),
      collapse = ";"
    ),
    SourceSymbolSource = paste(
      sort(unique(SourceSymbolSource)),
      collapse = ";"
    ),
    .groups = "drop"
  )

uniprot_annotation <- all_uniprot |>
  transmute(
    BaseAccession = base_accession(Entry),
    PrimaryGeneSymbol = str_trim(
      str_split_fixed(
        coalesce(`Gene Names (primary)`, ""),
        ";",
        2
      )[, 1]
    ),
    GeneNamesFallback = str_trim(
      str_split_fixed(
        coalesce(`Gene Names`, ""),
        "[;[:space:]]+",
        2
      )[, 1]
    ),
    UniProtProteinName = coalesce(`Protein names`, ""),
    ReviewedStatus = ifelse(
      tolower(coalesce(Reviewed, "")) == "reviewed",
      "reviewed",
      "non-reviewed"
    )
  ) |>
  mutate(
    UniProtGeneSymbol = ifelse(
      nzchar(PrimaryGeneSymbol),
      PrimaryGeneSymbol,
      GeneNamesFallback
    ),
    UniProtGeneSymbolSource = ifelse(
      nzchar(PrimaryGeneSymbol),
      "UniProt primary gene name",
      "UniProt gene-names fallback"
    )
  ) |>
  filter(
    !is.na(BaseAccession),
    nzchar(BaseAccession)
  ) |>
  select(
    BaseAccession,
    UniProtGeneSymbol,
    UniProtProteinName,
    ReviewedStatus,
    UniProtGeneSymbolSource
  ) |>
  distinct(BaseAccession, .keep_all = TRUE)
reviewed_annotation <- full_join(
  uniprot_annotation,
  source_annotation,
  by = "BaseAccession"
) |>
  mutate(
    GeneSymbolAudit = ifelse(
      !is.na(UniProtGeneSymbol) & nzchar(UniProtGeneSymbol),
      UniProtGeneSymbol,
      SourceGeneSymbol
    ),
    ProteinNameAudit = ifelse(
      !is.na(UniProtProteinName) & nzchar(UniProtProteinName),
      UniProtProteinName,
      SourceProteinName
    ),
    ReviewedStatus = ifelse(
      !is.na(ReviewedStatus),
      ReviewedStatus,
      "source-only"
    ),
    GeneSymbolAuditSource = case_when(
      !is.na(UniProtGeneSymbol) & nzchar(UniProtGeneSymbol) ~
        paste0("UniProt ", UniProtGeneSymbolSource),
      !is.na(SourceGeneSymbol) & nzchar(SourceGeneSymbol) ~
        SourceSymbolSource,
      !is.na(UniProtProteinName) & nzchar(UniProtProteinName) ~
        "UniProt protein name only; no gene symbol in UniProt/source",
      TRUE ~ "No matching human UniProt or source annotation; audit only"
    ),
    AnnotationMappingSource = case_when(
      !is.na(UniProtGeneSymbol) & nzchar(UniProtGeneSymbol) ~ paste0(
        "UniProtKB Homo sapiens ",
        ReviewedStatus,
        " ",
        UniProtGeneSymbolSource,
        "; audit only"
      ),
      !is.na(SourceGeneSymbol) & nzchar(SourceGeneSymbol) ~
        paste0(SourceSymbolSource, "; UniProt ID not matched in 2026-08-06 snapshot"),
      !is.na(UniProtProteinName) & nzchar(UniProtProteinName) ~ paste0(
        "UniProtKB Homo sapiens ",
        ReviewedStatus,
        " protein name only; audit only"
      ),
      TRUE ~ "No matching human UniProt or source annotation; audit only"
    )
  ) |>
  filter(
    !is.na(BaseAccession),
    nzchar(BaseAccession),
    (
      (!is.na(GeneSymbolAudit) & nzchar(GeneSymbolAudit)) |
        (!is.na(ProteinNameAudit) & nzchar(ProteinNameAudit))
    )
  ) |>
  select(
    BaseAccession,
    GeneSymbolAudit,
    ProteinNameAudit,
    ReviewedStatus,
    GeneSymbolAuditSource,
    AnnotationMappingSource
  ) |>
  distinct(BaseAccession, .keep_all = TRUE)

kla <- read_csv(kla_path, show_col_types = FALSE) |>
  mutate(
    SampleGroupKey = paste(PXD, SampleGroup, sep = "__"),
    Category = unname(category_by_group[SampleGroupKey]),
    BaseAccession = base_accession(BaseAccession)
  ) |>
  filter(
    !is.na(Category),
    !is.na(BaseAccession),
    nzchar(BaseAccession)
  ) |>
  semi_join(
    paired_stats |>
      select(SampleGroupKey),
    by = "SampleGroupKey"
  ) |>
  distinct(PXD, SampleGroup, Category, BaseAccession)

reference <- read_csv(reference_path, show_col_types = FALSE) |>
  transmute(
    PXD,
    SampleGroup,
    SampleGroupKey = paste(PXD, SampleGroup, sep = "__"),
    Category = unname(category_by_group[SampleGroupKey]),
    MappedBaseAccessions
  ) |>
  filter(!is.na(Category)) |>
  separate_rows(MappedBaseAccessions, sep = ";") |>
  mutate(BaseAccession = base_accession(MappedBaseAccessions)) |>
  filter(!is.na(BaseAccession), nzchar(BaseAccession)) |>
  distinct(PXD, SampleGroup, Category, BaseAccession)
reference_group_keys <- unique(paste(reference$PXD, reference$SampleGroup, sep = "__"))
expected_reference_group_keys <- unique(paired_stats$SampleGroupKey)
if (
  length(reference_group_keys) != 33 ||
    length(setdiff(reference_group_keys, expected_reference_group_keys)) ||
    length(setdiff(expected_reference_group_keys, reference_group_keys))
  ) {
  stop(
    "Reference Venn membership does not match the 33 exact quantitative ",
    "reference groups"
  )
}
venn_sample_scope <- stats |>
  transmute(
    PXD,
    SampleGroup,
    SampleGroupKey,
    Category,
    KlaIncludedInVenn = SampleGroupKey %in% paired_stats$SampleGroupKey,
    ReferenceIncludedInVenn = SampleGroupKey %in% expected_reference_group_keys,
    ExclusionReason = ifelse(
      KlaIncludedInVenn,
      "",
      "没有完全匹配且可审计的普通全蛋白参照；不进入当前Kla/参照Venn"
    )
  )
write_csv(
  venn_sample_scope,
  file.path(table_root, "venn_sample_group_scope.csv")
)
write_csv(
  paired_stats |>
    select(-SampleGroupKey),
  file.path(
    project_root,
    "results", "tables",
    "cell_type_kla_vs_reference_ddr_statistics_accession_only_paired_33.csv"
  )
)
paired_stats_zh <- read_csv(
  file.path(
    project_root,
    "results", "tables",
    "cell_type_kla_vs_reference_ddr_statistics_accession_only_zh.csv"
  ),
  show_col_types = FALSE
) |>
  mutate(SampleGroupKey = paste(`乳酸化PXD`, `样本组`, sep = "__")) |>
  semi_join(paired_stats |> select(SampleGroupKey), by = "SampleGroupKey") |>
  select(-SampleGroupKey)
write_csv(
  paired_stats_zh,
  file.path(
    project_root,
    "results", "tables",
    "cell_type_kla_vs_reference_ddr_statistics_accession_only_paired_33_zh.csv"
  )
)

make_sets <- function(data, ddr_only = FALSE) {
  if (ddr_only) data <- data |> filter(BaseAccession %in% go_ddr_ids)
  sets <- lapply(
    category_order,
    function(category) sort(unique(
      data$BaseAccession[data$Category == category]
    ))
  )
  names(sets) <- category_order
  sets
}

region_name <- function(present) {
  active <- category_order[present]
  if (!length(active)) return("outside_all_sets")
  if (length(active) == 1) return(paste0(active, "_only"))
  if (length(active) == length(category_order)) return("all_four")
  paste0(paste(active, collapse = "_and_"), "_only")
}

build_membership <- function(sets) {
  all_ids <- sort(unique(unlist(sets)))
  membership <- data.frame(
    BaseAccession = all_ids,
    stringsAsFactors = FALSE
  )
  for (category in category_order) {
    membership[[paste0("In_", category)]] <-
      membership$BaseAccession %in% sets[[category]]
  }
  membership$Region <- apply(
    membership[paste0("In_", category_order)],
    1,
    function(values) region_name(as.logical(values))
  )
  membership
}

write_analysis_tables <- function(sets, membership, analysis_name) {
  output_dir <- file.path(table_root, analysis_name)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  set_counts <- data.frame(
    Category = category_order,
    CategoryZh = unname(category_labels_zh[category_order]),
    CategoryEn = unname(category_labels_en[category_order]),
    ProteinCount = vapply(sets, length, integer(1)),
    stringsAsFactors = FALSE
  )
  region_counts <- membership |>
    count(Region, name = "ProteinCount") |>
    right_join(
      data.frame(
        Region = c(
          paste0(category_order, "_only"),
          combn(category_order, 2, paste, collapse = "_and_") |>
            paste0("_only"),
          combn(category_order, 3, paste, collapse = "_and_") |>
            paste0("_only"),
          "all_four"
        ),
        stringsAsFactors = FALSE
      ),
      by = "Region"
    ) |>
    mutate(ProteinCount = replace_na(ProteinCount, 0L)) |>
    arrange(match(Region, c(
      paste0(category_order, "_only"),
      combn(category_order, 2, paste, collapse = "_and_") |>
        paste0("_only"),
      combn(category_order, 3, paste, collapse = "_and_") |>
        paste0("_only"),
      "all_four"
    )))

  membership_annotated <- membership |>
    left_join(reviewed_annotation, by = "BaseAccession") |>
    mutate(
      GeneSymbolAudit = ifelse(
        !is.na(GeneSymbolAudit) & nzchar(GeneSymbolAudit),
        GeneSymbolAudit
        ,
        ifelse(
          !is.na(ProteinNameAudit) & nzchar(ProteinNameAudit),
          "[no gene symbol in UniProt/source]",
          "[unmapped]"
        )
      ),
      ProteinNameAudit = ifelse(
        is.na(ProteinNameAudit) | !nzchar(ProteinNameAudit),
        "[unmapped]",
        ProteinNameAudit
      ),
      ReviewedStatus = ifelse(
        is.na(ReviewedStatus),
        "unmapped",
        ReviewedStatus
      ),
      GeneSymbolAuditSource = ifelse(
        is.na(GeneSymbolAuditSource),
        "No matching human UniProt annotation; audit only",
        GeneSymbolAuditSource
      ),
      AnnotationMappingSource = ifelse(
        is.na(AnnotationMappingSource),
        "No matching human UniProt annotation; audit only",
        AnnotationMappingSource
      )
    ) |>
    select(
      BaseAccession,
      GeneSymbolAudit,
      ProteinNameAudit,
      ReviewedStatus,
      GeneSymbolAuditSource,
      AnnotationMappingSource,
      everything()
    )

  write_csv(set_counts, file.path(output_dir, "set_counts.csv"))
  write_csv(region_counts, file.path(output_dir, "region_counts.csv"))
  write_csv(
    membership_annotated,
    file.path(output_dir, "membership.csv")
  )
  write_csv(
    membership_annotated |>
      select(
        BaseAccession,
        GeneSymbolAudit,
        ProteinNameAudit,
        ReviewedStatus,
        GeneSymbolAuditSource,
        AnnotationMappingSource
      ) |>
      distinct(),
    file.path(output_dir, "id_annotation_mapping.csv")
  )
  invisible(set_counts)
}

draw_area_proportional <- function(
  sets,
  analysis_name,
  title_zh,
  title_en
) {
  fit <- euler(sets)
  membership <- build_membership(sets)
  write_analysis_tables(sets, membership, analysis_name)

  fit_regions <- data.frame(
    Region = names(fit$original.values),
    Original = as.numeric(fit$original.values),
    Fitted = as.numeric(fit$fitted.values),
    Residual = as.numeric(fit$residuals),
    stringsAsFactors = FALSE
  )
  write_csv(
    fit_regions,
    file.path(table_root, analysis_name, "euler_fit_regions.csv")
  )
  write_csv(
    data.frame(
      Analysis = analysis_name,
      Stress = fit$stress,
      stringsAsFactors = FALSE
    ),
    file.path(table_root, analysis_name, "euler_fit_summary.csv")
  )

  colors <- c(
    "#E69F00",
    "#D55E00",
    "#009E73",
    "#CC79A7"
  )
  render <- function(language = c("zh", "en")) {
    language <- match.arg(language)
    labels <- if (language == "zh") {
      unname(category_labels_zh[category_order])
    } else {
      unname(category_labels_en[category_order])
    }
    title <- if (language == "zh") title_zh else title_en
    font_family <- if (language == "zh") "PingFang SC" else "Arial"
    scope_title <- if (language == "zh") {
      "严格配对33组（9/2/9/13；排除4组）"
    } else {
      "33 strict pairs (9/2/9/13; 4 unpaired groups excluded)"
    }
    visible_title <- paste(title, scope_title, sep = " | ")
    subtitle <- if (language == "zh") {
      paste0(
        "严格配对范围：33组（正常组织9、癌症组织2、正常细胞9、癌症细胞13）；",
        "4组无严格参照已排除\n",
        "按 UniProt BaseAccession 去重；椭圆面积按集合数量比例拟合"
      )
    } else {
      paste0(
        "Strict paired scope: 33 groups (normal tissue 9, cancer tissue 2, ",
        "normal cells 9, cancer cells 13); 4 groups without exact references excluded\n",
        "Deduplicated by UniProt BaseAccession; ellipse areas are fitted proportional to set sizes"
      )
    }
    output_stems <- file.path(
      figure_root,
      c(
        paste0(analysis_name, "_", language),
        paste0(analysis_name, "_33groups_", language)
      )
    )
    for (output_stem in output_stems) {
      png(
        paste0(output_stem, ".png"),
        width = 2400,
        height = 2300,
        res = 300,
        type = "cairo"
      )
      par(family = font_family, mar = c(1, 1, 8, 1))
      diagram <- plot(
        fit,
        labels = FALSE,
        legend = list(
          labels = labels,
          side = "bottom",
          nrow = 1,
          ncol = 4,
          byrow = TRUE,
          cex = 0.9
        ),
        quantities = list(cex = 0.9),
        fills = list(fill = colors, alpha = 0.52),
        edges = list(col = "#4B4B4B", lwd = 1.2),
        main = list(label = visible_title, cex = 0.85),
        sub = subtitle,
        sub.cex = 0.78,
        quantities.cex = 1.0
      )
      print(diagram)
      dev.off()

      cairo_pdf(
        paste0(output_stem, ".pdf"),
        width = 8.5,
        height = 8.0,
        family = font_family
      )
      par(family = font_family, mar = c(1, 1, 8, 1))
      diagram <- plot(
        fit,
        labels = FALSE,
        legend = list(
          labels = labels,
          side = "bottom",
          nrow = 1,
          ncol = 4,
          byrow = TRUE,
          cex = 0.9
        ),
        quantities = list(cex = 0.9),
        fills = list(fill = colors, alpha = 0.52),
        edges = list(col = "#4B4B4B", lwd = 1.2),
        main = list(label = visible_title, cex = 0.85),
        sub = subtitle,
        sub.cex = 0.78,
        quantities.cex = 1.0
      )
      print(diagram)
      dev.off()
    }
    invisible(NULL)
  }
  render("zh")
  render("en")
  invisible(fit)
}

draw_area_proportional(
  make_sets(kla, ddr_only = FALSE),
  "all_kla_four_class_venn",
  "四类组织/细胞中的全部乳酸化蛋白（Kla）",
  "All Kla proteins across four tissue and cell categories"
)
draw_area_proportional(
  make_sets(kla, ddr_only = TRUE),
  "kla_ddr_four_class_venn",
  "四类组织/细胞中的乳酸化 DDR 蛋白",
  "Kla and DDR proteins across four tissue and cell categories"
)
draw_area_proportional(
  make_sets(reference, ddr_only = FALSE),
  "reference_proteome_four_class_venn",
  "四类组织/细胞中的普通全蛋白组蛋白",
  "Reference whole-proteome proteins across four tissue and cell categories"
)
draw_area_proportional(
  make_sets(reference, ddr_only = TRUE),
  "reference_proteome_ddr_four_class_venn",
  "四类组织/细胞中的普通全蛋白组 DDR 蛋白",
  "Reference whole-proteome DDR proteins across four tissue and cell categories"
)

message(
  "Generated four area-proportional Venn/Euler analyses in ",
  figure_root
)
