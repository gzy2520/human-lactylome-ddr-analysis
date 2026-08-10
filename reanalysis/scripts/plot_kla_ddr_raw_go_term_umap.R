#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(uwot)
})

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
script_path <- if (length(file_arg) > 0L) {
  normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE)
} else {
  normalizePath("reanalysis/scripts/plot_kla_ddr_raw_go_term_umap.R", mustWork = TRUE)
}
project_root <- normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = TRUE)

kla_membership_path <- file.path(
  project_root,
  "reanalysis/intermediate/expanded_ddr_by_accession/kla_proteins_by_sample_group.csv"
)
exclusions_path <- file.path(
  project_root,
  "reanalysis/config/final_sample_group_exclusions.csv"
)
four_class_path <- file.path(
  project_root,
  "reanalysis/config/four_class_sample_grouping.csv"
)
go_path <- file.path(
  project_root,
  "data/annotations/GO-repair+damage(human).tsv"
)

table_dir <- file.path(
  project_root,
  "reanalysis/results/tables/kla_ddr_raw_go_umap_33groups"
)
figure_dir <- file.path(
  project_root,
  "reanalysis/results/figures/kla_ddr_raw_go_umap_33groups"
)
report_path <- file.path(
  project_root,
  "reanalysis/reports/UMAP_RAW_GO_33GROUP_DATA_SCOPE.md"
)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

required_paths <- c(
  kla_membership_path,
  exclusions_path,
  four_class_path,
  go_path
)
missing_paths <- required_paths[!file.exists(required_paths)]
if (length(missing_paths) > 0L) {
  stop("Missing required input(s): ", paste(missing_paths, collapse = "; "))
}

stop_if_false <- function(condition, message) {
  if (!isTRUE(condition)) {
    stop(message, call. = FALSE)
  }
}

group_key <- function(pxd, sample_group) {
  paste(trimws(as.character(pxd)), trimws(as.character(sample_group)), sep = "\r")
}

base_accession <- function(x) {
  x <- trimws(as.character(x))
  x <- sub("^(sp|tr)\\|", "", x)
  x <- sub("\\|.*$", "", x)
  x <- sub("^.*:", "", x)
  sub("-[0-9]+$", "", x)
}

kla_membership <- fread(kla_membership_path)
exclusions <- fread(exclusions_path)
four_class <- fread(four_class_path)
go_raw <- fread(go_path, sep = "\t", encoding = "UTF-8")
setnames(go_raw, sub("^\\ufeff", "", names(go_raw)))

stop_if_false(
  all(c("PXD", "SampleGroup", "BaseAccession", "IsDdr") %in% names(kla_membership)),
  "Kla membership table is missing required columns."
)
stop_if_false(
  all(c("PXD", "SampleGroup", "ExclusionStatus") %in% names(exclusions)),
  "Exclusion table is missing required columns."
)
stop_if_false(
  all(c("PXD", "SampleGroup", "Category") %in% names(four_class)),
  "Four-class scope table is missing required columns."
)
stop_if_false(
  all(
    c(
      "GENE PRODUCT DB", "GENE PRODUCT ID", "QUALIFIER",
      "GO TERM", "GO NAME", "TAXON ID"
    ) %in% names(go_raw)
  ),
  "GO annotation table is missing required columns."
)

source_groups <- unique(kla_membership[, .(PXD, SampleGroup)])
setorder(source_groups, PXD, SampleGroup)

removed_four <- exclusions[
  ExclusionStatus == "removed_from_final_33group_scope",
  .(PXD, SampleGroup, ExclusionStatus, DecisionDate, Reason)
]
stop_if_false(nrow(removed_four) == 4L, "Expected exactly four final-scope removals.")
stop_if_false(
  uniqueN(removed_four, by = c("PXD", "SampleGroup")) == 4L,
  "The four final-scope removals are not unique by PXD + SampleGroup."
)

source_keys <- group_key(source_groups$PXD, source_groups$SampleGroup)
removed_keys <- group_key(removed_four$PXD, removed_four$SampleGroup)
stop_if_false(
  all(removed_keys %chin% source_keys),
  "At least one of the four removed groups is absent from the 37-group source table."
)

kla_membership[, ScopeKey := group_key(PXD, SampleGroup)]
kla_33 <- kla_membership[!ScopeKey %chin% removed_keys]
scope_33 <- unique(kla_33[, .(PXD, SampleGroup)])
setorder(scope_33, PXD, SampleGroup)

four_class[, ScopeKey := group_key(PXD, SampleGroup)]
scope_33[, ScopeKey := group_key(PXD, SampleGroup)]
stop_if_false(
  uniqueN(source_groups, by = c("PXD", "SampleGroup")) == 37L,
  "The source membership table is no longer a 37-group table."
)
stop_if_false(
  uniqueN(scope_33, by = c("PXD", "SampleGroup")) == 33L,
  "Removing the four groups did not yield exactly 33 groups."
)
stop_if_false(
  uniqueN(four_class, by = c("PXD", "SampleGroup")) == 33L,
  "The formal four-class configuration does not contain exactly 33 unique groups."
)
stop_if_false(
  setequal(scope_33$ScopeKey, four_class$ScopeKey),
  "The filtered membership groups do not exactly match the formal 33-group configuration."
)

class_expected <- c(
  normal_cells = 9L,
  cancer_cells = 13L,
  normal_tissue = 9L,
  cancer_tissue = 2L
)
class_observed <- four_class[, .N, by = Category]
stop_if_false(
  setequal(class_observed$Category, names(class_expected)) &&
    all(
      class_observed$N[
        match(names(class_expected), class_observed$Category)
      ] == unname(class_expected)
    ),
  "The formal 33-group class counts are not 9/13/9/2."
)

project_scope <- merge(
  scope_33[, .(PXD, SampleGroup)],
  four_class[, .(PXD, SampleGroup, Category)],
  by = c("PXD", "SampleGroup"),
  all.x = TRUE,
  sort = TRUE
)
project_scope[, IncludedInCurrentProjectScope := TRUE]

removed_audit <- copy(removed_four)
removed_audit[, PresentIn37GroupSource := group_key(PXD, SampleGroup) %chin% source_keys]
removed_counts <- kla_membership[
  ScopeKey %chin% removed_keys,
  .(
    SourceMembershipRows = .N,
    UniqueKlaProteins = uniqueN(BaseAccession),
    UniqueKlaDdrProteins = uniqueN(BaseAccession[IsDdr %in% TRUE])
  ),
  by = .(PXD, SampleGroup)
]
removed_audit <- merge(
  removed_audit,
  removed_counts,
  by = c("PXD", "SampleGroup"),
  all.x = TRUE,
  sort = TRUE
)

kla_ddr_membership <- unique(
  kla_33[IsDdr %in% TRUE, .(PXD, SampleGroup, BaseAccession)]
)
kla_ddr_membership <- kla_ddr_membership[nzchar(trimws(BaseAccession))]
kla_ddr_proteins <- sort(unique(kla_ddr_membership$BaseAccession))
stop_if_false(
  length(kla_ddr_proteins) == 507L,
  "The 33-group Kla intersect DDR union does not contain exactly 507 BaseAccessions."
)

go_raw[, BaseAccession := base_accession(`GENE PRODUCT ID`)]
qualifier_has_not <- grepl(
  "(^|[|,;[:space:]])NOT($|[|,;[:space:]])",
  as.character(go_raw$QUALIFIER),
  ignore.case = TRUE,
  perl = TRUE
)
qualifier_has_not[is.na(qualifier_has_not)] <- FALSE

go_filtered <- go_raw[
  `GENE PRODUCT DB` == "UniProtKB" &
    as.character(`TAXON ID`) == "9606" &
    !qualifier_has_not &
    nzchar(BaseAccession) &
    nzchar(trimws(as.character(`GO TERM`)))
]

protein_go_long <- unique(
  go_filtered[
    BaseAccession %chin% kla_ddr_proteins,
    .(
      BaseAccession,
      GO_TERM = trimws(as.character(`GO TERM`)),
      GO_NAME = trimws(as.character(`GO NAME`))
    )
  ]
)
setorder(protein_go_long, BaseAccession, GO_TERM, GO_NAME)

stop_if_false(
  uniqueN(protein_go_long$BaseAccession) == 507L,
  "At least one of the 507 Kla-DDR proteins lacks a qualifying raw GO term."
)
stop_if_false(
  uniqueN(protein_go_long, by = c("BaseAccession", "GO_TERM")) == nrow(protein_go_long),
  "Protein-GO term pairs were not reduced to one binary hit."
)

term_names_per_id <- protein_go_long[, uniqueN(GO_NAME), by = GO_TERM]
stop_if_false(
  all(term_names_per_id$V1 == 1L),
  "At least one GO term has conflicting GO names in the filtered annotation."
)

go_terms <- sort(unique(protein_go_long$GO_TERM))
stop_if_false(length(go_terms) == 66L, "Expected 66 raw GO terms in the 33-group intersection.")
stop_if_false(nrow(protein_go_long) == 1029L, "Expected 1,029 unique protein-GO term pairs.")

protein_index <- match(protein_go_long$BaseAccession, kla_ddr_proteins)
term_index <- match(protein_go_long$GO_TERM, go_terms)
umap_matrix <- matrix(
  0,
  nrow = length(kla_ddr_proteins),
  ncol = length(go_terms),
  dimnames = list(kla_ddr_proteins, go_terms)
)
umap_matrix[cbind(protein_index, term_index)] <- 1

stop_if_false(all(umap_matrix %in% c(0, 1)), "UMAP input matrix is not binary.")
stop_if_false(all(rowSums(umap_matrix) >= 1L), "At least one protein has no GO-term feature.")
stop_if_false(
  sum(umap_matrix) == nrow(protein_go_long),
  "Binary-matrix hits do not match the protein-GO long table."
)

protein_patterns <- data.table(
  BaseAccession = kla_ddr_proteins,
  RawGOTermPattern = apply(
    umap_matrix,
    1L,
    function(x) paste(go_terms[x == 1], collapse = ";")
  )
)
pattern_summary <- protein_patterns[
  ,
  .(
    NumberOfProteins = .N,
    BaseAccessions = paste(sort(BaseAccession), collapse = ";")
  ),
  by = RawGOTermPattern
]
setorder(pattern_summary, -NumberOfProteins, RawGOTermPattern)
pattern_summary[, PatternID := sprintf("PATTERN_%03d", seq_len(.N))]
setcolorder(
  pattern_summary,
  c("PatternID", "RawGOTermPattern", "NumberOfProteins", "BaseAccessions")
)
protein_patterns <- merge(
  protein_patterns,
  pattern_summary[, .(PatternID, RawGOTermPattern)],
  by = "RawGOTermPattern",
  all.x = TRUE,
  sort = FALSE
)
stop_if_false(nrow(pattern_summary) == 205L, "Expected 205 distinct raw-GO membership patterns.")
stop_if_false(
  pattern_summary$NumberOfProteins[[1L]] == 102L &&
    pattern_summary$RawGOTermPattern[[1L]] == "GO:0006974",
  "The largest repeated raw-GO pattern is no longer 102 proteins annotated only to GO:0006974."
)

# Keep the broad visual geometry of the previous Kla-DDR UMAP
# (30 neighbors, min_dist 0.8, spread 3.0), but use cosine distance because
# this input is a sparse binary raw-GO-term membership matrix. Random
# initialization plus a fixed seed avoids an excessively elongated placement
# of disconnected annotation components while preserving reproducibility.
umap_seed <- 25L
umap_neighbors <- 30L
umap_min_dist <- 0.8
umap_spread <- 3.0
umap_epochs <- 500L

set.seed(umap_seed)
umap_coordinates_matrix <- uwot::umap(
  X = umap_matrix,
  n_neighbors = umap_neighbors,
  n_components = 2L,
  metric = "cosine",
  n_epochs = umap_epochs,
  scale = FALSE,
  init = "random",
  spread = umap_spread,
  min_dist = umap_min_dist,
  fast_sgd = FALSE,
  n_threads = 1L,
  n_sgd_threads = 1L,
  seed = umap_seed,
  verbose = FALSE
)

stop_if_false(
  identical(dim(umap_coordinates_matrix), c(507L, 2L)),
  "UMAP did not return a 507 x 2 coordinate matrix."
)
stop_if_false(all(is.finite(umap_coordinates_matrix)), "UMAP coordinates contain non-finite values.")

umap_coordinates <- data.table(
  BaseAccession = kla_ddr_proteins,
  UMAP_1 = umap_coordinates_matrix[, 1L],
  UMAP_2 = umap_coordinates_matrix[, 2L]
)

protein_scope <- kla_ddr_membership[
  ,
  .(
    NumberOfIncludedGroups = uniqueN(group_key(PXD, SampleGroup)),
    IncludedPXD = paste(sort(unique(PXD)), collapse = ";"),
    IncludedSampleGroups = paste(
      sort(unique(paste(PXD, SampleGroup, sep = "::"))),
      collapse = ";"
    )
  ),
  by = BaseAccession
]
go_counts <- protein_go_long[, .(NumberOfRawGOTerms = uniqueN(GO_TERM)), by = BaseAccession]
protein_scope <- merge(protein_scope, go_counts, by = "BaseAccession", all.x = TRUE, sort = TRUE)
protein_scope <- merge(
  protein_scope,
  protein_patterns[, .(BaseAccession, PatternID)],
  by = "BaseAccession",
  all.x = TRUE,
  sort = TRUE
)

go_term_summary <- protein_go_long[
  ,
  .(
    GO_NAME = first(GO_NAME),
    NumberOfKlaDdrProteins = uniqueN(BaseAccession)
  ),
  by = GO_TERM
]
setorder(go_term_summary, -NumberOfKlaDdrProteins, GO_TERM)

binary_output <- data.table(BaseAccession = kla_ddr_proteins)
binary_output <- cbind(binary_output, as.data.table(umap_matrix))

scope_summary <- data.table(
  Item = c(
    "Source sample groups",
    "Removed unusable groups",
    "Final included sample groups",
    "Normal-cell groups",
    "Cancer-cell groups",
    "Normal-tissue groups",
    "Cancer-tissue groups",
    "Unique 33-group Kla-DDR BaseAccessions",
    "Qualifying raw GO terms",
    "Unique protein-GO term pairs",
    "Distinct raw-GO membership patterns",
    "Largest identical raw-GO pattern",
    "Proteins lacking a qualifying GO term"
  ),
  Value = c(
    37L,
    4L,
    33L,
    class_expected[["normal_cells"]],
    class_expected[["cancer_cells"]],
    class_expected[["normal_tissue"]],
    class_expected[["cancer_tissue"]],
    length(kla_ddr_proteins),
    length(go_terms),
    nrow(protein_go_long),
    nrow(pattern_summary),
    pattern_summary$NumberOfProteins[[1L]],
    length(setdiff(kla_ddr_proteins, protein_go_long$BaseAccession))
  )
)

parameters <- data.table(
  Parameter = c(
    "AnalysisUnit",
    "ProteinAnalysisKey",
    "UMAPInput",
    "GOFeatureEncoding",
    "GOCategoryAggregation",
    "GeneSymbolUsedAsAnalysisKey",
    "DistanceMetric",
    "NNeighbors",
    "MinDist",
    "Spread",
    "NEpochs",
    "Initialization",
    "ScaleFeatures",
    "RandomSeed",
    "NearestNeighborThreads",
    "SGDThreads",
    "RVersion",
    "uwotVersion",
    "ggplot2Version",
    "data.tableVersion"
  ),
  Value = c(
    "one unique protein per point",
    "isoform-stripped UniProt BaseAccession",
    "protein x raw GO term binary matrix only",
    "presence=1; absence=0",
    "none",
    "FALSE",
    "cosine",
    umap_neighbors,
    umap_min_dist,
    umap_spread,
    umap_epochs,
    "random",
    "FALSE",
    umap_seed,
    1L,
    1L,
    paste(R.version$major, R.version$minor, sep = "."),
    as.character(packageVersion("uwot")),
    as.character(packageVersion("ggplot2")),
    as.character(packageVersion("data.table"))
  )
)

input_file_audit <- data.table(
  InputRole = c(
    "37-group Kla membership source",
    "central exclusion decisions",
    "formal 33-group four-class scope",
    "human GO annotation"
  ),
  Path = sub(
    paste0("^", gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", project_root), "/?"),
    "",
    required_paths
  ),
  MD5 = unname(tools::md5sum(required_paths))
)

fwrite(project_scope, file.path(table_dir, "project_sample_group_scope_33.csv"))
fwrite(removed_audit, file.path(table_dir, "removed_sample_groups_4.csv"))
fwrite(scope_summary, file.path(table_dir, "scope_summary.csv"))
fwrite(protein_scope, file.path(table_dir, "kla_ddr_protein_scope.csv"))
fwrite(protein_go_long, file.path(table_dir, "protein_raw_go_term_long.csv"))
fwrite(go_term_summary, file.path(table_dir, "raw_go_term_summary.csv"))
fwrite(pattern_summary, file.path(table_dir, "raw_go_pattern_summary.csv"))
fwrite(binary_output, file.path(table_dir, "protein_raw_go_term_binary_matrix.csv"))
fwrite(umap_coordinates, file.path(table_dir, "umap_coordinates_fixed.csv"))
fwrite(parameters, file.path(table_dir, "umap_parameters.csv"))
fwrite(input_file_audit, file.path(table_dir, "input_file_audit.csv"))

saveRDS(
  list(
    coordinates = umap_coordinates,
    input_matrix = umap_matrix,
    protein_go_long = protein_go_long,
    pattern_summary = pattern_summary,
    parameters = parameters
  ),
  file.path(table_dir, "umap_raw_go_33groups_bundle.rds"),
  version = 3
)

writeLines(
  capture.output(sessionInfo()),
  file.path(table_dir, "session_info.txt"),
  useBytes = TRUE
)

umap_plot <- ggplot(umap_coordinates, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(
    color = "#2C7FB8",
    size = 2.15,
    alpha = 0.90,
    shape = 16
  ) +
  coord_equal() +
  labs(
    title = "Kla-DDR proteins embedded by raw GO-term membership",
    subtitle = "33 sample groups | 507 UniProt BaseAccessions | 66 ungrouped GO terms",
    x = "UMAP 1",
    y = "UMAP 2",
    caption = paste0(
      "One point per protein; input: 507 × 66 binary matrix; ",
      "cosine distance; n_neighbors = 30; seed = 25."
    )
  ) +
  theme_classic(base_size = 11, base_family = "sans") +
  theme(
    axis.line = element_line(linewidth = 0.45, colour = "#333333"),
    axis.ticks = element_line(linewidth = 0.4, colour = "#333333"),
    axis.text = element_text(colour = "#333333"),
    plot.title = element_text(face = "bold", size = 13, hjust = 0),
    plot.subtitle = element_text(size = 10, colour = "#4D4D4D"),
    plot.caption = element_text(size = 8.5, colour = "#5A5A5A", hjust = 0),
    plot.margin = margin(10, 13, 10, 10)
  )

png_path <- file.path(figure_dir, "kla_ddr_raw_go_umap_33groups_solid_blue.png")
pdf_path <- file.path(figure_dir, "kla_ddr_raw_go_umap_33groups_solid_blue.pdf")
svg_path <- file.path(figure_dir, "kla_ddr_raw_go_umap_33groups_solid_blue.svg")

ggsave(
  png_path,
  umap_plot,
  width = 8.2,
  height = 6.6,
  units = "in",
  dpi = 600,
  bg = "white",
  limitsize = TRUE
)
ggsave(
  pdf_path,
  umap_plot,
  width = 8.2,
  height = 6.6,
  units = "in",
  bg = "white",
  limitsize = TRUE
)
grDevices::svg(svg_path, width = 8.2, height = 6.6, onefile = FALSE, bg = "white")
print(umap_plot)
grDevices::dev.off()

removed_lines <- paste0(
  seq_len(nrow(removed_audit)),
  ". `",
  removed_audit$PXD,
  " / ",
  removed_audit$SampleGroup,
  "`"
)

report_lines <- c(
  "# 33组项目范围与Kla∩DDR原始GO-term UMAP数据审计",
  "",
  "## 当前项目正式范围",
  "",
  "- 原始Kla成员表含37个可量化样本组。",
  "- 本次及后续正式分析按`PXD + SampleGroup`精确排除4个不可用组，得到33组。",
  "- 正式33组与`four_class_sample_grouping.csv`逐组完全一致：正常细胞9组、癌细胞13组、正常组织9组、癌组织2组。",
  "- `PXD037371`的3个临床组此前已因TMT通道无法可靠映射而排除；它们不在上述37组成员源表中，也不计入本次“37减4等于33”。",
  "",
  "本次从37组中排除的4组为：",
  "",
  removed_lines,
  "",
  "## 本次UMAP选择的数据",
  "",
  "- 分析对象：33组中所有`IsDdr == TRUE`的Kla蛋白并集。",
  "- 唯一分析键：去除isoform后缀的UniProt `BaseAccession`；不以Gene Symbol匹配、去重或建模。",
  "- 一个点代表一个唯一`BaseAccession`，共507个点。",
  "- GO来源：`data/annotations/GO-repair+damage(human).tsv`。",
  "- GO过滤：仅保留`GENE PRODUCT DB == UniProtKB`、人源`TAXON ID == 9606`、不含`NOT`限定词且GO term非空的记录。",
  "- 同一蛋白与同一GO term的重复记录折叠为一个二值命中。",
  "- 507个蛋白全部至少有1个合格GO term；最终为66个原始GO term、1,029个唯一蛋白–GO term配对。",
  "- 507个蛋白形成205种不同的原始GO成员模式；最大的一组为102个仅命中`GO:0006974`的蛋白。",
  "- UMAP的唯一输入为507 × 66的“蛋白 × 原始GO term”0/1矩阵。",
  "- 不对GO term做HR、NHEJ、BER等类别归并；样本类别、蛋白出现组数、表达/强度和Gene Symbol均不进入UMAP特征矩阵。",
  "",
  "## UMAP参数与后续使用",
  "",
  "- R包：`uwot`。",
  "- 距离：cosine，适用于当前稀疏二值成员矩阵。",
  "- 参数：`n_neighbors = 30`、`min_dist = 0.8`、`spread = 3.0`、`n_epochs = 500`。",
  "- 使用随机初始化并固定随机种子25；最近邻与SGD均单线程，以保证当前环境内可重复。",
  "- 当前图中所有点为同一种蓝色，不编码通路。",
  "- 后续获得通路归属后，应直接按`BaseAccession`合并到`umap_coordinates_fixed.csv`并在固定坐标上将点替换为饼图；不得重新拟合UMAP。",
  "- 具有完全相同GO成员模式的蛋白在输入空间中不可区分；它们在二维图内的小幅分散不代表额外生物学差异。",
  "",
  "## 输出边界",
  "",
  "- 本流程只新建`kla_ddr_raw_go_umap_33groups`目录及本审计文档。",
  "- 未重跑、未覆盖既有Venn图、DDR柱状图或`previous_umap`。",
  "- UMAP是基于GO注释相似性的二维可视化；坐标轴及点间全局距离不应解释为定量生物学效应。",
  "",
  "## 关键文件",
  "",
  "- 正式33组：`reanalysis/results/tables/kla_ddr_raw_go_umap_33groups/project_sample_group_scope_33.csv`",
  "- 4个排除组：`reanalysis/results/tables/kla_ddr_raw_go_umap_33groups/removed_sample_groups_4.csv`",
  "- 507个蛋白：`reanalysis/results/tables/kla_ddr_raw_go_umap_33groups/kla_ddr_protein_scope.csv`",
  "- 原始GO长表：`reanalysis/results/tables/kla_ddr_raw_go_umap_33groups/protein_raw_go_term_long.csv`",
  "- 原始GO模式汇总：`reanalysis/results/tables/kla_ddr_raw_go_umap_33groups/raw_go_pattern_summary.csv`",
  "- 507 × 66二值矩阵：`reanalysis/results/tables/kla_ddr_raw_go_umap_33groups/protein_raw_go_term_binary_matrix.csv`",
  "- 固定坐标：`reanalysis/results/tables/kla_ddr_raw_go_umap_33groups/umap_coordinates_fixed.csv`",
  "- UMAP参数：`reanalysis/results/tables/kla_ddr_raw_go_umap_33groups/umap_parameters.csv`",
  "- 单色图：`reanalysis/results/figures/kla_ddr_raw_go_umap_33groups/kla_ddr_raw_go_umap_33groups_solid_blue.{png,pdf,svg}`"
)
writeLines(report_lines, report_path, useBytes = TRUE)

message("Formal sample groups: ", nrow(project_scope))
message("Kla-DDR BaseAccessions: ", length(kla_ddr_proteins))
message("Raw GO terms: ", length(go_terms))
message("Protein-GO pairs: ", nrow(protein_go_long))
message("Fixed coordinates: ", file.path(table_dir, "umap_coordinates_fixed.csv"))
message("Figure: ", png_path)
message("Scope report: ", report_path)
