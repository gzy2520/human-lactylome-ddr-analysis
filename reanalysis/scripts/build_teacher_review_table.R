#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args)) normalizePath(args[[1]]) else normalizePath(".")
table_dir <- file.path(project_root, "reanalysis", "results", "tables")
config_dir <- file.path(project_root, "reanalysis", "config")
output_csv <- file.path(table_dir, "kla_and_reference_teacher_review_zh.csv")
focus_csv <- file.path(table_dir, "kla_and_reference_teacher_review_focus_zh.csv")

read_project_csv <- function(path) {
  read.csv(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    na.strings = c("", "NA")
  )
}

pairing <- read_project_csv(
  file.path(table_dir, "lactylome_and_reference_proteome_pairing_zh.csv")
)
statistics <- read_project_csv(
  file.path(
    table_dir,
    "cell_type_kla_vs_reference_ddr_statistics_accession_only_zh.csv"
  )
)
audit <- read_project_csv(
  file.path(table_dir, "strict_reference_material_identity_audit_zh.csv")
)

key <- function(data) paste(data$乳酸化PXD, data$样本组, sep = "\r")
if (nrow(statistics) != 33 || anyDuplicated(key(statistics))) {
  stop("Teacher review table requires exactly 33 unique Kla PXD + sample groups")
}
if (anyDuplicated(key(pairing)) || anyDuplicated(key(audit))) {
  stop("Pairing and material-audit inputs must have unique PXD + sample-group keys")
}
if (!setequal(key(statistics), key(audit))) {
  stop("The 33-row statistics and strict material audit do not have the same scope")
}

pairing_33 <- pairing[match(key(statistics), key(pairing)), , drop = FALSE]
audit_33 <- audit[match(key(statistics), key(audit)), , drop = FALSE]
if (any(is.na(pairing_33$乳酸化PXD)) || any(is.na(audit_33$乳酸化PXD))) {
  stop("Not every final sample group has pairing and audit metadata")
}

translate_state <- function(value) {
  labels <- c(
    baseline_reference_only =
      "同一材料的基线参照；并非Kla实验的同处理或同时间点",
    independent_cohort_same_phenotype =
      "独立队列的同一组织/疾病表型；不是同一供体",
    same_study_matched_samples = "同一研究中的匹配样本",
    not_available = "无可用的严格材料匹配定量参照",
    same_study_matched_conditions = "同一研究中的匹配实验条件",
    exact_same_models = "相同细胞或实验模型",
    exact_same_samples = "同一批实际样本",
    same_cell_line_independent_experiment = "同一细胞系，但属于独立实验",
    same_study_matched_material = "同一研究中的匹配生物材料",
    same_study_states_present_method_pending =
      "同一研究包含对应状态，但普通全蛋白定量方法仍需老师确认"
  )
  translated <- unname(labels[value])
  translated[is.na(translated)] <- value[is.na(translated)]
  translated
}

translate_decision <- function(value) {
  ifelse(
    value == "included_exact_material_identity",
    "是",
    ifelse(value == "excluded_no_exact_material_reference", "否", value)
  )
}

translate_granularity <- function(value) {
  labels <- c(
    adjacent_tissue = "癌旁组织",
    adjacent_tissue_sample_group = "癌旁组织样本组",
    biospecimen = "具体生物样本",
    cell_line = "细胞系",
    cell_line_baseline = "细胞系基线",
    cell_line_with_conditions = "带实验条件的细胞系",
    cell_line_with_genotypes = "带基因型的细胞系",
    cell_line_with_treatment = "带处理的细胞系",
    cell_model_with_treatment = "带处理的细胞模型",
    disease_tissue = "疾病组织",
    disease_tissue_sample_group = "疾病组织样本组",
    disease_tissue_subgroups = "疾病组织亚组",
    named_cell_models = "明确命名的细胞模型",
    primary_cell_culture_baseline = "原代细胞培养基线",
    primary_cell_culture_with_conditions = "带条件的原代细胞培养",
    same_biospecimens = "同一批实际生物样本",
    tissue_sample_group = "组织样本组",
    tissue_unspecified_subregion = "组织，未细分解剖亚区",
    whole_organ_unspecified_subregion = "完整器官，未细分解剖亚区",
    `未纳入` = "未纳入"
  )
  translated <- unname(labels[value])
  translated[is.na(translated)] <- value[is.na(translated)]
  translated
}

translate_reference_strategy <- function(value) {
  labels <- c(
    external_exact_adjacent_tissue =
      "外部独立队列的精确癌旁组织普通全蛋白",
    external_exact_disease_tissue =
      "外部独立队列的精确疾病组织普通全蛋白",
    same_biospecimen_non_ptm_proteome =
      "同一生物样本的非PTM普通全蛋白",
    same_study_non_enrichment =
      "同研究非富集普通全蛋白",
    same_study_non_ptm_baseline =
      "同研究非PTM基线普通全蛋白",
    same_study_same_biospecimen_non_ptm_proteome =
      "同研究同一生物样本的非PTM普通全蛋白"
  )
  translated <- unname(labels[value])
  translated[is.na(translated)] <- value[is.na(translated)]
  translated
}

clean_text <- function(value, fallback) {
  value <- trimws(ifelse(is.na(value), "", value))
  value[value == ""] <- fallback
  value
}

included <- statistics$纳入严格配对分析 %in%
  c(TRUE, "TRUE", "True", 1, "1")
material_match <- audit_33$材料身份是否严格匹配 %in%
  c(TRUE, "TRUE", "True", 1, "1")
if (!identical(included, material_match)) {
  stop("Statistics inclusion disagrees with strict material-identity audit")
}

status_note <- translate_state(audit_33$实验状态匹配情况)
original_note <- clean_text(audit_33$原注意事项, "")
attention <- ifelse(
  included,
  paste0(status_note, ifelse(original_note == "", "", paste0("；", original_note))),
  clean_text(audit_33$决定依据, original_note)
)
attention <- gsub("；；+", "；", attention)

review <- data.frame(
  序号 = seq_len(nrow(statistics)),
  四分类 = statistics$四分类,
  乳酸化PXD = statistics$乳酸化PXD,
  样本组 = statistics$样本组,
  Kla材料 = audit_33$Kla材料,
  Kla取材粒度 = translate_granularity(audit_33$Kla取材粒度),
  Kla实际使用文件 = statistics$乳酸化证据文件,
  Kla证据形式 = pairing_33$乳酸化数据判定,
  Kla获取状态 = pairing_33$乳酸化获取状态,
  Kla蛋白数 = as.integer(statistics$乳酸化蛋白ID数),
  Kla_DDR蛋白数 = as.integer(statistics$乳酸化DDR蛋白ID数),
  Kla_DDR占比 = as.numeric(statistics$乳酸化DDR占比),
  普通全蛋白PXD = ifelse(
    included,
    statistics$常规蛋白组PXD,
    "未纳入：无完全匹配且可定量参照"
  ),
  参照材料 = ifelse(
    included,
    audit_33$参照材料,
    "未纳入：无完全匹配且可定量参照"
  ),
  参照取材粒度 = ifelse(
    included,
    translate_granularity(audit_33$参照取材粒度),
    "未纳入"
  ),
  普通全蛋白实际使用文件 = ifelse(
    included,
    audit_33$普通全蛋白文件,
    "未纳入"
  ),
  实际读取样本_工作表_列 = ifelse(
    included,
    audit_33$实际读取子集,
    "未纳入"
  ),
  普通全蛋白证据形式 = ifelse(
    included,
    translate_reference_strategy(pairing_33$常规蛋白组策略),
    "未纳入"
  ),
  普通全蛋白获取状态 = ifelse(
    included,
    pairing_33$常规蛋白组获取状态,
    "未纳入"
  ),
  普通全蛋白蛋白数 = as.integer(statistics$常规蛋白ID数),
  普通全蛋白DDR数 = as.integer(statistics$常规DDR蛋白ID数),
  普通全蛋白DDR占比 = as.numeric(statistics$常规DDR占比),
  材料身份严格匹配 = ifelse(material_match, "是", "否"),
  实验状态匹配情况 = status_note,
  是否纳入配对分析 = translate_decision(audit_33$分析决定),
  参照关系说明 = clean_text(statistics$参照匹配说明, "未记录"),
  需要老师注意的问题 = clean_text(attention, "无额外重大问题"),
  论文DOI = clean_text(
    sub("\\.$", "", trimws(pairing_33$论文DOI)),
    "元数据未记录或尚未核实"
  ),
  ProteomeXchange链接 = pairing_33$数据集链接,
  check.names = FALSE
)

category_levels <- c(
  "正常/非肿瘤组织",
  "癌症组织",
  "正常/非肿瘤细胞",
  "癌症细胞"
)
category_counts <- table(factor(review$四分类, levels = category_levels))
if (!identical(as.integer(category_counts), c(9L, 2L, 9L, 13L))) {
  stop("Unexpected four-class group counts: ", paste(category_counts, collapse = ", "))
}
if (sum(review$是否纳入配对分析 == "是") != 33 ||
    sum(review$是否纳入配对分析 == "否") != 0) {
  stop("Teacher review table must contain exactly 33 included rows")
}

write.csv(review, output_csv, row.names = FALSE, na = "")

focus <- review |>
  filter(
    是否纳入配对分析 == "否" |
      grepl(
        "基线参照|独立队列|独立实验|仍需老师确认",
        实验状态匹配情况
      )
  ) |>
  select(
    序号, 四分类, 乳酸化PXD, 样本组, Kla材料,
    普通全蛋白PXD, 参照材料, 实际读取样本_工作表_列,
    实验状态匹配情况, 是否纳入配对分析, 需要老师注意的问题
  )
write.csv(focus, focus_csv, row.names = FALSE, na = "")

message(
  "Wrote teacher review table: ", nrow(review), " rows; ",
  sum(review$是否纳入配对分析 == "是"), " included; ",
  sum(review$是否纳入配对分析 == "否"), " excluded."
)
