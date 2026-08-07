#!/usr/bin/env Rscript
# encoding: utf-8
# build_workbooks.R — 合并 7 个 Node/Python 交付工具的全 R 重写（Task 10）
# 参照（只读，未修改）：
#   build_cell_type_statistics_workbook.mjs
#   build_hippocampus_review_md.mjs
#   build_reference_ddr_comparison_workbook.mjs
#   build_reference_proteome_selection_workbook.mjs
#   build_venn_combined_workbook.mjs
#   create_bilingual_figure_legends_docx.js
#   build_project_metadata.py
# 输出规格见 docs/superpowers/plans/task10_output_spec.md
#
# 用法: Rscript build_workbooks.R <project_root> [--stage workbooks|metadata]
#   --stage workbooks : xlsx（3 个）+ md（1 个）+ 附属 csv（4 个）+ docx（1 个）
#   --stage metadata  : per-PXD metadata csv + project_file_inventory.csv + reconstructed manifest
#   默认（无 --stage）: 两者都运行（先 workbooks 后 metadata）

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("usage: Rscript build_workbooks.R <project_root> [--stage workbooks|metadata]")
}
project_root <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
stage <- "all"
if (length(args) >= 3 && args[[2]] == "--stage") {
  stage <- args[[3]]
}
if (!stage %in% c("all", "workbooks", "metadata")) {
  stop("--stage must be one of: workbooks | metadata")
}

suppressPackageStartupMessages({
  library(writexl)
  library(officer)
  library(jsonlite)
  library(data.table)
  library(Rmpfr)
})

p <- function(...) file.path(project_root, ...)

# ============================================================
# 通用工具
# ============================================================

write_utf8 <- function(path, text) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  con <- file(path, open = "wb")
  on.exit(close(con))
  writeBin(charToRaw(enc2utf8(text)), con)
}

# python csv.QUOTE_MINIMAL 字段转义
csv_field <- function(v) {
  s <- if (is.null(v) || is.na(v)) "" else as.character(v)
  if (grepl("[,\"\r\n]", s)) paste0('"', gsub('"', '""', s, fixed = TRUE), '"') else s
}

# .mjs 风格 csv：BOM + LF + QUOTE_MINIMAL + 无结尾换行
write_mjs_csv <- function(path, header, rows, bom = TRUE) {
  all_rows <- c(list(header), rows)
  lines <- vapply(all_rows, function(r) {
    paste(vapply(r, csv_field, character(1)), collapse = ",")
  }, character(1))
  text <- paste0(lines, collapse = "\n")
  if (bom) text <- paste0("\uFEFF", text)
  write_utf8(path, text)
}

# .py 风格 csv：无 BOM + LF + QUOTE_MINIMAL + 每行（含末行）以 \n 结尾
write_py_csv <- function(path, fields, rows) {
  lines <- c(
    paste(vapply(fields, csv_field, character(1)), collapse = ","),
    vapply(rows, function(r) {
      paste(vapply(fields, function(f) csv_field(r[[f]]), character(1)), collapse = ",")
    }, character(1))
  )
  write_utf8(path, paste0(paste0(lines, collapse = "\n"), "\n"))
}

# 读 CSV（utf-8-sig），全部按字符读，保留原始字符串（不解析数字、不转 NA）
read_csv_chr <- function(path) {
  utils::read.csv(path, fileEncoding = "UTF-8-BOM", colClasses = "character",
                  check.names = FALSE, stringsAsFactors = FALSE,
                  na.strings = "!!!__NO_SUCH_NA_SENTINEL__!!!")
}

# JS String(value ?? "未报告") + |→\|、\r\n→<br>、\n→<br>
# （注意：JS 的 ?? 只对 null/undefined 生效，空字符串保持原样）
js_cell <- function(v) {
  if (is.null(v) || length(v) == 0 || is.na(v)) v <- "未报告"
  s <- as.character(v)
  s <- gsub("|", "\\|", s, fixed = TRUE)
  s <- gsub("\r\n", "<br>", s, fixed = TRUE)
  s <- gsub("\n", "<br>", s, fixed = TRUE)
  s
}

# JS value.split(";") → trim → filter(Boolean)
js_split_trim <- function(v) {
  if (is.null(v) || length(v) == 0 || is.na(v)) return(character(0))
  parts <- strsplit(as.character(v), ";", fixed = TRUE)[[1]]
  parts <- trimws(parts)
  parts[parts != ""]
}

# Python Path.suffix 语义（小写由调用方处理）
py_suffix <- function(name) {
  n <- nchar(name)
  dots <- gregexpr(".", name, fixed = TRUE)[[1]]
  if (length(dots) == 1 && dots[1] == -1L) return("")
  i <- dots[[length(dots)]]
  if (i > 1 && i < n) substr(name, i, n) else ""
}

# JS Number→String 最短回读格式（ECMAScript Number::toString 语义）。
# 用于 .mjs 版 xlsx 中公式缓存值的文本一致：writexl 对 double 只保留 16 位小数序列化，
# 会把 0.036806883365200764 写成 0.03680688336520076；而 .mjs 引擎存的是最短回读串。
# 因此计算列以字符串形式写入，逐字符复现 .mjs 缓存值。
# 回读判定：十进制串先经 mpfr 精确解析（128 位），再按最近偶数舍入回 double
# （mpfr 的 double 转换是正确舍入的；R 自带 strtod 在个别边界值有 1 ulp 偏差，
# 已实测 as.numeric 对 "0.9555274736986343"/"0.03680688336520077" 判读与 glibc/JS 不同）。
str_roundtrips <- function(s, x) {
  as.numeric(Rmpfr::mpfr(s, precBits = 128)) == x
}

js_num_str <- function(x) {
  if (is.nan(x)) return("NaN")
  if (is.infinite(x)) return(if (x > 0) "Infinity" else "-Infinity")
  if (x == 0) return("0")
  if (x < 0) return(paste0("-", js_num_str(-x)))
  if (x >= 1e21) {
    e <- floor(log10(x))
    m <- x / 10^e
    ms <- ""
    for (d in 1:17) {
      cand <- sprintf(paste0("%.", d, "g"), m)
      if (str_roundtrips(cand, m)) { ms <- cand; break }
    }
    if (ms == "") ms <- sprintf("%.17g", m)
    return(paste0(ms, "e+", sprintf("%02d", e)))
  }
  if (x == floor(x)) return(sprintf("%.0f", x))
  for (d in 1:17) {
    s <- sprintf(paste0("%.", d, "g"), x)
    if (str_roundtrips(s, x)) return(s)
  }
  sprintf("%.17g", x)
}

# 公式计算值 → JS 等价文本（向量化）
js_fraction_str <- function(num, den) {
  vapply(seq_along(num), function(i) {
    if (is.na(den[[i]]) || den[[i]] == 0) "0" else js_num_str(num[[i]] / den[[i]])
  }, "")
}

# ============================================================
# 工具 1: build_cell_type_statistics_workbook.mjs
#   → cell_type_kla_ddr_statistics.xlsx
# ============================================================
build_cell_type_statistics_xlsx <- function() {
  stats <- read_csv_chr(p("reanalysis/results/tables/cell_type_kla_ddr_statistics.csv"))
  total <- as.numeric(stats$TotalKlaProteins)
  ddr <- as.numeric(stats$KlaGoDdrProteins)
  # D 列 = .mjs 公式 =IFERROR(C/B,0) 的缓存值（JS 最短回读文本）
  df <- data.frame(
    `细胞系或组织` = stats$CellOrTissueType,
    `总 Kla 蛋白数` = total,
    `与 DDR 交集蛋白数` = ddr,
    `DDR 交集 / 总数` = js_fraction_str(ddr, total),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  out <- p("reanalysis/results/tables/cell_type_kla_ddr_statistics.xlsx")
  writexl::write_xlsx(list(`按细胞类型统计` = df), out)
  cat("workbook:", basename(out), "rows:", nrow(df), "\n")
  out
}

# ============================================================
# 工具 2: build_hippocampus_review_md.mjs
#   → HUMAN_HIPPOCAMPUS_25_DATASET_REVIEW_TABLE.md
# ============================================================
build_hippocampus_review_md <- function() {
  source_path <- p("reanalysis/config/human_hippocampus_25_datasets.json")
  output_path <- p("reanalysis/reports/HUMAN_HIPPOCAMPUS_25_DATASET_REVIEW_TABLE.md")
  rows <- fromJSON(source_path, simplifyVector = FALSE)
  pxds <- vapply(rows, function(r) r$pxd, "")
  if (length(rows) != 25 || length(unique(pxds)) != 25) {
    stop("Expected exactly 25 unique PXD records.")
  }

  pxd_link <- function(pxd) sprintf("[%s](https://proteomecentral.proteomexchange.org/?pxid=%s)", pxd, pxd)

  doi_links <- function(v) {
    if (is.null(v) || is.na(v) || v == "" || v == "未报告") return("未报告")
    parts <- js_split_trim(v)
    paste0(vapply(parts, function(doi) sprintf("[%s](https://doi.org/%s)", doi, doi), ""), collapse = "<br>")
  }

  pmid_links <- function(v) {
    if (is.null(v) || is.na(v) || v == "" || v == "未报告") return("")
    parts <- js_split_trim(v)
    paste0(vapply(parts, function(pmid) sprintf("[PMID %s](https://pubmed.ncbi.nlm.nih.gov/%s/)", pmid, pmid), ""), collapse = "<br>")
  }

  headers <- c(
    "PXD", "数据集标题", "年份/发表", "真实材料与模型", "是否原生人海马", "人供者数",
    "性别", "年龄", "疾病/对照", "研究者额外操作", "样本处理与MS方法", "论文",
    "提交者/PI", "重复或关联关系", "Kla适用性", "等级与纳入建议", "证据来源/备注"
  )

  lines <- c(
    "# ProteomeCentral 人海马 25 份数据审核大表",
    "",
    "生成日期：2026-07-30",
    "",
    "检索条件：[ProteomeCentral：hippocampus + Homo sapiens](https://proteomecentral.proteomexchange.org/ui?view=datasets&search=hippocampus&species=Homo%20sapiens)",
    "",
    "> 审核提醒：数据库物种字段不能直接等同于真实样本。当前 25 条中，14 条含原生人海马，1 条为含人尸检海马的混合模型，10 条并非原生人海马。普通蛋白组鉴定也不能直接作为 Kla 位点证据。",
    "",
    paste0("| ", paste(headers, collapse = " | "), " |"),
    paste0("| ", paste(rep("---", length(headers)), collapse = " | "), " |")
  )

  row_list <- list()
  for (row in rows) {
    article <- c(
      doi_links(row$doi),
      pmid_links(row$pmid),
      if (!is.null(row$first_author) && !is.na(row$first_author) &&
          row$first_author != "" && row$first_author != "未报告") {
        paste0("第一作者：", js_cell(row$first_author))
      } else ""
    )
    article <- article[article != ""]

    model <- paste0(c(js_cell(row$actual_material), paste0("模型：", js_cell(row$model_type))), collapse = "<br>")
    methods <- paste0(c(paste0("处理：", js_cell(row$sample_processing)),
                        paste0("MS：", js_cell(row$ms_method))), collapse = "<br>")
    grade <- paste0(c(paste0("等级：", js_cell(row$recommendation_grade)),
                      js_cell(row$inclusion_advice)), collapse = "<br>")
    evidence <- paste0(c(js_cell(row$evidence_level),
                         paste0("人口学来源：", js_cell(row$demographic_source)),
                         js_cell(row$notes)), collapse = "<br>")

    dup <- if (is.null(row$duplicate_related) || is.na(row$duplicate_related) ||
               row$duplicate_related == "") "无已知重复登记" else js_cell(row$duplicate_related)

    values <- c(
      pxd_link(row$pxd),
      js_cell(row$title),
      paste0(js_cell(row$paper_year), "<br>", js_cell(row$publication_status),
             "<br>数据库发布：", js_cell(row$announce_date)),
      model,
      js_cell(row$native_human_hippocampus),
      js_cell(row$human_donor_n),
      js_cell(row$sex),
      js_cell(row$age),
      js_cell(row$disease_group),
      js_cell(row$extra_operation),
      methods,
      paste0(article, collapse = "<br>"),
      js_cell(row$submitter_pi),
      dup,
      js_cell(row$kla_suitability),
      grade,
      evidence
    )
    lines <- c(lines, paste0("| ", paste(values, collapse = " | "), " |"))
  }

  lines <- c(
    lines,
    "",
    "## 等级说明",
    "",
    "- A：正常或基线原生人海马。",
    "- B：原生人海马，但包含疾病、用药、保存方法或其他处理。",
    "- C：混合模型，包含人尸检海马但还混有体外或动物模型。",
    "- D：人源体外模型，如 iPSC、类器官或海马球体。",
    "- E：小鼠海马、其他人细胞或关键词误命中。",
    "- R：旧数据重分析或与其他 PXD 明确重叠。",
    "",
    "## 重点复核项",
    "",
    "- `PXD010543` 与 `PXD010544` 是同一批 4 TLE + 4 control 生物样本的两个 TMT 技术批次。",
    "- `PXD000395` 与 `PXD000950` 存在 control PE 原始文件重叠。",
    "- `PXD062981` 有相关论文，但论文 Data availability 未明确列出该 PXD，关联仍需复核。",
    "- `PXD050470` 是本批数据中唯一直接以正常人海马 Kla 图谱为主要目标的数据。",
    ""
  )

  write_utf8(output_path, paste0(lines, collapse = "\n"))
  cat("report:", basename(output_path), "rows:", length(rows), "\n")
  output_path
}

# ============================================================
# 工具 3: build_reference_ddr_comparison_workbook.mjs
#   → cell_type_kla_vs_reference_ddr_statistics.xlsx
#     + cell_type_kla_vs_reference_ddr_statistics_zh.csv
#     + cell_type_reference_control_information_zh.csv
# ============================================================

# chineseInfo（逐字照抄 .mjs）
chinese_info <- list(
  "Human hippocampus" = list(
    name = "人海马组织", matchType = "组织精确匹配",
    subset = "74名神经学正常的CA1海马供者",
    baseline = "神经学正常的死后CA1海马组织；未进行PTM富集",
    sampleCount = "74名供者，74个raw文件",
    acquisition = "非标记LC-MS/MS；Q Exactive HF",
    search = "Spectronaut Pulsar X；1% FDR；至少3条肽段用于定量",
    ptm = "否",
    rationale = "人CA1海马组织精确匹配；74名神经学正常供者；年龄范围接近Kla海马样本；无PTM富集且有论文发布的处理矩阵。",
    caveat = "供者年龄为66-104岁，材料为CA1；适合作为老年海马背景，但不能代表年轻海马或整个海马组织。",
    countDetail = "74名供者共2,092个论文发布的基因/蛋白特征；单供者2,018-2,092个，中位数2,090个。",
    completeness = "ProteomeXchange为PARTIAL，但74个raw文件和论文处理矩阵均可获得。",
    statistics = "可在控制年龄、性别和死后间隔后进行供者层面统计。"
  ),
  HEK293T = list(
    name = "HEK293T", matchType = "细胞系精确匹配",
    subset = "Control_HEK293T_lys；主分析排除Control_HEK293T_std_H002标准QC",
    baseline = "未处理HEK293T裂解物过程控制；未进行PTM富集",
    sampleCount = "401个裂解物运行；另有663个标准QC运行",
    acquisition = "DIA/SWATH；TripleTOF 6600",
    search = "DIA-NN检索，随后使用maxLFQ定量",
    ptm = "否",
    rationale = "细胞系精确匹配；未处理裂解物；标准化DIA平台且运行数多；标准QC材料被明确排除。",
    caveat = "这些样本是过程控制而非独立生物学重复，只能作为蛋白检出背景，不能把401个运行当作401个生物重复。",
    countDetail = "401个Control_HEK293T_lys运行中共检出6,441个高置信BaseAccession；单运行1,987-4,759个，中位数3,512个。",
    completeness = "PRIDE为PARTIAL；Figshare提供最终矩阵、映射和肽段计数文件。",
    statistics = "不可将过程控制运行作为独立生物重复进行差异检验。"
  ),
  "HK-2" = list(
    name = "HK-2", matchType = "细胞系精确匹配",
    subset = "HK-2_Control_1、HK-2_Control_3和HK-2_Control_4",
    baseline = "KSFM加2% FBS培养24小时的未处理对照；未进行PTM富集",
    sampleCount = "3个未处理对照；项目共9个raw文件",
    acquisition = "DirectDIA非标记LC-MS/MS；Orbitrap Exploris 240",
    search = "Spectronaut 19 directDIA；跨运行归一化；未插补",
    ptm = "否",
    rationale = "细胞系精确匹配；有3个明确的未处理对照；蛋白组完整度高；提供Spectronaut蛋白组报告。",
    caveat = "只纳入Control_1、Control_3和Control_4，不能混入Cu或TTM处理组；保留仓库原始编号。",
    countDetail = "3个对照的4,791个定量蛋白组展开后得到4,897个唯一BaseAccession；单对照蛋白组数为4,514-4,670。",
    completeness = "ProteomeXchange为PARTIAL，但9个raw文件和关键Spectronaut报告可获得。",
    statistics = "可对3个未处理对照进行组内描述，但不能混入Cu或TTM样本。"
  ),
  MCF10A = list(
    name = "MCF10A", matchType = "细胞系精确匹配",
    subset = "未处理MCF-10A基线细胞裂解物",
    baseline = "汇合状态的未处理MCF-10A细胞；未进行PTM富集",
    sampleCount = "10个肽段IEF组分×2次进样，共20个raw文件",
    acquisition = "DDA；LTQ Orbitrap Velos；240分钟梯度",
    search = "MaxQuant 1.4.1.2 / Andromeda",
    ptm = "否；肽段等电聚焦属于分析分级，不是PTM富集",
    rationale = "未处理MCF-10A精确匹配；分级深度较高；仓库提供MaxQuant evidence，适合作为检出背景。",
    caveat = "平台较旧且进行了深度分级；组分和重复进样属于技术深度，不能作为生物重复，也不宜与现代DIA绝对强度直接比较。",
    countDetail = "20个MCF-10A组分/进样的evidence中共有4,839个有效leading-razor BaseAccession；仓库未提供proteinGroups.txt。",
    completeness = "ProteomeXchange为PARTIAL；包含大量raw文件和一个合并MaxQuant检索结果包。",
    statistics = "不可将IEF组分和重复进样作为独立生物重复。"
  ),
  MCF7 = list(
    name = "MCF7", matchType = "细胞系精确匹配",
    subset = "MCF7 / SIDM00148",
    baseline = "未处理基线培养；未进行PTM富集",
    sampleCount = "6次DIA采集",
    acquisition = "DIA/SWATH；TripleTOF 6600",
    search = "DIA-NN检索，随后使用maxLFQ定量",
    ptm = "否",
    rationale = "未处理MCF7精确匹配；来自统一的泛癌DIA图谱；具有高置信蛋白矩阵和明确的细胞系映射。",
    caveat = "Kla MCF7数据包含DCA、缺氧、鱼藤酮和同位素示踪条件；该对照只能作为未处理检出背景，不是配对实验对照。",
    countDetail = "MCF7/SIDM00148平均细胞系行中有3,099个非缺失高置信BaseAccession。",
    completeness = "PRIDE为PARTIAL；Figshare提供完整分析矩阵。",
    statistics = "6次采集主要反映技术采集，不能直接作为独立生物重复。"
  ),
  HCT116 = list(
    name = "HCT116", matchType = "细胞系精确匹配",
    subset = "HCT-116 / SIDM00783",
    baseline = "未处理基线培养；未进行PTM富集",
    sampleCount = "6次DIA采集",
    acquisition = "DIA/SWATH；TripleTOF 6600",
    search = "DIA-NN检索，随后使用maxLFQ定量",
    ptm = "否",
    rationale = "未处理HCT-116精确匹配；统一DIA平台；可作为独立于Kla亚细胞分级实验的全细胞检出背景。",
    caveat = "PXD053474含亚细胞组分；全细胞参考不能用于推断某蛋白在特定亚细胞组分中缺失。",
    countDetail = "HCT-116/SIDM00783平均细胞系行中有4,010个非缺失高置信BaseAccession。",
    completeness = "PRIDE为PARTIAL；Figshare提供完整分析矩阵。",
    statistics = "6次采集主要反映技术采集，不能直接作为独立生物重复。"
  ),
  "T-ALL" = list(
    name = "T-ALL（TALL-1替代）", matchType = "疾病类别替代，非精确细胞系匹配",
    subset = "主替代：TALL-1 / SIDM00370；敏感性替代：Jurkat / SIDM01016",
    baseline = "未处理T淋巴母细胞白血病细胞系基线；未进行PTM富集",
    sampleCount = "TALL-1为6次DIA采集；Jurkat为6次DIA采集",
    acquisition = "DIA/SWATH；TripleTOF 6600",
    search = "DIA-NN检索，随后使用maxLFQ定量",
    ptm = "否",
    rationale = "未找到独立、未处理的TALL-104全蛋白组；因此TALL-1作为主疾病类别替代，Jurkat单独用于敏感性分析。",
    caveat = "TALL-1和Jurkat均不是TALL-104。PXD028488的TALL-104非富集组经历乳酸暴露，也不是独立正常基线。",
    countDetail = "TALL-1主分析3,383个蛋白；Jurkat 3,363个；并集3,835个；交集2,911个。",
    completeness = "PRIDE为PARTIAL；Figshare提供完整分析矩阵。",
    statistics = "只能作为替代背景，不可写成TALL-104精确对照。"
  ),
  "MDA-MB-468" = list(
    name = "MDA-MB-468", matchType = "细胞系精确匹配",
    subset = "MDA-MB-468 / SIDM00628",
    baseline = "未处理基线培养；未进行PTM富集",
    sampleCount = "6次DIA采集",
    acquisition = "DIA/SWATH；TripleTOF 6600",
    search = "DIA-NN检索，随后使用maxLFQ定量",
    ptm = "否",
    rationale = "未处理MDA-MB-468精确匹配；来自统一的泛癌DIA图谱和高置信蛋白矩阵。",
    caveat = "只用于同细胞系的蛋白检出背景，不与MaxQuant Kla富集数据直接比较绝对强度。",
    countDetail = "MDA-MB-468/SIDM00628平均细胞系行中有3,760个非缺失高置信BaseAccession。",
    completeness = "PRIDE为PARTIAL；Figshare提供完整分析矩阵。",
    statistics = "6次采集主要反映技术采集，不能直接作为独立生物重复。"
  ),
  "T-47D" = list(
    name = "T-47D", matchType = "细胞系精确匹配",
    subset = "T47D / SIDM00097",
    baseline = "未处理基线培养；未进行PTM富集",
    sampleCount = "6次DIA采集",
    acquisition = "DIA/SWATH；TripleTOF 6600",
    search = "DIA-NN检索，随后使用maxLFQ定量",
    ptm = "否",
    rationale = "未处理T47D精确匹配；来自统一泛癌DIA图谱；已将来源标签T47D明确映射到项目标签T-47D。",
    caveat = "来源中的名称为T47D而非T-47D；只作检出背景。",
    countDetail = "T47D/SIDM00097平均细胞系行中有3,516个非缺失高置信BaseAccession。",
    completeness = "PRIDE为PARTIAL；Figshare提供完整分析矩阵。",
    statistics = "6次采集主要反映技术采集，不能直接作为独立生物重复。"
  ),
  RKO = list(
    name = "RKO", matchType = "细胞系精确匹配",
    subset = "RKO / SIDM01090",
    baseline = "未处理基线培养；未进行PTM富集",
    sampleCount = "6次DIA采集",
    acquisition = "DIA/SWATH；TripleTOF 6600",
    search = "DIA-NN检索，随后使用maxLFQ定量",
    ptm = "否",
    rationale = "未处理RKO精确匹配；来自统一泛癌DIA图谱；作为Kla WT/KO实验的共同检出背景。",
    caveat = "外部RKO蛋白组不能替代GSK3B WT与KO的组内比较。",
    countDetail = "RKO/SIDM01090平均细胞系行中有3,615个非缺失高置信BaseAccession。",
    completeness = "PRIDE为PARTIAL；Figshare提供完整分析矩阵。",
    statistics = "6次采集主要反映技术采集，不能直接作为独立生物重复。"
  )
)

methods_rows <- list(
  c("项目", "说明"),
  c("DDR定义", "使用 data/annotations/GO-repair+damage(human).tsv；仅保留人源且排除 qualifier 含 NOT 的注释。"),
  c("匹配顺序", "UniProt BaseAccession 优先，来源数据具有 GeneSymbol 时再作辅助匹配。"),
  c("主分母", "每个参考对照实际纳入的唯一 BaseAccession；海马体为论文发布的 2,092 个 gene/protein feature。"),
  c("PXD030304", "使用 6,692 蛋白高置信矩阵。HEK293T 由 401 个 Control_HEK293T_lys 运行的 peptide-count 并集计算，标准QC排除。"),
  c("T-ALL", "10类主表使用 TALL-1 作为 TALL-104 的疾病类别替代；Jurkat、并集和交集结果单独保留。"),
  c("统计解释", "本结果是描述性蛋白集合比例，不把技术进样当作生物重复，不进行跨平台丰度显著性检验。"),
  c("未检出解释", "常规蛋白组或 Kla 数据未检出均不能直接解释为蛋白不存在或不发生乳酸化。")
)

build_reference_ddr_comparison <- function() {
  tables_dir <- p("reanalysis/results/tables")
  comparison <- read_csv_chr(file.path(tables_dir, "cell_type_kla_vs_reference_ddr_statistics.csv"))
  tall <- read_csv_chr(file.path(tables_dir, "tall104_surrogate_ddr_sensitivity.csv"))
  control_info <- read_csv_chr(file.path(tables_dir, "cell_type_reference_control_information.csv"))
  if (nrow(comparison) != 10) stop("Expected 10 comparison rows, found ", nrow(comparison))
  if (nrow(control_info) != 10) stop("Expected 10 control-information rows, found ", nrow(control_info))

  zh <- lapply(comparison$CellOrTissueType, function(ct) {
    ci <- chinese_info[[ct]]
    if (is.null(ci)) stop("Missing chineseInfo for ", ct)
    ci
  })
  names(zh) <- comparison$CellOrTissueType

  # ---- sheet 1: DDR占比对照（16 列）----
  ref_d <- as.numeric(comparison$ReferenceProteinCount)
  ref_e <- as.numeric(comparison$ReferenceDdrProteinCount)
  kla_g <- as.numeric(comparison$TotalKlaProteins)
  kla_h <- as.numeric(comparison$KlaGoDdrProteins)
  # F/I/J/K 列 = .mjs 公式缓存值（JS 最短回读文本），J、P 与 .mjs 公式同式
  f_txt <- js_fraction_str(ref_e, ref_d)
  i_txt <- js_fraction_str(kla_h, kla_g)
  f_dbl <- ifelse(ref_d == 0, 0, ref_e / ref_d)
  i_dbl <- ifelse(kla_g == 0, 0, kla_h / kla_g)
  j_dbl <- i_dbl - f_dbl
  j_txt <- vapply(j_dbl, js_num_str, "")
  k_txt <- vapply(ifelse(f_dbl == 0, 0, i_dbl / f_dbl), js_num_str, "")
  p_desc <- ifelse(abs(j_dbl) < 0.002, "基本相当",
                   ifelse(j_dbl > 0, "Kla中DDR占比更高", "参考组中DDR占比更高"))
  unit_label <- ifelse(comparison$ReferenceCountUnit == "BaseAccession",
                       "UniProt基础登录号", "论文发布的基因/蛋白特征")

  df_comparison <- data.frame(
    `细胞系或组织` = vapply(zh, `[[`, "", "name"),
    `参考PXD` = comparison$ReferencePXD,
    `参考样本子集` = vapply(zh, `[[`, "", "subset"),
    `参考蛋白数` = ref_d,
    `参考DDR蛋白数` = ref_e,
    `参考DDR占比` = f_txt,
    `Kla蛋白数` = kla_g,
    `Kla∩DDR蛋白数` = kla_h,
    `Kla DDR占比` = i_txt,
    `Kla-参考百分点差` = j_txt,
    `Kla/参考占比倍数` = k_txt,
    `参考计数单位` = unit_label,
    `登录号匹配数` = as.numeric(comparison$BaseAccessionMatches),
    `基因符号辅助匹配数` = as.numeric(comparison$GeneSymbolFallbackMatches),
    `蛋白数计数口径` = vapply(zh, `[[`, "", "countDetail"),
    `描述性结论` = p_desc,
    check.names = FALSE, stringsAsFactors = FALSE
  )

  # ---- sheet 2: 对照选择信息（23 列）----
  ctl <- lapply(control_info[["细胞或组织"]], function(ct) {
    ci <- chinese_info[[ct]]
    if (is.null(ci)) stop("Missing chineseInfo for ", ct)
    ci
  })
  names(ctl) <- control_info[["细胞或组织"]]

  df_control <- data.frame(
    `细胞或组织` = vapply(ctl, `[[`, "", "name"),
    `参考PXD` = control_info[["参考PXD"]],
    `年份` = as.numeric(control_info[["年份"]]),
    `匹配类型` = vapply(ctl, `[[`, "", "matchType"),
    `适用等级` = control_info[["适用等级"]],
    `选用样本或子集` = vapply(ctl, `[[`, "", "subset"),
    `正常或基线条件` = vapply(ctl, `[[`, "", "baseline"),
    `样本或采集数` = vapply(ctl, `[[`, "", "sampleCount"),
    `对照蛋白数` = as.numeric(control_info[["对照蛋白数"]]),
    `DDR蛋白数` = as.numeric(control_info[["DDR蛋白数"]]),
    # DDR占比 为小数：writexl 数值序列化会丢 1 ulp，用原串（JS 规范形）以文本写入
    `DDR占比` = control_info[["DDR占比"]],
    `采集方式或仪器` = vapply(ctl, `[[`, "", "acquisition"),
    `检索与定量` = vapply(ctl, `[[`, "", "search"),
    `是否PTM富集` = vapply(ctl, `[[`, "", "ptm"),
    `选择理由` = vapply(ctl, `[[`, "", "rationale"),
    `主要限制` = vapply(ctl, `[[`, "", "caveat"),
    `蛋白数计数口径` = vapply(ctl, `[[`, "", "countDetail"),
    `实际分析来源文件` = control_info[["实际分析来源文件"]],
    `仓库与分析完整度` = vapply(ctl, `[[`, "", "completeness"),
    `可作统计差异` = vapply(ctl, `[[`, "", "statistics"),
    `数据集URL` = control_info[["数据集URL"]],
    `论文URL` = control_info[["论文URL"]],
    `处理数据URL` = control_info[["处理数据URL"]],
    check.names = FALSE, stringsAsFactors = FALSE
  )

  # ---- sheet 3: TALL替代敏感性（4 列）----
  tall_b <- as.numeric(tall$ProteinCount)
  tall_c <- as.numeric(tall$DdrProteinCount)
  df_tall <- data.frame(
    `参考集合` = tall$TAllReferenceSet,
    `蛋白数` = tall_b,
    `DDR蛋白数` = tall_c,
    `DDR占比` = js_fraction_str(tall_c, tall_b),
    check.names = FALSE, stringsAsFactors = FALSE
  )

  # ---- sheet 4: 方法与解释（2 列）----
  df_methods <- data.frame(
    `项目` = vapply(methods_rows[-1], `[[`, "", 1),
    `说明` = vapply(methods_rows[-1], `[[`, "", 2),
    check.names = FALSE, stringsAsFactors = FALSE
  )

  out_xlsx <- file.path(tables_dir, "cell_type_kla_vs_reference_ddr_statistics.xlsx")
  writexl::write_xlsx(list(
    `DDR占比对照` = df_comparison,
    `对照选择信息` = df_control,
    `TALL替代敏感性` = df_tall,
    `方法与解释` = df_methods
  ), out_xlsx)

  # ---- zh csv 输出（数值 = CSV 原串透传；已实测全部为 JS 规范形）----
  comparison_zh_headers <- c(
    "细胞或组织", "参考PXD", "参考样本子集", "参考蛋白数", "参考DDR蛋白数", "参考DDR占比",
    "Kla蛋白数", "Kla与DDR交集蛋白数", "Kla DDR占比", "Kla减参考百分点差",
    "Kla与参考占比倍数", "选择理由", "主要限制"
  )
  comparison_zh_rows <- lapply(seq_len(nrow(comparison)), function(i) {
    row <- comparison[i, ]
    ci <- zh[[row$CellOrTissueType]]
    c(ci$name, row$ReferencePXD, ci$subset,
      row$ReferenceProteinCount, row$ReferenceDdrProteinCount, row$ReferenceDdrFraction,
      row$TotalKlaProteins, row$KlaGoDdrProteins, row$KlaGoDdrFraction,
      row$DdrFractionPercentagePointDifference, row$KlaToReferenceDdrFractionRatio,
      ci$rationale, ci$caveat)
  })

  control_zh_headers <- c(
    "细胞或组织", "参考PXD", "年份", "匹配类型", "适用等级", "选用样本或子集",
    "正常或基线条件", "样本或采集数", "对照蛋白数", "DDR蛋白数", "DDR占比",
    "采集方式或仪器", "检索与定量", "是否PTM富集", "选择理由", "主要限制",
    "蛋白数计数口径", "实际分析来源文件", "仓库与分析完整度", "可作统计差异",
    "数据集URL", "论文URL", "处理数据URL"
  )
  control_zh_rows <- lapply(seq_len(nrow(control_info)), function(i) {
    row <- control_info[i, ]
    ci <- ctl[[row[["细胞或组织"]]]]
    c(ci$name, row[["参考PXD"]], row[["年份"]], ci$matchType, row[["适用等级"]],
      ci$subset, ci$baseline, ci$sampleCount,
      row[["对照蛋白数"]], row[["DDR蛋白数"]], row[["DDR占比"]],
      ci$acquisition, ci$search, ci$ptm, ci$rationale, ci$caveat, ci$countDetail,
      row[["实际分析来源文件"]], ci$completeness, ci$statistics,
      row[["数据集URL"]], row[["论文URL"]], row[["处理数据URL"]])
  })

  out_zh_1 <- file.path(tables_dir, "cell_type_kla_vs_reference_ddr_statistics_zh.csv")
  out_zh_2 <- file.path(tables_dir, "cell_type_reference_control_information_zh.csv")
  write_mjs_csv(out_zh_1, comparison_zh_headers, comparison_zh_rows, bom = TRUE)
  write_mjs_csv(out_zh_2, control_zh_headers, control_zh_rows, bom = TRUE)

  cat("workbook:", basename(out_xlsx), "rows:", nrow(comparison), "\n")
  c(out_xlsx, out_zh_1, out_zh_2)
}

# ============================================================
# 工具 4: build_reference_proteome_selection_workbook.mjs
#   → cell_type_reference_proteome_selection.xlsx
#     + cell_type_reference_proteome_selection.csv
#     + reference_proteome_pxd030304_sample_audit.csv
# ============================================================

main_columns <- list(
  c("cell_type", "当前细胞/组织类型", 20),
  c("exact_identity", "准确模型身份", 38),
  c("kla_pxd", "当前Kla来源PXD", 26),
  c("total_kla_proteins", "Kla蛋白数", 13),
  c("kla_go_ddr_proteins", "Kla∩DDR蛋白数", 16),
  c("kla_go_ddr_fraction", "Kla∩DDR/Kla", 15),
  c("reference_pxd", "推荐常规蛋白组PXD", 19),
  c("reference_subset", "应使用的样本子集", 48),
  c("match_type", "匹配类型", 26),
  c("reference_year", "参考年份", 12),
  c("baseline_definition", "正常/基线定义", 44),
  c("selection_rationale", "选择理由", 60),
  c("reference_sample_count", "参考样本/采集数", 34),
  c("reference_protein_count_main", "对照蛋白数（主计数）", 18),
  c("reference_protein_count_detail", "蛋白数计数口径", 58),
  c("reference_protein_count_source", "蛋白数来源文件", 52),
  c("proteome_depth", "蛋白组深度", 42),
  c("acquisition", "采集方式/仪器", 34),
  c("search_quantification", "检索与定量", 40),
  c("ptm_enrichment", "是否PTM富集", 15),
  c("raw_data_status", "原始数据情况", 42),
  c("processed_data_status", "处理结果情况", 48),
  c("repository_completeness", "仓库与分析完整度", 48),
  c("recommended_file", "优先读取文件", 52),
  c("download_priority", "下载优先级", 44),
  c("suitability_grade", "适用等级", 12),
  c("use_for_detection_background", "可作检出背景", 18),
  c("use_for_statistical_differential", "可作统计差异", 34),
  c("main_caveat", "主要限制", 56),
  c("backup_reference", "备用数据", 48),
  c("dataset_url", "数据集URL", 48),
  c("paper_url", "论文URL", 44),
  c("processed_data_url", "处理数据URL", 48),
  c("selection_status", "选择状态", 34),
  c("decision_date", "核验日期", 14)
)

audit_rows <- list(
  c("HEK293T", "Control_HEK293T_lys", "Control_HEK293T", "401", "主要检出背景", "未处理裂解物过程控制；不能当作401个生物重复"),
  c("HEK293T", "Control_HEK293T_std_H002", "Control_HEK293T", "663", "仅用于技术稳定性检查", "标准QC，不进入主要蛋白集合"),
  c("MCF7", "MCF7", "SIDM00148", "6", "主要对照", "准确细胞系"),
  c("HCT116", "HCT-116", "SIDM00783", "6", "主要对照", "名称需映射 HCT116 -> HCT-116"),
  c("MDA-MB-468", "MDA-MB-468", "SIDM00628", "6", "主要对照", "准确细胞系"),
  c("T-47D", "T47D", "SIDM00097", "6", "主要对照", "名称需映射 T-47D -> T47D"),
  c("RKO", "RKO", "SIDM01090", "6", "主要对照", "准确细胞系"),
  c("T-ALL", "TALL-1", "SIDM00370", "6", "主要替代对照", "与TALL-104不是同一细胞系"),
  c("T-ALL", "Jurkat", "SIDM01016", "6", "敏感性替代对照", "与TALL-104不是同一细胞系")
)

dataset_rows <- list(
  c("PXD030304",
    "HEK293T、MCF7、HCT116、T-ALL替代、MDA-MB-468、T-47D、RKO",
    "PXD030304为PARTIAL；论文配套Figshare提供最终矩阵和映射",
    "6,981次采集（仓库说明）；分析矩阵README列出6,864次运行",
    "6,692蛋白高置信矩阵；8,498蛋白敏感性矩阵",
    "mapping_file_averaged / mapping_file_replicates / protein_matrix_6692 / protein_matrix_8498",
    "先下载处理矩阵，不下载全部原始DIA",
    "https://proteomecentral.proteomexchange.org/?pxid=PXD030304",
    "https://doi.org/10.6084/m9.figshare.19345397"),
  c("PXD043880",
    "Human hippocampus",
    "PARTIAL，但74个raw和74x2,092蛋白处理矩阵可用",
    "74名神经学正常供者",
    "2,092个蛋白特征",
    "13024_2023_650_MOESM1_ESM.xlsx / Source Data Proteins",
    "先使用补充矩阵；统一重检索时再下载raw",
    "https://proteomecentral.proteomexchange.org/?pxid=PXD043880",
    "https://doi.org/10.1186/s13024-023-00650-3"),
  c("PXD072220",
    "HK-2",
    "PARTIAL，但9个raw和Spectronaut关键报告可用",
    "3 control + 3 copper + 3 TTM",
    "4,933蛋白；89-90%完整度",
    "HK-2_Spectronaut-report_PG_Quantity.txt",
    "只读取3个Control样本",
    "https://proteomecentral.proteomexchange.org/?pxid=PXD072220",
    "https://doi.org/10.1152/ajpcell.00311.2026"),
  c("PXD002400",
    "MCF10A；MCF7备用",
    "PARTIAL；106个raw和一个261.6MB MaxQuant结果包",
    "MCF10A 10个IEF组分x2次进样=20个raw",
    "深度分级DDA蛋白组",
    "msms.zip",
    "先下载MaxQuant结果包",
    "https://proteomecentral.proteomexchange.org/?pxid=PXD002400",
    "https://www.ebi.ac.uk/pride/archive/projects/PXD002400"),
  c("PXD027472",
    "HEK293T备用",
    "jPOST数据；whole-cell lysate和crude membrane两组",
    "样本数需在下载前复核",
    "全细胞与膜蛋白组",
    "whole-cell-lysate arm",
    "只使用全细胞裂解物，排除膜富集组",
    "https://proteomecentral.proteomexchange.org/?pxid=PXD027472",
    "https://doi.org/10.1016/j.mcpro.2022.100206"),
  c("PXD028488",
    "TALL-104同细胞系技术敏感性检查",
    "现有Kla研究本身，不是独立外部基线",
    "TALL-104 enriched/non-enriched",
    "乳酸暴露后的非富集蛋白组",
    "TALL-Nonenrichment",
    "不能作为正常主对照；只作同细胞系敏感性检查",
    "https://proteomecentral.proteomexchange.org/?pxid=PXD028488",
    "https://doi.org/10.1038/s41592-022-01523-1")
)

rule_rows <- list(
  c("正常的定义", "同一细胞系/组织、未处理或使用明确的 untreated/control 子集；癌细胞系的“正常”不代表非癌，只代表该细胞系的基线状态。"),
  c("无偏的定义", "不做Kla、磷酸化、乙酰化、泛素化、免疫沉淀、BioID或膜蛋白富集。允许为了提高深度进行普通肽段分级。"),
  c("主比较层级", "仅在蛋白层面比较。常规蛋白组没有Kla位点信息，不能用于判断某个赖氨酸位点是否未乳酸化。"),
  c("主背景集合", "每种模型的参考蛋白组中通过质量控制并被检出的BaseAccession集合。"),
  c("PXD030304主规则", "使用6,692蛋白矩阵作为主结果；8,498蛋白矩阵作为敏感性结果。"),
  c("Kla未检出解释", "常规蛋白组检测到、Kla表未出现的蛋白只能称为“未在当前Kla实验中鉴定”，不能称为“不乳酸化”。"),
  c("参考组未检出解释", "Kla蛋白不在参考蛋白组中可能来自低丰度、批次、仪器、搜索库或富集增敏，不应直接删除。"),
  c("UniProt规则", "继续使用去isoform后缀的BaseAccession优先匹配，GeneSymbol仅作辅助。"),
  c("TALL-104规则", "TALL-1和Jurkat结果必须分别输出；两者共同支持时才标记为T-ALL替代背景稳定检出。"),
  c("禁止的统计", "不能把技术进样、IEF组分或HEK293T QC运行数当作独立生物重复进行显著性检验。"),
  c("建议输出", "每个模型保存reference_detected、Kla∩reference、Kla-reference、DDR归一化比例和映射失败表。")
)

build_reference_proteome_selection <- function() {
  tables_dir <- p("reanalysis/results/tables")
  selection <- fromJSON(p("reanalysis/config/reference_proteome_selection.json"), simplifyVector = FALSE)
  stats <- read_csv_chr(file.path(tables_dir, "cell_type_kla_ddr_statistics.csv"))
  stats_by_type <- split(stats, stats$CellOrTissueType)

  if (length(selection) != 10) stop("Expected 10 selected cell/tissue types, found ", length(selection))

  keys <- vapply(main_columns, `[[`, "", 1)
  headers <- vapply(main_columns, `[[`, "", 2)

  rows <- lapply(selection, function(row) {
    current <- stats_by_type[[row$cell_type]]
    if (is.null(current)) stop("Missing current Kla statistics for ", row$cell_type)
    out <- list()
    for (key in keys) {
      out[[key]] <- switch(key,
        total_kla_proteins = current$TotalKlaProteins,
        kla_go_ddr_proteins = current$KlaGoDdrProteins,
        kla_go_ddr_fraction = current$KlaGoDdrFraction,
        decision_date = "2026-07-30",
        if (is.null(row[[key]])) "" else row[[key]]
      )
    }
    out
  })

  # ---- csv 输出（35 列，无辅助列）----
  csv_rows <- lapply(rows, function(r) {
    vapply(keys, function(k) {
      v <- r[[k]]
      if (is.null(v) || is.na(v)) "" else as.character(v)
    }, "")
  })
  out_csv <- file.path(tables_dir, "cell_type_reference_proteome_selection.csv")
  write_mjs_csv(out_csv, headers, csv_rows, bom = TRUE)

  audit_headers <- c("当前类型", "PXD030304标签", "SIDM/项目标识", "运行数", "用途", "说明")
  out_audit <- file.path(tables_dir, "reference_proteome_pxd030304_sample_audit.csv")
  write_mjs_csv(out_audit, audit_headers, audit_rows, bom = TRUE)

  # ---- xlsx ----
  # sheet 1: 10类对照选择（35 主列 + 辅助列）
  num_cols <- c("total_kla_proteins", "kla_go_ddr_proteins", "reference_year",
                "reference_protein_count_main")
  main_data <- lapply(keys, function(k) {
    vals <- vapply(rows, function(r) {
      v <- r[[k]]
      if (is.null(v) || is.na(v)) return(NA_character_)
      as.character(v)
    }, "")
    if (k %in% num_cols) as.numeric(vals) else vals
  })
  names(main_data) <- headers
  # F 列（Kla∩DDR/Kla）= .mjs 公式 =E2/D2 的缓存值（JS 最短回读文本）
  total <- as.numeric(main_data[["Kla蛋白数"]])
  ddr <- as.numeric(main_data[["Kla∩DDR蛋白数"]])
  main_data[["Kla∩DDR/Kla"]] <- js_fraction_str(ddr, total)
  # 辅助列：参考PXD 首现
  pxd_seq <- main_data[["推荐常规蛋白组PXD"]]
  helper <- as.numeric(!duplicated(pxd_seq))
  main_data[["QC辅助：PXD首次出现"]] <- helper
  df_main <- as.data.frame(main_data, stringsAsFactors = FALSE, check.names = FALSE)

  # sheet 2: PXD030304样本核验
  df_audit <- data.frame(
    `当前类型` = vapply(audit_rows, `[[`, "", 1),
    `PXD030304标签` = vapply(audit_rows, `[[`, "", 2),
    `SIDM/项目标识` = vapply(audit_rows, `[[`, "", 3),
    `运行数` = as.numeric(vapply(audit_rows, `[[`, "", 4)),
    `用途` = vapply(audit_rows, `[[`, "", 5),
    `说明` = vapply(audit_rows, `[[`, "", 6),
    check.names = FALSE, stringsAsFactors = FALSE
  )

  # sheet 3: 数据集与比较规则
  # writexl 的表头行 = data.frame 列名 → 对应 JS 第 1 行（A1=选择摘要 B1=数量，C1:I1 为空）；
  # 数据行 df 第 i 行 ↔ JS 第 i+1 行：1-5 摘要、6-7 空、8 数据集表头、9-14 数据集、
  # 15-16 空、17 规则表头、18-28 规则（共 28 数据行）
  m <- matrix(NA_character_, nrow = 28, ncol = 9)
  m[1, 1] <- "当前细胞/组织类型"
  m[2, 1] <- "精确匹配"
  m[3, 1] <- "替代匹配"
  m[4, 1] <- "A级"
  m[5, 1] <- "主要参考PXD"
  exact <- sum(vapply(rows, function(r) r$match_type == "Exact cell-line match", logical(1)))
  tissue <- sum(vapply(rows, function(r) r$match_type == "Exact tissue match", logical(1)))
  surrogate <- sum(vapply(rows, function(r) r$match_type == "Disease-matched surrogate, not an exact cell-line match", logical(1)))
  grade_a <- sum(vapply(rows, function(r) r$suitability_grade == "A", logical(1)))
  m[1, 2] <- as.character(length(rows))
  m[2, 2] <- as.character(exact + tissue)
  m[3, 2] <- as.character(surrogate)
  m[4, 2] <- as.character(grade_a)
  m[5, 2] <- as.character(sum(helper))
  dataset_headers <- c("PXD", "覆盖模型", "完整度说明", "样本/采集规模", "蛋白组深度",
                       "优先文件", "下载策略", "数据集URL", "处理数据/论文URL")
  m[8, ] <- dataset_headers
  for (i in seq_along(dataset_rows)) m[8 + i, ] <- dataset_rows[[i]]
  m[17, 1:2] <- c("规则", "执行说明")
  for (i in seq_along(rule_rows)) m[17 + i, 1:2] <- rule_rows[[i]]
  df_rules <- as.data.frame(m, stringsAsFactors = FALSE)
  # 表头行须与 JS 版第 1 行一致（A1=选择摘要 B1=数量，C1:I1 为空）
  names(df_rules) <- c("选择摘要", "数量", "", "", "", "", "", "", "")

  out_xlsx <- file.path(tables_dir, "cell_type_reference_proteome_selection.xlsx")
  writexl::write_xlsx(list(
    `10类对照选择` = df_main,
    `PXD030304样本核验` = df_audit,
    `数据集与比较规则` = df_rules
  ), out_xlsx)

  cat("workbook:", basename(out_xlsx), "rows:", length(rows), "\n")
  c(out_xlsx, out_csv, out_audit)
}

# ============================================================
# 工具 5: build_venn_combined_workbook.mjs
#   → venn_combined_tables.xlsx
# ============================================================
build_venn_combined_xlsx <- function() {
  tables_dir <- p("reanalysis/results/tables")
  sources <- list(
    list(sheetName = "Kla_unique", fileName = "all_kla_three_groups_combined.csv"),
    list(sheetName = "Kla_non_dedup", fileName = "all_kla_three_groups_combined_non_deduplicated.csv"),
    list(sheetName = "Kla_DDR_unique", fileName = "kla_go_ddr_three_groups_combined.csv"),
    list(sheetName = "Kla_DDR_non_dedup", fileName = "kla_go_ddr_three_groups_combined_non_deduplicated.csv")
  )
  sheets <- list()
  counts <- integer(4)
  for (i in seq_along(sources)) {
    df <- read_csv_chr(file.path(tables_dir, sources[[i]]$fileName))
    counts[[i]] <- nrow(df)
    sheets[[sources[[i]]$sheetName]] <- df
  }

  meaning <- c(
    "All Kla proteins; one row per BaseAccession",
    "Three Kla groups appended; intersection proteins repeat by SourceCategory",
    "Kla intersected with GO repair/damage; one row per BaseAccession",
    "Three Kla-DDR groups appended; intersection proteins repeat by SourceCategory"
  )
  # 表头行（df 列名）↔ JS 第 3 行；df 第 i 行 ↔ JS 第 i+2 行：
  # df 1-4 = 四个 sheet 行（JS 4-7）、5-6 空（JS 8-9）、7 = Source CSV files（JS 10）、8-11 = 文件名（JS 11-14）。
  # （JS 第 1 行合并标题为呈现层，省略；JS 第 2 行为空行，省略）
  readme <- matrix(NA_character_, nrow = 11, ncol = 3)
  for (i in seq_len(4)) {
    readme[i, ] <- c(sources[[i]]$sheetName, as.character(counts[[i]]), meaning[[i]])
  }
  readme[7, 1] <- "Source CSV files"
  for (i in seq_len(4)) readme[7 + i, 1] <- sources[[i]]$fileName
  df_readme <- as.data.frame(readme, stringsAsFactors = FALSE)
  names(df_readme) <- c("Sheet", "Rows", "Meaning")
  sheets[["README"]] <- df_readme

  out <- file.path(tables_dir, "venn_combined_tables.xlsx")
  writexl::write_xlsx(sheets, out)
  cat("workbook:", basename(out), "rows:", paste(counts, collapse = "/"), "\n")
  out
}

# ============================================================
# 工具 6: create_bilingual_figure_legends_docx.js
#   → Kla_Venn_figure_legends_bilingual.docx
# ============================================================

figure_legends <- list(
  list(
    section = "Figure X / 图 X",
    img = "all_kla_three_group_venn.png",
    width = 500, height = 449,
    alt = "Venn diagram of all Kla proteins",
    en_title = "Figure X. Distribution of lysine-lactylated proteins among hippocampal tissue, immortalized models, and tumor cell lines. ",
    en_body = "Venn diagram showing the overlap of 3,112 lysine-lactylated (Kla) proteins identified across human hippocampal tissue, immortalized/non-tumor cell models (HEK293T, HK-2, and MCF10A), and tumor cell lines (MCF7, HCT116, T-ALL, MDA-MB-468, T-47D, and RKO). Proteins were collapsed to unique UniProt base accessions after removal of isoform suffixes. Numbers indicate the counts of unique proteins in each mutually exclusive region.",
    zh_title = "图 X. 人海马组织、永生化模型与肿瘤细胞系中赖氨酸乳酰化蛋白的分布。",
    zh_body = "维恩图展示在人海马组织、永生化/非肿瘤细胞模型（HEK293T、HK-2 和 MCF10A）及肿瘤细胞系（MCF7、HCT116、T-ALL、MDA-MB-468、T-47D 和 RKO）中鉴定到的 3,112 个赖氨酸乳酰化（Kla）蛋白的重叠关系。去除异构体后缀后，以唯一 UniProt 基础登录号（BaseAccession）作为蛋白计数单位。数字表示各互斥区域内的唯一蛋白数。"
  ),
  list(
    section = "Figure Y / 图 Y",
    img = "kla_go_ddr_three_group_venn.png",
    width = 478, height = 466,
    alt = "Venn diagram of Kla proteins associated with DNA repair and damage responses",
    en_title = "Figure Y. Distribution of Kla proteins associated with DNA repair and DNA damage responses among hippocampal tissue, immortalized models, and tumor cell lines. ",
    en_body = "Venn diagram showing the overlap of 275 Kla proteins annotated to DNA repair- or DNA damage response-related Gene Ontology (GO) biological processes across human hippocampal tissue, immortalized/non-tumor cell models (HEK293T, HK-2, and MCF10A), and tumor cell lines (MCF7, HCT116, T-ALL, MDA-MB-468, T-47D, and RKO). Kla proteins were matched to human GO annotations primarily by UniProt base accession, with gene-symbol matching used only when accession-level matching was unavailable; annotations carrying the qualifier “NOT” were excluded. Numbers indicate the counts of unique proteins in each mutually exclusive region.",
    zh_title = "图 Y. 人海马组织、永生化模型与肿瘤细胞系中 DNA 修复和 DNA 损伤应答相关 Kla 蛋白的分布。",
    zh_body = "维恩图展示在人海马组织、永生化/非肿瘤细胞模型（HEK293T、HK-2 和 MCF10A）及肿瘤细胞系（MCF7、HCT116、T-ALL、MDA-MB-468、T-47D 和 RKO）中鉴定到的 275 个 Kla 蛋白的重叠关系；这些蛋白被注释至 DNA 修复或 DNA 损伤应答相关的基因本体（GO）生物学过程。Kla 蛋白优先通过 UniProt 基础登录号与人源 GO 注释匹配，仅在登录号无法匹配时使用基因符号进行辅助匹配；带有“NOT”限定词的注释被排除。数字表示各互斥区域内的唯一蛋白数。"
  )
)

# 页脚页码字段：officer 0.7.3 原生支持（prop_section(footer_default = ...)）
footer_content <- block_list(
  fpar(
    run_word_field(field = "PAGE",
                   prop = fp_text(color = "777777", font.size = 9,
                                  font.family = "Times New Roman",
                                  eastasia.family = "SimSun")),
    fp_p = fp_par(text.align = "center")
  )
)

build_figure_legends_docx <- function() {
  figures_dir <- p("reanalysis/results/figures")
  out <- p("reanalysis/results/Kla_Venn_figure_legends_bilingual.docx")

  doc <- read_docx()
  doc <- set_doc_properties(
    doc,
    title = "Bilingual figure legends for Kla Venn diagrams",
    creator = "Codex",
    description = "Publication-ready English and Chinese legends for two Kla Venn diagrams."
  )

  heading_font <- fp_text(bold = TRUE, font.size = 16, font.family = "Arial",
                          eastasia.family = "SimHei")
  chinese_font <- function(size, ...) fp_text(..., font.size = size,
                                              font.family = "Times New Roman",
                                              eastasia.family = "SimSun")
  english_font <- function(size, ...) fp_text(..., font.size = size,
                                              font.family = "Times New Roman",
                                              eastasia.family = "SimSun")

  # 标题
  doc <- body_add_fpar(doc, fpar(
    ftext("Bilingual Figure Legends / 双语图注", prop = heading_font),
    fp_p = fp_par(text.align = "center")
  ))
  # 注释
  doc <- body_add_fpar(doc, fpar(
    ftext("注：Figure X/Y 与图 X/Y 为占位符，使用时请替换为稿件中的实际图号。",
          prop = chinese_font(9, italic = TRUE, color = "666666")),
    fp_p = fp_par(text.align = "center")
  ))

  for (i in seq_along(figure_legends)) {
    fl <- figure_legends[[i]]
    img_path <- file.path(figures_dir, fl$img)
    if (!file.exists(img_path)) stop("Missing figure image: ", img_path)
    # 小节（图 Y 前分页）
    par_children <- list()
    if (i > 1) par_children[[length(par_children) + 1]] <- run_pagebreak()
    par_children[[length(par_children) + 1]] <- ftext(fl$section, prop = fp_text(
      bold = TRUE, font.size = 12, font.family = "Arial", eastasia.family = "SimHei"))
    doc <- body_add_fpar(doc, do.call(fpar, c(par_children, list(fp_p = fp_par(text.align = "left")))))
    # 图片（px → in: 1px = 9525 EMU, 1 in = 914400 EMU）
    doc <- body_add_fpar(doc, fpar(
      external_img(img_path, width = fl$width * 9525 / 914400,
                   height = fl$height * 9525 / 914400, alt = fl$alt),
      fp_p = fp_par(text.align = "center")
    ))
    # 图注 EN / ZH
    doc <- body_add_fpar(doc, fpar(
      ftext(fl$en_title, prop = english_font(10.5, bold = TRUE)),
      ftext(fl$en_body, prop = english_font(10.5)),
      fp_p = fp_par(text.align = "justify")
    ))
    doc <- body_add_fpar(doc, fpar(
      ftext(fl$zh_title, prop = chinese_font(10.5, bold = TRUE)),
      ftext(fl$zh_body, prop = chinese_font(10.5)),
      fp_p = fp_par(text.align = "justify")
    ))
  }

  # 页面与页脚（须在全部内容之后设置默认节属性）
  doc <- body_set_default_section(doc, prop_section(
    page_size = page_size(width = 11906 / 1440, height = 16838 / 1440),
    page_margins = page_mar(top = 1134 / 1440, right = 1418 / 1440,
                            bottom = 1134 / 1440, left = 1418 / 1440),
    footer_default = footer_content
  ))

  print(doc, target = out)
  cat("docx:", basename(out), "\n")
  out
}

# ============================================================
# 工具 7: build_project_metadata.py（metadata 段）
# ============================================================

expected_pxds <- c("PXD014870", "PXD028488", "PXD038880", "PXD050470",
                   "PXD050906", "PXD053474", "PXD060185", "PXD078013", "PXD078736")
generated_metadata_files <- c("dataset_metadata.csv", "file_inventory.csv")
ephemeral_suffixes <- ".pyc"
inventory_fields <- c("PXD", "Area", "RelativePath", "CurrentPath", "FileName",
                      "Extension", "SizeBytes")
manifest_fields <- c("current_path", "size_bytes", "reconstructed_original_path",
                     "recorded_timestamp_utc", "recorded_action", "recorded_reason",
                     "provenance_evidence", "reconstruction_confidence", "notes")

read_dataset_config <- function() {
  path <- p("reanalysis/config/datasets.csv")
  rows <- read_csv_chr(path)
  fields <- names(rows)
  if (!"PXD" %in% fields) stop("Missing PXD column in ", path)
  configured <- unique(rows$PXD)
  if (!setequal(configured, expected_pxds) || nrow(rows) != length(expected_pxds)) {
    stop("datasets.csv must contain exactly the expected nine unique PXD rows; found ",
         paste(sort(configured), collapse = ", "))
  }
  list(fields = fields, rows = rows)
}

inventory_dataset <- function(pxd) {
  pxd_dir <- p("data", pxd)
  if (!dir.exists(pxd_dir)) stop("Missing dataset directory: ", pxd_dir)
  all_paths <- list.files(pxd_dir, recursive = TRUE, all.files = TRUE, full.names = TRUE)
  rel_paths <- list.files(pxd_dir, recursive = TRUE, all.files = TRUE, full.names = FALSE)
  keep <- !file.info(all_paths)$isdir
  all_paths <- all_paths[keep]
  rel_paths <- rel_paths[keep]
  ord <- order(vapply(rel_paths, function(x) x, ""), method = "radix")
  all_paths <- all_paths[ord]
  rel_paths <- rel_paths[ord]

  out <- list()
  for (i in seq_along(rel_paths)) {
    rp <- rel_paths[[i]]
    name <- basename(rp)
    if (name == ".DS_Store" || tolower(py_suffix(name)) %in% ephemeral_suffixes) next
    if (dirname(rp) == "metadata" && name %in% generated_metadata_files) next
    parts <- strsplit(rp, "/", fixed = TRUE)[[1]]
    area <- if (length(parts) > 1) parts[[1]] else "root"
    out[[length(out) + 1]] <- list(
      PXD = pxd,
      Area = area,
      RelativePath = rp,
      CurrentPath = file.path("data", pxd, rp),
      FileName = name,
      Extension = tolower(py_suffix(name)),
      SizeBytes = as.character(file.info(all_paths[[i]])$size)
    )
  }
  out
}

load_prior_manifest <- function() {
  path <- p("archive/migration_manifest_2026-07-21.csv")
  rows <- read_csv_chr(path)
  by_destination <- list()
  for (i in seq_len(nrow(rows))) {
    dest <- rows$destination[[i]]
    if (!is.null(by_destination[[dest]])) stop("Ambiguous destination in prior manifest: ", dest)
    by_destination[[dest]] <- list(
      timestamp_utc = rows$timestamp_utc[[i]],
      action = rows$action[[i]],
      source = rows$source[[i]],
      reason = rows$reason[[i]]
    )
  }
  by_destination
}

exact_manifest_provenance <- function(current_path, by_destination) {
  path <- current_path
  chain <- list()
  seen <- character(0)
  while (!is.null(by_destination[[path]])) {
    if (path %in% seen) stop("Migration cycle involving ", path)
    seen <- c(seen, path)
    record <- by_destination[[path]]
    chain[[length(chain) + 1]] <- record
    path <- record$source
  }
  if (!length(chain)) return(NULL)
  chrono <- rev(chain)
  list(
    reconstructed_original_path = path,
    recorded_timestamp_utc = paste(vapply(chrono, `[[`, "", "timestamp_utc"), collapse = "; "),
    recorded_action = paste(vapply(chrono, `[[`, "", "action"), collapse = "; "),
    recorded_reason = paste(vapply(chrono, `[[`, "", "reason"), collapse = " | "),
    provenance_evidence = "archive/migration_manifest_2026-07-21.csv",
    reconstruction_confidence = "high",
    notes = "Exact source/destination chain recorded in the prior manifest."
  )
}

inferred_readme_provenance <- function(current_path) {
  archive_prefix <- "archive/2026-07-21_pre_restructure/"
  if (startsWith(current_path, archive_prefix)) {
    relative <- substring(current_path, nchar("archive/") + 1)
    return(list(
      reconstructed_original_path = paste0("99_archive/", relative),
      recorded_timestamp_utc = "",
      recorded_action = "inferred_prefix_rename",
      recorded_reason = "Prior README identifies 99_archive as the historical archive area.",
      provenance_evidence = "archive/README_before_final_2026-07-21.md; archive/README.md",
      reconstruction_confidence = "medium",
      notes = "Directory-prefix reconstruction from README evidence; not an exact file-level move record."
    ))
  }
  reanalysis_prefix <- "archive/reanalysis_v1_2026-07-21/"
  if (startsWith(current_path, reanalysis_prefix)) {
    relative <- substring(current_path, nchar(reanalysis_prefix) + 1)
    return(list(
      reconstructed_original_path = paste0("03_reanalysis/", relative),
      recorded_timestamp_utc = "",
      recorded_action = "inferred_archive_move",
      recorded_reason = "Prior README identifies 03_reanalysis as the first-stage reanalysis area.",
      provenance_evidence = "archive/README_before_final_2026-07-21.md; archive/reanalysis_v1_2026-07-21/README.md",
      reconstruction_confidence = "medium",
      notes = "Directory-prefix reconstruction from README evidence; not an exact file-level move record."
    ))
  }
  NULL
}

reconstruct_migration_manifest <- function(by_destination) {
  out_path <- p("archive/migration_manifest_reconstructed_2026-07-22.csv")
  all_paths <- list.files(project_root, recursive = TRUE, all.files = TRUE, full.names = TRUE)
  rel_paths <- list.files(project_root, recursive = TRUE, all.files = TRUE, full.names = FALSE)
  keep <- !file.info(all_paths)$isdir
  all_paths <- all_paths[keep]
  rel_paths <- rel_paths[keep]
  ord <- order(vapply(rel_paths, function(x) x, ""), method = "radix")
  all_paths <- all_paths[ord]
  rel_paths <- rel_paths[ord]

  rows <- list()
  for (i in seq_along(rel_paths)) {
    rp <- rel_paths[[i]]
    name <- basename(rp)
    if (name == ".DS_Store" || tolower(py_suffix(name)) %in% ephemeral_suffixes) next
    current_path <- rp
    if (file.path(project_root, rp) == out_path) next
    provenance <- exact_manifest_provenance(current_path, by_destination)
    if (is.null(provenance)) provenance <- inferred_readme_provenance(current_path)
    if (is.null(provenance)) {
      provenance <- list(
        reconstructed_original_path = "",
        recorded_timestamp_utc = "",
        recorded_action = "",
        recorded_reason = "",
        provenance_evidence = "Current filesystem snapshot on 2026-07-22",
        reconstruction_confidence = "unknown",
        notes = "Current path is known; no supported original path was recovered."
      )
    }
    rows[[length(rows) + 1]] <- c(
      list(current_path = current_path,
           size_bytes = as.character(file.info(all_paths[[i]])$size)),
      provenance
    )
  }
  write_py_csv(out_path, manifest_fields, rows)
  out_path
}

build_project_metadata <- function() {
  cfg <- read_dataset_config()
  project_rows <- list()
  pxd_order <- sort(unique(cfg$rows$PXD), method = "radix")
  for (pxd in pxd_order) {
    metadata_row <- cfg$rows[cfg$rows$PXD == pxd, , drop = FALSE]
    metadata_dir <- p("data", pxd, "metadata")
    write_py_csv(file.path(metadata_dir, "dataset_metadata.csv"), cfg$fields,
                 list(as.list(metadata_row)))
    inventory_rows <- inventory_dataset(pxd)
    write_py_csv(file.path(metadata_dir, "file_inventory.csv"), inventory_fields, inventory_rows)
    project_rows <- c(project_rows, inventory_rows)
  }

  project_inventory_path <- p("reanalysis/reports/project_file_inventory.csv")
  write_py_csv(project_inventory_path, inventory_fields, project_rows)

  by_destination <- load_prior_manifest()
  manifest_path <- reconstruct_migration_manifest(by_destination)

  cat("Datasets:", length(pxd_order), "\n")
  cat("Inventoried files:", length(project_rows), "\n")
  cat("Project inventory:", project_inventory_path, "\n")
  cat("Reconstructed manifest:", manifest_path, "\n")
  project_inventory_path
}

# ============================================================
# main
# ============================================================

if (stage %in% c("all", "workbooks")) {
  build_cell_type_statistics_xlsx()
  build_hippocampus_review_md()
  build_reference_ddr_comparison()
  build_reference_proteome_selection()
  build_venn_combined_xlsx()
  build_figure_legends_docx()
}
if (stage %in% c("all", "metadata")) {
  build_project_metadata()
}

cat("build_workbooks.R done (stage:", stage, ")\n")
