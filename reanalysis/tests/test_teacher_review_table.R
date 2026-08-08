#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args)) normalizePath(args[[1]]) else normalizePath(".")
table_dir <- file.path(project_root, "reanalysis", "results", "tables")

read_table <- function(name) {
  read.csv(
    file.path(table_dir, name),
    check.names = FALSE,
    stringsAsFactors = FALSE,
    na.strings = c("", "NA")
  )
}

review <- read_table("kla_and_reference_teacher_review_zh.csv")
statistics <- read_table(
  "cell_type_kla_vs_reference_ddr_statistics_accession_only_zh.csv"
)

key <- function(data) paste(data$乳酸化PXD, data$样本组, sep = "\r")
stopifnot(nrow(review) == 37)
stopifnot(!anyDuplicated(key(review)))
stopifnot(identical(key(review), key(statistics)))
stopifnot(sum(review$是否纳入配对分析 == "是") == 33)
stopifnot(sum(review$是否纳入配对分析 == "否") == 4)
stopifnot(!any(review$乳酸化PXD == "PXD037371"))
stopifnot(all(review$Kla蛋白数 == statistics$乳酸化蛋白ID数))
stopifnot(all(review$Kla_DDR蛋白数 == statistics$乳酸化DDR蛋白ID数))
stopifnot(all.equal(
  review$Kla_DDR占比,
  statistics$乳酸化DDR占比,
  check.attributes = FALSE
))

included <- review$是否纳入配对分析 == "是"
stopifnot(all(
  review$普通全蛋白蛋白数[included] ==
    statistics$常规蛋白ID数[included]
))
stopifnot(all(
  review$普通全蛋白DDR数[included] ==
    statistics$常规DDR蛋白ID数[included]
))
stopifnot(all.equal(
  review$普通全蛋白DDR占比[included],
  statistics$常规DDR占比[included],
  check.attributes = FALSE
))
stopifnot(all(review$材料身份严格匹配[included] == "是"))
stopifnot(all(review$材料身份严格匹配[!included] == "否"))

split_paths <- function(value) {
  paths <- trimws(unlist(strsplit(value, ";", fixed = TRUE)))
  paths[nzchar(paths)]
}
assert_paths_exist <- function(values, label) {
  paths <- unique(unlist(lapply(values, split_paths)))
  missing <- paths[!file.exists(file.path(project_root, paths))]
  if (length(missing)) {
    stop(label, " has missing paths: ", paste(missing, collapse = ", "))
  }
}
assert_paths_exist(review$Kla实际使用文件, "Kla evidence")
assert_paths_exist(
  review$普通全蛋白实际使用文件[included],
  "Reference evidence"
)

hippocampus <- review[
  review$乳酸化PXD == "PXD050470" &
    review$样本组 == "human hippocampus",
]
stopifnot(nrow(hippocampus) == 1)
stopifnot(hippocampus$普通全蛋白PXD == "PXD050470")
stopifnot(
  hippocampus$普通全蛋白实际使用文件 ==
    "data/PXD050470/supplementary/prca2331-sup-0006-tables4.xlsx"
)
stopifnot(hippocampus$实际读取样本_工作表_列 == "H072;H081;H0187")
stopifnot(hippocampus$普通全蛋白蛋白数 == 6082)
stopifnot(hippocampus$普通全蛋白DDR数 == 219)

hcc <- review[review$乳酸化PXD == "PXD075377", ]
stopifnot(grepl(
  "^CISs sheet",
  hcc$实际读取样本_工作表_列[hcc$样本组 == "HCC"]
))
stopifnot(grepl(
  "^ANTs sheet",
  hcc$实际读取样本_工作表_列[hcc$样本组 == "adjacent liver"]
))

scar <- review[review$乳酸化PXD == "PXD046800", ]
stopifnot(
  scar$实际读取样本_工作表_列[scar$样本组 == "hypertrophic scar"] ==
    "HSP1;HSP2;HSP3;HSP4"
)
stopifnot(
  scar$实际读取样本_工作表_列[scar$样本组 == "adjacent skin"] ==
    "4 adjacent skin samples"
)

prostate <- review[review$乳酸化PXD == "PXD066054", ]
stopifnot(grepl(
  "^NAT1 to NAT5",
  prostate$实际读取样本_工作表_列[prostate$样本组 == "BPH"]
))
stopifnot(grepl(
  "^PCa1 to PCa5",
  prostate$实际读取样本_工作表_列[prostate$样本组 == "prostate cancer"]
))

hk2 <- review[
  review$乳酸化PXD %in% c("PXD058534", "PXD078736"),
]
stopifnot(nrow(hk2) == 2)
stopifnot(all(hk2$普通全蛋白PXD == "PXD072220"))
stopifnot(all(grepl("amostra1", hk2$实际读取样本_工作表_列)))

huvec <- review[
  review$乳酸化PXD == "PXD073311" &
    review$样本组 == "HUVEC control and Pg infection",
]
stopifnot(nrow(huvec) == 1)
stopifnot(huvec$普通全蛋白PXD == "PXD073311")
stopifnot(grepl("report.pg_matrix.tsv$", huvec$普通全蛋白实际使用文件))
stopifnot(grepl("A0h_1", huvec$实际读取样本_工作表_列))

mcf7 <- review[
  review$乳酸化PXD == "PXD014870" &
    grepl("MCF7", review$样本组),
]
stopifnot(nrow(mcf7) >= 1)
stopifnot(all(mcf7$普通全蛋白PXD == "PXD030304"))
stopifnot(all(statistics$GeneSymbol回退数 == 0))

message("Teacher review table tests passed.")
