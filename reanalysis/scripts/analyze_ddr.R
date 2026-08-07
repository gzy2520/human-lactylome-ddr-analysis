#!/usr/bin/env Rscript
# analyze_ddr.R - Kla reanalysis main pipeline (R port of run_pipeline.py)
#
# Usage: Rscript analyze_ddr.R <project_root> [--stage pipeline]
#
# This file is the single entry point for the kla reanalysis pipeline.
# It currently implements the `pipeline` stage (Python equivalent:
# run_pipeline.py + extractors.py + common.py + audit_target_sources.py).
# Later tasks append the expanded / reference / audit stages as additional
# branches of run_stage(). Keep run_stage() as the dispatcher so future
# stages slot in without restructuring.

suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
  library(readxl)
  library(ggVennDiagram)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

# ---------------------------------------------------------------------------
# Constants (mirror common.py / run_pipeline.py / extractors.py)
# ---------------------------------------------------------------------------

SITE_COLUMNS <- c(
  "PXD", "DOI", "SampleName", "CellOrTissueType", "ExperimentalGroup",
  "Replicate", "EnrichmentStatus", "AcquisitionMode", "Accession",
  "BaseAccession", "GeneSymbol", "ProteinName", "KlaSite", "ModifiedPeptide",
  "LocalizationProb", "PEP", "Score", "DiagnosticIon",
  "DiagnosticIonIntensity", "ClassI", "SourceFile", "SourceRow", "SiteID",
  "EvidenceMode", "SourceConfidence", "PrimaryIncluded", "InclusionStatus",
  "ExclusionReason"
)
NUMERIC_SITE_COLUMNS <- c("LocalizationProb", "PEP", "Score", "DiagnosticIonIntensity")
LOGICAL_SITE_COLUMNS <- c("ClassI", "PrimaryIncluded")

INCLUDED_PXDS <- c(
  "PXD014870", "PXD028488", "PXD050470", "PXD053474",
  "PXD060185", "PXD078013", "PXD078736"
)
EXCLUDED_PXDS <- c("PXD038880", "PXD050906")
CATEGORIES <- c("hippocampus_tissue", "normal_immortalized_cell_lines", "tumor_cell_lines")
REGION_ORDER <- c(
  "hippocampus_only", "normal_only", "tumor_only",
  "hippocampus_and_normal_only", "hippocampus_and_tumor_only",
  "normal_and_tumor_only", "all_three"
)

# pandas read_csv(dtype=str) default NA values
PANDAS_NA <- c(
  "", "NA", "NaN", "N/A", "N/A N/A", "#N/A", "#N/A N/A", "#NA",
  "-1.#IND", "-1.#QNAN", "-NaN", "-nan", "1.#IND", "1.#QNAN",
  "<NA>", "None", "n/a", "nan", "null", "NULL"
)

# ---------------------------------------------------------------------------
# lib (shared R utilities, do not redefine)
# ---------------------------------------------------------------------------

# Locate the scripts/lib directory. When run via Rscript the --file=
# argument points at this script; when sourced from a test script the
# working directory (project root or reanalysis/) is used as fallback.
find_lib_dir <- function() {
  args <- commandArgs(FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  script_dir <- if (length(file_arg)) dirname(sub("^--file=", "", file_arg[1])) else ""
  candidates <- c(
    file.path(script_dir, "lib"),
    file.path(getwd(), "lib"),
    file.path(getwd(), "scripts", "lib"),
    file.path(getwd(), "reanalysis", "scripts", "lib"),
    file.path(dirname(getwd()), "scripts", "lib")
  )
  found <- candidates[vapply(candidates, file.exists, logical(1))]
  if (!length(found)) stop("Cannot locate lib/ directory for analyze_ddr.R")
  found[1]
}
lib_dir <- find_lib_dir()
source(file.path(lib_dir, "accession_utils.R"))
source(file.path(lib_dir, "io_utils.R"))
# base_accession comes from lib/accession_utils.R
# relative_path comes from lib/io_utils.R

# ---------------------------------------------------------------------------
# common.py translation
# ---------------------------------------------------------------------------

clean_text <- function(value) {
  if (length(value) == 0) return("")
  if (is.na(value)) return("")
  text <- trimws(as.character(value))
  if (tolower(text) %in% c("nan", "none")) "" else text
}

number <- function(value) {
  text <- clean_text(value)
  if (!nzchar(text)) return(NA_real_)
  suppressWarnings(as.numeric(text))
}

integer_value <- function(value) {
  parsed <- number(value)
  if (is.na(parsed)) NA_integer_ else as.integer(parsed)
}

is_true <- function(value) {
  tolower(clean_text(value)) %in% c("+", "1", "true", "yes", "y")
}

split_tokens <- function(value, pattern = ";") {
  text <- clean_text(value)
  if (!nzchar(text)) return(character(0))
  tokens <- strsplit(text, pattern, perl = TRUE)[[1]]
  tokens <- trimws(tokens)
  tokens[nzchar(tokens)]
}

normalize_accession <- function(value) {
  text <- clean_text(value)
  text <- gsub("^(?:REV__|CON__)+", "", text, perl = TRUE)
  text <- gsub("^(?:sp|tr)\\|", "", text, perl = TRUE)
  if (grepl("|", text, fixed = TRUE)) {
    text <- sub("\\|.*$", "", text)
  }
  trimws(text)
}

accession_candidates <- function(value) {
  raw <- clean_text(value)
  candidates <- character(0)
  if (nzchar(raw)) {
    for (token in strsplit(raw, "[;:]", perl = TRUE)[[1]]) {
      accession <- normalize_accession(token)
      if (nzchar(accession) && !(accession %in% candidates)) {
        candidates <- c(candidates, accession)
      }
    }
  }
  if (!length(candidates)) {
    accession <- normalize_accession(raw)
    if (nzchar(accession)) candidates <- accession
  }
  candidates
}

unique_join <- function(values) {
  tokens <- character(0)
  for (value in values) {
    tokens <- c(tokens, split_tokens(value))
  }
  tokens <- unique(tokens)
  if (!length(tokens)) "" else paste(sort(tokens, method = "radix"), collapse = ";")
}

best_annotation <- function(values) {
  cleaned <- vapply(values, clean_text, character(1))
  cleaned <- cleaned[nzchar(cleaned)]
  if (!length(cleaned)) return("")
  counts <- table(cleaned)
  counts_vec <- as.integer(counts)
  ord <- order(-counts_vec, -nchar(names(counts_vec)), names(counts_vec), method = "radix")
  names(counts_vec)[ord[1]]
}

annotation_from_description <- function(value) {
  description <- clean_text(value)
  gene <- ""
  m <- regexec("(?:^|\\s)GN=([^\\s]+)", description, perl = TRUE)
  mm <- regmatches(description, m)[[1]]
  if (length(mm) >= 2) gene <- mm[2]
  protein_name <- sub("\\s+OS=.*$", "", description)
  c(gene, protein_name)
}

annotation_from_fasta <- function(value, accession) {
  for (header in split_tokens(value)) {
    parts <- strsplit(header, " ", fixed = TRUE)[[1]]
    token <- parts[1]
    if (base_accession(token) != base_accession(accession)) next
    rest <- if (length(parts) > 1) sub("^[^ ]+ ", "", header) else ""
    return(annotation_from_description(rest))
  }
  c("", "")
}

blank_site <- function(...) {
  row <- list(
    PXD = "", DOI = "", SampleName = "", CellOrTissueType = "",
    ExperimentalGroup = "", Replicate = "", EnrichmentStatus = "",
    AcquisitionMode = "", Accession = "", BaseAccession = "", GeneSymbol = "",
    ProteinName = "", KlaSite = "", ModifiedPeptide = "",
    LocalizationProb = NA_real_, PEP = NA_real_, Score = NA_real_,
    DiagnosticIon = "", DiagnosticIonIntensity = NA_real_, ClassI = FALSE,
    SourceFile = "", SourceRow = NA_integer_, SiteID = "", EvidenceMode = "",
    SourceConfidence = "", PrimaryIncluded = TRUE, InclusionStatus = "included",
    ExclusionReason = ""
  )
  updates <- list(...)
  row[names(updates)] <- updates
  row
}

strip_peptide_modifications <- function(value) {
  peptide <- clean_text(value)
  m <- regexec("^[A-Z]\\.([A-Z].*)\\.[A-Z]$", peptide, perl = TRUE)
  mm <- regmatches(peptide, m)[[1]]
  if (length(mm) >= 2) peptide <- mm[2]
  gsub("\\([^)]*\\)|\\[[^]]*\\]", "", peptide, perl = TRUE)
}

lactyl_positions_from_ascore <- function(value) {
  text <- clean_text(value)
  m <- gregexpr("K(\\d+)\\s*:\\s*(?:Lac|Lactyl(?:ation)?)", text,
                perl = TRUE, ignore.case = TRUE)
  hits <- regmatches(text, m)[[1]]
  if (!length(hits) || !nzchar(hits[1])) return(integer(0))
  positions <- suppressWarnings(as.integer(sub("^K(\\d+).*$", "\\1", hits, perl = TRUE)))
  sort(unique(positions[!is.na(positions)]))
}

lactyl_positions_from_peptide <- function(value) {
  peptide <- clean_text(value)
  positions <- integer(0)
  residue_position <- 0
  index <- 1
  len <- nchar(peptide)
  while (index <= len) {
    char <- substr(peptide, index, index)
    if (char >= "A" && char <= "Z") {
      residue_position <- residue_position + 1
      if (char == "K" && index + 1 <= len &&
          substr(peptide, index + 1, index + 1) %in% c("(", "[")) {
        closer <- if (substr(peptide, index + 1, index + 1) == "(") ")" else "]"
        end <- regexpr(closer, substr(peptide, index + 2, len), fixed = TRUE)
        if (end != -1) {
          end <- end + index + 1
          modification <- substr(peptide, index + 2, end - 1)
          if (grepl("72\\.02|\\bLac\\b|Lactyl", modification,
                    perl = TRUE, ignore.case = TRUE)) {
            positions <- c(positions, residue_position)
          }
          index <- end
        }
      }
    }
    index <- index + 1
  }
  sort(unique(positions))
}

parse_probability_values <- function(value) {
  text <- clean_text(value)
  m <- gregexpr("\\(([-+]?\\d*\\.?\\d+(?:[eE][-+]?\\d+)?)\\)", text, perl = TRUE)
  hits <- regmatches(text, m)[[1]]
  if (!length(hits) || !nzchar(hits[1])) return(numeric(0))
  nums <- suppressWarnings(as.numeric(sub("^\\((.*)\\)$", "\\1", hits, perl = TRUE)))
  nums[!is.na(nums)]
}

apply_annotation_supplement <- function(sites, path) {
  if (!file.exists(path) || !nrow(sites)) return(sites)
  supplement <- read_dtype_str(path, sep = "\t", keep_default_na = FALSE)
  gene_map <- setNames(supplement[["Gene Names (primary)"]], supplement$Entry)
  name_map <- setNames(supplement[["Protein names"]], supplement$Entry)
  missing_gene <- !nzchar(trimws(ifelse(is.na(sites$GeneSymbol), "", sites$GeneSymbol)))
  missing_name <- !nzchar(trimws(ifelse(is.na(sites$ProteinName), "", sites$ProteinName)))
  if (any(missing_gene)) {
    key <- sites$BaseAccession[missing_gene]
    mapped <- gene_map[key]
    sites$GeneSymbol[missing_gene] <- ifelse(is.na(mapped), "", mapped)
  }
  if (any(missing_name)) {
    key <- sites$BaseAccession[missing_name]
    mapped <- name_map[key]
    sites$ProteinName[missing_name] <- ifelse(is.na(mapped), "", mapped)
  }
  sites
}

read_go_annotations <- function(path) {
  raw <- read_dtype_str(path, sep = "\t", file_encoding = "UTF-8-BOM")
  raw$BaseAccession <- vapply(raw[["GENE PRODUCT ID"]],
                              function(v) base_accession(v), character(1))
  taxon <- suppressWarnings(as.numeric(raw[["TAXON ID"]]))
  taxon[is.na(taxon)] <- 9606
  raw$ExcludedNOT <- grepl("(?:^|\\|)NOT(?:\\||$)", raw[["QUALIFIER"]],
                           ignore.case = TRUE, perl = TRUE)
  retained <- raw[taxon == 9606 & !raw$ExcludedNOT & raw$BaseAccession != "", , drop = FALSE]
  accs <- sort(unique(retained$BaseAccession), method = "radix")
  rows <- lapply(accs, function(acc) {
    g <- retained[retained$BaseAccession == acc, , drop = FALSE]
    list(
      BaseAccession = acc,
      GOSymbol = best_annotation(g[["SYMBOL"]]),
      GOTerms = unique_join(g[["GO TERM"]]),
      GONames = unique_join(g[["GO NAME"]]),
      GOEvidenceCodes = unique_join(g[["GO EVIDENCE CODE"]]),
      GOReferences = unique_join(g[["REFERENCE"]]),
      GOAnnotationCount = nrow(g)
    )
  })
  summary_frame <- frame_from_rows(
    rows,
    columns = c(
      "BaseAccession", "GOSymbol", "GOTerms", "GONames", "GOEvidenceCodes",
      "GOReferences", "GOAnnotationCount"
    ),
    numeric_cols = "GOAnnotationCount"
  )
  list(retained = summary_frame, raw = raw)
}

# ---------------------------------------------------------------------------
# I/O helpers (pandas-compatible reading/writing)
# ---------------------------------------------------------------------------

read_dtype_str <- function(path, sep = ",", keep_default_na = TRUE,
                           file_encoding = "") {
  na <- if (keep_default_na) PANDAS_NA else "\u0001NO_NA\u0001"
  args <- list(
    file = path, sep = sep, header = TRUE, check.names = FALSE,
    stringsAsFactors = FALSE, colClasses = "character",
    na.strings = na, quote = "\"", comment.char = "",
    fill = TRUE, blank.lines.skip = FALSE
  )
  if (nzchar(file_encoding)) args$fileEncoding <- file_encoding
  suppressWarnings(do.call(read.table, args))
}

# Python-style repr() for doubles: shortest decimal that round-trips,
# fixed notation when -5 < exponent < 16, else exponential with "e+NN".
#
# R limitations worked around here:
#  - R's sprintf/formatReal is only correctly rounded to 15 significant
#    digits, so digits 16-17 are generated with exact decimal-string
#    arithmetic on the binary value (from sprintf("%a")).
#  - R's as.numeric() reads at most 15 significant digits, so the
#    round-trip test is exact only for candidates with <= 15 digits; for
#    16-17 digit candidates an exact interval test is used instead
#    (parse(D) == v iff D lies in [v - ulp/2, v + ulp/2] with round-half-even
#    boundary rules, checked via cross-multiplied decimal-string integers).
#
# The result reproduces CPython's repr() byte for byte (verified against a
# 12.5k-value corpus and against all floats in the baseline outputs).

# ---- exact non-negative integer arithmetic on decimal strings ----

str_trim_leading_zeros <- function(s) {
  s <- sub("^0+", "", s)
  if (!nzchar(s)) "0" else s
}

str_cmp <- function(a, b) {
  # compare two non-negative integer strings; -1, 0, +1
  a <- str_trim_leading_zeros(a)
  b <- str_trim_leading_zeros(b)
  la <- nchar(a)
  lb <- nchar(b)
  if (la != lb) return(if (la < lb) -1L else 1L)
  if (a == b) return(0L)
  if (a < b) -1L else 1L
}

str_add_one <- function(s) {
  digits <- strsplit(s, "", fixed = TRUE)[[1]]
  carry <- 1L
  for (i in rev(seq_along(digits))) {
    cur <- as.integer(digits[i]) + carry
    if (cur == 10L) {
      digits[i] <- "0"
      carry <- 1L
    } else {
      digits[i] <- as.character(cur)
      carry <- 0L
      break
    }
  }
  if (carry) paste0("1", paste(digits, collapse = "")) else paste(digits, collapse = "")
}

str_sub_one <- function(s) {
  digits <- strsplit(s, "", fixed = TRUE)[[1]]
  carry <- 1L
  for (i in rev(seq_along(digits))) {
    cur <- as.integer(digits[i]) - carry
    if (cur < 0L) {
      digits[i] <- "9"
      carry <- 1L
    } else {
      digits[i] <- as.character(cur)
      carry <- 0L
      break
    }
  }
  paste(digits, collapse = "")
}

str_add <- function(a, b) {
  da <- rev(strsplit(a, "", fixed = TRUE)[[1]])
  db <- rev(strsplit(b, "", fixed = TRUE)[[1]])
  n <- max(length(da), length(db))
  da <- c(da, rep("0", n - length(da)))
  db <- c(db, rep("0", n - length(db)))
  carry <- 0L
  out <- integer(n)
  for (i in seq_len(n)) {
    cur <- as.integer(da[i]) + as.integer(db[i]) + carry
    out[i] <- cur %% 10L
    carry <- cur %/% 10L
  }
  if (carry) out <- c(out, carry)
  str_trim_leading_zeros(paste(rev(out), collapse = ""))
}

str_mul <- function(a, b) {
  # schoolbook multiply of two non-negative decimal strings
  if (a == "0" || b == "0") return("0")
  da <- rev(strsplit(a, "", fixed = TRUE)[[1]])
  db <- rev(strsplit(b, "", fixed = TRUE)[[1]])
  na <- length(da)
  nb <- length(db)
  res <- integer(na + nb)
  for (i in seq_len(na)) {
    ai <- as.integer(da[i])
    if (ai == 0L) next
    carry <- 0L
    for (j in seq_len(nb)) {
      cur <- res[i + j - 1L] + ai * as.integer(db[j]) + carry
      res[i + j - 1L] <- cur %% 10L
      carry <- cur %/% 10L
    }
    res[i + nb] <- res[i + nb] + carry
  }
  for (i in seq_len(length(res) - 1L)) {
    if (res[i] >= 10L) {
      res[i + 1L] <- res[i + 1L] + res[i] %/% 10L
      res[i] <- res[i] %% 10L
    }
  }
  str_trim_leading_zeros(paste(rev(res), collapse = ""))
}

# global precomputed powers of two and five as decimal strings
P2_STR <- vector("list", 1075)
P5_STR <- vector("list", 1075)
P2_STR[[1]] <- "1"
P5_STR[[1]] <- "1"
{
  s2 <- "1"
  s5 <- "1"
  for (k in 1:1074) {
    s2 <- str_add(s2, s2)
    s5 <- str_mul(s5, "5")
    P2_STR[[k + 1]] <- s2
    P5_STR[[k + 1]] <- s5
  }
}

# global precomputed powers of ten as doubles (for exponent estimation)
POW10_DOUBLE <- setNames(10^(-324:309), as.character(-324:309))

# max k with 10^k <= m * 2^e using exact decimal-string comparison
exp10_exact <- function(m_str, e) {
  # log10(m) < 16, so the true exponent is within [floor(e*log10 2),
  # floor(e*log10 2) + 16]; search from one above that range downward.
  k_hi <- floor(e * 0.301029995663981195) + 17L
  m2e <- str_mul(m_str, P2_STR[[max(e, 0L) + 1L]])
  p2ne <- P2_STR[[max(-e, 0L) + 1L]]
  for (k in k_hi:(k_hi - 20L)) {
    lhs <- p2ne
    if (k > 0) lhs <- paste0(lhs, strrep("0", k))
    rhs <- m2e
    if (k < 0) rhs <- paste0(rhs, strrep("0", -k))
    if (str_cmp(lhs, rhs) <= 0) return(k)
  }
  k_hi - 20L
}

# exp10 for the fast path: log10 estimate + exact-double adjustment loop
exp10_fast <- function(av) {
  k <- floor(log10(av))
  while (k < 308 && av >= POW10_DOUBLE[[as.character(k + 1L)]]) k <- k + 1L
  while (k > -323 && av < POW10_DOUBLE[[as.character(k)]]) k <- k - 1L
  # exact-power-of-ten boundary guard: if av is exactly a rounded 10^k
  # double, decide with exact arithmetic (double(10^k) may sit below the
  # true power of ten, in which case Python's exponent is k-1).
  if (av == POW10_DOUBLE[[as.character(k)]]) {
    h <- parse_hex_float(av)
    k <- exp10_exact(h$m, h$e)
  }
  k
}

parse_hex_float <- function(v) {
  # parse sprintf("%a") of a positive double -> list(m, e, m_even)
  # value = m * 2^e with m integer (2^52 + frac for normals, frac for
  # subnormals); m_even = parity of m.
  h <- sprintf("%a", v)
  h <- sub("^-", "", h)
  h <- sub("^0x", "", h)
  parts <- strsplit(h, "p", fixed = TRUE)[[1]]
  mant <- parts[[1]]
  E <- as.numeric(parts[[2]])
  mf <- strsplit(mant, ".", fixed = TRUE)[[1]]
  intpart <- mf[[1]]
  frac <- if (length(mf) > 1) mf[[2]] else ""
  hexdigits <- paste0(intpart, frac)
  hexval <- "0"
  for (ch in strsplit(hexdigits, "", fixed = TRUE)[[1]]) {
    hexval <- str_add(str_mul(hexval, "16"), as.character(strtoi(ch, 16L)))
  }
  m_even <- strtoi(substr(hexdigits, nchar(hexdigits), nchar(hexdigits)), 16L) %% 2L == 0L
  # For normals the leading "1" of %a is the implicit bit, so hexval is the
  # full integer mantissa (2^52 + fraction). For subnormals hexval is the
  # fraction itself. In both cases value = m * 2^(E-52).
  m <- hexval
  list(m = m, e = as.integer(E - 52L), m_even = m_even)
}

# round-half-even division: N / (2^a * 10^b) rounded to nearest integer
round_half_even_div <- function(N_str, a, b) {
  # N = m * 2^max(e,0) * 10^max(s,0); divide by 2^a * 10^b
  r1_any <- FALSE
  r1 <- "0"
  q1 <- N_str
  if (b > 0) {
    L <- nchar(N_str)
    if (L <= b) {
      r1_any <- any(strsplit(N_str, "", fixed = TRUE)[[1]] != "0")
      q1 <- "0"
      r1 <- paste0(strrep("0", b - L), N_str)
    } else {
      q1 <- substr(N_str, 1, L - b)
      r1 <- substr(N_str, L - b + 1, L)
      r1_any <- any(strsplit(r1, "", fixed = TRUE)[[1]] != "0")
    }
  }
  q2 <- q1
  bits <- NULL
  if (a > 0) {
    bits <- integer(a)
    for (i in seq_len(a)) {
      h <- str_halve(q2)
      bits[i] <- h$rem
      q2 <- h$q
    }
  }
  q_odd <- as.integer(substr(q2, nchar(q2), nchar(q2))) %% 2L == 1L
  round_up <- FALSE
  if (a == 0 && b == 0) {
    round_up <- FALSE
  } else if (a == 0) {
    half_str <- paste0("5", strrep("0", b - 1))
    r1_pad <- if (nchar(r1) < b) paste0(strrep("0", b - nchar(r1)), r1) else r1
    cmp <- str_cmp(r1_pad, half_str)
    round_up <- cmp > 0 || (cmp == 0 && q_odd)
  } else {
    top <- bits[a]
    lower_any <- if (a > 1) any(bits[1:(a - 1)] == 1L) else FALSE
    r_gt <- top == 1L && lower_any
    r_gt2 <- top == 1L && !lower_any && r1_any
    r_eq <- top == 1L && !lower_any && !r1_any
    round_up <- r_gt || r_gt2 || (r_eq && q_odd)
  }
  if (round_up) str_add_one(q2) else q2
}

# exact interval test: does the d-digit decimal R_d * 10^(-s) parse
# (round-to-nearest, half-even) to the double m * 2^e?
in_interval <- function(R_d, s, m_str, e, m_even) {
  A <- str_mul(str_mul("2", R_d), P2_STR[[max(-e, 0L) + 1L]])
  if (s < 0) A <- paste0(A, strrep("0", -s))
  m2 <- str_mul("2", m_str)
  Blow <- str_mul(str_sub_one(m2), P2_STR[[max(e, 0L) + 1L]])
  Bhigh <- str_mul(str_add_one(m2), P2_STR[[max(e, 0L) + 1L]])
  if (s > 0) {
    Blow <- paste0(Blow, strrep("0", s))
    Bhigh <- paste0(Bhigh, strrep("0", s))
  }
  cmp_low <- str_cmp(A, Blow)
  cmp_high <- str_cmp(A, Bhigh)
  (cmp_low > 0 || (cmp_low == 0 && m_even)) &&
    (cmp_high < 0 || (cmp_high == 0 && !m_even))
}

format_digits <- function(digits, exp10) {
  if (exp10 >= 16 || exp10 <= -5) {
    mantissa <- if (nchar(digits) == 1) {
      digits
    } else {
      paste0(substr(digits, 1, 1), ".", substr(digits, 2, nchar(digits)))
    }
    paste0(mantissa, sprintf("e%+03d", exp10))
  } else if (exp10 >= 0) {
    if (exp10 + 1 >= nchar(digits)) {
      paste0(digits, strrep("0", exp10 + 1 - nchar(digits)), ".0")
    } else {
      paste0(substr(digits, 1, exp10 + 1), ".",
             substr(digits, exp10 + 2, nchar(digits)))
    }
  } else {
    paste0("0.", strrep("0", -exp10 - 1), digits)
  }
}

py_repr_one <- function(v) {
  if (is.na(v)) return("")
  if (is.nan(v)) return("nan")
  if (is.infinite(v)) return(if (v > 0) "inf" else "-inf")
  if (v == 0) return("0.0")
  sign <- if (v < 0) "-" else ""
  av <- abs(v)
  exp10 <- exp10_fast(av)
  digits <- NULL
  for (d in 1:15) {
    cand <- sprintf("%.*e", d - 1, av)
    if (suppressWarnings(as.numeric(cand)) == av) {
      digits <- gsub("\\.", "", sub("e[+-][0-9]+$", "", cand))
      break
    }
  }
  if (is.null(digits)) {
    h <- parse_hex_float(av)
    exp10 <- exp10_exact(h$m, h$e)
    for (d in 16:17) {
      s <- d - 1L - exp10
      N_str <- str_mul(h$m, P2_STR[[max(h$e, 0L) + 1L]])
      if (s > 0) N_str <- paste0(N_str, strrep("0", s))
      R_d <- round_half_even_div(N_str, max(-h$e, 0L), max(-s, 0L))
      if (in_interval(R_d, s, h$m, h$e, h$m_even)) {
        digits <- R_d
        break
      }
    }
  }
  if (is.null(digits)) digits <- "1" # unreachable: 17 digits always round-trip
  paste0(sign, format_digits(digits, exp10))
}

PY_REPR_CACHE <- new.env(hash = TRUE, parent = emptyenv())

py_repr <- function(x) {
  vapply(x, function(v) {
    key <- sprintf("%a", v)
    if (exists(key, envir = PY_REPR_CACHE, inherits = FALSE)) {
      return(get(key, envir = PY_REPR_CACHE))
    }
    out <- py_repr_one(v)
    assign(key, out, envir = PY_REPR_CACHE)
    out
  }, character(1))
}

# pandas to_csv(index=False, encoding="utf-8-sig") with QUOTE_MINIMAL quoting:
# BOM prefix, quote only fields containing ',' '"' '\r' or '\n', double quotes
# inside quoted fields, NA/NaN written as empty strings, booleans as
# True/False, integers plain, floats via Python-repr.
format_char_field <- function(x) {
  if (is.na(x)) return("")
  if (grepl("[,\"\r\n]", x)) {
    paste0("\"", gsub("\"", "\"\"", x), "\"")
  } else {
    x
  }
}

format_column <- function(col) {
  if (is.logical(col)) {
    return(ifelse(is.na(col), "", ifelse(col, "True", "False")))
  }
  if (is.integer(col)) {
    return(ifelse(is.na(col), "", as.character(col)))
  }
  if (is.numeric(col)) {
    return(py_repr(col))
  }
  vapply(col, format_char_field, character(1))
}

write_csv <- function(frame, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  con <- file(path, open = "wt", encoding = "UTF-8")
  on.exit(close(con))
  writeChar("﻿", con, eos = NULL)
  nms <- names(frame)
  writeLines(paste(nms, collapse = ","), con)
  if (!nrow(frame)) return(invisible(NULL))
  cols <- lapply(frame, format_column)
  n <- nrow(frame)
  lines <- character(n)
  for (i in seq_len(n)) {
    lines[i] <- paste(vapply(cols, function(c) c[[i]], character(1)),
                      collapse = ",")
  }
  writeLines(lines, con)
  invisible(NULL)
}

# Build a data.frame from a list of named rows (all rows share the same keys).
frame_from_rows <- function(rows, columns, numeric_cols = character(0),
                            logical_cols = character(0),
                            integer_cols = character(0)) {
  out <- vector("list", length(columns))
  names(out) <- columns
  for (col in columns) {
    vals <- lapply(rows, function(r) r[[col]])
    if (col %in% numeric_cols) {
      out[[col]] <- as.numeric(unlist(vals))
    } else if (col %in% integer_cols) {
      out[[col]] <- as.integer(unlist(vals))
    } else if (col %in% logical_cols) {
      out[[col]] <- as.logical(unlist(vals))
    } else {
      out[[col]] <- unlist(vals)
    }
    if (is.null(out[[col]])) out[[col]] <- rep(NA, length(rows))
  }
  as.data.frame(out, stringsAsFactors = FALSE)
}

# pandas concat(frames, ignore_index=TRUE, sort=FALSE): union of columns in
# first-appearance order; missing cells NA-typed per column class.
rbind_frames <- function(frames) {
  frames <- frames[vapply(frames, nrow, integer(1)) > 0]
  if (!length(frames)) return(NULL)
  all_cols <- character(0)
  for (fr in frames) {
    for (cn in names(fr)) {
      if (!(cn %in% all_cols)) all_cols <- c(all_cols, cn)
    }
  }
  parts <- lapply(frames, function(fr) {
    out <- lapply(all_cols, function(cn) {
      if (cn %in% names(fr)) return(fr[[cn]])
      if (cn %in% NUMERIC_SITE_COLUMNS) return(rep(NA_real_, nrow(fr)))
      if (cn %in% c(LOGICAL_SITE_COLUMNS, "Sensitivity075Pass",
                    "PXD053_DDA", "PXD053_DIA", "PXD053_S3")) {
        return(rep(NA, nrow(fr)))
      }
      rep(NA_character_, nrow(fr))
    })
    names(out) <- all_cols
    as.data.frame(out, stringsAsFactors = FALSE)
  })
  do.call(rbind, parts)
}

# ---------------------------------------------------------------------------
# extractors.py translation
# ---------------------------------------------------------------------------

exclusion_row <- function(pxd, source, reason, count = 1, detail = "") {
  list(PXD = pxd, Source = source, Reason = reason, Count = count, Detail = detail)
}

log_frame <- function(logs) {
  frame_from_rows(
    logs,
    columns = c("PXD", "Source", "Reason", "Count", "Detail"),
    integer_cols = "Count"
  )
}

dataframe <- function(rows) {
  if (!length(rows)) {
    out <- as.data.frame(
      setNames(lapply(SITE_COLUMNS, function(c) character(0)), SITE_COLUMNS),
      stringsAsFactors = FALSE
    )
    out$LocalizationProb <- numeric(0)
    out$PEP <- numeric(0)
    out$Score <- numeric(0)
    out$DiagnosticIonIntensity <- numeric(0)
    out$ClassI <- logical(0)
    out$PrimaryIncluded <- logical(0)
    out$SourceRow <- integer(0)
    return(out)
  }
  frame <- frame_from_rows(
    rows, columns = SITE_COLUMNS,
    numeric_cols = NUMERIC_SITE_COLUMNS,
    logical_cols = LOGICAL_SITE_COLUMNS,
    integer_cols = "SourceRow"
  )
  for (col in SITE_COLUMNS) {
    if (!(col %in% names(frame))) frame[[col]] <- ""
  }
  frame
}

valid_maxquant_site_ids <- function(path) {
  table <- read_dtype_str(path, sep = "\t")
  keep <- rep(TRUE, nrow(table))
  for (cn in c("Reverse", "Potential contaminant", "Contaminant")) {
    if (cn %in% names(table)) keep <- keep & !is_true(table[[cn]])
  }
  n_la <- if ("Lactyl (K)" %in% names(table)) {
    vapply(table[["Lactyl (K)"]], number, numeric(1))
  } else {
    rep(NA_real_, nrow(table))
  }
  keep <- keep & !is.na(n_la) & n_la > 0
  ids_col <- if ("Lactyl (K) site IDs" %in% names(table)) {
    table[["Lactyl (K) site IDs"]]
  } else {
    rep("", nrow(table))
  }
  valid <- unique(unlist(lapply(ids_col[keep], split_tokens)))
  list(ids = valid, rows = list())
}

extract_pxd014870 <- function(data_root, project_root, doi) {
  root <- file.path(data_root, "PXD014870", "search_results")
  source_map <- list(
    "MCF7_DCA_SILAC_Kla_IP" = c("MCF7_DCA", "DCA"),
    "MCF7_Hypoxia_SILAC_Kla_IP" = c("MCF7_Hypoxia", "Hypoxia"),
    "MCF7_Rotenone_SILAC_Kla_IP" = c("MCF7_Rotenone", "Rotenone"),
    "MCF7_U_13C6_Glucose_SILAC_Kla_IP" = c("MCF7_U13C6_Glucose", "U-13C6 glucose")
  )
  rows <- list()
  logs <- list()
  for (directory_name in names(source_map)) {
    source_token <- source_map[[directory_name]][1]
    group <- source_map[[directory_name]][2]
    txt <- file.path(root, directory_name, "txt")
    site_path <- file.path(txt, "Lactyl (K)Sites.txt")
    peptide_path <- file.path(txt, "modificationSpecificPeptides.txt")
    if (!file.exists(site_path) || !file.exists(peptide_path)) {
      logs[[length(logs) + 1]] <- exclusion_row(
        "PXD014870", directory_name, "missing_required_search_table"
      )
      next
    }
    valid_ids <- valid_maxquant_site_ids(peptide_path)$ids
    sites <- read_dtype_str(site_path, sep = "\t")
    n <- nrow(sites)
    bad_cols <- rep(FALSE, n)
    for (cn in c("Reverse", "Potential contaminant", "Contaminant")) {
      if (cn %in% names(sites)) bad_cols <- bad_cols | is_true(sites[[cn]])
    }
    id_col <- if ("id" %in% names(sites)) sites$id else rep("", n)
    amino_col <- if ("Amino acid" %in% names(sites)) sites$`Amino acid` else rep("K", n)
    protein_col <- if ("Protein" %in% names(sites)) {
      sites$Protein
    } else if ("Leading proteins" %in% names(sites)) {
      sites$`Leading proteins`
    } else {
      rep("", n)
    }
    position_col <- if ("Position" %in% names(sites)) {
      sites$Position
    } else if ("Positions within proteins" %in% names(sites)) {
      sites$`Positions within proteins`
    } else {
      rep("", n)
    }
    fasta_col <- if ("Fasta headers" %in% names(sites)) sites$`Fasta headers` else rep("", n)
    gene_col <- if ("Gene names" %in% names(sites)) sites$`Gene names` else rep("", n)
    pname_col <- if ("Protein names" %in% names(sites)) sites$`Protein names` else rep("", n)
    loc_col <- if ("Localization prob" %in% names(sites)) sites$`Localization prob` else rep("", n)
    rawfile_col <- if ("Best localization raw file" %in% names(sites)) {
      sites$`Best localization raw file`
    } else {
      rep("", n)
    }
    diag_col <- if ("Diagnostic peak" %in% names(sites)) sites$`Diagnostic peak` else rep("", n)
    modseq_col <- if ("Modified sequence" %in% names(sites)) sites$`Modified sequence` else rep("", n)
    pep_col <- if ("PEP" %in% names(sites)) sites$PEP else rep("", n)
    score_col <- if ("Score" %in% names(sites)) sites$Score else rep("", n)
    reason_counts <- list()
    for (i in seq_len(n)) {
      reason <- ""
      if (bad_cols[i]) reason <- "reverse_or_contaminant"
      site_id <- clean_text(id_col[i])
      if (!nzchar(reason) &&
          (!nzchar(site_id) || !(site_id %in% valid_ids))) {
        reason <- "no_valid_kla_site_id_in_modificationSpecificPeptides"
      }
      if (!nzchar(reason) && toupper(clean_text(amino_col[i])) != "K") {
        reason <- "non_lysine_site"
      }
      accession <- normalize_accession(protein_col[i])
      position <- integer_value(position_col[i])
      if (!nzchar(reason) && (!nzchar(accession) || is.na(position))) {
        reason <- "missing_accession_or_kla_position"
      }
      if (nzchar(reason)) {
        reason_counts[[reason]] <- (reason_counts[[reason]] %||% 0L) + 1L
        next
      }
      ann <- annotation_from_fasta(fasta_col[i], accession)
      gene <- ann[1]
      protein_name <- ann[2]
      if (!nzchar(gene)) {
        gt <- split_tokens(gene_col[i])
        gene <- if (length(gt)) gt[1] else ""
      }
      if (!nzchar(protein_name)) {
        pt <- split_tokens(pname_col[i])
        protein_name <- if (length(pt)) pt[1] else ""
      }
      localization <- number(loc_col[i])
      raw_file <- clean_text(rawfile_col[i])
      if (!nzchar(raw_file)) raw_file <- source_token
      diagnostic <- clean_text(diag_col[i])
      rows[[length(rows) + 1]] <- blank_site(
        PXD = "PXD014870", DOI = doi, SampleName = raw_file,
        CellOrTissueType = "MCF7", ExperimentalGroup = group,
        Replicate = raw_file, EnrichmentStatus = "Enriched",
        AcquisitionMode = "DDA", Accession = accession,
        BaseAccession = base_accession(accession), GeneSymbol = gene,
        ProteinName = protein_name, KlaSite = paste0("K", position),
        ModifiedPeptide = clean_text(modseq_col[i]),
        LocalizationProb = localization, PEP = number(pep_col[i]),
        Score = number(score_col[i]), DiagnosticIon = diagnostic,
        ClassI = !is.na(localization) && localization >= 0.75,
        SourceFile = relative_path(site_path, project_root),
        SourceRow = i + 1L, SiteID = site_id,
        EvidenceMode = "maxquant_site_table_plus_modified_peptides",
        SourceConfidence = if (!is.na(localization) && localization >= 0.75) {
          "high_localized"
        } else {
          "standard_author_search_result"
        }
      )
    }
    for (reason in sort(names(reason_counts), method = "radix")) {
      logs[[length(logs) + 1]] <- exclusion_row(
        "PXD014870", relative_path(site_path, project_root),
        reason, reason_counts[[reason]]
      )
    }
  }
  frame <- dataframe(rows)
  if (nrow(frame) > 0) {
    lp <- suppressWarnings(as.numeric(frame$LocalizationProb))
    frame$Sensitivity075Pass <- !is.na(lp) & lp >= 0.75
  }
  list(frame, log_frame(logs))
}

peaks_annotations <- function(proteins) {
  annotations <- list()
  acc_col <- if ("Accession" %in% names(proteins)) proteins$Accession else rep("", nrow(proteins))
  desc_col <- if ("Description" %in% names(proteins)) {
    proteins$Description
  } else {
    rep("", nrow(proteins))
  }
  for (i in seq_len(nrow(proteins))) {
    accession <- normalize_accession(acc_col[i])
    ann <- annotation_from_description(desc_col[i])
    if (nzchar(accession) && is.null(annotations[[accession]])) {
      annotations[[accession]] <- ann
    }
  }
  annotations
}

peaks_start_map <- function(protein_peptides) {
  mapping <- list()
  acc_col <- if ("Protein Accession" %in% names(protein_peptides)) {
    protein_peptides$`Protein Accession`
  } else {
    rep("", nrow(protein_peptides))
  }
  pep_col <- if ("Peptide" %in% names(protein_peptides)) {
    protein_peptides$Peptide
  } else {
    rep("", nrow(protein_peptides))
  }
  start_col <- if ("Start" %in% names(protein_peptides)) {
    protein_peptides$Start
  } else {
    rep("", nrow(protein_peptides))
  }
  for (i in seq_len(nrow(protein_peptides))) {
    accession <- normalize_accession(acc_col[i])
    peptide <- strip_peptide_modifications(pep_col[i])
    start <- integer_value(start_col[i])
    if (nzchar(accession) && nzchar(peptide) && !is.na(start)) {
      mapping[[peptide]] <- c(mapping[[peptide]], list(list(accession, start)))
    }
  }
  mapping
}

matching_accession_mappings <- function(matches, raw_accession) {
  candidates <- unique(vapply(
    accession_candidates(raw_accession),
    function(a) base_accession(a), character(1)
  ))
  if (!length(candidates)) return(matches)
  Filter(function(m) base_accession(m[[1]]) %in% candidates, matches)
}

site_position_pairs <- function(site_ids, positions) {
  if (length(site_ids) != length(positions)) return(list())
  Map(list, site_ids, positions)
}

normalize_scan <- function(value) {
  text <- clean_text(value)
  sub("^.*:", "", text)
}

read_marker_156 <- function(path) {
  markers <- list()
  if (!file.exists(path)) return(markers)
  table <- suppressWarnings(read.csv(
    path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE,
    na.strings = PANDAS_NA, comment.char = "", fill = TRUE,
    blank.lines.skip = FALSE
  ))
  nms <- names(table)
  required <- c("Peptide", "scan", "Source File", "ion", "theo m/z",
                "ion relative intensity(%)")
  if (!all(required %in% nms)) return(markers)
  for (i in seq_len(nrow(table))) {
    ion <- clean_text(table$ion[i])
    if (tolower(ion) != "marker(156)") next
    theoretical <- number(table$`theo m/z`[i])
    if (is.na(theoretical) || abs(theoretical - 156.1025) > 0.01) next
    key <- paste(
      clean_text(table$Peptide[i]),
      basename(clean_text(table$`Source File`[i])),
      normalize_scan(table$scan[i]),
      sep = "\u0001"
    )
    intensity <- number(table$`ion relative intensity(%)`[i])
    if (is.na(intensity)) intensity <- 0
    markers[[key]] <- max(markers[[key]] %||% 0, intensity)
  }
  markers
}

pxd028_directory_metadata <- function(directory) {
  name <- basename(directory)
  parent <- basename(dirname(directory))
  enriched <- startsWith(parent, "Enrichment")
  species <- if (startsWith(name, "BV2") || startsWith(name, "RAW-")) "mouse" else "human"
  cell_type <- if (startsWith(name, "HEK293T")) {
    "HEK293T"
  } else if (startsWith(name, "HCT116")) {
    "HCT116"
  } else if (startsWith(name, "TALL")) {
    "T-ALL"
  } else if (startsWith(name, "BV2")) {
    "BV2"
  } else {
    "RAW264.7"
  }
  aggregate <- grepl("all HCD", name, fixed = TRUE)
  hcd <- regexpr("HCD[0-9]+", name, perl = TRUE)
  hcd_label <- if (hcd != -1) {
    regmatches(name, hcd)
  } else if (aggregate) {
    "all HCD aggregate"
  } else {
    "not specified"
  }
  list(
    Directory = name,
    SearchCollection = parent,
    Species = species,
    CellType = cell_type,
    EnrichmentStatus = if (enriched) "Enriched" else "Unenriched",
    HCDCondition = hcd_label,
    AggregateDirectory = aggregate
  )
}

extract_pxd028488 <- function(data_root, project_root, doi) {
  root <- file.path(data_root, "PXD028488", "search_results")
  rows <- list()
  logs <- list()
  directory_audit <- list()
  for (collection in c("Enrichment-Search files", "Nonenrichment-Search files")) {
    dirs <- sort(
      list.files(file.path(root, collection), pattern = "-Search files$",
                 full.names = TRUE, sort = FALSE),
      method = "radix"
    )
    for (directory in dirs) {
      meta <- pxd028_directory_metadata(directory)
      psm_path <- file.path(directory, "DB search psm.csv")
      meta$PSMAvailable <- file.exists(psm_path)
      if (meta$Species != "human") {
        meta$PrimaryStatus <- "excluded"
        meta$Reason <- "non_human_cell_line"
        meta$KlaEvidenceRows <- 0L
        directory_audit[[length(directory_audit) + 1]] <- meta
        logs[[length(logs) + 1]] <- exclusion_row(
          "PXD028488", relative_path(directory, project_root),
          "non_human_directory", detail = meta$CellType
        )
        next
      }
      if (!file.exists(psm_path)) {
        meta$PrimaryStatus <- "excluded"
        meta$Reason <- "missing_DB_search_psm.csv"
        meta$KlaEvidenceRows <- 0L
        directory_audit[[length(directory_audit) + 1]] <- meta
        logs[[length(logs) + 1]] <- exclusion_row(
          "PXD028488", relative_path(directory, project_root),
          "missing_required_psm_table"
        )
        next
      }
      validation_only <- isTRUE(meta$AggregateDirectory)
      meta$PrimaryStatus <- if (validation_only) "validation_only" else "included"
      meta$Reason <- if (validation_only) {
        "aggregate_of_individual_HCD_runs"
      } else {
        "human_analyzable_search_directory"
      }

      psm <- read_dtype_str(psm_path)
      pp_path <- file.path(directory, "protein-peptides.csv")
      proteins_path <- file.path(directory, "proteins.csv")
      protein_peptides <- read_dtype_str(pp_path)
      proteins <- read_dtype_str(proteins_path)
      annotations <- peaks_annotations(proteins)
      starts <- peaks_start_map(protein_peptides)
      markers <- read_marker_156(file.path(directory, "PSM ions.csv"))

      counters <- c(
        score_below_20 = 0L, no_kla_modification = 0L, no_position_mapping = 0L,
        accession_mapping_mismatch = 0L, reverse_or_contaminant = 0L
      )

      acc_col <- if ("Accession" %in% names(psm)) psm$Accession else rep("", nrow(psm))
      score_col <- if ("-10lgP" %in% names(psm)) psm$`-10lgP` else rep("", nrow(psm))
      pep_col <- if ("Peptide" %in% names(psm)) psm$Peptide else rep("", nrow(psm))
      ascore_col <- if ("AScore" %in% names(psm)) psm$AScore else rep("", nrow(psm))
      ptm_col <- if ("PTM" %in% names(psm)) psm$PTM else rep("", nrow(psm))
      src_col <- if ("Source File" %in% names(psm)) psm$`Source File` else rep("", nrow(psm))
      scan_col <- if ("Scan" %in% names(psm)) psm$Scan else rep("", nrow(psm))

      before <- length(rows)
      for (i in seq_len(nrow(psm))) {
        raw_accession <- clean_text(acc_col[i])
        if (grepl("DECOY|REV|REVERSED|CON_|CONTAM|CRAP", raw_accession,
                  ignore.case = TRUE, perl = TRUE)) {
          counters["reverse_or_contaminant"] <- counters["reverse_or_contaminant"] + 1L
          next
        }
        score <- number(score_col[i])
        if (is.na(score) || score < 20) {
          counters["score_below_20"] <- counters["score_below_20"] + 1L
          next
        }
        peptide_raw <- clean_text(pep_col[i])
        positions <- lactyl_positions_from_ascore(ascore_col[i])
        if (!length(positions)) {
          positions <- lactyl_positions_from_peptide(peptide_raw)
        }
        ptm_has_lac <- grepl("(?:^|;)\\s*(?:Lac|Lactyl(?:ation)?)\\b",
                             clean_text(ptm_col[i]), ignore.case = TRUE, perl = TRUE)
        if (!length(positions) && !ptm_has_lac) {
          counters["no_kla_modification"] <- counters["no_kla_modification"] + 1L
          next
        }
        if (!length(positions)) {
          counters["no_position_mapping"] <- counters["no_position_mapping"] + 1L
          next
        }
        peptide_clean <- strip_peptide_modifications(peptide_raw)
        matches <- starts[[peptide_clean]] %||% list()
        mapped_matches <- matching_accession_mappings(matches, raw_accession)
        if (length(matches) && !length(mapped_matches)) {
          counters["accession_mapping_mismatch"] <- counters["accession_mapping_mismatch"] + 1L
        }
        matches <- mapped_matches
        if (!length(matches)) {
          counters["no_position_mapping"] <- counters["no_position_mapping"] + 1L
          next
        }
        source_file <- basename(clean_text(src_col[i]))
        marker_key <- paste(
          peptide_raw, source_file, normalize_scan(scan_col[i]),
          sep = "\u0001"
        )
        marker_intensity <- markers[[marker_key]]
        for (m in matches) {
          accession <- m[[1]]
          start <- m[[2]]
          ann <- annotations[[accession]] %||% c("", "")
          for (peptide_position in positions) {
            protein_position <- start + peptide_position - 1
            rows[[length(rows) + 1]] <- blank_site(
              PXD = "PXD028488", DOI = doi,
              SampleName = if (nzchar(source_file)) source_file else meta$Directory,
              CellOrTissueType = meta$CellType,
              ExperimentalGroup = paste(meta$EnrichmentStatus, meta$HCDCondition),
              Replicate = source_file,
              EnrichmentStatus = meta$EnrichmentStatus,
              AcquisitionMode = "DDA",
              Accession = accession,
              BaseAccession = base_accession(accession),
              GeneSymbol = ann[1],
              ProteinName = ann[2],
              KlaSite = paste0("K", protein_position),
              ModifiedPeptide = peptide_raw,
              Score = score,
              DiagnosticIon = if (!is.null(marker_intensity)) {
                "CycIm m/z 156.1025"
              } else {
                "not detected in exported ion assignment"
              },
              DiagnosticIonIntensity = marker_intensity,
              SourceFile = relative_path(psm_path, project_root),
              SourceRow = i + 1L,
              SiteID = paste0(
                meta$Directory, ":", i + 1L, ":", accession,
                ":K", protein_position
              ),
              EvidenceMode = "PEAKS_PSM_plus_protein_peptides",
              SourceConfidence = if (!is.null(marker_intensity)) {
                "high_diagnostic_supported"
              } else {
                "standard_search_result"
              },
              PrimaryIncluded = !validation_only,
              InclusionStatus = if (validation_only) "validation_only" else "included",
              ExclusionReason = if (validation_only) "aggregate_HCD_duplicate" else ""
            )
          }
        }
      }
      evidence_count <- length(rows) - before
      meta$KlaEvidenceRows <- as.integer(evidence_count)
      meta$DiagnosticSupportedRows <- if (evidence_count > 0) {
        sum(vapply(
          rows[(before + 1):length(rows)],
          function(r) !is.na(r$DiagnosticIonIntensity), logical(1)
        ))
      } else {
        0L
      }
      directory_audit[[length(directory_audit) + 1]] <- meta
      for (reason in names(counters)) {
        if (counters[[reason]]) {
          logs[[length(logs) + 1]] <- exclusion_row(
            "PXD028488", relative_path(psm_path, project_root),
            reason, counters[[reason]]
          )
        }
      }
    }
  }
  audit_frame <- frame_from_rows(
    directory_audit,
    columns = c(
      "Directory", "SearchCollection", "Species", "CellType",
      "EnrichmentStatus", "HCDCondition", "AggregateDirectory", "PSMAvailable",
      "PrimaryStatus", "Reason", "KlaEvidenceRows", "DiagnosticSupportedRows"
    ),
    logical_cols = c("AggregateDirectory", "PSMAvailable"),
    integer_cols = "KlaEvidenceRows",
    numeric_cols = "DiagnosticSupportedRows"
  )
  list(dataframe(rows), log_frame(logs), audit_frame)
}

extract_pxd050470 <- function(data_root, project_root, doi) {
  root <- file.path(data_root, "PXD050470", "supplementary")
  sources <- list(
    list("prca2331-sup-0005-tables3.xlsx", 12L, "Proteins accession",
         "Positions within proteins", "Gene names", "Table S3"),
    list("prca2331-sup-0014-tables12.xlsx", 10L, "Proteins accession",
         "Positions within proteins", "Gene name", "Table S12")
  )
  rows <- list()
  logs <- list()
  for (src in sources) {
    filename <- src[[1]]
    header <- src[[2]]
    accession_col <- src[[3]]
    position_col <- src[[4]]
    gene_col <- src[[5]]
    table_name <- src[[6]]
    path <- file.path(root, filename)
    if (!file.exists(path)) {
      logs[[length(logs) + 1]] <- exclusion_row(
        "PXD050470", filename, "missing_supplementary_table"
      )
      next
    }
    table <- suppressWarnings(read_excel(
      path, sheet = "Sheet1", skip = header, col_names = TRUE,
      col_types = "text"
    ))
    names(table) <- vapply(names(table), clean_text, character(1))
    n <- nrow(table)
    intensity_columns <- names(table)[vapply(
      names(table),
      function(cn) startsWith(clean_text(cn), "Intensity_"), logical(1)
    )]
    acc_vec <- if (accession_col %in% names(table)) {
      table[[accession_col]]
    } else {
      rep("", n)
    }
    pos_vec <- if (position_col %in% names(table)) {
      table[[position_col]]
    } else {
      rep("", n)
    }
    gene_vec <- if (gene_col %in% names(table)) {
      table[[gene_col]]
    } else {
      rep("", n)
    }
    loc_vec <- if ("Localization probability" %in% names(table)) {
      table[["Localization probability"]]
    } else {
      rep("", n)
    }
    modseq_vec <- if ("Modified sequence" %in% names(table)) {
      table[["Modified sequence"]]
    } else {
      rep("", n)
    }
    pep_vec <- if ("PEP" %in% names(table)) table$PEP else rep("", n)
    score_vec <- if ("Score" %in% names(table)) table$Score else rep("", n)
    for (i in seq_len(n)) {
      accession <- normalize_accession(acc_vec[i])
      position <- integer_value(pos_vec[i])
      if (!nzchar(accession) || is.na(position)) next
      sample_values <- intensity_columns[vapply(
        intensity_columns,
        function(cn) nzchar(clean_text(table[[cn]][i])), logical(1)
      )]
      if (!length(sample_values)) sample_values <- "author_table_aggregate"
      for (sample_column in sample_values) {
        sample <- if (sample_column == "author_table_aggregate") {
          sample_column
        } else {
          sub("^Intensity_", "", sample_column)
        }
        localization <- number(loc_vec[i])
        rows[[length(rows) + 1]] <- blank_site(
          PXD = "PXD050470", DOI = doi, SampleName = sample,
          CellOrTissueType = "Human hippocampus",
          ExperimentalGroup = "Physiological hippocampus",
          Replicate = sample, EnrichmentStatus = "Enriched",
          AcquisitionMode = "DDA", Accession = accession,
          BaseAccession = base_accession(accession),
          GeneSymbol = clean_text(gene_vec[i]), ProteinName = "",
          KlaSite = paste0("K", position),
          ModifiedPeptide = clean_text(modseq_vec[i]),
          LocalizationProb = localization, PEP = number(pep_vec[i]),
          Score = number(score_vec[i]), ClassI = TRUE,
          SourceFile = relative_path(path, project_root),
          SourceRow = i + header + 1L,
          SiteID = paste0(table_name, ":", i + header + 1L, ":",
                          accession, ":K", position),
          EvidenceMode = "author_supplementary_table",
          SourceConfidence = "high_author_class_I"
        )
      }
    }
  }
  raw_dir <- file.path(data_root, "PXD050470", "raw")
  if (length(list.files(raw_dir, all.files = TRUE, no.. = TRUE)) == 0) {
    logs[[length(logs) + 1]] <- exclusion_row(
      "PXD050470", relative_path(raw_dir, project_root),
      "repository_raw_files_not_present_locally", 27L,
      "Raw files were not re-searched; S3 and S12 are the analysis evidence."
    )
  }
  list(dataframe(rows), log_frame(logs))
}

extract_pxd053474_dda <- function(data_root, project_root, doi) {
  root <- file.path(data_root, "PXD053474", "search_results", "extracted")
  rows <- list()
  logs <- list()
  paths <- sort(
    list.files(root, pattern = "^peptide\\.csv$", recursive = TRUE,
               full.names = TRUE, sort = FALSE),
    method = "radix"
  )
  for (peptide_path in paths) {
    directory <- dirname(peptide_path)
    pp_path <- file.path(directory, "protein-peptides.csv")
    proteins_path <- file.path(directory, "proteins.csv")
    if (!file.exists(pp_path) || !file.exists(proteins_path)) {
      logs[[length(logs) + 1]] <- exclusion_row(
        "PXD053474", relative_path(directory, project_root),
        "missing_DDA_mapping_table"
      )
      next
    }
    peptide_table <- read_dtype_str(peptide_path)
    protein_peptides <- read_dtype_str(pp_path)
    proteins <- read_dtype_str(proteins_path)
    starts <- peaks_start_map(protein_peptides)
    annotations <- peaks_annotations(proteins)
    enriched <- grepl("Enriched-DDA", peptide_path, fixed = TRUE)
    compartment <- sub("-Search files$", "", basename(directory))
    counters <- c(
      score_below_20 = 0L, no_kla_modification = 0L, no_position_mapping = 0L,
      accession_mapping_mismatch = 0L
    )
    acc_col <- if ("Accession" %in% names(peptide_table)) {
      peptide_table$Accession
    } else {
      rep("", nrow(peptide_table))
    }
    score_col <- if ("-10lgP" %in% names(peptide_table)) {
      peptide_table$`-10lgP`
    } else {
      rep("", nrow(peptide_table))
    }
    pep_col <- if ("Peptide" %in% names(peptide_table)) {
      peptide_table$Peptide
    } else {
      rep("", nrow(peptide_table))
    }
    ascore_col <- if ("AScore" %in% names(peptide_table)) {
      peptide_table$AScore
    } else {
      rep("", nrow(peptide_table))
    }
    ptm_col <- if ("PTM" %in% names(peptide_table)) {
      peptide_table$PTM
    } else {
      rep("", nrow(peptide_table))
    }
    src_col <- if ("Source File" %in% names(peptide_table)) {
      peptide_table$`Source File`
    } else {
      rep("", nrow(peptide_table))
    }
    for (i in seq_len(nrow(peptide_table))) {
      score <- number(score_col[i])
      if (is.na(score) || score < 20) {
        counters["score_below_20"] <- counters["score_below_20"] + 1L
        next
      }
      peptide_raw <- clean_text(pep_col[i])
      positions <- lactyl_positions_from_ascore(ascore_col[i])
      if (!length(positions)) {
        positions <- lactyl_positions_from_peptide(peptide_raw)
      }
      ptm_has_lac <- grepl("(?:^|;)\\s*(?:Lac|Lactyl(?:ation)?)\\b",
                           clean_text(ptm_col[i]), ignore.case = TRUE, perl = TRUE)
      if (!length(positions) && !ptm_has_lac) {
        counters["no_kla_modification"] <- counters["no_kla_modification"] + 1L
        next
      }
      if (!length(positions)) {
        counters["no_position_mapping"] <- counters["no_position_mapping"] + 1L
        next
      }
      peptide_clean <- strip_peptide_modifications(peptide_raw)
      matches <- starts[[peptide_clean]] %||% list()
      mapped_matches <- matching_accession_mappings(matches, acc_col[i])
      if (length(matches) && !length(mapped_matches)) {
        counters["accession_mapping_mismatch"] <- counters["accession_mapping_mismatch"] + 1L
      }
      matches <- mapped_matches
      if (!length(matches)) {
        counters["no_position_mapping"] <- counters["no_position_mapping"] + 1L
        next
      }
      sample <- basename(clean_text(src_col[i]))
      if (!nzchar(sample)) sample <- compartment
      for (m in matches) {
        accession <- m[[1]]
        start <- m[[2]]
        ann <- annotations[[accession]] %||% c("", "")
        for (peptide_position in positions) {
          position <- start + peptide_position - 1
          rows[[length(rows) + 1]] <- blank_site(
            PXD = "PXD053474", DOI = doi, SampleName = sample,
            CellOrTissueType = "HCT116", ExperimentalGroup = compartment,
            Replicate = sample,
            EnrichmentStatus = if (enriched) "Enriched" else "Unenriched",
            AcquisitionMode = "DDA", Accession = accession,
            BaseAccession = base_accession(accession), GeneSymbol = ann[1],
            ProteinName = ann[2], KlaSite = paste0("K", position),
            ModifiedPeptide = peptide_raw, Score = score,
            DiagnosticIon = "not available in deposited PXD053474 DDA tables",
            SourceFile = relative_path(peptide_path, project_root),
            SourceRow = i + 1L,
            SiteID = paste0("DDA:", basename(directory), ":", i + 1L, ":",
                            accession, ":K", position),
            EvidenceMode = "PEAKS_DDA_peptide_plus_protein_peptides",
            SourceConfidence = "standard_search_result"
          )
        }
      }
    }
    for (reason in names(counters)) {
      if (counters[[reason]]) {
        logs[[length(logs) + 1]] <- exclusion_row(
          "PXD053474", relative_path(peptide_path, project_root),
          reason, counters[[reason]]
        )
      }
    }
  }
  list(dataframe(rows), log_frame(logs))
}

run_prefixes <- function(columns) {
  prefixes <- list()
  suffix <- ".PG.IsIdentified"
  for (column in columns) {
    if (!endsWith(column, suffix)) next
    prefix <- substr(column, 1, nchar(column) - nchar(suffix))
    sample <- sub("^\\[\\d+\\]\\s*", "", prefix, perl = TRUE)
    prefixes[[length(prefixes) + 1]] <- c(prefix, sample)
  }
  prefixes
}

extract_pxd053474_dia <- function(data_root, project_root, doi) {
  root <- file.path(data_root, "PXD053474", "search_results", "extracted")
  rows <- list()
  logs <- list()
  paths <- sort(
    list.files(root, pattern = "ptm-site.*Report\\.tsv$", recursive = TRUE,
               full.names = TRUE, sort = FALSE),
    method = "radix"
  )
  for (path in paths) {
    table <- read_dtype_str(path, sep = "\t")
    n <- nrow(table)
    prefixes <- run_prefixes(names(table))
    enriched <- grepl("Enriched-DIA", path, fixed = TRUE)
    compartment <- basename(dirname(path))
    title_col <- if ("PTM.ModificationTitle" %in% names(table)) {
      table[["PTM.ModificationTitle"]]
    } else {
      rep("", n)
    }
    aa_col <- if ("PTM.SiteAA" %in% names(table)) table[["PTM.SiteAA"]] else rep("", n)
    org_col <- if ("PG.Organisms" %in% names(table)) {
      table[["PG.Organisms"]]
    } else {
      rep("", n)
    }
    protid_col <- if ("PTM.ProteinId" %in% names(table)) {
      table[["PTM.ProteinId"]]
    } else {
      rep("", n)
    }
    loc_col <- if ("PTM.SiteLocation" %in% names(table)) {
      table[["PTM.SiteLocation"]]
    } else {
      rep("", n)
    }
    genes_col <- if ("PG.Genes" %in% names(table)) table[["PG.Genes"]] else rep("", n)
    desc_col <- if ("PG.ProteinDescriptions" %in% names(table)) {
      table[["PG.ProteinDescriptions"]]
    } else {
      rep("", n)
    }
    flank_col <- if ("PTM.FlankingRegion" %in% names(table)) {
      table[["PTM.FlankingRegion"]]
    } else {
      rep("", n)
    }
    ckey_col <- if ("PTM.CollapseKey" %in% names(table)) {
      table[["PTM.CollapseKey"]]
    } else {
      rep("", n)
    }
    pg_pep_col <- if ("PG.PEP" %in% names(table)) table[["PG.PEP"]] else rep("", n)
    pg_cscore_col <- if ("PG.Cscore" %in% names(table)) {
      table[["PG.Cscore"]]
    } else {
      rep("", n)
    }
    prefix_isid <- lapply(prefixes, function(p) table[[paste0(p[1], ".PG.IsIdentified")]])
    prefix_loc <- lapply(prefixes, function(p) {
      cn <- paste0(p[1], ".PTM.SiteProbability")
      if (cn %in% names(table)) table[[cn]] else rep("", n)
    })
    prefix_group <- lapply(prefixes, function(p) {
      cn <- paste0(p[1], ".PTM.Group")
      if (cn %in% names(table)) table[[cn]] else rep("", n)
    })
    prefix_pep <- lapply(prefixes, function(p) {
      cn <- paste0(p[1], ".PG.PEP (Run-Wise)")
      if (cn %in% names(table)) table[[cn]] else rep("", n)
    })
    prefix_score <- lapply(prefixes, function(p) {
      cn <- paste0(p[1], ".PG.Cscore (Run-Wise)")
      if (cn %in% names(table)) table[[cn]] else rep("", n)
    })
    no_lac <- 0L
    for (i in seq_len(n)) {
      if (tolower(clean_text(title_col[i])) != "lac" ||
          toupper(clean_text(aa_col[i])) != "K") {
        no_lac <- no_lac + 1L
        next
      }
      org <- clean_text(org_col[i])
      if (nzchar(org) && !grepl("Homo sapiens", org, fixed = TRUE)) next
      accession <- normalize_accession(protid_col[i])
      position <- integer_value(loc_col[i])
      if (!nzchar(accession) || is.na(position)) next
      genes <- split_tokens(genes_col[i])
      names <- split_tokens(desc_col[i])
      emitted <- FALSE
      for (k in seq_along(prefixes)) {
        if (!is_true(prefix_isid[[k]][i])) next
        emitted <- TRUE
        localization <- number(prefix_loc[[k]][i])
        group_value <- clean_text(prefix_group[[k]][i])
        if (!nzchar(group_value)) group_value <- clean_text(flank_col[i])
        site_id <- clean_text(ckey_col[i])
        if (!nzchar(site_id)) {
          site_id <- paste0("DIA:", basename(dirname(path)), ":", i + 1L)
        }
        rows[[length(rows) + 1]] <- blank_site(
          PXD = "PXD053474", DOI = doi, SampleName = prefixes[[k]][2],
          CellOrTissueType = "HCT116", ExperimentalGroup = compartment,
          Replicate = prefixes[[k]][2],
          EnrichmentStatus = if (enriched) "Enriched" else "Unenriched",
          AcquisitionMode = "DIA", Accession = accession,
          BaseAccession = base_accession(accession),
          GeneSymbol = if (length(genes)) genes[1] else "",
          ProteinName = if (length(names)) names[1] else "",
          KlaSite = paste0("K", position),
          ModifiedPeptide = group_value,
          LocalizationProb = localization,
          PEP = number(prefix_pep[[k]][i]),
          Score = number(prefix_score[[k]][i]),
          ClassI = !is.na(localization) && localization >= 0.75,
          SourceFile = relative_path(path, project_root),
          SourceRow = i + 1L,
          SiteID = site_id,
          EvidenceMode = "DIA_ptm_site_report",
          SourceConfidence = if (!is.na(localization) && localization >= 0.75) {
            "high_localized"
          } else {
            "standard_search_result"
          }
        )
      }
      if (!emitted) {
        site_id <- clean_text(ckey_col[i])
        if (!nzchar(site_id)) {
          site_id <- paste0("DIA:", basename(dirname(path)), ":", i + 1L)
        }
        rows[[length(rows) + 1]] <- blank_site(
          PXD = "PXD053474", DOI = doi,
          SampleName = paste0(compartment, "_experiment_wide"),
          CellOrTissueType = "HCT116", ExperimentalGroup = compartment,
          Replicate = "experiment-wide",
          EnrichmentStatus = if (enriched) "Enriched" else "Unenriched",
          AcquisitionMode = "DIA", Accession = accession,
          BaseAccession = base_accession(accession),
          GeneSymbol = if (length(genes)) genes[1] else "",
          ProteinName = if (length(names)) names[1] else "",
          KlaSite = paste0("K", position),
          ModifiedPeptide = clean_text(flank_col[i]),
          PEP = number(pg_pep_col[i]),
          Score = number(pg_cscore_col[i]),
          SourceFile = relative_path(path, project_root),
          SourceRow = i + 1L,
          SiteID = site_id,
          EvidenceMode = "DIA_ptm_site_report_experiment_wide",
          SourceConfidence = "standard_search_result"
        )
      }
    }
    if (no_lac) {
      logs[[length(logs) + 1]] <- exclusion_row(
        "PXD053474", relative_path(path, project_root),
        "non_lactyl_PTM_rows", no_lac
      )
    }
  }
  list(dataframe(rows), log_frame(logs))
}

extract_pxd053474_supplementary <- function(data_root, project_root, doi) {
  path <- file.path(data_root, "PXD053474", "supplementary", "js4c00366_si_003.xlsx")
  rows <- list()
  logs <- list()
  if (!file.exists(path)) {
    return(list(
      dataframe(list()),
      log_frame(list(exclusion_row("PXD053474", path, "missing_author_S3")))
    ))
  }
  for (spec in list(c("Subcellular", "Accession (553)"), c("Whole-cell lysates", "Accession"))) {
    sheet <- spec[1]
    accession_column <- spec[2]
    table <- suppressWarnings(read_excel(
      path, sheet = sheet, skip = 1, col_names = TRUE, col_types = "text"
    ))
    n <- nrow(table)
    acc_vec <- if (accession_column %in% names(table)) {
      table[[accession_column]]
    } else {
      rep("", n)
    }
    sites_vec <- if ("Sites" %in% names(table)) table$Sites else rep("", n)
    pep_vec <- if ("Peptide" %in% names(table)) table$Peptide else rep("", n)
    id_vec <- if ("ID" %in% names(table)) table$ID else rep("", n)
    for (i in seq_len(n)) {
      accession <- normalize_accession(acc_vec[i])
      sites_text <- clean_text(sites_vec[i])
      positions <- character(0)
      if (nzchar(sites_text)) {
        positions <- unlist(regmatches(
          sites_text,
          gregexpr("\\d+", sites_text, perl = TRUE)
        ))
      }
      if (!nzchar(accession) || !length(positions)) next
      for (position_raw in positions) {
        position <- suppressWarnings(as.integer(position_raw))
        site_id <- clean_text(id_vec[i])
        if (!nzchar(site_id)) {
          site_id <- paste0("S3:", sheet, ":", i + 2L, ":", accession,
                            ":K", position)
        }
        rows[[length(rows) + 1]] <- blank_site(
          PXD = "PXD053474", DOI = doi,
          SampleName = paste0("Author_S3_", sheet),
          CellOrTissueType = "HCT116", ExperimentalGroup = sheet,
          Replicate = "author table aggregate",
          EnrichmentStatus = "Enriched",
          AcquisitionMode = "DDA/DIA author union",
          Accession = accession,
          BaseAccession = base_accession(accession),
          GeneSymbol = "", ProteinName = "",
          KlaSite = paste0("K", position),
          ModifiedPeptide = clean_text(pep_vec[i]),
          ClassI = TRUE,
          SourceFile = relative_path(path, project_root),
          SourceRow = i + 2L,
          SiteID = site_id,
          EvidenceMode = "author_supplementary_table",
          SourceConfidence = "high_author_reported"
        )
      }
    }
  }
  list(dataframe(rows), log_frame(logs))
}

reconcile_pxd053474 <- function(dda, dia, supplementary) {
  evidence <- rbind_frames(list(dda, dia, supplementary))
  dda_keys <- unique(paste0(dda$BaseAccession, "\u0001", dda$KlaSite))
  dia_keys <- unique(paste0(dia$BaseAccession, "\u0001", dia$KlaSite))
  s3_keys <- unique(paste0(supplementary$BaseAccession, "\u0001", supplementary$KlaSite))
  all_keys <- sort(unique(c(dda_keys, dia_keys, s3_keys)), method = "radix")
  support_map <- list()
  comparison_rows <- list()
  for (key in all_keys) {
    parts <- strsplit(key, "\u0001", fixed = TRUE)[[1]]
    acc <- parts[1]
    site <- parts[2]
    in_dda <- key %in% dda_keys
    in_dia <- key %in% dia_keys
    in_s3 <- key %in% s3_keys
    in_search <- in_dda || in_dia
    comparison <- if (in_search && in_s3) {
      "search_and_supplementary"
    } else if (in_search) {
      "search_only"
    } else {
      "supplementary_only"
    }
    primary <- in_s3 || (in_dda && in_dia)
    confidence <- if (in_s3 && in_search) {
      "high_search_and_author_confirmed"
    } else if (in_s3) {
      "high_author_reported"
    } else if (in_dda && in_dia) {
      "high_cross_mode_search_support"
    } else {
      "moderate_single_mode_search_only"
    }
    support_map[[key]] <- list(
      PXD053_DDA = in_dda, PXD053_DIA = in_dia, PXD053_S3 = in_s3,
      PXD053Comparison = comparison, PrimaryIncluded = primary,
      SourceConfidence = confidence,
      InclusionStatus = if (primary) "included" else "audit_only",
      ExclusionReason = if (primary) {
        ""
      } else {
        "single_mode_search_only_not_in_author_S3"
      }
    )
    comparison_rows[[length(comparison_rows) + 1]] <- list(
      BaseAccession = acc, KlaSite = site, PXD053_DDA = in_dda,
      PXD053_DIA = in_dia, PXD053_S3 = in_s3,
      PXD053Comparison = comparison, PrimaryIncluded = primary,
      SourceConfidence = confidence,
      InclusionStatus = if (primary) "included" else "audit_only",
      ExclusionReason = if (primary) {
        ""
      } else {
        "single_mode_search_only_not_in_author_S3"
      }
    )
  }
  evidence$PXD053_DDA <- FALSE
  evidence$PXD053_DIA <- FALSE
  evidence$PXD053_S3 <- FALSE
  evidence$PXD053Comparison <- ""
  key_of <- paste0(evidence$BaseAccession, "\u0001", evidence$KlaSite)
  for (i in seq_len(nrow(evidence))) {
    sup <- support_map[[key_of[i]]]
    evidence$PXD053_DDA[i] <- sup$PXD053_DDA
    evidence$PXD053_DIA[i] <- sup$PXD053_DIA
    evidence$PXD053_S3[i] <- sup$PXD053_S3
    evidence$PXD053Comparison[i] <- sup$PXD053Comparison
    evidence$PrimaryIncluded[i] <- sup$PrimaryIncluded
    evidence$SourceConfidence[i] <- sup$SourceConfidence
    evidence$InclusionStatus[i] <- sup$InclusionStatus
    evidence$ExclusionReason[i] <- sup$ExclusionReason
  }
  comparison_frame <- frame_from_rows(
    comparison_rows,
    columns = c(
      "BaseAccession", "KlaSite", "PXD053_DDA", "PXD053_DIA", "PXD053_S3",
      "PXD053Comparison", "PrimaryIncluded", "SourceConfidence",
      "InclusionStatus", "ExclusionReason"
    ),
    logical_cols = c("PXD053_DDA", "PXD053_DIA", "PXD053_S3", "PrimaryIncluded")
  )
  list(evidence, comparison_frame)
}

extract_maxquant_site_table <- function(pxd, doi, path, project_root,
                                        cell_default, sample_map,
                                        evidence_mode) {
  table <- read_dtype_str(path, sep = "\t")
  n <- nrow(table)
  rows <- list()
  logs <- list()
  reason_counts <- list()
  labels <- sort(
    unique(sub("^Identification type ", "",
               names(table)[startsWith(names(table), "Identification type ")]),
           fixed = FALSE),
    method = "radix"
  )
  bad_cols <- rep(FALSE, n)
  for (cn in c("Reverse", "Potential contaminant", "Contaminant")) {
    if (cn %in% names(table)) bad_cols <- bad_cols | is_true(table[[cn]])
  }
  amino_col <- if ("Amino acid" %in% names(table)) table$`Amino acid` else rep("K", n)
  protein_col <- if ("Protein" %in% names(table)) {
    table$Protein
  } else if ("Leading proteins" %in% names(table)) {
    table$`Leading proteins`
  } else {
    rep("", n)
  }
  position_col <- if ("Position" %in% names(table)) {
    table$Position
  } else if ("Positions within proteins" %in% names(table)) {
    table$`Positions within proteins`
  } else {
    rep("", n)
  }
  fasta_col <- if ("Fasta headers" %in% names(table)) table$`Fasta headers` else rep("", n)
  gene_col <- if ("Gene names" %in% names(table)) table$`Gene names` else rep("", n)
  pname_col <- if ("Protein names" %in% names(table)) table$`Protein names` else rep("", n)
  id_col <- if ("id" %in% names(table)) table$id else rep("", n)
  diag_col <- if ("Diagnostic peak" %in% names(table)) table$`Diagnostic peak` else rep("", n)
  la_prob_col <- if ("La(K) Probabilities" %in% names(table)) {
    table[["La(K) Probabilities"]]
  } else if ("La (K) Probabilities" %in% names(table)) {
    table[["La (K) Probabilities"]]
  } else if ("Modified sequence" %in% names(table)) {
    table[["Modified sequence"]]
  } else {
    rep("", n)
  }
  loc_base_col <- if ("Localization prob" %in% names(table)) {
    table[["Localization prob"]]
  } else {
    rep("", n)
  }
  pep_base_col <- if ("PEP" %in% names(table)) table$PEP else rep("", n)
  score_base_col <- if ("Score" %in% names(table)) table$Score else rep("", n)
  ident_cols <- lapply(labels, function(lb) {
    cn <- paste0("Identification type ", lb)
    if (cn %in% names(table)) table[[cn]] else rep("", n)
  })
  loc_token_cols <- lapply(labels, function(lb) {
    cn <- paste0("Localization prob ", lb)
    if (cn %in% names(table)) table[[cn]] else NULL
  })
  pep_token_cols <- lapply(labels, function(lb) {
    cn <- paste0("PEP ", lb)
    if (cn %in% names(table)) table[[cn]] else NULL
  })
  score_token_cols <- lapply(labels, function(lb) {
    cn <- paste0("Score ", lb)
    if (cn %in% names(table)) table[[cn]] else NULL
  })
  for (i in seq_len(n)) {
    reason <- ""
    if (bad_cols[i]) reason <- "reverse_or_contaminant"
    if (!nzchar(reason) && toupper(clean_text(amino_col[i])) != "K") {
      reason <- "non_lysine_site"
    }
    accession <- normalize_accession(protein_col[i])
    position <- integer_value(position_col[i])
    if (!nzchar(reason) && (!nzchar(accession) || is.na(position))) {
      reason <- "missing_accession_or_kla_position"
    }
    if (nzchar(reason)) {
      reason_counts[[reason]] <- (reason_counts[[reason]] %||% 0L) + 1L
      next
    }
    ann <- annotation_from_fasta(fasta_col[i], accession)
    gene <- ann[1]
    protein_name <- ann[2]
    if (!nzchar(gene)) {
      gt <- split_tokens(gene_col[i])
      gene <- if (length(gt)) gt[1] else ""
    }
    if (!nzchar(protein_name)) {
      pt <- split_tokens(pname_col[i])
      protein_name <- if (length(pt)) pt[1] else ""
    }
    detected <- character(0)
    for (k in seq_along(labels)) {
      if (nzchar(clean_text(ident_cols[[k]][i]))) detected <- c(detected, labels[k])
    }
    if (!length(detected)) {
      if (length(labels)) {
        reason_counts[["site_not_identified_in_any_sample"]] <-
          (reason_counts[["site_not_identified_in_any_sample"]] %||% 0L) + 1L
        next
      }
      detected <- cell_default
    }
    for (token in detected) {
      mapped <- sample_map[[token]] %||% list()
      k <- match(token, labels)
      if (is.na(k)) k <- NULL
      loc_val <- if (!is.null(k) && !is.null(loc_token_cols[[k]])) {
        loc_token_cols[[k]][i]
      } else {
        loc_base_col[i]
      }
      pep_val <- if (!is.null(k) && !is.null(pep_token_cols[[k]])) {
        pep_token_cols[[k]][i]
      } else {
        pep_base_col[i]
      }
      score_val <- if (!is.null(k) && !is.null(score_token_cols[[k]])) {
        score_token_cols[[k]][i]
      } else {
        score_base_col[i]
      }
      localization <- number(loc_val)
      site_id <- clean_text(id_col[i])
      if (!nzchar(site_id)) {
        site_id <- paste0(pxd, ":", i + 1L, ":", accession, ":K", position)
      }
      rows[[length(rows) + 1]] <- blank_site(
        PXD = pxd, DOI = doi,
        SampleName = mapped$SampleName %||% token,
        CellOrTissueType = mapped$CellType %||% cell_default,
        ExperimentalGroup = mapped$ExperimentalGroup %||% token,
        Replicate = mapped$Replicate %||% token,
        EnrichmentStatus = mapped$EnrichmentStatus %||% "Enriched",
        AcquisitionMode = mapped$AcquisitionMode %||% "DDA",
        Accession = accession, BaseAccession = base_accession(accession),
        GeneSymbol = gene, ProteinName = protein_name,
        KlaSite = paste0("K", position),
        ModifiedPeptide = clean_text(la_prob_col[i]),
        LocalizationProb = localization,
        PEP = number(pep_val), Score = number(score_val),
        DiagnosticIon = clean_text(diag_col[i]),
        ClassI = !is.na(localization) && localization >= 0.75,
        SourceFile = relative_path(path, project_root),
        SourceRow = i + 1L, SiteID = site_id,
        EvidenceMode = evidence_mode,
        SourceConfidence = if (!is.na(localization) && localization >= 0.75) {
          "high_localized"
        } else {
          "standard_author_search_result"
        }
      )
    }
  }
  for (reason in names(reason_counts)) {
    logs[[length(logs) + 1]] <- exclusion_row(
      pxd, relative_path(path, project_root), reason, reason_counts[[reason]]
    )
  }
  list(dataframe(rows), log_frame(logs))
}

extract_pxd078013 <- function(data_root, project_root, doi, sample_map) {
  root <- file.path(data_root, "PXD078013", "search_results")
  evidence_path <- file.path(root, "evidence.txt")
  proteins_path <- file.path(root, "proteinGroups.txt")
  evidence <- read_dtype_str(evidence_path, sep = "\t")
  proteins <- read_dtype_str(proteins_path, sep = "\t")
  rows <- list()
  logs <- list()
  reason_counts <- list()
  n <- nrow(evidence)
  keep_rc <- rep(TRUE, n)
  for (cn in c("Reverse", "Potential contaminant")) {
    if (cn %in% names(evidence)) keep_rc <- keep_rc & !is_true(evidence[[cn]])
  }
  n_la <- if ("La (K)" %in% names(evidence)) {
    vapply(evidence[["La (K)"]], number, numeric(1))
  } else {
    rep(NA_real_, n)
  }
  keep_la <- !is.na(n_la) & n_la > 0
  ids_col <- if ("La (K) site IDs" %in% names(evidence)) {
    evidence[["La (K) site IDs"]]
  } else {
    rep("", n)
  }
  site_ids_rows <- lapply(seq_len(n), function(i) split_tokens(ids_col[i]))
  has_ids <- lengths(site_ids_rows) > 0
  exp_col <- if ("Experiment" %in% names(evidence)) {
    evidence$Experiment
  } else if ("Raw file" %in% names(evidence)) {
    evidence$`Raw file`
  } else {
    rep("", n)
  }
  modseq_col <- if ("Modified sequence" %in% names(evidence)) {
    evidence[["Modified sequence"]]
  } else {
    rep("", n)
  }
  pep_col <- if ("PEP" %in% names(evidence)) evidence$PEP else rep("", n)
  score_col <- if ("Score" %in% names(evidence)) evidence$Score else rep("", n)
  probs_col <- if ("La (K) Probabilities" %in% names(evidence)) {
    evidence[["La (K) Probabilities"]]
  } else {
    rep("", n)
  }
  for (i in seq_len(n)) {
    reason <- ""
    if (!keep_rc[i]) reason <- "reverse_or_contaminant"
    if (!nzchar(reason) && !keep_la[i]) reason <- "evidence_without_positive_LaK"
    if (!nzchar(reason) && !has_ids[i]) reason <- "positive_LaK_without_site_ID"
    if (nzchar(reason)) {
      reason_counts[[reason]] <- (reason_counts[[reason]] %||% 0L) + 1L
    }
  }
  sel <- which(keep_rc & keep_la & has_ids)
  expanded_site <- unlist(site_ids_rows[sel])
  expanded_idx <- rep(sel, lengths(site_ids_rows[sel]))
  evidence_by_site <- split(expanded_idx, expanded_site)
  locprob_vec <- vapply(seq_len(n), function(i) {
    p <- parse_probability_values(probs_col[i])
    if (length(p)) max(p) else NA_real_
  }, numeric(1))
  pep_vec <- vapply(pep_col, number, numeric(1))
  score_vec <- vapply(score_col, number, numeric(1))

  np <- nrow(proteins)
  keep_p <- rep(TRUE, np)
  for (cn in c("Reverse", "Potential contaminant", "Only identified by site")) {
    if (cn %in% names(proteins)) keep_p <- keep_p & !is_true(proteins[[cn]])
  }
  p_site_ids <- lapply(seq_len(np), function(i) {
    split_tokens(if ("La (K) site IDs" %in% names(proteins)) {
      proteins[["La (K) site IDs"]][i]
    } else {
      ""
    })
  })
  p_positions <- lapply(seq_len(np), function(i) {
    split_tokens(if ("La (K) site positions" %in% names(proteins)) {
      proteins[["La (K) site positions"]][i]
    } else {
      ""
    })
  })
  p_has <- lengths(p_site_ids) > 0 & lengths(p_positions) > 0
  majority_col <- if ("Majority protein IDs" %in% names(proteins)) {
    proteins[["Majority protein IDs"]]
  } else {
    rep("", np)
  }
  protein_ids_col <- if ("Protein IDs" %in% names(proteins)) {
    proteins[["Protein IDs"]]
  } else {
    rep("", np)
  }
  fasta_col <- if ("Fasta headers" %in% names(proteins)) {
    proteins[["Fasta headers"]]
  } else {
    rep("", np)
  }
  gene_col <- if ("Gene names" %in% names(proteins)) proteins[["Gene names"]] else rep("", np)
  pname_col <- if ("Protein names" %in% names(proteins)) {
    proteins[["Protein names"]]
  } else {
    rep("", np)
  }
  for (pi in seq_len(np)) {
    if (!keep_p[pi]) next
    site_ids <- p_site_ids[[pi]]
    positions_raw <- p_positions[[pi]]
    if (!length(site_ids) || !length(positions_raw)) next
    pairs <- site_position_pairs(site_ids, positions_raw)
    if (!length(pairs)) {
      reason_counts[["site_id_position_count_mismatch"]] <-
        (reason_counts[["site_id_position_count_mismatch"]] %||% 0L) + 1L
      next
    }
    accessions <- split_tokens(majority_col[pi])
    if (!length(accessions)) accessions <- split_tokens(protein_ids_col[pi])
    accession <- if (length(accessions)) normalize_accession(accessions[1]) else ""
    if (!nzchar(accession)) next
    ann <- annotation_from_fasta(fasta_col[pi], accession)
    gene <- ann[1]
    protein_name <- ann[2]
    if (!nzchar(gene)) {
      gt <- split_tokens(gene_col[pi])
      gene <- if (length(gt)) gt[1] else ""
    }
    if (!nzchar(protein_name)) {
      pt <- split_tokens(pname_col[pi])
      protein_name <- if (length(pt)) pt[1] else ""
    }
    for (pair in pairs) {
      position <- integer_value(pair[[2]])
      if (is.na(position)) next
      ev_idxs <- evidence_by_site[[pair[[1]]]]
      if (is.null(ev_idxs)) {
        reason_counts[["proteinGroups_site_without_positive_evidence"]] <-
          (reason_counts[["proteinGroups_site_without_positive_evidence"]] %||% 0L) + 1L
        next
      }
      for (ev_idx in ev_idxs) {
        token <- clean_text(exp_col[ev_idx])
        mapped <- sample_map[[token]] %||% list()
        rows[[length(rows) + 1]] <- blank_site(
          PXD = "PXD078013", DOI = doi,
          SampleName = mapped$SampleName %||% token,
          CellOrTissueType = "RKO",
          ExperimentalGroup = mapped$ExperimentalGroup %||% token,
          Replicate = mapped$Replicate %||% token,
          EnrichmentStatus = "Enriched", AcquisitionMode = "DDA",
          Accession = accession, BaseAccession = base_accession(accession),
          GeneSymbol = gene, ProteinName = protein_name,
          KlaSite = paste0("K", position),
          ModifiedPeptide = modseq_col[ev_idx],
          LocalizationProb = locprob_vec[ev_idx],
          PEP = pep_vec[ev_idx], Score = score_vec[ev_idx],
          ClassI = !is.na(locprob_vec[ev_idx]) && locprob_vec[ev_idx] >= 0.75,
          SourceFile = paste0(
            relative_path(proteins_path, project_root), ";",
            relative_path(evidence_path, project_root)
          ),
          SourceRow = paste0("proteinGroups:", pi + 1L, ";evidence:", ev_idx + 1L),
          SiteID = pair[[1]],
          EvidenceMode = "derived_proteinGroups_plus_Kla_positive_evidence",
          SourceConfidence = "moderate_two_table_derived_site"
        )
      }
    }
  }
  for (reason in names(reason_counts)) {
    logs[[length(logs) + 1]] <- exclusion_row(
      "PXD078013", relative_path(evidence_path, project_root),
      reason, reason_counts[[reason]]
    )
  }
  list(dataframe(rows), log_frame(logs))
}

# ---------------------------------------------------------------------------
# run_pipeline.py logic translation
# ---------------------------------------------------------------------------

write_csv_std <- function(frame, path) write_csv(frame, path)

as_bool <- function(series) {
  if (is.logical(series)) {
    series[is.na(series)] <- FALSE
    return(series)
  }
  tolower(ifelse(is.na(series), "", series)) %in% c("true", "1", "yes", "+")
}

sample_mapping <- function(config, pxd) {
  subset <- config[config$PXD == pxd, , drop = FALSE]
  out <- list()
  for (i in seq_len(nrow(subset))) {
    out[[subset$SourceToken[i]]] <- as.list(subset[i, , drop = FALSE])
  }
  out
}

ordered_sites <- function(values) {
  tokens <- character(0)
  for (value in values) {
    t <- strsplit(clean_text(value), ";", fixed = TRUE)[[1]]
    tokens <- c(tokens, t[nzchar(t)])
  }
  tokens <- unique(tokens)
  if (!length(tokens)) return("")
  key1 <- vapply(tokens, function(t) {
    if (startsWith(t, "K") && grepl("^[0-9]+$", substring(t, 2), perl = TRUE)) {
      as.numeric(substring(t, 2))
    } else {
      1e9
    }
  }, numeric(1))
  paste(tokens[order(key1, tokens, method = "radix")], collapse = ";")
}

max_number <- function(values) {
  numeric_values <- suppressWarnings(as.numeric(values))
  if (all(is.na(numeric_values))) return(NA_real_)
  max(numeric_values, na.rm = TRUE)
}

AGGREGATE_COLUMNS <- c(
  "BaseAccession", "GeneSymbol", "ProteinName", "KlaSites", "PXD", "Sample",
  "CellType", "ExperimentalGroup", "Category", "SourceFile", "EvidenceMode",
  "LocalizationProb", "SourceConfidence", "EvidenceRows"
)

aggregate_evidence <- function(frame) {
  if (!nrow(frame)) {
    out <- as.data.frame(
      setNames(lapply(AGGREGATE_COLUMNS, function(c) character(0)),
               AGGREGATE_COLUMNS),
      stringsAsFactors = FALSE
    )
    out$LocalizationProb <- numeric(0)
    out$EvidenceRows <- integer(0)
    return(out)
  }
  by_acc <- split(seq_len(nrow(frame)), frame$BaseAccession)
  accs <- sort(names(by_acc), method = "radix")
  has_category <- "Category" %in% names(frame)
  rows <- lapply(accs, function(acc) {
    idx <- by_acc[[acc]]
    g <- frame[idx, , drop = FALSE]
    list(
      BaseAccession = acc,
      GeneSymbol = best_annotation(g$GeneSymbol),
      ProteinName = best_annotation(g$ProteinName),
      KlaSites = ordered_sites(g$KlaSite),
      PXD = unique_join(g$PXD),
      Sample = unique_join(g$SampleName),
      CellType = unique_join(g$CellOrTissueType),
      ExperimentalGroup = unique_join(g$ExperimentalGroup),
      Category = if (has_category) unique_join(g$Category) else "",
      SourceFile = unique_join(g$SourceFile),
      EvidenceMode = unique_join(g$EvidenceMode),
      LocalizationProb = max_number(g$LocalizationProb),
      SourceConfidence = unique_join(g$SourceConfidence),
      EvidenceRows = length(idx)
    )
  })
  frame_from_rows(
    rows, columns = AGGREGATE_COLUMNS,
    numeric_cols = "LocalizationProb", integer_cols = "EvidenceRows"
  )
}

cell_type_statistics <- function(evidence, go_accessions, grouping) {
  rows <- list()
  category_order <- setNames(seq_along(CATEGORIES) - 1L, CATEGORIES)
  ordered_grouping <- grouping[
    order(category_order[grouping$teacher_requested_grouping], method = "radix"),
    , drop = FALSE
  ]
  for (i in seq_len(nrow(ordered_grouping))) {
    cell_type <- ordered_grouping$CellType[i]
    subset <- evidence[evidence$CellOrTissueType == cell_type, , drop = FALSE]
    if (!nrow(subset)) next
    accessions <- unique(subset$BaseAccession)
    ddr_accessions <- intersect(accessions, go_accessions)
    total <- length(accessions)
    rows[[length(rows) + 1]] <- list(
      CellOrTissueType = cell_type,
      TotalKlaProteins = total,
      KlaGoDdrProteins = length(ddr_accessions),
      KlaGoDdrFraction = if (total) length(ddr_accessions) / total else 0
    )
  }
  frame_from_rows(
    rows,
    columns = c("CellOrTissueType", "TotalKlaProteins", "KlaGoDdrProteins",
                "KlaGoDdrFraction"),
    integer_cols = c("TotalKlaProteins", "KlaGoDdrProteins"),
    numeric_cols = "KlaGoDdrFraction"
  )
}

extract_all <- function(project_root) {
  data_root <- file.path(project_root, "data")
  config_root <- file.path(project_root, "reanalysis", "config")
  dataset_config <- read_dtype_str(file.path(config_root, "datasets.csv"))
  sample_config <- read_dtype_str(file.path(config_root, "sample_map.csv"))
  doi <- setNames(dataset_config$DOI, dataset_config$PXD)
  frames <- list()
  logs <- list()
  audits <- list()

  r <- extract_pxd014870(data_root, project_root, doi[["PXD014870"]])
  frames[[length(frames) + 1]] <- r[[1]]
  logs[[length(logs) + 1]] <- r[[2]]

  r <- extract_pxd028488(data_root, project_root, doi[["PXD028488"]])
  frames[[length(frames) + 1]] <- r[[1]]
  logs[[length(logs) + 1]] <- r[[2]]
  audits$pxd028_directory_audit <- r[[3]]

  r <- extract_pxd050470(data_root, project_root, doi[["PXD050470"]])
  frames[[length(frames) + 1]] <- r[[1]]
  logs[[length(logs) + 1]] <- r[[2]]

  dda <- extract_pxd053474_dda(data_root, project_root, doi[["PXD053474"]])
  dia <- extract_pxd053474_dia(data_root, project_root, doi[["PXD053474"]])
  supplementary <- extract_pxd053474_supplementary(
    data_root, project_root, doi[["PXD053474"]]
  )
  reconciled <- reconcile_pxd053474(dda[[1]], dia[[1]], supplementary[[1]])
  frames[[length(frames) + 1]] <- reconciled[[1]]
  logs[[length(logs) + 1]] <- dda[[2]]
  logs[[length(logs) + 1]] <- dia[[2]]
  logs[[length(logs) + 1]] <- supplementary[[2]]
  audits$pxd053_comparison <- reconciled[[2]]
  audits$pxd053_dda <- dda[[1]]
  audits$pxd053_dia <- dia[[1]]
  audits$pxd053_supplementary <- supplementary[[1]]

  r <- extract_maxquant_site_table(
    "PXD060185", doi[["PXD060185"]],
    file.path(data_root, "PXD060185", "search_results", "RESULT", "combined",
              "txt", "La (K)Sites.txt"),
    project_root, "breast cell line",
    sample_mapping(sample_config, "PXD060185"),
    "maxquant_LaK_site_table"
  )
  frames[[length(frames) + 1]] <- r[[1]]
  logs[[length(logs) + 1]] <- r[[2]]

  r <- extract_pxd078013(
    data_root, project_root, doi[["PXD078013"]],
    sample_mapping(sample_config, "PXD078013")
  )
  frames[[length(frames) + 1]] <- r[[1]]
  logs[[length(logs) + 1]] <- r[[2]]

  r <- extract_maxquant_site_table(
    "PXD078736", doi[["PXD078736"]],
    file.path(data_root, "PXD078736", "search_results", "txt", "La(K)Sites.txt"),
    project_root, "HK-2",
    sample_mapping(sample_config, "PXD078736"),
    "maxquant_LaK_site_table"
  )
  frames[[length(frames) + 1]] <- r[[1]]
  logs[[length(logs) + 1]] <- r[[2]]

  all_evidence <- rbind_frames(frames)
  all_evidence <- apply_annotation_supplement(
    all_evidence,
    file.path(config_root, "uniprot_annotation_supplement.tsv")
  )
  all_evidence$PrimaryIncluded <- as_bool(all_evidence$PrimaryIncluded)
  all_evidence$ClassI <- as_bool(all_evidence$ClassI)
  ord <- order(
    all_evidence$PXD, all_evidence$CellOrTissueType, all_evidence$SampleName,
    all_evidence$BaseAccession, all_evidence$KlaSite, all_evidence$SourceFile,
    all_evidence$SourceRow, method = "radix", na.last = TRUE
  )
  all_evidence <- all_evidence[ord, , drop = FALSE]
  rownames(all_evidence) <- NULL
  nonempty_logs <- logs[vapply(logs, nrow, integer(1)) > 0]
  exclusion_log <- if (length(nonempty_logs)) {
    do.call(rbind, nonempty_logs)
  } else {
    log_frame(list())
  }
  list(all_evidence, exclusion_log, audits)
}

attach_go <- function(proteins, go_summary, go_raw) {
  go_columns <- c(
    "GOSymbol", "GOTerms", "GONames", "GOEvidenceCodes", "GOReferences",
    "GOAnnotationCount"
  )
  merged <- proteins
  idx <- match(merged$BaseAccession, go_summary$BaseAccession)
  hit <- !is.na(idx)
  for (col in go_columns) {
    merged[[col]] <- if (col == "GOAnnotationCount") {
      rep(NA_real_, nrow(merged))
    } else {
      rep(NA_character_, nrow(merged))
    }
  }
  for (col in go_columns) {
    if (any(hit)) merged[[col]][hit] <- go_summary[[col]][idx[hit]]
  }
  merged$GOMatchMode <- ifelse(
    !is.na(merged$GONames) & nzchar(merged$GONames),
    "BaseAccession", "unmatched"
  )

  taxon <- suppressWarnings(as.numeric(go_raw[["TAXON ID"]]))
  taxon[is.na(taxon)] <- 9606
  retained <- go_raw[
    !go_raw$ExcludedNOT & taxon == 9606 & go_raw$BaseAccession != "",
    , drop = FALSE
  ]
  retained$SymbolKey <- toupper(trimws(retained$SYMBOL))
  gene_rows <- list()
  symbols_all <- sort(
    unique(retained$SymbolKey[retained$SymbolKey != ""]),
    method = "radix"
  )
  for (symbol in symbols_all) {
    g <- retained[retained$SymbolKey == symbol, , drop = FALSE]
    gene_rows[[symbol]] <- list(
      SymbolKey = symbol,
      GOSymbol = best_annotation(g$SYMBOL),
      GOTerms = unique_join(g[["GO TERM"]]),
      GONames = unique_join(g[["GO NAME"]]),
      GOEvidenceCodes = unique_join(g[["GO EVIDENCE CODE"]]),
      GOReferences = unique_join(g[["REFERENCE"]]),
      GOAnnotationCount = nrow(g)
    )
  }
  for (i in which(merged$GOMatchMode == "unmatched")) {
    symbols <- toupper(split_tokens(merged$GeneSymbol[i]))
    match_row <- NULL
    for (s in symbols) {
      if (!is.null(gene_rows[[s]])) {
        match_row <- gene_rows[[s]]
        break
      }
    }
    if (is.null(match_row)) next
    for (col in go_columns) merged[[col]][i] <- match_row[[col]]
    merged$GOMatchMode[i] <- "GeneSymbol_fallback"
  }
  merged
}

venn_regions <- function(sets) {
  h <- sets[["hippocampus_tissue"]]
  n <- sets[["normal_immortalized_cell_lines"]]
  t <- sets[["tumor_cell_lines"]]
  list(
    hippocampus_only = setdiff(h, c(n, t)),
    normal_only = setdiff(n, c(h, t)),
    tumor_only = setdiff(t, c(h, n)),
    hippocampus_and_normal_only = setdiff(intersect(h, n), t),
    hippocampus_and_tumor_only = setdiff(intersect(h, t), n),
    normal_and_tumor_only = setdiff(intersect(n, t), h),
    all_three = intersect(intersect(h, n), t)
  )
}

plot_venn <- function(regions, output_stem, title) {
  region_order <- names(regions)
  sets_map <- list(
    hippocampus_only = "Hippocampus tissue",
    normal_only = "Immortalized models",
    tumor_only = "Tumor cell lines",
    hippocampus_and_normal_only = c("Hippocampus tissue", "Immortalized models"),
    hippocampus_and_tumor_only = c("Hippocampus tissue", "Tumor cell lines"),
    normal_and_tumor_only = c("Immortalized models", "Tumor cell lines"),
    all_three = c("Hippocampus tissue", "Immortalized models", "Tumor cell lines")
  )
  value <- character(0)
  group <- character(0)
  for (rn in region_order) {
    for (v in regions[[rn]]) {
      value <- c(value, v)
      group <- c(group, sets_map[[rn]])
    }
  }
  items <- data.frame(value = value, Group = group, stringsAsFactors = FALSE)
  p <- ggVennDiagram(
    items,
    label = "count",
    label_alpha = 0,
    label_size = 10,
    set_size = 4,
    category.names = c(
      "Hippocampus tissue", "Immortalized models", "Tumor cell lines"
    )
  ) +
    scale_fill_gradient(low = "#F7F3E8", high = "#2A6F97") +
    ggplot2::labs(title = title) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 13, hjust = 0.5),
      legend.position = "none"
    )
  dir.create(dirname(output_stem), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(
    paste0(output_stem, ".png"), p, width = 7.2, height = 6.2, dpi = 300,
    bg = "white"
  )
  grDevices::cairo_pdf(
    file = paste0(output_stem, ".pdf"), width = 7.2, height = 6.2
  )
  print(p)
  grDevices::dev.off()
  grDevices::svg(
    file = paste0(output_stem, ".svg"), width = 7.2, height = 6.2
  )
  print(p)
  grDevices::dev.off()
  invisible(NULL)
}

write_group_outputs <- function(evidence, grouping, scheme, analysis_name,
                                table_root, figure_root) {
  mapping <- setNames(grouping[[scheme]], grouping$CellType)
  assigned <- evidence
  assigned$Category <- unname(mapping[assigned$CellOrTissueType])
  assigned$Category[is.na(assigned$Category)] <- ""
  unmapped <- unique(assigned[assigned$Category == "",
                              c("PXD", "CellOrTissueType", "SampleName",
                                "BaseAccession"), drop = FALSE])
  if (nrow(unmapped)) {
    stop(sprintf(
      "Unmapped cell types for %s: %s",
      scheme,
      paste(sort(unique(unmapped$CellOrTissueType)), collapse = ", ")
    ))
  }
  sets <- lapply(CATEGORIES, function(category) {
    unique(assigned$BaseAccession[assigned$Category == category])
  })
  names(sets) <- CATEGORIES
  regions <- venn_regions(sets)
  output_dir <- file.path(table_root, "venn_regions", scheme, analysis_name)
  counts <- list()
  membership_rows <- list()
  for (region_name in REGION_ORDER) {
    accessions <- regions[[region_name]]
    region_frame <- aggregate_evidence(
      assigned[assigned$BaseAccession %in% accessions, , drop = FALSE]
    )
    write_csv(region_frame, file.path(output_dir, paste0(region_name, ".csv")))
    counts[[length(counts) + 1]] <- list(
      Region = region_name, ProteinCount = length(accessions)
    )
    for (accession in sort(accessions, method = "radix")) {
      membership_rows[[length(membership_rows) + 1]] <- list(
        BaseAccession = accession, Region = region_name
      )
    }
  }
  count_frame <- frame_from_rows(
    counts, columns = c("Region", "ProteinCount"), integer_cols = "ProteinCount"
  )
  write_csv(count_frame, file.path(output_dir, "venn_region_counts.csv"))
  membership_frame <- frame_from_rows(
    membership_rows, columns = c("BaseAccession", "Region")
  )
  write_csv(membership_frame, file.path(output_dir, "venn_membership.csv"))
  combined <- aggregate_evidence(assigned)
  region_by_accession <- list()
  for (region_name in REGION_ORDER) {
    for (acc in regions[[region_name]]) region_by_accession[[acc]] <- region_name
  }
  membership_columns <- c(
    InHippocampusTissue = "hippocampus_tissue",
    InNormalImmortalizedCellLines = "normal_immortalized_cell_lines",
    InTumorCellLines = "tumor_cell_lines"
  )
  for (col in names(membership_columns)) {
    category <- membership_columns[[col]]
    combined[[col]] <- ifelse(
      combined$BaseAccession %in% sets[[category]], "Yes", "No"
    )
  }
  combined$DetectedGroupCount <- as.integer(rowSums(
    vapply(names(membership_columns), function(col) {
      combined[[col]] == "Yes"
    }, logical(nrow(combined)))
  ))
  combined$VennRegion <- vapply(
    combined$BaseAccession,
    function(acc) region_by_accession[[acc]] %||% "", character(1)
  )
  combined$AnalysisSet <- analysis_name
  leading_columns <- c(
    "BaseAccession", "GeneSymbol", "ProteinName", "KlaSites",
    "InHippocampusTissue", "InNormalImmortalizedCellLines", "InTumorCellLines",
    "DetectedGroupCount", "VennRegion", "AnalysisSet"
  )
  combined <- combined[
    , c(leading_columns,
        names(combined)[!(names(combined) %in% leading_columns)]),
    drop = FALSE
  ]
  write_csv(
    combined, file.path(output_dir, paste0(analysis_name, "_three_groups_combined.csv"))
  )
  category_frames <- list()
  for (category in CATEGORIES) {
    category_frame <- aggregate_evidence(
      assigned[assigned$BaseAccession %in% sets[[category]], , drop = FALSE]
    )
    write_csv(category_frame, file.path(output_dir, paste0(category, "_all.csv")))
    category_frame <- cbind(
      category_frame[, 1:4, drop = FALSE],
      data.frame(SourceCategory = category, stringsAsFactors = FALSE),
      category_frame[, 5:ncol(category_frame), drop = FALSE]
    )
    category_frames[[length(category_frames) + 1]] <- category_frame
  }
  non_deduplicated <- do.call(rbind, category_frames)
  for (col in names(membership_columns)) {
    category <- membership_columns[[col]]
    non_deduplicated[[col]] <- ifelse(
      non_deduplicated$BaseAccession %in% sets[[category]], "Yes", "No"
    )
  }
  non_deduplicated$DetectedGroupCount <- as.integer(rowSums(
    vapply(names(membership_columns), function(col) {
      non_deduplicated[[col]] == "Yes"
    }, logical(nrow(non_deduplicated)))
  ))
  non_deduplicated$VennRegion <- vapply(
    non_deduplicated$BaseAccession,
    function(acc) region_by_accession[[acc]] %||% "", character(1)
  )
  non_deduplicated$AnalysisSet <- analysis_name
  leading_columns_non_dedup <- c(
    leading_columns[1:4], "SourceCategory", leading_columns[5:length(leading_columns)]
  )
  non_deduplicated <- non_deduplicated[
    , c(leading_columns_non_dedup,
        names(non_deduplicated)[!(names(non_deduplicated) %in% leading_columns_non_dedup)]),
    drop = FALSE
  ]
  write_csv(
    non_deduplicated,
    file.path(
      output_dir,
      paste0(analysis_name, "_three_groups_combined_non_deduplicated.csv")
    )
  )
  plot_venn(
    regions,
    file.path(figure_root, scheme, paste0(analysis_name, "_three_group_venn")),
    if (analysis_name == "all_kla") {
      "All Kla proteins"
    } else {
      "Kla and GO repair/damage proteins"
    }
  )
  list(count_frame, regions, combined, non_deduplicated)
}

target_trace <- function(all_evidence) {
  targets <- list(
    MRE11 = list(genes = c("MRE11"), accessions = c("P49959")),
    `XLF/NHEJ1` = list(genes = c("NHEJ1", "XLF"), accessions = c("Q9H9Q4")),
    `NBS1/NBN` = list(genes = c("NBN", "NBS1"), accessions = c("O60934"))
  )
  rows <- list()
  for (pxd in INCLUDED_PXDS) {
    dataset <- all_evidence[all_evidence$PXD == pxd, , drop = FALSE]
    for (target in names(targets)) {
      aliases <- targets[[target]]
      gene_match <- vapply(
        dataset$GeneSymbol,
        function(value) {
          tokens <- toupper(trimws(strsplit(clean_text(value), ";", fixed = TRUE)[[1]]))
          any(tokens %in% aliases$genes)
        },
        logical(1)
      )
      accession_match <- dataset$BaseAccession %in% aliases$accessions
      matched <- dataset[gene_match | accession_match, , drop = FALSE]
      primary <- matched[matched$PrimaryIncluded, , drop = FALSE]
      status <- if (nrow(primary)) {
        "present_in_primary_kla_evidence"
      } else if (nrow(matched)) {
        "present_only_in_audit_kla_evidence"
      } else {
        "not_present_in_extracted_kla_evidence"
      }
      rows[[length(rows) + 1]] <- list(
        PXD = pxd,
        Target = target,
        TargetGenes = paste(sort(aliases$genes), collapse = ";"),
        TargetAccessions = paste(sort(aliases$accessions), collapse = ";"),
        Status = status,
        PrimaryEvidenceRows = nrow(primary),
        AllEvidenceRows = nrow(matched),
        KlaSites = if (nrow(matched)) ordered_sites(matched$KlaSite) else "",
        Samples = if (nrow(matched)) unique_join(matched$SampleName) else "",
        SourceFiles = if (nrow(matched)) unique_join(matched$SourceFile) else ""
      )
    }
  }
  frame_from_rows(
    rows,
    columns = c(
      "PXD", "Target", "TargetGenes", "TargetAccessions", "Status",
      "PrimaryEvidenceRows", "AllEvidenceRows", "KlaSites", "Samples",
      "SourceFiles"
    ),
    integer_cols = c("PrimaryEvidenceRows", "AllEvidenceRows")
  )
}

regression_outputs <- function(new_proteins, project_root) {
  old_path <- file.path(
    project_root,
    "archive/reanalysis_v1_2026-07-21/outputs/01_kla_extraction",
    "kla_proteins_by_dataset_cell_type_group.csv"
  )
  empty_detail <- frame_from_rows(
    list(),
    columns = c("PXD", "BaseAccession", "RegressionStatus", "Reason")
  )
  empty_summary <- frame_from_rows(
    list(),
    columns = c("PXD", "RegressionStatus", "ProteinCount"),
    integer_cols = "ProteinCount"
  )
  if (!file.exists(old_path)) return(list(empty_detail, empty_summary))
  old <- read_dtype_str(old_path)
  old_keys <- paste0(old$PXD, "\u0001", old$BaseAccession)
  new_keys <- paste0(new_proteins$PXD, "\u0001", new_proteins$BaseAccession)
  all_keys <- sort(unique(c(old_keys, new_keys)), method = "radix")
  reason_map <- c(
    PXD028488 = "expanded_all_human_PEAKS_directories",
    PXD050470 = "corrected_supplementary_header_and_S3_S12_union",
    PXD053474 = "four_search_suites_plus_reconciled_S3",
    PXD078013 = "strict_evidence_plus_proteinGroups_site_link"
  )
  rows <- lapply(all_keys, function(key) {
    parts <- strsplit(key, "\u0001", fixed = TRUE)[[1]]
    pxd <- parts[1]
    accession <- parts[2]
    in_old <- key %in% old_keys
    in_new <- key %in% new_keys
    if (in_old && in_new) {
      status <- "shared"
      reason <- "present_in_both"
    } else if (in_new) {
      status <- "new_only"
      reason <- unname(reason_map[pxd])
      if (is.na(reason)) reason <- "revised_sample_level_extraction"
    } else {
      status <- "old_only"
      reason <- "removed_by_revised_traceable_evidence_rule_or_accession_normalization"
    }
    list(PXD = pxd, BaseAccession = accession, RegressionStatus = status,
         Reason = reason)
  })
  detail <- frame_from_rows(
    rows,
    columns = c("PXD", "BaseAccession", "RegressionStatus", "Reason")
  )
  summary_rows <- list()
  for (pxd in sort(unique(detail$PXD), method = "radix")) {
    for (status in sort(unique(detail$RegressionStatus), method = "radix")) {
      cnt <- sum(detail$PXD == pxd & detail$RegressionStatus == status)
      if (cnt) {
        summary_rows[[length(summary_rows) + 1]] <- list(
          PXD = pxd, RegressionStatus = status, ProteinCount = cnt
        )
      }
    }
  }
  summary <- frame_from_rows(
    summary_rows,
    columns = c("PXD", "RegressionStatus", "ProteinCount"),
    integer_cols = "ProteinCount"
  )
  list(detail, summary)
}

# ---------------------------------------------------------------------------
# audit_target_sources.py translation (target source audit)
# ---------------------------------------------------------------------------

AUDIT_TARGETS <- list(
  MRE11 = list(accession = "P49959", genes = c("MRE11")),
  `XLF/NHEJ1` = list(accession = "Q9H9Q4", genes = c("NHEJ1", "XLF")),
  `NBS1/NBN` = list(accession = "O60934", genes = c("NBN", "NBS1"))
)

target_pattern <- function(accession, genes) {
  aliases <- unique(c(accession, genes))
  aliases <- aliases[order(-nchar(aliases), aliases)]
  paste0(
    "(?<![A-Za-z0-9])(?:", paste(aliases, collapse = "|"),
    ")(?![A-Za-z0-9])"
  )
}

target_mask <- function(table, pattern) {
  likely <- grep(
    "accession|protein|gene|description|fasta", names(table),
    ignore.case = TRUE, value = TRUE
  )
  columns <- if (length(likely)) likely else names(table)
  mask <- rep(FALSE, nrow(table))
  for (cn in columns) {
    m <- grepl(pattern, ifelse(is.na(table[[cn]]), "", table[[cn]]),
               ignore.case = TRUE, perl = TRUE)
    mask <- mask | m
  }
  mask
}

lactyl_mask <- function(table) {
  mask <- rep(FALSE, nrow(table))
  for (cn in names(table)) {
    if (grepl("peptide|ptm|ascore|modification|modified", cn,
              ignore.case = TRUE, perl = TRUE)) {
      m <- grepl(
        "Lactyl|\\bLac\\b|K\\(\\+?72\\.02|K\\(Lac|La\\s*\\(K\\)",
        ifelse(is.na(table[[cn]]), "", table[[cn]]),
        ignore.case = TRUE, perl = TRUE
      )
      mask <- mask | m
    }
    if (grepl("La\\s*\\(K\\)", cn, ignore.case = TRUE, perl = TRUE)) {
      numeric_values <- suppressWarnings(as.numeric(table[[cn]]))
      numeric_values[is.na(numeric_values)] <- 0
      mask <- mask | numeric_values > 0
      if (grepl("site\\s*IDs?", cn, ignore.case = TRUE, perl = TRUE)) {
        mask <- mask | nzchar(trimws(ifelse(is.na(table[[cn]]), "", table[[cn]])))
      }
    }
  }
  if ("PTM.ModificationTitle" %in% names(table)) {
    mask <- mask | tolower(ifelse(is.na(table[["PTM.ModificationTitle"]]), "",
                                  table[["PTM.ModificationTitle"]])) == "lac"
  }
  mask
}

read_table_audit <- function(path) {
  sep <- if (tolower(tools::file_ext(path)) %in% c("txt", "tsv")) "\t" else ","
  read_dtype_str(path, sep = sep)
}

scan_text_file <- function(path, project_root, rows, kla_defining) {
  table <- read_table_audit(path)
  lactyl <- if (kla_defining) lactyl_mask(table) else rep(FALSE, nrow(table))
  for (target in names(AUDIT_TARGETS)) {
    t <- AUDIT_TARGETS[[target]]
    matches <- target_mask(table, target_pattern(t$accession, t$genes))
    if (!any(matches)) next
    rows[[length(rows) + 1]] <- list(
      Target = target,
      SourceFile = relative_path(path, project_root),
      TargetRows = sum(matches),
      KlaTargetRows = if (kla_defining) sum(matches & lactyl) else 0L,
      SourceRole = if (kla_defining) {
        "kla_candidate_table"
      } else {
        "protein_or_peptide_support_table"
      }
    )
  }
  invisible(rows)
}

scan_excel_file <- function(path, project_root, rows) {
  sheets <- excel_sheets(path)
  for (sheet in sheets) {
    table <- suppressWarnings(read_excel(
      path, sheet = sheet, col_names = FALSE, col_types = "text"
    ))
    if (!nrow(table)) next
    for (target in names(AUDIT_TARGETS)) {
      t <- AUDIT_TARGETS[[target]]
      pattern <- target_pattern(t$accession, t$genes)
      matches <- rep(FALSE, nrow(table))
      for (cn in names(table)) {
        m <- grepl(ifelse(is.na(table[[cn]]), "", table[[cn]]),
                   pattern, ignore.case = TRUE, perl = TRUE)
        matches <- matches | m
      }
      if (!any(matches)) next
      rows[[length(rows) + 1]] <- list(
        Target = target,
        SourceFile = paste0(relative_path(path, project_root), "#", sheet),
        TargetRows = sum(matches),
        KlaTargetRows = sum(matches),
        SourceRole = "author_kla_supplementary_table"
      )
    }
  }
  invisible(rows)
}

source_files_audit <- function(project_root, pxd) {
  root <- file.path(project_root, "data", pxd)
  files <- list()
  append_paths <- function(pattern, kla) {
    paths <- list.files(root, pattern = pattern, recursive = TRUE,
                        full.names = TRUE, sort = FALSE)
    for (p in paths) {
      files[[length(files) + 1]] <<- list(path = p, kla = kla)
    }
  }
  if (pxd == "PXD014870") {
    for (p in list.files(root, pattern = "^Lactyl \\(K\\)Sites\\.txt$",
                         recursive = TRUE, full.names = TRUE, sort = FALSE)) {
      files[[length(files) + 1]] <- list(
        path = p, kla = grepl("Sites", basename(p))
      )
    }
    append_paths("^modificationSpecificPeptides\\.txt$", TRUE)
    append_paths("^proteinGroups\\.txt$", FALSE)
  } else if (pxd == "PXD028488") {
    append_paths("^DB search psm\\.csv$", TRUE)
    append_paths("^protein-peptides\\.csv$", FALSE)
    append_paths("^proteins\\.csv$", FALSE)
  } else if (pxd == "PXD053474") {
    append_paths("^peptide\\.csv$", TRUE)
    append_paths("ptm-site.*Report\\.tsv$", TRUE)
    append_paths("^protein-peptides\\.csv$", FALSE)
    append_paths("^proteins\\.csv$", FALSE)
  }
  keys <- unique(vapply(files, function(f) f$path, character(1)))
  keys <- sort(keys, method = "radix")
  lapply(keys, function(k) {
    list(path = k, kla = files[[match(k, vapply(files, function(f) f$path, character(1)))]][["kla"]])
  })
}

build_target_source_audit <- function(project_root, extracted_evidence) {
  scan_rows <- list()
  scanned_pxds <- c("PXD014870", "PXD028488", "PXD050470", "PXD053474")
  excel_files <- list(
    PXD050470 = sort(
      list.files(file.path(project_root, "data", "PXD050470", "supplementary"),
                 pattern = "\\.xlsx$", full.names = TRUE, sort = FALSE),
      method = "radix"
    ),
    PXD053474 = list(
      file.path(project_root, "data", "PXD053474", "supplementary",
                "js4c00366_si_003.xlsx")
    )
  )
  for (pxd in scanned_pxds) {
    local_rows <- list()
    for (f in source_files_audit(project_root, pxd)) {
      local_rows <- scan_text_file(f$path, project_root, local_rows, f$kla)
    }
    for (p in excel_files[[pxd]]) {
      if (file.exists(p)) local_rows <- scan_excel_file(p, project_root, local_rows)
    }
    if (length(local_rows)) {
      for (k in seq_along(local_rows)) local_rows[[k]]$PXD <- pxd
    }
    scan_rows <- c(scan_rows, local_rows)
  }
  scan <- if (length(scan_rows)) {
    frame_from_rows(
      scan_rows,
      columns = c("Target", "SourceFile", "TargetRows", "KlaTargetRows",
                  "SourceRole", "PXD"),
      integer_cols = c("TargetRows", "KlaTargetRows")
    )
  } else {
    frame_from_rows(
      list(),
      columns = c("Target", "SourceFile", "TargetRows", "KlaTargetRows",
                  "SourceRole", "PXD"),
      integer_cols = c("TargetRows", "KlaTargetRows")
    )
  }
  output_rows <- list()
  for (pxd in c("PXD014870", "PXD028488", "PXD050470", "PXD053474",
                "PXD060185", "PXD078013", "PXD078736")) {
    for (target in names(AUDIT_TARGETS)) {
      t <- AUDIT_TARGETS[[target]]
      extracted <- extracted_evidence[
        extracted_evidence$PXD == pxd &
          (extracted_evidence$BaseAccession == t$accession |
             vapply(extracted_evidence$GeneSymbol, function(value) {
               tokens <- toupper(trimws(strsplit(clean_text(value), ";", fixed = TRUE)[[1]]))
               any(tokens %in% t$genes)
             }, logical(1))),
        , drop = FALSE
      ]
      source <- if (nrow(scan)) {
        scan[scan$PXD == pxd & scan$Target == target, , drop = FALSE]
      } else {
        scan
      }
      target_rows <- if (nrow(source)) sum(source$TargetRows) else 0L
      kla_rows <- if (nrow(source)) sum(source$KlaTargetRows) else 0L
      conclusion <- if (nrow(extracted)) {
        "present_in_extracted_primary_kla_evidence"
      } else if (kla_rows) {
        "kla_candidate_rows_present_but_not_selected_review_required"
      } else if (target_rows) {
        "protein_or_unmodified_peptide_present_but_no_kla_evidence"
      } else if (pxd %in% scanned_pxds) {
        "not_found_in_available_author_search_or_kla_tables"
      } else {
        "not_present_in_extracted_kla_evidence_source_not_deep_scanned"
      }
      output_rows[[length(output_rows) + 1]] <- list(
        PXD = pxd,
        Target = target,
        Accession = t$accession,
        GeneAliases = paste(sort(t$genes), collapse = ";"),
        SourceTargetRows = target_rows,
        SourceKlaTargetRows = kla_rows,
        ExtractedPrimaryRows = nrow(extracted),
        Conclusion = conclusion,
        SourceFiles = if (nrow(source)) unique_join(source$SourceFile) else "",
        Interpretation = if (conclusion ==
                             "protein_or_unmodified_peptide_present_but_no_kla_evidence") {
          "The target was detected at protein/unmodified-peptide level, but the available author tables do not assign Kla to it."
        } else {
          ""
        }
      )
    }
  }
  frame_from_rows(
    output_rows,
    columns = c(
      "PXD", "Target", "Accession", "GeneAliases", "SourceTargetRows",
      "SourceKlaTargetRows", "ExtractedPrimaryRows", "Conclusion",
      "SourceFiles", "Interpretation"
    ),
    integer_cols = c("SourceTargetRows", "SourceKlaTargetRows",
                     "ExtractedPrimaryRows")
  )
}

# ---------------------------------------------------------------------------
# Pipeline stage
# ---------------------------------------------------------------------------

run_pipeline_stage <- function(project_root) {
  reanalysis <- file.path(project_root, "reanalysis")
  intermediate <- file.path(reanalysis, "intermediate")
  tables <- file.path(reanalysis, "results", "tables")
  figures <- file.path(reanalysis, "results", "figures")
  reports <- file.path(reanalysis, "reports")
  logs_dir <- file.path(reanalysis, "logs")
  for (path in c(intermediate, tables, figures, reports, logs_dir)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }

  extracted <- extract_all(project_root)
  all_evidence <- extracted[[1]]
  exclusion_log <- extracted[[2]]
  audits <- extracted[[3]]
  primary <- all_evidence[all_evidence$PrimaryIncluded, , drop = FALSE]
  if (length(intersect(unique(primary$PXD), EXCLUDED_PXDS))) {
    stop("Excluded PXD entered primary evidence")
  }
  missing_pxd <- setdiff(INCLUDED_PXDS, unique(primary$PXD))
  if (length(missing_pxd)) {
    stop(sprintf(
      "Included PXD produced no primary evidence: %s",
      paste(sort(missing_pxd), collapse = ", ")
    ))
  }

  write_csv(all_evidence,
            file.path(intermediate, "kla_by_dataset", "all_included_and_audit_kla_evidence.csv"))
  write_csv(primary,
            file.path(intermediate, "kla_by_dataset", "all_primary_sample_level_kla_sites.csv"))
  write_csv(exclusion_log, file.path(tables, "exclusion_log.csv"))
  for (pxd in INCLUDED_PXDS) {
    dataset <- primary[primary$PXD == pxd, , drop = FALSE]
    write_csv(dataset,
              file.path(intermediate, "kla_by_dataset", paste0(pxd, "_sample_level_kla_sites.csv")))
    write_csv(aggregate_evidence(dataset),
              file.path(tables, "per_pxd", paste0(pxd, "_unique_kla_proteins.csv")))
  }
  lp014 <- suppressWarnings(as.numeric(primary$LocalizationProb))
  pxd014_sensitivity <- primary[
    primary$PXD == "PXD014870" & !is.na(lp014) & lp014 >= 0.75,
    , drop = FALSE
  ]
  write_csv(pxd014_sensitivity,
            file.path(intermediate, "kla_by_dataset", "PXD014870_sensitivity_localization_0.75.csv"))

  directory_audit <- audits$pxd028_directory_audit
  directory_audit$OldWorkflowDirectory <- directory_audit$Directory %in% c(
    "HCT116-Enrichment-Search files",
    "TALL-NALAC-Search files",
    "HEK293T-Enrichment-all HCD-Search files"
  )
  write_csv(directory_audit, file.path(tables, "pxd028488", "directory_audit.csv"))
  write_csv(directory_audit[directory_audit$PrimaryStatus == "included", , drop = FALSE],
            file.path(tables, "pxd028488", "included_directories.csv"))
  write_csv(directory_audit[directory_audit$PrimaryStatus != "included", , drop = FALSE],
            file.path(tables, "pxd028488", "excluded_directories_and_reasons.csv"))
  write_csv(directory_audit,
            file.path(tables, "pxd028488", "old_vs_new_directory_coverage.csv"))
  write_csv(
    all_evidence[
      all_evidence$PXD == "PXD028488" &
        !is.na(all_evidence$DiagnosticIonIntensity),
      , drop = FALSE
    ],
    file.path(tables, "pxd028488", "diagnostic_ion_156_supported_evidence.csv")
  )

  comparison <- audits$pxd053_comparison
  comparison$ComparisonCategory <- comparison$PXD053Comparison
  comparison$ComparisonCategory[comparison$ComparisonCategory == "search_and_supplementary"] <- "consistent"
  write_csv(comparison,
            file.path(tables, "pxd053474", "search_vs_supplementary_all.csv"))
  for (category in c("consistent", "search_only", "supplementary_only")) {
    write_csv(comparison[comparison$ComparisonCategory == category, , drop = FALSE],
              file.path(tables, "pxd053474", paste0(category, ".csv")))
  }
  write_csv(comparison[0, , drop = FALSE],
            file.path(tables, "pxd053474", "inconsistent.csv"))
  write_csv(
    all_evidence[
      all_evidence$PXD == "PXD053474" & !all_evidence$PrimaryIncluded,
      , drop = FALSE
    ],
    file.path(tables, "pxd053474", "single_mode_search_only_audit_evidence.csv")
  )

  go <- read_go_annotations(
    file.path(project_root, "data", "annotations", "GO-repair+damage(human).tsv")
  )
  go_summary <- go$retained
  go_raw <- go$raw
  proteins <- aggregate_evidence(primary)
  proteins_by_pxd <- do.call(rbind, lapply(INCLUDED_PXDS, function(pxd) {
    aggregate_evidence(primary[primary$PXD == pxd, , drop = FALSE])
  }))
  proteins_go <- attach_go(proteins_by_pxd, go_summary, go_raw)
  matched <- proteins_go[proteins_go$GOMatchMode != "unmatched", , drop = FALSE]
  unmatched <- proteins_go[proteins_go$GOMatchMode == "unmatched", , drop = FALSE]
  write_csv(proteins, file.path(tables, "all_unique_kla_proteins.csv"))
  write_csv(matched,
            file.path(intermediate, "go_intersection", "all_pxd_kla_go_ddr_proteins.csv"))
  write_csv(unmatched, file.path(tables, "go_unmatched_kla_proteins.csv"))
  write_csv(
    unmatched[
      unmatched$BaseAccession == "" |
        !nzchar(trimws(ifelse(is.na(unmatched$GeneSymbol), "", unmatched$GeneSymbol))),
      , drop = FALSE
    ],
    file.path(tables, "accession_gene_mapping_failures.csv")
  )
  for (pxd in INCLUDED_PXDS) {
    write_csv(matched[matched$PXD == pxd, , drop = FALSE],
              file.path(intermediate, "go_intersection", paste0(pxd, "_kla_go_ddr.csv")))
  }

  grouping <- read_dtype_str(file.path(reanalysis, "config", "grouping_schemes.csv"))
  review_needed <- grouping[
    tolower(ifelse(is.na(grouping$review_needed), "", grouping$review_needed)) == "yes" |
      ifelse(is.na(grouping$classification_warning), "", grouping$classification_warning) != "",
    , drop = FALSE
  ]
  write_csv(review_needed, file.path(tables, "classification_review_needed.csv"))
  go_accessions <- unique(matched$BaseAccession)
  go_evidence <- primary[primary$BaseAccession %in% go_accessions, , drop = FALSE]
  write_csv(
    cell_type_statistics(primary, go_accessions, grouping),
    file.path(tables, "cell_type_kla_ddr_statistics.csv")
  )
  all_counts <- list()
  region_cache <- list()
  combined_cache <- list()
  non_deduplicated_cache <- list()
  for (scheme in c("teacher_requested_grouping", "biologically_conventional_grouping")) {
    for (analysis_name in c("all_kla", "kla_go_ddr")) {
      evidence <- if (analysis_name == "all_kla") primary else go_evidence
      result <- write_group_outputs(evidence, grouping, scheme, analysis_name,
                                    tables, figures)
      count_frame <- result[[1]]
      count_frame <- cbind(
        data.frame(GroupingScheme = scheme, Analysis = analysis_name,
                   stringsAsFactors = FALSE),
        count_frame
      )
      all_counts[[length(all_counts) + 1]] <- count_frame
      region_cache[[paste(scheme, analysis_name)]] <- result[[2]]
      combined_cache[[paste(scheme, analysis_name)]] <- result[[3]]
      non_deduplicated_cache[[paste(scheme, analysis_name)]] <- result[[4]]
    }
  }
  write_csv(do.call(rbind, all_counts), file.path(tables, "venn_all_schemes_counts.csv"))
  write_csv(combined_cache[["teacher_requested_grouping all_kla"]],
            file.path(tables, "all_kla_three_groups_combined.csv"))
  write_csv(combined_cache[["teacher_requested_grouping kla_go_ddr"]],
            file.path(tables, "kla_go_ddr_three_groups_combined.csv"))
  write_csv(non_deduplicated_cache[["teacher_requested_grouping all_kla"]],
            file.path(tables, "all_kla_three_groups_combined_non_deduplicated.csv"))
  write_csv(non_deduplicated_cache[["teacher_requested_grouping kla_go_ddr"]],
            file.path(tables, "kla_go_ddr_three_groups_combined_non_deduplicated.csv"))

  teacher_all <- region_cache[["teacher_requested_grouping all_kla"]][["tumor_only"]]
  teacher_go <- region_cache[["teacher_requested_grouping kla_go_ddr"]][["tumor_only"]]
  grouping_map <- setNames(grouping$teacher_requested_grouping, grouping$CellType)
  teacher_primary <- primary
  teacher_primary$Category <- unname(grouping_map[teacher_primary$CellOrTissueType])
  teacher_primary$Category[is.na(teacher_primary$Category)] <- ""
  write_csv(
    aggregate_evidence(
      teacher_primary[teacher_primary$BaseAccession %in% teacher_all, , drop = FALSE]
    ),
    file.path(tables, "tumor_specific_kla_proteins.csv")
  )
  write_csv(
    aggregate_evidence(
      teacher_primary[teacher_primary$BaseAccession %in% teacher_go, , drop = FALSE]
    ),
    file.path(tables, "tumor_specific_kla_ddr_proteins.csv")
  )

  trace <- target_trace(all_evidence)
  write_csv(trace,
            file.path(tables, "target_protein_evidence_trace_MRE11_XLF_NBS1.csv"))
  write_csv(
    build_target_source_audit(project_root, primary),
    file.path(tables, "target_protein_source_level_audit_MRE11_XLF_NBS1.csv")
  )
  regression <- regression_outputs(proteins_by_pxd, project_root)
  write_csv(regression[[1]], file.path(tables, "regression_old_vs_new_detail.csv"))
  write_csv(regression[[2]], file.path(tables, "regression_old_vs_new_summary.csv"))

  dataset_config <- read_dtype_str(file.path(reanalysis, "config", "datasets.csv"))
  summary_rows <- list()
  for (i in seq_len(nrow(dataset_config))) {
    metadata <- as.list(dataset_config[i, , drop = FALSE])
    pxd <- metadata$PXD
    subset <- primary[primary$PXD == pxd, , drop = FALSE]
    go_subset <- matched[matched$PXD == pxd, , drop = FALSE]
    site_keys <- unique(paste0(subset$BaseAccession, "\u0001", subset$KlaSite))
    summary_rows[[length(summary_rows) + 1]] <- c(
      metadata,
      list(
        PrimarySampleLevelRows = nrow(subset),
        UniqueKlaSites = length(site_keys),
        UniqueKlaProteins = length(unique(subset$BaseAccession)),
        KlaGoDdrProteins = length(unique(go_subset$BaseAccession))
      )
    )
  }
  summary_frame <- frame_from_rows(
    summary_rows,
    columns = c(names(dataset_config), "PrimarySampleLevelRows",
                "UniqueKlaSites", "UniqueKlaProteins", "KlaGoDdrProteins"),
    integer_cols = c("PrimarySampleLevelRows", "UniqueKlaSites",
                     "UniqueKlaProteins", "KlaGoDdrProteins")
  )
  write_csv(summary_frame, file.path(tables, "dataset_analysis_summary.csv"))

  env_frame <- data.frame(
    Component = c("R", "Platform", "data.table", "ggplot2"),
    Version = c(
      R.version.string,
      R.version$platform,
      as.character(packageVersion("data.table")),
      as.character(packageVersion("ggplot2"))
    ),
    stringsAsFactors = FALSE
  )
  write_csv(env_frame, file.path(logs_dir, "software_environment.csv"))

  cat(sprintf("Primary sample-level rows: %s\n",
              format(nrow(primary), big.mark = ",")))
  cat(sprintf("Unique Kla proteins: %s\n",
              format(length(unique(proteins$BaseAccession)), big.mark = ",")))
  cat(sprintf("Kla GO-DDR proteins: %s\n",
              format(length(unique(matched$BaseAccession)), big.mark = ",")))
  cat("Results: ", tables, "\n", sep = "")
  invisible(0L)
}

# ---------------------------------------------------------------------------
# Stage dispatcher: later tasks append expanded / reference / audit stages
# ---------------------------------------------------------------------------

run_stage <- function(stage, project_root) {
  switch(stage,
    pipeline = run_pipeline_stage(project_root),
    stop(sprintf("Unknown stage: %s", stage))
  )
}

if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  project_root <- if (length(args)) normalizePath(args[[1]]) else normalizePath(".")
  stage <- "pipeline"
  if ("--stage" %in% args) {
    stage <- args[match("--stage", args) + 1]
  }
  invisible(run_stage(stage, project_root))
}
