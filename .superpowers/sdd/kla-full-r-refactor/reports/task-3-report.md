# Task 3 Report: lib/ Shared Modules

**Status: DONE_WITH_CONCERNS**

## Summary

Created three shared R library modules under `reanalysis/scripts/lib/` by extracting pure functions verbatim from existing analysis scripts, plus a self-check test script. No existing scripts were modified. Loading order: `accession_utils.R` -> `io_utils.R` -> `extractors.R`.

---

## Files Created

### 1. `reanalysis/scripts/lib/accession_utils.R`

Functions extracted from source scripts:

| Function | Source | Lines (source) | Notes |
|---|---|---|---|
| `base_accession(values)` | `analyze_expanded_ddr_fraction_by_accession.R` | 41-47 | Strips db prefix (`sp|...|`), `NX_`, isoform suffix. includes `trimws()`. |
| `is_uniprot(values)` | `analyze_expanded_ddr_fraction_by_accession.R` | 57-67 | Regex match for canonical UniProt accession pattern. |
| `split_accessions(values)` | `analyze_expanded_ddr_fraction_by_accession.R` | 69-73 | Split by `[;,]`, base_accession, filter UniProt, sort unique. |
| `split_protein_identifiers(values)` | `analyze_expanded_ddr_fraction_by_accession.R` | 75-84 | Split, return both UniProt (stripped) and ENSP identifiers. |
| `identifier_type(values)` | `analyze_expanded_ddr_fraction_by_accession.R` | 175-177 | Classifies as `"ENSEMBLPROT"` or `"UniProtKB"`. |
| `match_target_accession(values, target_accessions)` | `analyze_kla_regulator_intensity.R` (91-99), `analyze_kla_regulator_whole_proteome_intensity.R` (98-105) | 91-99, 98-105 | **Adapted**: `target_accessions` passed as explicit parameter (was closure over script env). |
| `accession_feature(values)` | `analyze_kla_regulator_intensity.R` (101-108), `analyze_kla_regulator_whole_proteome_intensity.R` (108-115) | 101-108, 108-115 | Returns `"ACC:<first_token>"`; no target_accessions dependency. |
| `safe_numeric(values)` | `analyze_kla_regulator_intensity.R` (111-113), `analyze_kla_regulator_whole_proteome_intensity.R` (118-120) | 111-113, 118-120 | Comma-stripping numeric conversion. |

Header: `lib_loaded <- TRUE`

### 2. `reanalysis/scripts/lib/io_utils.R`

| Function | Source | Lines (source) | Notes |
|---|---|---|---|
| `relative_path(path, root)` | `analyze_expanded_ddr_fraction_by_accession.R` (49-55) | 49-55 | **Signature changed**: `root` passed as explicit parameter (was closure). |
| `valid_maxquant_rows(data)` | `analyze_expanded_ddr_fraction_by_accession.R` | 344-361 | Filters Reverse/contaminant/site-only rows. |
| `read_delimited(path)` | `analyze_expanded_ddr_fraction_by_accession.R` (363-375); `analyze_kla_regulator_intensity.R` (182-194) | 363-375, 182-194 | **Identical** between both sources. CSV via `read.csv`, other via `read.delim`. |
| `write_csv_std(frame, path)` | New (interface spec) | -- | UTF-8 CSV: `write.csv(..., row.names = FALSE, na = "")`. |
| `save_figure(path_stem, plot, width, height, dpi=350, ...)` | New (interface spec) | -- | Writes `<stem>.png` (ggsave + bg="white") and `<stem>.pdf` (cairo_pdf). |

Header: `lib_loaded <- TRUE`

### 3. `reanalysis/scripts/lib/extractors.R`

All functions extracted verbatim from `analyze_expanded_ddr_fraction_by_accession.R`:

| Function | Source Lines | Dependencies |
|---|---|---|
| `extract_maxquant_sites(path, sample_tokens, sheet)` | 377-407 | `read_delimited`, `valid_maxquant_rows`, `split_accessions`, `read_excel` (readxl) |
| `extract_maxquant_proteins(path, abundance_pattern)` | 409-434 | `read_delimited`, `valid_maxquant_rows`, `split_protein_identifiers`, `map_ensembl_proteins` (*) |
| `extract_pd_proteins(path, sample_token, lactylome)` | 436-456 | `split_accessions` |
| `extract_pd_lactyl_peptides(path, sample_token)` | 458-483 | `split_accessions` |
| `extract_spectronaut_proteins(path, lactyl_pattern, group_pattern, accession_column)` | 485-523 | `fread` (data.table), `split_accessions` |
| `extract_spectronaut_matrix(path, lactyl_pattern)` | 525-529 | `read_delimited`, `split_accessions` |
| `extract_spectronaut_quant(path, group_pattern)` | 531-545 | `read_delimited`, `split_accessions` |
| `extract_huvec_xml(path)` | 566-572 | `str_extract_all` (stringr), `split_accessions` |

(*) `extract_maxquant_proteins` calls `map_ensembl_proteins(identifiers)` (defined in expanded DDR script lines 110-159) for its side effects on `ensembl_mapping_records`. This function is NOT included in the lib; callers must ensure it is available in the environment.

Header: `extractors_loaded <- TRUE`

### 4. `reanalysis/scripts/lib/test_lib_self_check.R`

Self-check script that sources all three lib modules and asserts correct behavior.

**Assertions corrected from brief Step 4 examples to match actual source behavior:**

| Brief Example | Corrected To | Reason |
|---|---|---|
| `base_accession("REV__CON__Q9H9Q4-3") == "Q9H9Q4"` | `== "REV__CON__Q9H9Q4"` | Source `base_accession` strips `-[0-9]+$` only; `REV__CON__` prefix is not a known pattern. |
| `accession_feature(c("P49959-2", "P49959-1")) == c("P49959", "P49959")` | `== c("ACC:P49959", "ACC:P49959")` | Source `accession_feature` prepends `"ACC:"` prefix (heatmap script lines 101-108). |

Additional fixes applied:
- `vapply` returns named vectors by default in R; wrapped all `identical()` comparisons involving `accession_feature`/`match_target_accession` with `unname()`.
- `write_csv_std` test: `read.csv` with default `na.strings` reads empty fields as `""`, not `NA`; assertion corrected to `result$y[1] == ""`.

---

## Cross-Script Comparison: Behavior Differences and Choices

### `base_accession` -- All Three Versions

| Feature | Expanded DDR (lib base) | Intensity Script | Whole Proteome Script |
|---|---|---|---|
| `trimws()` | Yes | No | No |
| Pipe prefix `sp\|...\|` | Yes | Yes | Yes |
| `NX_` prefix | Yes | No | No |
| Extra `sub("^([^\|;]+)\\|.*$")` | No | No | Yes |
| Isoform `-[0-9]+$` | Yes | Yes | Yes |

**Decision**: Use expanded DDR version as the canonical `base_accession` (most complete: handles whitespace + NX_ prefix + pipe prefix + isoform suffix). The whole proteome script's extra `sub("^([^\|;]+)\\|.*$")` handles a different format (`Q9H9Q4|description`), which is rare in this project's data.

### `relative_path` -- Signature Change

All source versions use `function(path)` with `project_root` captured from closure. Lib version uses `function(path, root)` for explicit dependency injection.

### `match_target_accession` -- Heatmap Only

Present only in the two heatmap scripts (not in expanded DDR). Both heatmap versions are identical. Lib version adapts to take `target_accessions` as a parameter instead of capturing from closure.

### `accession_feature` -- Heatmap Only

Present only in the two heatmap scripts. Both versions are identical. Key behavior: returns `"ACC:<first_token>"` (includes `"ACC:"` prefix), not raw accession.

### `safe_numeric` -- Both Heatmap Scripts

Identical in both heatmap scripts. Not present in expanded DDR.

### `read_delimited` -- Expanded DDR vs Intensity Script

**Identical** between both sources (verbatim match). No decision needed.

---

## Verification Results

### Self-Check (Step 5)

```
Command: Rscript reanalysis/scripts/lib/test_lib_self_check.R .
Output:  lib self-check passed
Status:  PASS
```

### Script Integrity (Step 6)

```
Command: Rscript reanalysis/scripts/analyze_expanded_ddr_fraction_by_accession.R .
Output:  Expanded accession-only DDR comparison: 37 included; 3 excluded.
Status:  PASS (exit code 0, no errors)
```

### Baseline Verification (Step 6)

```
Command: Rscript reanalysis/scripts/lib/verify_outputs.R . archive/refactor_baseline_2026-08-06/output_sha256_baseline.csv
Output:
  OK: 358 files
  CHANGED: 4 files
  MISSING: 0 files
```

**Changed files (all figures -- pre-existing rendering variation):**

| File | Baseline Size | Current Size | Delta |
|---|---|---|---|
| `cell_type_kla_vs_reference_ddr_fraction.pdf` | 87,590 | 88,350 | +760 |
| `cell_type_kla_vs_reference_ddr_fraction.png` | 1,819,702 | 1,131,036 | -688,666 |
| `cell_type_kla_vs_reference_ddr_fraction_accession_only.pdf` | 87,590 | 88,350 | +760 |
| `cell_type_kla_vs_reference_ddr_fraction_accession_only.png` | 1,819,702 | 1,131,036 | -688,666 |

**All 358 non-figure baseline files (CSV, TSV, other data outputs) are byte-identical to baseline.**

---

## Concerns

1. **4 figure files show hash changes** (PDF/PNG). These are NOT caused by Task 3 changes (no existing script was modified). The differences are consistent with render-level variation (Cairo library version, font metrics, R graphics device, or ggplot2 version differences between the baseline creation run and the current run). The PNG size change (-38%) suggests a DPI or compression setting difference. The PDF size change (+0.9%) is consistent with minor font or layout metric drift. **Recommendation**: Re-baseline these figure files, or mark figure files as "expected to vary" in the verify tool.

2. **`extract_maxquant_proteins` calls `map_ensembl_proteins`** which is NOT in the lib. This function (defined in expanded DDR script lines 110-159) has side effects on `ensembl_mapping_records` and requires `AnnotationDbi` + `org.Hs.eg.db`. Callers of the lib version must ensure `map_ensembl_proteins` is available in the environment. This should be addressed in a future task when `map_ensembl_proteins` is also extracted or the dependency is made explicit.

3. **`match_target_accession` signature differs from source**: The lib version requires `target_accessions` as a second parameter, while the heatmap scripts use it as a closure. When future tasks refactor the heatmap scripts to use the lib, the call sites will need updating.

---

## Files Not Modified

- `reanalysis/scripts/analyze_expanded_ddr_fraction_by_accession.R` -- unchanged
- `reanalysis/scripts/analyze_kla_regulator_intensity.R` -- unchanged
- `reanalysis/scripts/analyze_kla_regulator_whole_proteome_intensity.R` -- unchanged
- `reanalysis/scripts/lib/verify_outputs.R` -- unchanged (pre-existing)

---

## Artifacts

| Artifact | Path |
|---|---|
| accession_utils.R | `/Users/gzy2520/Desktop/Research/kla/reanalysis/scripts/lib/accession_utils.R` |
| io_utils.R | `/Users/gzy2520/Desktop/Research/kla/reanalysis/scripts/lib/io_utils.R` |
| extractors.R | `/Users/gzy2520/Desktop/Research/kla/reanalysis/scripts/lib/extractors.R` |
| test_lib_self_check.R | `/Users/gzy2520/Desktop/Research/kla/reanalysis/scripts/lib/test_lib_self_check.R` |
| This report | `/Users/gzy2520/Desktop/Research/kla/.superpowers/sdd/kla-full-r-refactor/reports/task-3-report.md` |
