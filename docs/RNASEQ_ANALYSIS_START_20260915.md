# RNA analysis start: 2026-09-15

Teacher review now governs inclusion: see [review](../audit/20260915_teacher_review/REPORT.md) and `config/rnaseq_teacher_review_gate_20260915.csv`. Existing matrix QC is technical preparation only. DepMap model profiles and the two ENCODE files do not satisfy the preferred three-independent-replicate primary analysis criterion. Old GSE269418 HK-2 references are excluded from the primary candidate pending replacement validation.

Server root: `/home/user/gzy/kla31-rnaseq-20260914`.

Download acceptance at 05:28 UTC: strict FASTQ 12/12 complete, 12 MD5 matches, zero failures; GDC 1017/1017 and DepMap 3/3 verified by download status. This does not promote any biological reference to AnalysisReady.

## Running preparation

- `workflow/start_rnaseq_qc_20260915.R`: source-local DepMap matrices for seven independent ACH models; extract numeric Entrez IDs, exclude all columns belonging to duplicated Entrez IDs with a full audit, retain release units. ENCODE ENCSR000COZ keeps versioned human ENSG rows and separately saves expected_count, TPM and FPKM. Per-sample missingness, negative values, detected genes and within-source Spearman correlations are reported. Inputs receive SHA256 hashes.
- `workflow/check_selected_fastq_integrity_20260915.R`: full gzip integrity traversal of all 12 selected files and structural/name pairing checks of the first 10,000 pairs per run. This is not full adapter, sequence quality or mapping QC.

Outputs: `outputs/20260915_matrix_qc` and `outputs/20260915_fastq_qc` on the server. Logs: `logs/matrix_qc_20260915.log` and `logs/fastq_qc_20260915.log`. Successful stages write `COMPLETE.txt`; absence is not success.

## Remaining analysis

Install isolated quantification tools and a declared reference/index before quantifying the two FASTQ cohorts. Prepare the remaining source-local matrices including recount3 acquisition; resolve GDC sample/diagnosis ambiguity and donor duplicates, GTEx eligibility, GEO column identity and annotation. Keep shared HCT116/HK2 references explicit. Do not directly concatenate source abundance units. Biological eligibility, stable protein-gene mapping and cross-source analysis remain pending.
