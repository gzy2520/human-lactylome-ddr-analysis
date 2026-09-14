# 2026-09-14 RNA source-routing server snapshot

`source_native_archive_probe_status.tsv` records the completed GEO archive probes that changed the routing decision:

- GSE203529 contains three WT untreated HEK293T RSEM gene-result files.
- GSE132714 contains the 18 selected BPH per-sample, versioned-Ensembl RPKM files.

`ena_fastq_active_manifest.tsv` is the deliberately reduced FASTQ scope. It contains only 12 MD5-resolved files: six for KLA31_01 pathological rotator-cuff tendon and six for KLA31_22 exact PC-3M DMSO controls. It supersedes the earlier 137-file exploratory ENA queue while retaining already completed files and partial files for continuation.

`ena_fastq_download_status.tsv` is an append-only status snapshot. Use `workflow/check_server_rnaseq_download_20260914.sh` for the current live state; it deduplicates repeated status records by target path and reports the active manifest rather than the superseded exploratory scope.

These are acquisition/provenance records. They do not make any group analysis-ready and contain no derived expression matrix.
