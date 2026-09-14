# GEO / GTEx / ENCODE server-download snapshot

`download_status.tsv` is the append-only server history. `download_status_latest.tsv` keeps the latest state for each target path; `GSE163787_series_matrix.txt.gz` is marked as superseded because GEO publishes platform-specific matrix files, which were recovered as the GPL20301 and GPL21103 files. The selected TALL-104 sample is in the GPL20301 matrix; the second platform file remains source metadata for species/platform exclusion.

The GEO raw-read manifest and its ENA download queue are separate. Raw files remain on the server and are not copied into Git.
