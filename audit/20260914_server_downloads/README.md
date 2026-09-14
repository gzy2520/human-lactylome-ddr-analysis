# 2026-09-14 server RNA download snapshots

This directory contains small provenance and status snapshots copied from the isolated server job at `192.168.3.45:/home/user/gzy/kla31-rnaseq-20260914`. Raw RNA files remain on the server and are intentionally not copied into Git.

The `depmap/` snapshot records the fixed DepMap Public 24Q4 (2024-12-16) manifest, the official no-captcha release catalog, the completed MD5 status rows, and the structure check. The expression matrix, `Model.csv`, and `OmicsProfiles.csv` all passed their published MD5 checks on the server. The structure check found all seven target ACH rows and 17 duplicate Entrez header IDs; duplicate resolution remains a downstream QC step.

To inspect the live queues from macOS, run `workflow/check_server_rnaseq_download_20260914.sh`; use `workflow/check_server_rnaseq_download_20260914.sh --watch 60` for a 60-second refresh interval.
