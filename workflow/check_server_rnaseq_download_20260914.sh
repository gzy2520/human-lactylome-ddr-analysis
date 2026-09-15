#!/usr/bin/env bash
# Inspect the isolated server-side RNA download queues without exposing signed URLs.
set -euo pipefail

remote_host="${REMOTE_HOST:-user@100.121.229.123}"
watch_seconds=""
if [[ "${1:-}" == "--watch" ]]; then
  watch_seconds="${2:-60}"
  if ! [[ "$watch_seconds" =~ ^[1-9][0-9]*$ ]]; then
    echo "Usage: $0 [--watch SECONDS]" >&2
    exit 2
  fi
elif [[ -n "${1:-}" ]]; then
  echo "Usage: $0 [--watch SECONDS]" >&2
  exit 2
fi

run_report() {
  ssh -T -o BatchMode=yes -o ConnectTimeout=10 -o HostKeyAlias=192.168.3.45 "$remote_host" 'bash -s' <<'REMOTE_REPORT'
set -u
server_root="/home/user/gzy/kla31-rnaseq-20260914"
cd "$server_root"

python3 - <<'PY'
from collections import Counter
from pathlib import Path
import csv
import re
import time

root = Path("/home/user/gzy/kla31-rnaseq-20260914")
print("server", root)
print("checked_utc", time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()))

jobs = [
    ("processed_geo_gtex_encode", "download_20260914"),
    ("gdc_star_counts", "gdc_download_20260914"),
    ("geo_raw_resolution", "geo_raw_resolution_20260914"),
    ("ena_fastq", "ena_fastq_download_20260914"),
    ("depmap", "depmap_download_20260914"),
    ("gse163787_recovery", "gse163787_recovery_20260914"),
]
print("jobs")
for label, name in jobs:
    pid_file = root / "logs" / f"{name}.pid"
    pid = pid_file.read_text().strip() if pid_file.exists() else ""
    running = bool(pid) and Path(f"/proc/{pid}").exists()
    print(f"  {label}: {'RUNNING' if running else 'STOPPED'} pid={pid or '-'}")

def rows_from(relative):
    path = root / relative
    if not path.exists():
        return None
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))

manifest = rows_from("config/gdc_star_counts_manifest_20260914.tsv")
gdc = rows_from("metadata/gdc_star_counts_download_status.tsv")
if manifest is not None:
    complete = [r for r in (gdc or []) if r.get("State") == "complete"]
    failed = [r for r in (gdc or []) if r.get("State") == "failed"]
    matches = [r for r in complete if r.get("ExpectedMD5") == r.get("ObservedMD5")]
    print(f"gdc: manifest={len(manifest)} status={len(gdc or [])} complete={len(complete)} failed={len(failed)} md5_match={len(matches)}")

depmap = rows_from("metadata/depmap_download_status.tsv")
if depmap is not None:
    complete = [r for r in depmap if r.get("State") == "complete"]
    failed = [r for r in depmap if r.get("State") == "failed"]
    matches = [r for r in complete if r.get("ExpectedMD5") == r.get("ObservedMD5")]
    print(f"depmap: manifest=3 status={len(depmap)} complete={len(complete)} failed={len(failed)} md5_match={len(matches)}")

geo = rows_from("metadata/download_status.tsv")
if geo is not None:
    latest = {row.get("RelativeTarget", ""): row for row in geo}
    for target, row in latest.items():
        if row.get("State") == "failed" and target.endswith("GSE163787_series_matrix.txt.gz"):
            row["State"] = "superseded_by_platform_specific_matrices"
    counts = Counter(row.get("State", "") for row in latest.values())
    print("processed_geo_gtex_encode_latest:", " ".join(f"{key}={value}" for key, value in sorted(counts.items())))

raw_manifest_path = "metadata/ena_fastq_active_manifest.tsv"
raw_manifest = rows_from(raw_manifest_path)
if raw_manifest is None:
    raw_manifest_path = "metadata/geo_selected_raw_fastq_manifest_20260914.tsv"
    raw_manifest = rows_from(raw_manifest_path)
if raw_manifest is not None:
    resolved = [r for r in raw_manifest if r.get("ResolutionState") == "resolved"]
    print(f"ena_manifest: file={raw_manifest_path} rows={len(raw_manifest)} resolved={len(resolved)}")
ena = rows_from("metadata/ena_fastq_download_status.tsv")
if ena is not None:
    latest = {r.get("RelativeTarget", ""): r for r in ena}
    complete = [r for r in latest.values() if r.get("State") == "complete"]
    failed = [r for r in latest.values() if r.get("State") == "failed"]
    matches = [r for r in complete if r.get("ExpectedMD5") == r.get("ObservedMD5")]
    print(f"ena_status: records={len(ena)} unique_files={len(latest)} complete={len(complete)} failed={len(failed)} md5_match={len(matches)}")
else:
    print("ena_status: not created yet (the current file may still be downloading)")

def directory_size(relative):
    directory = root / relative
    file_count = 0
    byte_count = 0
    if directory.exists():
        for item in directory.rglob("*"):
            if item.is_file():
                file_count += 1
                byte_count += item.stat().st_size
    print(f"size {relative}: files={file_count} bytes={byte_count}")

for relative in ["raw/gdc_star_counts", "raw/ena_fastq", "raw/depmap", "raw/geo", "raw/gtex_v8", "raw/encode"]:
    directory_size(relative)

print("active_partial_files")
for partial in sorted((item for item in root.glob("raw/**/*.part") if item.stat().st_size > 0), key=lambda item: item.stat().st_size, reverse=True)[:5]:
    print(f"  {partial.relative_to(root)} bytes={partial.stat().st_size}")

print("recent_progress")
for name in ["download_20260914", "gdc_download_20260914", "ena_fastq_download_20260914", "depmap_download_20260914"]:
    log = root / "logs" / f"{name}.log"
    if not log.exists():
        continue
    lines = log.read_text(errors="replace").splitlines()
    meaningful = [line for line in lines if line.strip()]
    if meaningful:
        last = re.sub(r"https?://\\S+", "<url>", meaningful[-1])
        print(f"  {name}: {last[:240]}")
PY
REMOTE_REPORT
}

if [[ -n "$watch_seconds" ]]; then
  while true; do
    run_report
    echo "next_check_seconds=$watch_seconds"
    sleep "$watch_seconds"
  done
else
  run_report
fi
