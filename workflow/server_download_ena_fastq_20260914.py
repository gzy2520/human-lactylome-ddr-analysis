#!/usr/bin/env python3
"""Download the exact ENA FASTQ files in a resolved GEO manifest.

Only rows with ResolutionState=resolved are eligible. The manifest is treated as an
external candidate contract; no sample is promoted to AnalysisReady by this script.
"""
import csv
import hashlib
import os
import subprocess
import sys
import time
from pathlib import Path
from urllib.parse import urlparse


def checksum(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def append_status(path: Path, row: dict, target: str, state: str, file: Path | None, detail: str) -> None:
    exists = path.exists()
    observed = checksum(file) if file and file.exists() and state == "complete" else ""
    size = str(file.stat().st_size) if file and file.exists() else ""
    record = {
        "TimestampUTC": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "TaskID": row["TaskID"], "GroupIDs": row["GroupIDs"], "GSE": row["GSE"], "SelectedAccession": row["SelectedAccession"],
        "SRR": row["SRR"], "ENAFastqURL": row["ENAFastqURL"], "RelativeTarget": target, "State": state,
        "Bytes": size, "ExpectedBytes": row["ExpectedBytes"], "ExpectedMD5": row["ExpectedMD5"], "ObservedMD5": observed, "Detail": detail,
    }
    with path.open("a", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(record), delimiter="\t", quoting=csv.QUOTE_ALL)
        if not exists:
            writer.writeheader()
        writer.writerow(record)


def main(root: Path, manifest_path: Path) -> None:
    with manifest_path.open(newline="") as handle:
        rows = [row for row in csv.DictReader(handle, delimiter="\t") if row["ResolutionState"] == "resolved"]
    if not rows:
        raise RuntimeError("The manifest contains no resolved FASTQ rows; refusing to start a download.")
    urls = [row["ENAFastqURL"] for row in rows]
    if len(urls) != len(set(urls)):
        raise RuntimeError("Resolved manifest contains duplicate ENA URLs.")
    out_root = root / "raw" / "ena_fastq"
    status_path = root / "metadata" / "ena_fastq_download_status.tsv"
    out_root.mkdir(parents=True, exist_ok=True)
    for index, row in enumerate(rows, start=1):
        url = row["ENAFastqURL"]
        name = Path(urlparse(url).path).name
        if not name:
            raise RuntimeError(f"Could not derive a file name from {url}")
        target = Path("ena_fastq") / row["GSE"] / f"{row['SRR']}_{name}"
        destination = root / "raw" / target
        destination.parent.mkdir(parents=True, exist_ok=True)
        if destination.exists() and checksum(destination) == row["ExpectedMD5"].lower():
            append_status(status_path, row, str(target), "complete", destination, "already present; ENA MD5 verified")
            continue
        if destination.exists():
            destination.unlink()
        partial = Path(str(destination) + ".part")
        expected_md5 = row["ExpectedMD5"].lower()
        expected_size = int(row["ExpectedBytes"]) if row["ExpectedBytes"].isdigit() else None
        if (
            partial.exists()
            and expected_size is not None
            and partial.stat().st_size >= expected_size
            and checksum(partial) != expected_md5
        ):
            quarantined = Path(f"{partial}.checksum-mismatch-{time.strftime('%Y%m%dT%H%M%SZ', time.gmtime())}")
            partial.replace(quarantined)
        command = ["wget", "--continue", "--timeout=60", "--read-timeout=60", "--tries=5", "--waitretry=5", "--retry-connrefused", "--no-verbose", "--show-progress", "-O", str(partial), url]
        code = subprocess.run(command, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT, check=False).returncode
        if code == 0 and partial.exists() and checksum(partial) == expected_md5:
            if expected_size is not None and partial.stat().st_size != expected_size:
                append_status(status_path, row, str(target), "failed", partial, "MD5 matched but ENA byte count differed; retained .part")
                continue
            partial.replace(destination)
            append_status(status_path, row, str(target), "complete", destination, "downloaded; ENA MD5 and byte count verified")
        else:
            append_status(status_path, row, str(target), "failed", partial if partial.exists() else None, f"wget exit {code}; retained .part for resumption or MD5 mismatch")
        if index % 10 == 0:
            print(f"ENA_FASTQ_PROGRESS {index}/{len(rows)}", flush=True)
    (root / "metadata" / "ena_fastq_README.txt").write_text(
        "Only resolved GEO/SRA accessions were downloaded. Each complete file has the ENA MD5 and byte count recorded in ena_fastq_download_status.tsv.\n"
        "The files remain candidate raw material and do not make any registered group AnalysisReady.\n",
        encoding="utf-8",
    )
    print(f"ENA_FASTQ_DOWNLOAD_STAGE_COMPLETE: {len(rows)} resolved FASTQ files inspected.")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit("Usage: server_download_ena_fastq_20260914.py <server_root> <resolved_manifest_tsv>")
    main(Path(sys.argv[1]), Path(sys.argv[2]))
