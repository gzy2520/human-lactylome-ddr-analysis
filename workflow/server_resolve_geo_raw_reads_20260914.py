#!/usr/bin/env python3
"""Resolve selected GEO accessions to SRA runs and ENA FASTQ URLs on the server.

Python is used only for its standard-library HTTP and JSON support because this server's
R installation does not include jsonlite. This is a download manifest generator, not an
analysis script: it never uses gene symbols or changes any AnalysisReady field.
"""
import csv
import json
import re
import sys
import time
import urllib.parse
import urllib.request
from collections import Counter
from pathlib import Path


def fetch_text(url: str) -> str:
    request = urllib.request.Request(url, headers={"User-Agent": "kla31-rnaseq-manifest/20260914"})
    last_error = None
    for wait in (0, 2, 5, 10):
        if wait:
            time.sleep(wait)
        try:
            with urllib.request.urlopen(request, timeout=60) as response:
                return response.read().decode("utf-8")
        except Exception as error:  # recorded per selected accession below
            last_error = error
    raise RuntimeError(str(last_error))


def expand_selected(value: str) -> list[str]:
    output = []
    for item in value.split(";"):
        item = item.strip()
        match = re.fullmatch(r"(GSM)(\d+)-(GSM)(\d+)", item)
        if match and match.group(1) == match.group(3):
            output.extend(f"{match.group(1)}{number}" for number in range(int(match.group(2)), int(match.group(4)) + 1))
        else:
            output.append(item)
    return output


def sra_ids(accession: str) -> list[str]:
    term = urllib.parse.quote(f"{accession}[All Fields]", safe="")
    answer = json.loads(fetch_text(
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&retmode=json&retmax=100&term=" + term
    ))
    return answer.get("esearchresult", {}).get("idlist", [])


def runinfo(ids: list[str]) -> list[dict]:
    if not ids:
        return []
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?" + urllib.parse.urlencode({
        "db": "sra", "rettype": "runinfo", "retmode": "text", "id": ",".join(ids)
    })
    return list(csv.DictReader(fetch_text(url).splitlines()))


def ena_fastq(run: str) -> list[dict]:
    url = "https://www.ebi.ac.uk/ena/portal/api/filereport?" + urllib.parse.urlencode({
        "accession": run, "result": "read_run",
        "fields": "run_accession,experiment_accession,sample_accession,fastq_ftp,fastq_md5,fastq_bytes,library_layout",
        "format": "tsv", "download": "true"
    })
    return list(csv.DictReader(fetch_text(url).splitlines(), delimiter="\t"))


def unresolved(task: dict, accession: str, state: str, detail: str, ids: list[str] | None = None, run: dict | None = None) -> dict:
    run = run or {}
    return {
        "TaskID": task["TaskID"], "GroupIDs": task["GroupIDs"], "GSE": task["Accession"],
        "SelectedAccession": accession, "SRAQueryIDs": ";".join(ids or []), "SRR": run.get("Run", ""),
        "Experiment": run.get("Experiment", ""), "RunInfoSample": run.get("SampleName", ""),
        "LibraryName": run.get("LibraryName", ""), "LibraryLayout": run.get("LibraryLayout", ""),
        "ENAFastqURL": "", "ExpectedMD5": "", "ExpectedBytes": "", "ResolutionState": state, "Detail": detail
    }


def main(specification: Path, output: Path) -> None:
    with specification.open(newline="") as handle:
        tasks = [row for row in csv.DictReader(handle, delimiter="\t") if row["Source"] == "GEO"]
    rows: list[dict] = []
    for task in tasks:
        for accession in expand_selected(task["SelectedSamples"]):
            print(f"RESOLVING {task['TaskID']} {task['Accession']} {accession}", flush=True)
            try:
                ids = sra_ids(accession)
                time.sleep(0.4)
                runs = runinfo(ids)
            except Exception as error:
                rows.append(unresolved(task, accession, "query_failed", f"SRA query failed: {error}"))
                continue
            if not runs:
                rows.append(unresolved(task, accession, "unresolved", "No SRA RunInfo record returned for selected accession.", ids))
                continue
            for run in runs:
                time.sleep(0.4)
                try:
                    ena = ena_fastq(run["Run"])
                except Exception as error:
                    rows.append(unresolved(task, accession, "ena_query_failed", f"ENA query failed: {error}", ids, run))
                    continue
                if not ena or not ena[0].get("fastq_ftp"):
                    rows.append(unresolved(task, accession, "no_ena_fastq", "SRA run resolved but ENA did not return FASTQ metadata.", ids, run))
                    continue
                record = ena[0]
                urls = record["fastq_ftp"].split(";")
                md5s = record["fastq_md5"].split(";")
                sizes = record["fastq_bytes"].split(";")
                if not (len(urls) == len(md5s) == len(sizes)):
                    raise RuntimeError(f"ENA FASTQ fields have unequal lengths for {run['Run']}")
                for url, md5, size in zip(urls, md5s, sizes):
                    rows.append({
                        "TaskID": task["TaskID"], "GroupIDs": task["GroupIDs"], "GSE": task["Accession"],
                        "SelectedAccession": accession, "SRAQueryIDs": ";".join(ids), "SRR": run["Run"],
                        "Experiment": run.get("Experiment", ""), "RunInfoSample": run.get("SampleName", ""),
                        "LibraryName": run.get("LibraryName", ""), "LibraryLayout": record.get("library_layout") or run.get("LibraryLayout", ""),
                        "ENAFastqURL": "https://" + url, "ExpectedMD5": md5.lower(), "ExpectedBytes": size,
                        "ResolutionState": "resolved", "Detail": "SRA search was initiated with the selected GEO/SRA accession."
                    })
            print(f"RESOLVED_ROWS {task['TaskID']} {accession} total={len(rows)}", flush=True)
    fields = ["TaskID", "GroupIDs", "GSE", "SelectedAccession", "SRAQueryIDs", "SRR", "Experiment", "RunInfoSample", "LibraryName", "LibraryLayout", "ENAFastqURL", "ExpectedMD5", "ExpectedBytes", "ResolutionState", "Detail"]
    resolved = [row for row in rows if row["ResolutionState"] == "resolved"]
    if any(not re.fullmatch(r"SRR\d+", row["SRR"]) for row in resolved):
        raise RuntimeError("Resolved row contains an invalid SRR accession")
    if any(not re.fullmatch(r"[0-9a-f]{32}", row["ExpectedMD5"]) for row in resolved):
        raise RuntimeError("Resolved row contains an invalid ENA MD5")
    # A technical replicate or shared sample record can legitimately return the same
    # FASTQ for more than one selected accession. Download it once while retaining
    # every selected accession in the merged provenance fields.
    merged = {}
    for row in resolved:
        key = row["ENAFastqURL"]
        if key not in merged:
            merged[key] = row.copy()
            continue
        prior = merged[key]
        for field in ("TaskID", "GroupIDs", "SelectedAccession", "SRAQueryIDs"):
            values = [item for item in (prior[field] + ";" + row[field]).split(";") if item]
            prior[field] = ";".join(dict.fromkeys(values))
        prior["Detail"] = prior["Detail"] + " Duplicate URL relation merged with another selected accession."
    unresolved_rows = [row for row in rows if row["ResolutionState"] != "resolved"]
    rows = unresolved_rows + list(merged.values())
    resolved = list(merged.values())
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + ".tmp")
    with temporary.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", quoting=csv.QUOTE_ALL)
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(output)
    summary = Counter((row["TaskID"], row["GSE"], row["ResolutionState"]) for row in rows)
    with output.with_name(output.stem + "_summary.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", quoting=csv.QUOTE_ALL)
        writer.writerow(["TaskID", "GSE", "ResolutionState", "RowCount"])
        for (task_id, gse, state), count in sorted(summary.items()):
            writer.writerow([task_id, gse, state, count])
    print(f"GEO_RAW_RESOLUTION_PASS: {len(resolved)} FASTQ files resolved; {len(rows) - len(resolved)} non-resolved rows retained explicitly.")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit("Usage: server_resolve_geo_raw_reads_20260914.py <download_spec_tsv> <output_tsv>")
    main(Path(sys.argv[1]), Path(sys.argv[2]))
