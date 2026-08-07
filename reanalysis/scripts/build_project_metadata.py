#!/usr/bin/env python3
from __future__ import annotations

import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data"
DATASETS_CONFIG = ROOT / "reanalysis" / "config" / "datasets.csv"
PROJECT_INVENTORY = ROOT / "reanalysis" / "reports" / "project_file_inventory.csv"
PRIOR_MANIFEST = ROOT / "archive" / "migration_manifest_2026-07-21.csv"
RECONSTRUCTED_MANIFEST = ROOT / "archive" / "migration_manifest_reconstructed_2026-07-22.csv"

EXPECTED_PXDS = {
    "PXD014870",
    "PXD028488",
    "PXD038880",
    "PXD050470",
    "PXD050906",
    "PXD053474",
    "PXD060185",
    "PXD078013",
    "PXD078736",
}
GENERATED_METADATA_FILES = {"dataset_metadata.csv", "file_inventory.csv"}
EPHEMERAL_SUFFIXES = {".pyc"}
INVENTORY_FIELDS = [
    "PXD",
    "Area",
    "RelativePath",
    "CurrentPath",
    "FileName",
    "Extension",
    "SizeBytes",
]
MANIFEST_FIELDS = [
    "current_path",
    "size_bytes",
    "reconstructed_original_path",
    "recorded_timestamp_utc",
    "recorded_action",
    "recorded_reason",
    "provenance_evidence",
    "reconstruction_confidence",
    "notes",
]


def read_dataset_config() -> tuple[list[str], list[dict[str, str]]]:
    with DATASETS_CONFIG.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        fields = reader.fieldnames
        rows = list(reader)
    if not fields or "PXD" not in fields:
        raise ValueError(f"Missing PXD column in {DATASETS_CONFIG}")

    configured = {row["PXD"] for row in rows}
    if configured != EXPECTED_PXDS or len(rows) != len(EXPECTED_PXDS):
        raise ValueError(
            "datasets.csv must contain exactly the expected nine unique PXD rows; "
            f"found {sorted(configured)}"
        )
    return fields, rows


def write_csv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def inventory_dataset(pxd: str) -> list[dict[str, object]]:
    pxd_dir = DATA / pxd
    if not pxd_dir.is_dir():
        raise FileNotFoundError(f"Missing dataset directory: {pxd_dir}")

    rows: list[dict[str, object]] = []
    paths = sorted(pxd_dir.rglob("*"), key=lambda path: path.as_posix())
    for path in paths:
        if not path.is_file() or path.name == ".DS_Store" or path.suffix in EPHEMERAL_SUFFIXES:
            continue
        if path.parent == pxd_dir / "metadata" and path.name in GENERATED_METADATA_FILES:
            continue

        relative = path.relative_to(pxd_dir)
        area = relative.parts[0] if len(relative.parts) > 1 else "root"
        rows.append(
            {
                "PXD": pxd,
                "Area": area,
                "RelativePath": relative.as_posix(),
                "CurrentPath": path.relative_to(ROOT).as_posix(),
                "FileName": path.name,
                "Extension": path.suffix.lower(),
                "SizeBytes": path.stat().st_size,
            }
        )
    return rows


def load_prior_manifest() -> dict[str, dict[str, str]]:
    with PRIOR_MANIFEST.open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    by_destination: dict[str, dict[str, str]] = {}
    for row in rows:
        destination = row["destination"]
        if destination in by_destination:
            raise ValueError(f"Ambiguous destination in prior manifest: {destination}")
        by_destination[destination] = row
    return by_destination


def exact_manifest_provenance(
    current_path: str, by_destination: dict[str, dict[str, str]]
) -> dict[str, str] | None:
    path = current_path
    chain: list[dict[str, str]] = []
    seen: set[str] = set()
    while path in by_destination:
        if path in seen:
            raise ValueError(f"Migration cycle involving {path}")
        seen.add(path)
        record = by_destination[path]
        chain.append(record)
        path = record["source"]

    if not chain:
        return None

    chronological = list(reversed(chain))
    return {
        "reconstructed_original_path": path,
        "recorded_timestamp_utc": "; ".join(row["timestamp_utc"] for row in chronological),
        "recorded_action": "; ".join(row["action"] for row in chronological),
        "recorded_reason": " | ".join(row["reason"] for row in chronological),
        "provenance_evidence": "archive/migration_manifest_2026-07-21.csv",
        "reconstruction_confidence": "high",
        "notes": "Exact source/destination chain recorded in the prior manifest.",
    }


def inferred_readme_provenance(current_path: str) -> dict[str, str] | None:
    archive_prefix = "archive/2026-07-21_pre_restructure/"
    if current_path.startswith(archive_prefix):
        relative = current_path.removeprefix("archive/")
        return {
            "reconstructed_original_path": f"99_archive/{relative}",
            "recorded_timestamp_utc": "",
            "recorded_action": "inferred_prefix_rename",
            "recorded_reason": "Prior README identifies 99_archive as the historical archive area.",
            "provenance_evidence": (
                "archive/README_before_final_2026-07-21.md; archive/README.md"
            ),
            "reconstruction_confidence": "medium",
            "notes": "Directory-prefix reconstruction from README evidence; not an exact file-level move record.",
        }

    reanalysis_prefix = "archive/reanalysis_v1_2026-07-21/"
    if current_path.startswith(reanalysis_prefix):
        relative = current_path.removeprefix(reanalysis_prefix)
        return {
            "reconstructed_original_path": f"03_reanalysis/{relative}",
            "recorded_timestamp_utc": "",
            "recorded_action": "inferred_archive_move",
            "recorded_reason": "Prior README identifies 03_reanalysis as the first-stage reanalysis area.",
            "provenance_evidence": (
                "archive/README_before_final_2026-07-21.md; "
                "archive/reanalysis_v1_2026-07-21/README.md"
            ),
            "reconstruction_confidence": "medium",
            "notes": "Directory-prefix reconstruction from README evidence; not an exact file-level move record.",
        }
    return None


def reconstruct_migration_manifest() -> None:
    by_destination = load_prior_manifest()
    rows: list[dict[str, object]] = []
    paths = sorted(ROOT.rglob("*"), key=lambda path: path.as_posix())
    for path in paths:
        if (
            not path.is_file()
            or path.name == ".DS_Store"
            or path.suffix in EPHEMERAL_SUFFIXES
            or path == RECONSTRUCTED_MANIFEST
        ):
            continue
        current_path = path.relative_to(ROOT).as_posix()
        provenance = exact_manifest_provenance(current_path, by_destination)
        provenance = provenance or inferred_readme_provenance(current_path)
        if provenance is None:
            provenance = {
                "reconstructed_original_path": "",
                "recorded_timestamp_utc": "",
                "recorded_action": "",
                "recorded_reason": "",
                "provenance_evidence": "Current filesystem snapshot on 2026-07-22",
                "reconstruction_confidence": "unknown",
                "notes": "Current path is known; no supported original path was recovered.",
            }
        rows.append(
            {
                "current_path": current_path,
                "size_bytes": path.stat().st_size,
                **provenance,
            }
        )
    write_csv(RECONSTRUCTED_MANIFEST, MANIFEST_FIELDS, rows)


def main() -> None:
    metadata_fields, metadata_rows = read_dataset_config()
    project_rows: list[dict[str, object]] = []
    for metadata in sorted(metadata_rows, key=lambda row: row["PXD"]):
        pxd = metadata["PXD"]
        metadata_dir = DATA / pxd / "metadata"
        write_csv(metadata_dir / "dataset_metadata.csv", metadata_fields, [metadata])
        inventory_rows = inventory_dataset(pxd)
        write_csv(metadata_dir / "file_inventory.csv", INVENTORY_FIELDS, inventory_rows)
        project_rows.extend(inventory_rows)

    write_csv(PROJECT_INVENTORY, INVENTORY_FIELDS, project_rows)
    reconstruct_migration_manifest()

    print(f"Datasets: {len(metadata_rows)}")
    print(f"Inventoried files: {len(project_rows)}")
    print(f"Project inventory: {PROJECT_INVENTORY.relative_to(ROOT)}")
    print(f"Reconstructed manifest: {RECONSTRUCTED_MANIFEST.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
