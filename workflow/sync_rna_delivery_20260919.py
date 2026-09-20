#!/usr/bin/env python3
"""Explicit, reversible desktop publication of a reviewed RNA release.

Usage: sync_rna_delivery_20260919.py <release_dir> <desktop_kla_dir>
No protein figures or historical source matrices are changed.
"""
import csv
import json
from pathlib import Path
import shutil
import sys

source = Path(sys.argv[1]).resolve()
target = Path(sys.argv[2]).resolve()
repo = Path(__file__).resolve().parents[1]
for name in ("validation.csv", "visual_qa.csv"):
    with (source / name).open() as f:
        rows = list(csv.DictReader(f))
    assert rows and all(r["Pass"].upper() == "TRUE" for r in rows), name
assert (source / "release_sha256.csv").is_file()
assert target.is_dir()
archive = target / "archive" / "rna_before_20260919"
if archive.exists():
    raise RuntimeError(f"Archive already exists; inspect before repeating publication: {archive}")
archive.mkdir(parents=True)
moves = []
seen = set()
for name in ("RNA", "rna"):
    old = target / name
    if not old.exists():
        continue
    ident = (old.stat().st_dev, old.stat().st_ino)
    if ident in seen:
        continue
    seen.add(ident)
    saved = archive / name
    old.rename(saved)
    moves.append({"from": str(old), "to": str(saved)})
top_archive = archive / "top_level"
top_archive.mkdir()
for old in sorted(target.glob("Figure_*RNA*")):
    if old.is_file():
        saved = top_archive / old.name
        old.rename(saved)
        moves.append({"from": str(old), "to": str(saved)})
shutil.copytree(source, target / "RNA")
(target / "RNA" / "SOURCE_LOCATION.txt").write_text(
    f"Authoritative release: {source}\nSHA256 paths in release_sha256.csv are relative to {repo}\n")
for p in sorted((source / "core").glob("*.png")):
    shutil.copy2(p, target / p.name)
old_results = repo / "results" / "rna_teacher_questions"
saved_results = repo / "results" / "archive" / "rna_teacher_questions_before_20260919"
if old_results.exists() or old_results.is_symlink():
    if saved_results.exists():
        raise RuntimeError(f"Local archive already exists: {saved_results}")
    saved_results.parent.mkdir(parents=True, exist_ok=True)
    old_results.rename(saved_results)
    moves.append({"from": str(old_results), "to": str(saved_results)})
old_results.symlink_to(source / "core", target_is_directory=True)
historical = repo / "outputs" / "20260917_delivery_31group" / "rna"
(historical / "HISTORICAL_NOT_CURRENT.md").write_text(
    f"# Historical RNA delivery\n\nThis directory predates the 2026-09-19 correction and contains mixed versions.\n"
    f"Use {source} or config/rna_current_release.csv. Historical files are retained for audit.\n")
log = repo / "audit" / "20260919_rna_repair" / "delivery_sync.json"
log.write_text(json.dumps({"release": str(source), "desktop": str(target / "RNA"),
                           "moves": moves}, indent=2) + "\n")
print(f"PUBLISHED {target / 'RNA'}; preserved {len(moves)} old paths")
