#!/usr/bin/env python3
"""Build the seven-study core Kla evidence table used by the R workflow."""

from __future__ import annotations

import argparse
from pathlib import Path
import platform
import sys

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from python.utils.common import apply_annotation_supplement
from python.utils.extractors import (
    extract_maxquant_site_table,
    extract_pxd014870,
    extract_pxd028488,
    extract_pxd050470,
    extract_pxd053474_dda,
    extract_pxd053474_dia,
    extract_pxd053474_supplementary,
    extract_pxd078013,
    reconcile_pxd053474,
)

INCLUDED_PXDS = (
    "PXD028488",
    "PXD050470",
    "PXD053474",
    "PXD060185",
    "PXD078013",
    "PXD078736",
)
EXCLUDED_PXDS = {"PXD038880", "PXD050906"}


def write_csv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, index=False, encoding="utf-8-sig")


def as_bool(series: pd.Series) -> pd.Series:
    if series.dtype == bool:
        return series.fillna(False)
    return (
        series.fillna("")
        .astype(str)
        .str.casefold()
        .isin({"true", "1", "yes", "+"})
    )


def sample_mapping(config: pd.DataFrame, pxd: str) -> dict[str, dict[str, str]]:
    subset = config[config["PXD"].eq(pxd)]
    return {row["SourceToken"]: row.to_dict() for _, row in subset.iterrows()}


def extract_all(project_root: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    data_root = project_root / "data"
    config_root = project_root / "config"
    dataset_config = pd.read_csv(
        config_root / "datasets.csv", dtype=str
    ).fillna("")
    sample_config = pd.read_csv(
        config_root / "sample_map.csv", dtype=str
    ).fillna("")
    doi = dataset_config.set_index("PXD")["DOI"].to_dict()
    frames: list[pd.DataFrame] = []
    logs: list[pd.DataFrame] = []

    sites, log, _ = extract_pxd028488(
        data_root, project_root, doi["PXD028488"]
    )
    frames.append(sites)
    logs.append(log)

    sites, log = extract_pxd050470(
        data_root, project_root, doi["PXD050470"]
    )
    frames.append(sites)
    logs.append(log)

    dda, dda_log = extract_pxd053474_dda(
        data_root, project_root, doi["PXD053474"]
    )
    dia, dia_log = extract_pxd053474_dia(
        data_root, project_root, doi["PXD053474"]
    )
    supplementary, supplementary_log = extract_pxd053474_supplementary(
        data_root, project_root, doi["PXD053474"]
    )
    sites, _ = reconcile_pxd053474(dda, dia, supplementary)
    frames.append(sites)
    logs.extend((dda_log, dia_log, supplementary_log))

    sites, log = extract_maxquant_site_table(
        "PXD060185",
        doi["PXD060185"],
        data_root
        / "PXD060185/search_results/RESULT/combined/txt/La (K)Sites.txt",
        project_root,
        "breast cell line",
        sample_mapping(sample_config, "PXD060185"),
        "maxquant_LaK_site_table",
    )
    frames.append(sites)
    logs.append(log)

    sites, log = extract_pxd078013(
        data_root,
        project_root,
        doi["PXD078013"],
        sample_mapping(sample_config, "PXD078013"),
    )
    frames.append(sites)
    logs.append(log)

    sites, log = extract_maxquant_site_table(
        "PXD078736",
        doi["PXD078736"],
        data_root / "PXD078736/search_results/txt/La(K)Sites.txt",
        project_root,
        "HK-2",
        sample_mapping(sample_config, "PXD078736"),
        "maxquant_LaK_site_table",
    )
    frames.append(sites)
    logs.append(log)

    evidence = pd.concat(frames, ignore_index=True, sort=False)
    evidence = apply_annotation_supplement(
        evidence, config_root / "uniprot_annotation_supplement.tsv"
    )
    evidence["PrimaryIncluded"] = as_bool(evidence["PrimaryIncluded"])
    evidence["ClassI"] = as_bool(evidence["ClassI"])
    evidence = evidence.sort_values(
        [
            "PXD",
            "CellOrTissueType",
            "SampleName",
            "BaseAccession",
            "KlaSite",
            "SourceFile",
            "SourceRow",
        ],
        kind="stable",
    ).reset_index(drop=True)
    exclusion_log = pd.concat(
        [frame for frame in logs if not frame.empty],
        ignore_index=True,
    )
    return evidence, exclusion_log


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Build the stable core Kla evidence inputs."
    )
    parser.add_argument("--project-root", type=Path, default=PROJECT_ROOT)
    args = parser.parse_args()
    project_root = args.project_root.resolve()

    evidence, exclusion_log = extract_all(project_root)
    primary = evidence[evidence["PrimaryIncluded"]].copy()
    if set(primary["PXD"].unique()) & EXCLUDED_PXDS:
        raise AssertionError("Excluded PXD entered primary evidence")
    missing = set(INCLUDED_PXDS) - set(primary["PXD"].unique())
    if missing:
        raise AssertionError(f"Included PXD produced no evidence: {sorted(missing)}")

    intermediate = project_root / "work/intermediate/kla_by_dataset"
    tables = project_root / "results/tables"
    write_csv(
        evidence,
        intermediate / "all_included_and_audit_kla_evidence.csv",
    )
    write_csv(
        primary,
        intermediate / "all_primary_sample_level_kla_sites.csv",
    )
    write_csv(exclusion_log, tables / "core_kla_exclusion_log.csv")

    for pxd in INCLUDED_PXDS:
        write_csv(
            primary[primary["PXD"].eq(pxd)],
            intermediate / f"{pxd}_sample_level_kla_sites.csv",
        )

    environment = pd.DataFrame(
        [
            {"Component": "Python", "Version": platform.python_version()},
            {"Component": "pandas", "Version": pd.__version__},
        ]
    )
    write_csv(environment, project_root / "work/logs/python_environment.csv")
    print(f"Primary sample-level rows: {len(primary):,}")
    print(
        "Unique core Kla proteins: "
        f"{primary['BaseAccession'].nunique():,}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
