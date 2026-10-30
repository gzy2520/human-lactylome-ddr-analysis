#!/usr/bin/env python3
"""Extract the archived point values used by the supplied figure screenshots."""
import argparse
import csv
import hashlib
import json
import math
import shutil
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CATEGORIES = ["normal_tissue", "cancer_tissue", "normal_cells", "cancer_cells"]
LABELS = ["non-tumor tissues", "tumor tissues", "normal cell lines", "cancer cell lines"]
PATHWAYS = [("d", "BER"), ("e", "NER"), ("f", "MMR"), ("g", "FA"),
            ("h", "HR"), ("i", "AEJ"), ("j", "NHEJ")]


def column_name(number):
    label = ""
    while number:
        number, remainder = divmod(number - 1, 26)
        label = chr(65 + remainder) + label
    return label


def load(path):
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def validate_counts(rows, expected):
    observed = Counter(row["Category"] for row in rows)
    assert [observed[category] for category in CATEGORIES] == expected, observed


def ordered(rows):
    return sorted(rows, key=lambda row: (int(row["RowOrder"]), row.get("PXD", ""),
                                        row.get("SampleGroup", ""), row.get("SampleID", "")))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-root", type=Path, default=ROOT / ".worktrees/kla31-rnaseq-20260912")
    parser.add_argument("--input-dir", type=Path, help="Read the three retained source CSVs from this directory.")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "outputs/01a18b74-d839-7a53-b3c3-3c45f301073e")
    args = parser.parse_args()
    frozen = args.source_root / "corrected_final_result_20260905_sample_only/Corrected_Data/sample_only_inputs"
    sources = {
        "fraction": frozen / "figure1_sample_boxplot_values.csv",
        "pathways": frozen / "figure1_pathway_summary_sample_boxplot_values.csv",
        "proliferation": args.source_root / "outputs/20260917_delivery_31group/rna/figure1b_hallmark_proliferation_31groups.csv",
    }
    if args.input_dir:
        sources = {key: args.input_dir / path.name for key, path in sources.items()}
    fraction, pathways, proliferation = [load(sources[key]) for key in ("fraction", "pathways", "proliferation")]
    assert len(fraction) == 272 and all(row["ObservationType"] == "sample" for row in fraction)
    assert len({tuple(row[key] for key in ("PXD", "SampleGroup", "Dataset", "SampleID")) for row in fraction}) == 272
    validate_counts([row for row in fraction if row["Dataset"] == "Whole proteome"], [60, 107, 13, 20])
    validate_counts([row for row in fraction if row["Dataset"] == "Lactylome (Kla)"], [31, 11, 9, 21])
    validate_counts(proliferation, [9, 3, 7, 12])
    assert len({row["GroupID"] for row in proliferation}) == 31
    pathway_samples = [row for row in pathways if row["ObservationType"] == "sample"]
    assert len(pathway_samples) == 504
    assert len({tuple(row[key] for key in ("PXD", "SampleGroup", "SampleID", "Pathway")) for row in pathway_samples}) == 504

    panels, details = [], []
    maximum_error = 0.0

    def record(panel, cell, series, category, index, row, value, unit, numerator=None, denominator=None):
        nonlocal maximum_error
        assert math.isfinite(value)
        if numerator is not None:
            error = abs(value - 100 * numerator / denominator)
            maximum_error = max(maximum_error, error)
            assert error < 1e-10
        details.append({
            "Sheet": panel, "Cell": cell, "Series": series, "Category": category,
            "Replicate": index, "PXD": row.get("PXD", ""), "SampleGroup": row.get("SampleGroup", ""),
            "SampleID": row.get("SampleID", ""), "GroupID": row.get("GroupID", ""),
            "ReferenceKey": row.get("ReferenceKey", ""), "ReferencePXD": row.get("ReferencePXD", ""),
            "ReferenceLabel": row.get("ReferenceLabelEn", ""), "RNA_sample_N": int(row["N"]) if row.get("N") else None,
            "Numerator": numerator, "Denominator": denominator, "Value": value, "Unit": unit,
            "ObservationType": row.get("ObservationType", "GroupID"), "SourceMode": row.get("SourceMode", ""),
            "SourceFile": row.get("SourceFile", ""),
        })

    def grouped(name, rows, series_specs, source_key):
        maximum_n = max(sum(row["Category"] == category and ("Dataset" not in row or row["Dataset"] == dataset)
                            for row in rows) for category in CATEGORIES for _, dataset, _, _, _ in series_specs)
        matrix = [["Series"] + [label for label in LABELS for _ in range(maximum_n)]]
        counts = {}
        for row_index, (series, dataset, value_key, numerator_key, denominator_key) in enumerate(series_specs, 2):
            output = [series]
            counts[series] = []
            for category_index, category in enumerate(CATEGORIES):
                subset = ordered([row for row in rows if row["Category"] == category and row.get("Dataset") == dataset])
                counts[series].append(len(subset))
                for replicate, row in enumerate(subset, 1):
                    value = float(row[value_key]) * (100 if source_key == "pathways" else 1)
                    output.append(value)
                    cell = f"{column_name(1 + category_index * maximum_n + replicate)}{row_index}"
                    record(name, cell, series, category, replicate, row, value, "%",
                           int(row[numerator_key]), int(row[denominator_key]))
                output.extend([None] * (maximum_n - len(subset)))
            matrix.append(output)
        panels.append({"name": name, "type": "Grouped", "replicates": maximum_n, "matrix": matrix,
                       "counts": counts, "unit": "%", "source": source_key})

    grouped("a_DDR比例", fraction, [
        ("Whole proteome", "Whole proteome", "DdrFractionPercentage", "DdrProteinCount", "ProteinCount"),
        ("Lactylome (Kla)", "Lactylome (Kla)", "DdrFractionPercentage", "DdrProteinCount", "ProteinCount"),
    ], "fraction")

    b_matrix = [LABELS] + [[None] * 4 for _ in range(12)]
    for category_index, category in enumerate(CATEGORIES):
        for replicate, row in enumerate(ordered([row for row in proliferation if row["Category"] == category]), 1):
            value = float(row["ProliferationMedian"])
            b_matrix[replicate][category_index] = value
            record("b_增殖得分", f"{column_name(category_index + 1)}{replicate + 1}", "Proliferation score",
                   category, replicate, row, value, "log2 score")
    panels.append({"name": "b_增殖得分", "type": "Column", "replicates": 1, "matrix": b_matrix,
                   "counts": {"Proliferation score": [9, 3, 7, 12]}, "unit": "log2 score", "source": "proliferation"})

    for letter, pathway in PATHWAYS:
        subset = [row for row in pathway_samples if row["Pathway"] == pathway]
        validate_counts(subset, [31, 11, 9, 21])
        grouped(f"{letter}_{pathway}", subset, [
            ("Pro", "Lactylome (Kla)", "PositiveFraction", "PositiveProteinCount", "KlaDdrProteinCount"),
            ("Inh", "Lactylome (Kla)", "NegativeFraction", "NegativeProteinCount", "KlaDdrProteinCount"),
        ], "pathways")

    assert len(details) == 1311
    assert len({(row["Sheet"], row["Cell"]) for row in details}) == 1311
    args.output_dir.mkdir(parents=True, exist_ok=True)
    snapshots = args.output_dir / "source_snapshot"
    snapshots.mkdir(exist_ok=True)
    for source in sources.values():
        target = snapshots / source.name
        if source.resolve() != target.resolve():
            shutil.copy2(source, target)
    with (args.output_dir / "逐点数据与样本编号.csv").open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(details[0]))
        writer.writeheader()
        writer.writerows(details)
    payload = {"panels": panels, "details": details,
               "sources": [{"key": key, "path": str(path.resolve()), "file": path.name,
                            "sha256": hashlib.sha256(path.read_bytes()).hexdigest()} for key, path in sources.items()],
               "validation": {"point_count": len(details), "maximum_percentage_error": maximum_error}}
    (args.output_dir / "prism_points.json").write_text(json.dumps(payload, ensure_ascii=False, allow_nan=False), encoding="utf-8")
    print(json.dumps({"panels": {panel["name"]: panel["counts"] for panel in panels}, "validation": payload["validation"]}, ensure_ascii=False))


if __name__ == "__main__":
    main()
