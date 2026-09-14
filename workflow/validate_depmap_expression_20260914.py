#!/usr/bin/env python3
"""Validate the downloaded DepMap candidate structure without using gene symbols as keys."""
import csv
import re
import sys
from collections import Counter
from pathlib import Path

TARGETS = {
    "ACH-000019", "ACH-000971", "ACH-000739", "ACH-000681",
    "ACH-000849", "ACH-000147", "ACH-000943",
}
ENTREZ_HEADER = re.compile(r"^.+ \((\d+)\)$")


def main(root: Path, output: Path) -> None:
    base = root / "raw" / "depmap" / "DepMap_Public_24Q4"
    expression = base / "OmicsExpressionProteinCodingGenesTPMLogp1.csv"
    model = base / "Model.csv"
    profiles = base / "OmicsProfiles.csv"
    for path in (expression, model, profiles):
        if not path.exists():
            raise FileNotFoundError(path)

    with expression.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.reader(handle)
        header = next(reader)
        expression_rows = 0
        target_rows = set()
        for row in reader:
            expression_rows += 1
            if row and row[0] in TARGETS:
                target_rows.add(row[0])
    entrez_headers = [match.group(1) for value in header[1:] if (match := ENTREZ_HEADER.match(value))]
    if len(entrez_headers) != len(header) - 1:
        raise ValueError("expression headers do not provide numeric Entrez keys")

    with model.open(newline="", encoding="utf-8-sig") as handle:
        model_rows = list(csv.DictReader(handle))
    model_targets = {row.get("ModelID") for row in model_rows} & TARGETS

    with profiles.open(newline="", encoding="utf-8-sig") as handle:
        profile_rows = list(csv.DictReader(handle))
    target_profiles = [row for row in profile_rows if row.get("ModelID") in TARGETS]
    rna_profiles = [row for row in target_profiles if row.get("Datatype", "").lower() == "rna"]
    if target_rows != TARGETS or model_targets != TARGETS or len(rna_profiles) != len(TARGETS):
        raise ValueError("not all target model IDs have expression, Model and RNA profile rows")

    duplicate_entrez = len(entrez_headers) - len(set(entrez_headers))
    records = [
        ("ExpressionRows", str(expression_rows), "streamed rows in the downloaded model-keyed matrix"),
        ("ExpressionColumns", str(len(header) - 1), "gene columns excluding the leading ModelID column"),
        ("EntrezGeneHeaders", str(len(entrez_headers)), "numeric Entrez IDs parsed from headers; symbols are not analytical keys"),
        ("UniqueEntrezGeneIDs", str(len(set(entrez_headers))), "unique numeric Entrez IDs before downstream duplicate resolution"),
        ("DuplicateEntrezGeneHeaders", str(duplicate_entrez), "duplicate stable IDs require deterministic downstream resolution"),
        ("TargetModelRows", str(len(target_rows)), ";".join(sorted(target_rows))),
        ("ModelMetadataRows", str(len(model_targets)), ";".join(sorted(model_targets))),
        ("TargetRNAProfiles", str(len(rna_profiles)), "one Datatype=rna profile for each target model"),
        ("Validation", "PASS_WITH_MAPPING_CAVEAT", "all seven target ACH IDs resolved; duplicate Entrez columns remain for downstream QC"),
    ]
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("Metric", "Value", "Detail"))
        writer.writerows(records)
    print("DEPMAP_STRUCTURE_VALIDATION_PASS", len(TARGETS), "target models", len(entrez_headers), "Entrez gene columns")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit("Usage: validate_depmap_expression_20260914.py <server_root> <output_tsv>")
    main(Path(sys.argv[1]), Path(sys.argv[2]))
