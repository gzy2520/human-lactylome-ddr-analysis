from __future__ import annotations

from pathlib import Path
import re

import pandas as pd

from common import clean_text, relative_path, unique_join


TARGETS = {
    "MRE11": ("P49959", {"MRE11"}),
    "XLF/NHEJ1": ("Q9H9Q4", {"NHEJ1", "XLF"}),
    "NBS1/NBN": ("O60934", {"NBN", "NBS1"}),
}


def target_pattern(accession: str, genes: set[str]) -> re.Pattern[str]:
    aliases = "|".join(re.escape(value) for value in sorted({accession, *genes}, key=len, reverse=True))
    return re.compile(rf"(?<![A-Za-z0-9])(?:{aliases})(?![A-Za-z0-9])", re.I)


def target_mask(table: pd.DataFrame, pattern: re.Pattern[str]) -> pd.Series:
    mask = pd.Series(False, index=table.index)
    likely = [
        column
        for column in table.columns
        if re.search(r"accession|protein|gene|description|fasta", str(column), re.I)
    ]
    for column in likely or list(table.columns):
        mask |= table[column].fillna("").astype(str).str.contains(pattern, na=False)
    return mask


def lactyl_mask(table: pd.DataFrame) -> pd.Series:
    mask = pd.Series(False, index=table.index)
    for column in table.columns:
        if re.search(r"peptide|ptm|ascore|modification|modified", str(column), re.I):
            mask |= table[column].fillna("").astype(str).str.contains(
                r"Lactyl|\bLac\b|K\(\+?72\.02|K\(Lac|La\s*\(K\)",
                case=False,
                regex=True,
                na=False,
            )
        if re.search(r"La\s*\(K\)", str(column), re.I):
            numeric = pd.to_numeric(table[column], errors="coerce")
            mask |= numeric.fillna(0).gt(0)
            if re.search(r"site\s*IDs?", str(column), re.I):
                mask |= table[column].fillna("").astype(str).str.strip().ne("")
    if "PTM.ModificationTitle" in table:
        mask |= table["PTM.ModificationTitle"].fillna("").str.casefold().eq("lac")
    return mask


def read_table(path: Path) -> pd.DataFrame:
    separator = "\t" if path.suffix.casefold() in {".txt", ".tsv"} else ","
    return pd.read_csv(path, sep=separator, dtype=str, low_memory=False).fillna("")


def scan_text_file(
    path: Path,
    project_root: Path,
    rows: list[dict[str, object]],
    kla_defining: bool,
) -> None:
    table = read_table(path)
    lactyl = lactyl_mask(table) if kla_defining else pd.Series(False, index=table.index)
    for target, (accession, genes) in TARGETS.items():
        matches = target_mask(table, target_pattern(accession, genes))
        if not matches.any():
            continue
        rows.append(
            {
                "Target": target,
                "SourceFile": relative_path(path, project_root),
                "TargetRows": int(matches.sum()),
                "KlaTargetRows": int((matches & lactyl).sum()) if kla_defining else 0,
                "SourceRole": "kla_candidate_table" if kla_defining else "protein_or_peptide_support_table",
            }
        )


def scan_excel_file(path: Path, project_root: Path, rows: list[dict[str, object]]) -> None:
    for sheet in pd.ExcelFile(path).sheet_names:
        table = pd.read_excel(path, sheet_name=sheet, header=None, dtype=str).fillna("")
        if table.empty:
            continue
        for target, (accession, genes) in TARGETS.items():
            pattern = target_pattern(accession, genes)
            matches = pd.Series(False, index=table.index)
            for column in table.columns:
                matches |= table[column].astype(str).str.contains(pattern, na=False)
            if matches.any():
                rows.append(
                    {
                        "Target": target,
                        "SourceFile": f"{relative_path(path, project_root)}#{sheet}",
                        "TargetRows": int(matches.sum()),
                        "KlaTargetRows": int(matches.sum()),
                        "SourceRole": "author_kla_supplementary_table",
                    }
                )


def source_files(project_root: Path, pxd: str) -> list[tuple[Path, bool]]:
    root = project_root / "data" / pxd
    files: list[tuple[Path, bool]] = []
    if pxd == "PXD014870":
        files.extend((path, "Sites" in path.name or path.name == "modificationSpecificPeptides.txt") for path in root.rglob("Lactyl (K)Sites.txt"))
        files.extend((path, True) for path in root.rglob("modificationSpecificPeptides.txt"))
        files.extend((path, False) for path in root.rglob("proteinGroups.txt"))
    elif pxd == "PXD028488":
        files.extend((path, True) for path in root.rglob("DB search psm.csv"))
        files.extend((path, False) for path in root.rglob("protein-peptides.csv"))
        files.extend((path, False) for path in root.rglob("proteins.csv"))
    elif pxd == "PXD053474":
        files.extend((path, True) for path in root.rglob("peptide.csv"))
        files.extend((path, True) for path in root.rglob("*ptm-site*Report.tsv"))
        files.extend((path, False) for path in root.rglob("protein-peptides.csv"))
        files.extend((path, False) for path in root.rglob("proteins.csv"))
    return sorted(set(files), key=lambda item: str(item[0]))


def build_target_source_audit(project_root: Path, extracted_evidence: pd.DataFrame) -> pd.DataFrame:
    scan_rows: list[dict[str, object]] = []
    scanned_pxds = ["PXD014870", "PXD028488", "PXD050470", "PXD053474"]
    excel_files = {
        "PXD050470": sorted((project_root / "data/PXD050470/supplementary").glob("*.xlsx")),
        "PXD053474": [project_root / "data/PXD053474/supplementary/js4c00366_si_003.xlsx"],
    }
    for pxd in scanned_pxds:
        local_rows: list[dict[str, object]] = []
        for path, kla_defining in source_files(project_root, pxd):
            scan_text_file(path, project_root, local_rows, kla_defining)
        for path in excel_files.get(pxd, []):
            if path.exists():
                scan_excel_file(path, project_root, local_rows)
        for row in local_rows:
            row["PXD"] = pxd
        scan_rows.extend(local_rows)

    scan = pd.DataFrame(scan_rows)
    output_rows = []
    for pxd in ["PXD014870", "PXD028488", "PXD050470", "PXD053474", "PXD060185", "PXD078013", "PXD078736"]:
        for target, (accession, genes) in TARGETS.items():
            extracted = extracted_evidence[
                extracted_evidence["PXD"].eq(pxd)
                & (
                    extracted_evidence["BaseAccession"].eq(accession)
                    | extracted_evidence["GeneSymbol"].fillna("").map(
                        lambda value: bool({token.strip().upper() for token in str(value).split(";")} & genes)
                    )
                )
            ]
            source = scan[(scan["PXD"].eq(pxd)) & (scan["Target"].eq(target))] if not scan.empty else scan
            target_rows = int(source["TargetRows"].sum()) if not source.empty else 0
            kla_rows = int(source["KlaTargetRows"].sum()) if not source.empty else 0
            if not extracted.empty:
                conclusion = "present_in_extracted_primary_kla_evidence"
            elif kla_rows:
                conclusion = "kla_candidate_rows_present_but_not_selected_review_required"
            elif target_rows:
                conclusion = "protein_or_unmodified_peptide_present_but_no_kla_evidence"
            elif pxd in scanned_pxds:
                conclusion = "not_found_in_available_author_search_or_kla_tables"
            else:
                conclusion = "not_present_in_extracted_kla_evidence_source_not_deep_scanned"
            output_rows.append(
                {
                    "PXD": pxd,
                    "Target": target,
                    "Accession": accession,
                    "GeneAliases": ";".join(sorted(genes)),
                    "SourceTargetRows": target_rows,
                    "SourceKlaTargetRows": kla_rows,
                    "ExtractedPrimaryRows": len(extracted),
                    "Conclusion": conclusion,
                    "SourceFiles": unique_join(source["SourceFile"]) if not source.empty else "",
                    "Interpretation": (
                        "The target was detected at protein/unmodified-peptide level, but the available author tables do not assign Kla to it."
                        if conclusion == "protein_or_unmodified_peptide_present_but_no_kla_evidence"
                        else ""
                    ),
                }
            )
    return pd.DataFrame(output_rows)
