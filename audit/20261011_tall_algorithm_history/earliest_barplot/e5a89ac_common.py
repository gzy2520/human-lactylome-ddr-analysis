from __future__ import annotations

from collections import Counter
from pathlib import Path
import re
from typing import Iterable

import pandas as pd


SITE_COLUMNS = [
    "PXD",
    "DOI",
    "SampleName",
    "CellOrTissueType",
    "ExperimentalGroup",
    "Replicate",
    "EnrichmentStatus",
    "AcquisitionMode",
    "Accession",
    "BaseAccession",
    "GeneSymbol",
    "ProteinName",
    "KlaSite",
    "ModifiedPeptide",
    "LocalizationProb",
    "PEP",
    "Score",
    "DiagnosticIon",
    "DiagnosticIonIntensity",
    "ClassI",
    "SourceFile",
    "SourceRow",
    "SiteID",
    "EvidenceMode",
    "SourceConfidence",
    "PrimaryIncluded",
    "InclusionStatus",
    "ExclusionReason",
]


def clean_text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    text = str(value).strip()
    return "" if text.casefold() in {"nan", "none"} else text


def number(value: object) -> float | None:
    parsed = pd.to_numeric(pd.Series([clean_text(value)]), errors="coerce").iloc[0]
    return None if pd.isna(parsed) else float(parsed)


def integer(value: object) -> int | None:
    parsed = number(value)
    return None if parsed is None else int(parsed)


def is_true(value: object) -> bool:
    return clean_text(value).casefold() in {"+", "1", "true", "yes", "y"}


def split_tokens(value: object, pattern: str = r";") -> list[str]:
    return [token.strip() for token in re.split(pattern, clean_text(value)) if token.strip()]


def normalize_accession(value: object) -> str:
    text = clean_text(value)
    text = re.sub(r"^(?:REV__|CON__)+", "", text)
    text = re.sub(r"^(?:sp|tr)\|", "", text)
    if "|" in text:
        text = text.split("|", 1)[0]
    return text.strip()


def base_accession(value: object) -> str:
    return re.sub(r"-\d+$", "", normalize_accession(value))


def accession_candidates(value: object) -> list[str]:
    raw = clean_text(value)
    candidates: list[str] = []
    for token in re.split(r"[;:]", raw):
        accession = normalize_accession(token)
        if accession and accession not in candidates:
            candidates.append(accession)
    if not candidates:
        accession = normalize_accession(raw)
        if accession:
            candidates.append(accession)
    return candidates


def unique_join(values: Iterable[object]) -> str:
    tokens: set[str] = set()
    for value in values:
        tokens.update(split_tokens(value))
    return ";".join(sorted(tokens))


def best_annotation(values: Iterable[object]) -> str:
    cleaned = [clean_text(value) for value in values if clean_text(value)]
    if not cleaned:
        return ""
    counts = Counter(cleaned)
    return sorted(counts, key=lambda value: (-counts[value], -len(value), value))[0]


def annotation_from_description(value: object) -> tuple[str, str]:
    description = clean_text(value)
    gene_match = re.search(r"(?:^|\s)GN=([^\s]+)", description)
    gene = gene_match.group(1) if gene_match else ""
    protein_name = re.sub(r"\s+OS=.*$", "", description).strip()
    return gene, protein_name


def annotation_from_fasta(value: object, accession: str) -> tuple[str, str]:
    for header in split_tokens(value):
        token = header.split(maxsplit=1)[0]
        if base_accession(token) != base_accession(accession):
            continue
        return annotation_from_description(header.split(maxsplit=1)[1] if " " in header else "")
    return "", ""


def blank_site(**updates: object) -> dict[str, object]:
    row = {column: "" for column in SITE_COLUMNS}
    row.update(
        {
            "LocalizationProb": None,
            "PEP": None,
            "Score": None,
            "DiagnosticIonIntensity": None,
            "ClassI": False,
            "PrimaryIncluded": True,
            "InclusionStatus": "included",
            "ExclusionReason": "",
        }
    )
    row.update(updates)
    return row


def strip_peptide_modifications(value: object) -> str:
    peptide = clean_text(value)
    match = re.match(r"^[A-Z]\.([A-Z].*)\.[A-Z]$", peptide)
    if match:
        peptide = match.group(1)
    return re.sub(r"\([^)]*\)|\[[^]]*\]", "", peptide)


def lactyl_positions_from_ascore(value: object) -> list[int]:
    return sorted(
        {
            int(position)
            for position in re.findall(
                r"K(\d+)\s*:\s*(?:Lac|Lactyl(?:ation)?)", clean_text(value), re.I
            )
        }
    )


def lactyl_positions_from_peptide(value: object) -> list[int]:
    peptide = clean_text(value)
    positions: list[int] = []
    residue_position = 0
    index = 0
    while index < len(peptide):
        char = peptide[index]
        if "A" <= char <= "Z":
            residue_position += 1
            if char == "K" and index + 1 < len(peptide) and peptide[index + 1] in "([":
                closer = ")" if peptide[index + 1] == "(" else "]"
                end = peptide.find(closer, index + 2)
                if end != -1:
                    modification = peptide[index + 2 : end]
                    if re.search(r"72\.02|\bLac\b|Lactyl", modification, re.I):
                        positions.append(residue_position)
                    index = end
        index += 1
    return sorted(set(positions))


def parse_probability_values(value: object) -> list[float]:
    probabilities: list[float] = []
    for match in re.findall(r"\(([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)\)", clean_text(value)):
        try:
            probabilities.append(float(match))
        except ValueError:
            continue
    return probabilities


def apply_annotation_supplement(sites: pd.DataFrame, path: Path) -> pd.DataFrame:
    if not path.exists() or sites.empty:
        return sites
    supplement = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    gene_map = supplement.set_index("Entry")["Gene Names (primary)"].to_dict()
    name_map = supplement.set_index("Entry")["Protein names"].to_dict()
    sites = sites.copy()
    missing_gene = sites["GeneSymbol"].fillna("").str.strip().eq("")
    missing_name = sites["ProteinName"].fillna("").str.strip().eq("")
    sites.loc[missing_gene, "GeneSymbol"] = sites.loc[missing_gene, "BaseAccession"].map(gene_map).fillna("")
    sites.loc[missing_name, "ProteinName"] = sites.loc[missing_name, "BaseAccession"].map(name_map).fillna("")
    return sites


def read_go_annotations(path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    raw = pd.read_csv(path, sep="\t", dtype=str, low_memory=False, encoding="utf-8-sig").fillna("")
    raw["BaseAccession"] = raw["GENE PRODUCT ID"].map(base_accession)
    taxon = pd.to_numeric(raw["TAXON ID"], errors="coerce")
    raw["ExcludedNOT"] = raw["QUALIFIER"].str.contains(r"(?:^|\|)NOT(?:\||$)", case=False, regex=True)
    retained = raw[(taxon.isna() | taxon.eq(9606)) & ~raw["ExcludedNOT"] & raw["BaseAccession"].ne("")].copy()
    rows = []
    for accession, group in retained.groupby("BaseAccession", sort=True):
        rows.append(
            {
                "BaseAccession": accession,
                "GOSymbol": best_annotation(group["SYMBOL"]),
                "GOTerms": unique_join(group["GO TERM"]),
                "GONames": unique_join(group["GO NAME"]),
                "GOEvidenceCodes": unique_join(group["GO EVIDENCE CODE"]),
                "GOReferences": unique_join(group["REFERENCE"]),
                "GOAnnotationCount": len(group),
            }
        )
    return pd.DataFrame(rows), raw


def relative_path(path: Path, root: Path) -> str:
    try:
        return str(path.resolve().relative_to(root.resolve()))
    except ValueError:
        return str(path.resolve())
