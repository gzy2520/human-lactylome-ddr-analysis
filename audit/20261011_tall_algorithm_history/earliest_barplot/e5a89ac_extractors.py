from __future__ import annotations

import csv
from pathlib import Path
import re
from typing import Iterable

import pandas as pd

from common import (
    SITE_COLUMNS,
    accession_candidates,
    annotation_from_description,
    annotation_from_fasta,
    base_accession,
    blank_site,
    clean_text,
    integer,
    is_true,
    lactyl_positions_from_ascore,
    lactyl_positions_from_peptide,
    normalize_accession,
    number,
    parse_probability_values,
    relative_path,
    split_tokens,
    strip_peptide_modifications,
)


def exclusion_row(
    pxd: str,
    source: str,
    reason: str,
    count: int = 1,
    detail: str = "",
) -> dict[str, object]:
    return {"PXD": pxd, "Source": source, "Reason": reason, "Count": count, "Detail": detail}


def dataframe(rows: list[dict[str, object]]) -> pd.DataFrame:
    if not rows:
        return pd.DataFrame(columns=SITE_COLUMNS)
    frame = pd.DataFrame(rows)
    for column in SITE_COLUMNS:
        if column not in frame.columns:
            frame[column] = ""
    return frame


def valid_maxquant_site_ids(path: Path) -> tuple[set[str], list[dict[str, object]]]:
    table = pd.read_csv(path, sep="\t", dtype=str, low_memory=False).fillna("")
    valid: set[str] = set()
    logs: list[dict[str, object]] = []
    for _, row in table.iterrows():
        if any(is_true(row.get(column, "")) for column in ("Reverse", "Potential contaminant", "Contaminant")):
            continue
        if (number(row.get("Lactyl (K)", "")) or 0) <= 0:
            continue
        valid.update(split_tokens(row.get("Lactyl (K) site IDs", "")))
    return valid, logs


def extract_pxd014870(
    data_root: Path,
    project_root: Path,
    doi: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    root = data_root / "PXD014870/search_results"
    source_map = {
        "MCF7_DCA_SILAC_Kla_IP": ("MCF7_DCA", "DCA"),
        "MCF7_Hypoxia_SILAC_Kla_IP": ("MCF7_Hypoxia", "Hypoxia"),
        "MCF7_Rotenone_SILAC_Kla_IP": ("MCF7_Rotenone", "Rotenone"),
        "MCF7_U_13C6_Glucose_Kla_IP": ("MCF7_U13C6_Glucose", "U-13C6 glucose"),
    }
    rows: list[dict[str, object]] = []
    logs: list[dict[str, object]] = []
    for directory_name, (source_token, group) in source_map.items():
        txt = root / directory_name / "txt"
        site_path = txt / "Lactyl (K)Sites.txt"
        peptide_path = txt / "modificationSpecificPeptides.txt"
        if not site_path.exists() or not peptide_path.exists():
            logs.append(exclusion_row("PXD014870", directory_name, "missing_required_search_table"))
            continue
        valid_ids, _ = valid_maxquant_site_ids(peptide_path)
        sites = pd.read_csv(site_path, sep="\t", dtype=str, low_memory=False).fillna("")
        reason_counts: dict[str, int] = {}
        for index, row in sites.iterrows():
            reason = ""
            if any(is_true(row.get(column, "")) for column in ("Reverse", "Potential contaminant", "Contaminant")):
                reason = "reverse_or_contaminant"
            site_id = clean_text(row.get("id", ""))
            if not reason and (not site_id or site_id not in valid_ids):
                reason = "no_valid_kla_site_id_in_modificationSpecificPeptides"
            if not reason and clean_text(row.get("Amino acid", "K")).upper() != "K":
                reason = "non_lysine_site"
            accession = normalize_accession(row.get("Protein", row.get("Leading proteins", "")))
            position = integer(row.get("Position", row.get("Positions within proteins", "")))
            if not reason and (not accession or position is None):
                reason = "missing_accession_or_kla_position"
            if reason:
                reason_counts[reason] = reason_counts.get(reason, 0) + 1
                continue
            gene, protein_name = annotation_from_fasta(row.get("Fasta headers", ""), accession)
            if not gene:
                gene = split_tokens(row.get("Gene names", ""))[0] if split_tokens(row.get("Gene names", "")) else ""
            if not protein_name:
                protein_name = split_tokens(row.get("Protein names", ""))[0] if split_tokens(row.get("Protein names", "")) else ""
            localization = number(row.get("Localization prob", ""))
            raw_file = clean_text(row.get("Best localization raw file", "")) or source_token
            diagnostic = clean_text(row.get("Diagnostic peak", ""))
            rows.append(
                blank_site(
                    PXD="PXD014870",
                    DOI=doi,
                    SampleName=raw_file,
                    CellOrTissueType="MCF7",
                    ExperimentalGroup=group,
                    Replicate=raw_file,
                    EnrichmentStatus="Enriched",
                    AcquisitionMode="DDA",
                    Accession=accession,
                    BaseAccession=base_accession(accession),
                    GeneSymbol=gene,
                    ProteinName=protein_name,
                    KlaSite=f"K{position}",
                    ModifiedPeptide=clean_text(row.get("Modified sequence", "")),
                    LocalizationProb=localization,
                    PEP=number(row.get("PEP", "")),
                    Score=number(row.get("Score", "")),
                    DiagnosticIon=diagnostic,
                    ClassI=localization is not None and localization >= 0.75,
                    SourceFile=relative_path(site_path, project_root),
                    SourceRow=index + 2,
                    SiteID=site_id,
                    EvidenceMode="maxquant_site_table_plus_modified_peptides",
                    SourceConfidence="high_localized" if localization is not None and localization >= 0.75 else "standard_author_search_result",
                )
            )
        for reason, count in sorted(reason_counts.items()):
            logs.append(exclusion_row("PXD014870", relative_path(site_path, project_root), reason, count))
    frame = dataframe(rows)
    if not frame.empty:
        frame["Sensitivity075Pass"] = pd.to_numeric(frame["LocalizationProb"], errors="coerce").ge(0.75)
    return frame, pd.DataFrame(logs)


def peaks_annotations(proteins: pd.DataFrame) -> dict[str, tuple[str, str]]:
    annotations: dict[str, tuple[str, str]] = {}
    for _, row in proteins.iterrows():
        accession = normalize_accession(row.get("Accession", ""))
        gene, name = annotation_from_description(row.get("Description", ""))
        if accession and accession not in annotations:
            annotations[accession] = (gene, name)
    return annotations


def peaks_start_map(protein_peptides: pd.DataFrame) -> dict[str, list[tuple[str, int]]]:
    mapping: dict[str, list[tuple[str, int]]] = {}
    for _, row in protein_peptides.iterrows():
        accession = normalize_accession(row.get("Protein Accession", ""))
        peptide = strip_peptide_modifications(row.get("Peptide", ""))
        start = integer(row.get("Start", ""))
        if accession and peptide and start is not None:
            mapping.setdefault(peptide, []).append((accession, start))
    return mapping


def matching_accession_mappings(
    matches: list[tuple[str, int]], raw_accession: object
) -> list[tuple[str, int]]:
    candidates = {base_accession(item) for item in accession_candidates(raw_accession)}
    if not candidates:
        return matches
    return [item for item in matches if base_accession(item[0]) in candidates]


def site_position_pairs(site_ids: list[str], positions: list[str]) -> list[tuple[str, str]]:
    if len(site_ids) != len(positions):
        return []
    return list(zip(site_ids, positions))


def normalize_scan(value: object) -> str:
    text = clean_text(value)
    return text.rsplit(":", 1)[-1]


def read_marker_156(path: Path) -> dict[tuple[str, str, str], float]:
    markers: dict[tuple[str, str, str], float] = {}
    if not path.exists():
        return markers
    with path.open(errors="replace", newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader, [])
        positions = {name: index for index, name in enumerate(header[:18])}
        required = {"Peptide", "scan", "Source File", "ion", "theo m/z", "ion relative intensity(%)"}
        if not required.issubset(positions):
            return markers
        for raw in reader:
            if len(raw) < 18:
                raw += [""] * (18 - len(raw))
            if clean_text(raw[positions["ion"]]).casefold() != "marker(156)":
                continue
            theoretical = number(raw[positions["theo m/z"]])
            if theoretical is None or abs(theoretical - 156.1025) > 0.01:
                continue
            key = (
                clean_text(raw[positions["Peptide"]]),
                Path(clean_text(raw[positions["Source File"]])).name,
                normalize_scan(raw[positions["scan"]]),
            )
            intensity = number(raw[positions["ion relative intensity(%)"]]) or 0.0
            markers[key] = max(markers.get(key, 0.0), intensity)
    return markers


def pxd028_directory_metadata(directory: Path) -> dict[str, object]:
    name = directory.name
    parent = directory.parent.name
    enriched = parent.startswith("Enrichment")
    species = "mouse" if name.startswith("BV2") or name.startswith("RAW-") else "human"
    if name.startswith("HEK293T"):
        cell_type = "HEK293T"
    elif name.startswith("HCT116"):
        cell_type = "HCT116"
    elif name.startswith("TALL"):
        cell_type = "T-ALL"
    elif name.startswith("BV2"):
        cell_type = "BV2"
    else:
        cell_type = "RAW264.7"
    aggregate = "all HCD" in name
    hcd = re.search(r"HCD(\d+)", name)
    hcd_label = f"HCD{hcd.group(1)}" if hcd else ("all HCD aggregate" if aggregate else "not specified")
    return {
        "Directory": name,
        "SearchCollection": parent,
        "Species": species,
        "CellType": cell_type,
        "EnrichmentStatus": "Enriched" if enriched else "Unenriched",
        "HCDCondition": hcd_label,
        "AggregateDirectory": aggregate,
    }


def extract_pxd028488(
    data_root: Path,
    project_root: Path,
    doi: str,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    root = data_root / "PXD028488/search_results"
    rows: list[dict[str, object]] = []
    logs: list[dict[str, object]] = []
    directory_audit: list[dict[str, object]] = []
    for collection in ("Enrichment-Search files", "Nonenrichment-Search files"):
        for directory in sorted((root / collection).glob("*-Search files")):
            meta = pxd028_directory_metadata(directory)
            psm_path = directory / "DB search psm.csv"
            meta["PSMAvailable"] = psm_path.exists()
            if meta["Species"] != "human":
                meta.update({"PrimaryStatus": "excluded", "Reason": "non_human_cell_line", "KlaEvidenceRows": 0})
                directory_audit.append(meta)
                logs.append(exclusion_row("PXD028488", relative_path(directory, project_root), "non_human_directory", detail=str(meta["CellType"])))
                continue
            if not psm_path.exists():
                meta.update({"PrimaryStatus": "excluded", "Reason": "missing_DB_search_psm.csv", "KlaEvidenceRows": 0})
                directory_audit.append(meta)
                logs.append(exclusion_row("PXD028488", relative_path(directory, project_root), "missing_required_psm_table"))
                continue
            validation_only = bool(meta["AggregateDirectory"])
            meta["PrimaryStatus"] = "validation_only" if validation_only else "included"
            meta["Reason"] = "aggregate_of_individual_HCD_runs" if validation_only else "human_analyzable_search_directory"

            psm = pd.read_csv(psm_path, dtype=str, low_memory=False).fillna("")
            protein_peptides = pd.read_csv(directory / "protein-peptides.csv", dtype=str, low_memory=False).fillna("")
            proteins = pd.read_csv(directory / "proteins.csv", dtype=str, low_memory=False).fillna("")
            annotations = peaks_annotations(proteins)
            starts = peaks_start_map(protein_peptides)
            markers = read_marker_156(directory / "PSM ions.csv")
            counters = {
                "score_below_20": 0,
                "no_kla_modification": 0,
                "no_position_mapping": 0,
                "accession_mapping_mismatch": 0,
                "reverse_or_contaminant": 0,
            }
            before = len(rows)

            for index, row in psm.iterrows():
                raw_accession = clean_text(row.get("Accession", ""))
                if re.search(r"DECOY|REV|REVERSED|CON_|CONTAM|CRAP", raw_accession, re.I):
                    counters["reverse_or_contaminant"] += 1
                    continue
                score = number(row.get("-10lgP", ""))
                if score is None or score < 20:
                    counters["score_below_20"] += 1
                    continue
                peptide_raw = clean_text(row.get("Peptide", ""))
                positions = lactyl_positions_from_ascore(row.get("AScore", ""))
                if not positions:
                    positions = lactyl_positions_from_peptide(peptide_raw)
                ptm_has_lac = bool(re.search(r"(?:^|;)\s*(?:Lac|Lactyl(?:ation)?)\b", clean_text(row.get("PTM", "")), re.I))
                if not positions and not ptm_has_lac:
                    counters["no_kla_modification"] += 1
                    continue
                if not positions:
                    counters["no_position_mapping"] += 1
                    continue
                peptide_clean = strip_peptide_modifications(peptide_raw)
                matches = starts.get(peptide_clean, [])
                mapped_matches = matching_accession_mappings(matches, raw_accession)
                if matches and not mapped_matches:
                    counters["accession_mapping_mismatch"] += 1
                matches = mapped_matches
                if not matches:
                    counters["no_position_mapping"] += 1
                    continue
                source_file = Path(clean_text(row.get("Source File", ""))).name
                marker_key = (peptide_raw, source_file, normalize_scan(row.get("Scan", "")))
                marker_intensity = markers.get(marker_key)
                for accession, start in matches:
                    gene, protein_name = annotations.get(accession, ("", ""))
                    for peptide_position in positions:
                        protein_position = start + peptide_position - 1
                        rows.append(
                            blank_site(
                                PXD="PXD028488",
                                DOI=doi,
                                SampleName=source_file or str(meta["Directory"]),
                                CellOrTissueType=meta["CellType"],
                                ExperimentalGroup=f"{meta['EnrichmentStatus']} {meta['HCDCondition']}",
                                Replicate=source_file,
                                EnrichmentStatus=meta["EnrichmentStatus"],
                                AcquisitionMode="DDA",
                                Accession=accession,
                                BaseAccession=base_accession(accession),
                                GeneSymbol=gene,
                                ProteinName=protein_name,
                                KlaSite=f"K{protein_position}",
                                ModifiedPeptide=peptide_raw,
                                Score=score,
                                DiagnosticIon="CycIm m/z 156.1025" if marker_intensity is not None else "not detected in exported ion assignment",
                                DiagnosticIonIntensity=marker_intensity,
                                SourceFile=relative_path(psm_path, project_root),
                                SourceRow=index + 2,
                                SiteID=f"{meta['Directory']}:{index + 2}:{accession}:K{protein_position}",
                                EvidenceMode="PEAKS_PSM_plus_protein_peptides",
                                SourceConfidence="high_diagnostic_supported" if marker_intensity is not None else "standard_search_result",
                                PrimaryIncluded=not validation_only,
                                InclusionStatus="validation_only" if validation_only else "included",
                                ExclusionReason="aggregate_HCD_duplicate" if validation_only else "",
                            )
                        )
            evidence_count = len(rows) - before
            meta["KlaEvidenceRows"] = evidence_count
            meta["DiagnosticSupportedRows"] = sum(
                1 for item in rows[before:] if item["DiagnosticIonIntensity"] not in (None, "")
            )
            directory_audit.append(meta)
            for reason, count in counters.items():
                if count:
                    logs.append(exclusion_row("PXD028488", relative_path(psm_path, project_root), reason, count))
    return dataframe(rows), pd.DataFrame(logs), pd.DataFrame(directory_audit)


def extract_pxd050470(
    data_root: Path,
    project_root: Path,
    doi: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    root = data_root / "PXD050470/supplementary"
    sources = (
        ("prca2331-sup-0005-tables3.xlsx", 12, "Proteins accession", "Positions within proteins", "Gene names", "Table S3"),
        ("prca2331-sup-0014-tables12.xlsx", 10, "Proteins accession", "Positions within proteins", "Gene name", "Table S12"),
    )
    rows: list[dict[str, object]] = []
    logs: list[dict[str, object]] = []
    for filename, header, accession_col, position_col, gene_col, table_name in sources:
        path = root / filename
        if not path.exists():
            logs.append(exclusion_row("PXD050470", filename, "missing_supplementary_table"))
            continue
        table = pd.read_excel(path, sheet_name="Sheet1", header=header, dtype=str).fillna("")
        table.columns = [clean_text(column) for column in table.columns]
        intensity_columns = [column for column in table.columns if clean_text(column).startswith("Intensity_")]
        for index, row in table.iterrows():
            accession = normalize_accession(row.get(accession_col, ""))
            position = integer(row.get(position_col, ""))
            if not accession or position is None:
                continue
            sample_values = [column for column in intensity_columns if clean_text(row.get(column, ""))]
            if not sample_values:
                sample_values = ["author_table_aggregate"]
            for sample_column in sample_values:
                sample = sample_column.removeprefix("Intensity_") if sample_column != "author_table_aggregate" else sample_column
                localization = number(row.get("Localization probability", ""))
                rows.append(
                    blank_site(
                        PXD="PXD050470",
                        DOI=doi,
                        SampleName=sample,
                        CellOrTissueType="Human hippocampus",
                        ExperimentalGroup="Physiological hippocampus",
                        Replicate=sample,
                        EnrichmentStatus="Enriched",
                        AcquisitionMode="DDA",
                        Accession=accession,
                        BaseAccession=base_accession(accession),
                        GeneSymbol=clean_text(row.get(gene_col, "")),
                        ProteinName="",
                        KlaSite=f"K{position}",
                        ModifiedPeptide=clean_text(row.get("Modified sequence", "")),
                        LocalizationProb=localization,
                        PEP=number(row.get("PEP", "")),
                        Score=number(row.get("Score", "")),
                        ClassI=True,
                        SourceFile=relative_path(path, project_root),
                        SourceRow=index + header + 2,
                        SiteID=f"{table_name}:{index + header + 2}:{accession}:K{position}",
                        EvidenceMode="author_supplementary_table",
                        SourceConfidence="high_author_class_I",
                    )
                )
    if not any((data_root / "PXD050470/raw").iterdir()):
        logs.append(
            exclusion_row(
                "PXD050470",
                relative_path(data_root / "PXD050470/raw", project_root),
                "repository_raw_files_not_present_locally",
                27,
                "Raw files were not re-searched; S3 and S12 are the analysis evidence.",
            )
        )
    return dataframe(rows), pd.DataFrame(logs)


def extract_pxd053474_dda(
    data_root: Path,
    project_root: Path,
    doi: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    root = data_root / "PXD053474/search_results/extracted"
    rows: list[dict[str, object]] = []
    logs: list[dict[str, object]] = []
    for peptide_path in sorted(root.rglob("peptide.csv")):
        directory = peptide_path.parent
        pp_path = directory / "protein-peptides.csv"
        proteins_path = directory / "proteins.csv"
        if not pp_path.exists() or not proteins_path.exists():
            logs.append(exclusion_row("PXD053474", relative_path(directory, project_root), "missing_DDA_mapping_table"))
            continue
        peptide_table = pd.read_csv(peptide_path, dtype=str, low_memory=False).fillna("")
        protein_peptides = pd.read_csv(pp_path, dtype=str, low_memory=False).fillna("")
        proteins = pd.read_csv(proteins_path, dtype=str, low_memory=False).fillna("")
        starts = peaks_start_map(protein_peptides)
        annotations = peaks_annotations(proteins)
        enriched = "Enriched-DDA" in str(peptide_path)
        compartment = directory.name.replace("-Search files", "")
        counters = {
            "score_below_20": 0,
            "no_kla_modification": 0,
            "no_position_mapping": 0,
            "accession_mapping_mismatch": 0,
        }
        for index, row in peptide_table.iterrows():
            score = number(row.get("-10lgP", ""))
            if score is None or score < 20:
                counters["score_below_20"] += 1
                continue
            peptide_raw = clean_text(row.get("Peptide", ""))
            positions = lactyl_positions_from_ascore(row.get("AScore", "")) or lactyl_positions_from_peptide(peptide_raw)
            ptm_has_lac = bool(re.search(r"(?:^|;)\s*(?:Lac|Lactyl(?:ation)?)\b", clean_text(row.get("PTM", "")), re.I))
            if not positions and not ptm_has_lac:
                counters["no_kla_modification"] += 1
                continue
            if not positions:
                counters["no_position_mapping"] += 1
                continue
            peptide_clean = strip_peptide_modifications(peptide_raw)
            matches = starts.get(peptide_clean, [])
            mapped_matches = matching_accession_mappings(matches, row.get("Accession", ""))
            if matches and not mapped_matches:
                counters["accession_mapping_mismatch"] += 1
            matches = mapped_matches
            if not matches:
                counters["no_position_mapping"] += 1
                continue
            sample = Path(clean_text(row.get("Source File", ""))).name or compartment
            for accession, start in matches:
                gene, protein_name = annotations.get(accession, ("", ""))
                for peptide_position in positions:
                    position = start + peptide_position - 1
                    rows.append(
                        blank_site(
                            PXD="PXD053474",
                            DOI=doi,
                            SampleName=sample,
                            CellOrTissueType="HCT116",
                            ExperimentalGroup=compartment,
                            Replicate=sample,
                            EnrichmentStatus="Enriched" if enriched else "Unenriched",
                            AcquisitionMode="DDA",
                            Accession=accession,
                            BaseAccession=base_accession(accession),
                            GeneSymbol=gene,
                            ProteinName=protein_name,
                            KlaSite=f"K{position}",
                            ModifiedPeptide=peptide_raw,
                            Score=score,
                            DiagnosticIon="not available in deposited PXD053474 DDA tables",
                            SourceFile=relative_path(peptide_path, project_root),
                            SourceRow=index + 2,
                            SiteID=f"DDA:{directory.name}:{index + 2}:{accession}:K{position}",
                            EvidenceMode="PEAKS_DDA_peptide_plus_protein_peptides",
                            SourceConfidence="standard_search_result",
                        )
                    )
        for reason, count in counters.items():
            if count:
                logs.append(exclusion_row("PXD053474", relative_path(peptide_path, project_root), reason, count))
    return dataframe(rows), pd.DataFrame(logs)


def run_prefixes(columns: Iterable[str]) -> list[tuple[str, str]]:
    prefixes: list[tuple[str, str]] = []
    suffix = ".PG.IsIdentified"
    for column in columns:
        if not column.endswith(suffix):
            continue
        prefix = column[: -len(suffix)]
        sample = re.sub(r"^\[\d+\]\s*", "", prefix)
        prefixes.append((prefix, sample))
    return prefixes


def extract_pxd053474_dia(
    data_root: Path,
    project_root: Path,
    doi: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    root = data_root / "PXD053474/search_results/extracted"
    rows: list[dict[str, object]] = []
    logs: list[dict[str, object]] = []
    for path in sorted(root.rglob("*ptm-site*Report.tsv")):
        table = pd.read_csv(path, sep="\t", dtype=str, low_memory=False).fillna("")
        prefixes = run_prefixes(table.columns)
        enriched = "Enriched-DIA" in str(path)
        compartment = path.parent.name
        no_lac = 0
        for index, row in table.iterrows():
            if clean_text(row.get("PTM.ModificationTitle", "")).casefold() != "lac" or clean_text(row.get("PTM.SiteAA", "")).upper() != "K":
                no_lac += 1
                continue
            if clean_text(row.get("PG.Organisms", "")) and "Homo sapiens" not in clean_text(row.get("PG.Organisms", "")):
                continue
            accession = normalize_accession(row.get("PTM.ProteinId", ""))
            position = integer(row.get("PTM.SiteLocation", ""))
            if not accession or position is None:
                continue
            genes = split_tokens(row.get("PG.Genes", ""))
            names = split_tokens(row.get("PG.ProteinDescriptions", ""))
            emitted = False
            for prefix, sample in prefixes:
                if not is_true(row.get(f"{prefix}.PG.IsIdentified", "")):
                    continue
                emitted = True
                localization = number(row.get(f"{prefix}.PTM.SiteProbability", ""))
                rows.append(
                    blank_site(
                        PXD="PXD053474",
                        DOI=doi,
                        SampleName=sample,
                        CellOrTissueType="HCT116",
                        ExperimentalGroup=compartment,
                        Replicate=sample,
                        EnrichmentStatus="Enriched" if enriched else "Unenriched",
                        AcquisitionMode="DIA",
                        Accession=accession,
                        BaseAccession=base_accession(accession),
                        GeneSymbol=genes[0] if genes else "",
                        ProteinName=names[0] if names else "",
                        KlaSite=f"K{position}",
                        ModifiedPeptide=clean_text(row.get(f"{prefix}.PTM.Group", row.get("PTM.FlankingRegion", ""))),
                        LocalizationProb=localization,
                        PEP=number(row.get(f"{prefix}.PG.PEP (Run-Wise)", row.get("PG.PEP", ""))),
                        Score=number(row.get(f"{prefix}.PG.Cscore (Run-Wise)", row.get("PG.Cscore", ""))),
                        ClassI=localization is not None and localization >= 0.75,
                        SourceFile=relative_path(path, project_root),
                        SourceRow=index + 2,
                        SiteID=clean_text(row.get("PTM.CollapseKey", "")) or f"DIA:{path.parent.name}:{index + 2}",
                        EvidenceMode="DIA_ptm_site_report",
                        SourceConfidence="high_localized" if localization is not None and localization >= 0.75 else "standard_search_result",
                    )
                )
            if not emitted:
                rows.append(
                    blank_site(
                        PXD="PXD053474",
                        DOI=doi,
                        SampleName=f"{compartment}_experiment_wide",
                        CellOrTissueType="HCT116",
                        ExperimentalGroup=compartment,
                        Replicate="experiment-wide",
                        EnrichmentStatus="Enriched" if enriched else "Unenriched",
                        AcquisitionMode="DIA",
                        Accession=accession,
                        BaseAccession=base_accession(accession),
                        GeneSymbol=genes[0] if genes else "",
                        ProteinName=names[0] if names else "",
                        KlaSite=f"K{position}",
                        ModifiedPeptide=clean_text(row.get("PTM.FlankingRegion", "")),
                        PEP=number(row.get("PG.PEP", "")),
                        Score=number(row.get("PG.Cscore", "")),
                        SourceFile=relative_path(path, project_root),
                        SourceRow=index + 2,
                        SiteID=clean_text(row.get("PTM.CollapseKey", "")) or f"DIA:{path.parent.name}:{index + 2}",
                        EvidenceMode="DIA_ptm_site_report_experiment_wide",
                        SourceConfidence="standard_search_result",
                    )
                )
        if no_lac:
            logs.append(exclusion_row("PXD053474", relative_path(path, project_root), "non_lactyl_PTM_rows", no_lac))
    return dataframe(rows), pd.DataFrame(logs)


def extract_pxd053474_supplementary(
    data_root: Path,
    project_root: Path,
    doi: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    path = data_root / "PXD053474/supplementary/js4c00366_si_003.xlsx"
    rows: list[dict[str, object]] = []
    logs: list[dict[str, object]] = []
    if not path.exists():
        return dataframe(rows), pd.DataFrame([exclusion_row("PXD053474", str(path), "missing_author_S3")])
    for sheet, accession_column in (("Subcellular", "Accession (553)"), ("Whole-cell lysates", "Accession")):
        table = pd.read_excel(path, sheet_name=sheet, header=1, dtype=str).fillna("")
        for index, row in table.iterrows():
            accession = normalize_accession(row.get(accession_column, ""))
            positions = [int(value) for value in re.findall(r"(?:K)?(\d+)", clean_text(row.get("Sites", "")), re.I)]
            if not accession or not positions:
                continue
            for position in positions:
                rows.append(
                    blank_site(
                        PXD="PXD053474",
                        DOI=doi,
                        SampleName=f"Author_S3_{sheet}",
                        CellOrTissueType="HCT116",
                        ExperimentalGroup=sheet,
                        Replicate="author table aggregate",
                        EnrichmentStatus="Enriched",
                        AcquisitionMode="DDA/DIA author union",
                        Accession=accession,
                        BaseAccession=base_accession(accession),
                        GeneSymbol="",
                        ProteinName="",
                        KlaSite=f"K{position}",
                        ModifiedPeptide=clean_text(row.get("Peptide", "")),
                        ClassI=True,
                        SourceFile=relative_path(path, project_root),
                        SourceRow=index + 3,
                        SiteID=clean_text(row.get("ID", "")) or f"S3:{sheet}:{index + 3}:{accession}:K{position}",
                        EvidenceMode="author_supplementary_table",
                        SourceConfidence="high_author_reported",
                    )
                )
    return dataframe(rows), pd.DataFrame(logs)


def reconcile_pxd053474(
    dda: pd.DataFrame,
    dia: pd.DataFrame,
    supplementary: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    evidence = pd.concat([dda, dia, supplementary], ignore_index=True)
    dda_keys = set(zip(dda["BaseAccession"], dda["KlaSite"]))
    dia_keys = set(zip(dia["BaseAccession"], dia["KlaSite"]))
    s3_keys = set(zip(supplementary["BaseAccession"], supplementary["KlaSite"]))
    all_keys = sorted(dda_keys | dia_keys | s3_keys)
    comparison_rows = []
    support_map: dict[tuple[str, str], dict[str, object]] = {}
    for key in all_keys:
        in_dda = key in dda_keys
        in_dia = key in dia_keys
        in_s3 = key in s3_keys
        in_search = in_dda or in_dia
        if in_search and in_s3:
            comparison = "search_and_supplementary"
        elif in_search:
            comparison = "search_only"
        else:
            comparison = "supplementary_only"
        primary = in_s3 or (in_dda and in_dia)
        if in_s3 and in_search:
            confidence = "high_search_and_author_confirmed"
        elif in_s3:
            confidence = "high_author_reported"
        elif in_dda and in_dia:
            confidence = "high_cross_mode_search_support"
        else:
            confidence = "moderate_single_mode_search_only"
        support_map[key] = {
            "PXD053_DDA": in_dda,
            "PXD053_DIA": in_dia,
            "PXD053_S3": in_s3,
            "PXD053Comparison": comparison,
            "PrimaryIncluded": primary,
            "SourceConfidence": confidence,
            "InclusionStatus": "included" if primary else "audit_only",
            "ExclusionReason": "" if primary else "single_mode_search_only_not_in_author_S3",
        }
        comparison_rows.append(
            {
                "BaseAccession": key[0],
                "KlaSite": key[1],
                **support_map[key],
            }
        )
    for column in ("PXD053_DDA", "PXD053_DIA", "PXD053_S3"):
        evidence[column] = False
    evidence["PXD053Comparison"] = ""
    for index, row in evidence.iterrows():
        support = support_map[(row["BaseAccession"], row["KlaSite"])]
        for key, value in support.items():
            evidence.at[index, key] = value
    return evidence, pd.DataFrame(comparison_rows)


def extract_maxquant_site_table(
    pxd: str,
    doi: str,
    path: Path,
    project_root: Path,
    cell_default: str,
    sample_map: dict[str, dict[str, str]],
    evidence_mode: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    table = pd.read_csv(path, sep="\t", dtype=str, low_memory=False).fillna("")
    rows: list[dict[str, object]] = []
    logs: list[dict[str, object]] = []
    reason_counts: dict[str, int] = {}
    labels = sorted(
        {
            column.removeprefix("Identification type ")
            for column in table.columns
            if column.startswith("Identification type ")
        }
    )
    for index, row in table.iterrows():
        reason = ""
        if any(is_true(row.get(column, "")) for column in ("Reverse", "Potential contaminant", "Contaminant")):
            reason = "reverse_or_contaminant"
        if not reason and clean_text(row.get("Amino acid", "K")).upper() != "K":
            reason = "non_lysine_site"
        accession = normalize_accession(row.get("Protein", row.get("Leading proteins", "")))
        position = integer(row.get("Position", row.get("Positions within proteins", "")))
        if not reason and (not accession or position is None):
            reason = "missing_accession_or_kla_position"
        if reason:
            reason_counts[reason] = reason_counts.get(reason, 0) + 1
            continue
        gene, protein_name = annotation_from_fasta(row.get("Fasta headers", ""), accession)
        if not gene:
            genes = split_tokens(row.get("Gene names", ""))
            gene = genes[0] if genes else ""
        if not protein_name:
            names = split_tokens(row.get("Protein names", ""))
            protein_name = names[0] if names else ""
        detected = [
            label
            for label in labels
            if clean_text(row.get(f"Identification type {label}", ""))
        ]
        if not detected:
            if labels:
                reason_counts["site_not_identified_in_any_sample"] = reason_counts.get(
                    "site_not_identified_in_any_sample", 0
                ) + 1
                continue
            detected = [cell_default]
        for token in detected:
            mapped = sample_map.get(token, {})
            localization = number(
                row.get(f"Localization prob {token}", row.get("Localization prob", ""))
            )
            rows.append(
                blank_site(
                    PXD=pxd,
                    DOI=doi,
                    SampleName=mapped.get("SampleName", token),
                    CellOrTissueType=mapped.get("CellType", cell_default),
                    ExperimentalGroup=mapped.get("ExperimentalGroup", token),
                    Replicate=mapped.get("Replicate", token),
                    EnrichmentStatus=mapped.get("EnrichmentStatus", "Enriched"),
                    AcquisitionMode=mapped.get("AcquisitionMode", "DDA"),
                    Accession=accession,
                    BaseAccession=base_accession(accession),
                    GeneSymbol=gene,
                    ProteinName=protein_name,
                    KlaSite=f"K{position}",
                    ModifiedPeptide=clean_text(row.get("La(K) Probabilities", row.get("La (K) Probabilities", row.get("Modified sequence", "")))),
                    LocalizationProb=localization,
                    PEP=number(row.get(f"PEP {token}", row.get("PEP", ""))),
                    Score=number(row.get(f"Score {token}", row.get("Score", ""))),
                    DiagnosticIon=clean_text(row.get("Diagnostic peak", "")),
                    ClassI=localization is not None and localization >= 0.75,
                    SourceFile=relative_path(path, project_root),
                    SourceRow=index + 2,
                    SiteID=clean_text(row.get("id", "")) or f"{pxd}:{index + 2}:{accession}:K{position}",
                    EvidenceMode=evidence_mode,
                    SourceConfidence="high_localized" if localization is not None and localization >= 0.75 else "standard_author_search_result",
                )
            )
    for reason, count in reason_counts.items():
        logs.append(exclusion_row(pxd, relative_path(path, project_root), reason, count))
    return dataframe(rows), pd.DataFrame(logs)


def extract_pxd078013(
    data_root: Path,
    project_root: Path,
    doi: str,
    sample_map: dict[str, dict[str, str]],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    root = data_root / "PXD078013/search_results"
    evidence_path = root / "evidence.txt"
    proteins_path = root / "proteinGroups.txt"
    evidence = pd.read_csv(evidence_path, sep="\t", dtype=str, low_memory=False).fillna("")
    proteins = pd.read_csv(proteins_path, sep="\t", dtype=str, low_memory=False).fillna("")
    evidence_by_site: dict[str, list[dict[str, object]]] = {}
    logs: list[dict[str, object]] = []
    reason_counts: dict[str, int] = {}
    for index, row in evidence.iterrows():
        if any(is_true(row.get(column, "")) for column in ("Reverse", "Potential contaminant")):
            reason_counts["reverse_or_contaminant"] = reason_counts.get("reverse_or_contaminant", 0) + 1
            continue
        if (number(row.get("La (K)", "")) or 0) <= 0:
            reason_counts["evidence_without_positive_LaK"] = reason_counts.get("evidence_without_positive_LaK", 0) + 1
            continue
        site_ids = split_tokens(row.get("La (K) site IDs", ""))
        if not site_ids:
            reason_counts["positive_LaK_without_site_ID"] = reason_counts.get("positive_LaK_without_site_ID", 0) + 1
            continue
        probabilities = parse_probability_values(row.get("La (K) Probabilities", ""))
        entry = {
            "EvidenceRow": index + 2,
            "SampleToken": clean_text(row.get("Experiment", row.get("Raw file", ""))),
            "ModifiedPeptide": clean_text(row.get("Modified sequence", "")),
            "LocalizationProb": max(probabilities) if probabilities else None,
            "PEP": number(row.get("PEP", "")),
            "Score": number(row.get("Score", "")),
        }
        for site_id in site_ids:
            evidence_by_site.setdefault(site_id, []).append(entry)
    rows: list[dict[str, object]] = []
    for protein_index, row in proteins.iterrows():
        if any(is_true(row.get(column, "")) for column in ("Reverse", "Potential contaminant", "Only identified by site")):
            continue
        site_ids = split_tokens(row.get("La (K) site IDs", ""))
        positions = split_tokens(row.get("La (K) site positions", ""))
        if not site_ids or not positions:
            continue
        pairs = site_position_pairs(site_ids, positions)
        if not pairs:
            reason_counts["site_id_position_count_mismatch"] = reason_counts.get(
                "site_id_position_count_mismatch", 0
            ) + 1
            continue
        accessions = split_tokens(row.get("Majority protein IDs", "")) or split_tokens(row.get("Protein IDs", ""))
        accession = normalize_accession(accessions[0] if accessions else "")
        if not accession:
            continue
        gene, protein_name = annotation_from_fasta(row.get("Fasta headers", ""), accession)
        if not gene:
            genes = split_tokens(row.get("Gene names", ""))
            gene = genes[0] if genes else ""
        if not protein_name:
            names = split_tokens(row.get("Protein names", ""))
            protein_name = names[0] if names else ""
        for site_id, position_raw in pairs:
            position = integer(position_raw)
            if position is None:
                continue
            entries = evidence_by_site.get(site_id, [])
            if not entries:
                reason_counts["proteinGroups_site_without_positive_evidence"] = reason_counts.get("proteinGroups_site_without_positive_evidence", 0) + 1
                continue
            for entry in entries:
                token = clean_text(entry["SampleToken"])
                mapped = sample_map.get(token, {})
                rows.append(
                    blank_site(
                        PXD="PXD078013",
                        DOI=doi,
                        SampleName=mapped.get("SampleName", token),
                        CellOrTissueType="RKO",
                        ExperimentalGroup=mapped.get("ExperimentalGroup", token),
                        Replicate=mapped.get("Replicate", token),
                        EnrichmentStatus="Enriched",
                        AcquisitionMode="DDA",
                        Accession=accession,
                        BaseAccession=base_accession(accession),
                        GeneSymbol=gene,
                        ProteinName=protein_name,
                        KlaSite=f"K{position}",
                        ModifiedPeptide=entry["ModifiedPeptide"],
                        LocalizationProb=entry["LocalizationProb"],
                        PEP=entry["PEP"],
                        Score=entry["Score"],
                        ClassI=entry["LocalizationProb"] is not None and float(entry["LocalizationProb"]) >= 0.75,
                        SourceFile=f"{relative_path(proteins_path, project_root)};{relative_path(evidence_path, project_root)}",
                        SourceRow=f"proteinGroups:{protein_index + 2};evidence:{entry['EvidenceRow']}",
                        SiteID=site_id,
                        EvidenceMode="derived_proteinGroups_plus_Kla_positive_evidence",
                        SourceConfidence="moderate_two_table_derived_site",
                    )
                )
    for reason, count in reason_counts.items():
        logs.append(exclusion_row("PXD078013", relative_path(evidence_path, project_root), reason, count))
    return dataframe(rows), pd.DataFrame(logs)
