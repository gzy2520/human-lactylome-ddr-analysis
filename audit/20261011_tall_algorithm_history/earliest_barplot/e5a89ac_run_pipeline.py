#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path
import platform
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib_venn import venn3
import pandas as pd

from common import (
    apply_annotation_supplement,
    best_annotation,
    clean_text,
    read_go_annotations,
    unique_join,
)
from extractors import (
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
from audit_target_sources import build_target_source_audit


INCLUDED_PXDS = [
    "PXD014870",
    "PXD028488",
    "PXD050470",
    "PXD053474",
    "PXD060185",
    "PXD078013",
    "PXD078736",
]
EXCLUDED_PXDS = ["PXD038880", "PXD050906"]
CATEGORIES = [
    "hippocampus_tissue",
    "normal_immortalized_cell_lines",
    "tumor_cell_lines",
]
REGION_ORDER = [
    "hippocampus_only",
    "normal_only",
    "tumor_only",
    "hippocampus_and_normal_only",
    "hippocampus_and_tumor_only",
    "normal_and_tumor_only",
    "all_three",
]


def write_csv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, index=False, encoding="utf-8-sig")


def as_bool(series: pd.Series) -> pd.Series:
    if series.dtype == bool:
        return series.fillna(False)
    return series.fillna("").astype(str).str.casefold().isin({"true", "1", "yes", "+"})


def sample_mapping(config: pd.DataFrame, pxd: str) -> dict[str, dict[str, str]]:
    subset = config[config["PXD"].eq(pxd)]
    return {row["SourceToken"]: row.to_dict() for _, row in subset.iterrows()}


def ordered_sites(values: pd.Series) -> str:
    sites = set()
    for value in values:
        sites.update(token for token in clean_text(value).split(";") if token)
    return ";".join(sorted(sites, key=lambda value: (int(value[1:]) if value.startswith("K") and value[1:].isdigit() else 10**9, value)))


def max_number(values: pd.Series) -> float | None:
    numeric = pd.to_numeric(values, errors="coerce").dropna()
    return None if numeric.empty else float(numeric.max())


def aggregate_evidence(frame: pd.DataFrame) -> pd.DataFrame:
    columns = [
        "BaseAccession",
        "GeneSymbol",
        "ProteinName",
        "KlaSites",
        "PXD",
        "Sample",
        "CellType",
        "ExperimentalGroup",
        "Category",
        "SourceFile",
        "EvidenceMode",
        "LocalizationProb",
        "SourceConfidence",
        "EvidenceRows",
    ]
    if frame.empty:
        return pd.DataFrame(columns=columns)
    rows = []
    for accession, group in frame.groupby("BaseAccession", sort=True):
        rows.append(
            {
                "BaseAccession": accession,
                "GeneSymbol": best_annotation(group["GeneSymbol"]),
                "ProteinName": best_annotation(group["ProteinName"]),
                "KlaSites": ordered_sites(group["KlaSite"]),
                "PXD": unique_join(group["PXD"]),
                "Sample": unique_join(group["SampleName"]),
                "CellType": unique_join(group["CellOrTissueType"]),
                "ExperimentalGroup": unique_join(group["ExperimentalGroup"]),
                "Category": unique_join(group["Category"]) if "Category" in group else "",
                "SourceFile": unique_join(group["SourceFile"]),
                "EvidenceMode": unique_join(group["EvidenceMode"]),
                "LocalizationProb": max_number(group["LocalizationProb"]),
                "SourceConfidence": unique_join(group["SourceConfidence"]),
                "EvidenceRows": len(group),
            }
        )
    return pd.DataFrame(rows, columns=columns)


def cell_type_statistics(
    evidence: pd.DataFrame,
    go_accessions: set[str],
    grouping: pd.DataFrame,
) -> pd.DataFrame:
    rows = []
    category_order = {category: index for index, category in enumerate(CATEGORIES)}
    ordered_grouping = (
        grouping.assign(
            _category_order=grouping["teacher_requested_grouping"].map(category_order)
        )
        .sort_values("_category_order", kind="stable")
    )
    for cell_type in ordered_grouping["CellType"]:
        subset = evidence[evidence["CellOrTissueType"].eq(cell_type)]
        if subset.empty:
            continue
        accessions = set(subset["BaseAccession"])
        ddr_accessions = accessions & go_accessions
        total = len(accessions)
        rows.append(
            {
                "CellOrTissueType": cell_type,
                "TotalKlaProteins": total,
                "KlaGoDdrProteins": len(ddr_accessions),
                "KlaGoDdrFraction": len(ddr_accessions) / total if total else 0.0,
            }
        )
    return pd.DataFrame(rows)


def extract_all(project_root: Path) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, pd.DataFrame]]:
    data_root = project_root / "data"
    config_root = project_root / "reanalysis/config"
    dataset_config = pd.read_csv(config_root / "datasets.csv", dtype=str).fillna("")
    sample_config = pd.read_csv(config_root / "sample_map.csv", dtype=str).fillna("")
    doi = dataset_config.set_index("PXD")["DOI"].to_dict()
    frames: list[pd.DataFrame] = []
    logs: list[pd.DataFrame] = []
    audits: dict[str, pd.DataFrame] = {}

    sites, log = extract_pxd014870(data_root, project_root, doi["PXD014870"])
    frames.append(sites)
    logs.append(log)

    sites, log, directory_audit = extract_pxd028488(data_root, project_root, doi["PXD028488"])
    frames.append(sites)
    logs.append(log)
    audits["pxd028_directory_audit"] = directory_audit

    sites, log = extract_pxd050470(data_root, project_root, doi["PXD050470"])
    frames.append(sites)
    logs.append(log)

    dda, dda_log = extract_pxd053474_dda(data_root, project_root, doi["PXD053474"])
    dia, dia_log = extract_pxd053474_dia(data_root, project_root, doi["PXD053474"])
    supplementary, supplementary_log = extract_pxd053474_supplementary(
        data_root, project_root, doi["PXD053474"]
    )
    sites, comparison = reconcile_pxd053474(dda, dia, supplementary)
    frames.append(sites)
    logs.extend([dda_log, dia_log, supplementary_log])
    audits["pxd053_comparison"] = comparison
    audits["pxd053_dda"] = dda
    audits["pxd053_dia"] = dia
    audits["pxd053_supplementary"] = supplementary

    sites, log = extract_maxquant_site_table(
        "PXD060185",
        doi["PXD060185"],
        data_root / "PXD060185/search_results/RESULT/combined/txt/La (K)Sites.txt",
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

    all_evidence = pd.concat(frames, ignore_index=True, sort=False)
    all_evidence = apply_annotation_supplement(
        all_evidence, config_root / "uniprot_annotation_supplement.tsv"
    )
    all_evidence["PrimaryIncluded"] = as_bool(all_evidence["PrimaryIncluded"])
    all_evidence["ClassI"] = as_bool(all_evidence["ClassI"])
    all_evidence = all_evidence.sort_values(
        ["PXD", "CellOrTissueType", "SampleName", "BaseAccession", "KlaSite", "SourceFile", "SourceRow"],
        kind="stable",
    ).reset_index(drop=True)
    exclusion_log = pd.concat([item for item in logs if not item.empty], ignore_index=True)
    return all_evidence, exclusion_log, audits


def attach_go(proteins: pd.DataFrame, go_summary: pd.DataFrame, go_raw: pd.DataFrame) -> pd.DataFrame:
    go_columns = [
        "GOSymbol",
        "GOTerms",
        "GONames",
        "GOEvidenceCodes",
        "GOReferences",
        "GOAnnotationCount",
    ]
    merged = proteins.merge(go_summary, on="BaseAccession", how="left")
    for column in go_columns:
        if column not in merged:
            merged[column] = ""
    merged["GOMatchMode"] = merged["GONames"].fillna("").ne("").map(
        {True: "BaseAccession", False: "unmatched"}
    )

    taxon = pd.to_numeric(go_raw["TAXON ID"], errors="coerce")
    retained = go_raw[
        (~go_raw["ExcludedNOT"])
        & (taxon.isna() | taxon.eq(9606))
        & go_raw["BaseAccession"].ne("")
    ].copy()
    retained["SymbolKey"] = retained["SYMBOL"].str.upper().str.strip()
    gene_rows = []
    for symbol, group in retained[retained["SymbolKey"].ne("")].groupby("SymbolKey", sort=True):
        gene_rows.append(
            {
                "SymbolKey": symbol,
                "GOSymbol": best_annotation(group["SYMBOL"]),
                "GOTerms": unique_join(group["GO TERM"]),
                "GONames": unique_join(group["GO NAME"]),
                "GOEvidenceCodes": unique_join(group["GO EVIDENCE CODE"]),
                "GOReferences": unique_join(group["REFERENCE"]),
                "GOAnnotationCount": len(group),
            }
        )
    gene_map = {row["SymbolKey"]: row for row in gene_rows}
    for index in merged.index[merged["GOMatchMode"].eq("unmatched")]:
        symbols = [token.strip().upper() for token in clean_text(merged.at[index, "GeneSymbol"]).split(";") if token.strip()]
        match = next((gene_map[symbol] for symbol in symbols if symbol in gene_map), None)
        if match is None:
            continue
        for column in go_columns:
            merged.at[index, column] = match[column]
        merged.at[index, "GOMatchMode"] = "GeneSymbol_fallback"
    return merged


def venn_regions(sets: dict[str, set[str]]) -> dict[str, set[str]]:
    h = sets["hippocampus_tissue"]
    n = sets["normal_immortalized_cell_lines"]
    t = sets["tumor_cell_lines"]
    return {
        "hippocampus_only": h - n - t,
        "normal_only": n - h - t,
        "tumor_only": t - h - n,
        "hippocampus_and_normal_only": (h & n) - t,
        "hippocampus_and_tumor_only": (h & t) - n,
        "normal_and_tumor_only": (n & t) - h,
        "all_three": h & n & t,
    }


def plot_venn(regions: dict[str, set[str]], output_stem: Path, title: str) -> None:
    counts = {name: len(values) for name, values in regions.items()}
    fig, axis = plt.subplots(figsize=(7.2, 6.2))
    diagram = venn3(
        subsets=(
            counts["hippocampus_only"],
            counts["normal_only"],
            counts["hippocampus_and_normal_only"],
            counts["tumor_only"],
            counts["hippocampus_and_tumor_only"],
            counts["normal_and_tumor_only"],
            counts["all_three"],
        ),
        set_labels=("Hippocampus tissue", "Immortalized models", "Tumor cell lines"),
        set_colors=("#2A6F97", "#E9C46A", "#C14953"),
        alpha=0.62,
        ax=axis,
    )
    if diagram.set_labels:
        for label in diagram.set_labels:
            if label is not None:
                label.set_fontsize(10)
    if diagram.subset_labels:
        for label in diagram.subset_labels:
            if label is not None:
                label.set_fontsize(10)
    axis.set_title(title, fontsize=13, pad=16)
    output_stem.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_stem.with_suffix(".png"), dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(output_stem.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
    fig.savefig(output_stem.with_suffix(".svg"), bbox_inches="tight", facecolor="white")
    plt.close(fig)


def write_group_outputs(
    evidence: pd.DataFrame,
    grouping: pd.DataFrame,
    scheme: str,
    analysis_name: str,
    table_root: Path,
    figure_root: Path,
) -> tuple[pd.DataFrame, dict[str, set[str]], pd.DataFrame, pd.DataFrame]:
    mapping = grouping.set_index("CellType")[scheme].to_dict()
    assigned = evidence.copy()
    assigned["Category"] = assigned["CellOrTissueType"].map(mapping).fillna("")
    unmapped = assigned[assigned["Category"].eq("")][
        ["PXD", "CellOrTissueType", "SampleName", "BaseAccession"]
    ].drop_duplicates()
    if not unmapped.empty:
        raise ValueError(f"Unmapped cell types for {scheme}: {sorted(unmapped['CellOrTissueType'].unique())}")
    sets = {
        category: set(assigned.loc[assigned["Category"].eq(category), "BaseAccession"])
        for category in CATEGORIES
    }
    regions = venn_regions(sets)
    output_dir = table_root / "venn_regions" / scheme / analysis_name
    counts = []
    membership_rows = []
    for region_name in REGION_ORDER:
        accessions = regions[region_name]
        region_frame = aggregate_evidence(assigned[assigned["BaseAccession"].isin(accessions)])
        write_csv(region_frame, output_dir / f"{region_name}.csv")
        counts.append({"Region": region_name, "ProteinCount": len(accessions)})
        membership_rows.extend(
            {"BaseAccession": accession, "Region": region_name} for accession in sorted(accessions)
        )
    count_frame = pd.DataFrame(counts)
    write_csv(count_frame, output_dir / "venn_region_counts.csv")
    write_csv(pd.DataFrame(membership_rows), output_dir / "venn_membership.csv")
    combined = aggregate_evidence(assigned)
    region_by_accession = {
        accession: region_name
        for region_name, accessions in regions.items()
        for accession in accessions
    }
    membership_columns = {
        "InHippocampusTissue": "hippocampus_tissue",
        "InNormalImmortalizedCellLines": "normal_immortalized_cell_lines",
        "InTumorCellLines": "tumor_cell_lines",
    }
    for column, category in membership_columns.items():
        combined[column] = combined["BaseAccession"].isin(sets[category]).map(
            {True: "Yes", False: "No"}
        )
    combined["DetectedGroupCount"] = combined[list(membership_columns)].eq("Yes").sum(axis=1)
    combined["VennRegion"] = combined["BaseAccession"].map(region_by_accession).fillna("")
    combined["AnalysisSet"] = analysis_name
    leading_columns = [
        "BaseAccession",
        "GeneSymbol",
        "ProteinName",
        "KlaSites",
        "InHippocampusTissue",
        "InNormalImmortalizedCellLines",
        "InTumorCellLines",
        "DetectedGroupCount",
        "VennRegion",
        "AnalysisSet",
    ]
    combined = combined[leading_columns + [column for column in combined if column not in leading_columns]]
    write_csv(combined, output_dir / f"{analysis_name}_three_groups_combined.csv")
    category_frames = []
    for category in CATEGORIES:
        category_frame = aggregate_evidence(
            assigned[assigned["BaseAccession"].isin(sets[category])]
        )
        write_csv(category_frame, output_dir / f"{category}_all.csv")
        category_frame.insert(4, "SourceCategory", category)
        category_frames.append(category_frame)
    non_deduplicated = pd.concat(category_frames, ignore_index=True)
    for column, category in membership_columns.items():
        non_deduplicated[column] = non_deduplicated["BaseAccession"].isin(
            sets[category]
        ).map({True: "Yes", False: "No"})
    non_deduplicated["DetectedGroupCount"] = non_deduplicated[
        list(membership_columns)
    ].eq("Yes").sum(axis=1)
    non_deduplicated["VennRegion"] = non_deduplicated["BaseAccession"].map(
        region_by_accession
    ).fillna("")
    non_deduplicated["AnalysisSet"] = analysis_name
    non_deduplicated_leading = leading_columns[:4] + ["SourceCategory"] + leading_columns[4:]
    non_deduplicated = non_deduplicated[
        non_deduplicated_leading
        + [column for column in non_deduplicated if column not in non_deduplicated_leading]
    ]
    write_csv(
        non_deduplicated,
        output_dir / f"{analysis_name}_three_groups_combined_non_deduplicated.csv",
    )
    plot_venn(
        regions,
        figure_root / scheme / f"{analysis_name}_three_group_venn",
        "All Kla proteins"
        if analysis_name == "all_kla"
        else "Kla and GO repair/damage proteins",
    )
    return count_frame, regions, combined, non_deduplicated


def target_trace(all_evidence: pd.DataFrame) -> pd.DataFrame:
    targets = {
        "MRE11": {"genes": {"MRE11"}, "accessions": {"P49959"}},
        "XLF/NHEJ1": {"genes": {"NHEJ1", "XLF"}, "accessions": {"Q9H9Q4"}},
        "NBS1/NBN": {"genes": {"NBN", "NBS1"}, "accessions": {"O60934"}},
    }
    rows = []
    for pxd in INCLUDED_PXDS:
        dataset = all_evidence[all_evidence["PXD"].eq(pxd)]
        for target, aliases in targets.items():
            gene_match = dataset["GeneSymbol"].fillna("").map(
                lambda value: bool({token.strip().upper() for token in str(value).split(";")} & aliases["genes"])
            )
            accession_match = dataset["BaseAccession"].isin(aliases["accessions"])
            matched = dataset[gene_match | accession_match]
            primary = matched[matched["PrimaryIncluded"]]
            if not primary.empty:
                status = "present_in_primary_kla_evidence"
            elif not matched.empty:
                status = "present_only_in_audit_kla_evidence"
            else:
                status = "not_present_in_extracted_kla_evidence"
            rows.append(
                {
                    "PXD": pxd,
                    "Target": target,
                    "TargetGenes": ";".join(sorted(aliases["genes"])),
                    "TargetAccessions": ";".join(sorted(aliases["accessions"])),
                    "Status": status,
                    "PrimaryEvidenceRows": len(primary),
                    "AllEvidenceRows": len(matched),
                    "KlaSites": ordered_sites(matched["KlaSite"]) if not matched.empty else "",
                    "Samples": unique_join(matched["SampleName"]) if not matched.empty else "",
                    "SourceFiles": unique_join(matched["SourceFile"]) if not matched.empty else "",
                }
            )
    return pd.DataFrame(rows)


def regression_outputs(new_proteins: pd.DataFrame, project_root: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    old_path = project_root / "archive/reanalysis_v1_2026-07-21/outputs/01_kla_extraction/kla_proteins_by_dataset_cell_type_group.csv"
    if not old_path.exists():
        return pd.DataFrame(), pd.DataFrame()
    old = pd.read_csv(old_path, dtype=str).fillna("")
    old_keys = set(zip(old["PXD"], old["BaseAccession"]))
    new_keys = set(zip(new_proteins["PXD"], new_proteins["BaseAccession"]))
    rows = []
    for key in sorted(old_keys | new_keys):
        if key in old_keys and key in new_keys:
            status = "shared"
            reason = "present_in_both"
        elif key in new_keys:
            status = "new_only"
            reason = {
                "PXD028488": "expanded_all_human_PEAKS_directories",
                "PXD050470": "corrected_supplementary_header_and_S3_S12_union",
                "PXD053474": "four_search_suites_plus_reconciled_S3",
                "PXD078013": "strict_evidence_plus_proteinGroups_site_link",
            }.get(key[0], "revised_sample_level_extraction")
        else:
            status = "old_only"
            reason = "removed_by_revised_traceable_evidence_rule_or_accession_normalization"
        rows.append({"PXD": key[0], "BaseAccession": key[1], "RegressionStatus": status, "Reason": reason})
    detail = pd.DataFrame(rows)
    summary = detail.groupby(["PXD", "RegressionStatus"], as_index=False).size().rename(columns={"size": "ProteinCount"})
    return detail, summary


def main() -> int:
    parser = argparse.ArgumentParser(description="Run the traceable Kla reanalysis pipeline.")
    parser.add_argument("--project-root", type=Path, default=Path(__file__).resolve().parents[2])
    args = parser.parse_args()
    project_root = args.project_root.resolve()
    reanalysis = project_root / "reanalysis"
    intermediate = reanalysis / "intermediate"
    tables = reanalysis / "results/tables"
    figures = reanalysis / "results/figures"
    reports = reanalysis / "reports"
    logs_dir = reanalysis / "logs"
    for path in (intermediate, tables, figures, reports, logs_dir):
        path.mkdir(parents=True, exist_ok=True)

    all_evidence, exclusion_log, audits = extract_all(project_root)
    primary = all_evidence[all_evidence["PrimaryIncluded"]].copy()
    if set(primary["PXD"].unique()) & set(EXCLUDED_PXDS):
        raise AssertionError("Excluded PXD entered primary evidence")
    missing_pxd = set(INCLUDED_PXDS) - set(primary["PXD"].unique())
    if missing_pxd:
        raise AssertionError(f"Included PXD produced no primary evidence: {sorted(missing_pxd)}")

    write_csv(all_evidence, intermediate / "kla_by_dataset/all_included_and_audit_kla_evidence.csv")
    write_csv(primary, intermediate / "kla_by_dataset/all_primary_sample_level_kla_sites.csv")
    write_csv(exclusion_log, tables / "exclusion_log.csv")
    for pxd in INCLUDED_PXDS:
        dataset = primary[primary["PXD"].eq(pxd)].copy()
        write_csv(dataset, intermediate / f"kla_by_dataset/{pxd}_sample_level_kla_sites.csv")
        write_csv(aggregate_evidence(dataset), tables / f"per_pxd/{pxd}_unique_kla_proteins.csv")
    pxd014_sensitivity = primary[
        primary["PXD"].eq("PXD014870")
        & pd.to_numeric(primary["LocalizationProb"], errors="coerce").ge(0.75)
    ]
    write_csv(
        pxd014_sensitivity,
        intermediate / "kla_by_dataset/PXD014870_sensitivity_localization_0.75.csv",
    )

    directory_audit = audits["pxd028_directory_audit"].copy()
    directory_audit["OldWorkflowDirectory"] = directory_audit["Directory"].isin(
        {
            "HCT116-Enrichment-Search files",
            "TALL-NALAC-Search files",
            "HEK293T-Enrichment-all HCD-Search files",
        }
    )
    write_csv(directory_audit, tables / "pxd028488/directory_audit.csv")
    write_csv(directory_audit[directory_audit["PrimaryStatus"].eq("included")], tables / "pxd028488/included_directories.csv")
    write_csv(directory_audit[directory_audit["PrimaryStatus"].ne("included")], tables / "pxd028488/excluded_directories_and_reasons.csv")
    write_csv(directory_audit, tables / "pxd028488/old_vs_new_directory_coverage.csv")
    write_csv(
        all_evidence[
            all_evidence["PXD"].eq("PXD028488")
            & all_evidence["DiagnosticIonIntensity"].notna()
        ],
        tables / "pxd028488/diagnostic_ion_156_supported_evidence.csv",
    )

    comparison = audits["pxd053_comparison"].copy()
    comparison["ComparisonCategory"] = comparison["PXD053Comparison"].map(
        {
            "search_and_supplementary": "consistent",
            "search_only": "search_only",
            "supplementary_only": "supplementary_only",
        }
    )
    write_csv(comparison, tables / "pxd053474/search_vs_supplementary_all.csv")
    for category in ("consistent", "search_only", "supplementary_only"):
        write_csv(comparison[comparison["ComparisonCategory"].eq(category)], tables / f"pxd053474/{category}.csv")
    write_csv(comparison.iloc[0:0], tables / "pxd053474/inconsistent.csv")
    write_csv(
        all_evidence[all_evidence["PXD"].eq("PXD053474") & ~all_evidence["PrimaryIncluded"]],
        tables / "pxd053474/single_mode_search_only_audit_evidence.csv",
    )

    go_summary, go_raw = read_go_annotations(project_root / "data/annotations/GO-repair+damage(human).tsv")
    proteins = aggregate_evidence(primary)
    proteins_by_pxd = pd.concat(
        [aggregate_evidence(primary[primary["PXD"].eq(pxd)]) for pxd in INCLUDED_PXDS],
        ignore_index=True,
    )
    proteins_go = attach_go(proteins_by_pxd, go_summary, go_raw)
    matched = proteins_go[proteins_go["GOMatchMode"].ne("unmatched")].copy()
    unmatched = proteins_go[proteins_go["GOMatchMode"].eq("unmatched")].copy()
    write_csv(proteins, tables / "all_unique_kla_proteins.csv")
    write_csv(matched, intermediate / "go_intersection/all_pxd_kla_go_ddr_proteins.csv")
    write_csv(unmatched, tables / "go_unmatched_kla_proteins.csv")
    write_csv(
        unmatched[
            unmatched["BaseAccession"].eq("") | unmatched["GeneSymbol"].fillna("").eq("")
        ],
        tables / "accession_gene_mapping_failures.csv",
    )
    for pxd in INCLUDED_PXDS:
        write_csv(matched[matched["PXD"].eq(pxd)], intermediate / f"go_intersection/{pxd}_kla_go_ddr.csv")

    grouping = pd.read_csv(reanalysis / "config/grouping_schemes.csv", dtype=str).fillna("")
    review_needed = grouping[
        grouping["review_needed"].str.casefold().eq("yes")
        | grouping["classification_warning"].ne("")
    ]
    write_csv(review_needed, tables / "classification_review_needed.csv")
    go_accessions = set(matched["BaseAccession"])
    go_evidence = primary[primary["BaseAccession"].isin(go_accessions)].copy()
    write_csv(
        cell_type_statistics(primary, go_accessions, grouping),
        tables / "cell_type_kla_ddr_statistics.csv",
    )
    all_counts = []
    region_cache: dict[tuple[str, str], dict[str, set[str]]] = {}
    combined_cache: dict[tuple[str, str], pd.DataFrame] = {}
    non_deduplicated_cache: dict[tuple[str, str], pd.DataFrame] = {}
    for scheme in ("teacher_requested_grouping", "biologically_conventional_grouping"):
        for analysis_name, evidence in (("all_kla", primary), ("kla_go_ddr", go_evidence)):
            count_frame, regions, combined, non_deduplicated = write_group_outputs(
                evidence, grouping, scheme, analysis_name, tables, figures
            )
            count_frame.insert(0, "Analysis", analysis_name)
            count_frame.insert(0, "GroupingScheme", scheme)
            all_counts.append(count_frame)
            region_cache[(scheme, analysis_name)] = regions
            combined_cache[(scheme, analysis_name)] = combined
            non_deduplicated_cache[(scheme, analysis_name)] = non_deduplicated
    write_csv(pd.concat(all_counts, ignore_index=True), tables / "venn_all_schemes_counts.csv")
    write_csv(
        combined_cache[("teacher_requested_grouping", "all_kla")],
        tables / "all_kla_three_groups_combined.csv",
    )
    write_csv(
        combined_cache[("teacher_requested_grouping", "kla_go_ddr")],
        tables / "kla_go_ddr_three_groups_combined.csv",
    )
    write_csv(
        non_deduplicated_cache[("teacher_requested_grouping", "all_kla")],
        tables / "all_kla_three_groups_combined_non_deduplicated.csv",
    )
    write_csv(
        non_deduplicated_cache[("teacher_requested_grouping", "kla_go_ddr")],
        tables / "kla_go_ddr_three_groups_combined_non_deduplicated.csv",
    )

    teacher_all = region_cache[("teacher_requested_grouping", "all_kla")]["tumor_only"]
    teacher_go = region_cache[("teacher_requested_grouping", "kla_go_ddr")]["tumor_only"]
    grouping_map = grouping.set_index("CellType")["teacher_requested_grouping"].to_dict()
    teacher_primary = primary.copy()
    teacher_primary["Category"] = teacher_primary["CellOrTissueType"].map(grouping_map).fillna("")
    write_csv(
        aggregate_evidence(teacher_primary[teacher_primary["BaseAccession"].isin(teacher_all)]),
        tables / "tumor_specific_kla_proteins.csv",
    )
    write_csv(
        aggregate_evidence(teacher_primary[teacher_primary["BaseAccession"].isin(teacher_go)]),
        tables / "tumor_specific_kla_ddr_proteins.csv",
    )

    trace = target_trace(all_evidence)
    write_csv(trace, tables / "target_protein_evidence_trace_MRE11_XLF_NBS1.csv")
    write_csv(
        build_target_source_audit(project_root, primary),
        tables / "target_protein_source_level_audit_MRE11_XLF_NBS1.csv",
    )
    regression_detail, regression_summary = regression_outputs(proteins_by_pxd, project_root)
    write_csv(regression_detail, tables / "regression_old_vs_new_detail.csv")
    write_csv(regression_summary, tables / "regression_old_vs_new_summary.csv")

    dataset_config = pd.read_csv(reanalysis / "config/datasets.csv", dtype=str).fillna("")
    summary_rows = []
    for _, metadata in dataset_config.iterrows():
        pxd = metadata["PXD"]
        subset = primary[primary["PXD"].eq(pxd)]
        go_subset = matched[matched["PXD"].eq(pxd)] if not matched.empty else matched
        summary_rows.append(
            {
                **metadata.to_dict(),
                "PrimarySampleLevelRows": len(subset),
                "UniqueKlaSites": subset[["BaseAccession", "KlaSite"]].drop_duplicates().shape[0],
                "UniqueKlaProteins": subset["BaseAccession"].nunique(),
                "KlaGoDdrProteins": go_subset["BaseAccession"].nunique(),
            }
        )
    write_csv(pd.DataFrame(summary_rows), tables / "dataset_analysis_summary.csv")

    environment = pd.DataFrame(
        [
            {"Component": "Python", "Version": platform.python_version()},
            {"Component": "Platform", "Version": platform.platform()},
            {"Component": "pandas", "Version": pd.__version__},
            {"Component": "matplotlib", "Version": matplotlib.__version__},
        ]
    )
    write_csv(environment, logs_dir / "software_environment.csv")
    print(f"Primary sample-level rows: {len(primary):,}")
    print(f"Unique Kla proteins: {proteins['BaseAccession'].nunique():,}")
    print(f"Kla GO-DDR proteins: {matched['BaseAccession'].nunique():,}")
    print(f"Results: {tables}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
