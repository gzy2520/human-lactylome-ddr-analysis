from __future__ import annotations

from collections import Counter, defaultdict
import csv
import hashlib
import io
import json
from pathlib import Path
import re
import statistics
import zipfile

import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# The workflow invokes this file by path, so Python places the script's
# directory (rather than the repository root) on sys.path. Add the project
# root before importing the shared utility module.
import sys

DEFAULT_PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(DEFAULT_PROJECT_ROOT))

from python.utils.reference_proteome_utils import (
    base_accession,
    best_annotation,
    clean_text,
    read_go_annotations,
    unique_join,
)

parser = argparse.ArgumentParser()
parser.add_argument("--project-root", type=Path, default=DEFAULT_PROJECT_ROOT)
PROJECT_ROOT = parser.parse_known_args()[0].project_root.resolve()
DATA_ROOT = PROJECT_ROOT / "data"
TABLES = PROJECT_ROOT / "results/tables"
FIGURES = PROJECT_ROOT / "results/figures"
INTERMEDIATE = PROJECT_ROOT / "work/intermediate/reference_proteomes"
CONFIG_PATH = PROJECT_ROOT / "config/reference_proteome_selection.json"
GO_PATH = DATA_ROOT / "annotations/GO-repair+damage(human).tsv"
KLA_STATS_PATH = TABLES / "cell_type_kla_ddr_statistics.csv"

MODEL_ORDER = [
    "Human hippocampus",
    "HEK293T",
    "HK-2",
    "MCF10A",
    "MCF7",
    "HCT116",
    "T-ALL",
    "MDA-MB-468",
    "T-47D",
    "RKO",
]

UNIPROT_ACCESSION = re.compile(
    r"^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9])(?:-\d+)?$"
)


def write_csv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, index=False, encoding="utf-8-sig")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def record(
    model: str,
    pxd: str,
    subset: str,
    source_file: Path,
    source_id: str,
    accession: str = "",
    gene: str = "",
    protein_name: str = "",
    evidence_count: int = 1,
    count_unit: str = "BaseAccession",
) -> dict[str, object]:
    accession = base_accession(accession)
    key = accession if accession else f"GENEFEATURE:{source_id}"
    return {
        "CellOrTissueType": model,
        "ReferencePXD": pxd,
        "ReferenceSubset": subset,
        "ReferenceProteinKey": key,
        "BaseAccession": accession,
        "GeneSymbol": clean_text(gene),
        "ProteinName": clean_text(protein_name),
        "EvidenceCount": int(evidence_count),
        "CountUnit": count_unit,
        "SourceFile": str(source_file.relative_to(PROJECT_ROOT)),
        "SourceIdentifier": clean_text(source_id),
    }


def load_pxd030304() -> tuple[list[dict[str, object]], pd.DataFrame]:
    source_root = DATA_ROOT / "PXD030304/search_results"
    averaged_path = source_root / "ProCan-DepMapSanger_protein_matrix_6692_averaged.txt"
    mapping_path = source_root / "ProCan-DepMapSanger_mapping_file_replicates.txt"
    peptide_path = source_root / "ProCan-DepMapSanger_peptide_counts_per_protein_per_sample.txt"

    target_rows = {
        "SIDM00148;MCF7": "MCF7",
        "SIDM00783;HCT-116": "HCT116",
        "SIDM00628;MDA-MB-468": "MDA-MB-468",
        "SIDM00097;T47D": "T-47D",
        "SIDM01090;RKO": "RKO",
        "SIDM00370;TALL-1": "T-ALL",
        "SIDM01016;Jurkat": "T-ALL_Jurkat_sensitivity",
    }
    subset_labels = {
        "MCF7": "MCF7 / SIDM00148",
        "HCT116": "HCT-116 / SIDM00783",
        "MDA-MB-468": "MDA-MB-468 / SIDM00628",
        "T-47D": "T47D / SIDM00097",
        "RKO": "RKO / SIDM01090",
        "T-ALL": "TALL-1 / SIDM00370",
        "T-ALL_Jurkat_sensitivity": "Jurkat / SIDM01016",
    }

    records: list[dict[str, object]] = []
    sensitivity_rows: list[dict[str, object]] = []
    with averaged_path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        proteins = [(item, base_accession(item.split(";", 1)[0])) for item in header[1:]]
        main_accessions = {accession for _, accession in proteins if accession}
        for row in reader:
            if not row or row[0] not in target_rows:
                continue
            model = target_rows[row[0]]
            target = sensitivity_rows if model == "T-ALL_Jurkat_sensitivity" else records
            for index, value in enumerate(row[1:]):
                if clean_text(value) == "":
                    continue
                source_id, accession = proteins[index]
                target.append(
                    record(
                        model,
                        "PXD030304",
                        subset_labels[model],
                        averaged_path,
                        source_id,
                        accession=accession,
                    )
                )

    selected_runs: set[str] = set()
    with mapping_path.open(encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["Cell_line"] == "Control_HEK293T_lys":
                selected_runs.add(row["Automatic_MS_filename"])

    hek_support: dict[str, set[str]] = defaultdict(set)
    with peptide_path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        accessions = [base_accession(item.split(";", 1)[0]) for item in header[1:]]
        keep = [accession in main_accessions for accession in accessions]
        for row in reader:
            if not row or row[0] not in selected_runs:
                continue
            run = row[0]
            for index, value in enumerate(row[1:]):
                if not keep[index] or clean_text(value) == "":
                    continue
                try:
                    present = float(value) > 0
                except ValueError:
                    present = False
                if present:
                    hek_support[accessions[index]].add(run)

    for accession, runs in sorted(hek_support.items()):
        records.append(
            record(
                "HEK293T",
                "PXD030304",
                "Control_HEK293T_lys; standard-QC excluded",
                peptide_path,
                accession,
                accession=accession,
                evidence_count=len(runs),
            )
        )

    sensitivity = pd.DataFrame(sensitivity_rows)
    return records, sensitivity


def load_pxd072220() -> list[dict[str, object]]:
    path = DATA_ROOT / "PXD072220/search_results/HK-2_Spectronaut-report_PG_Quantity.txt"
    support: dict[str, set[str]] = defaultdict(set)
    genes: dict[str, Counter[str]] = defaultdict(Counter)
    with path.open(encoding="utf-8-sig", newline="") as handle:
        handle.readline()
        handle.readline()
        reader = csv.DictReader(handle, delimiter="\t")
        control_columns = [
            column for column in reader.fieldnames or [] if column.endswith(".PG.Log2Quantity")
        ][:3]
        for row in reader:
            present_columns = [
                column
                for column in control_columns
                if clean_text(row.get(column)) not in {"", "NaN", "nan", "NA"}
            ]
            if not present_columns:
                continue
            accessions = [base_accession(item) for item in clean_text(row["PG.ProteinGroups"]).split(";")]
            gene_tokens = [item.strip() for item in clean_text(row["PG.Genes"]).split(";")]
            for index, accession in enumerate(accessions):
                if not accession:
                    continue
                support[accession].update(present_columns)
                if index < len(gene_tokens) and gene_tokens[index]:
                    genes[accession][gene_tokens[index]] += 1

    return [
        record(
            "HK-2",
            "PXD072220",
            "HK-2 Control_1, Control_3 and Control_4",
            path,
            accession,
            accession=accession,
            gene=genes[accession].most_common(1)[0][0] if genes[accession] else "",
            evidence_count=len(columns),
        )
        for accession, columns in sorted(support.items())
    ]


def load_pxd002400() -> list[dict[str, object]]:
    path = DATA_ROOT / "PXD002400/search_results/msms.zip"
    support: dict[str, set[str]] = defaultdict(set)
    genes: dict[str, Counter[str]] = defaultdict(Counter)
    names: dict[str, Counter[str]] = defaultdict(Counter)

    with zipfile.ZipFile(path) as archive, archive.open("evidence.txt") as raw:
        handle = io.TextIOWrapper(raw, encoding="utf-8-sig", newline="")
        for row in csv.DictReader(handle, delimiter="\t"):
            raw_file = clean_text(row.get("Raw file"))
            if not raw_file.startswith("Toni_20111109_FB_10A_"):
                continue
            if clean_text(row.get("Reverse")) == "+" or clean_text(row.get("Contaminant")) == "+":
                continue
            accession = clean_text(row.get("Leading Razor Protein"))
            if not UNIPROT_ACCESSION.fullmatch(accession):
                continue
            accession = base_accession(accession)
            support[accession].add(raw_file)
            gene = clean_text(row.get("Gene Names"))
            protein_name = clean_text(row.get("Protein Names"))
            if gene:
                genes[accession][gene] += 1
            if protein_name:
                names[accession][protein_name] += 1

    return [
        record(
            "MCF10A",
            "PXD002400",
            "MCF-10A baseline; 10 IEF fractions x 2 injections",
            path,
            accession,
            accession=accession,
            gene=genes[accession].most_common(1)[0][0] if genes[accession] else "",
            protein_name=names[accession].most_common(1)[0][0] if names[accession] else "",
            evidence_count=len(raw_files),
        )
        for accession, raw_files in sorted(support.items())
    ]


def load_pxd043880() -> list[dict[str, object]]:
    path = DATA_ROOT / "PXD043880/supplementary/13024_2023_650_MOESM1_ESM.xlsx"
    frame = pd.read_excel(path, sheet_name="Source Data Proteins", header=1)
    protein_columns = list(frame.columns[8:])
    records: list[dict[str, object]] = []
    for column in protein_columns:
        gene_feature = clean_text(column)
        if not gene_feature:
            continue
        symbols = ";".join(token for token in gene_feature.split("_") if token)
        evidence_count = int(pd.to_numeric(frame[column], errors="coerce").notna().sum())
        records.append(
            record(
                "Human hippocampus",
                "PXD043880",
                "74 neurologically normal CA1 hippocampus donors",
                path,
                gene_feature,
                gene=symbols,
                evidence_count=evidence_count,
                count_unit="published_gene/protein_feature",
            )
        )
    return records


def build_go_maps() -> tuple[dict[str, dict[str, object]], dict[str, dict[str, object]]]:
    go_summary, go_raw = read_go_annotations(GO_PATH)
    accession_map = go_summary.set_index("BaseAccession").to_dict("index")

    taxon = pd.to_numeric(go_raw["TAXON ID"], errors="coerce")
    retained = go_raw[
        (~go_raw["ExcludedNOT"])
        & (taxon.isna() | taxon.eq(9606))
        & go_raw["BaseAccession"].ne("")
    ].copy()
    retained["SymbolKey"] = retained["SYMBOL"].str.upper().str.strip()
    gene_map: dict[str, dict[str, object]] = {}
    for symbol, group in retained[retained["SymbolKey"].ne("")].groupby("SymbolKey", sort=True):
        gene_map[symbol] = {
            "GOSymbol": best_annotation(group["SYMBOL"]),
            "GOTerms": unique_join(group["GO TERM"]),
            "GONames": unique_join(group["GO NAME"]),
            "GOEvidenceCodes": unique_join(group["GO EVIDENCE CODE"]),
            "GOReferences": unique_join(group["REFERENCE"]),
            "GOAnnotationCount": len(group),
        }
    return accession_map, gene_map


def annotate_go(frame: pd.DataFrame) -> pd.DataFrame:
    accession_map, gene_map = build_go_maps()
    annotated = frame.copy()
    go_columns = [
        "GOSymbol",
        "GOTerms",
        "GONames",
        "GOEvidenceCodes",
        "GOReferences",
        "GOAnnotationCount",
        "GOMatchMode",
    ]
    for column in go_columns:
        annotated[column] = "" if column != "GOAnnotationCount" else 0

    for index, row in annotated.iterrows():
        accession = clean_text(row["BaseAccession"])
        if accession in accession_map:
            match = accession_map[accession]
            for column in go_columns[:-1]:
                annotated.at[index, column] = match.get(column, "")
            annotated.at[index, "GOMatchMode"] = "BaseAccession"
            continue

        symbols = [
            token.strip().upper()
            for token in re.split(r"[;,\s]+", clean_text(row["GeneSymbol"]))
            if token.strip()
        ]
        matches = [gene_map[symbol] for symbol in symbols if symbol in gene_map]
        if not matches:
            annotated.at[index, "GOMatchMode"] = "unmatched"
            continue
        annotated.at[index, "GOSymbol"] = unique_join(item["GOSymbol"] for item in matches)
        annotated.at[index, "GOTerms"] = unique_join(item["GOTerms"] for item in matches)
        annotated.at[index, "GONames"] = unique_join(item["GONames"] for item in matches)
        annotated.at[index, "GOEvidenceCodes"] = unique_join(
            item["GOEvidenceCodes"] for item in matches
        )
        annotated.at[index, "GOReferences"] = unique_join(item["GOReferences"] for item in matches)
        annotated.at[index, "GOAnnotationCount"] = sum(
            int(item["GOAnnotationCount"]) for item in matches
        )
        annotated.at[index, "GOMatchMode"] = "GeneSymbol_fallback"
    return annotated


def build_reference_statistics(
    annotated: pd.DataFrame, config: list[dict[str, object]]
) -> pd.DataFrame:
    config_by_model = {row["cell_type"]: row for row in config}
    rows: list[dict[str, object]] = []
    for model in MODEL_ORDER:
        subset = annotated[annotated["CellOrTissueType"].eq(model)]
        ddr = subset[subset["GOMatchMode"].ne("unmatched")]
        total = len(subset)
        expected = int(config_by_model[model]["reference_protein_count_main"])
        if total != expected:
            raise AssertionError(f"{model}: extracted {total} proteins/features, expected {expected}")
        rows.append(
            {
                "CellOrTissueType": model,
                "ReferencePXD": config_by_model[model]["reference_pxd"],
                "ReferenceSubset": config_by_model[model]["reference_subset"],
                "ReferenceProteinCount": total,
                "ReferenceDdrProteinCount": len(ddr),
                "ReferenceDdrFraction": len(ddr) / total if total else 0.0,
                "BaseAccessionMatches": int(ddr["GOMatchMode"].eq("BaseAccession").sum()),
                "GeneSymbolFallbackMatches": int(
                    ddr["GOMatchMode"].eq("GeneSymbol_fallback").sum()
                ),
                "ReferenceCountUnit": unique_join(subset["CountUnit"]),
                "ProteinCountDetail": config_by_model[model]["reference_protein_count_detail"],
            }
        )
    return pd.DataFrame(rows)


def build_control_information(
    statistics_frame: pd.DataFrame, config: list[dict[str, object]]
) -> pd.DataFrame:
    statistics_by_model = statistics_frame.set_index("CellOrTissueType").to_dict("index")
    rows = []
    for item in config:
        model = str(item["cell_type"])
        stats = statistics_by_model[model]
        rows.append(
            {
                "细胞或组织": model,
                "参考PXD": item["reference_pxd"],
                "年份": int(item["reference_year"]),
                "匹配类型": item["match_type"],
                "适用等级": item["suitability_grade"],
                "选用样本或子集": item["reference_subset"],
                "正常或基线条件": item["baseline_definition"],
                "样本或采集数": item["reference_sample_count"],
                "对照蛋白数": int(stats["ReferenceProteinCount"]),
                "DDR蛋白数": int(stats["ReferenceDdrProteinCount"]),
                "DDR占比": float(stats["ReferenceDdrFraction"]),
                "采集方式或仪器": item["acquisition"],
                "检索与定量": item["search_quantification"],
                "是否PTM富集": item["ptm_enrichment"],
                "选择理由": item["selection_rationale"],
                "主要限制": item["main_caveat"],
                "蛋白数计数口径": item["reference_protein_count_detail"],
                "实际分析来源文件": item["reference_protein_count_source"],
                "仓库与分析完整度": item["repository_completeness"],
                "可作统计差异": item["use_for_statistical_differential"],
                "数据集URL": item["dataset_url"],
                "论文URL": item["paper_url"],
                "处理数据URL": item["processed_data_url"],
            }
        )
    return pd.DataFrame(rows)


def build_tall_sensitivity(
    primary: pd.DataFrame, jurkat: pd.DataFrame
) -> pd.DataFrame:
    primary_set = set(primary["BaseAccession"])
    jurkat_set = set(jurkat["BaseAccession"])
    go_accessions = set(
        pd.concat([primary, jurkat], ignore_index=True)
        .loc[lambda frame: frame["GOMatchMode"].eq("BaseAccession"), "BaseAccession"]
    )
    rows = []
    for label, proteins in [
        ("TALL-1_primary", primary_set),
        ("Jurkat_sensitivity", jurkat_set),
        ("TALL-1_Jurkat_union", primary_set | jurkat_set),
        ("TALL-1_Jurkat_intersection", primary_set & jurkat_set),
    ]:
        ddr = proteins & go_accessions
        rows.append(
            {
                "TAllReferenceSet": label,
                "ProteinCount": len(proteins),
                "DdrProteinCount": len(ddr),
                "DdrFraction": len(ddr) / len(proteins) if proteins else 0.0,
            }
        )
    return pd.DataFrame(rows)


def plot_comparison(frame: pd.DataFrame) -> None:
    plot = frame.set_index("CellOrTissueType").loc[MODEL_ORDER].reset_index()
    display_names = {
        "Human hippocampus": "人海马组织",
        "HEK293T": "HEK293T",
        "HK-2": "HK-2",
        "MCF10A": "MCF10A",
        "MCF7": "MCF7",
        "HCT116": "HCT116",
        "T-ALL": "T-ALL（TALL-1参照）",
        "MDA-MB-468": "MDA-MB-468",
        "T-47D": "T-47D",
        "RKO": "RKO",
    }
    plot["DisplayName"] = plot["CellOrTissueType"].map(display_names)
    y = np.arange(len(plot))
    height = 0.32

    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": [
                "Arial Unicode MS",
                "Hiragino Sans GB",
                "PingFang SC",
                "DejaVu Sans",
            ],
            "axes.unicode_minus": False,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )
    fig, ax = plt.subplots(figsize=(11.4, 7.2))
    reference_bars = ax.barh(
        y - height / 2,
        plot["ReferenceDdrFraction"] * 100,
        height,
        label="常规全蛋白组参照",
        color="#4E79A7",
        edgecolor="white",
        linewidth=0.5,
    )
    kla_bars = ax.barh(
        y + height / 2,
        plot["KlaGoDdrFraction"] * 100,
        height,
        label="乳酸化蛋白组（Kla）",
        color="#F28E2B",
        edgecolor="white",
        linewidth=0.5,
    )
    ax.set_yticks(y, plot["DisplayName"])
    ax.invert_yaxis()
    ax.set_xlabel("DDR 注释蛋白占比（%）", fontsize=11)
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines["left"].set_color("#8A8A8A")
    ax.spines["bottom"].set_color("#8A8A8A")
    ax.grid(axis="x", color="#D9DDE3", linewidth=0.7, alpha=0.9)
    ax.set_axisbelow(True)
    ax.tick_params(axis="both", labelsize=10, colors="#30343B")
    ax.legend(
        frameon=False,
        loc="lower right",
        fontsize=10,
        ncols=2,
        bbox_to_anchor=(1, 1.005),
        borderaxespad=0,
    )

    maximum = max(
        plot["ReferenceDdrFraction"].max(),
        plot["KlaGoDdrFraction"].max(),
    ) * 100
    label_offset = 0.12
    for bar, count, total, fraction in zip(
        reference_bars,
        plot["ReferenceDdrProteinCount"],
        plot["ReferenceProteinCount"],
        plot["ReferenceDdrFraction"],
    ):
        ax.text(
            bar.get_width() + label_offset,
            bar.get_y() + bar.get_height() / 2,
            f"{int(count):,}/{int(total):,}（{fraction:.1%}）",
            va="center",
            ha="left",
            fontsize=8.5,
            color="#2F3B49",
        )
    for bar, count, total, fraction in zip(
        kla_bars,
        plot["KlaGoDdrProteins"],
        plot["TotalKlaProteins"],
        plot["KlaGoDdrFraction"],
    ):
        ax.text(
            bar.get_width() + label_offset,
            bar.get_y() + bar.get_height() / 2,
            f"{int(count):,}/{int(total):,}（{fraction:.1%}）",
            va="center",
            ha="left",
            fontsize=8.5,
            color="#4A3420",
        )
    ax.set_xlim(0, maximum + 4.2)
    ax.set_xticks(np.arange(0, 17, 2))
    ax.margins(y=0.04)

    fig.suptitle(
        "各细胞系与组织中 DDR 蛋白占比",
        x=0.12,
        y=0.975,
        ha="left",
        fontsize=17,
        fontweight="bold",
        color="#20252B",
    )
    fig.text(
        0.12,
        0.925,
        "常规全蛋白组参照与乳酸化蛋白组比较；柱端为 DDR 蛋白数/总蛋白数（占比）",
        ha="left",
        fontsize=10.5,
        color="#525A64",
    )
    fig.text(
        0.12,
        0.018,
        (
            "注：仅纳入同时具有完整 Kla 蛋白集合、DDR 交集和可计数常规蛋白组的 10 类模型。"
            "T-ALL 使用 TALL-1 作为主参照；海马参照为 74 名神经学正常 CA1 供者。"
            "该图为描述性集合比较，不代表蛋白丰度或统计学显著性。"
        ),
        fontsize=8.2,
        color="#5C626A",
        ha="left",
    )
    fig.tight_layout(rect=(0.095, 0.055, 0.995, 0.91))
    FIGURES.mkdir(parents=True, exist_ok=True)
    for stem in (
        "cell_type_kla_vs_reference_ddr_fraction",
        "cell_type_kla_vs_reference_ddr_fraction_v2",
    ):
        fig.savefig(
            FIGURES / f"{stem}.png",
            dpi=350,
            bbox_inches="tight",
            facecolor="white",
        )
        fig.savefig(
            FIGURES / f"{stem}.pdf",
            bbox_inches="tight",
            facecolor="white",
        )
    plt.close(fig)


def build_source_manifest() -> pd.DataFrame:
    sources = [
        DATA_ROOT
        / "PXD030304/search_results/ProCan-DepMapSanger_protein_matrix_6692_averaged.txt",
        DATA_ROOT
        / "PXD030304/search_results/ProCan-DepMapSanger_mapping_file_replicates.txt",
        DATA_ROOT
        / "PXD030304/search_results/ProCan-DepMapSanger_peptide_counts_per_protein_per_sample.txt",
        DATA_ROOT / "PXD043880/supplementary/13024_2023_650_MOESM1_ESM.xlsx",
        DATA_ROOT / "PXD072220/search_results/HK-2_Spectronaut-report_PG_Quantity.txt",
        DATA_ROOT / "PXD072220/search_results/HK-2_Spectronaut-report_Summary.txt",
        DATA_ROOT / "PXD002400/search_results/msms.zip",
        GO_PATH,
    ]
    return pd.DataFrame(
        [
            {
                "Path": str(path.relative_to(PROJECT_ROOT)),
                "SizeBytes": path.stat().st_size,
                "SHA256": sha256_file(path),
            }
            for path in sources
        ]
    )


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    INTERMEDIATE.mkdir(parents=True, exist_ok=True)
    config = json.loads(CONFIG_PATH.read_text(encoding="utf-8"))

    pxd030304, jurkat_raw = load_pxd030304()
    records = [
        *load_pxd043880(),
        *pxd030304,
        *load_pxd072220(),
        *load_pxd002400(),
    ]
    reference = pd.DataFrame(records).sort_values(
        ["CellOrTissueType", "ReferenceProteinKey"], kind="stable"
    )
    annotated = annotate_go(reference)
    jurkat = annotate_go(jurkat_raw)

    statistics_frame = build_reference_statistics(annotated, config)
    control_information = build_control_information(statistics_frame, config)
    comparison = pd.DataFrame()
    if KLA_STATS_PATH.exists():
        kla = pd.read_csv(KLA_STATS_PATH, encoding="utf-8-sig")
        comparison = statistics_frame.merge(
            kla, on="CellOrTissueType", how="left", validate="one_to_one"
        )
        comparison["DdrFractionPercentagePointDifference"] = (
            comparison["KlaGoDdrFraction"] - comparison["ReferenceDdrFraction"]
        )
        comparison["KlaToReferenceDdrFractionRatio"] = (
            comparison["KlaGoDdrFraction"] / comparison["ReferenceDdrFraction"]
        )

    matched = annotated[annotated["GOMatchMode"].ne("unmatched")].copy()
    unmatched = annotated[annotated["GOMatchMode"].eq("unmatched")].copy()
    tall_primary = annotated[annotated["CellOrTissueType"].eq("T-ALL")]
    tall_sensitivity = build_tall_sensitivity(tall_primary, jurkat)

    write_csv(annotated, TABLES / "reference_proteome_all_proteins.csv")
    write_csv(matched, TABLES / "reference_proteome_ddr_proteins.csv")
    write_csv(unmatched, TABLES / "reference_proteome_go_unmatched.csv")
    write_csv(statistics_frame, TABLES / "cell_type_reference_ddr_statistics.csv")
    write_csv(
        control_information,
        TABLES / "cell_type_reference_control_information.csv",
    )
    if not comparison.empty:
        write_csv(comparison, TABLES / "cell_type_kla_vs_reference_ddr_statistics.csv")
    write_csv(tall_sensitivity, TABLES / "tall104_surrogate_ddr_sensitivity.csv")
    write_csv(build_source_manifest(), TABLES / "reference_proteome_source_manifest.csv")

    for model in MODEL_ORDER:
        slug = re.sub(r"[^a-z0-9]+", "_", model.casefold()).strip("_")
        subset = annotated[annotated["CellOrTissueType"].eq(model)]
        write_csv(subset, INTERMEDIATE / f"{slug}_reference_proteins.csv")
        write_csv(
            subset[subset["GOMatchMode"].ne("unmatched")],
            INTERMEDIATE / f"{slug}_reference_ddr_proteins.csv",
        )

    if not comparison.empty:
        plot_comparison(comparison)

    if len(statistics_frame) != 10 or statistics_frame["ReferenceDdrFraction"].isna().any():
        raise AssertionError("Reference DDR statistics are incomplete")
    if not (
        statistics_frame["ReferenceDdrProteinCount"]
        <= statistics_frame["ReferenceProteinCount"]
    ).all():
        raise AssertionError("DDR counts exceed reference protein counts")

    if not comparison.empty:
        if len(comparison) != len(statistics_frame) or comparison["ReferenceDdrFraction"].isna().any():
            raise AssertionError("Kla-to-reference DDR comparison is incomplete")
        print(
            comparison[
                [
                    "CellOrTissueType",
                    "ReferenceProteinCount",
                    "ReferenceDdrProteinCount",
                    "ReferenceDdrFraction",
                    "KlaGoDdrFraction",
                    "KlaToReferenceDdrFractionRatio",
                ]
            ].to_string(index=False)
        )


if __name__ == "__main__":
    main()
