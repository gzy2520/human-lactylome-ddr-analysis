#!/usr/bin/env python3
"""Read-only audit of the frozen Kla publication tables.

This script only reads the frozen table package and writes audit artefacts under
the current workspace. It does not modify the publication files.
"""

from __future__ import annotations

import hashlib
import re
from pathlib import Path

import numpy as np
import openpyxl
import pandas as pd


WORKSPACE = Path(__file__).resolve().parents[2]
TABLES = Path("/Users/gzy2520/Desktop/Research/kla/results/final_figures_and_tables/Tables")
FINAL = TABLES.parent
SOURCE = Path("/Users/gzy2520/Desktop/Research/kla/data/publication_input/human_ddr_go_annotations.tsv")
FIGURE_AUDIT = WORKSPACE / "audits/20260905_final_result"
OUT = Path(__file__).resolve().parent
OUT.mkdir(parents=True, exist_ok=True)


def clean(v):
    if pd.isna(v):
        return ""
    return str(v).strip()


def save(name: str, rows: list[dict]) -> pd.DataFrame:
    frame = pd.DataFrame(rows)
    frame.to_csv(OUT / name, index=False)
    return frame


def xlsx_sheets(path: Path) -> dict[str, pd.DataFrame]:
    return {sheet: pd.read_excel(path, sheet_name=sheet) for sheet in pd.ExcelFile(path).sheet_names}


def split_accessions(value) -> list[str]:
    if pd.isna(value):
        return []
    return [x.strip() for x in re.split(r"[;,]", str(value)) if x.strip() and x.strip().lower() not in {"nan", "na", "none"}]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    # Inventory, missingness, duplicate rows, and formula presence.
    inventory = []
    for path in sorted(TABLES.glob("*.xlsx")):
        workbook = openpyxl.load_workbook(path, read_only=True, data_only=False)
        sheet_frames = xlsx_sheets(path)
        for sheet, frame in sheet_frames.items():
            ws = workbook[sheet]
            formulas = sum(
                1
                for row in ws.iter_rows()
                for cell in row
                if isinstance(cell.value, str) and cell.value.startswith("=")
            )
            inventory.append(
                {
                    "Workbook": path.name,
                    "Sheet": sheet,
                    "Rows": len(frame),
                    "Columns": len(frame.columns),
                    "NA_cells": int(frame.isna().sum().sum()),
                    "Exact_duplicate_rows": int(frame.duplicated().sum()),
                    "Formula_cells": formulas,
                }
            )
        workbook.close()
    save("workbook_inventory.csv", inventory)

    # Table 14: stable accession key, categorical score states, and score arithmetic.
    weights = {"BER": 1, "NER": 2, "MMR": 3, "FA": 4, "HR": 5, "AEJ": 6, "NHEJ": 7}
    t14_csv = pd.read_csv(TABLES / "Table_14_new_ESCC_Kla_DDR_proteins.csv")
    t14_xlsx = pd.read_excel(TABLES / "Table_14_new_ESCC_Kla_DDR_proteins.xlsx")
    t14_rows = []
    for label, frame in [("csv", t14_csv), ("xlsx", t14_xlsx)]:
        score = sum(weight * frame[col] for col, weight in weights.items())
        t14_rows.append(
            {
                "Check": f"Table14_{label}",
                "Rows": len(frame),
                "Duplicate_BaseAccession": int(frame["BaseAccession"].duplicated().sum()),
                "Invalid_state_cells": int((~frame[list(weights)].isin([-1, 0, 1])).sum().sum()),
                "SignedScore_mismatches": int((score != frame["SignedScore"]).sum()),
                "Status": "PASS" if len(frame) == 14 and frame["BaseAccession"].is_unique and (score == frame["SignedScore"]).all() else "FAIL",
            }
        )
    t14_rows.append(
        {
            "Check": "Table14_csv_xlsx_identity",
            "Rows": len(t14_csv),
            "Duplicate_BaseAccession": "",
            "Invalid_state_cells": "",
            "SignedScore_mismatches": "",
            "Status": "PASS" if t14_csv.fillna("").astype(str).equals(t14_xlsx.fillna("").astype(str)) else "FAIL",
        }
    )
    save("table14_checks.csv", t14_rows)

    # Table S1 group and membership reconciliation.
    s1 = xlsx_sheets(TABLES / "Table_S1_Kla_Data.xlsx")
    g = s1["Group_Summary"]
    k = s1["Kla_Protein_Membership"]
    kd = s1["Kla_DDR_Membership"]
    group_keys = ["PXD", "SampleGroup"]
    kla_counts = k.groupby(group_keys).agg(KlaProteinCount_membership=("BaseAccession", "nunique")).reset_index()
    ddr_counts = kd.groupby(group_keys).agg(KlaDdrProteinCount_membership=("BaseAccession", "nunique")).reset_index()
    reconciled = g.merge(kla_counts, on=group_keys, how="outer").merge(ddr_counts, on=group_keys, how="outer")
    reconciled["KlaCount_diff"] = reconciled["KlaProteinCount"] - reconciled["KlaProteinCount_membership"]
    reconciled["KlaDDRCount_diff"] = reconciled["KlaDdrProteinCount"] - reconciled["KlaDdrProteinCount_membership"]
    reconciled.to_csv(OUT / "s1_membership_reconciliation.csv", index=False)
    s1_checks = [
        {
            "Check": "S1_group_registry",
            "Observed": len(g),
            "Expected": 31,
            "Status": "PASS" if len(g) == 31 and g[group_keys].drop_duplicates().shape[0] == 31 else "FAIL",
            "Note": "31 unique (PXD, SampleGroup); category counts: " + "; ".join(f"{k}={v}" for k, v in g["Category"].value_counts().sort_index().items()),
        },
        {
            "Check": "S1_membership_counts",
            "Observed": int(((reconciled["KlaCount_diff"] != 0) | (reconciled["KlaDDRCount_diff"] != 0)).sum()),
            "Expected": 0,
            "Status": "PASS" if ((reconciled["KlaCount_diff"] == 0) & (reconciled["KlaDDRCount_diff"] == 0)).all() else "FAIL",
            "Note": "Group_Summary counts equal unique BaseAccession counts in membership sheets",
        },
        {
            "Check": "S1_membership_integrity",
            "Observed": int(k.duplicated(group_keys + ["BaseAccession"]).sum() + kd.duplicated(group_keys + ["BaseAccession"]).sum() + (kd["IsDdr"] == False).sum()),
            "Expected": 0,
            "Status": "PASS" if k.duplicated(group_keys + ["BaseAccession"]).sum() == 0 and kd.duplicated(group_keys + ["BaseAccession"]).sum() == 0 and (kd["IsDdr"] == False).sum() == 0 else "FAIL",
            "Note": "No duplicate group/accession keys; Kla_DDR_Membership is all IsDdr=TRUE",
        },
        {
            "Check": "S1_paired_flag",
            "Observed": int((g["PairedAnalysisIncluded"] != True).sum()),
            "Expected": 0,
            "Status": "PASS" if (g["PairedAnalysisIncluded"] == True).all() else "WARN",
            "Note": "All 31 rows are marked paired; this is a metadata flag, not proof of independent biological pairing",
        },
        {
            "Check": "S1_pair_reuse",
            "Observed": f"{g['PXD'].nunique()} Kla PXD; {g['ReferencePXD'].nunique()} reference PXD",
            "Expected": "31 independent pairs",
            "Status": "WARN" if g["ReferencePXD"].nunique() < len(g) else "PASS",
            "Note": "These are 31 group rows, not 31 independent dataset pairs; PXD030304 is reused for 8 groups and several other references are reused",
        },
        {
            "Check": "S1_roworder_contiguity",
            "Observed": ",".join(map(str, sorted(set(range(1, len(g) + 1)) - set(g["RowOrder"].astype(int))))),
            "Expected": "none",
            "Status": "WARN" if set(range(1, len(g) + 1)) - set(g["RowOrder"].astype(int)) else "PASS",
            "Note": "RowOrder preserves legacy category ordering and is non-contiguous; document if used as a display index",
        },
    ]
    save("s1_checks.csv", s1_checks)

    # Table S2 mapping and summary reconciliation.
    s2 = xlsx_sheets(TABLES / "Table_S2_Reference_Data.xlsx")
    rg = s2["Reference_Group_Summary"]
    r = s2["Reference_Protein_Membership"]
    rd = s2["Reference_DDR_Membership"]
    rg_merge = g[group_keys + ["ReferencePXD", "ReferenceProteinCount", "ReferenceDdrProteinCount", "ReferenceDdrFraction", "MatchMode"]].merge(
        rg[group_keys + ["ReferencePXD", "ReferenceProteinCount", "ReferenceDdrProteinCount", "ReferenceDdrFraction", "MatchMode"]],
        on=group_keys,
        suffixes=("_S1", "_S2"),
        how="outer",
    )
    comparable = ["ReferencePXD", "ReferenceProteinCount", "ReferenceDdrProteinCount", "ReferenceDdrFraction", "MatchMode"]
    for col in comparable:
        a, b = f"{col}_S1", f"{col}_S2"
        rg_merge[f"{col}_mismatch"] = rg_merge[a].fillna("<NA>").astype(str) != rg_merge[b].fillna("<NA>").astype(str)
    rg_merge.to_csv(OUT / "s1_s2_summary_reconciliation.csv", index=False)

    s2_rows = []
    for key, frame in r.groupby(group_keys):
        mapped = {acc for value in frame["MappedBaseAccessions"] for acc in split_accessions(value)}
        mapped_ddr = {acc for value in frame.loc[frame["IsDdr"] == True, "MappedBaseAccessions"] for acc in split_accessions(value)}
        summary = rg.set_index(group_keys).loc[key]
        unmapped = frame["MappedBaseAccessions"].isna() | frame["MappedBaseAccessions"].astype(str).str.strip().isin(["", "nan"])
        s2_rows.append(
            {
                "PXD": key[0],
                "SampleGroup": key[1],
                "SourceRows": len(frame),
                "UniqueMappedBaseAccessions": len(mapped),
                "SummaryReferenceProteinCount": int(summary["ReferenceProteinCount"]),
                "UniqueMappedDDRBaseAccessions": len(mapped_ddr),
                "SummaryReferenceDdrProteinCount": int(summary["ReferenceDdrProteinCount"]),
                "UnmappedRows": int(unmapped.sum()),
                "UnmappedDDRRows": int((unmapped & (frame["IsDdr"] == True)).sum()),
                "Status": "WARN" if unmapped.any() else "PASS",
            }
        )
    s2_frame = pd.DataFrame(s2_rows)
    s2_frame["Mapped_vs_summary_diff"] = s2_frame["UniqueMappedBaseAccessions"] - s2_frame["SummaryReferenceProteinCount"]
    s2_frame["MappedDDR_vs_summary_diff"] = s2_frame["UniqueMappedDDRBaseAccessions"] - s2_frame["SummaryReferenceDdrProteinCount"]
    s2_frame.to_csv(OUT / "s2_mapping_reconciliation.csv", index=False)
    s2_checks = [
        {
            "Check": "S1_S2_summary_fields",
            "Observed": int(sum(rg_merge[f"{c}_mismatch"].sum() for c in comparable)),
            "Expected": 0,
            "Status": "PASS" if int(sum(rg_merge[f"{c}_mismatch"].sum() for c in comparable)) == 0 else "FAIL",
            "Note": "Reference fields agree between Table S1 Group_Summary and Table S2 Reference_Group_Summary",
        },
        {
            "Check": "S2_membership_keys",
            "Observed": int(r.duplicated(["PXD", "SampleGroup", "ReferencePXD", "SourceProteinID"]).sum() + rd.duplicated(["PXD", "SampleGroup", "ReferencePXD", "SourceProteinID"]).sum() + (rd["IsDdr"] == False).sum()),
            "Expected": 0,
            "Status": "PASS" if r.duplicated(["PXD", "SampleGroup", "ReferencePXD", "SourceProteinID"]).sum() == 0 and rd.duplicated(["PXD", "SampleGroup", "ReferencePXD", "SourceProteinID"]).sum() == 0 and (rd["IsDdr"] == False).sum() == 0 else "FAIL",
            "Note": "No duplicate source-protein keys; Reference_DDR_Membership is all IsDdr=TRUE",
        },
        {
            "Check": "S2_unmapped_reference_ids",
            "Observed": int(s2_frame["UnmappedRows"].sum()),
            "Expected": 0,
            "Status": "WARN" if s2_frame["UnmappedRows"].sum() else "PASS",
            "Note": "Unmapped ENSEMBLPROT rows remain in normal human lung (595) and normal pregnancy placenta (546); summary counts include them, while BaseAccession-based membership cannot",
        },
    ]
    save("s2_checks.csv", s2_checks)

    # Table S3 exact source reconciliation.
    s3 = pd.read_excel(TABLES / "Table_S3_Human_DDR_GO_Annotations.xlsx")
    src = pd.read_csv(SOURCE, sep="\t")
    s3_mismatches = 0
    for col in src.columns:
        s3_mismatches += int((s3[col].fillna("").astype(str).str.strip() != src[col].fillna("").astype(str).str.strip()).sum())
    save(
        "s3_source_reconciliation.csv",
        [
            {
                "TableRows": len(s3),
                "SourceRows": len(src),
                "TableColumns": len(s3.columns),
                "SourceColumns": len(src.columns),
                "CellMismatchesAfter_NA_Whitespace_Normalization": s3_mismatches,
                "TableUniqueGeneProductIDs": s3["GENE PRODUCT ID"].nunique(),
                "SourceUniqueGeneProductIDs": src["GENE PRODUCT ID"].nunique(),
                "ExactDuplicateRows": int(s3.duplicated().sum()),
                "Status": "PASS" if s3_mismatches == 0 and len(s3) == len(src) and list(s3.columns) == list(src.columns) else "FAIL",
                "Note": "The workbook exactly matches the current 5,111-row TSV; the package README's 6,707-row description is stale",
            }
        ],
    )

    # Table S4: stable accession membership, region labels, set counts.
    s4_file = TABLES / "Table_S4_Venn_Membership.xlsx"
    flags = ["In_normal_tissue", "In_cancer_tissue", "In_cancer_cells", "In_normal_cells"]
    s4_checks = []
    set_counts = pd.read_excel(s4_file, sheet_name="Set_Counts")
    region_counts = pd.read_excel(s4_file, sheet_name="Region_Counts")

    def expected_region(row):
        present = [flag.replace("In_", "") for flag in flags if bool(row[flag])]
        if len(present) == 4:
            return "all_four"
        if not present:
            return "none"
        return "_and_".join(present) + "_only"

    for sheet, analysis in [("AllKla_Members", "AllKla"), ("KlaDDR_Members", "KlaDDR"), ("Reference_Members", "Reference"), ("ReferenceDDR_Members", "ReferenceDDR")]:
        frame = pd.read_excel(s4_file, sheet_name=sheet)
        region_mismatch = int((frame.apply(expected_region, axis=1) != frame["Region"]).sum())
        set_mismatch = 0
        for flag in flags:
            expected = int(frame[flag].sum())
            listed = set_counts.loc[(set_counts["Analysis"] == analysis) & (set_counts["Category"] == flag.replace("In_", "")), "ProteinCount"]
            set_mismatch += int(len(listed) != 1 or expected != int(listed.iloc[0]))
        listed_region_total = int(region_counts.loc[region_counts["Analysis"] == analysis, "ProteinCount"].sum())
        s4_checks.append(
            {
                "Sheet": sheet,
                "Rows": len(frame),
                "Duplicate_BaseAccession": int(frame["BaseAccession"].duplicated().sum()),
                "Invalid_flag_cells": int((~frame[flags].isin([True, False])).sum().sum()),
                "Region_label_mismatches": region_mismatch,
                "Set_Count_mismatches": set_mismatch,
                "Region_total_mismatch": int(listed_region_total != len(frame)),
                "Status": "PASS" if frame["BaseAccession"].is_unique and region_mismatch == 0 and set_mismatch == 0 and listed_region_total == len(frame) else "FAIL",
            }
        )
    save("s4_checks.csv", s4_checks)

    # Table S5: signed score arithmetic and consistency with S4 Kla-DDR sets.
    s5_file = TABLES / "Table_S5_Pathway_Protein_Ranking.xlsx"
    s5_checks = []
    for sheet, category in [("NonTumorTissues", "normal_tissue"), ("TumorTissues", "cancer_tissue"), ("CancerCellLines", "cancer_cells"), ("NormalCellLines", "normal_cells")]:
        frame = pd.read_excel(s5_file, sheet_name=sheet)
        score = sum(weight * frame[col] for col, weight in weights.items())
        s4 = pd.read_excel(s4_file, sheet_name="KlaDDR_Members")
        expected_ids = set(s4.loc[s4[f"In_{category}"], "BaseAccession"])
        observed_ids = set(frame["BaseAccession"])
        s5_checks.append(
            {
                "Sheet": sheet,
                "Rows": len(frame),
                "Duplicate_BaseAccession": int(frame["BaseAccession"].duplicated().sum()),
                "Invalid_state_cells": int((~frame[list(weights)].isin([-1, 0, 1])).sum().sum()),
                "SignedScore_mismatches": int((score != frame["SignedScore"]).sum()),
                "Missing_from_S4": len(observed_ids - expected_ids),
                "Extra_vs_S4": len(expected_ids - observed_ids),
                "Status": "PASS" if frame["BaseAccession"].is_unique and (score == frame["SignedScore"]).all() and observed_ids == expected_ids else "FAIL",
            }
        )
    save("s5_checks.csv", s5_checks)

    # Table S6: regulator record duplication and reference provenance.
    s6 = xlsx_sheets(TABLES / "Table_S6_Lactylation_Regulators.xlsx")
    ann = s6["Regulator_Annotations"]
    mapping = s6["Regulator_ID_Mapping"]
    duplicate_record = ann[ann.duplicated(["Role", "GeneSymbol", "BaseAccession"], keep=False)].copy()
    missing_refs = ann[ann["References"].isna()].copy()
    # A newline-delimited list is readable; a second URL glued directly to the
    # previous DOI is not.  Detect the latter without rejecting valid lists.
    malformed_multi_url = ann["References"].fillna("").astype(str).map(
        lambda x: x.count("https://") > 1 and bool(re.search(r"\Shttps://", x))
    )
    save("s6_duplicate_records.csv", duplicate_record)
    save("s6_missing_references.csv", missing_refs)
    save("s6_checks.csv", [
        {
            "Check": "S6_annotation_rows",
            "Observed": len(ann),
            "Expected_documented_unique_records": 49,
            "Status": "WARN" if len(ann) != 49 or ann[["Role", "GeneSymbol", "BaseAccession"]].drop_duplicates().shape[0] != 49 else "PASS",
            "Note": "50 rows, 49 unique (Role, GeneSymbol, BaseAccession) records, 48 unique BaseAccessions; HDAC5/Q9UQL6 is duplicated in the same role (records 23 and 24), while HDAC8/Q9BY41 has two roles",
        },
        {
            "Check": "S6_mapping_coverage",
            "Observed": len(set(ann["BaseAccession"]) - set(mapping["BaseAccession"])) + len(set(mapping["BaseAccession"]) - set(ann["BaseAccession"])),
            "Expected": 0,
            "Status": "PASS" if set(ann["BaseAccession"]) == set(mapping["BaseAccession"]) else "FAIL",
            "Note": "Regulator_ID_Mapping covers all 48 unique annotation accessions",
        },
        {
            "Check": "S6_missing_references",
            "Observed": len(missing_refs),
            "Expected": 0,
            "Status": "WARN" if len(missing_refs) else "PASS",
            "Note": "12 annotation rows have no reference; this is a provenance gap unless explicitly marked as unavailable",
        },
        {
            "Check": "S6_reference_delimiter",
            "Observed": int(malformed_multi_url.sum()),
            "Expected": 0,
            "Status": "WARN" if malformed_multi_url.any() else "PASS",
            "Note": "At least one cell concatenates multiple DOI URLs without a delimiter (PARK7 record 42); newline-separated multiple URLs are acceptable",
        },
    ])

    # Statistics tables: denominator and observation-level audit.
    stat_dir = TABLES / "statistical_tests"
    stats_rows = []
    f1_omni = pd.read_csv(stat_dir / "figure1_category_omnibus_anova.csv")
    f1_one = pd.read_csv(stat_dir / "figure1_category_one_way_anova.csv")
    pathway = pd.read_csv(stat_dir / "pathway_summary_two_way_anova.csv")
    mean_median = pd.read_csv(stat_dir / "figure1_category_boxplot_mean_median.csv")
    stats_rows.extend([
        {
            "Table": "figure1_category_boxplot_mean_median.csv",
            "Rows": len(mean_median),
            "N_or_NPoint": "not included",
            "NonSample_issue": "Cannot audit denominator from this table alone",
            "Status": "WARN",
            "Note": "Add N, sample/aggregate/pool/condition/model classification, and source version to make the summary auditable",
        },
        {
            "Table": "figure1_category_omnibus_anova.csv",
            "Rows": len(f1_omni),
            "N_or_NPoint": ",".join(map(str, f1_omni["N"].astype(int).unique())),
            "NonSample_issue": "N=310 includes 38 non-sample observations; sample-only expanded input is 272",
            "Status": "FAIL",
            "Note": "Comparison says controlling for modality, but the default sequential two-way ANOVA category term is order-dependent; this is a method/label issue",
        },
        {
            "Table": "figure1_category_one_way_anova.csv",
            "Rows": len(f1_one),
            "N_or_NPoint": ",".join(map(str, f1_one["N"].astype(int))),
            "NonSample_issue": "Category totals include 38 non-sample observations; sample-only totals would be 91, 118, 41, 22",
            "Status": "FAIL",
            "Note": "Dataset-factor tests inherit the mixed observation unit; unused NA columns are harmless but should be documented or removed",
        },
        {
            "Table": "pathway_summary_two_way_anova.csv",
            "Rows": len(pathway),
            "N_or_NPoint": f"N={int(pathway['N'].iloc[0])}; NPoint={int(pathway['NPoint'].iloc[0])}",
            "NonSample_issue": "NPoint=98 includes 26 non-sample source observations; sample-only expanded input is 72",
            "Status": "FAIL",
            "Note": "Pro/Inh are two measurements from one source observation but are fit as independent rows; the current p/q values are not valid paired inference",
        },
    ])
    save("statistics_table_checks.csv", stats_rows)

    # Cross-artifact mismatch: formal figures in the frozen package are old while tables/statistics are expanded.
    figure_rows = []
    if (FIGURE_AUDIT / "per_figure_audit.csv").exists():
        figure_audit = pd.read_csv(FIGURE_AUDIT / "per_figure_audit.csv")
        selected = figure_audit[figure_audit["File"].str.contains(r"Figure_1a|Figure_2[c-i]_DDR_pathway_summary", regex=True, na=False)]
        for _, row in selected.iterrows():
            figure_rows.append(
                {
                    "FigureFile": row["File"],
                    "FigureScope": row["Scope"],
                    "FigureUnit": row["Unit"],
                    "FigureAssessment": row["Assessment"],
                    "TableEvidence": "Figure 1 omnibus N=310; pathway NPoint=98",
                    "Status": "FAIL" if row["Scope"] == "baseline30" else "WARN",
                    "Note": row["Evidence"],
                }
            )
    save("cross_artifact_mismatch.csv", figure_rows)

    # Hashes provide a simple read-only preservation trace for the frozen table files.
    hash_rows = []
    for path in sorted(TABLES.rglob("*")):
        if path.is_file() and path.name != ".DS_Store":
            hash_rows.append({"File": str(path.relative_to(TABLES)), "SHA256": sha256(path), "Bytes": path.stat().st_size})
    save("table_hashes.csv", hash_rows)

    print(f"Wrote table audit artefacts to {OUT}")


if __name__ == "__main__":
    main()
