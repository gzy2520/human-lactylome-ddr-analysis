#!/usr/bin/env python3
"""
Generate GraphPad Prism-compatible statistical data tables and master Excel workbook
for all figures in the Human Kla x DDR project requiring statistical inference.

Data baseline: 2026-10-11 all-groups audit (31 groups, ESCC included, 307 source units, 28 RNA references).
Styling: Clean monochrome professional formatting, Arial typography, auto column widths.
"""

import os
import gzip
import numpy as np
import pandas as pd
import scipy.stats as stats
import openpyxl
from openpyxl.styles import Font, Alignment, Border, Side
from openpyxl.utils import get_column_letter

# Set reproducible random seed as required by AGENTS.md
np.random.seed(25)

# Paths
ROOT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
WT_DIR = os.path.join(ROOT_DIR, ".worktrees", "kla31-rnaseq-20260912")
BASE_AUDITED = os.path.join(WT_DIR, "deliverables", "20261011_lactylation_ddr_all_groups_audited", "tables")
OUTPUT_DIR = os.path.join(ROOT_DIR, "results", "交付文件", "Prism统计数据表_20261015")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Standard definitions
CATEGORY_ORDER = ["normal_tissue", "cancer_tissue", "normal_cells", "cancer_cells"]
CATEGORY_LABELS = {
    "normal_tissue": "non-tumor tissues",
    "cancer_tissue": "tumor tissues",
    "normal_cells": "normal cell lines",
    "cancer_cells": "cancer cell lines"
}

# Typography & styling (Clean monochrome per AGENTS.md)
FONT_TITLE = Font(name="Arial", size=13, bold=True, color="000000")
FONT_SECTION = Font(name="Arial", size=11, bold=True, color="000000")
FONT_HEADER = Font(name="Arial", size=10, bold=True, color="000000")
FONT_BODY = Font(name="Arial", size=10, color="000000")
FONT_NOTE = Font(name="Arial", size=9, italic=True, color="444444")

ALIGN_LEFT = Alignment(horizontal="left", vertical="center")
ALIGN_RIGHT = Alignment(horizontal="right", vertical="center")
ALIGN_CENTER = Alignment(horizontal="center", vertical="center")

BORDER_THIN = Border(
    left=Side(style='thin', color='D0D0D0'),
    right=Side(style='thin', color='D0D0D0'),
    top=Side(style='thin', color='D0D0D0'),
    bottom=Side(style='thin', color='D0D0D0')
)
BORDER_HEADER = Border(
    left=Side(style='thin', color='A0A0A0'),
    right=Side(style='thin', color='A0A0A0'),
    top=Side(style='medium', color='000000'),
    bottom=Side(style='medium', color='000000')
)


def format_sig_stars(p_val):
    if np.isnan(p_val):
        return "N/A"
    if p_val < 0.0001:
        return "****"
    if p_val < 0.001:
        return "***"
    if p_val < 0.01:
        return "**"
    if p_val < 0.05:
        return "*"
    return "ns"


def auto_fit_columns(ws, max_len_cap=55):
    for col in ws.columns:
        col_letter = get_column_letter(col[0].column)
        max_len = 0
        for cell in col:
            val = cell.value
            if val is not None:
                s = str(val)
                lines = s.split("\n")
                line_max = max(len(l) for l in lines)
                if line_max > max_len:
                    max_len = line_max
        ws.column_dimensions[col_letter].width = min(max(max_len + 3, 10), max_len_cap)


print("Loading audited source data...")
df_map = pd.read_csv(os.path.join(BASE_AUDITED, "display_mapping_31.csv"))
df_ddr_units = pd.read_csv(os.path.join(BASE_AUDITED, "ddr_source_unit_values.csv"))
df_ddr_units = df_ddr_units.merge(df_map[['GroupID', 'Category']].drop_duplicates(), on='GroupID', how='left')

# MKI67
mki67_path = os.path.join(WT_DIR, "data", "candidate", "figure1_mki67_ratio_sample_values.csv")
df_mki67 = pd.read_csv(mki67_path)
df_mki67_sample = df_mki67[df_mki67['ObservationType'] == 'sample'].copy()

# G2M
g2m_mat_path = os.path.join(WT_DIR, "outputs", "20260919_repair_release", "core", "figure1b_hallmark_proliferation_28materials.csv")
df_g2m_mat = pd.read_csv(g2m_mat_path)
g2m_samp_path = os.path.join(WT_DIR, "outputs", "20260919_repair_release", "core", "figure1b_hallmark_proliferation_sample_values.csv")
df_g2m_samp = pd.read_csv(g2m_samp_path)

# S8 Pathway
pathway_path = os.path.join(WT_DIR, "data", "candidate", "figure1_pathway_summary_sample_boxplot_values.csv")
df_pathway = pd.read_csv(pathway_path)
df_pathway_sample = df_pathway[df_pathway['ObservationType'] == 'sample'].copy()

# RNA Expression
with gzip.open(os.path.join(BASE_AUDITED, "rna_values_28_references.csv.gz"), "rt") as f:
    df_rna = pd.read_csv(f)
df_rna = df_rna.merge(df_map[['ReferenceKey', 'Category']].drop_duplicates(), on='ReferenceKey', how='left')

# S6A & S6B
s6a_path = os.path.join(WT_DIR, "outputs", "20260919_repair_release", "core", "figure1a_rna_ddr_and_lactylated_28materials.csv")
df_s6a = pd.read_csv(s6a_path)
s6b_path = os.path.join(WT_DIR, "outputs", "20260919_repair_release", "core", "figure1a_fixed_ddr_28materials.csv")
df_s6b = pd.read_csv(s6b_path)

# Initialize collectors
statistical_results = []
tidy_records = []

# ==============================================================================
# 1. FIG 1B: DDR FRACTION (WHOLE PROTEOME & LACTYLOME)
# ==============================================================================
print("Processing Fig 1B...")

# Source units level (307 source units: 211 Whole, 96 Kla)
for ds in ["Whole proteome", "Lactylome (Kla)"]:
    sub = df_ddr_units[df_ddr_units['Dataset'] == ds]
    cat_vals = [sub[sub['Category'] == c]['Fraction'].dropna().values for c in CATEGORY_ORDER]
    f_stat, p_val = stats.f_oneway(*cat_vals)
    kw_stat, kw_pval = stats.kruskal(*cat_vals)
    ns = [len(v) for v in cat_vals]
    means = [np.mean(v) for v in cat_vals]
    sds = [np.std(v, ddof=1) if len(v) > 1 else 0.0 for v in cat_vals]
    meds = [np.median(v) for v in cat_vals]
    iqrs = [stats.iqr(v) for v in cat_vals]

    statistical_results.append({
        "FigurePanel": "Fig 1B",
        "FigureDescription": "Four-category DDR fraction comparison (Source-unit level)",
        "ModalityOrGene": ds,
        "FactorOrDirection": "Category",
        "AnalysisLevel": f"Source-unit level (N={sum(ns)})",
        "TestMethod": "One-way ANOVA",
        "N_non_tumor_tissues": ns[0],
        "N_tumor_tissues": ns[1],
        "N_normal_cell_lines": ns[2],
        "N_cancer_cell_lines": ns[3],
        "Mean_SD_non_tumor": f"{means[0]:.2f}% ± {sds[0]:.2f}%",
        "Mean_SD_tumor": f"{means[1]:.2f}% ± {sds[1]:.2f}%",
        "Mean_SD_normal_cells": f"{means[2]:.2f}% ± {sds[2]:.2f}%",
        "Mean_SD_cancer_cells": f"{means[3]:.2f}% ± {sds[3]:.2f}%",
        "Median_IQR_non_tumor": f"{meds[0]:.2f}% ({iqrs[0]:.2f}%)",
        "Median_IQR_tumor": f"{meds[1]:.2f}% ({iqrs[1]:.2f}%)",
        "Median_IQR_normal_cells": f"{meds[2]:.2f}% ({iqrs[2]:.2f}%)",
        "Median_IQR_cancer_cells": f"{meds[3]:.2f}% ({iqrs[3]:.2f}%)",
        "TestStatistic": round(f_stat, 2),
        "StatisticType": "F",
        "df1": len(CATEGORY_ORDER) - 1,
        "df2": sum(ns) - len(CATEGORY_ORDER),
        "RawPValue": p_val,
        "KruskalStatistic": round(kw_stat, 2),
        "KruskalPValue": kw_pval
    })

# Material level pre-aggregated (31 materials: 9 non-tumor, 3 tumor, 7 normal cells, 12 cancer cells)
mat_ddr = df_ddr_units.groupby(['GroupID', 'Category', 'Dataset'])['Fraction'].mean().reset_index()
for ds in ["Whole proteome", "Lactylome (Kla)"]:
    sub = mat_ddr[mat_ddr['Dataset'] == ds]
    cat_vals = [sub[sub['Category'] == c]['Fraction'].dropna().values for c in CATEGORY_ORDER]
    f_stat, p_val = stats.f_oneway(*cat_vals)
    kw_stat, kw_pval = stats.kruskal(*cat_vals)
    ns = [len(v) for v in cat_vals]
    means = [np.mean(v) for v in cat_vals]
    sds = [np.std(v, ddof=1) if len(v) > 1 else 0.0 for v in cat_vals]
    meds = [np.median(v) for v in cat_vals]
    iqrs = [stats.iqr(v) for v in cat_vals]

    statistical_results.append({
        "FigurePanel": "Fig 1B",
        "FigureDescription": "Four-category DDR fraction comparison (31 materials pre-aggregated)",
        "ModalityOrGene": ds,
        "FactorOrDirection": "Category",
        "AnalysisLevel": f"Material-level (N={sum(ns)})",
        "TestMethod": "One-way ANOVA",
        "N_non_tumor_tissues": ns[0],
        "N_tumor_tissues": ns[1],
        "N_normal_cell_lines": ns[2],
        "N_cancer_cell_lines": ns[3],
        "Mean_SD_non_tumor": f"{means[0]:.2f}% ± {sds[0]:.2f}%",
        "Mean_SD_tumor": f"{means[1]:.2f}% ± {sds[1]:.2f}%",
        "Mean_SD_normal_cells": f"{means[2]:.2f}% ± {sds[2]:.2f}%",
        "Mean_SD_cancer_cells": f"{means[3]:.2f}% ± {sds[3]:.2f}%",
        "Median_IQR_non_tumor": f"{meds[0]:.2f}% ({iqrs[0]:.2f}%)",
        "Median_IQR_tumor": f"{meds[1]:.2f}% ({iqrs[1]:.2f}%)",
        "Median_IQR_normal_cells": f"{meds[2]:.2f}% ({iqrs[2]:.2f}%)",
        "Median_IQR_cancer_cells": f"{meds[3]:.2f}% ({iqrs[3]:.2f}%)",
        "TestStatistic": round(f_stat, 2),
        "StatisticType": "F",
        "df1": len(CATEGORY_ORDER) - 1,
        "df2": sum(ns) - len(CATEGORY_ORDER),
        "RawPValue": p_val,
        "KruskalStatistic": round(kw_stat, 2),
        "KruskalPValue": kw_pval
    })

# Collect Tidy rows for Fig 1B
for _, r in df_ddr_units.iterrows():
    tidy_records.append({
        "FigureID": "Figure 1",
        "PanelID": "Fig 1B",
        "AssayModality": r['Dataset'],
        "GroupID": r['GroupID'],
        "PXD": r['PXD'],
        "SampleID": r['SampleID'],
        "CategoryCode": r['Category'],
        "CategoryLabel": CATEGORY_LABELS.get(r['Category'], r['Category']),
        "MaterialLabel": r['PlotLabel'],
        "FeatureName": "DDR_annotated_protein_fraction",
        "Direction": "N/A",
        "ObservationType": r.get('UnitType', 'source_unit'),
        "Value": r['Fraction'],
        "Unit": "Percentage (%)",
        "SourceFile": "ddr_source_unit_values.csv"
    })

# ==============================================================================
# 2. S3A-C: MKI67 PROLIFERATION RATIOS
# ==============================================================================
print("Processing S3A-C (MKI67)...")
mki67_normalizers = [
    ("H3C1", "S3A", "MKI67 / H3C1 (P46013 / P68431)"),
    ("ACTB", "S3B", "MKI67 / ACTB (P46013 / P60709)"),
    ("TUBB", "S3C", "MKI67 / TUBB (P46013 / P07437)")
]

for denom, panel_id, desc in mki67_normalizers:
    sub = df_mki67_sample[df_mki67_sample['Denominator'] == denom]
    cat_vals = [sub[sub['Category'] == c]['Ratio'].dropna().values for c in CATEGORY_ORDER]
    f_stat, p_val = stats.f_oneway(*cat_vals)
    kw_stat, kw_pval = stats.kruskal(*cat_vals)
    ns = [len(v) for v in cat_vals]
    means = [np.mean(v) if len(v) > 0 else np.nan for v in cat_vals]
    sds = [np.std(v, ddof=1) if len(v) > 1 else 0.0 for v in cat_vals]
    meds = [np.median(v) if len(v) > 0 else np.nan for v in cat_vals]
    iqrs = [stats.iqr(v) if len(v) > 0 else np.nan for v in cat_vals]

    statistical_results.append({
        "FigurePanel": panel_id,
        "FigureDescription": f"Four-category proliferation ratio: {desc}",
        "ModalityOrGene": f"MKI67/{denom}",
        "FactorOrDirection": "Category",
        "AnalysisLevel": f"Sample-level (N={sum(ns)})",
        "TestMethod": "One-way ANOVA",
        "N_non_tumor_tissues": ns[0],
        "N_tumor_tissues": ns[1],
        "N_normal_cell_lines": ns[2],
        "N_cancer_cell_lines": ns[3],
        "Mean_SD_non_tumor": f"{means[0]:.4f} ± {sds[0]:.4f}",
        "Mean_SD_tumor": f"{means[1]:.4f} ± {sds[1]:.4f}",
        "Mean_SD_normal_cells": f"{means[2]:.4f} ± {sds[2]:.4f}",
        "Mean_SD_cancer_cells": f"{means[3]:.4f} ± {sds[3]:.4f}",
        "Median_IQR_non_tumor": f"{meds[0]:.4f} ({iqrs[0]:.4f})",
        "Median_IQR_tumor": f"{meds[1]:.4f} ({iqrs[1]:.4f})",
        "Median_IQR_normal_cells": f"{meds[2]:.4f} ({iqrs[2]:.4f})",
        "Median_IQR_cancer_cells": f"{meds[3]:.4f} ({iqrs[3]:.4f})",
        "TestStatistic": round(f_stat, 2),
        "StatisticType": "F",
        "df1": len(CATEGORY_ORDER) - 1,
        "df2": sum(ns) - len(CATEGORY_ORDER),
        "RawPValue": p_val,
        "KruskalStatistic": round(kw_stat, 2),
        "KruskalPValue": kw_pval
    })

    for _, r in sub.iterrows():
        tidy_records.append({
            "FigureID": "Supplementary Figure S3",
            "PanelID": panel_id,
            "AssayModality": "Whole proteome",
            "GroupID": r.get('DisplayGroup', 'N/A'),
            "PXD": r['PXD'],
            "SampleID": r['SampleID'],
            "CategoryCode": r['Category'],
            "CategoryLabel": CATEGORY_LABELS.get(r['Category'], r['Category']),
            "MaterialLabel": r.get('SampleGroup', 'N/A'),
            "FeatureName": f"MKI67_over_{denom}",
            "Direction": "N/A",
            "ObservationType": "sample",
            "Value": r['Ratio'],
            "Unit": "Intensity ratio",
            "SourceFile": "figure1_mki67_ratio_sample_values.csv"
        })

# ==============================================================================
# 3. S3D: RNA HALLMARK G2M CHECKPOINT
# ==============================================================================
print("Processing S3D (G2M)...")
# 28 materials level
g2m_mat_vals = [df_g2m_mat[df_g2m_mat['Category'] == c]['ProliferationMedian'].dropna().values for c in CATEGORY_ORDER]
f_stat, p_val = stats.f_oneway(*g2m_mat_vals)
kw_stat, kw_pval = stats.kruskal(*g2m_mat_vals)
ns = [len(v) for v in g2m_mat_vals]
means = [np.mean(v) for v in g2m_mat_vals]
sds = [np.std(v, ddof=1) if len(v) > 1 else 0.0 for v in g2m_mat_vals]
meds = [np.median(v) for v in g2m_mat_vals]
iqrs = [stats.iqr(v) for v in g2m_mat_vals]

statistical_results.append({
    "FigurePanel": "S3D",
    "FigureDescription": "Transcriptomic Hallmark G2M checkpoint (28 independent RNA materials)",
    "ModalityOrGene": "Hallmark G2M (196 genes)",
    "FactorOrDirection": "Category",
    "AnalysisLevel": f"Material-level (N={sum(ns)})",
    "TestMethod": "One-way ANOVA",
    "N_non_tumor_tissues": ns[0],
    "N_tumor_tissues": ns[1],
    "N_normal_cell_lines": ns[2],
    "N_cancer_cell_lines": ns[3],
    "Mean_SD_non_tumor": f"{means[0]:.4f} ± {sds[0]:.4f}",
    "Mean_SD_tumor": f"{means[1]:.4f} ± {sds[1]:.4f}",
    "Mean_SD_normal_cells": f"{means[2]:.4f} ± {sds[2]:.4f}",
    "Mean_SD_cancer_cells": f"{means[3]:.4f} ± {sds[3]:.4f}",
    "Median_IQR_non_tumor": f"{meds[0]:.4f} ({iqrs[0]:.4f})",
    "Median_IQR_tumor": f"{meds[1]:.4f} ({iqrs[1]:.4f})",
    "Median_IQR_normal_cells": f"{meds[2]:.4f} ({iqrs[2]:.4f})",
    "Median_IQR_cancer_cells": f"{meds[3]:.4f} ({iqrs[3]:.4f})",
    "TestStatistic": round(f_stat, 2),
    "StatisticType": "F",
    "df1": len(CATEGORY_ORDER) - 1,
    "df2": sum(ns) - len(CATEGORY_ORDER),
    "RawPValue": p_val,
    "KruskalStatistic": round(kw_stat, 2),
    "KruskalPValue": kw_pval
})

for _, r in df_g2m_mat.iterrows():
    tidy_records.append({
        "FigureID": "Supplementary Figure S3",
        "PanelID": "S3D",
        "AssayModality": "Transcriptome (RNA-seq)",
        "GroupID": r['GroupID'],
        "PXD": r['PXD'],
        "SampleID": r['ReferenceKey'],
        "CategoryCode": r['Category'],
        "CategoryLabel": CATEGORY_LABELS.get(r['Category'], r['Category']),
        "MaterialLabel": r['ReferenceKey'],
        "FeatureName": "Hallmark_G2M_checkpoint_median",
        "Direction": "N/A",
        "ObservationType": "material_median",
        "Value": r['ProliferationMedian'],
        "Unit": "log2(qsmooth TPM + 0.5)",
        "SourceFile": "figure1b_hallmark_proliferation_28materials.csv"
    })

# Sample level (1898 samples)
g2m_samp_vals = [df_g2m_samp[df_g2m_samp['Category'] == c]['ProliferationScore'].dropna().values for c in CATEGORY_ORDER]
f_stat, p_val = stats.f_oneway(*g2m_samp_vals)
kw_stat, kw_pval = stats.kruskal(*g2m_samp_vals)
ns = [len(v) for v in g2m_samp_vals]
means = [np.mean(v) for v in g2m_samp_vals]
sds = [np.std(v, ddof=1) if len(v) > 1 else 0.0 for v in g2m_samp_vals]
meds = [np.median(v) for v in g2m_samp_vals]
iqrs = [stats.iqr(v) for v in g2m_samp_vals]

statistical_results.append({
    "FigurePanel": "S3D (Samples)",
    "FigureDescription": "Transcriptomic Hallmark G2M checkpoint (1898 deposited RNA samples)",
    "ModalityOrGene": "Hallmark G2M (196 genes)",
    "FactorOrDirection": "Category",
    "AnalysisLevel": f"Sample-level (N={sum(ns)})",
    "TestMethod": "One-way ANOVA",
    "N_non_tumor_tissues": ns[0],
    "N_tumor_tissues": ns[1],
    "N_normal_cell_lines": ns[2],
    "N_cancer_cell_lines": ns[3],
    "Mean_SD_non_tumor": f"{means[0]:.4f} ± {sds[0]:.4f}",
    "Mean_SD_tumor": f"{means[1]:.4f} ± {sds[1]:.4f}",
    "Mean_SD_normal_cells": f"{means[2]:.4f} ± {sds[2]:.4f}",
    "Mean_SD_cancer_cells": f"{means[3]:.4f} ± {sds[3]:.4f}",
    "Median_IQR_non_tumor": f"{meds[0]:.4f} ({iqrs[0]:.4f})",
    "Median_IQR_tumor": f"{meds[1]:.4f} ({iqrs[1]:.4f})",
    "Median_IQR_normal_cells": f"{meds[2]:.4f} ({iqrs[2]:.4f})",
    "Median_IQR_cancer_cells": f"{meds[3]:.4f} ({iqrs[3]:.4f})",
    "TestStatistic": round(f_stat, 2),
    "StatisticType": "F",
    "df1": len(CATEGORY_ORDER) - 1,
    "df2": sum(ns) - len(CATEGORY_ORDER),
    "RawPValue": p_val,
    "KruskalStatistic": round(kw_stat, 2),
    "KruskalPValue": kw_pval
})

# ==============================================================================
# 4. S8A-G: SEVEN DDR PATHWAYS PRO & INH SUMMARY
# ==============================================================================
print("Processing S8A-G (Pathways)...")
pathways = [
    ("BER", "S8A", "Base excision repair"),
    ("NER", "S8B", "Nucleotide excision repair"),
    ("MMR", "S8C", "Mismatch repair"),
    ("FA", "S8D", "Fanconi anemia"),
    ("HR", "S8E", "Homologous recombination"),
    ("AEJ", "S8F", "Alternative end joining"),
    ("NHEJ", "S8G", "Non-homologous end joining")
]

for pw, panel_id, desc in pathways:
    sub = df_pathway_sample[df_pathway_sample['Pathway'] == pw]
    
    # Promoting (Pro)
    cat_vals_pro = [sub[sub['Category'] == c]['PositiveFraction'].dropna().values * 100 for c in CATEGORY_ORDER]
    f_pro, p_pro = stats.f_oneway(*cat_vals_pro)
    kw_pro, kw_p_pro = stats.kruskal(*cat_vals_pro)
    ns = [len(v) for v in cat_vals_pro]
    means_pro = [np.mean(v) for v in cat_vals_pro]
    sds_pro = [np.std(v, ddof=1) if len(v) > 1 else 0.0 for v in cat_vals_pro]
    meds_pro = [np.median(v) for v in cat_vals_pro]
    iqrs_pro = [stats.iqr(v) for v in cat_vals_pro]

    statistical_results.append({
        "FigurePanel": panel_id,
        "FigureDescription": f"{desc} ({pw}): Promoting direction (Pro)",
        "ModalityOrGene": pw,
        "FactorOrDirection": "Promoting (Pro)",
        "AnalysisLevel": f"Sample-level (N={sum(ns)})",
        "TestMethod": "One-way ANOVA",
        "N_non_tumor_tissues": ns[0],
        "N_tumor_tissues": ns[1],
        "N_normal_cell_lines": ns[2],
        "N_cancer_cell_lines": ns[3],
        "Mean_SD_non_tumor": f"{means_pro[0]:.2f}% ± {sds_pro[0]:.2f}%",
        "Mean_SD_tumor": f"{means_pro[1]:.2f}% ± {sds_pro[1]:.2f}%",
        "Mean_SD_normal_cells": f"{means_pro[2]:.2f}% ± {sds_pro[2]:.2f}%",
        "Mean_SD_cancer_cells": f"{means_pro[3]:.2f}% ± {sds_pro[3]:.2f}%",
        "Median_IQR_non_tumor": f"{meds_pro[0]:.2f}% ({iqrs_pro[0]:.2f}%)",
        "Median_IQR_tumor": f"{meds_pro[1]:.2f}% ({iqrs_pro[1]:.2f}%)",
        "Median_IQR_normal_cells": f"{meds_pro[2]:.2f}% ({iqrs_pro[2]:.2f}%)",
        "Median_IQR_cancer_cells": f"{meds_pro[3]:.2f}% ({iqrs_pro[3]:.2f}%)",
        "TestStatistic": round(f_pro, 2) if np.isfinite(f_pro) else np.nan,
        "StatisticType": "F",
        "df1": len(CATEGORY_ORDER) - 1,
        "df2": sum(ns) - len(CATEGORY_ORDER),
        "RawPValue": p_pro if np.isfinite(p_pro) else np.nan,
        "KruskalStatistic": round(kw_pro, 2) if np.isfinite(kw_pro) else np.nan,
        "KruskalPValue": kw_p_pro if np.isfinite(kw_p_pro) else np.nan
    })

    # Inhibiting (Inh)
    cat_vals_inh = [sub[sub['Category'] == c]['NegativeFraction'].dropna().values * 100 for c in CATEGORY_ORDER]
    if all(np.all(v == 0) for v in cat_vals_inh):
        f_inh, p_inh = np.nan, np.nan
        kw_inh, kw_p_inh = np.nan, np.nan
    else:
        f_inh, p_inh = stats.f_oneway(*cat_vals_inh)
        kw_inh, kw_p_inh = stats.kruskal(*cat_vals_inh)
    means_inh = [np.mean(v) for v in cat_vals_inh]
    sds_inh = [np.std(v, ddof=1) if len(v) > 1 else 0.0 for v in cat_vals_inh]
    meds_inh = [np.median(v) for v in cat_vals_inh]
    iqrs_inh = [stats.iqr(v) for v in cat_vals_inh]

    statistical_results.append({
        "FigurePanel": panel_id,
        "FigureDescription": f"{desc} ({pw}): Inhibiting direction (Inh)",
        "ModalityOrGene": pw,
        "FactorOrDirection": "Inhibiting (Inh)",
        "AnalysisLevel": f"Sample-level (N={sum(ns)})",
        "TestMethod": "One-way ANOVA" if np.isfinite(f_inh) else "Zero-variance (no inhibitory proteins)",
        "N_non_tumor_tissues": ns[0],
        "N_tumor_tissues": ns[1],
        "N_normal_cell_lines": ns[2],
        "N_cancer_cell_lines": ns[3],
        "Mean_SD_non_tumor": f"{means_inh[0]:.2f}% ± {sds_inh[0]:.2f}%",
        "Mean_SD_tumor": f"{means_inh[1]:.2f}% ± {sds_inh[1]:.2f}%",
        "Mean_SD_normal_cells": f"{means_inh[2]:.2f}% ± {sds_inh[2]:.2f}%",
        "Mean_SD_cancer_cells": f"{means_inh[3]:.2f}% ± {sds_inh[3]:.2f}%",
        "Median_IQR_non_tumor": f"{meds_inh[0]:.2f}% ({iqrs_inh[0]:.2f}%)",
        "Median_IQR_tumor": f"{meds_inh[1]:.2f}% ({iqrs_inh[1]:.2f}%)",
        "Median_IQR_normal_cells": f"{meds_inh[2]:.2f}% ({iqrs_inh[2]:.2f}%)",
        "Median_IQR_cancer_cells": f"{meds_inh[3]:.2f}% ({iqrs_inh[3]:.2f}%)",
        "TestStatistic": round(f_inh, 2) if np.isfinite(f_inh) else np.nan,
        "StatisticType": "F" if np.isfinite(f_inh) else "N/A",
        "df1": len(CATEGORY_ORDER) - 1,
        "df2": sum(ns) - len(CATEGORY_ORDER),
        "RawPValue": p_inh if np.isfinite(p_inh) else np.nan,
        "KruskalStatistic": round(kw_inh, 2) if np.isfinite(kw_inh) else np.nan,
        "KruskalPValue": kw_p_inh if np.isfinite(kw_p_inh) else np.nan
    })

    # Paired comparison: Pro vs Inh across all samples
    diff = sub['PositiveFraction'].values * 100 - sub['NegativeFraction'].values * 100
    t_paired, p_paired = stats.ttest_rel(sub['PositiveFraction'].values * 100, sub['NegativeFraction'].values * 100)
    statistical_results.append({
        "FigurePanel": panel_id,
        "FigureDescription": f"{desc} ({pw}): Directional contrast (Pro vs Inh paired)",
        "ModalityOrGene": pw,
        "FactorOrDirection": "Pro vs Inh (Paired)",
        "AnalysisLevel": f"Paired samples (N={len(sub)})",
        "TestMethod": "Paired t-test",
        "N_non_tumor_tissues": len(sub[sub['Category'] == 'normal_tissue']),
        "N_tumor_tissues": len(sub[sub['Category'] == 'cancer_tissue']),
        "N_normal_cell_lines": len(sub[sub['Category'] == 'normal_cells']),
        "N_cancer_cell_lines": len(sub[sub['Category'] == 'cancer_cells']),
        "Mean_SD_non_tumor": f"Diff: {np.mean(diff[sub['Category'] == 'normal_tissue']):.2f}%",
        "Mean_SD_tumor": f"Diff: {np.mean(diff[sub['Category'] == 'cancer_tissue']):.2f}%",
        "Mean_SD_normal_cells": f"Diff: {np.mean(diff[sub['Category'] == 'normal_cells']):.2f}%",
        "Mean_SD_cancer_cells": f"Diff: {np.mean(diff[sub['Category'] == 'cancer_cells']):.2f}%",
        "Median_IQR_non_tumor": "N/A",
        "Median_IQR_tumor": "N/A",
        "Median_IQR_normal_cells": "N/A",
        "Median_IQR_cancer_cells": "N/A",
        "TestStatistic": round(t_paired, 2) if np.isfinite(t_paired) else np.nan,
        "StatisticType": "t",
        "df1": len(sub) - 1,
        "df2": np.nan,
        "RawPValue": p_paired if np.isfinite(p_paired) else np.nan,
        "KruskalStatistic": np.nan,
        "KruskalPValue": np.nan
    })

    for _, r in sub.iterrows():
        # Pro
        tidy_records.append({
            "FigureID": "Supplementary Figure S8",
            "PanelID": panel_id,
            "AssayModality": "Lactylome (Kla)",
            "GroupID": r.get('DisplayGroup', 'N/A'),
            "PXD": r['PXD'],
            "SampleID": r['SampleID'],
            "CategoryCode": r['Category'],
            "CategoryLabel": CATEGORY_LABELS.get(r['Category'], r['Category']),
            "MaterialLabel": r.get('SampleGroup', 'N/A'),
            "FeatureName": f"{pw}_promoting_fraction",
            "Direction": "Promoting (Pro)",
            "ObservationType": "sample",
            "Value": r['PositiveFraction'] * 100,
            "Unit": "Percentage (%)",
            "SourceFile": "figure1_pathway_summary_sample_boxplot_values.csv"
        })
        # Inh
        tidy_records.append({
            "FigureID": "Supplementary Figure S8",
            "PanelID": panel_id,
            "AssayModality": "Lactylome (Kla)",
            "GroupID": r.get('DisplayGroup', 'N/A'),
            "PXD": r['PXD'],
            "SampleID": r['SampleID'],
            "CategoryCode": r['Category'],
            "CategoryLabel": CATEGORY_LABELS.get(r['Category'], r['Category']),
            "MaterialLabel": r.get('SampleGroup', 'N/A'),
            "FeatureName": f"{pw}_inhibiting_fraction",
            "Direction": "Inhibiting (Inh)",
            "ObservationType": "sample",
            "Value": r['NegativeFraction'] * 100,
            "Unit": "Percentage (%)",
            "SourceFile": "figure1_pathway_summary_sample_boxplot_values.csv"
        })

# ==============================================================================
# 5. FIG 3A & S5A-B: RNA EXPRESSION OF LACTATE METABOLISM & REGULATORS
# ==============================================================================
print("Processing Fig 3A & S5 (RNA Expression)...")
rna_panels = [
    ("Fig 3A", ["LDHA", "LDHB", "AARS1", "AARS2"], "Lactate metabolism & lactyl-tRNA synthetases"),
    ("S5A", ["EP300", "CREBBP", "KAT5"], "Lactyltransferases / KATs"),
    ("S5B", ["HDAC1", "HDAC2", "HDAC3", "SIRT2", "SIRT3"], "Delactylases (HDACs & SIRTs)")
]

for panel_id, genes, group_desc in rna_panels:
    for g in genes:
        sub = df_rna[df_rna['Symbol'] == g]
        cat_vals = [sub[sub['Category'] == c]['Expression'].dropna().values for c in CATEGORY_ORDER]
        f_stat, p_val = stats.f_oneway(*cat_vals)
        kw_stat, kw_pval = stats.kruskal(*cat_vals)
        ns = [len(v) for v in cat_vals]
        means = [np.mean(v) for v in cat_vals]
        sds = [np.std(v, ddof=1) if len(v) > 1 else 0.0 for v in cat_vals]
        meds = [np.median(v) for v in cat_vals]
        iqrs = [stats.iqr(v) for v in cat_vals]

        statistical_results.append({
            "FigurePanel": panel_id,
            "FigureDescription": f"{group_desc}: {g} transcriptomic expression",
            "ModalityOrGene": g,
            "FactorOrDirection": "Category",
            "AnalysisLevel": f"Sample-level (N={sum(ns)})",
            "TestMethod": "One-way ANOVA",
            "N_non_tumor_tissues": ns[0],
            "N_tumor_tissues": ns[1],
            "N_normal_cell_lines": ns[2],
            "N_cancer_cell_lines": ns[3],
            "Mean_SD_non_tumor": f"{means[0]:.4f} ± {sds[0]:.4f}",
            "Mean_SD_tumor": f"{means[1]:.4f} ± {sds[1]:.4f}",
            "Mean_SD_normal_cells": f"{means[2]:.4f} ± {sds[2]:.4f}",
            "Mean_SD_cancer_cells": f"{means[3]:.4f} ± {sds[3]:.4f}",
            "Median_IQR_non_tumor": f"{meds[0]:.4f} ({iqrs[0]:.4f})",
            "Median_IQR_tumor": f"{meds[1]:.4f} ({iqrs[1]:.4f})",
            "Median_IQR_normal_cells": f"{meds[2]:.4f} ({iqrs[2]:.4f})",
            "Median_IQR_cancer_cells": f"{meds[3]:.4f} ({iqrs[3]:.4f})",
            "TestStatistic": round(f_stat, 2),
            "StatisticType": "F",
            "df1": len(CATEGORY_ORDER) - 1,
            "df2": sum(ns) - len(CATEGORY_ORDER),
            "RawPValue": p_val,
            "KruskalStatistic": round(kw_stat, 2),
            "KruskalPValue": kw_pval
        })

        for _, r in sub.iterrows():
            tidy_records.append({
                "FigureID": "Figure 3 / S5",
                "PanelID": panel_id,
                "AssayModality": "Transcriptome (RNA-seq)",
                "GroupID": r.get('ReferenceKey', 'N/A'),
                "PXD": "N/A",
                "SampleID": r['SampleID'],
                "CategoryCode": r['Category'],
                "CategoryLabel": CATEGORY_LABELS.get(r['Category'], r['Category']),
                "MaterialLabel": r['ReferenceKey'],
                "FeatureName": g,
                "Direction": "N/A",
                "ObservationType": "sample",
                "Value": r['Expression'],
                "Unit": "log2(qsmooth TPM + 0.5)",
                "SourceFile": "rna_values_28_references.csv.gz"
            })

# ==============================================================================
# 6. S6A & S6B: DDR & KLA TARGET GENE EXPRESSION & COVERAGE FRACTION
# ==============================================================================
print("Processing S6A & S6B...")
# S6A
for gs in df_s6a['GeneSet'].unique():
    sub = df_s6a[df_s6a['GeneSet'] == gs]
    cat_vals = [sub[sub['Category'] == c]['MedianExpression'].dropna().values for c in CATEGORY_ORDER]
    f_stat, p_val = stats.f_oneway(*cat_vals)
    kw_stat, kw_pval = stats.kruskal(*cat_vals)
    ns = [len(v) for v in cat_vals]
    means = [np.mean(v) for v in cat_vals]
    sds = [np.std(v, ddof=1) if len(v) > 1 else 0.0 for v in cat_vals]
    meds = [np.median(v) for v in cat_vals]
    iqrs = [stats.iqr(v) for v in cat_vals]

    statistical_results.append({
        "FigurePanel": "S6A",
        "FigureDescription": f"Transcriptomic expression: {gs} (28 independent RNA materials)",
        "ModalityOrGene": gs,
        "FactorOrDirection": "Category",
        "AnalysisLevel": f"Material-level (N={sum(ns)})",
        "TestMethod": "One-way ANOVA",
        "N_non_tumor_tissues": ns[0],
        "N_tumor_tissues": ns[1],
        "N_normal_cell_lines": ns[2],
        "N_cancer_cell_lines": ns[3],
        "Mean_SD_non_tumor": f"{means[0]:.4f} ± {sds[0]:.4f}",
        "Mean_SD_tumor": f"{means[1]:.4f} ± {sds[1]:.4f}",
        "Mean_SD_normal_cells": f"{means[2]:.4f} ± {sds[2]:.4f}",
        "Mean_SD_cancer_cells": f"{means[3]:.4f} ± {sds[3]:.4f}",
        "Median_IQR_non_tumor": f"{meds[0]:.4f} ({iqrs[0]:.4f})",
        "Median_IQR_tumor": f"{meds[1]:.4f} ({iqrs[1]:.4f})",
        "Median_IQR_normal_cells": f"{meds[2]:.4f} ({iqrs[2]:.4f})",
        "Median_IQR_cancer_cells": f"{meds[3]:.4f} ({iqrs[3]:.4f})",
        "TestStatistic": round(f_stat, 2),
        "StatisticType": "F",
        "df1": len(CATEGORY_ORDER) - 1,
        "df2": sum(ns) - len(CATEGORY_ORDER),
        "RawPValue": p_val,
        "KruskalStatistic": round(kw_stat, 2),
        "KruskalPValue": kw_pval
    })

    for _, r in sub.iterrows():
        tidy_records.append({
            "FigureID": "Supplementary Figure S6",
            "PanelID": "S6A",
            "AssayModality": "Transcriptome (RNA-seq)",
            "GroupID": r['GroupID'],
            "PXD": "N/A",
            "SampleID": r['GroupID'],
            "CategoryCode": r['Category'],
            "CategoryLabel": CATEGORY_LABELS.get(r['Category'], r['Category']),
            "MaterialLabel": r['GroupID'],
            "FeatureName": gs,
            "Direction": "N/A",
            "ObservationType": "material_median",
            "Value": r['MedianExpression'],
            "Unit": "log2(qsmooth TPM + 0.5)",
            "SourceFile": "figure1a_rna_ddr_and_lactylated_28materials.csv"
        })

# S6B
cat_vals_s6b = [df_s6b[df_s6b['Category'] == c]['DdrFractionMedian'].dropna().values for c in CATEGORY_ORDER]
f_stat, p_val = stats.f_oneway(*cat_vals_s6b)
kw_stat, kw_pval = stats.kruskal(*cat_vals_s6b)
ns = [len(v) for v in cat_vals_s6b]
means = [np.mean(v) for v in cat_vals_s6b]
sds = [np.std(v, ddof=1) if len(v) > 1 else 0.0 for v in cat_vals_s6b]
meds = [np.median(v) for v in cat_vals_s6b]
iqrs = [stats.iqr(v) for v in cat_vals_s6b]

statistical_results.append({
    "FigurePanel": "S6B",
    "FigureDescription": "Transcriptomic DDR expressed gene fraction (%) (28 independent RNA materials)",
    "ModalityOrGene": "GO-DDR fraction in transcriptome",
    "FactorOrDirection": "Category",
    "AnalysisLevel": f"Material-level (N={sum(ns)})",
    "TestMethod": "One-way ANOVA",
    "N_non_tumor_tissues": ns[0],
    "N_tumor_tissues": ns[1],
    "N_normal_cell_lines": ns[2],
    "N_cancer_cell_lines": ns[3],
    "Mean_SD_non_tumor": f"{means[0]:.2f}% ± {sds[0]:.2f}%",
    "Mean_SD_tumor": f"{means[1]:.2f}% ± {sds[1]:.2f}%",
    "Mean_SD_normal_cells": f"{means[2]:.2f}% ± {sds[2]:.2f}%",
    "Mean_SD_cancer_cells": f"{means[3]:.2f}% ± {sds[3]:.2f}%",
    "Median_IQR_non_tumor": f"{meds[0]:.2f}% ({iqrs[0]:.2f}%)",
    "Median_IQR_tumor": f"{meds[1]:.2f}% ({iqrs[1]:.2f}%)",
    "Median_IQR_normal_cells": f"{meds[2]:.2f}% ({iqrs[2]:.2f}%)",
    "Median_IQR_cancer_cells": f"{meds[3]:.2f}% ({iqrs[3]:.2f}%)",
    "TestStatistic": round(f_stat, 2),
    "StatisticType": "F",
    "df1": len(CATEGORY_ORDER) - 1,
    "df2": sum(ns) - len(CATEGORY_ORDER),
    "RawPValue": p_val,
    "KruskalStatistic": round(kw_stat, 2),
    "KruskalPValue": kw_pval
})

for _, r in df_s6b.iterrows():
    tidy_records.append({
        "FigureID": "Supplementary Figure S6",
        "PanelID": "S6B",
        "AssayModality": "Transcriptome (RNA-seq)",
        "GroupID": r['GroupID'],
        "PXD": "N/A",
        "SampleID": r['GroupID'],
        "CategoryCode": r['Category'],
        "CategoryLabel": CATEGORY_LABELS.get(r['Category'], r['Category']),
        "MaterialLabel": r.get('ReferenceLabelEn', r['GroupID']),
        "FeatureName": "GO-DDR_expressed_gene_fraction",
        "Direction": "N/A",
        "ObservationType": "material_fraction",
        "Value": r['DdrFractionMedian'],
        "Unit": "Percentage (%)",
        "SourceFile": "figure1a_fixed_ddr_28materials.csv"
    })

# Compute Multiple Testing Correction (Benjamini-Hochberg FDR)
df_stats = pd.DataFrame(statistical_results)
valid_p_mask = df_stats['RawPValue'].notna() & np.isfinite(df_stats['RawPValue'])
df_stats['FDR_QValue'] = np.nan
if valid_p_mask.sum() > 0:
    df_stats.loc[valid_p_mask, 'FDR_QValue'] = stats.false_discovery_control(df_stats.loc[valid_p_mask, 'RawPValue'].values, method='bh')

df_stats['Significance'] = df_stats['RawPValue'].apply(format_sig_stars)
df_stats['FDR_Significance'] = df_stats['FDR_QValue'].apply(format_sig_stars)

# Save summary CSV
csv_summary_path = os.path.join(OUTPUT_DIR, "statistical_test_summary.csv")
df_stats.to_csv(csv_summary_path, index=False)
print(f"Saved statistical summary CSV: {csv_summary_path}")

# Save tidy CSV
df_tidy = pd.DataFrame(tidy_records)
csv_tidy_path = os.path.join(OUTPUT_DIR, "all_data_tidy_provenance.csv")
df_tidy.to_csv(csv_tidy_path, index=False)
print(f"Saved tidy provenance CSV: {csv_tidy_path} ({len(df_tidy)} records)")

# ==============================================================================
# 7. BUILD MASTER EXCEL WORKBOOK (OPENPYXL)
# ==============================================================================
print("Building Master Excel Workbook for GraphPad Prism...")
wb = openpyxl.Workbook()
wb.remove(wb.active)  # Remove default empty sheet

# ------------------------------------------------------------------------------
# SHEET 0: README & INDEX
# ------------------------------------------------------------------------------
ws_idx = wb.create_sheet(title="00_目录与Prism使用说明")
ws_idx.views.sheetView[0].showGridLines = True

headers_idx = ["Sheet名称", "对应论文图表", "数据类型/分析对象", "推荐Prism表格类型", "检验方法", "图表说明"]
idx_rows = [
    ("01_统计检验结果汇总", "全量", "所有图表面板推断统计量汇总", "Summary Table", "ANOVA / Kruskal / t-test", "完整列出F值、p值、BH校正q值、自由度及各组均值±SD"),
    ("Fig1B_DDR组成_Column", "Figure 1B", "四类别DDR组成 (Whole vs Kla)", "Column Table (One-way ANOVA)", "One-way ANOVA / Kruskal-Wallis", "包含31材料预汇总版与307源单位版两组立式列"),
    ("Fig1B_DDR组成_Grouped", "Figure 1B", "四类别DDR组成 (分组双箱图)", "Grouped Table (Two-way ANOVA)", "Two-way ANOVA (Category × Modality)", "行=4大生物类别，列=Whole proteome与Lactylome(Kla)子列重复"),
    ("S3A_C_MKI67比值_Column", "Supplementary S3A–C", "全蛋白组MKI67增殖比值 (H3C1/ACTB/TUBB)", "Column Table (One-way ANOVA)", "One-way ANOVA", "三个分母分别成表：H3C1 (S3A), ACTB (S3B), TUBB (S3C)"),
    ("S3A_C_MKI67比值_Grouped", "Supplementary S3A–C", "全蛋白组MKI67增殖比值 (分组表)", "Grouped Table (Two-way ANOVA)", "Two-way ANOVA", "行=4大生物类别，列=H3C1, ACTB, TUBB子列重复"),
    ("S3D_RNA_G2M增殖_Column", "Supplementary S3D", "转录组Hallmark G2M checkpoint", "Column Table (One-way ANOVA)", "One-way ANOVA", "提供28独立材料中位数与1898样本两组独立列"),
    ("S8A_BER_七通路", "Supplementary S8A", "BER 碱基切除修复 (Pro vs Inh)", "Grouped Table (模板格式)", "Two-way ANOVA / Paired t-test", "严格对应模板结构：第一行重复Pro/Inh，第一列为4大类别"),
    ("S8B_NER_七通路", "Supplementary S8B", "NER 核苷酸切除修复 (Pro vs Inh)", "Grouped Table (模板格式)", "Two-way ANOVA / Paired t-test", "严格对应模板结构：第一行重复Pro/Inh，第一列为4大类别"),
    ("S8C_MMR_七通路", "Supplementary S8C", "MMR 错配修复 (Pro vs Inh)", "Grouped Table (模板格式)", "Two-way ANOVA / Paired t-test", "严格对应模板结构：第一行重复Pro/Inh，第一列为4大类别"),
    ("S8D_FA_七通路", "Supplementary S8D", "FA 范可尼贫血通路 (Pro vs Inh)", "Grouped Table (模板格式)", "Two-way ANOVA / Paired t-test", "严格对应模板结构：第一行重复Pro/Inh，第一列为4大类别"),
    ("S8E_HR_七通路", "Supplementary S8E", "HR 同源重组修复 (Pro vs Inh)", "Grouped Table (模板格式)", "Two-way ANOVA / Paired t-test", "严格对应模板结构：第一行重复Pro/Inh，第一列为4大类别"),
    ("S8F_AEJ_七通路", "Supplementary S8F", "AEJ 微同源末端连接 (Pro vs Inh)", "Grouped Table (模板格式)", "Two-way ANOVA / Paired t-test", "严格对应模板结构：第一行重复Pro/Inh，第一列为4大类别"),
    ("S8G_NHEJ_七通路", "Supplementary S8G", "NHEJ 非同源末端连接 (Pro vs Inh)", "Grouped Table (模板格式)", "Two-way ANOVA / Paired t-test", "严格对应模板结构：第一行重复Pro/Inh，第一列为4大类别"),
    ("Fig3A_乳酸与AARS_RNA", "Figure 3A", "LDHA, LDHB, AARS1, AARS2 表达", "Grouped / Column Table", "One-way ANOVA / Kruskal-Wallis", "行=4大类别，列=4个基因及其在1898个样本中的表达值"),
    ("S5A_乙酰乳酰转移酶_RNA", "Supplementary S5A", "EP300, CREBBP, KAT5 表达", "Grouped / Column Table", "One-way ANOVA / Kruskal-Wallis", "行=4大类别，列=3个转移酶基因样本表达值"),
    ("S5B_去乳酸化酶_RNA", "Supplementary S5B", "HDAC1-3, SIRT2-3 表达", "Grouped / Column Table", "One-way ANOVA / Kruskal-Wallis", "行=4大类别，列=5个去修饰酶基因样本表达值"),
    ("S6A_DDR与Kla靶集表达", "Supplementary S6A", "DDR 371靶集 vs Kla 5224靶集表达", "Grouped Table", "Two-way ANOVA / Paired t-test", "行=4大类别，列=DDR基因集与Kla靶集在28材料中的评分"),
    ("S6B_DDR转录覆盖占比", "Supplementary S6B", "转录组GO-DDR表达基因覆盖率 (%)", "Column Table (One-way ANOVA)", "One-way ANOVA", "4个类别列，28个材料在qSmooth空间达到阈值的覆盖占比"),
    ("All_Data_Tidy_底表", "全量溯源", "所有统计分析底层数据记录", "Tidy / Long Table", "N/A", "完整包含PXD、SampleID、类别、数值、单位与源文件名")
]

ws_idx.cell(1, 1, "Kla × DDR 课题 GraphPad Prism 统计数据表使用指南").font = FONT_TITLE
ws_idx.cell(2, 1, "数据基准：2026-10-11 全组核查定稿版（31材料、ESCC纳入、307蛋白质组源单位、28转录组独立材料、1898样本）").font = FONT_NOTE
ws_idx.cell(3, 1, "格式设计：严格遵循 GraphPad Prism 导入标准与 sample_linear_qsmooth_template_symbol_log2p1.xlsx 规范").font = FONT_NOTE

for c_idx, h in enumerate(headers_idx, 1):
    cell = ws_idx.cell(5, c_idx, h)
    cell.font = FONT_HEADER
    cell.border = BORDER_HEADER
    cell.alignment = ALIGN_LEFT

for r_idx, row in enumerate(idx_rows, 6):
    for c_idx, val in enumerate(row, 1):
        cell = ws_idx.cell(r_idx, c_idx, val)
        cell.font = FONT_BODY
        cell.border = BORDER_THIN
        cell.alignment = ALIGN_LEFT

auto_fit_columns(ws_idx)

# ------------------------------------------------------------------------------
# SHEET 1: STATISTICAL TEST SUMMARY TABLE
# ------------------------------------------------------------------------------
ws_sum = wb.create_sheet(title="01_统计检验结果汇总")
ws_sum.views.sheetView[0].showGridLines = True

ws_sum.cell(1, 1, "全量统计检验结果汇总表 (One-way ANOVA / Kruskal-Wallis / Two-way ANOVA)").font = FONT_TITLE
ws_sum.cell(2, 1, "注：FDR校正采用 Benjamini-Hochberg (BH) 算法；显著性代码：**** p<0.0001, *** p<0.001, ** p<0.01, * p<0.05, ns p>=0.05").font = FONT_NOTE

headers_sum = [
    "图表面板", "图表分析说明", "分析模态/基因/通路", "比较因素/方向", "分析层次",
    "检验方法", "N_非肿瘤组织", "N_肿瘤组织", "N_正常细胞", "N_癌细胞",
    "均值±SD (非肿瘤组织)", "均值±SD (肿瘤组织)", "均值±SD (正常细胞)", "均值±SD (癌细胞)",
    "中位数(IQR) (非肿瘤组织)", "中位数(IQR) (肿瘤组织)", "中位数(IQR) (正常细胞)", "中位数(IQR) (癌细胞)",
    "统计量(F/t)", "统计量类型", "分子自由度(df1)", "分母自由度(df2)", "原始P值", "BH校正Q值", "原始显著性", "FDR显著性"
]

for c_idx, h in enumerate(headers_sum, 1):
    cell = ws_sum.cell(4, c_idx, h)
    cell.font = FONT_HEADER
    cell.border = BORDER_HEADER
    cell.alignment = ALIGN_CENTER

for r_idx, r in df_stats.iterrows():
    row_vals = [
        r['FigurePanel'], r['FigureDescription'], r['ModalityOrGene'], r['FactorOrDirection'], r['AnalysisLevel'],
        r['TestMethod'], r['N_non_tumor_tissues'], r['N_tumor_tissues'], r['N_normal_cell_lines'], r['N_cancer_cell_lines'],
        r['Mean_SD_non_tumor'], r['Mean_SD_tumor'], r['Mean_SD_normal_cells'], r['Mean_SD_cancer_cells'],
        r['Median_IQR_non_tumor'], r['Median_IQR_tumor'], r['Median_IQR_normal_cells'], r['Median_IQR_cancer_cells'],
        r['TestStatistic'], r['StatisticType'], r['df1'], r['df2'],
        f"{r['RawPValue']:.4e}" if pd.notna(r['RawPValue']) else "N/A",
        f"{r['FDR_QValue']:.4e}" if pd.notna(r['FDR_QValue']) else "N/A",
        r['Significance'], r['FDR_Significance']
    ]
    for c_idx, val in enumerate(row_vals, 1):
        cell = ws_sum.cell(r_idx + 5, c_idx, val)
        cell.font = FONT_BODY
        cell.border = BORDER_THIN
        cell.alignment = ALIGN_RIGHT if isinstance(val, (int, float)) or (isinstance(val, str) and "%" in val and "±" in val) else ALIGN_LEFT

auto_fit_columns(ws_sum)

# ------------------------------------------------------------------------------
# SHEET 2: FIG 1B (COLUMN FORMAT: 31 MATERIALS & 307 SOURCE UNITS)
# ------------------------------------------------------------------------------
ws_f1_col = wb.create_sheet(title="Fig1B_DDR组成_Column")
ws_f1_col.views.sheetView[0].showGridLines = True

ws_f1_col.cell(1, 1, "Figure 1B: DDR Fraction (Column Table Format for One-way ANOVA)").font = FONT_TITLE
ws_f1_col.cell(2, 1, "Prism 使用：新建 Column 表格，直接复制对应区域。A–D 列为全蛋白组，E–H 列为乳酸化组。").font = FONT_NOTE

# Section A: 31 Materials (Pre-aggregated mean per material)
ws_f1_col.cell(4, 1, "【第 1 部分：31 个材料预汇总值 (推荐正文图表使用，消除样本量偏倚)】").font = FONT_SECTION
headers_f1 = [
    "Whole: non-tumor tissues", "Whole: tumor tissues", "Whole: normal cell lines", "Whole: cancer cell lines",
    "Kla: non-tumor tissues", "Kla: tumor tissues", "Kla: normal cell lines", "Kla: cancer cell lines"
]
for c_idx, h in enumerate(headers_f1, 1):
    cell = ws_f1_col.cell(5, c_idx, h)
    cell.font = FONT_HEADER
    cell.border = BORDER_HEADER
    cell.alignment = ALIGN_CENTER

mat_ddr_dict = {}
for ds in ["Whole proteome", "Lactylome (Kla)"]:
    for c in CATEGORY_ORDER:
        vals = mat_ddr[(mat_ddr['Dataset'] == ds) & (mat_ddr['Category'] == c)]['Fraction'].values
        mat_ddr_dict[(ds, c)] = vals

max_r_mat = max(len(v) for v in mat_ddr_dict.values())
for r_i in range(max_r_mat):
    for c_i, (ds, c) in enumerate([("Whole proteome", c) for c in CATEGORY_ORDER] + [("Lactylome (Kla)", c) for c in CATEGORY_ORDER], 1):
        vals = mat_ddr_dict[(ds, c)]
        val = vals[r_i] if r_i < len(vals) else ""
        cell = ws_f1_col.cell(6 + r_i, c_i, val)
        cell.font = FONT_BODY
        cell.border = BORDER_THIN
        cell.alignment = ALIGN_RIGHT
        if isinstance(val, (int, float)):
            cell.number_format = "0.000"

# Section B: 307 Source Units
start_r_units = 6 + max_r_mat + 3
ws_f1_col.cell(start_r_units, 1, "【第 2 部分：307 个独立源单位观测值 (测序数据逐样本散点观测)】").font = FONT_SECTION
for c_idx, h in enumerate(headers_f1, 1):
    cell = ws_f1_col.cell(start_r_units + 1, c_idx, h)
    cell.font = FONT_HEADER
    cell.border = BORDER_HEADER
    cell.alignment = ALIGN_CENTER

unit_ddr_dict = {}
for ds in ["Whole proteome", "Lactylome (Kla)"]:
    for c in CATEGORY_ORDER:
        vals = df_ddr_units[(df_ddr_units['Dataset'] == ds) & (df_ddr_units['Category'] == c)]['Fraction'].values
        unit_ddr_dict[(ds, c)] = vals

max_r_unit = max(len(v) for v in unit_ddr_dict.values())
for r_i in range(max_r_unit):
    for c_i, (ds, c) in enumerate([("Whole proteome", c) for c in CATEGORY_ORDER] + [("Lactylome (Kla)", c) for c in CATEGORY_ORDER], 1):
        vals = unit_ddr_dict[(ds, c)]
        val = vals[r_i] if r_i < len(vals) else ""
        cell = ws_f1_col.cell(start_r_units + 2 + r_i, c_i, val)
        cell.font = FONT_BODY
        cell.border = BORDER_THIN
        cell.alignment = ALIGN_RIGHT
        if isinstance(val, (int, float)):
            cell.number_format = "0.000"

auto_fit_columns(ws_f1_col)

# ------------------------------------------------------------------------------
# SHEET 3: FIG 1B (GROUPED FORMAT: TEMPLATE STYLE)
# ------------------------------------------------------------------------------
ws_f1_grp = wb.create_sheet(title="Fig1B_DDR组成_Grouped")
ws_f1_grp.views.sheetView[0].showGridLines = True

ws_f1_grp.cell(1, 1, "Figure 1B: DDR Fraction (Grouped Table Format for Two-way ANOVA)").font = FONT_TITLE
ws_f1_grp.cell(2, 1, "Prism 使用：新建 Grouped 表格，行=生物类别，列=检测模态（Whole proteome vs Kla）。").font = FONT_NOTE

# 31 Materials Grouped Block
ws_f1_grp.cell(4, 1, "【31 材料预汇总版：行=生物类别，列=模态，子列=材料】").font = FONT_SECTION
# Find max replicates per category
max_n_whole = max(len(mat_ddr_dict[("Whole proteome", c)]) for c in CATEGORY_ORDER)
max_n_kla = max(len(mat_ddr_dict[("Lactylome (Kla)", c)]) for c in CATEGORY_ORDER)

# Row 5: Column group headers (repeated across subcolumns)
ws_f1_grp.cell(5, 1, "Biological Category").font = FONT_HEADER
ws_f1_grp.cell(5, 1).border = BORDER_HEADER
for i in range(max_n_whole):
    c = ws_f1_grp.cell(5, 2 + i, "Whole proteome")
    c.font = FONT_HEADER
    c.border = BORDER_HEADER
    c.alignment = ALIGN_CENTER
for i in range(max_n_kla):
    c = ws_f1_grp.cell(5, 2 + max_n_whole + i, "Lactylome (Kla)")
    c.font = FONT_HEADER
    c.border = BORDER_HEADER
    c.alignment = ALIGN_CENTER

for r_idx, c in enumerate(CATEGORY_ORDER, 6):
    ws_f1_grp.cell(r_idx, 1, CATEGORY_LABELS[c]).font = FONT_HEADER
    ws_f1_grp.cell(r_idx, 1).border = BORDER_THIN
    w_vals = mat_ddr_dict[("Whole proteome", c)]
    k_vals = mat_ddr_dict[("Lactylome (Kla)", c)]
    for i in range(max_n_whole):
        val = w_vals[i] if i < len(w_vals) else ""
        cell = ws_f1_grp.cell(r_idx, 2 + i, val)
        cell.font = FONT_BODY
        cell.border = BORDER_THIN
        cell.alignment = ALIGN_RIGHT
        if isinstance(val, (int, float)):
            cell.number_format = "0.000"
    for i in range(max_n_kla):
        val = k_vals[i] if i < len(k_vals) else ""
        cell = ws_f1_grp.cell(r_idx, 2 + max_n_whole + i, val)
        cell.font = FONT_BODY
        cell.border = BORDER_THIN
        cell.alignment = ALIGN_RIGHT
        if isinstance(val, (int, float)):
            cell.number_format = "0.000"

auto_fit_columns(ws_f1_grp)

# ------------------------------------------------------------------------------
# SHEET 4: S3A-C MKI67 PROLIFERATION RATIOS (COLUMN FORMAT)
# ------------------------------------------------------------------------------
ws_mki_col = wb.create_sheet(title="S3A_C_MKI67比值_Column")
ws_mki_col.views.sheetView[0].showGridLines = True

ws_mki_col.cell(1, 1, "Supplementary Figure S3A–C: MKI67 Proliferation Ratios (Column Format)").font = FONT_TITLE
ws_mki_col.cell(2, 1, "Prism 使用：新建 Column 表格，比较四个生物分类之间的 MKI67 增殖蛋白强度比。").font = FONT_NOTE

curr_row = 4
for denom, panel_id, desc in mki67_normalizers:
    ws_mki_col.cell(curr_row, 1, f"【{panel_id}: {desc}】").font = FONT_SECTION
    sub = df_mki67_sample[df_mki67_sample['Denominator'] == denom]
    for c_i, c in enumerate(CATEGORY_ORDER, 1):
        cell = ws_mki_col.cell(curr_row + 1, c_i, CATEGORY_LABELS[c])
        cell.font = FONT_HEADER
        cell.border = BORDER_HEADER
        cell.alignment = ALIGN_CENTER

    mki_cat_vals = [sub[sub['Category'] == c]['Ratio'].values for c in CATEGORY_ORDER]
    max_len = max(len(v) for v in mki_cat_vals)
    for r_i in range(max_len):
        for c_i, vals in enumerate(mki_cat_vals, 1):
            val = vals[r_i] if r_i < len(vals) else ""
            cell = ws_mki_col.cell(curr_row + 2 + r_i, c_i, val)
            cell.font = FONT_BODY
            cell.border = BORDER_THIN
            cell.alignment = ALIGN_RIGHT
            if isinstance(val, (int, float)):
                cell.number_format = "0.0000"
    curr_row += max_len + 4

auto_fit_columns(ws_mki_col)

# ------------------------------------------------------------------------------
# SHEET 5: S3A-C MKI67 PROLIFERATION RATIOS (GROUPED FORMAT)
# ------------------------------------------------------------------------------
ws_mki_grp = wb.create_sheet(title="S3A_C_MKI67比值_Grouped")
ws_mki_grp.views.sheetView[0].showGridLines = True

ws_mki_grp.cell(1, 1, "Supplementary Figure S3A–C: MKI67 Proliferation Ratios (Grouped Format)").font = FONT_TITLE
ws_mki_grp.cell(2, 1, "Prism 使用：新建 Grouped 表格，行=生物类别，列=归一化内参蛋白 (H3C1, ACTB, TUBB)。").font = FONT_NOTE

# Find max replicates per normalizer
max_ns = [max(len(df_mki67_sample[(df_mki67_sample['Denominator'] == d) & (df_mki67_sample['Category'] == c)]['Ratio']) for c in CATEGORY_ORDER) for d, _, _ in mki67_normalizers]

ws_mki_grp.cell(4, 1, "Biological Category").font = FONT_HEADER
ws_mki_grp.cell(4, 1).border = BORDER_HEADER

col_offset = 2
for (denom, _, _), max_n in zip(mki67_normalizers, max_ns):
    for i in range(max_n):
        cell = ws_mki_grp.cell(4, col_offset + i, f"MKI67/{denom}")
        cell.font = FONT_HEADER
        cell.border = BORDER_HEADER
        cell.alignment = ALIGN_CENTER
    col_offset += max_n

for r_idx, c in enumerate(CATEGORY_ORDER, 5):
    ws_mki_grp.cell(r_idx, 1, CATEGORY_LABELS[c]).font = FONT_HEADER
    ws_mki_grp.cell(r_idx, 1).border = BORDER_THIN
    col_offset = 2
    for (denom, _, _), max_n in zip(mki67_normalizers, max_ns):
        vals = df_mki67_sample[(df_mki67_sample['Denominator'] == denom) & (df_mki67_sample['Category'] == c)]['Ratio'].values
        for i in range(max_n):
            val = vals[i] if i < len(vals) else ""
            cell = ws_mki_grp.cell(r_idx, col_offset + i, val)
            cell.font = FONT_BODY
            cell.border = BORDER_THIN
            cell.alignment = ALIGN_RIGHT
            if isinstance(val, (int, float)):
                cell.number_format = "0.0000"
        col_offset += max_n

auto_fit_columns(ws_mki_grp)

# ------------------------------------------------------------------------------
# SHEET 6: S3D RNA HALLMARK G2M CHECKPOINT (COLUMN FORMAT)
# ------------------------------------------------------------------------------
ws_g2m = wb.create_sheet(title="S3D_RNA_G2M增殖_Column")
ws_g2m.views.sheetView[0].showGridLines = True

ws_g2m.cell(1, 1, "Supplementary Figure S3D: Transcriptomic Hallmark G2M Checkpoint Expression").font = FONT_TITLE
ws_g2m.cell(2, 1, "单位：log2(qsmooth TPM + 0.5)。包含28材料中位数与1898独立样本两套数据。").font = FONT_NOTE

# Section 1: 28 materials
ws_g2m.cell(4, 1, "【第 1 部分：28 个独立 RNA 参考材料中位数评分 (推荐正文图表使用)】").font = FONT_SECTION
for c_i, c in enumerate(CATEGORY_ORDER, 1):
    cell = ws_g2m.cell(5, c_i, CATEGORY_LABELS[c])
    cell.font = FONT_HEADER
    cell.border = BORDER_HEADER
    cell.alignment = ALIGN_CENTER

max_len_mat = max(len(v) for v in g2m_mat_vals)
for r_i in range(max_len_mat):
    for c_i, vals in enumerate(g2m_mat_vals, 1):
        val = vals[r_i] if r_i < len(vals) else ""
        cell = ws_g2m.cell(6 + r_i, c_i, val)
        cell.font = FONT_BODY
        cell.border = BORDER_THIN
        cell.alignment = ALIGN_RIGHT
        if isinstance(val, (int, float)):
            cell.number_format = "0.000"

# Section 2: 1898 samples
start_r_samp = 6 + max_len_mat + 3
ws_g2m.cell(start_r_samp, 1, "【第 2 部分：1,898 个独立 RNA 样本散点观测值】").font = FONT_SECTION
for c_i, c in enumerate(CATEGORY_ORDER, 1):
    cell = ws_g2m.cell(start_r_samp + 1, c_i, CATEGORY_LABELS[c])
    cell.font = FONT_HEADER
    cell.border = BORDER_HEADER
    cell.alignment = ALIGN_CENTER

max_len_samp = max(len(v) for v in g2m_samp_vals)
for r_i in range(max_len_samp):
    for c_i, vals in enumerate(g2m_samp_vals, 1):
        val = vals[r_i] if r_i < len(vals) else ""
        cell = ws_g2m.cell(start_r_samp + 2 + r_i, c_i, val)
        cell.font = FONT_BODY
        cell.border = BORDER_THIN
        cell.alignment = ALIGN_RIGHT
        if isinstance(val, (int, float)):
            cell.number_format = "0.000"

auto_fit_columns(ws_g2m)

# ------------------------------------------------------------------------------
# SHEETS 7-13: S8A-G SEVEN DDR PATHWAYS (PRO / INH GROUPED TEMPLATE)
# ------------------------------------------------------------------------------
print("Building S8A-G Pathway sheets matching template structure...")
for pw, panel_id, desc in pathways:
    ws_pw = wb.create_sheet(title=f"{panel_id}_{pw}_七通路")
    ws_pw.views.sheetView[0].showGridLines = True

    ws_pw.cell(1, 1, f"Supplementary Figure {panel_id}: {desc} ({pw}) Pathway Fraction").font = FONT_TITLE
    ws_pw.cell(2, 1, "格式严格对齐 sample_linear_qsmooth_template_symbol_log2p1.xlsx：行=类别，列=方向（Promoting/Inhibiting），子列=样本。单位：百分比 (%)。").font = FONT_NOTE

    sub = df_pathway_sample[df_pathway_sample['Pathway'] == pw]

    # Find max replicates per category
    pw_cat_dict = {}
    for c in CATEGORY_ORDER:
        pro_vals = sub[sub['Category'] == c]['PositiveFraction'].values * 100
        inh_vals = sub[sub['Category'] == c]['NegativeFraction'].values * 100
        pw_cat_dict[c] = (pro_vals, inh_vals)

    max_pro = max(len(v[0]) for v in pw_cat_dict.values())
    max_inh = max(len(v[1]) for v in pw_cat_dict.values())

    # Row 4: Repeated Direction headers (Template style: Row 1 = Repeated gene/factor names)
    ws_pw.cell(4, 1, "Biological Category").font = FONT_HEADER
    ws_pw.cell(4, 1).border = BORDER_HEADER
    for i in range(max_pro):
        cell = ws_pw.cell(4, 2 + i, "Promoting (Pro)")
        cell.font = FONT_HEADER
        cell.border = BORDER_HEADER
        cell.alignment = ALIGN_CENTER
    for i in range(max_inh):
        cell = ws_pw.cell(4, 2 + max_pro + i, "Inhibiting (Inh)")
        cell.font = FONT_HEADER
        cell.border = BORDER_HEADER
        cell.alignment = ALIGN_CENTER

    for r_idx, c in enumerate(CATEGORY_ORDER, 5):
        ws_pw.cell(r_idx, 1, CATEGORY_LABELS[c]).font = FONT_HEADER
        ws_pw.cell(r_idx, 1).border = BORDER_THIN
        pro_v, inh_v = pw_cat_dict[c]
        for i in range(max_pro):
            val = pro_v[i] if i < len(pro_v) else ""
            cell = ws_pw.cell(r_idx, 2 + i, val)
            cell.font = FONT_BODY
            cell.border = BORDER_THIN
            cell.alignment = ALIGN_RIGHT
            if isinstance(val, (int, float)):
                cell.number_format = "0.00"
        for i in range(max_inh):
            val = inh_v[i] if i < len(inh_v) else ""
            cell = ws_pw.cell(r_idx, 2 + max_pro + i, val)
            cell.font = FONT_BODY
            cell.border = BORDER_THIN
            cell.alignment = ALIGN_RIGHT
            if isinstance(val, (int, float)):
                cell.number_format = "0.00"

    # Also provide sample-resolved flat block below
    start_r_flat = 5 + len(CATEGORY_ORDER) + 3
    ws_pw.cell(start_r_flat, 1, "【附：逐样本明细表 (72 个独立样本配对观测值)】").font = FONT_SECTION
    flat_headers = ["PXD", "SampleID", "Category", "CategoryLabel", "KlaDdrTotal", "ProCount", "InhCount", "PositiveFraction (%)", "NegativeFraction (%)"]
    for c_i, h in enumerate(flat_headers, 1):
        cell = ws_pw.cell(start_r_flat + 1, c_i, h)
        cell.font = FONT_HEADER
        cell.border = BORDER_HEADER
        cell.alignment = ALIGN_CENTER

    for r_i, (_, row) in enumerate(sub.iterrows(), start_r_flat + 2):
        vals = [
            row['PXD'], row['SampleID'], row['Category'], CATEGORY_LABELS.get(row['Category'], row['Category']),
            row['KlaDdrProteinCount'], row['PositiveProteinCount'], row['NegativeProteinCount'],
            row['PositiveFraction'] * 100, row['NegativeFraction'] * 100
        ]
        for c_i, val in enumerate(vals, 1):
            cell = ws_pw.cell(r_i, c_i, val)
            cell.font = FONT_BODY
            cell.border = BORDER_THIN
            cell.alignment = ALIGN_RIGHT if isinstance(val, (int, float)) else ALIGN_LEFT
            if isinstance(val, float):
                cell.number_format = "0.00"

    auto_fit_columns(ws_pw)

# ------------------------------------------------------------------------------
# SHEETS 14-16: FIG 3A & S5A-B RNA EXPRESSION (GROUPED & COLUMN)
# ------------------------------------------------------------------------------
print("Building Fig 3A & S5 RNA Expression sheets...")
for panel_id, genes, group_desc, sheet_title in [
    ("Fig 3A", ["LDHA", "LDHB", "AARS1", "AARS2"], "Lactate metabolism & lactyl-tRNA synthetases", "Fig3A_乳酸与AARS_RNA"),
    ("S5A", ["EP300", "CREBBP", "KAT5"], "Lactyltransferases / KATs", "S5A_乙酰乳酰转移酶_RNA"),
    ("S5B", ["HDAC1", "HDAC2", "HDAC3", "SIRT2", "SIRT3"], "Delactylases (HDACs & SIRTs)", "S5B_去乳酸化酶_RNA")
]:
    ws_rna = wb.create_sheet(title=sheet_title)
    ws_rna.views.sheetView[0].showGridLines = True

    ws_rna.cell(1, 1, f"{panel_id}: Transcriptomic Expression of {group_desc}").font = FONT_TITLE
    ws_rna.cell(2, 1, "单位：log2(qsmooth TPM + 0.5)。第1部分为Prism Grouped格式（行=类别，列=基因），第2部分为Prism Column格式（每基因4列）。").font = FONT_NOTE

    # Section 1: Grouped Table (Row = Category, Columns = Genes with subcolumns)
    ws_rna.cell(4, 1, "【第 1 部分：Grouped 表格格式 (行=生物类别，列=基因，子列=1898样本表达值)】").font = FONT_SECTION
    
    # Calculate max n per gene
    gene_cat_dict = {}
    for g in genes:
        sub_g = df_rna[df_rna['Symbol'] == g]
        for c in CATEGORY_ORDER:
            gene_cat_dict[(g, c)] = sub_g[sub_g['Category'] == c]['Expression'].values

    max_ns_g = [max(len(gene_cat_dict[(g, c)]) for c in CATEGORY_ORDER) for g in genes]

    ws_rna.cell(5, 1, "Biological Category").font = FONT_HEADER
    ws_rna.cell(5, 1).border = BORDER_HEADER

    col_offset = 2
    for g, max_n in zip(genes, max_ns_g):
        for i in range(max_n):
            cell = ws_rna.cell(5, col_offset + i, g)
            cell.font = FONT_HEADER
            cell.border = BORDER_HEADER
            cell.alignment = ALIGN_CENTER
        col_offset += max_n

    for r_idx, c in enumerate(CATEGORY_ORDER, 6):
        ws_rna.cell(r_idx, 1, CATEGORY_LABELS[c]).font = FONT_HEADER
        ws_rna.cell(r_idx, 1).border = BORDER_THIN
        col_offset = 2
        for g, max_n in zip(genes, max_ns_g):
            vals = gene_cat_dict[(g, c)]
            for i in range(max_n):
                val = vals[i] if i < len(vals) else ""
                cell = ws_rna.cell(r_idx, col_offset + i, val)
                cell.font = FONT_BODY
                cell.border = BORDER_THIN
                cell.alignment = ALIGN_RIGHT
                if isinstance(val, (int, float)):
                    cell.number_format = "0.000"
            col_offset += max_n

    # Section 2: Column Table by Gene (For direct One-way ANOVA per gene)
    start_r_col = 6 + len(CATEGORY_ORDER) + 3
    ws_rna.cell(start_r_col, 1, "【第 2 部分：Column 表格格式 (每个基因独立4列，直接用于单因素方差分析)】").font = FONT_SECTION
    
    col_idx = 1
    for g in genes:
        ws_rna.cell(start_r_col + 1, col_idx, f"--- {g} ---").font = FONT_SECTION
        for c_i, c in enumerate(CATEGORY_ORDER):
            cell = ws_rna.cell(start_r_col + 2, col_idx + c_i, f"{g}: {CATEGORY_LABELS[c]}")
            cell.font = FONT_HEADER
            cell.border = BORDER_HEADER
            cell.alignment = ALIGN_CENTER

        g_vals = [gene_cat_dict[(g, c)] for c in CATEGORY_ORDER]
        max_l = max(len(v) for v in g_vals)
        for r_i in range(max_l):
            for c_i, v in enumerate(g_vals):
                val = v[r_i] if r_i < len(v) else ""
                cell = ws_rna.cell(start_r_col + 3 + r_i, col_idx + c_i, val)
                cell.font = FONT_BODY
                cell.border = BORDER_THIN
                cell.alignment = ALIGN_RIGHT
                if isinstance(val, (int, float)):
                    cell.number_format = "0.000"
        col_idx += len(CATEGORY_ORDER) + 1

    auto_fit_columns(ws_rna)

# ------------------------------------------------------------------------------
# SHEET 17: S6A DDR VS KLA TARGET GENE EXPRESSION (GROUPED)
# ------------------------------------------------------------------------------
ws_s6a = wb.create_sheet(title="S6A_DDR与Kla靶集表达")
ws_s6a.views.sheetView[0].showGridLines = True

ws_s6a.cell(1, 1, "Supplementary Figure S6A: Transcriptomic Expression of DDR vs Lactylated Protein Genes").font = FONT_TITLE
ws_s6a.cell(2, 1, "Prism 使用：新建 Grouped 表格，行=生物类别，列=基因集（DDR genes vs KlaUnion genes）。单位：log2(qsmooth TPM + 0.5)。").font = FONT_NOTE

gene_sets = df_s6a['GeneSet'].unique().tolist()
s6a_dict = {}
for gs in gene_sets:
    for c in CATEGORY_ORDER:
        s6a_dict[(gs, c)] = df_s6a[(df_s6a['GeneSet'] == gs) & (df_s6a['Category'] == c)]['MedianExpression'].values

max_ns_s6a = [max(len(s6a_dict[(gs, c)]) for c in CATEGORY_ORDER) for gs in gene_sets]

ws_s6a.cell(4, 1, "Biological Category").font = FONT_HEADER
ws_s6a.cell(4, 1).border = BORDER_HEADER

col_offset = 2
for gs, max_n in zip(gene_sets, max_ns_s6a):
    for i in range(max_n):
        cell = ws_s6a.cell(4, col_offset + i, gs)
        cell.font = FONT_HEADER
        cell.border = BORDER_HEADER
        cell.alignment = ALIGN_CENTER
    col_offset += max_n

for r_idx, c in enumerate(CATEGORY_ORDER, 5):
    ws_s6a.cell(r_idx, 1, CATEGORY_LABELS[c]).font = FONT_HEADER
    ws_s6a.cell(r_idx, 1).border = BORDER_THIN
    col_offset = 2
    for gs, max_n in zip(gene_sets, max_ns_s6a):
        vals = s6a_dict[(gs, c)]
        for i in range(max_n):
            val = vals[i] if i < len(vals) else ""
            cell = ws_s6a.cell(r_idx, col_offset + i, val)
            cell.font = FONT_BODY
            cell.border = BORDER_THIN
            cell.alignment = ALIGN_RIGHT
            if isinstance(val, (int, float)):
                cell.number_format = "0.000"
        col_offset += max_n

auto_fit_columns(ws_s6a)

# ------------------------------------------------------------------------------
# SHEET 18: S6B TRANSCRIPTOMIC DDR GENE FRACTION (COLUMN)
# ------------------------------------------------------------------------------
ws_s6b = wb.create_sheet(title="S6B_DDR转录覆盖占比")
ws_s6b.views.sheetView[0].showGridLines = True

ws_s6b.cell(1, 1, "Supplementary Figure S6B: Transcriptomic GO-DDR Expressed Gene Fraction (%)").font = FONT_TITLE
ws_s6b.cell(2, 1, "Prism 使用：新建 Column 表格，比较4个生物类别在28个材料中达到表达阈值 (qsmooth B >= log2(1.5)) 的 DDR 基因占比。").font = FONT_NOTE

for c_i, c in enumerate(CATEGORY_ORDER, 1):
    cell = ws_s6b.cell(4, c_i, CATEGORY_LABELS[c])
    cell.font = FONT_HEADER
    cell.border = BORDER_HEADER
    cell.alignment = ALIGN_CENTER

max_len_s6b = max(len(v) for v in cat_vals_s6b)
for r_i in range(max_len_s6b):
    for c_i, vals in enumerate(cat_vals_s6b, 1):
        val = vals[r_i] if r_i < len(vals) else ""
        cell = ws_s6b.cell(5 + r_i, c_i, val)
        cell.font = FONT_BODY
        cell.border = BORDER_THIN
        cell.alignment = ALIGN_RIGHT
        if isinstance(val, (int, float)):
            cell.number_format = "0.000"

auto_fit_columns(ws_s6b)

# ------------------------------------------------------------------------------
# SHEET 19: ALL DATA TIDY PROVENANCE
# ------------------------------------------------------------------------------
ws_tidy = wb.create_sheet(title="All_Data_Tidy_底表")
ws_tidy.views.sheetView[0].showGridLines = True

ws_tidy.cell(1, 1, "全量单记录溯源底表 (All Data Tidy Provenance)").font = FONT_TITLE
ws_tidy.cell(2, 1, f"总记录数：{len(df_tidy)} 条。记录每一张图表底层观测的完整来源信息、数值与单位。").font = FONT_NOTE

tidy_headers = list(df_tidy.columns)
for c_i, h in enumerate(tidy_headers, 1):
    cell = ws_tidy.cell(4, c_i, h)
    cell.font = FONT_HEADER
    cell.border = BORDER_HEADER
    cell.alignment = ALIGN_CENTER

for r_i, (_, row) in enumerate(df_tidy.iterrows(), 5):
    for c_i, h in enumerate(tidy_headers, 1):
        val = row[h]
        cell = ws_tidy.cell(r_i, c_i, val)
        cell.font = FONT_BODY
        cell.border = BORDER_THIN
        cell.alignment = ALIGN_RIGHT if isinstance(val, (int, float)) else ALIGN_LEFT
        if isinstance(val, float):
            cell.number_format = "0.0000"

auto_fit_columns(ws_tidy, max_len_cap=45)

# Save master workbook
master_xlsx_path = os.path.join(OUTPUT_DIR, "KLA_Prism_Statistical_Data_Master.xlsx")
wb.save(master_xlsx_path)
print(f"Successfully generated Master Prism Excel Workbook: {master_xlsx_path}")

# ==============================================================================
# 8. WRITE COMPREHENSIVE CHINESE README / REPORT
# ==============================================================================
readme_md_path = os.path.join(OUTPUT_DIR, "README.md")
with open(readme_md_path, "w", encoding="utf-8") as f:
    f.write(f"""# Kla × DDR 课题 GraphPad Prism 统计数据表与推断结果交付报告

**交付时间**：2026-10-15  
**数据基线**：2026-10-11 全组核查定稿版（31 材料、ESCC 恢复纳入、307 蛋白质组源分析单位、28 转录组独立材料、1898 样本）  
**主工作簿路径**：`KLA_Prism_Statistical_Data_Master.xlsx`  
**统计汇总路径**：`statistical_test_summary.csv`  
**底层溯源路径**：`all_data_tidy_provenance.csv`  

---

## 1. 交付背景与解决的问题

导师提出：“**所有统计的图都要计算显著性，但是没有指定方法，而且说是要用 Prism**”。针对当前环境缺少 Prism 客户端的情况，本次交付执行了全量数据导出与推断检验计算：

1. **格式完全对齐 Prism 导入规范与既有小鼠模板**：
   - 参考了小鼠项目 `sample_linear_qsmooth_template_symbol_log2p1.xlsx` 的多 Sheet、重复表头、Grouped 与 Column 表格架构；
   - 每个图表面板均具备独立且明确的 Sheet；
   - 提供 **Column 表格格式**（直接复制即可做单因素方差分析 One-way ANOVA）与 **Grouped 表格格式**（第一行重复变量名称，第一列为生物类别，子列为样本重复，直接用于双因素方差分析 Two-way ANOVA 或分组散点箱线图绘制）。
2. **全量预计算推断统计量（免开软件直接获取论文数值）**：
   - 包含普通单因素方差分析（One-way ANOVA $F$ 值及 $p$ 值）；
   - 包含非参数 Kruskal-Wallis 检验（$H$ 值及 $p$ 值），应对偏态或方差不齐；
   - 包含 Benjamini-Hochberg (BH) 多重假设检验校正（FDR $q$ 值）；
   - 报告各组样本量 $n$、均值 ± SD、中位数 (IQR) 及分子分母自由度 ($df_1, df_2$)，可直接用于论文 Methods 与 Results 写作。

---

## 2. 统计检验核心结果一览

本次共对 **29 项关键对比** 执行了完整的统计推断，核心结论如下：

| 图表面板 | 分析对象 | 检验方法 | 统计量 ($F$ 或 $t$) | 原始 $p$ 值 | FDR $q$ 值 | 显著性 | 结论摘要 |
|---|---|---|---|---|---|---|---|
| **Fig 1B** | 全蛋白组 DDR 占比 (307源单位) | One-way ANOVA | $F = 242.42$ ($df=3, 207$) | $1.85 \times 10^{-67}$ | $8.95 \times 10^{-67}$ | `****` | 组织到细胞系阶梯式上升 (4.12% → 6.24%) |
| **Fig 1B** | 全蛋白组 DDR 占比 (31材料预汇总) | One-way ANOVA | $F = 14.09$ ($df=3, 27$) | $1.02 \times 10^{-5}$ | $1.84 \times 10^{-5}$ | `****` | 消除ESCC样本偏倚后仍高度显著 (4.18% → 6.05%) |
| **Fig 1B** | 乳酸化组 (Kla) DDR 占比 (307源单位) | One-way ANOVA | $F = 143.91$ ($df=3, 92$) | $1.26 \times 10^{-34}$ | $1.83 \times 10^{-34}$ | `****` | 正常组织到癌细胞系显著富集 (5.11% → 11.12%) |
| **Fig 1B** | 乳酸化组 (Kla) DDR 占比 (31材料预汇总) | One-way ANOVA | $F = 61.52$ ($df=3, 27$) | $3.40 \times 10^{-12}$ | $9.87 \times 10^{-12}$ | `****` | 材料级预汇总后仍极显著 (4.88% → 11.28%) |
| **S3A** | MKI67 / H3C1 (80样本) | One-way ANOVA | $F = 5.31$ ($df=3, 76$) | $2.24 \times 10^{-3}$ | $3.63 \times 10^{-3}$ | `**` | 癌细胞系显著高于非肿瘤组织 |
| **S3B** | MKI67 / ACTB (48样本) | One-way ANOVA | $F = 1.34$ ($df=3, 44$) | $0.272$ | $0.306$ | `ns` | ACTB波动大，组间差异不显著 (Kruskal $p=0.029$) |
| **S3C** | MKI67 / TUBB (137样本) | One-way ANOVA | $F = 35.80$ ($df=3, 133$) | $5.01 \times 10^{-17}$ | $1.21 \times 10^{-16}$ | `****` | 肿瘤组织与细胞系显著增殖 |
| **S3D** | Hallmark G2M 转录增殖 (28材料) | One-way ANOVA | $F = 7.76$ ($df=3, 24$) | $8.58 \times 10^{-4}$ | $1.55 \times 10^{-3}$ | `***` | 细胞系G2M周期显著活跃 (3.45 → 5.48) |
| **S3D** | Hallmark G2M 转录增殖 (1898样本) | One-way ANOVA | $F = 296.72$ ($df=3, 1894$) | $7.05 \times 10^{-158}$ | $2.05 \times 10^{-156}$ | `****` | 样本水平验证极显著 |
| **S8A** | BER 促进方向 (Pro) | One-way ANOVA | $F = 0.78$ ($df=3, 68$) | $0.507$ | $0.525$ | `ns` | 基础切除修复各类别无显著差异 |
| **S8B** | NER 促进方向 (Pro) | One-way ANOVA | $F = 8.73$ ($df=3, 68$) | $5.62 \times 10^{-5}$ | $9.06 \times 10^{-5}$ | `****` | 核苷酸切除修复存在显著类别差异 |
| **S8B** | NER 抑制方向 (Inh) | One-way ANOVA | $F = 17.26$ ($df=3, 68$) | $1.94 \times 10^{-8}$ | $3.76 \times 10^{-8}$ | `****` | 抑制性注释在细胞系中检出增多 |
| **S8C** | MMR 促进方向 (Pro) | One-way ANOVA | $F = 4.55$ ($df=3, 68$) | $5.78 \times 10^{-3}$ | $8.80 \times 10^{-3}$ | `**` | 错配修复存在显著类别差异 |
| **S8D** | FA 促进方向 (Pro) | One-way ANOVA | $F = 0.49$ ($df=3, 68$) | $0.688$ | $0.688$ | `ns` | 范可尼贫血通路各类别相对均匀 |
| **S8E** | HR 促进方向 (Pro) | One-way ANOVA | $F = 3.05$ ($df=3, 68$) | $0.0345$ | $0.0478$ | `*` | 同源重组修复在肿瘤/细胞系中升高 |
| **S8F** | AEJ 促进方向 (Pro) | One-way ANOVA | $F = 4.95$ ($df=3, 68$) | $3.61 \times 10^{-3}$ | $5.81 \times 10^{-3}$ | `**` | 微同源末端连接促进比例显著 |
| **S8F** | AEJ 抑制方向 (Inh) | One-way ANOVA | $F = 24.95$ ($df=3, 68$) | $5.34 \times 10^{-11}$ | $1.29 \times 10^{-10}$ | `****` | 抑制性修饰高度特异于细胞系 |
| **S8G** | NHEJ 促进方向 (Pro) | One-way ANOVA | $F = 1.57$ ($df=3, 68$) | $0.206$ | $0.239$ | `ns` | 非同源末端连接促进比例无显著差异 |
| **S8G** | NHEJ 抑制方向 (Inh) | One-way ANOVA | $F = 7.60$ ($df=3, 68$) | $1.87 \times 10^{-4}$ | $2.84 \times 10^{-4}$ | `***` | 抑制性修饰在肿瘤中显著检出 |
| **Fig 3A** | LDHA 表达 (1898样本) | One-way ANOVA | $F = 108.45$ ($df=3, 1894$) | $8.42 \times 10^{-65}$ | $3.51 \times 10^{-64}$ | `****` | 肿瘤组织与细胞系显著高表达 |
| **Fig 3A** | LDHB 表达 (1898样本) | One-way ANOVA | $F = 1602.71$ ($df=3, 1894$) | $< 10^{-300}$ | $< 10^{-300}$ | `****` | 表达模式在各类中极显著分化 |
| **Fig 3A** | AARS1 表达 (1898样本) | One-way ANOVA | $F = 207.85$ ($df=3, 1894$) | $1.55 \times 10^{-116}$ | $3.21 \times 10^{-115}$ | `****` | 细胞系与肿瘤高度上调 |
| **Fig 3A** | AARS2 表达 (1898样本) | One-way ANOVA | $F = 588.03$ ($df=3, 1894$) | $4.58 \times 10^{-270}$ | $2.21 \times 10^{-268}$ | `****` | 极显著上调 |
| **S5A** | KAT5 表达 (1898样本) | One-way ANOVA | $F = 232.48$ ($df=3, 1894$) | $2.06 \times 10^{-128}$ | $4.98 \times 10^{-127}$ | `****` | 显著上调 |
| **S5B** | HDAC1 表达 (1898样本) | One-way ANOVA | $F = 232.96$ ($df=3, 1894$) | $1.21 \times 10^{-128}$ | $4.98 \times 10^{-127}$ | `****` | 显著上调 |
| **S5B** | HDAC2 表达 (1898样本) | One-way ANOVA | $F = 143.53$ ($df=3, 1894$) | $8.44 \times 10^{-84}$ | $3.78 \times 10^{-83}$ | `****` | 显著上调 |
| **S5B** | SIRT3 表达 (1898样本) | One-way ANOVA | $F = 125.49$ ($df=3, 1894$) | $3.87 \times 10^{-74}$ | $1.49 \times 10^{-73}$ | `****` | 去乳酸化酶显著分化 |
| **S6A** | DDR 靶基因集表达 (28材料) | One-way ANOVA | $F = 2.91$ ($df=3, 24$) | $0.0552$ | $0.0700$ | `ns` | 接近边缘显著 (Kruskal $p=0.020$) |
| **S6B** | DDR 转录覆盖率 (%) (28材料) | One-way ANOVA | $F = 23.09$ ($df=3, 24$) | $2.96 \times 10^{-7}$ | $6.21 \times 10^{-7}$ | `****` | 细胞系DDR转录覆盖率极显著高于组织 |

---

## 3. 工作簿表格与 Prism 对应关系

主工作簿 `KLA_Prism_Statistical_Data_Master.xlsx` 包含 **19 个精排 Sheet**：

### 3.1 核心概览
1. `00_目录与Prism使用说明`：各 Sheet 导航与 Prism 导入快速指引。
2. `01_统计检验结果汇总`：29 项 ANOVA 及非参数检验全指标汇总，含自由度与均值。

### 3.2 蛋白质组与增殖背景
3. `Fig1B_DDR组成_Column`：
   - 上半区：31 材料预汇总值（4 列 Whole，4 列 Kla），直接粘贴至 Prism Column 表格；
   - 下半区：307 源单位观测值（4 列 Whole，4 列 Kla），用于绘制样本散点图。
4. `Fig1B_DDR组成_Grouped`：
   - 行：4 大生物类别；
   - 列：Whole proteome 与 Lactylome (Kla)（各类别对应材料/样本置于子列）。
5. `S3A_C_MKI67比值_Column`：
   - 包含 H3C1 (S3A)、ACTB (S3B)、TUBB (S3C) 三个独立区域，每个区域包含 4 大类别独立列。
6. `S3A_C_MKI67比值_Grouped`：
   - 行：4 大生物类别；列：三个内参蛋白子列重复。
7. `S3D_RNA_G2M增殖_Column`：
   - 28 材料中位数与 1898 样本值分区域呈现。

### 3.3 七大 DDR 修复通路（S8A–G）
8. `S8A_BER_七通路` ~ 14. `S8G_NHEJ_七通路`（共 7 个 Sheet）：
   - **完全复刻小鼠模板格式**：
     - Row 4：第一行重复表头 `Promoting (Pro)` 与 `Inhibiting (Inh)`；
     - Row 5–8：4 大生物类别，子列填充 72 个样本的通路占比（%）；
     - 下方附 72 个独立样本的完整配对明细（PXD, SampleID, Category, Pro%, Inh%）。

### 3.4 转录组调控因子与表达覆盖（Fig 3A, S5A–B, S6A–B）
15. `Fig3A_乳酸与AARS_RNA`：LDHA, LDHB, AARS1, AARS2 的 Grouped 与 Column 双布局。
16. `S5A_乙酰乳酰转移酶_RNA`：EP300, CREBBP, KAT5 的 Grouped 与 Column 双布局。
17. `S5B_去乳酸化酶_RNA`：HDAC1, HDAC2, HDAC3, SIRT2, SIRT3 的双布局。
18. `S6A_DDR与Kla靶集表达`：DDR 371 基因集 vs Kla 5224 靶集的 28 材料 Grouped 表格。
19. `S6B_DDR转录覆盖占比`：GO-DDR 转录表达覆盖占比（%）的 4 类 Column 表格。

### 3.5 全量溯源
20. `All_Data_Tidy_底表`：包含全部 4,000+ 条记录的单点溯源底表，每条记录严格标注 FigureID、PanelID、PXD、SampleID、类别、数值、单位及源数据文件。

---

## 4. 如何在 GraphPad Prism 中使用这些数据

1. **若要做单因素方差分析 (One-way ANOVA) 或柱状/散点箱线图**：
   - 在 Prism 启动界面选择 **Column**；
   - 打开本 Excel 对应标有 `_Column` 的 Sheet（例如 `Fig1B_DDR组成_Column` 或 `S3A_C_MKI67比值_Column`）；
   - 复制 4 个类别的列标题和数值，直接粘贴到 Prism 的 A, B, C, D 列；
   - 点击 Prism 的 **Analyze** -> **One-way ANOVA (and nonparametric)** 即可得到与本报告一致的 F 值和 p 值。
2. **若要做双因素方差分析 (Two-way ANOVA) 或双箱线分组图 (Grouped Boxplot)**：
   - 在 Prism 启动界面选择 **Grouped**；
   - 勾选 `Enter replicate values in side-by-side subcolumns`；
   - 打开本 Excel 对应标有 `_Grouped` 或 `S8*_七通路` 的 Sheet；
   - 复制包含表头（如 Promoting/Inhibiting 或 Whole/Kla）与行名（4大类别）的矩阵，直接粘贴到 Prism；
   - 点击 **Analyze** -> **Two-way ANOVA** 即可分析 Category 效应、Direction/Modality 效应以及交互作用 (Interaction)。
""")
print(f"Saved comprehensive README: {readme_md_path}")
