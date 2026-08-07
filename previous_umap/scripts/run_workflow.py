#!/usr/bin/env python3
from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
import re

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "input"
RESULTS = ROOT / "results"
PROJECT_ARCHIVE = ROOT.parent / "archive" / "previous_umap"

CATEGORY_CODES = ["HR", "NHEJ", "BER", "NER", "MMR", "TLS", "DRR", "CP", "Other"]
MECHANISM_CODES = ["HR", "NHEJ", "BER", "NER", "MMR", "TLS", "DRR"]
CATEGORY_LABELS = {
    "HR": "Homologous recombination repair",
    "NHEJ": "Non-homologous end joining",
    "BER": "Base excision / single-strand break repair",
    "NER": "Nucleotide excision repair",
    "MMR": "Mismatch repair",
    "TLS": "Translesion synthesis / DNA damage tolerance",
    "DRR": "Direct reversal repair",
    "CP": "DDR signaling / generic DNA repair",
    "Other": "Other DNA damage-related repair terms",
}

EVIDENCE_CODE_WEIGHT = {
    "EXP": 3.0,
    "IDA": 3.0,
    "IPI": 3.0,
    "IMP": 3.0,
    "IGI": 3.0,
    "IEP": 3.0,
    "HTP": 3.0,
    "HDA": 3.0,
    "HMP": 3.0,
    "HGI": 3.0,
    "HEP": 3.0,
    "TAS": 2.0,
    "NAS": 2.0,
    "IC": 2.0,
    "ISS": 1.5,
    "ISO": 1.5,
    "ISA": 1.5,
    "ISM": 1.5,
    "IBA": 1.5,
    "IBD": 1.5,
    "IKR": 1.5,
    "IRD": 1.5,
    "RCA": 1.5,
    "IEA": 1.0,
}

# Curated GO-to-category overrides. Rules are intentionally conservative:
# only explicit pathway terms enter concrete repair mechanisms; generic DNA
# repair and DDR signaling enter CP. Damage-type, chromatin, telomere,
# mitochondrial, and cross-link terms enter Other unless the term itself
# specifies one of the seven mechanisms.
GO_TERM_TO_CATEGORY = {
    "GO:0000724": "HR",  # double-strand break repair via homologous recombination
    "GO:0000725": "HR",  # recombinational repair
    "GO:0000730": "HR",  # DNA recombinase assembly
    "GO:0010792": "HR",  # DSB processing involved in single-strand annealing
    "GO:0045002": "HR",  # DSB repair via single-strand annealing
    "GO:0045003": "HR",  # synthesis-dependent strand annealing
    "GO:0006303": "NHEJ",  # DSB repair via nonhomologous end joining
    "GO:0031848": "NHEJ",  # protection from NHEJ at telomere
    "GO:0097680": "NHEJ",  # classical NHEJ
    "GO:0097681": "NHEJ",  # alternative NHEJ
    "GO:0000012": "BER",  # single strand break repair
    "GO:0006284": "BER",  # base-excision repair
    "GO:0006285": "BER",  # base-excision repair, AP site formation
    "GO:0006287": "BER",  # base-excision repair, gap-filling
    "GO:0097510": "BER",  # base-excision repair, AP site formation via deaminated base removal
    "GO:0097698": "BER",  # telomere maintenance via base-excision repair
    "GO:0006283": "NER",  # transcription-coupled NER
    "GO:0006289": "NER",  # nucleotide-excision repair
    "GO:0006297": "NER",  # nucleotide-excision repair, DNA gap filling
    "GO:0000720": "NER",  # pyrimidine dimer repair by NER
    "GO:0070914": "NER",  # UV-damage excision repair
    "GO:1904161": "NER",  # DNA damage recognition, global genome NER
    "GO:1901255": "NER",  # nucleotide-excision repair, DNA incision
    "GO:0000710": "MMR",  # meiotic mismatch repair
    "GO:0006298": "MMR",  # mismatch repair
    "GO:0070716": "MMR",  # mismatch repair involved in meiotic recombination
    "GO:0006301": "TLS",  # DNA damage tolerance
    "GO:0019985": "TLS",  # translesion synthesis
    "GO:0042276": "TLS",  # error-prone TLS
    "GO:0070987": "TLS",  # error-free TLS
    "GO:0000719": "DRR",  # photoreactive repair
    "GO:0000077": "CP",  # DNA damage checkpoint signaling
    "GO:0000729": "CP",  # DNA double-strand break processing
    "GO:0000731": "Other",  # DNA synthesis involved in DNA repair
    "GO:0006281": "CP",  # DNA repair
    "GO:0006302": "CP",  # double-strand break repair
    "GO:0006307": "Other",  # DNA alkylation repair, damage-type term
    "GO:0006974": "CP",  # DNA damage response
    "GO:0008630": "CP",  # apoptotic signaling in response to DNA damage
    "GO:0030330": "CP",  # p53-mediated DDR signaling
    "GO:0031571": "CP",  # mitotic G1 DNA damage checkpoint signaling
    "GO:0031573": "CP",  # mitotic intra-S DNA damage checkpoint signaling
    "GO:0036297": "Other",  # interstrand cross-link repair, composite repair program
    "GO:0036298": "Other",  # recombinational interstrand cross-link repair
    "GO:0042770": "CP",  # signal transduction in response to DNA damage
    "GO:0042771": "CP",  # p53 apoptotic signaling in response to DNA damage
    "GO:0042772": "CP",  # DNA damage response, signal transduction resulting in transcription
    "GO:0043247": "Other",  # telomere maintenance in response to DNA damage
    "GO:0043504": "Other",  # mitochondrial DNA repair
    "GO:0044773": "CP",  # mitotic DNA damage checkpoint signaling
    "GO:0106300": "Other",  # protein-DNA covalent cross-linking repair
    "GO:0140861": "Other",  # DNA repair-dependent chromatin remodeling
    "GO:1990414": "Other",  # replication-born DSB repair via sister chromatid exchange
}

GO_TERM_SPECIFICITY = {
    # Broad HR/recombination terms are less specific than explicit DSB-HR terms.
    "GO:0000725": 1.0,
    # Exact mechanism terms used as primary tie-break evidence.
    "GO:0000724": 3.0,
    "GO:0006303": 3.0,
    "GO:0097680": 3.0,
    "GO:0097681": 3.0,
    "GO:0000012": 3.0,
    "GO:0006284": 3.0,
    "GO:0006289": 3.0,
    "GO:0006298": 3.0,
    "GO:0019985": 3.0,
    "GO:0006301": 3.0,
    "GO:0000719": 3.0,
}

GO_NAME_RULES = [
    ("HR", re.compile(r"\bhomologous recombination\b|\bhomology[- ]directed repair\b|\bsynthesis-dependent strand annealing\b|\bsingle-strand annealing\b", re.I)),
    ("NHEJ", re.compile(r"\bnon[- ]?homologous end joining\b|\bend joining\b", re.I)),
    ("BER", re.compile(r"\bbase[- ]excision repair\b|\bsingle strand break repair\b", re.I)),
    ("NER", re.compile(r"\bnucleotide[- ]excision repair\b|\buv[- ]damage excision repair\b|\bpyrimidine dimer repair\b", re.I)),
    ("MMR", re.compile(r"\bmismatch repair\b", re.I)),
    ("TLS", re.compile(r"\btranslesion synthesis\b|\bdna damage tolerance\b|\bpost[- ]?replication repair\b", re.I)),
    ("DRR", re.compile(r"\bphotoreactive repair\b|\bdirect reversal\b|\bdna direct repair\b", re.I)),
    ("Other", re.compile(r"\bcross-link repair\b|\balkylation repair\b|\bmitochondrial dna repair\b|\btelomere maintenance in response to dna damage\b|\bdna repair-dependent chromatin remodeling\b|\bdna synthesis involved in dna repair\b", re.I)),
    ("CP", re.compile(r"\bdna damage response\b|\bdna repair\b|\bdouble-strand break repair\b|\bcheckpoint signaling\b|\bsignal transduction in response to dna damage\b", re.I)),
]

ALLOWED_MSIGDB_COLLECTIONS = {
    "C5:GO:BP",
    "C2:CP:REACTOME",
    "C2:CP:KEGG_LEGACY",
    "C2:CP:KEGG_MEDICUS",
    "C2:CP:WIKIPATHWAYS",
    "ARCHIVED:C5_GO:BP",
    "ARCHIVED:C2_CP:REACTOME",
}

MSIGDB_RULES = [
    ("HR", re.compile(r"HOMOLOG(?:OUS)?_RECOMBINATION|HOMOLOGY_DIRECTED_REPAIR|SINGLE_STRAND_ANNEALING|SYNTHESIS_DEPENDENT_STRAND_ANNEALING", re.I)),
    ("NHEJ", re.compile(r"NON_?HOMOLOGOUS_END_JOINING|NONHOMOLOGOUS_END_JOINING", re.I)),
    ("BER", re.compile(r"BASE_EXCISION_REPAIR|SINGLE_STRAND_BREAK_REPAIR", re.I)),
    ("NER", re.compile(r"NUCLEOTIDE_EXCISION_REPAIR|UV_DAMAGE_EXCISION_REPAIR|PYRIMIDINE_DIMER_REPAIR", re.I)),
    ("MMR", re.compile(r"MISMATCH_REPAIR", re.I)),
    ("TLS", re.compile(r"TRANSLESION_SYNTHESIS|DNA_DAMAGE_TOLERANCE|POSTREPLICATION_REPAIR|POST_REPLICATION_REPAIR", re.I)),
    ("DRR", re.compile(r"DNA_DIRECT_REPAIR|DIRECT_REVERSAL|PHOTOREACTIVE_REPAIR", re.I)),
    ("Other", re.compile(r"CROSS_LINK_REPAIR|ALKYLATION_REPAIR|MITOCHONDRIAL_DNA_REPAIR|TELOMERE_MAINTENANCE_IN_RESPONSE_TO_DNA_DAMAGE|DNA_REPAIR_DEPENDENT_CHROMATIN_REMODELING", re.I)),
    ("CP", re.compile(r"(^|_)DNA_REPAIR($|_)|(^|_)DOUBLE_STRAND_BREAK_REPAIR($|_)|DNA_DAMAGE_RESPONSE|DNA_DAMAGE_CHECKPOINT|DNA_INTEGRITY_CHECKPOINT|CHECKPOINT_SIGNALING|SIGNAL_TRANSDUCTION_IN_RESPONSE_TO_DNA_DAMAGE", re.I)),
]


@dataclass
class ClassificationResult:
    evidence: pd.DataFrame
    matrix: pd.DataFrame
    primary: pd.DataFrame
    pie: pd.DataFrame
    counts: pd.DataFrame


def base_acc(value: object) -> str:
    text = "" if pd.isna(value) else str(value).strip()
    text = re.sub(r"^(sp|tr)\|", "", text)
    if "|" in text:
        text = text.split("|", 1)[0]
    text = re.sub(r"^.*?:", "", text)
    text = re.sub(r"-\d+$", "", text)
    return text.strip()


def first_non_empty(values: pd.Series) -> str:
    for value in values:
        if not pd.isna(value) and str(value).strip():
            return str(value).strip()
    return ""


def read_article_table(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    symbol_col = next((c for c in ("GeneSymbol", "Symbol") if c in df.columns), None)
    uniprot_col = next((c for c in ("UniProtKB", "UniProt") if c in df.columns), None)
    if symbol_col is None or uniprot_col is None:
        raise ValueError(f"{path} needs GeneSymbol/Symbol and UniProtKB/UniProt columns")

    out = pd.DataFrame(
        {
            "Source": path.stem,
            "Symbol": df[symbol_col].astype(str).str.strip().str.upper(),
            "UniProtKB": df[uniprot_col].astype(str).str.strip(),
            "ProteinName": df["ProteinName"].astype(str) if "ProteinName" in df.columns else "",
            "KlaSites": df["KlaSites"].astype(str) if "KlaSites" in df.columns else "",
        }
    )
    out["BaseAcc"] = out["UniProtKB"].map(base_acc)
    out = out[(out["Symbol"] != "") & (out["BaseAcc"] != "")]
    return out.drop_duplicates(["Source", "Symbol", "BaseAcc", "KlaSites"])


def build_kla_gene_table(article_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    files = sorted(article_dir.glob("PXD*.csv"))
    if not files:
        raise FileNotFoundError(f"No PXD*.csv files found in {article_dir}")

    provenance = pd.concat([read_article_table(path) for path in files], ignore_index=True)
    provenance = provenance.sort_values(["Source", "Symbol", "BaseAcc"]).reset_index(drop=True)

    grouped = provenance.groupby("BaseAcc", as_index=False)
    unique_rows = []
    for base_acc_value, group in grouped:
        symbols = sorted(set(group["Symbol"]))
        sources = sorted(set(group["Source"]))
        sites = sorted({str(x).strip() for x in group["KlaSites"] if not pd.isna(x) and str(x).strip() not in {"", "nan", "NA"}})
        unique_rows.append(
            {
                "Symbol": symbols[0],
                "UniProtKB": first_non_empty(group["UniProtKB"]) or base_acc_value,
                "AliasSymbols": ";".join(symbols),
                "Sources": ";".join(sources),
                "KlaSites": ";".join(sites),
                "Source_Count": len(sources),
            }
        )

    genes = pd.DataFrame(unique_rows).sort_values("Symbol").reset_index(drop=True)
    return provenance, genes


def classify_go_name(go_name: object, go_term: object | None = None) -> str | None:
    term = "" if go_term is None or pd.isna(go_term) else str(go_term).strip()
    if term in GO_TERM_TO_CATEGORY:
        return GO_TERM_TO_CATEGORY[term]
    name = "" if pd.isna(go_name) else str(go_name).strip()
    for category, pattern in GO_NAME_RULES:
        if pattern.search(name):
            return category
    return None


def classify_msigdb_name(name: object, collection: object | None = None) -> str | None:
    pathway = "" if pd.isna(name) else str(name).strip().upper()
    coll = "" if collection is None or pd.isna(collection) else str(collection).strip()
    if coll and coll not in ALLOWED_MSIGDB_COLLECTIONS:
        return None
    for category, pattern in MSIGDB_RULES:
        if pattern.search(pathway):
            return category
    return None


def _empty_evidence() -> pd.DataFrame:
    return pd.DataFrame(
        columns=["Symbol", "BaseAcc", "Category", "Evidence_ID", "Evidence_Name", "Evidence_Code", "Reference", "Source"]
    )


def load_go_evidence(go_path: Path, genes: pd.DataFrame) -> pd.DataFrame:
    go = pd.read_csv(go_path, sep="\t")
    required = {"GENE PRODUCT ID", "SYMBOL", "QUALIFIER", "GO TERM", "GO NAME", "GO EVIDENCE CODE", "TAXON ID"}
    missing = required.difference(go.columns)
    if missing:
        raise ValueError(f"GO file missing columns: {sorted(missing)}")

    gene_symbols = set(genes["Symbol"])
    gene_accs = set(genes["UniProtKB"].map(base_acc))
    work = go.copy()
    work["SymbolRaw"] = work["SYMBOL"].astype(str).str.strip().str.upper()
    work["BaseAcc"] = work["GENE PRODUCT ID"].map(base_acc)
    work["Qualifier"] = work["QUALIFIER"].fillna("").astype(str)
    work["Taxon"] = pd.to_numeric(work["TAXON ID"], errors="coerce")
    work["Category"] = [classify_go_name(name, term) for name, term in zip(work["GO NAME"], work["GO TERM"])]

    work = work[
        work["Category"].notna()
        & (work["BaseAcc"].isin(gene_accs) | work["SymbolRaw"].isin(gene_symbols))
        & (work["Taxon"].isna() | (work["Taxon"] == 9606))
        & ~work["Qualifier"].str.contains("NOT", case=False, na=False)
    ].copy()
    if work.empty:
        return _empty_evidence()

    acc_lookup = genes.assign(BaseAcc=genes["UniProtKB"].map(base_acc)).drop_duplicates("BaseAcc").set_index("BaseAcc")["Symbol"].to_dict()
    valid_symbols = set(genes["Symbol"])
    work["Symbol"] = work["BaseAcc"].map(acc_lookup).fillna(work["SymbolRaw"].where(work["SymbolRaw"].isin(valid_symbols)))

    return work[
        ["Symbol", "BaseAcc", "Category", "GO TERM", "GO NAME", "GO EVIDENCE CODE", "REFERENCE"]
    ].rename(
        columns={
            "GO TERM": "Evidence_ID",
            "GO NAME": "Evidence_Name",
            "GO EVIDENCE CODE": "Evidence_Code",
            "REFERENCE": "Reference",
        }
    ).assign(Source="GO")


def parse_msigdb_records(msigdb_path: Path) -> list[dict[str, str]]:
    records: list[dict[str, str]] = []
    current: dict[str, str] = {}
    with open(msigdb_path, encoding="utf-8", errors="ignore") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if not line or "\t" not in line:
                continue
            key, value = line.split("\t", 1)
            if key == "STANDARD_NAME":
                if current:
                    records.append(current)
                current = {"STANDARD_NAME": value}
            else:
                current[key] = value
    if current:
        records.append(current)
    return records


def load_msigdb_evidence(msigdb_path: Path, genes: pd.DataFrame) -> pd.DataFrame:
    if not msigdb_path.exists():
        return _empty_evidence()
    rows = []
    gene_symbols = set(genes["Symbol"])
    for record in parse_msigdb_records(msigdb_path):
        category = classify_msigdb_name(record.get("STANDARD_NAME", ""), record.get("COLLECTION", ""))
        if category is None:
            continue
        record_genes = {x.strip().upper() for x in record.get("GENE_SYMBOLS", "").split(",") if x.strip()}
        for symbol in sorted(gene_symbols.intersection(record_genes)):
            rows.append(
                {
                    "Symbol": symbol,
                    "BaseAcc": "",
                    "Category": category,
                    "Evidence_ID": record.get("STANDARD_NAME", ""),
                    "Evidence_Name": record.get("STANDARD_NAME", ""),
                    "Evidence_Code": record.get("COLLECTION", ""),
                    "Reference": record.get("PMID", ""),
                    "Source": "MSigDB",
                }
            )
    return pd.DataFrame(rows).drop_duplicates() if rows else _empty_evidence()


def build_matrix(genes: pd.DataFrame, evidence: pd.DataFrame) -> pd.DataFrame:
    matrix = genes[["Symbol", "UniProtKB"]].copy()
    matrix["BaseAcc"] = matrix["UniProtKB"].map(base_acc)
    for code in CATEGORY_CODES:
        matrix[code] = matrix["Symbol"].isin(set(evidence.loc[evidence["Category"] == code, "Symbol"])).astype(int)
    matrix["Mechanism_Count"] = matrix[MECHANISM_CODES].sum(axis=1)
    matrix["Any_Category_Count"] = matrix[CATEGORY_CODES].sum(axis=1)
    return matrix


def primary_category(counts: dict[str, int]) -> str:
    mechanisms = [code for code in MECHANISM_CODES if counts.get(code, 0) > 0]
    if mechanisms:
        return sorted(mechanisms, key=lambda c: (-counts.get(c, 0), MECHANISM_CODES.index(c)))[0]
    if counts.get("CP", 0) > 0:
        return "CP"
    if counts.get("Other", 0) > 0:
        return "Other"
    return "Unclassified"


def _evidence_weight(code: object) -> float:
    text = "" if pd.isna(code) else str(code).strip().upper()
    if text in ALLOWED_MSIGDB_COLLECTIONS:
        return 1.2
    return EVIDENCE_CODE_WEIGHT.get(text, 1.0)


def _term_specificity(evidence_id: object) -> float:
    text = "" if pd.isna(evidence_id) else str(evidence_id).strip()
    return GO_TERM_SPECIFICITY.get(text, 2.0)


def primary_mechanism_from_evidence(gene_evidence: pd.DataFrame) -> tuple[str | None, bool, str, str]:
    mechanism_evidence = gene_evidence[gene_evidence["Category"].isin(MECHANISM_CODES)].copy()
    if mechanism_evidence.empty:
        return None, False, "", ""

    rows = []
    for category, sub in mechanism_evidence.groupby("Category", observed=True):
        rows.append(
            {
                "Category": category,
                "Count": len(sub),
                "Specificity": sum(_term_specificity(x) for x in sub["Evidence_ID"]),
                "EvidenceWeight": sum(_evidence_weight(x) for x in sub["Evidence_Code"]),
                "Order": MECHANISM_CODES.index(category),
            }
        )
    scores = pd.DataFrame(rows)
    scores = scores.sort_values(["Count", "Specificity", "EvidenceWeight", "Order"], ascending=[False, False, False, True])
    top = scores.iloc[0]
    tied = scores[
        (scores["Count"] == top["Count"])
        & (scores["Specificity"] == top["Specificity"])
        & (scores["EvidenceWeight"] == top["EvidenceWeight"])
    ]
    tie_categories = ";".join(tied.sort_values("Order")["Category"].tolist()) if len(tied) > 1 else ""
    score_text = ";".join(
        f"{row.Category}:n={int(row.Count)},specificity={row.Specificity:g},evidence={row.EvidenceWeight:g}"
        for row in scores.itertuples(index=False)
    )
    return str(top["Category"]), len(tied) > 1, tie_categories, score_text


def proportions(counts: dict[str, int]) -> dict[str, float]:
    mechanism_total = sum(counts.get(code, 0) for code in MECHANISM_CODES)
    if mechanism_total > 0:
        return {code: (counts.get(code, 0) / mechanism_total if code in MECHANISM_CODES else 0.0) for code in CATEGORY_CODES}
    nonmechanism_total = counts.get("CP", 0) + counts.get("Other", 0)
    if nonmechanism_total > 0:
        return {
            code: (counts.get(code, 0) / nonmechanism_total if code in {"CP", "Other"} else 0.0)
            for code in CATEGORY_CODES
        }
    return {code: 0.0 for code in CATEGORY_CODES}


def classify_genes(genes: pd.DataFrame, go_path: Path, msigdb_path: Path | None = None) -> ClassificationResult:
    go_evidence = load_go_evidence(go_path, genes)
    evidence = go_evidence
    if msigdb_path is not None:
        evidence = pd.concat([evidence, load_msigdb_evidence(msigdb_path, genes)], ignore_index=True)

    evidence = evidence.drop_duplicates(["Symbol", "Category", "Evidence_ID", "Source"]).sort_values(
        ["Symbol", "Category", "Source", "Evidence_ID"]
    )
    matrix = build_matrix(genes, evidence)

    evidence_counts = evidence.groupby(["Symbol", "Category"], observed=True).size().unstack(fill_value=0)
    evidence_counts = evidence_counts.reindex(columns=CATEGORY_CODES, fill_value=0)
    primary_rows = []
    pie_rows = []
    for _, row in matrix.iterrows():
        symbol = row["Symbol"]
        counts = {code: int(evidence_counts.loc[symbol, code]) if symbol in evidence_counts.index else 0 for code in CATEGORY_CODES}
        gene_evidence = evidence[evidence["Symbol"] == symbol]
        primary_mechanism, is_primary_tie, tie_categories, score_text = primary_mechanism_from_evidence(gene_evidence)
        primary = primary_mechanism or primary_category(counts)
        mechanisms = [code for code in MECHANISM_CODES if counts[code] > 0]
        active = mechanisms if mechanisms else ([primary] if primary in {"CP", "Other"} else [])
        primary_rows.append(
            {
                "Symbol": symbol,
                "UniProtKB": row["UniProtKB"],
                "Primary_Category": primary,
                "Primary_Label": CATEGORY_LABELS.get(primary, "Unclassified"),
                "Mechanism_Count": len(mechanisms),
                "Is_Multimechanism": len(mechanisms) > 1,
                "Primary_Is_Tie": is_primary_tie,
                "Primary_Tie_Categories": tie_categories,
                "Primary_Score_Detail": score_text,
                "Active_Categories": ";".join(active),
                "Evidence_Count": sum(counts.values()),
            }
        )
        pie_rows.append({"Symbol": symbol, **proportions(counts)})

    counts_df = matrix[CATEGORY_CODES].sum().rename_axis("Category").reset_index(name="Gene_Count")
    counts_df["Label"] = counts_df["Category"].map(CATEGORY_LABELS)
    return ClassificationResult(
        evidence=evidence.reset_index(drop=True),
        matrix=matrix.reset_index(drop=True),
        primary=pd.DataFrame(primary_rows),
        pie=pd.DataFrame(pie_rows),
        counts=counts_df,
    )


def build_umap_features(genes: pd.DataFrame, evidence: pd.DataFrame) -> pd.DataFrame:
    term_matrix = (
        evidence.assign(value=1)
        .pivot_table(index="Symbol", columns="Evidence_ID", values="value", aggfunc="max", fill_value=0)
        .astype(float)
    )
    term_matrix = term_matrix.reindex(genes["Symbol"], fill_value=0)
    if term_matrix.shape[1] == 0:
        term_matrix["no_evidence"] = 0.0

    category_counts = (
        evidence.groupby(["Symbol", "Category"], observed=True)
        .size()
        .unstack(fill_value=0)
        .reindex(index=genes["Symbol"], columns=CATEGORY_CODES, fill_value=0)
        .astype(float)
    )
    category_counts = np.log1p(category_counts) * 1.5
    category_counts.columns = [f"CategoryCount_{col}" for col in category_counts.columns]

    features = pd.concat([term_matrix, category_counts], axis=1)
    features = features.reset_index().rename(columns={"index": "Symbol"})
    return features


def compute_embedding(pie: pd.DataFrame, features_df: pd.DataFrame) -> tuple[pd.DataFrame, str]:
    features = features_df.drop(columns=["Symbol"]).to_numpy(dtype=float)
    nonzero_sd = features.std(axis=0)
    nonzero_sd[nonzero_sd == 0] = 1.0
    features = features / nonzero_sd
    if features.shape[0] < 3:
        coords = np.zeros((features.shape[0], 2))
        method = "degenerate"
    else:
        try:
            import umap  # type: ignore

            reducer = umap.UMAP(
                n_neighbors=min(30, max(2, features.shape[0] - 1)),
                min_dist=0.8,
                spread=3.0,
                metric="correlation",
                random_state=25,
                n_jobs=1,
            )
            coords = reducer.fit_transform(features)
            method = "UMAP on GO-term plus category-count evidence matrix"
        except Exception as exc:
            centered = features - features.mean(axis=0, keepdims=True)
            _, _, vt = np.linalg.svd(centered, full_matrices=False)
            coords = centered @ vt[:2].T
            if coords.shape[1] == 1:
                coords = np.column_stack([coords[:, 0], np.zeros(coords.shape[0])])
            method = f"PCA fallback on GO-term plus category-count evidence matrix: {exc}"

    out = pie.copy()
    out["Strict_UMAP_1"] = coords[:, 0]
    out["Strict_UMAP_2"] = coords[:, 1]
    out["UMAP_1"] = coords[:, 0]
    out["UMAP_2"] = coords[:, 1]
    out["Embedding_Method"] = method
    out["Layout_Source"] = method
    return out, method


def add_legacy_coordinates(plot_data: pd.DataFrame, legacy_path: Path) -> pd.DataFrame:
    if not legacy_path.exists():
        return plot_data
    legacy = pd.read_csv(legacy_path)[["Symbol", "UMAP_1", "UMAP_2"]].rename(
        columns={"Symbol": "Legacy_Symbol", "UMAP_1": "Legacy_UMAP_1", "UMAP_2": "Legacy_UMAP_2"}
    )
    work = plot_data.copy()
    # H2AFX is the current HGNC symbol; the old figure used H2AX.
    alias_to_legacy = {"H2AFX": "H2AX"}
    work["Legacy_Symbol"] = work["Symbol"].replace(alias_to_legacy)
    work = work.merge(legacy, on="Legacy_Symbol", how="left")
    has_legacy = work["Legacy_UMAP_1"].notna() & work["Legacy_UMAP_2"].notna()
    work["UMAP_1"] = work["UMAP_1"].astype(float)
    work["UMAP_2"] = work["UMAP_2"].astype(float)
    work["Has_Legacy_Coordinates"] = has_legacy
    return work


def write_rule_table(path: Path) -> None:
    rule_labels = {
        **{code: "Concrete repair mechanism" for code in MECHANISM_CODES},
        "CP": "Generic DNA repair / DDR signaling",
        "Other": "Conservative bucket for DNA damage-related repair terms not forced into a mechanism",
    }
    rows = [
        {"Category": code, "Label": CATEGORY_LABELS[code], "Rule": rule_labels[code]}
        for code in CATEGORY_CODES
    ]
    for term, category in sorted(GO_TERM_TO_CATEGORY.items()):
        rows.append({"Category": category, "Label": CATEGORY_LABELS[category], "Rule": f"{term} curated GO override"})
    pd.DataFrame(rows).to_csv(path, index=False, encoding="utf-8-sig")


def write_symbol_conflicts(provenance: pd.DataFrame, path: Path) -> None:
    conflicts = (
        provenance.groupby("BaseAcc")
        .agg(
            Symbols=("Symbol", lambda x: ";".join(sorted(set(x)))),
            Source_Count=("Source", "nunique"),
            Sources=("Source", lambda x: ";".join(sorted(set(x)))),
        )
        .reset_index()
    )
    conflicts["Symbol_Count"] = conflicts["Symbols"].str.split(";").map(len)
    conflicts = conflicts.loc[conflicts["Symbol_Count"] > 1].sort_values(["Symbol_Count", "BaseAcc"], ascending=[False, True])
    conflicts.to_csv(path, index=False, encoding="utf-8-sig")


def main() -> None:
    parser = argparse.ArgumentParser(description="Build KLA-DDR gene classification from local GO repair/damage annotations.")
    parser.add_argument("--include-msigdb", action="store_true", help="Add MSigDB pathway evidence as sensitivity analysis. Default is GO-only.")
    args = parser.parse_args()

    RESULTS.mkdir(exist_ok=True)
    provenance, genes = build_kla_gene_table(DATA / "article_tables")
    result = classify_genes(
        genes,
        DATA / "GO-repair+damage(human).tsv",
        PROJECT_ARCHIVE / "msigdb_sensitivity" / "genesets.tsv" if args.include_msigdb else None,
    )
    umap_features = build_umap_features(genes, result.evidence)
    pie_with_coords, embedding_method = compute_embedding(result.pie, umap_features)
    umap_data = pie_with_coords.merge(result.primary, on="Symbol", how="left").merge(
        genes[["Symbol", "AliasSymbols", "Sources", "KlaSites", "Source_Count"]], on="Symbol", how="left"
    )
    umap_data = add_legacy_coordinates(
        umap_data, PROJECT_ARCHIVE / "legacy_layout" / "umap_8type_data.csv"
    )
    umap_data["Has_Legacy_Coordinates"] = umap_data.get("Has_Legacy_Coordinates", False)

    provenance.to_csv(RESULTS / "kla_ddr_source_provenance.csv", index=False, encoding="utf-8-sig")
    genes.to_csv(RESULTS / "kla_ddr_unique_go_repair_damage.csv", index=False, encoding="utf-8-sig")
    result.evidence.to_csv(RESULTS / "gene_repair_evidence_long.csv", index=False, encoding="utf-8-sig")
    result.matrix.to_csv(RESULTS / "gene_repair_category_matrix.csv", index=False, encoding="utf-8-sig")
    result.primary.to_csv(RESULTS / "gene_primary_category.csv", index=False, encoding="utf-8-sig")
    result.counts.to_csv(RESULTS / "repair_category_counts.csv", index=False, encoding="utf-8-sig")
    umap_features.to_csv(RESULTS / "umap_go_term_feature_matrix.csv", index=False, encoding="utf-8-sig")
    umap_data.to_csv(RESULTS / "umap_pie_data.csv", index=False, encoding="utf-8-sig")
    umap_data.loc[~umap_data.get("Has_Legacy_Coordinates", False)].to_csv(
        RESULTS / "genes_without_legacy_umap_coordinates.csv", index=False, encoding="utf-8-sig"
    )
    write_rule_table(RESULTS / "classification_rules.csv")
    write_symbol_conflicts(provenance, RESULTS / "gene_symbol_conflicts.csv")

    source_counts = result.evidence.groupby(["Source", "Category"], observed=True).size().reset_index(name="Evidence_Rows")
    source_counts.to_csv(RESULTS / "repair_category_evidence_source_counts.csv", index=False, encoding="utf-8-sig")
    evidence_code_counts = result.evidence.groupby(["Source", "Evidence_Code"], observed=True).size().reset_index(name="Evidence_Rows")
    evidence_code_counts.to_csv(RESULTS / "evidence_code_counts.csv", index=False, encoding="utf-8-sig")

    print(f"Candidate genes: {len(genes)}")
    print(f"Evidence rows: {len(result.evidence)}")
    print(f"Genes with legacy UMAP coordinates: {umap_data.get('Has_Legacy_Coordinates', pd.Series(False, index=umap_data.index)).sum()}")
    print(f"Embedding method: {embedding_method}")
    print(f"Main evidence source: {'GO + MSigDB' if args.include_msigdb else 'GO only'}")
    print(f"Results written to: {RESULTS}")


if __name__ == "__main__":
    main()
