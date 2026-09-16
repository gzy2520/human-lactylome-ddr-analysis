#!/usr/bin/env python3
"""Lay the DDR annotation over the smoothed 31-group reference profile.

The two sides are keyed differently - the DDR annotation by UniProt BaseAccession, the RNA
profile by Ensembl gene - so they are joined through the mapping built by
map_ddr_uniprot_to_ensembl_20260916.py, never through a gene symbol.

Three query tables come out, smallest to largest:

  Kla n DDR        the 401 proteins that are both lactylated and DDR. The paper's main line.
  reference DDR    the 836 DDR proteins detected in the reference proteome.
  DDR GO universe  every accession carrying a DDR-related GO annotation (2,719).

Values are the qsmooth-smoothed log2(TPM + 0.5) of the material-class reference profile, one
column per group row (31 columns; the rows that share a reference matrix are identical and are
flagged rather than hidden).

Usage: overlay_ddr_panel_20260916.py <qsmooth_dir> <mapping_dir> <out_dir>
"""
import csv
import gzip
import os
import re
import sys
import xml.etree.ElementTree as ET
import zipfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
M = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
NS = {"m": M[1:-1]}
PATHWAYS = ["BER", "NER", "MMR", "FA", "HR", "AEJ", "NHEJ"]


def read_s4_pathways(path):
    """BaseAccession -> {pathway: state} from every sheet of the S4 workbook.

    A protein can appear on several sheets; states are identical for a given accession by
    construction, so the first non-empty state wins.
    """
    out = {}
    with zipfile.ZipFile(path) as zf:
        if "xl/sharedStrings.xml" in zf.namelist():
            root = ET.fromstring(zf.read("xl/sharedStrings.xml"))
            strings = ["".join(t.text or "" for t in si.iter(M + "t"))
                       for si in root.findall("m:si", NS)]
        else:
            strings = []

        def col_index(ref):
            n = 0
            for ch in "".join(c for c in ref if c.isalpha()):
                n = n * 26 + (ord(ch) - 64)
            return n - 1

        sheets = sorted(n for n in zf.namelist() if re.match(r"xl/worksheets/sheet\d+\.xml$", n))
        for sheet in sheets:
            root = ET.fromstring(zf.read(sheet))
            header = None
            for row in root.iter(M + "row"):
                cells = {}
                for c in row.findall("m:c", NS):
                    v = c.find("m:v", NS)
                    if v is None or v.text is None:
                        continue
                    cells[col_index(c.get("r") or "")] = (
                        strings[int(v.text)] if c.get("t") == "s" else v.text)
                if not cells:
                    continue
                vals = [cells.get(i, "") for i in range(max(cells) + 1)]
                if header is None:
                    header = vals
                    continue
                if "BaseAccession" not in header:
                    continue
                rec = dict(zip(header, vals))
                acc = rec.get("BaseAccession", "").strip()
                if not acc:
                    continue
                out.setdefault(acc, {p: rec.get(p, "") for p in PATHWAYS})
    return out


def read_go_terms(path):
    """BaseAccession -> set of (GO term, GO name)."""
    out = {}
    with open(path, newline="", encoding="utf-8-sig") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            acc = (row.get("GENE PRODUCT ID") or "").strip()
            if not acc:
                continue
            out.setdefault(acc, set()).add(((row.get("GO TERM") or "").strip(),
                                            (row.get("GO NAME") or "").strip()))
    return out


def read_mapping(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return {r["BaseAccession"]: r for r in csv.DictReader(fh, delimiter="\t")}


def read_panel(expr_dir, name, accessions, mapping, columns, genes_index):
    """One row per accession that resolves to a gene present in the profile."""
    rows, unmapped_gene, missing_matrix = [], [], []
    for acc in sorted(accessions):
        m = mapping.get(acc)
        gene = ""
        if m and m["EnsemblGeneIDs"]:
            for g in m["EnsemblGeneIDs"].split(";"):
                if g in genes_index:
                    gene = g
                    break
            if not gene:
                unmapped_gene.append(acc)
                continue
        else:
            unmapped_gene.append(acc)
            continue
        rows.append((acc, gene, genes_index[gene]))
    return rows, unmapped_gene


def main(qsmooth_dir, mapping_dir, out_dir):
    os.makedirs(out_dir, exist_ok=True)

    # smoothed reference profile, one column per group row
    prof_path = os.path.join(qsmooth_dir, "matrices", "qsmooth_A_collapsed_log2tpm.tsv.gz")
    expansion = list(csv.DictReader(open(os.path.join(qsmooth_dir, "group_expansion_31.csv"),
                                         encoding="utf-8")))
    with gzip.open(prof_path, "rt", encoding="utf-8") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        per_ref = header[1:]
        ref_col = {r: i for i, r in enumerate(per_ref)}
        columns = [r["GroupID"] for r in expansion]
        col_index = [ref_col[r["ReferenceKey"]] for r in expansion]
        shared = [r["SharedReference"] == "True" for r in expansion]
        genes_index = {}
        for line in fh:
            p = line.rstrip("\n").split("\t")
            genes_index[p[0]] = [p[1 + i] for i in col_index]
    print(f"reference profile: {len(genes_index)} genes x {len(expansion)} group rows")

    mapping = read_mapping(os.path.join(mapping_dir, "ddr_uniprot_to_ensembl.tsv"))
    s4 = read_s4_pathways(os.path.join(ROOT, "data/publication_input",
                                       "Supplementary_Table_S4_Pathway_Protein_Ranking.xlsx"))
    go = read_go_terms(os.path.join(ROOT, "data/publication_input",
                                    "human_ddr_go_annotations.tsv"))

    def panel_accs(fname):
        with open(os.path.join(ROOT, "data/publication_input", fname),
                  newline="", encoding="utf-8-sig") as fh:
            return set(r["BaseAccession"] for r in csv.DictReader(fh))

    kla = panel_accs("venn_kla_ddr.csv")
    ref = panel_accs("venn_reference_ddr.csv")
    universe = set(go)

    summary = []
    for label, accs in [("kla_ddr", kla), ("reference_ddr", ref), ("ddr_go_universe", universe)]:
        rows, unmapped = read_panel(qsmooth_dir, label, accs, mapping, columns, genes_index)
        out_path = os.path.join(out_dir, f"{label}_expression_31groups.tsv.gz")
        with gzip.open(out_path, "wt", encoding="utf-8") as fh:
            fh.write("\t".join(["BaseAccession", "EnsemblGeneID"] + columns) + "\n")
            for acc, gene, vals in rows:
                fh.write("\t".join([acc, gene] + vals) + "\n")
        summary.append((label, len(accs), len(rows), len(unmapped)))
        print(f"{label:<18} accessions={len(accs):>5} with gene in profile={len(rows):>5} "
              f"unmapped/absent={len(unmapped):>5}")

    # annotation companion for the main panel
    with open(os.path.join(out_dir, "kla_ddr_annotation.csv"), "w", newline="",
              encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["BaseAccession", "EnsemblGeneID", "CrossCheck"] + PATHWAYS
                   + ["SignedScore", "GO terms"])
        for acc in sorted(kla):
            m = mapping.get(acc, {})
            gene = m.get("EnsemblGeneIDs", "")
            states = s4.get(acc, {p: "" for p in PATHWAYS})
            score = ""
            terms = ";".join(f"{t} {n}" for t, n in sorted(go.get(acc, set())))
            w.writerow([acc, gene, m.get("CrossCheck", "")] +
                       [states.get(p, "") for p in PATHWAYS] + [score, terms])

    with open(os.path.join(out_dir, "overlay_summary.csv"), "w") as fh:
        fh.write("panel,accessions,with_gene_in_profile,unmapped_or_absent\n")
        for label, n, kept, unmapped in summary:
            fh.write(f"{label},{n},{kept},{unmapped}\n")
        fh.write(f"profile_genes,{len(genes_index)},\n")
        fh.write(f"group_rows,{len(columns)},\n")
        fh.write(f"shared_reference_rows,{sum(shared)},\n")
    print("OVERLAY_DONE")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])
