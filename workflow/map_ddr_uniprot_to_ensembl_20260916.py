#!/usr/bin/env python3
"""Map the DDR panel's UniProt BaseAccessions to Ensembl gene identifiers.

The DDR annotation the project holds is keyed by UniProt BaseAccession (protein side: the GO
annotation table, the venn membership files, the S4 pathway ranking). The RNA matrices are
keyed by Ensembl gene ID. This builds the bridge between the two.

Source: UniProt's own HUMAN_9606_idmapping.dat.gz, read locally. The online idmapping API was
tried first and abandoned - a few hundred identifiers sit queued there for tens of minutes,
while the file answers the same question in seconds and does not depend on the service.

Two routes are taken so they can be cross-checked rather than trusted singly:

  route 1  UniProt -> Ensembl gene                                  (what the overlay needs)
  route 2  UniProt -> Entrez GeneID -> Ensembl gene via NCBI gene2ensembl

UniProt -> Ensembl_PRO was tried as the second route first and does not work here: this
release of the human idmapping file carries Ensembl and GeneID but no Ensembl_PRO for many of
these entries, so an empty protein mapping means "not carried", not "disagrees".

Disagreements between the routes are reported, not silently resolved.

Usage: map_ddr_uniprot_to_ensembl_20260916.py <out_dir> <HUMAN_9606_idmapping.dat.gz> <gene2ensembl.gz>
"""
import csv
import gzip
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# DB labels as they appear in the idmapping file
WANTED = {"Ensembl": "gene", "Ensembl_PRO": "protein", "GeneID": "entrez"}


def read_venn(path):
    with open(path, newline="", encoding="utf-8-sig") as fh:
        return list(csv.DictReader(fh))


def read_entrez_to_ensembl(path):
    """Entrez GeneID -> Ensembl gene, human only; genes reached by several GeneIDs are kept
    against each of them because the join direction here is GeneID -> gene."""
    out = {}
    with gzip.open(path, "rt", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 3 or p[0] != "9606" or not p[2]:
                continue
            out.setdefault(p[1], set()).add(p[2].split(".")[0])
    return out


def main(out_dir, idmap_path, gene2ensembl_path):
    os.makedirs(out_dir, exist_ok=True)

    kla = read_venn(os.path.join(ROOT, "data/publication_input/venn_kla_ddr.csv"))
    ref = read_venn(os.path.join(ROOT, "data/publication_input/venn_reference_ddr.csv"))
    panel = {}
    for row in kla:
        panel.setdefault(row["BaseAccession"], set()).add("KlaDdr")
    for row in ref:
        panel.setdefault(row["BaseAccession"], set()).add("ReferenceDdr")

    # the GO annotation table carries the wider DDR universe the panels are drawn from
    with open(os.path.join(ROOT, "data/publication_input/human_ddr_go_annotations.tsv"),
              newline="", encoding="utf-8-sig") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            acc = (row.get("GENE PRODUCT ID") or "").strip()
            if acc:
                panel.setdefault(acc, set()).add("DdrGoUniverse")

    # lactylation writer / eraser / reader enzymes. These are not DDR proteins, so they are
    # not in the DDR sets, but the regulator expression panel needs them mapped the same way.
    reg_path = os.path.join(ROOT, "data/publication_input", "regulator_kla_percentiles_31.csv")
    with open(reg_path, newline="", encoding="utf-8-sig") as fh:
        for row in csv.DictReader(fh):
            acc = (row.get("RegulatorBaseAccession") or "").strip()
            if acc:
                panel.setdefault(acc, set()).add("LactylationRegulator")

    # every lactylated protein, for the "genes behind the lactylated proteins" question
    with open(os.path.join(ROOT, "data/publication_input/venn_all_kla.csv"),
              newline="", encoding="utf-8-sig") as fh:
        for row in csv.DictReader(fh):
            acc = (row.get("BaseAccession") or "").strip()
            if acc:
                panel.setdefault(acc, set()).add("KlaUnion")

    ids = set(panel)
    n_go = sum(1 for a in panel if "DdrGoUniverse" in panel[a])
    n_reg = sum(1 for a in panel if "LactylationRegulator" in panel[a])
    print(f"accessions to map: {len(ids)} "
          f"(Kla n DDR {len(kla)}, reference DDR {len(ref)}, GO universe {n_go}, "
          f"lactylation regulators {n_reg}, Kla union {sum(1 for a in panel if 'KlaUnion' in panel[a])})")

    hits = {acc: {"gene": set(), "protein": set(), "entrez": set()} for acc in ids}
    with gzip.open(idmap_path, "rt", encoding="utf-8") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 3:
                continue
            acc, db, val = parts
            kind = WANTED.get(db)
            if kind and acc in hits:
                hits[acc][kind].add(val.split(".")[0])
    print(f"accessions with an Ensembl gene: {sum(1 for a in ids if hits[a]['gene'])}")

    # the project's cached ENSEMBLPROT -> UniProt mapping, inverted
    cached = {}
    with open(os.path.join(ROOT, "data/publication_input/reference_protein_membership_31.csv"),
              newline="", encoding="utf-8-sig") as fh:
        for r in csv.DictReader(fh):
            if r["IdentifierType"] != "ENSEMBLPROT":
                continue
            ensp = r["SourceProteinID"].split(".")[0]
            for acc in r["MappedBaseAccessions"].split(";"):
                if acc:
                    cached.setdefault(acc, set()).add(ensp)

    print("reading gene2ensembl for the second route ...")
    entrez2ensg = read_entrez_to_ensembl(gene2ensembl_path)
    print(f"  human GeneIDs with an Ensembl gene: {len(entrez2ensg)}")

    rows = []
    for acc in sorted(ids):
        genes = sorted(hits[acc]["gene"])
        entrez = sorted(hits[acc]["entrez"])
        via_entrez = set()
        for e in entrez:
            via_entrez |= entrez2ensg.get(e, set())
        via_entrez = sorted(via_entrez)
        cached_p = sorted(cached.get(acc, set()))

        if genes and via_entrez:
            status = "agree" if set(genes) & set(via_entrez) else "disagree"
        elif genes:
            status = "route2_unavailable_no_entrez"
        else:
            status = "no_ensembl_gene"

        rows.append({
            "BaseAccession": acc,
            "PanelMembership": ";".join(sorted(panel[acc])),
            "EnsemblGeneIDs": ";".join(genes),
            "NGeneIDs": len(genes),
            "AmbiguousGeneMapping": len(genes) > 1,
            "EntrezGeneIDs": ";".join(entrez),
            "EnsemblGeneViaEntrez": ";".join(via_entrez),
            "CachedEnsemblProteinIDs": ";".join(cached_p),
            "CrossCheck": status,
        })

    out_path = os.path.join(out_dir, "ddr_uniprot_to_ensembl.tsv")
    with open(out_path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    usable = [r for r in rows if r["EnsemblGeneIDs"]]
    unambiguous = [r for r in usable if not r["AmbiguousGeneMapping"]]
    agree = [r for r in rows if r["CrossCheck"] == "agree"]
    disagree = [r for r in rows if r["CrossCheck"] == "disagree"]
    r2_na = [r for r in rows if r["CrossCheck"] == "route2_unavailable_no_entrez"]
    print(f"mapped to >=1 Ensembl gene : {len(usable)}/{len(rows)}")
    print(f"unambiguous 1:1 mapping    : {len(unambiguous)}")
    print(f"many-to-one (several ENSG) : {len(usable) - len(unambiguous)}")
    print(f"routes agree               : {len(agree)}")
    print(f"routes disagree            : {len(disagree)}")
    print(f"route 2 unavailable        : {len(r2_na)}")
    for r in disagree[:10]:
        print(f"   DISAGREE {r['BaseAccession']}: direct={r['EnsemblGeneIDs']} "
              f"via_entrez={r['EnsemblGeneViaEntrez']}")

    with open(os.path.join(out_dir, "mapping_summary.csv"), "w") as fh:
        fh.write("metric,value\n")
        for k, v in [("panel_accessions", len(rows)),
                     ("mapped_to_ensembl_gene", len(usable)),
                     ("unambiguous_gene_mapping", len(unambiguous)),
                     ("multi_gene_mapping", len(usable) - len(unambiguous)),
                     ("no_ensembl_gene", len(rows) - len(usable)),
                     ("routes_agree", len(agree)),
                     ("routes_disagree", len(disagree)),
                     ("route2_unavailable", len(r2_na))]:
            fh.write(f"{k},{v}\n")
    print("MAPPING_DONE", out_path)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])
