# Lactate metabolism workflow

Started 2026-09-20, delivered 2026-09-21. Statistical methods in the existing release are unchanged.

## Reproduce from frozen configuration

From the repository root, run:

```bash
bash workflow/lactate_metabolism/run.sh outputs/new_lactate_reproduction
```

The output must not exist or must be empty. Input locations are resolved from `config/rna_current_release.csv`. R dependencies: data.table, ggplot2, ragg, jsonlite. Catalog regeneration additionally requires AnnotationDbi/org.Hs.eg.db, whose metadata is frozen. Cairo/Arial are used by the PDF device. `set.seed(25)` is set; no random inferential procedures are run.

Optionally set `KLA_QA_PYTHON` to a Python executable with Pillow and ensure Poppler `pdftoppm` is available to generate previews. Inspect previews manually; successful rendering is not visual approval.

## Annotation acquisition (explicit refresh, not part of routine reproduction)

The frozen download responses and parameters are versioned. To refresh, first create a separate configuration snapshot and update paths rather than overwriting an accepted catalog. The acquisition/build sequence used here was:

1. `fetch_go.py`: six terms, ontology descendants and all annotation pages.
2. `build_catalog.R`: initial stable-ID map and mechanism candidate resolution.
3. `fetch_go_mapping.py`: current UniProt Ensembl cross-references for GO protein products (ComplexPortal is not a UniProt accession).
4. `build_catalog.R`: current UniProt mapping preferred where available, installed org.Hs.eg.db fallback recorded.
5. `fetch_protein_evidence.py` then `extract_evidence.py`: current reviewed protein functions for mechanistic candidates.
6. `analyze.R`, `validate_catalog.R`, `report.R`, `render_qa.py`.

Candidate names are resolved only when constructing the frozen mapping; all expression joins use Ensembl IDs. GO annotations and curated mechanism members are separate. No missing-value imputation or module flux score is used.

`freeze_release.py` packages methods/evidence and writes source/output SHA256 inventories after numerical and visual checks; it does not grant visual approval. Numerical tables are reproducible; device metadata can change PDF bytes between runs.
