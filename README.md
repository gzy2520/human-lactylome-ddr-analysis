# Integrated human lysine lactylation proteomics and DNA-damage response

This repository contains the current publication workflow for a cross-study
human lysine lactylation (Kla) proteomics analysis. The active teacher-approved
scope contains 30 paired sample groups and 399 unique Kla-intersecting DDR
proteins.

The current publication outputs are:

- a DDR-fraction bar chart;
- two regulator percentile heatmaps;
- four exact-membership, four-category Venn diagrams;
- separate 4+1 direct-GO-term pathway matrices and summaries.

The Venn geometry is schematic and not area proportional. All displayed
numbers come from exact membership and region tables. The pathway matrices no
longer use the former manual ±1 workbook: direct BP, CC, and MF GO terms are
mapped to BER, NER, MMR, FA, HR, NHEJ, AEJ, or `Others`.

An optional supervisor-score preview is retained separately and does not
replace the direct-GO-term publication result:

```bash
Rscript workflow/run_pipeline.R revised_score_preview
```

It reads the revised 2026-08-16 workbook, filters its 507 scored proteins to
the current 399-protein membership, and writes to a dedicated preview
directory.

## Reproduce

```bash
python3 -m pip install -r requirements.txt
Rscript workflow/install_r_dependencies.R
Rscript workflow/run_pipeline.R selected_figures
Rscript workflow/record_environment.R .
Rscript workflow/build_manifest.R .
```

The analytical key is the isoform-stripped human UniProt `BaseAccession`.
Gene symbols and protein names are display/audit annotations only. Randomized
steps use seed 25.

Raw and processed PXD files remain under `data/` locally and are not tracked by
Git. Generated `results/` and historical `archive/` contents are also outside
Git. Archived 33-group, 507-protein, UMAP, t-SNE, PCA, Cytoscape, circular
matrix, and manual-score materials are historical only and must not be used to
describe the current analysis.

Current entry points:

- `PROJECT_INDEX.md`: authoritative project map and scope;
- `NEW_CHAT_PROJECT_PROMPT.md`: copy-paste project handoff;
- `manuscript/methods/METHODS_EN.md` and `METHODS_ZH.md`: active Methods;
- `docs/DATA_PROVENANCE.md`: current provenance and inclusion logic;
- `docs/GO_TERM_PATHWAY_SCORING_30GROUPS.md`: pathway-scoring contract;
- `docs/REVISED_SCORE_WORKBOOK_PREVIEW_20260816.md`: separate exploratory
  supervisor-score preview;
- `tests/validate_publication_contract.R`: executable result contract.
