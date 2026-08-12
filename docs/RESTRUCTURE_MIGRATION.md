# Publication restructure migration

Date: 2026-08-12

The former `reanalysis/` tree was promoted into a shallow publication layout:
`config/`, `python/`, `R/`, `workflow/`, `tests/`, `manuscript/`, `docs/`,
`results/`, and `work/`.

Current code retains only:

- the verified seven-PXD raw evidence parser;
- 40/37/33-group data preparation and reference audit;
- regulator intensity and DDR analyses;
- four-class Venn/Euler analysis;
- fixed 507-protein functional input preparation;
- BP-semantic UMAP, five-set UMAP/t-SNE/PCA, pathway-specific UMAP, and the
  final five-set separate linear pathway matrices.

Removed from the active tree:

- early three-group analytical endpoints;
- V1/V2 embedding variants and collision-only comparisons;
- circular pathway matrices;
- Cytoscape experiments;
- one-off workbook builders and download probes;
- incomplete all-R ports and temporary refactor planning files.

No raw PXD data were moved or deleted. Legacy source remains recoverable from
Git history at commit `8ba4052`; large local historical outputs remain under
the ignored `archive/` directory.
