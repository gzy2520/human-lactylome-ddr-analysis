# Previous UMAP methodology

## Scope

The retained primary coordinates are an actual unsupervised UMAP embedding.
They are not a supervised projection, a manually arranged layout, or the
archived legacy coordinates. The current results record the embedding method as
`UMAP on GO-term plus category-count evidence matrix`.

## Evidence matrix

`scripts/run_workflow.py` builds one feature row per Kla-associated gene from
the retained article tables and human GO repair/damage annotations. The matrix
combines:

- binary columns for individual GO evidence terms; and
- log-transformed evidence-row counts for HR, NHEJ, BER, NER, MMR, TLS, DRR,
  CP, and Other, multiplied by 1.5.

Each feature is divided by its standard deviation before embedding. Category
labels are derived from documented GO-term rules, but no category is supplied
to UMAP as a prediction target.

## UMAP parameters

The workflow uses `umap-learn` with correlation distance, `min_dist=0.8`,
`spread=3.0`, `random_state=25`, one job, and up to 30 neighbors (bounded by the
number of genes). Pie slices and primary-category colors are overlays added
after the coordinates are computed; they do not place the points.

The code contains a deterministic PCA fallback for environments where UMAP
cannot run. That fallback was not used for the retained primary results. Legacy
coordinates are loaded only into separately named `Legacy_UMAP_1` and
`Legacy_UMAP_2` columns for historical comparison and do not replace the
primary `UMAP_1` and `UMAP_2` values.

## Reproduction paths

Run the retained workflow from `previous_umap`:

```bash
python3 scripts/run_workflow.py
Rscript scripts/plot_umap.R
```

Active inputs are under `input/`; generated tables and figures go to
`results/` and `figures/`. Optional MSigDB sensitivity data and legacy-layout
coordinates remain under `../archive/previous_umap/` and are referenced there
by the workflow.
