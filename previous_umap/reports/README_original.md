# GO repair/damage KLA-DDR workflow

This folder is an isolated rerun of the lactylation-related DNA damage/repair
classification workflow. It does not depend on, overwrite, or move the older
root-level pipeline outputs.

## Inputs

- `data/article_tables/PXD*.csv`: lactylated protein/gene tables extracted from
  the source studies. `PXD053474` is handled with the same `GeneSymbol` +
  `UniProtKB` logic as the other source tables.
- `data/GO-repair+damage(human).tsv`: QuickGO/GO table downloaded for human DNA
  repair and DNA damage terms.
- `data/legacy_layout/umap_8type_data.csv`: old UMAP coordinates used only to
  keep the figure layout comparable with the previous `umap_8type.png`.
- `data/genesets.tsv`: MSigDB GMT-like metadata file. This is kept for optional
  sensitivity analysis, but the default workflow is GO-only.

## Classification logic

Default evidence source is GO-only. Evidence codes are not filtered, because the
current aim is a broad annotation map rather than a high-confidence functional
claim. `NOT` qualifiers are excluded.

Categories are intentionally conservative:

- `HR`: explicit homologous recombination / single-strand annealing terms
- `NHEJ`: explicit non-homologous end joining terms
- `BER`: explicit base excision or single-strand break repair terms
- `NER`: explicit nucleotide excision / UV-damage excision terms
- `MMR`: explicit mismatch repair terms
- `TLS`: explicit translesion synthesis or DNA damage tolerance terms
- `DRR`: direct reversal / photoreactive repair terms
- `CP`: checkpoint, DDR signaling, or generic DNA repair / DSB repair terms
- `Other`: DNA damage-related repair terms that should not be forced into the
  seven mechanisms, such as interstrand cross-link repair, DNA alkylation repair,
  protein-DNA cross-link repair, mitochondrial DNA repair, telomere maintenance
  in response to DNA damage, and repair-dependent chromatin remodeling.

The rule table written to `results/classification_rules.csv` records the curated
GO overrides used by the script.

The UMAP layout is computed from the GO-term evidence matrix, not from the
collapsed category proportions. This gives the embedding enough feature
resolution to separate genes that otherwise collapse into identical 8/9-class
profiles.

## Run

```bash
python3 run_workflow.py
Rscript plot_umap.R
```

Optional MSigDB sensitivity run:

```bash
python3 run_workflow.py --include-msigdb
Rscript plot_umap.R
```

The MSigDB option is not recommended for the main figure unless explicitly
explained, because broad pathway gene sets increase multi-mechanism assignments.

## Outputs

- `results/kla_ddr_unique_go_repair_damage.csv`: deduplicated candidate genes
- `results/gene_repair_evidence_long.csv`: long evidence table
- `results/gene_repair_category_matrix.csv`: gene-by-category matrix
- `results/gene_primary_category.csv`: primary category per gene
- `results/umap_pie_data.csv`: plot-ready table with old coordinates where
  available
- `results/genes_excluded_from_legacy_umap_plot.csv`: genes without old
  coordinates, kept out of the legacy-layout UMAP figure
- `figures/umap_repair_pies.png` and `.pdf`: pie UMAP
- `figures/umap_repair_primary.png` and `.pdf`: primary-category UMAP
