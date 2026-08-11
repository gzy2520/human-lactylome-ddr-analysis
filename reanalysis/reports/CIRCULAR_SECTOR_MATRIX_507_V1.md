# Circular seven-sector matrix for 507 Kla-DDR proteins

## Scope

- The figure contains all 507 fixed Kla-intersection-DDR proteins.
- Protein identity is the isoform-stripped UniProt `BaseAccession`.
- The seven sectors follow the fixed order
  `BER`, `NER`, `MMR`, `FA`, `HR`, `AEJ`, `NHEJ`.
- The two auxiliary non-pathway categories are excluded.
- Proteins are ordered from the inner to the outer radius by ascending signed
  weighted score, with `BaseAccession` as the tie-breaker.

## Visual encoding

- Each concentric arc is one protein.
- Each angular sector is one pathway.
- `+1` (promoting) is shown as the solid pathway color.
- `-1` (suppressing) is shown in dark charcoal (`#2F3437`).
- `0` (unassigned) is shown in very light gray.
- Existing pathway colors are reused for continuity with the UMAP, t-SNE, PCA,
  and Cytoscape figures.

The weighted score is used only as the protein-ordering key. It does not
determine sector color, arc width, pathway state, or pathway-state proportion.

The pathway sectors are deliberately separated by narrow white gaps. Individual
protein boundaries are not drawn for zero or positive cells because 507
concentric outlines would create a dark moire pattern.

## What density can and cannot show

The circular matrix can reveal descriptive concentration: for example, a dense
block of solid color means that many proteins in a particular score-ranked
interval have `+1` for that pathway. This can be reported as a higher local
proportion or concentration.

It is not, by itself, a statistical enrichment analysis. The rank is calculated
from the same seven pathway assignments and their coefficients, so the ordering
and the colored pattern are mathematically coupled. Calling a region "enriched"
would require a pre-defined region, an independent background, and an explicit
test such as Fisher's exact or a hypergeometric test with multiple-testing
correction.

For quantitative follow-up, the output table
`pathway_state_density_by_25_rank_bin.csv` reports `-1/0/+1` counts and
fractions in consecutive 25-protein rank bins. These bins are descriptive and
are not presented as hypothesis tests.

## Outputs

Circle-focused versions:

- `kla_ddr_circular_sector_matrix_507_en.{png,pdf,svg}`
- `kla_ddr_circular_sector_matrix_507_zh.{png,pdf,svg}`

Versions with a quantitative pathway-state summary panel:

- `kla_ddr_circular_sector_matrix_with_summary_507_en.{png,pdf,svg}`
- `kla_ddr_circular_sector_matrix_with_summary_507_zh.{png,pdf,svg}`

Linear unrolled versions:

- `kla_ddr_linear_pathway_matrix_507_en.{png,pdf,svg}`
- `kla_ddr_linear_pathway_matrix_507_zh.{png,pdf,svg}`

Linear versions with the quantitative pathway-state summary panel:

- `kla_ddr_linear_pathway_matrix_with_summary_507_en.{png,pdf,svg}`
- `kla_ddr_linear_pathway_matrix_with_summary_507_zh.{png,pdf,svg}`

The PNG files are exported at 600 dpi. PDF and SVG files are retained for
lossless enlargement and editing.
