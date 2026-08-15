# Final separate 4+1 linear pathway-state figures for the current 30-group analysis

## Final figure choice

The linear layout is selected as the final display. The earlier circular layout
is retained as a versioned comparison but is not part of this final 4+1 output.
The linear state matrix and pathway summary are exported as separate figures.
In addition, the four biological categories are not combined into one
multi-panel figure: each category is exported independently, and the all-protein
set is exported as a fifth independent figure.

## Five protein sets

1. tumor tissues, n = 178.
2. non-tumor tissues, n = 183.
3. cancer cell lines, n = 381.
4. normal cell lines, n = 292.
5. All current Kla-intersection-DDR proteins, n = 399.

The four biological categories reuse the fixed membership from the existing
30-group analysis and can overlap. Every set is keyed by isoform-stripped
UniProt `BaseAccession`.

## Independent ordering within each panel

Proteins are independently sorted within each of the five panels by the signed
weighted score:

`BER*1 + NER*2 + MMR*3 + FA*4 + HR*5 + AEJ*6 + NHEJ*7`

`BaseAccession` is the tie-breaker. The score is used only to determine protein
order. It does not determine color, cell width, pathway state, or summary
proportion.

## State encoding

- `+1` promoting: solid pathway color.
- `-1` suppressing: dark charcoal `#2F3437`.
- `0` unassigned: light gray `#F1F3F5`.

The same pathway colors and state colors are used across all five panels and
both figure types.

## Fully separate outputs

For each set key
`normal_tissue`, `cancer_tissue`, `cancer_cells`, `normal_cells`, and
`all_kla_ddr`, the independent linear matrix is:

- `kla_ddr_linear_pathway_matrix_<set>_en.{png,pdf,svg}`
- `kla_ddr_linear_pathway_matrix_<set>_zh.{png,pdf,svg}`

The corresponding independent summary is:

- `kla_ddr_pathway_state_summary_<set>_en.{png,pdf,svg}`
- `kla_ddr_pathway_state_summary_<set>_zh.{png,pdf,svg}`

All five summary panels use the same percentage axis. Counts and within-set
percentages are displayed for `-1` and `+1`, and the number of `0` proteins is
listed separately. These summaries are descriptive and are not statistical
enrichment tests.

This yields 10 figure stems per language and 60 files across PNG, PDF, and SVG.
PNG files are exported at 600 dpi. PDF and SVG files are retained for lossless
enlargement and editing.
