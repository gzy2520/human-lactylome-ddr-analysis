# Difficulty-weighted DDR score rank grids for the fixed 33-group analysis

## Scope and score definition

- The analysis is restricted to the fixed 507 Kla-intersection-DDR proteins,
  keyed by isoform-stripped UniProt `BaseAccession`.
- The revised source workbook is
  `data/identifier/260810乳酸化DDR基因评分表.xlsx`.
- The `评分表` sheet contains 507 active protein rows and seven trailing rows
  without `BaseAccession`; the blank rows are excluded.
- The signed weighted score is reproduced exactly as:

  `BER*1 + NER*2 + MMR*3 + FA*4 + HR*5 + AEJ*6 + NHEJ*7`

- `Chromatin Interaction` and
  `Others (Transcription, RNA processing and proteostasis)` are not used.
- The second figure row is `abs(signed weighted score)`. It is not the
  alternative component-wise calculation `sum(weight * abs(pathway score))`.

All 507 cached Excel `score` values exactly match the independently recalculated
seven-pathway formula.

## The 2x5 layout

Rows:

1. Original signed weighted score.
2. Absolute value of the signed weighted score.

Columns:

1. Normal/non-tumor tissues, n = 183.
2. Normal/non-tumor cells, n = 471.
3. Cancer tissues, n = 178.
4. Cancer cells, n = 383.
5. All Kla-DDR proteins, n = 507.

Category membership is the fixed four-category membership used by the existing
33-group Venn and embedding analyses. Categories overlap; each panel is the set
of proteins detected in that category.

## Ordering recommendation

The primary version ranks proteins independently within every panel by the
plotted score, with `BaseAccession` as the tie-breaker. This produces a rank
distribution curve and makes the score range, zero plateau and upper tail
directly readable. The horizontal axis should therefore be described as
`protein rank`, not as a continuous biological trajectory.

An alphabetical comparison is retained because it was specifically requested
by the teacher. It uses `GeneSymbol` only for display ordering, with
`BaseAccession` as the stable identity and tie-breaker. Alphabetical ordering
does not encode a biological or quantitative relationship and produces a
jagged connected line, so it is less suitable as the primary analytical
figure.

## Outputs

Primary ascending-score grids:

- `kla_ddr_weighted_score_ranked_2x5_en.{png,pdf,svg}`
- `kla_ddr_weighted_score_ranked_2x5_zh.{png,pdf,svg}`

Alphabetical comparison grids:

- `kla_ddr_weighted_score_alphabetical_2x5_en.{png,pdf,svg}`
- `kla_ddr_weighted_score_alphabetical_2x5_zh.{png,pdf,svg}`

All graphics are stored in:

`reanalysis/results/figures/kla_ddr_weighted_score_rank_grids_33groups_v1/`

Auditable score, ordering, summary, coefficient, manifest and input-hash tables
are stored in:

`reanalysis/results/tables/kla_ddr_weighted_score_rank_grids_33groups_v1/`

The PNG files are exported at 450 dpi. PDF and SVG files are retained for
vector editing and lossless enlargement.
