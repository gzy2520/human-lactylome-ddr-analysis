# Publication result contract

The repository is considered reproducible only when
`tests/validate_publication_contract.R` passes all of the following:

- 40 catalogued source groups and 37 quantifiable Kla groups.
- 33 paired groups and 30 unique whole-proteome display rows.
- Paired four-class counts of 9 normal tissues, 2 cancer tissues, 9 normal
  cells, and 13 cancer cells.
- 507 unique Kla-DDR BaseAccessions.
- Five protein sets of 183, 471, 178, 383, and 507 proteins.
- 507 by 3,008 BP semantic feature matrix.
- 1,175 signed nonzero pathway assignments: 1,108 promoting and 67 suppressing.
- 35 five-set by seven-pathway summary rows.
- Required bilingual heatmaps, DDR figures, pathway-specific UMAP, and
  five-set pathway matrices.

The contract checks stable identifiers and analytical counts. It does not
compare raster image bytes because font rendering may vary across operating
systems. Figure existence is supplemented by manual visual inspection of both
language variants.
