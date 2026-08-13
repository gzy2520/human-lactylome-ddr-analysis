# Publication result contract

The repository is considered reproducible only when
`tests/validate_publication_contract.R` passes all of the following:

- 40 catalogued source groups and 37 quantifiable Kla groups.
- 30 paired groups and 28 unique whole-proteome display rows.
- Paired four-class counts of 9 normal tissues, 2 cancer tissues, 7 normal
  cells, and 12 cancer cells.
- 399 unique Kla-DDR BaseAccessions.
- Five protein sets of 183, 292, 178, 381, and 399 proteins.
- 35 five-set by seven-pathway summary rows.
- A 4 x 4 set-count summary for the four Venn analyses, without intersection counts.
- Required bilingual heatmaps, DDR bar/Venn figures, and five-set linear
  pathway matrices. Embeddings and Cytoscape are outside this revision.

The contract checks stable identifiers and analytical counts. It does not
compare raster image bytes because font rendering may vary across operating
systems. Figure existence is supplemented by manual visual inspection of both
language variants.
