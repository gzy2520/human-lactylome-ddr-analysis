# RNA reference panels (2026-09-17)

Three sources meet on every (material x DDR protein) pair: the lactylome (Kla), the
matched un-enriched whole proteome (Ref), and the transcriptome (RNA). These panels record
the RNA side of that: the qsmooth-smoothed expression profile, and how the transcript level
relates to what each proteomics assay picks up.

RNA_1   Expression of the 371 Kla n DDR genes across the 31 material classes. Faceted by DDR
        pathway, genes ordered by mean expression within each block, columns ordered by
        material category, row-scaled so pattern rather than absolute level is visible.
RNA_2a  Joint state of all 11,501 pairs: the four Kla/Ref combinations, each split by whether
        the gene's transcript is expressed (TPM >= 1) or low (TPM < 1).
RNA_2b  The same pairs read from the RNA side: detection rate of each assay at expressed (TPM >= 1)
        and low (TPM < 1) transcript level.
RNA_2   Composite of 2a and 2b.

Scope: descriptive. The panels record joint states and detection rates; no statistical test
is applied and no mechanism is inferred.
The RNA reference is a material-class profile drawn from different studies than the proteome,
so groups are compared as classes, not as paired samples.
Expression is called at an absolute TPM >= 1, not a within-group rank, so the split differs
between materials; the current assisted directory records the threshold sensitivity.
