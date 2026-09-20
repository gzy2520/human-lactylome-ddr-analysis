# Frozen Hallmark G2M panel

This CSV is copied without changing values from the pre-repair accepted export:
`results/rna_teacher_questions/official_hallmark_g2m_checkpoint_genes_mapped.csv`.
It contains the 196 Ensembl genes already used in the current 18,332-gene analysis space, plus Entrez and display symbols. It is a mapped subset of HALLMARK_G2M_CHECKPOINT, not a claim that the complete upstream gene set contains only 196 genes.

The original export was produced by `msigdbr(species="Homo sapiens", collection="H")`. Its CSV did not record an upstream database version; this repair does not invent one. The fixed CSV and its hashes now define the exact analysis panel, avoiding silent upstream changes during plotting. The regression comparison confirms unchanged scores when the input matrices are unchanged.
