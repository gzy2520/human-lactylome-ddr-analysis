# NSC identity-resolution decision (2026-09-16)

## Established evidence

- KLA31_30 is registered as **neural stem cells**. Its proteome comparison uses the two cell models `ENSA` and `HMP1` (PXD070007; reference PXD069969).
- GSE119834 has nine bulk-RNA-seq neural-stem-cell libraries. The source annotation states that these RNA-seq cells were not treated.
- One RNA library is the same cell model as the proteome: `GSM3384849`, model `ENSA`.
- The other proteome model, `HMP1`, is absent from GSE119834. The remaining eight libraries are distinct NSC models; their names and the exact match flag are in `nsc_rna_vs_proteome_model_map.csv`.

## What this means for the current inclusion rules

GSE119834 passes the stated material category (neural stem cells) and no-perturbation requirement. It does **not** establish a two-model ENSA/HMP1 RNA match. Its nine libraries should be described as independent external NSC models, not technical replicates of ENSA or HMP1.

## Inclusion decision (recorded 2026-09-16)

The user approved category-matched NSC libraries when they have no extra treatment. Therefore all nine GSE119834 NSC libraries are selected as the RNA reference for KLA31_30.
The next gate is source-matrix QC: validate the downloaded author FPKM matrix, retain its RefSeq stable identifiers, verify all nine selected columns and finite numeric values, and then run source-local analysis. The ENSA/HMP1 distinction remains a provenance caveat only.

No gene symbols are used as analytical join keys in either path. The model labels above are sample-identity annotations only.
