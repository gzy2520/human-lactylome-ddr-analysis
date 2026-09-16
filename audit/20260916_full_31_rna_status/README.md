# Full 31-group RNA status (2026-09-16)

This is the master ledger for all 31 proteome/Kla group rows. It distinguishes the 13 rows already in source-local analysis from 18 rows with confirmed candidates queued for selection or matrix QC. Zero rows are currently blocked.
Queued means a verified eligible candidate source exists with matching unperturbed biological replicates (or accepted n=1 descriptive reference), pending stable-ID mapping, sample freezing, or local matrix QC.
The 13 started rows use material/cell identity matching and selected samples with no knockdown, overexpression or experimental drug. DMSO vehicle controls remain explicitly marked and separated from untreated sources.

ExpressionStatus records the outcome of the 2026-09-16 server-side expression build: every one of the 31 group rows now has a per-group log2(TPM + 0.5) matrix keyed by stable Ensembl gene identifiers, with cleaned counts archived alongside where the source provides them. ExpressionSamples and ExpressionGenes give the size of that matrix; three rows share the HCT116 matrix and two share the HK-2 matrix, so 31 rows map onto 28 distinct matrices.
The matrices are expression profiles only. No row was compared against another here; the cross-tissue comparison is a separate downstream step.
Limitations entries carrying the "[2026-09-16 expression build]" marker are caveats established while assembling those matrices.
