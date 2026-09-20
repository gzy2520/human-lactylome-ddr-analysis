# RNA release 2026-09-19

Current source conversion corrected; 28 distinct RNA references / 1,898 unique sample IDs / 31 proteome mappings.
Statistical models are unchanged at the user's request. Updated inputs are evaluated by the existing methods.

## Deliverables
* `core/`: 14 PDF/PNG pairs plus underlying tables; frozen 196-gene G2M panel, fixed 371-gene DDR panel, 5,224 Kla genes, 48 regulators (49 role entries).
* `reference/`: 4 PDF/PNG pairs plus three-source descriptive tables; 31 rows retain explicit shared-reference mapping.
* `source_changes.csv`: only BPH and TALL-104 source expression values change; all other 26 references remain numerically unchanged.
* `validation.csv`, `visual_qa.csv`: numerical/scope and rendered-output checks.
* `release_sha256.csv`: manifest of this release and its supporting matrices, metadata and scripts (paths relative to repository root).
* `normalization_summary.csv`: separate mean and maximum A/B differences; A and B weights reported separately.

Historical releases are retained for audit and are not merged into this directory.
Source/batch and donor/model limitations remain as documented in the Methods. Passing these checks is not statistical-method endorsement.
