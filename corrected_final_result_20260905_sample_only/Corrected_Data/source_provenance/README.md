# ESCC inclusion audit (2026-09-03)

- Selected representative group: `PXD064038 / MEC and NEC ESCC groups`. MEC and NEC are retained as six source observations under one cancer-tissue group.
- Excluded from this update: PXD048995 (KYSE30 histone-focused) and PXD063945 (neoadjuvant ESCC).
- Parsed PXD064038: 1239 Kla BaseAccessions; 92 Kla-DDR BaseAccessions.
- Non-lactylated control: iProX/PXD065830 Dataset1, sheet 2.a, using only 94 ESCC tumor T samples (the 24 N non-tumor samples are excluded). The union contains 8083 BaseAccessions and 420 DDR BaseAccessions (5.196%).
- PXD053809 is not selected for this reproducible run because its deposited processed XLS is encrypted; PXD010154 remains a healthy-esophagus background and is not the control.
- Added 72 new all-Kla union proteins and 2 new Kla-DDR union proteins.
- The expanded S4 tumor panel adds 14 cancer-tissue Kla-DDR proteins: existing seven-pathway states are retained where frozen S4 already curated the accession; otherwise broad GO-only accessions receive zero pathway states.
- Dataset-level boxplot inputs are regenerated from the expanded group unions: one point per PXD/sample-group and seven pathways for each of Kla and whole proteome.

The dated input and candidate output directories are isolated from the frozen publication input and approved result directories.
