# Data boundary

`publication_input/` is the only data directory read by the final workflow.
It is a frozen, publication-scoped snapshot with exactly 30 Kla sample groups.
The files are standardized analysis inputs rather than raw mass-spectrometry
files; each is listed with a checksum in `publication_input/INPUT_MANIFEST.csv`.

`Supplementary_Table_S4_Pathway_Protein_Ranking.xlsx` and
`Supplementary_Table_S5_Lactylation_Regulators.xlsx` are author-approved
release assets. The workflow reads S4 for pathway figures and copies both files
unchanged to the final supplementary-output directory.

The final workflow never scans the local PXD download cache. Public source
accessions and original evidence-file locators remain in
`publication_input/group_summary_30.csv` and Supplementary Table S1.
