# Scope of this input directory: **31 groups**

These files are the **31-group** version of the input set. The 31st group is
`PXD064038 / MEC and NEC ESCC groups`, added on 2026-09-03; the previous 30-group versions
were superseded on 2026-09-16.

## Why the `_30` file names are kept

`group_summary_30.csv`, `kla_protein_membership_30.csv` and the other `_30`-suffixed names
predate the ESCC inclusion and now describe 31 groups. The names are deliberately **not**
changed: they are hard-coded string literals in 16 scripts under `R/` and `workflow/`, and
references outside this repository cannot be checked. The suffix therefore no longer states
the group count — this file and `INPUT_MANIFEST.csv` are authoritative for the scope.

## What changed (30 groups -> 31 groups)

| File | 30 groups | 31 groups |
|---|---|---|
| `group_summary_30.csv` | 30 rows | 31 rows |
| `venn_kla_ddr.csv` (Kla ∩ DDR) | 399 | **401** |
| `venn_reference_ddr.csv` (reference DDR) | 758 | **836** |
| `venn_all_kla.csv` (Kla union) | 5742 | 5814 |
| `venn_reference.csv` (reference union) | 24266 | 24397 |
| `kla_protein_membership_30.csv`, `reference_protein_membership_30.csv` | 21 PXD | 22 PXD |
| `Supplementary_Table_S4_Pathway_Protein_Ranking.xlsx` | 183-protein tumour panel | 192-protein tumour panel |

Every promoted file is a strict **superset** of the version it replaced: nothing was removed
and no accession changed side. The two proteins added to the Kla ∩ DDR union are `P28066`
(PSMA5) and `Q8N122` (RPTOR).

## Provenance

* Source of the promoted files:
  `data/candidate/escc_inclusion_20260903_pxd065830_tumor_reference/publication_input/`
  (produced by `R/candidate/prepare_escc_inclusion_inputs.R`).
* The superseded 30-group files are preserved at
  `work/20260916_publication_input_30group_backup/`.
* `INPUT_MANIFEST.csv` carries the MD5 of every file here; all 13 entries were verified to
  match after the promotion.
* Three files were already identical between the two versions and were left untouched:
  `human_ddr_go_annotations.tsv`, `pathway_display.csv`,
  `Supplementary_Table_S5_Lactylation_Regulators.xlsx`.

## Consequence

Derived products built from the previous 30-group inputs — the figure and table outputs under
`corrected_final_result_20260905_sample_only/` and the audit at
`audits/20260905_final_result/` — correspond to the 30-group scope. They are not regenerated
by this promotion and must be rebuilt before they are used again.
