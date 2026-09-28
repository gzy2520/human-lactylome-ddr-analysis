# Scientific correction of the supplied DDR / DNA repair report

Entry point: `revise_doc.py` edits a copy of the user's Desktop Word document; the source document and earlier result directories are retained. `render_extended.R` and `render_top25.R` render corrected terminology from existing predictions and importance values, without refitting models. The UMAP explicitly uses full-data fitted predictions and is descriptive only.

Reproduction from the repository root (R with global package library, Python with python-docx):

1. `Rscript --vanilla workflow/external_material/prepare_rna.R outputs/20260928_external_escc_rna_corrected`
2. `Rscript --vanilla workflow/external_material/predict.R outputs/20260928_external_escc_rna_corrected/full_rna.rds outputs/20260928_report_revision/external newRNA`
3. `Rscript --vanilla workflow/report_revision/render_extended.R`
4. `Rscript --vanilla workflow/report_revision/render_top25.R`
5. Run `workflow/report_revision/revise_doc.py` with Python.

Steps 1–2 require fresh output directories; do not overwrite historical outputs. The previous external preparation used log2(TPM+1), while frozen training used log2(TPM+0.5). This correction aligns the external transformation to training. Raw-count reconstruction matched the independently transformed old matrix within 1.78e-15. Corrected MAEs (percentage points): DDR tumor 12.30183, adjacent 15.91094; DNA repair tumor 14.25926, adjacent 18.53333. The older external tables and figure 13 are superseded for reporting.

The original general report builder under `workflow/ddr_fraction/` is retained as historical code; use this revision entry point for the corrected supplied report. Model weights, labels and internal evaluation predictions are unchanged. External evaluation compares 46 RNA records from 23 paired donors with two shared material references per target; there is no patient-level RNA–proteome matching, and RNA pretreatment history is unverified.
