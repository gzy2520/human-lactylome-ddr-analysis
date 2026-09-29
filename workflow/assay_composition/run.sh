#!/bin/bash
set -euo pipefail
# Run from the repository root. Requires fresh labels/model/external output directories.
mkdir -p outputs/20260929_assay_composition/logs
Rscript --vanilla workflow/assay_composition/labels.R
Rscript --vanilla workflow/assay_composition/external_labels.R outputs/20260929_assay_composition/external_labels
for target in Kla_DDR Kla_DNA_repair Proteome_DDR Proteome_DNA_repair; do
 Rscript --vanilla workflow/assay_composition/train.R outputs/20260926_full_rna_inputs "outputs/20260929_assay_composition/models/$target" "$target" > "outputs/20260929_assay_composition/logs/$target.log" 2>&1
done
Rscript --vanilla workflow/assay_composition/predict_external.R outputs/20260928_external_escc_rna_corrected/full_rna.rds outputs/20260929_assay_composition/external/newRNA newRNA
Rscript --vanilla workflow/assay_composition/predict_external.R outputs/20260926_full_rna_inputs/full_rna.rds outputs/20260929_assay_composition/external/reused_TCGA reused_TCGA
Rscript --vanilla workflow/assay_composition/summarize.R
# Generate the report using Python with python-docx installed.
"${PYTHON_BIN:-python3}" workflow/assay_composition/report.py
