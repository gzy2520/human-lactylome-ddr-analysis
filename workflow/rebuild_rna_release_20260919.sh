#!/usr/bin/env bash
# Local downstream rebuild; the expression directory must already have passed server QC.
# Usage: bash workflow/rebuild_rna_release_20260919.sh <expression_dir> <release_label>
set -euo pipefail
cd "$(dirname "$0")/.."
expr_dir="${1:?expression directory required}"
label="${2:?new release label required}"
[[ "$label" =~ ^[0-9]{8}_[A-Za-z0-9_]+$ ]] || exit 2
qs_dir="outputs/${label}_qsmooth"
panel_dir="outputs/${label}_panel"
assisted_dir="outputs/${label}_assisted"
delivery_dir="outputs/${label}_release"
for target_dir in "$qs_dir" "$panel_dir" "$assisted_dir" "$delivery_dir"; do
  if [[ -e "$target_dir" ]]; then
    echo "Refusing existing release path: $target_dir" >&2
    exit 1
  fi
done
test -f "$expr_dir/verification_summary.csv"
mkdir -p "$delivery_dir"
Rscript --vanilla workflow/qsmooth_31group_20260916.R "$expr_dir" "$qs_dir" config/rnaseq_expression_extraction_contract_20260919.csv
python3 workflow/overlay_ddr_panel_20260916.py "$qs_dir" outputs/20260916_ddr_panel_31group "$panel_dir"
cp outputs/20260916_ddr_panel_31group/ddr_uniprot_to_ensembl.tsv "$panel_dir/"
cp outputs/20260916_ddr_panel_31group/mapping_summary.csv "$panel_dir/"
Rscript --vanilla workflow/explore_rna_assisted_ddr_20260917.R . "$assisted_dir" "$qs_dir" "$panel_dir"
Rscript --vanilla workflow/plot_rna_reference_31group_20260917.R . "$delivery_dir/reference" "$qs_dir" "$panel_dir" "$assisted_dir"
Rscript --vanilla workflow/render_teacher_questions_rna_exact.R . "$delivery_dir/core" "$qs_dir" "$expr_dir" "$panel_dir"
Rscript --vanilla workflow/validate_rna_release_20260919.R "$expr_dir" "$qs_dir" "$panel_dir" "$assisted_dir" "$delivery_dir"
echo "REBUILD_VALIDATED: $delivery_dir; publish only after visual QA."
