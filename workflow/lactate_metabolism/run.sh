#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/../.."
output_dir=${1:?Provide a NEW empty output directory}
Rscript --vanilla workflow/preflight_kla31.R
Rscript --vanilla workflow/lactate_metabolism/analyze.R "$output_dir"
Rscript --vanilla workflow/lactate_metabolism/validate_catalog.R "$output_dir"
Rscript --vanilla workflow/lactate_metabolism/report.R "$output_dir"
# Visual inspection remains necessary; this command renders, but does not approve, PDFs.
if [[ -n "${KLA_QA_PYTHON:-}" ]]; then
  "$KLA_QA_PYTHON" workflow/lactate_metabolism/render_qa.py "$output_dir" "$output_dir/visual_preview"
fi
Rscript --vanilla workflow/preflight_kla31.R
