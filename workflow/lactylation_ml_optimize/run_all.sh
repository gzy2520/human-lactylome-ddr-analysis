#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/../.."
if [ -e outputs/20260922_br_v5/oos_predictions.csv.gz ]; then
  echo 'Existing v5 predictions: refusing overwrite. Use an isolated copy for reproduction.' >&2
  exit 1
fi
mkdir -p audit/20260922_br_v5
for stage in tests run evaluate sensitivity tests verify_results report; do
  Rscript --vanilla "workflow/lactylation_ml_optimize/${stage}.R" > "audit/20260922_br_v5/${stage}.log" 2>&1
done
