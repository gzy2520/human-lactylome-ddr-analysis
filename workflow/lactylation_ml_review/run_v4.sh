#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/../.."
Rscript --vanilla workflow/preflight_kla31.R
Rscript --vanilla workflow/lactylation_ml_review/06_test_v4_weighting.R
Rscript --vanilla workflow/lactylation_ml_review/07_prepare_v4.R
Rscript --vanilla workflow/lactylation_ml_review/08_train_BR_v4.R
BR_V4_OUT=outputs/20260921_br_v4_weight_only Rscript --vanilla workflow/lactylation_ml_review/08_train_BR_v4.R
Rscript --vanilla workflow/lactylation_ml_review/09_sensitivity_v4.R
Rscript --vanilla workflow/lactylation_ml_review/10_acceptance_tests_v4.R
Rscript --vanilla workflow/lactylation_ml_review/10_additional_tests_v4.R
Rscript --vanilla workflow/lactylation_ml_review/11_figures_v4.R
Rscript --vanilla workflow/lactylation_ml_review/12_report_v4.R
