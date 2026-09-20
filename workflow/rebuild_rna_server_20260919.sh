#!/usr/bin/env bash
# Run ON THE SERVER, using an immutable copy of these four R scripts.
# Usage: bash rebuild_rna_server_20260919.sh <source_root> <fresh_output> <script_dir>
set -euo pipefail
source_root="${1:?source root required}"
output_dir="${2:?fresh output directory required}"
script_dir="${3:?versioned script directory required}"
contract="$source_root/config/rnaseq_expression_extraction_contract_20260919.csv"
Rscript --vanilla "$script_dir/build_31_expression_matrices_20260916.R" "$source_root" "$contract" "$output_dir"
Rscript --vanilla "$script_dir/finalize_31_expression_matrices_20260916.R" "$output_dir" "$contract"
Rscript --vanilla "$script_dir/verify_31_expression_matrices_20260916.R" "$output_dir"
Rscript --vanilla "$script_dir/fingerprint_rna_sources_20260919.R" "$source_root" "$contract" "$output_dir"
