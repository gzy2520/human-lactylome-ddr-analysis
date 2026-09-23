#!/usr/bin/env bash
# Fetch original GEO/annotation files into the ignored candidate-data area.
# Files are verified before use; the >100 MB author archive stays out of Git.
set -euo pipefail
root=${1:-data/candidate/rna_expansion}

fetch_one() {
  local relative=$1 expected=$2 url=$3 target="$root/$1" actual tmp
  mkdir -p "$(dirname "$target")"
  if [ ! -s "$target" ]; then
    tmp="$target.part.$$"
    trap 'rm -f "$tmp"' EXIT
    curl --fail --location --retry 3 --continue-at - --output "$tmp" "$url"
    mv "$tmp" "$target"
    trap - EXIT
  fi
  actual=$(shasum -a 256 "$target" | awk '{print $1}')
  if [ "$actual" != "$expected" ]; then
    printf 'SHA-256 mismatch: %s\nexpected=%s\nactual=%s\n' \
      "$target" "$expected" "$actual" >&2
    exit 1
  fi
  printf 'OK %s\n' "$relative"
}

fetch_one GSE178411/GSE178411_counts.txt.gz \
  19622a1b543d9b67481ca5bb13e35f73aafd28fe56f62dde3309e91c37ad0228 \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE178nnn/GSE178411/suppl/GSE178411_counts.txt.gz
fetch_one GSE178411/GSE178411_series_matrix.txt.gz \
  1d94d0e62bffb500da013e3be8d5c212f925162cf0e1b66f9e94e1721672e9f4 \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE178nnn/GSE178411/matrix/GSE178411_series_matrix.txt.gz
fetch_one GSE180836/GSE180836_RAW.tar \
  a7cb6e33131c7edfaab4877c3592d6dae0c021148da889b947617d9388ae7368 \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE180nnn/GSE180836/suppl/GSE180836_RAW.tar
fetch_one GSE180836/GSE180836_series_matrix.txt.gz \
  f026fd6f48664421941aebafb064f1b41d7604ed9a4cf4bf1661802022f8e6b6 \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE180nnn/GSE180836/matrix/GSE180836_series_matrix.txt.gz
fetch_one annotation/Homo_sapiens.GRCh38.111.gtf.gz \
  a52356765a41264e17a2076aff1abab703ec3ba239b78f88560756e85c169831 \
  https://ftp.ensembl.org/pub/release-111/gtf/homo_sapiens/Homo_sapiens.GRCh38.111.gtf.gz
fetch_one annotation/gene2ensembl.gz \
  51168dd97bcc57535aaab88fae221dbeaed1c62dab7aa4239bacad346f8a0043 \
  https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2ensembl.gz
