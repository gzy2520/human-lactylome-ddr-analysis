#!/usr/bin/env bash
# Run from repository root. Download only into a NEW directory, then checksum/inspect.
set -euo pipefail
out=${1:?Usage: fetch.sh NEW_OUTPUT_DIRECTORY}
if [[ -e "$out" ]]; then echo 'Fresh output directory required' >&2; exit 1; fi
mkdir -p "$out"
curl -fLsS --retry 2 --max-time 300 'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13195773/fullTextXML' -o "$out/PMC13195773.xml"
curl -fLsS --retry 2 --max-time 300 'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13195773/supplementaryFiles' -o "$out/PMC13195773_supp_complete.zip"
curl -fLsS --retry 2 --max-time 300 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE130078&format=file' -o "$out/GSE130078_download.tar"
curl -fLsS --retry 2 --max-time 300 'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz' -o "$out/gencode.v19.annotation.gtf.gz"
curl -fLsS --retry 2 --max-time 300 'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC6900598/fullTextXML' -o "$out/PMC6900598.xml"
curl -fLsS --retry 2 --max-time 300 'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC6900598/supplementaryFiles' -o "$out/PMC6900598_supp.zip"
unzip -t "$out/PMC13195773_supp_complete.zip"
unzip -t "$out/PMC6900598_supp.zip"
gzip -t "$out/gencode.v19.annotation.gtf.gz"
tar -tf "$out/GSE130078_download.tar"
