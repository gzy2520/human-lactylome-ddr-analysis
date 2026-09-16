#!/usr/bin/env bash
# Transcript -> gene map for the one group whose source matrix is transcript-level and was
# quantified against an older (GRCh37-era) transcriptome: GSE114691 placenta.
#
# NCBI's gene2ensembl table only links transcripts that have a RefSeq RNA counterpart, and it
# is keyed to current annotations, so it resolves only ~23% of that series' older ENST IDs.
# Adding the Ensembl GRCh37 transcript definitions raises the resolvable share substantially.
# Ensembl gene identifiers are stable across assemblies, so both maps land on the same ENSG
# space used everywhere else.
#
# Usage: server_prepare_grch37_transcript_map_20260916.sh <server_root>
set -euo pipefail

root=${1:?usage: $0 <server_root>}
ann=${root}/metadata/annotation
gtf=${ann}/Homo_sapiens.GRCh37.87.gtf.gz
mkdir -p "${ann}"

if [ ! -s "${gtf}" ]; then
  curl -L --fail --retry 5 --retry-delay 5 \
    -o "${gtf}.part" \
    https://ftp.ensembl.org/pub/grch37/release-87/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz
  mv "${gtf}.part" "${gtf}"
fi
echo "GTF37_BYTES=$(stat -c %s "${gtf}")"

zcat "${gtf}" | awk -F'\t' '$3=="transcript" {
    if (match($9, /gene_id "[^"]+"/)) {
      g = substr($9, RSTART+9, RLENGTH-10)
      if (match($9, /transcript_id "[^"]+"/)) {
        t = substr($9, RSTART+15, RLENGTH-16)
        print t"\t"g
      }
    }
  }' | sort -u > "${ann}/human_grch37_transcript_to_gene.tsv"

echo "GRCH37_TRANSCRIPTS=$(wc -l < "${ann}/human_grch37_transcript_to_gene.tsv")"
echo "GRCH37_MAP_DONE"
