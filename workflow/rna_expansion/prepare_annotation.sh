#!/usr/bin/env bash
# Build the same Ensembl 111 merged-exon length definition as the 31-group RNA
# pipeline, and freeze the human subset of NCBI gene2ensembl for Entrez mapping.
set -euo pipefail

gtf=${1:?usage: prepare_annotation.sh <Ensembl111.gtf.gz> <gene2ensembl.gz> <output_dir>}
gene2ensembl=${2:?usage: prepare_annotation.sh <Ensembl111.gtf.gz> <gene2ensembl.gz or -> <output_dir>}
out=${3:?usage: prepare_annotation.sh <Ensembl111.gtf.gz> <gene2ensembl.gz> <output_dir>}
mkdir -p "$out"
gzip -t "$gtf"
if [ "$gene2ensembl" != "-" ]; then gzip -t "$gene2ensembl"; fi

if [ ! -s "$out/human_gene_lengths_ensembl111.tsv" ]; then
  exons="$out/.human_ensembl111_exons.tmp.tsv"
  gzip -dc "$gtf" \
    | awk -F '\t' '$3 == "exon" && match($9, /gene_id "[^"]+"/) {
        gid = substr($9, RSTART + 9, RLENGTH - 10)
        print gid "\t" $4 "\t" $5
      }' > "$exons"
  LC_ALL=C sort -T "$out" -k1,1 -k2,2n "$exons" \
    | awk -F '\t' '
        $1 != g {
          if (g != "") print g "\t" len
          g = $1; cs = $2; ce = $3; len = $3 - $2 + 1; next
        }
        {
          if ($2 <= ce) { if ($3 > ce) { len += $3 - ce; ce = $3 } }
          else { len += $3 - $2 + 1; ce = $3 }
        }
        END { if (g != "") print g "\t" len }
      ' > "$out/human_gene_lengths_ensembl111.tsv"
  rm "$exons"
fi

if [ "$gene2ensembl" != "-" ]; then
  gzip -dc "$gene2ensembl" \
    | awk -F '\t' 'BEGIN { OFS = "\t" } NR == 1 || $1 == "9606" { print }' \
    > "$out/human_gene2ensembl.tsv"
fi

echo "length_rows=$(wc -l < "$out/human_gene_lengths_ensembl111.tsv")"
if [ "$gene2ensembl" != "-" ]; then
  echo "mapping_rows=$(wc -l < "$out/human_gene2ensembl.tsv")"
fi
