#!/usr/bin/env bash
# Stage 0 of the 31-group expression extraction: shared stable-ID annotation.
#
# Produces, once, on the server:
#   metadata/annotation/human_gene_lengths_ensembl111.tsv   gene_id(ENSG, version-stripped) <TAB> merged-exon length (bp)
#   metadata/annotation/human_symbol_to_entrez.tsv          GeneSymbol <TAB> GeneID  (NCBI gene_info, taxid 9606)
#   metadata/annotation/hgnc_complete_set.tsv.gz            HGNC complete set (symbol / prev_symbol / entrez_id)
#
# The length table is the single length reference used for every counts->TPM conversion in
# stage 2, so that all count-based groups are normalised identically. No gene symbol is used
# as an analysis identifier anywhere; the symbol table exists only to convert two author
# matrices (GSE171750, GSE283812) whose rows are symbols into official NCBI GeneIDs, and the
# HGNC set supplies the rename history those 2020-era matrices still depend on.
#
# Usage: server_prepare_gene_annotation_20260916.sh <server_root>
set -euo pipefail

root=${1:?usage: $0 <server_root>}
ann=${root}/metadata/annotation
gtf=${ann}/Homo_sapiens.GRCh38.111.gtf.gz
mkdir -p "${ann}"

# ---- 1. Ensembl release 111 GTF (only needed to build the length table) -------------------
if [ ! -s "${gtf}" ]; then
  curl -L --fail --retry 5 --retry-delay 5 \
    -o "${gtf}.part" \
    https://ftp.ensembl.org/pub/release-111/gtf/homo_sapiens/Homo_sapiens.GRCh38.111.gtf.gz
  mv "${gtf}.part" "${gtf}"
fi
echo "GTF_BYTES=$(stat -c %s "${gtf}")"

# ---- 2. Merged-exon gene length table ----------------------------------------------------
# Exons are merged per gene before summing, so a gene's length is the size of its exon union
# rather than the sum of overlapping exon records.
exons=${ann}/human_ensembl111_exons.tsv
zcat "${gtf}" \
  | awk -F'\t' '$3=="exon" {
      if (match($9, /gene_id "[^"]+"/)) {
        gid = substr($9, RSTART+9, RLENGTH-10)
        sub(/\.(UTR|exon|CDS).*/, "", gid)
        print gid"\t"$4"\t"$5
      }
    }' > "${exons}"

sort -k1,1 -k2,2n "${exons}" \
  | awk -F'\t' '
      $1 != g {
        if (g != "") print g"\t"len
        g = $1; cs = $2; ce = $3; len = $3 - $2 + 1; next
      }
      {
        if ($2 <= ce) { if ($3 > ce) { len += $3 - ce; ce = $3 } }
        else          { len += $3 - $2 + 1; ce = $3 }
      }
      END { if (g != "") print g"\t"len }
    ' > "${ann}/human_gene_lengths_ensembl111.tsv"
rm -f "${exons}"
echo "GENES_WITH_LENGTH=$(wc -l < "${ann}/human_gene_lengths_ensembl111.tsv")"

# ---- 3. Official NCBI symbol -> GeneID table (taxid 9606) --------------------------------
gene_info=${ann}/Homo_sapiens.gene_info.gz
if [ ! -s "${gene_info}" ]; then
  curl -L --fail --retry 5 --retry-delay 5 \
    -o "${gene_info}.part" \
    https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz
  mv "${gene_info}.part" "${gene_info}"
fi
zcat "${gene_info}" | awk -F'\t' 'BEGIN{OFS="\t"} $1=="9606" {print $3, $2}' > "${ann}/human_symbol_to_entrez.tsv"
echo "SYMBOL_MAP_ROWS=$(wc -l < "${ann}/human_symbol_to_entrez.tsv")"

# ---- 4. HGNC complete set (rename history for the symbol-routed author matrices) ----------
# gene_info above carries only today's official symbol, so a matrix written when AARS was
# still called AARS loses every gene that has been renamed since. HGNC is the nomenclature
# authority and records those renames in prev_symbol, which build_symbol_lookup() layers on
# top of the official table.
hgnc=${ann}/hgnc_complete_set.tsv.gz
if [ ! -s "${hgnc}" ]; then
  curl -L --fail --retry 5 --retry-delay 5 \
    -o "${hgnc}.part" \
    https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt
  gzip -c "${hgnc}.part" > "${hgnc}"
  rm -f "${hgnc}.part"
fi
echo "HGNC_ROWS=$(zcat "${hgnc}" | wc -l)"
echo "ANNOTATION_PREP_DONE"
