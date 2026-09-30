#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT="${ROOT}/data/raw/mappings"
BASE="https://ftp.ebi.ac.uk/pub/databases/msd/sifts/csv"
mkdir -p "${OUT}"
for NAME in pdb_chain_uniprot.csv pdb_chain_pfam.csv pdb_chain_cath_uniprot.csv; do
  curl -L --fail --retry 3 -o "${OUT}/${NAME}.part" "${BASE}/${NAME}"
  mv "${OUT}/${NAME}.part" "${OUT}/${NAME}"
done
curl -L --fail --retry 3 -o "${OUT}/cath-domain-list.txt.part" \
  "https://download.cathdb.info/cath/releases/latest-release/cath-classification-data/cath-domain-list.txt"
mv "${OUT}/cath-domain-list.txt.part" "${OUT}/cath-domain-list.txt"
(cd "${OUT}" && sha256sum pdb_chain_uniprot.csv pdb_chain_pfam.csv \
  pdb_chain_cath_uniprot.csv cath-domain-list.txt > SHA256SUMS)
{
  printf 'file\turl\tsifts_release_header\n'
  for NAME in pdb_chain_uniprot.csv pdb_chain_pfam.csv pdb_chain_cath_uniprot.csv; do
    HEADER=$(head -1 "${OUT}/${NAME}" | tr -d '\r\n\t')
    printf '%s\t%s/%s\t%s\n' "${NAME}" "${BASE}" "${NAME}" "${HEADER}"
  done
  HEADER=$(head -1 "${OUT}/cath-domain-list.txt" | tr -d '\r\n\t')
  printf '%s\t%s\t%s\n' "cath-domain-list.txt" \
    "https://download.cathdb.info/cath/releases/latest-release/cath-classification-data/cath-domain-list.txt" \
    "${HEADER}"
} > "${OUT}/sources.tsv"
