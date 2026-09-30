#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
STAGE="${1:?Usage: scripts/run.sh discover|download|verify|features|candidates|cluster|split|export|sequence-audit|repair|structure-audit}"
PYTHON="${PDB3DI_PYTHON:-python}"
FOLDSEEK="${FOLDSEEK_BIN:-foldseek}"
MMSEQS="${MMSEQS_BIN:-mmseqs}"
THREADS="${SLURM_CPUS_PER_TASK:-${THREADS:-16}}"
VERSION="pdb_v1"
RAW="${ROOT}/data/raw/structures"
INTERIM="${ROOT}/data/interim"
CANDIDATES="${INTERIM}/candidates"
MANIFEST="${ROOT}/manifests/${VERSION}"
PROCESSED="${ROOT}/data/processed/${VERSION}"
DB="${INTERIM}/foldseek/experimental"
mkdir -p "${INTERIM}" "${MANIFEST}" "${PROCESSED}" "${ROOT}/runs"
export PYTHONPATH="${ROOT}/src${PYTHONPATH:+:${PYTHONPATH}}"

case "${STAGE}" in
  discover)
    "${PYTHON}" -m pdb_3di.source discover --out-dir "${MANIFEST}"
    ;;
  download)
    DOWNLOAD_ARGS=()
    if [[ -n "${DOWNLOAD_LIMIT:-}" ]]; then DOWNLOAD_ARGS+=(--limit "${DOWNLOAD_LIMIT}"); fi
    DOWNLOAD_ARGS+=(--shard-index "${DOWNLOAD_SHARD_INDEX:-0}" --shard-count "${DOWNLOAD_SHARD_COUNT:-1}")
    DOWNLOAD_LOG="${MANIFEST}/downloads.tsv"
    if [[ "${DOWNLOAD_SHARD_COUNT:-1}" != "1" ]]; then
      DOWNLOAD_LOG="${MANIFEST}/downloads_${DOWNLOAD_SHARD_INDEX}.tsv"
    fi
    "${PYTHON}" -m pdb_3di.source download --ids "${MANIFEST}/source_ids.txt" \
      --raw-dir "${RAW}" --log "${DOWNLOAD_LOG}" \
      --workers "${DOWNLOAD_WORKERS:-8}" "${DOWNLOAD_ARGS[@]}"
    ;;
  verify)
    "${PYTHON}" -m pdb_3di.source verify --ids "${MANIFEST}/source_ids.txt" \
      --raw-dir "${RAW}" --log-dir "${MANIFEST}" \
      --report "${ROOT}/reports/${VERSION}/source_validation.json"
    ;;
  features)
    mkdir -p "$(dirname "${DB}")"
    "${FOLDSEEK}" createdb "${RAW}" "${DB}" --file-include '\.cif\.gz$' \
      --chain-name-mode 1 --threads "${THREADS}"
    "${FOLDSEEK}" convert2fasta "${DB}" "${INTERIM}/aa.fasta"
    "${FOLDSEEK}" lndb "${DB}_h" "${DB}_ss_h"
    "${FOLDSEEK}" convert2fasta "${DB}_ss" "${INTERIM}/3di.fasta"
    "${FOLDSEEK}" version > "${INTERIM}/foldseek/version.txt"
    ;;
  candidates)
    "${PYTHON}" -m pdb_3di.candidates --raw-dir "${RAW}" \
      --aa "${INTERIM}/aa.fasta" --di "${INTERIM}/3di.fasta" \
      --out-dir "${CANDIDATES}" --min-length 30 \
      --max-xray-resolution 3.0 --max-em-resolution 4.0 \
      --workers "${METADATA_WORKERS:-4}"
    ;;
  cluster)
    mkdir -p "${INTERIM}/mmseqs"
    rm -f "${INTERIM}/mmseqs/clusterdb"* "${INTERIM}/mmseqs/clusters.tsv"
    rm -rf "${INTERIM}/mmseqs/tmp"
    "${MMSEQS}" createdb "${CANDIDATES}/aa.filtered.fasta" "${INTERIM}/mmseqs/seqdb"
    "${MMSEQS}" cluster "${INTERIM}/mmseqs/seqdb" "${INTERIM}/mmseqs/clusterdb" \
      "${INTERIM}/mmseqs/tmp" --min-seq-id 0.30 -c 0.80 --cov-mode 0 \
      --alignment-mode 3 --cluster-mode 0 -s 7.5 --max-seqs 300 --threads "${THREADS}"
    "${MMSEQS}" createtsv "${INTERIM}/mmseqs/seqdb" "${INTERIM}/mmseqs/seqdb" \
      "${INTERIM}/mmseqs/clusterdb" "${INTERIM}/mmseqs/clusters.tsv"
    "${MMSEQS}" version > "${INTERIM}/mmseqs/version.txt"
    ;;
  split)
    rm -f "${MANIFEST}/repair_state.json"
    "${PYTHON}" -m pdb_3di.split_clusters --records "${CANDIDATES}/records.tsv" \
      --clusters "${INTERIM}/mmseqs/clusters.tsv" --out-dir "${MANIFEST}" \
      --seed 42 --valid-groups "${VALID_GROUPS:-1000}" \
      --test-groups "${TEST_GROUPS:-1000}" --eval-min-length 30 --eval-max-length 510
    ;;
  export)
    EXPORT_ARGS=()
    if [[ "${SKIP_STRUCTURES:-0}" == "1" ]]; then EXPORT_ARGS+=(--skip-structures); fi
    "${PYTHON}" -m pdb_3di.export --records "${CANDIDATES}/records.tsv" \
      --manifest-dir "${MANIFEST}" --raw-dir "${RAW}" --out-dir "${PROCESSED}" \
      "${EXPORT_ARGS[@]}"
    ;;
  sequence-audit)
    bash "${ROOT}/scripts/audit_sequences.sh"
    ;;
  repair)
    "${PYTHON}" -m pdb_3di.repair_sequence_conflicts --records "${CANDIDATES}/records.tsv" \
      --manifest-dir "${MANIFEST}" --audit-dir "${AUDIT_DIR:-${INTERIM}/sequence}" \
      --valid-count "${VALID_GROUPS:-1000}" --test-count "${TEST_GROUPS:-1000}"
    ;;
  structure-audit)
    bash "${ROOT}/scripts/audit_structures.sh"
    ;;
  *) echo "Unknown stage: ${STAGE}" >&2; exit 2 ;;
esac
