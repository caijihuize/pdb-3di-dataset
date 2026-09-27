#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
MMSEQS="${MMSEQS_BIN:-mmseqs}"
PYTHON="${PDB3DI_PYTHON:-python}"
THREADS="${SLURM_CPUS_PER_TASK:-${THREADS:-16}}"
IN="${ROOT}/data/processed/pdb_v1"
OUT="${AUDIT_DIR:-${ROOT}/data/interim/sequence}"
mkdir -p "${OUT}"
"${PYTHON}" "${ROOT}/src/pdb_3di/split_fasta.py" --input "${IN}/train.fasta" \
  --out-dir "${OUT}/target_chunks" --records-per-shard 50000
for QUERY in valid test; do
  "${MMSEQS}" createdb "${IN}/${QUERY}.fasta" "${OUT}/${QUERY}.db"
  : > "${OUT}/${QUERY}_train.tsv"
  SHARD=0
  for TARGET in "${OUT}"/target_chunks/target_*.fasta; do
    TARGET_DB="${OUT}/${QUERY}_target_${SHARD}.db"
    RESULT_DB="${OUT}/${QUERY}_target_${SHARD}.res"
    "${MMSEQS}" createdb "${TARGET}" "${TARGET_DB}"
    "${MMSEQS}" search "${OUT}/${QUERY}.db" "${TARGET_DB}" "${RESULT_DB}" \
      "${OUT}/tmp_${QUERY}_${SHARD}" --min-seq-id 0.30 -c 0.80 --cov-mode 0 \
      --alignment-mode 3 -s 7.5 --max-seqs 100000 --threads "${THREADS}"
    "${MMSEQS}" convertalis "${OUT}/${QUERY}.db" "${TARGET_DB}" "${RESULT_DB}" \
      "${OUT}/${QUERY}_target_${SHARD}.tsv" \
      --format-output query,target,fident,alnlen,qcov,tcov,evalue,bits
    cat "${OUT}/${QUERY}_target_${SHARD}.tsv" >> "${OUT}/${QUERY}_train.tsv"
    SHARD=$((SHARD + 1))
  done
done
"${MMSEQS}" createdb "${IN}/valid.fasta" "${OUT}/valid_audit.db"
"${MMSEQS}" createdb "${IN}/test.fasta" "${OUT}/test_audit.db"
"${MMSEQS}" search "${OUT}/valid_audit.db" "${OUT}/test_audit.db" \
  "${OUT}/valid_test.res" "${OUT}/tmp_valid_test" --min-seq-id 0.30 \
  -c 0.80 --cov-mode 0 --alignment-mode 3 -s 7.5 --max-seqs 500000 --threads "${THREADS}"
"${MMSEQS}" convertalis "${OUT}/valid_audit.db" "${OUT}/test_audit.db" \
  "${OUT}/valid_test.res" "${OUT}/valid_test.tsv" \
  --format-output query,target,fident,alnlen,qcov,tcov,evalue,bits
