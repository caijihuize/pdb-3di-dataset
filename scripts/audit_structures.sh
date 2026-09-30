#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
FOLDSEEK="${FOLDSEEK_BIN:-foldseek}"
PYTHON="${PDB3DI_PYTHON:-python}"
THREADS="${SLURM_CPUS_PER_TASK:-${THREADS:-16}}"
FSDB="${ROOT}/data/interim/foldseek/experimental"
OUT="${STRUCT_AUDIT_DIR:-${ROOT}/data/interim/structure}"
MANIFEST="${ROOT}/manifests/pdb_v1"
mkdir -p "${OUT}"
find "${OUT}" -mindepth 1 -depth -delete
for SPLIT in train valid test; do
  "${PYTHON}" "${ROOT}/src/pdb_3di/foldseek_keys.py" \
    --header-db "${FSDB}_h" --header-index "${FSDB}_h.index" \
    --ids "${MANIFEST}/${SPLIT}_ids.txt" --output "${OUT}/${SPLIT}.keys" \
    --mapping "${MANIFEST}/${SPLIT}_structure_ids.tsv" --split "${SPLIT}"
  "${FOLDSEEK}" createsubdb "${OUT}/${SPLIT}.keys" "${FSDB}" "${OUT}/${SPLIT}"
done
for QUERY in valid test; do
  "${FOLDSEEK}" search "${OUT}/${QUERY}" "${OUT}/train" \
    "${OUT}/${QUERY}_train.res" "${OUT}/tmp_${QUERY}" \
    --alignment-type 2 -a 1 -e 1e-3 -c 0.50 --cov-mode 0 \
    --max-seqs 10 --threads "${THREADS}"
  "${FOLDSEEK}" convertalis "${OUT}/${QUERY}" "${OUT}/train" \
    "${OUT}/${QUERY}_train.res" "${OUT}/${QUERY}_train.tsv" \
    --format-output query,target,fident,alnlen,qcov,tcov,evalue,bits,qtmscore,ttmscore
done
"${FOLDSEEK}" search "${OUT}/valid" "${OUT}/test" \
  "${OUT}/valid_test.res" "${OUT}/tmp_valid_test" \
  --alignment-type 2 -a 1 -e 1e-3 -c 0.50 --cov-mode 0 \
  --max-seqs 10 --threads "${THREADS}"
"${FOLDSEEK}" convertalis "${OUT}/valid" "${OUT}/test" \
  "${OUT}/valid_test.res" "${OUT}/valid_test.tsv" \
  --format-output query,target,fident,alnlen,qcov,tcov,evalue,bits,qtmscore,ttmscore
