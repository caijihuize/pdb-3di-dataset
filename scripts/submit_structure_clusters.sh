#!/usr/bin/env bash
#SBATCH --job-name=pdb3di_scluster
#SBATCH --partition=qgpu_3090
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --gres=gpu:1
#SBATCH --time=2-00:00:00
#SBATCH --output=runs/structure_cluster_%j.out
#SBATCH --error=runs/structure_cluster_%j.err
set -euo pipefail
ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
cd "${ROOT}"
PYTHON="${PDB3DI_PYTHON:-/hpcfs/fhome/caihuize/miniconda3/envs/ESM3_3Di_3090/bin/python}"
FOLDSEEK="${FOLDSEEK_BIN:-/hpcfs/fhome/caihuize/.local/bin/foldseek}"
export PYTHONPATH="${ROOT}/src${PYTHONPATH:+:${PYTHONPATH}}"
OUT=data/interim/structure_clusters
FSDB=data/interim/foldseek/experimental
mkdir -p "${OUT}" manifests/pdb_v2 runs
tail -n +2 data/interim/candidates/records.tsv | cut -f1 > "${OUT}/candidate_ids.txt"
"${PYTHON}" -m pdb_3di.foldseek_keys \
  --header-db "${FSDB}_h" --header-index "${FSDB}_h.index" \
  --ids "${OUT}/candidate_ids.txt" --output "${OUT}/candidate.keys" \
  --mapping "${OUT}/candidate_structure_ids.tsv" --split candidates
rm -f "${OUT}/candidates"* "${OUT}/clusters"* "${OUT}/clusters.tsv"
rm -rf "${OUT}/tmp"
"${FOLDSEEK}" createsubdb "${OUT}/candidate.keys" "${FSDB}" "${OUT}/candidates"
"${FOLDSEEK}" cluster "${OUT}/candidates" "${OUT}/clusters" "${OUT}/tmp" \
  -e 1e-3 -c 0.80 --cov-mode 0 --alignment-type 2 --cluster-mode 0 \
  --threads "${SLURM_CPUS_PER_TASK}"
"${FOLDSEEK}" createtsv "${OUT}/candidates" "${OUT}/candidates" \
  "${OUT}/clusters" "${OUT}/clusters.tsv"
"${FOLDSEEK}" version > "${OUT}/version.txt"
