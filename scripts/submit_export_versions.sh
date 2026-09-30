#!/usr/bin/env bash
#SBATCH --job-name=pdb3di_export_v2
#SBATCH --partition=qgpu_3090
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --gres=gpu:1
#SBATCH --time=1-00:00:00
#SBATCH --output=runs/export_versions_%j.out
#SBATCH --error=runs/export_versions_%j.err
set -euo pipefail
ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
cd "${ROOT}"
PYTHON="${PDB3DI_PYTHON:-/hpcfs/fhome/caihuize/miniconda3/envs/ESM3_3Di_3090/bin/python}"
export PYTHONPATH="${ROOT}/src${PYTHONPATH:+:${PYTHONPATH}}"
"${PYTHON}" -m pdb_3di.export_versions \
  --records data/interim/candidates/records.tsv \
  --annotations data/interim/annotations/record_annotations.tsv \
  --structure-key-map data/interim/structure_clusters/candidate_structure_ids.tsv \
  --manifest-root manifests --output-root data/processed
"${PYTHON}" -m pdb_3di.validate_versions --root "${ROOT}" \
  --output reports/pdb_v2/validation.json
