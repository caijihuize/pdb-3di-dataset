#!/usr/bin/env bash
#SBATCH --job-name=pdb3di_annotate
#SBATCH --partition=qgpu_3090
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --gres=gpu:1
#SBATCH --time=2-00:00:00
#SBATCH --output=runs/annotate_%j.out
#SBATCH --error=runs/annotate_%j.err
set -euo pipefail
ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
cd "${ROOT}"
PYTHON="${PDB3DI_PYTHON:-/hpcfs/fhome/caihuize/miniconda3/envs/ESM3_3Di_3090/bin/python}"
export PYTHONPATH="${ROOT}/src${PYTHONPATH:+:${PYTHONPATH}}"
mkdir -p data/interim/annotations reports/pdb_v1_1 runs
"${PYTHON}" -m pdb_3di.annotations \
  --records data/interim/candidates/records.tsv \
  --raw-dir data/raw/structures --sifts-dir data/raw/mappings \
  --output data/interim/annotations/record_annotations.tsv \
  --errors reports/pdb_v1_1/annotation_errors.tsv \
  --workers "${SLURM_CPUS_PER_TASK}"
