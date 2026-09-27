#!/usr/bin/env bash
#SBATCH --job-name=pdb3di_download
#SBATCH --partition=qfree
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=03:00:00
#SBATCH --array=0-15%2
#SBATCH --output=runs/download_%A_%a.out
#SBATCH --error=runs/download_%A_%a.err
set -euo pipefail
ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
cd "${ROOT}"
export PDB3DI_PYTHON="${PDB3DI_PYTHON:-/hpcfs/fhome/caihuize/miniconda3/envs/ESM3_3Di_3090/bin/python}"
export DOWNLOAD_SHARD_INDEX="${SLURM_ARRAY_TASK_ID}"
export DOWNLOAD_SHARD_COUNT=16
export DOWNLOAD_WORKERS=8
bash scripts/run.sh download
