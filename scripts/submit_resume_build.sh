#!/usr/bin/env bash
#SBATCH --job-name=pdb3di_cluster
#SBATCH --partition=qgpu_3090
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --gres=gpu:1
#SBATCH --time=2-00:00:00
#SBATCH --output=runs/cluster_%j.out
#SBATCH --error=runs/cluster_%j.err
set -euo pipefail
ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
cd "${ROOT}"
export PDB3DI_PYTHON="${PDB3DI_PYTHON:-/hpcfs/fhome/caihuize/miniconda3/envs/ESM3_3Di_3090/bin/python}"
export FOLDSEEK_BIN="${FOLDSEEK_BIN:-/hpcfs/fhome/caihuize/.local/bin/foldseek}"
export MMSEQS_BIN="${MMSEQS_BIN:-/hpcfs/fhome/caihuize/.local/bin/mmseqs}"
bash scripts/run.sh cluster
bash scripts/run.sh split
