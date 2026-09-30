#!/usr/bin/env bash
#SBATCH --job-name=pdb3di_hf_v2
#SBATCH --partition=qfree
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=03:00:00
#SBATCH --output=runs/publish_enhanced_hf_%j.out
#SBATCH --error=runs/publish_enhanced_hf_%j.err
set -euo pipefail
ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
cd "${ROOT}"
PYTHON="${PDB3DI_PYTHON:-/hpcfs/fhome/caihuize/miniconda3/envs/ESM3_3Di_3090/bin/python}"
HF_BIN="${HF_BIN:-/hpcfs/fhome/caihuize/.local/bin/hf}"
export PYTHONPATH="${ROOT}/src${PYTHONPATH:+:${PYTHONPATH}}"
"${PYTHON}" -m pdb_3di.prepare_enhanced_release \
  --root "${ROOT}" --output "${ROOT}/release/enhanced"
(cd release/enhanced && sha256sum -c SHA256SUMS)
"${HF_BIN}" auth whoami
"${HF_BIN}" upload caijihuize/pdb-3di-dataset release/enhanced versions \
  --repo-type dataset --commit-message "Publish validated enhanced PDB AA-3Di versions"
