#!/usr/bin/env bash
#SBATCH --job-name=pdb3di_publish_hf
#SBATCH --partition=qfree
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=03:00:00
#SBATCH --output=runs/publish_hf_%j.out
#SBATCH --error=runs/publish_hf_%j.err
set -euo pipefail
ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
cd "${ROOT}"
RELEASE="${ROOT}/release/pdb_v1"
HF_BIN="${HF_BIN:-/hpcfs/fhome/caihuize/.local/bin/hf}"
[[ -s "${ROOT}/reports/pdb_v1/validation.json" ]]
[[ -s "${ROOT}/reports/pdb_v1/summary.json" ]]
[[ -s "${RELEASE}/checksums/SHA256SUMS" ]]
(
  cd "${RELEASE}"
  sha256sum -c checksums/SHA256SUMS
)
"${HF_BIN}" auth whoami
"${HF_BIN}" repo create caijihuize/pdb-3di-dataset --repo-type dataset --exist-ok
"${HF_BIN}" upload caijihuize/pdb-3di-dataset "${RELEASE}" . \
  --repo-type dataset --commit-message "Publish validated experimental PDB AA-3Di v1"
