#!/usr/bin/env bash
#SBATCH --job-name=pdb3di_versions
#SBATCH --partition=qgpu_3090
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --gres=gpu:1
#SBATCH --time=1-00:00:00
#SBATCH --output=runs/versions_%j.out
#SBATCH --error=runs/versions_%j.err
set -euo pipefail
ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
cd "${ROOT}"
PYTHON="${PDB3DI_PYTHON:-/hpcfs/fhome/caihuize/miniconda3/envs/ESM3_3Di_3090/bin/python}"
export PYTHONPATH="${ROOT}/src${PYTHONPATH:+:${PYTHONPATH}}"
"${PYTHON}" -m pdb_3di.build_versions \
  --records data/interim/candidates/records.tsv \
  --annotations data/interim/annotations/record_annotations.tsv \
  --sequence-clusters data/interim/mmseqs/clusters.tsv \
  --structure-clusters data/interim/structure_clusters/clusters.tsv \
  --structure-key-map data/interim/structure_clusters/candidate_structure_ids.tsv \
  --cath-domain-list data/raw/mappings/cath-domain-list.txt \
  --v1-manifests manifests/pdb_v1 --output-root manifests \
  --seed 42 --valid-groups 1000 --test-groups 1000
