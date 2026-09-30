#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${ROOT}"
mkdir -p runs
ANNOTATE_JOB=$(sbatch --parsable scripts/submit_annotations.sh)
STRUCTURE_JOB=$(sbatch --parsable scripts/submit_structure_clusters.sh)
VERSIONS_JOB=$(sbatch --parsable --dependency="afterok:${ANNOTATE_JOB}:${STRUCTURE_JOB}" scripts/submit_build_versions.sh)
printf 'annotations_job\t%s\nstructure_clusters_job\t%s\nversions_job\t%s\n' \
  "${ANNOTATE_JOB}" "${STRUCTURE_JOB}" "${VERSIONS_JOB}" | tee runs/v2_pipeline_jobs.tsv
