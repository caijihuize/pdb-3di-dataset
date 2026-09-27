#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PYTHON="${PDB3DI_PYTHON:-python}"
export PYTHONPATH="${ROOT}/src${PYTHONPATH:+:${PYTHONPATH}}"
for ROUND in $(seq -w 1 20); do
  export AUDIT_DIR="${ROOT}/data/interim/sequence_round_${ROUND}"
  SKIP_STRUCTURES=1 bash "${ROOT}/scripts/run.sh" export
  bash "${ROOT}/scripts/run.sh" sequence-audit
  bash "${ROOT}/scripts/run.sh" repair
  STATUS="$("${PYTHON}" -c 'import json,sys; print(json.load(open(sys.argv[1]))["status"])' \
    "${ROOT}/manifests/pdb_v1/repair_state.json")"
  if [[ "${STATUS}" == "sequence_audit_passed" ]]; then break; fi
done
if [[ "${STATUS}" != "sequence_audit_passed" ]]; then
  echo "Sequence audit did not pass in 20 rounds" >&2
  exit 1
fi
bash "${ROOT}/scripts/run.sh" export
bash "${ROOT}/scripts/run.sh" structure-audit
"${PYTHON}" -m pdb_3di.validate --root "${ROOT}"
"${PYTHON}" -m pdb_3di.report --root "${ROOT}" \
  --audit-dir "${AUDIT_DIR}" --structure-dir "${ROOT}/data/interim/structure"
"${PYTHON}" -m pdb_3di.prepare_release --root "${ROOT}" \
  --out-dir "${ROOT}/release/pdb_v1"
