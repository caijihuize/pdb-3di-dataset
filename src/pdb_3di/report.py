"""Record source and tool provenance plus final audit counts."""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path


def command(args: list[str]) -> str:
    try:
        result = subprocess.run(args, capture_output=True, text=True, timeout=20, check=True)
        return (result.stdout or result.stderr).strip().splitlines()[0]
    except (OSError, subprocess.SubprocessError, IndexError):
        return "unavailable"


def audit(path: Path) -> dict:
    if not path.is_file():
        return {"status": "missing"}
    rows = []
    with path.open() as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if len(row) >= 8:
                rows.append(row)
    violations = [r for r in rows if float(r[2]) >= 0.30 and
                  float(r[4]) >= 0.80 and float(r[5]) >= 0.80]
    return {"reported_alignments": len(rows), "threshold_violations": len(violations)}


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--root", type=Path, required=True)
    p.add_argument("--audit-dir", type=Path, required=True)
    p.add_argument("--structure-dir", type=Path, required=True)
    a = p.parse_args()
    manifest = a.root / "manifests/pdb_v1"
    counts = {name: len((manifest / f"{name}_ids.txt").read_text().splitlines())
              for name in ("train", "valid", "test")}
    sequence = {name: audit(a.audit_dir / f"{name}.tsv")
                for name in ("valid_train", "test_train", "valid_test")}
    structural = {}
    for name in ("valid_train", "test_train", "valid_test"):
        path = a.structure_dir / f"{name}.tsv"
        if path.is_file():
            with path.open() as handle:
                rows = [line for line in handle if line.strip()]
            structural[name] = {"reported_neighbors": len(rows)}
        else:
            structural[name] = {"status": "missing"}
    report = {
        "dataset_version": "pdb_v1", "recorded_utc": datetime.now(timezone.utc).isoformat(),
        "source": json.loads((manifest / "source_query.json").read_text()),
        "split_counts": counts, "sequence_audit": sequence, "structure_audit": structural,
        "tools": {"foldseek": command(["foldseek", "version"]),
                  "mmseqs": command(["mmseqs", "version"]),
                  "git_revision": command(["git", "-C", str(a.root), "rev-parse", "HEAD"])},
    }
    target = a.root / "reports/pdb_v1/summary.json"
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(report, indent=2) + "\n")
    print(target)


if __name__ == "__main__":
    main()
