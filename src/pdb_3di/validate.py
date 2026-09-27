"""Validate exported rows, split IDs and structure references."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import pyarrow.parquet as pq


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--root", type=Path, required=True)
    a = p.parse_args()
    root = a.root
    manifest = root / "manifests/pdb_v1"
    processed = root / "data/processed/pdb_v1"
    with (root / "data/interim/candidates/records.tsv").open(newline="") as handle:
        records = {row["id"]: row for row in csv.DictReader(handle, delimiter="\t")}
    all_ids = set()
    checks = {"split_counts": {}, "status": "passed"}
    for split in ("train", "valid", "test"):
        ids = (manifest / f"{split}_ids.txt").read_text().splitlines()
        name = "validation" if split == "valid" else split
        table = pq.read_table(processed / f"{name}.parquet").to_pydict()
        assert table["id"] == ids, f"{split} row order mismatch"
        assert len(set(ids)) == len(ids), f"{split} duplicate IDs"
        assert not (set(ids) & all_ids), f"{split} overlaps another split"
        all_ids.update(ids)
        for i, key in enumerate(ids):
            row = records[key]
            assert table["sequence_aa"][i] == row["aa"], f"{key} AA mismatch"
            assert table["sequence_3di"][i] == row["3di"], f"{key} 3Di mismatch"
            assert len(row["aa"]) == len(row["3di"]), f"{key} length mismatch"
            if split != "train":
                source = processed / "structures" / name / f"{row['pdb_id'].lower()}.cif.gz"
                assert source.is_file(), f"{key} missing evaluation structure"
        checks["split_counts"][split] = len(ids)
    path = root / "reports/pdb_v1/validation.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(checks, indent=2) + "\n")
    print(path)


if __name__ == "__main__":
    main()
