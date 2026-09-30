"""Validate enhanced dataset annotations, splits, row order and mappings."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import pyarrow.parquet as pq

from pdb_3di.export_versions import SPLITS, VERSIONS


EXPECTED = {
    "pdb_v1_1": {"train": 134637, "valid": 1000, "test": 1000},
    "pdb_v2_sequence": {"valid": 1000, "test": 1000},
    "pdb_v2_structure": {"valid": 500, "test": 500},
}


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = {"status": "passed", "versions": {}, "mapping_checksums": {}}
    candidate_count = sum(1 for _ in (args.root / "data/interim/annotations/record_annotations.tsv").open()) - 1
    error_count = sum(1 for _ in (args.root / "reports/pdb_v1_1/annotation_errors.tsv").open()) - 1
    if candidate_count != 837946 or error_count != 0:
        raise RuntimeError(f"Annotation validation failed: records={candidate_count} errors={error_count}")
    for version in VERSIONS:
        manifest = args.root / "manifests" / version
        processed = args.root / "data/processed" / version
        seen = set()
        counts = {}
        for split in SPLITS:
            ids = (manifest / f"{split}_ids.txt").read_text().splitlines()
            if len(ids) != len(set(ids)) or seen.intersection(ids):
                raise RuntimeError(f"{version}/{split}: duplicate or cross-split IDs")
            seen.update(ids)
            counts[split] = len(ids)
            with (manifest / f"{split}_structure_ids.tsv").open(newline="") as handle:
                mapping = list(csv.DictReader(handle, delimiter="\t"))
            if [row["id"] for row in mapping] != ids or [int(row["row_index"]) for row in mapping] != list(range(len(ids))):
                raise RuntimeError(f"{version}/{split}: structure mapping order mismatch")
            parquet_name = "validation" if split == "valid" else split
            table_ids = pq.read_table(processed / f"{parquet_name}.parquet", columns=["id"])["id"].to_pylist()
            if table_ids != ids:
                raise RuntimeError(f"{version}/{split}: Parquet row order mismatch")
        for split, expected in EXPECTED.get(version, {}).items():
            if counts[split] != expected:
                raise RuntimeError(f"{version}/{split}: {counts[split]} != {expected}")
        report["versions"][version] = counts
    if ((args.root / "manifests/pdb_v2_sequence/train_ids.txt").read_bytes()
            == (args.root / "manifests/pdb_v2_structure/train_ids.txt").read_bytes()):
        raise RuntimeError("Sequence and structure isolated training splits are identical")
    sources = args.root / "data/raw/mappings/SHA256SUMS"
    for line in sources.read_text().splitlines():
        expected, filename = line.split(maxsplit=1)
        filename = filename.lstrip("* ")
        actual = digest(sources.parent / filename)
        if actual != expected:
            raise RuntimeError(f"Mapping checksum mismatch: {filename}")
        report["mapping_checksums"][filename] = actual
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, sort_keys=True))


if __name__ == "__main__":
    main()
