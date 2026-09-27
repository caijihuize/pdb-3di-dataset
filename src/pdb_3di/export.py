"""Export split Parquet, FASTA and evaluation mmCIF files."""

from __future__ import annotations

import argparse
import csv
import json
import shutil
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

from pdb_3di.source import location


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--records", type=Path, required=True)
    p.add_argument("--manifest-dir", type=Path, required=True)
    p.add_argument("--raw-dir", type=Path, required=True)
    p.add_argument("--out-dir", type=Path, required=True)
    p.add_argument("--skip-structures", action="store_true")
    a = p.parse_args()
    records = {}
    with a.records.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            records[row["id"]] = row
    a.out_dir.mkdir(parents=True, exist_ok=True)
    splits = {}
    for split in ("train", "valid", "test"):
        ids = (a.manifest_dir / f"{split}_ids.txt").read_text().splitlines()
        if len(ids) != len(set(ids)) or set(ids) - records.keys():
            raise ValueError(f"{split}: duplicate or missing record IDs")
        splits[split] = ids
        with (a.out_dir / f"{split}.fasta").open("w") as fa, \
             (a.out_dir / f"{split}_3di.fasta").open("w") as fd, \
             (a.manifest_dir / f"{split}_row_index.tsv").open("w") as index:
            index.write("row_index\tid\n")
            for i, key in enumerate(ids):
                fa.write(f">{key}\n{records[key]['aa']}\n")
                fd.write(f">{key}\n{records[key]['3di']}\n")
                index.write(f"{i}\t{key}\n")
        table = pa.table({"id": ids,
                          "sequence_3di": [records[x]["3di"] for x in ids],
                          "sequence_aa": [records[x]["aa"] for x in ids]})
        out_name = "validation" if split == "valid" else split
        pq.write_table(table, a.out_dir / f"{out_name}.parquet", compression="zstd")
    if len(set(sum(splits.values(), []))) != sum(map(len, splits.values())):
        raise ValueError("Cross-split duplicate ID")
    if not a.skip_structures:
        missing = []
        for split in ("valid", "test"):
            dest = a.out_dir / "structures" / ("validation" if split == "valid" else split)
            dest.mkdir(parents=True, exist_ok=True)
            for old in dest.glob("*.cif.gz"):
                old.unlink()
            for key in splits[split]:
                pdb_id = records[key]["pdb_id"]
                source = location(a.raw_dir, pdb_id)
                if source.is_file():
                    shutil.copy2(source, dest / source.name)
                else:
                    missing.append(key)
        (a.manifest_dir / "structure_missing_ids.txt").write_text(
            "".join(f"{x}\n" for x in missing))
        if missing:
            raise ValueError(f"Missing {len(missing)} evaluation mmCIF files")
    (a.out_dir / "schema.json").write_text(json.dumps({
        "columns": ["id", "sequence_3di", "sequence_aa"],
        "sample_id": "<PDB_ID>_<Foldseek chain ID>",
        "evaluation_structures": "full deposited mmCIF; target chain is in records.tsv",
        "split_counts": {key: len(ids) for key, ids in splits.items()},
    }, indent=2) + "\n")
    print({key: len(ids) for key, ids in splits.items()})


if __name__ == "__main__":
    main()
