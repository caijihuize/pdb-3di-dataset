"""Export enhanced manifests as annotated Parquet datasets and ID mappings."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq


VERSIONS = ("pdb_v1_1", "pdb_v2_sequence", "pdb_v2_structure", "pdb_v2_time")
SPLITS = ("train", "valid", "test")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--records", type=Path, required=True)
    parser.add_argument("--annotations", type=Path, required=True)
    parser.add_argument("--structure-key-map", type=Path, required=True)
    parser.add_argument("--manifest-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    args = parser.parse_args()
    selected = set()
    ids_by_version = {}
    for version in VERSIONS:
        ids_by_version[version] = {}
        for split in SPLITS:
            ids = (args.manifest_root / version / f"{split}_ids.txt").read_text().splitlines()
            ids_by_version[version][split] = ids
            selected.update(ids)
    records = {}
    with args.records.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["id"] in selected:
                records[row["id"]] = row
    annotations = {}
    with args.annotations.open(newline="") as handle:
        annotation_fields = None
        for row in csv.DictReader(handle, delimiter="\t"):
            annotation_fields = list(row)
            if row["id"] in selected:
                annotations[row["id"]] = row
    structure_keys = {}
    with args.structure_key_map.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["id"] in selected:
                structure_keys[row["id"]] = row
    missing = selected - records.keys() | selected - annotations.keys() | selected - structure_keys.keys()
    if missing:
        raise RuntimeError(f"Missing records, annotations or structure keys for {len(missing)} IDs")
    columns = ["id", "sequence_3di", "sequence_aa", "method", "method_category",
               "resolution", "deposited_length", "coverage", "length"] + [
        field for field in annotation_fields if field not in {"id", "pdb_id", "chain_id"}
    ]
    for version in VERSIONS:
        out = args.output_root / version
        out.mkdir(parents=True, exist_ok=True)
        for split, ids in ids_by_version[version].items():
            output_name = "validation" if split == "valid" else split
            data = {column: [] for column in columns}
            for record_id in ids:
                record, annotation = records[record_id], annotations[record_id]
                values = {"id": record_id, "sequence_3di": record["3di"],
                          "sequence_aa": record["aa"],
                          **{field: record[field] for field in columns if field in record},
                          **{field: annotation[field] for field in columns if field in annotation}}
                for column in columns:
                    data[column].append(values.get(column, ""))
            pq.write_table(pa.table(data), out / f"{output_name}.parquet", compression="zstd")
            mapping_path = args.manifest_root / version / f"{split}_structure_ids.tsv"
            with mapping_path.open("w", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(["row_index", "split", "id", "pdb_id", "model_id",
                                 "chain_id", "foldseek_db_key"])
                for index, record_id in enumerate(ids):
                    row = structure_keys[record_id]
                    writer.writerow([index, split, record_id, row["pdb_id"], row["model_id"],
                                     row["chain_id"], row["foldseek_db_key"]])
        (out / "schema.txt").write_text("\n".join(columns) + "\n")


if __name__ == "__main__":
    main()
