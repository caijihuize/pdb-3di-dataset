"""Build annotated v1.1 and sequence, structure and chronological v2 manifests."""

from __future__ import annotations

import argparse
import csv
import random
import shutil
from collections import defaultdict
from pathlib import Path


class UnionFind:
    def __init__(self, items):
        self.parent = {item: item for item in items}

    def find(self, item):
        parent = self.parent
        while parent[item] != item:
            parent[item] = parent[parent[item]]
            item = parent[item]
        return item

    def union(self, left, right):
        left, right = self.find(left), self.find(right)
        if left != right:
            self.parent[max(left, right)] = min(left, right)

    def copy(self):
        other = object.__new__(UnionFind)
        other.parent = self.parent.copy()
        return other


def read_tsv(path: Path) -> list[dict]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def unions_from_tsv(uf: UnionFind, path: Path, key_map: dict[str, str] | None = None) -> None:
    if not path.is_file():
        raise FileNotFoundError(path)
    with path.open() as handle:
        for line in handle:
            columns = line.rstrip("\n").split("\t")
            if len(columns) < 2:
                continue
            left, right = columns[:2]
            if key_map is not None:
                left, right = key_map.get(left, left), key_map.get(right, right)
            if left in uf.parent and right in uf.parent:
                uf.union(left, right)


def union_buckets(uf: UnionFind, buckets) -> None:
    for members in buckets:
        members = list(members)
        if members:
            first = members[0]
            for member in members[1:]:
                uf.union(first, member)


def groups(uf: UnionFind, ids) -> dict[str, list[str]]:
    output = defaultdict(list)
    for record_id in ids:
        output[uf.find(record_id)].append(record_id)
    return output


def rank(record: dict, annotation: dict) -> tuple:
    missing = int(annotation.get("internal_missing_residues") or 10**9)
    completeness = float(annotation.get("backbone_complete_fraction") or 0)
    resolution = float(record.get("resolution") or 999)
    return (missing, -completeness, resolution, record["id"])


def usable(record: dict, annotation: dict) -> bool:
    return (float(annotation.get("three_di_max_token_fraction") or 1) <= 0.95
            and float(annotation.get("backbone_complete_fraction") or 0) >= 0.90
            and int(annotation.get("internal_missing_residues") or 10**9) <= max(10, int(record["length"]) // 10))


def write_ids(directory: Path, split_ids: dict[str, list[str]]) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    for split, ids in split_ids.items():
        (directory / f"{split}_ids.txt").write_text("".join(f"{x}\n" for x in ids))


def split_groups(name: str, uf: UnionFind, records: dict, annotations: dict, output: Path,
                 seed: int, valid_count: int, test_count: int) -> None:
    grouped = groups(uf, records)
    eligible = []
    for root, members in grouped.items():
        choices = [x for x in members if usable(records[x], annotations[x])
                   and 30 <= int(records[x]["length"]) <= 510]
        if choices:
            eligible.append((root, min(choices, key=lambda x: rank(records[x], annotations[x]))))
    if len(eligible) < valid_count + test_count:
        raise RuntimeError(f"{name}: only {len(eligible)} eligible groups")
    eligible.sort()
    random.Random(seed).shuffle(eligible)
    test_roots = {root for root, _ in eligible[:test_count]}
    valid_roots = {root for root, _ in eligible[test_count:test_count + valid_count]}
    heldout = test_roots | valid_roots
    test = sorted(rep for _, rep in eligible[:test_count])
    valid = sorted(rep for _, rep in eligible[test_count:test_count + valid_count])
    best_by_sequence = {}
    for root, members in grouped.items():
        if root in heldout:
            continue
        for record_id in members:
            if not usable(records[record_id], annotations[record_id]):
                continue
            sequence = records[record_id]["aa"]
            previous = best_by_sequence.get(sequence)
            if previous is None or rank(records[record_id], annotations[record_id]) < rank(records[previous], annotations[previous]):
                best_by_sequence[sequence] = record_id
    train = sorted(best_by_sequence.values())
    split_ids = {"train": train, "valid": valid, "test": test}
    write_ids(output, split_ids)
    selected = set(train + valid + test)
    with (output / "cluster_membership.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["id", "group_id", "group_size", "sampling_weight", "split", "selected"])
        for root, members in sorted(grouped.items()):
            group_id = f"{name}:{root}"
            selected_split = "test" if root in test_roots else "valid" if root in valid_roots else "train"
            for record_id in sorted(members):
                writer.writerow([record_id, group_id, len(members), f"{1 / len(members):.8g}",
                                 selected_split, int(record_id in selected)])
    with (output / "split_counts.tsv").open("w") as handle:
        handle.write(f"candidates\t{len(records)}\nusable\t{sum(usable(records[x], annotations[x]) for x in records)}\n")
        handle.write(f"groups\t{len(grouped)}\ntrain\t{len(train)}\nvalid\t{len(valid)}\ntest\t{len(test)}\nseed\t{seed}\n")


def chronological(records: dict, annotations: dict, sequence_uf: UnionFind,
                  structure_uf: UnionFind, output: Path) -> None:
    canonical = {}
    for record_id, record in records.items():
        annotation = annotations[record_id]
        if not usable(record, annotation) or not annotation.get("release_date"):
            continue
        old = canonical.get(record["aa"])
        if old is None or rank(record, annotation) < rank(records[old], annotations[old]):
            canonical[record["aa"]] = record_id
    periods = {"train": [], "valid": [], "test": []}
    for record_id in canonical.values():
        date = annotations[record_id]["release_date"]
        if date <= "2024-12-31":
            periods["train"].append(record_id)
        elif date <= "2025-12-31":
            periods["valid"].append(record_id)
        elif date <= "2026-09-30":
            periods["test"].append(record_id)
    for ids in periods.values():
        ids.sort()
    write_ids(output, periods)
    train_sequence = {sequence_uf.find(x) for x in periods["train"]}
    train_structure = {structure_uf.find(x) for x in periods["train"]}
    with (output / "novelty.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["id", "split", "release_date", "sequence_novel", "structure_novel"])
        for split in ("valid", "test"):
            for record_id in periods[split]:
                writer.writerow([record_id, split, annotations[record_id]["release_date"],
                                 int(sequence_uf.find(record_id) not in train_sequence),
                                 int(structure_uf.find(record_id) not in train_structure)])
    with (output / "split_counts.tsv").open("w") as handle:
        for split, ids in periods.items():
            handle.write(f"{split}\t{len(ids)}\n")


def v11(source: Path, annotations_path: Path, output: Path) -> None:
    output.mkdir(parents=True, exist_ok=True)
    selected = set()
    for split in ("train", "valid", "test"):
        ids = (source / f"{split}_ids.txt").read_text().splitlines()
        selected.update(ids)
        shutil.copy2(source / f"{split}_ids.txt", output / f"{split}_ids.txt")
        mapping = source / f"{split}_structure_ids.tsv"
        if mapping.is_file():
            shutil.copy2(mapping, output / mapping.name)
    with annotations_path.open(newline="") as source_handle, \
         (output / "record_annotations.tsv").open("w", newline="") as target_handle:
        reader = csv.DictReader(source_handle, delimiter="\t")
        writer = csv.DictWriter(target_handle, fieldnames=reader.fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        count = 0
        for row in reader:
            if row["id"] in selected:
                writer.writerow(row)
                count += 1
    if count != len(selected):
        raise RuntimeError(f"pdb_v1.1 annotation mismatch: {count} for {len(selected)} IDs")
    (output / "VERSION.txt").write_text("pdb_v1.1\nSame AA/3Di rows and split IDs as pdb_v1; adds quality, date and identifier annotations.\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--records", type=Path, required=True)
    parser.add_argument("--annotations", type=Path, required=True)
    parser.add_argument("--sequence-clusters", type=Path, required=True)
    parser.add_argument("--structure-clusters", type=Path, required=True)
    parser.add_argument("--structure-key-map", type=Path, required=True)
    parser.add_argument("--v1-manifests", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--valid-groups", type=int, default=1000)
    parser.add_argument("--test-groups", type=int, default=1000)
    args = parser.parse_args()
    records = {row["id"]: row for row in read_tsv(args.records)}
    annotations = {row["id"]: row for row in read_tsv(args.annotations)}
    missing = records.keys() - annotations.keys()
    if missing:
        raise RuntimeError(f"Missing annotations for {len(missing)} records")
    sequence_uf = UnionFind(records)
    unions_from_tsv(sequence_uf, args.sequence_clusters)
    by_pdb, by_sequence, by_uniprot = defaultdict(list), defaultdict(list), defaultdict(list)
    for record_id, record in records.items():
        by_pdb[record["pdb_id"]].append(record_id)
        by_sequence[record["aa"]].append(record_id)
        accessions = annotations[record_id].get("sifts_uniprot_ids") or annotations[record_id].get("native_uniprot_ids")
        for accession in accessions.split(";"):
            if accession:
                by_uniprot[accession].append(record_id)
    union_buckets(sequence_uf, by_pdb.values())
    union_buckets(sequence_uf, by_sequence.values())
    union_buckets(sequence_uf, by_uniprot.values())
    structure_uf = sequence_uf.copy()
    key_map = {}
    with args.structure_key_map.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            key_map[row["foldseek_db_key"]] = row["id"]
    unions_from_tsv(structure_uf, args.structure_clusters, key_map)
    v11(args.v1_manifests, args.annotations, args.output_root / "pdb_v1_1")
    split_groups("sequence", sequence_uf, records, annotations,
                 args.output_root / "pdb_v2_sequence", args.seed, args.valid_groups, args.test_groups)
    split_groups("structure", structure_uf, records, annotations,
                 args.output_root / "pdb_v2_structure", args.seed, args.valid_groups, args.test_groups)
    chronological(records, annotations, sequence_uf, structure_uf,
                  args.output_root / "pdb_v2_time")


if __name__ == "__main__":
    main()
