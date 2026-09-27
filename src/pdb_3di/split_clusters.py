"""Merge MMseqs clusters by PDB entry/exact sequence and make deterministic splits."""

from __future__ import annotations

import argparse
import csv
import random
from collections import defaultdict
from pathlib import Path


class UnionFind:
    def __init__(self, items):
        self.parent = {x: x for x in items}

    def find(self, x):
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, a, b):
        a, b = self.find(a), self.find(b)
        if a != b:
            self.parent[max(a, b)] = min(a, b)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--records", required=True)
    p.add_argument("--clusters", required=True)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--valid-groups", type=int, default=1000)
    p.add_argument("--test-groups", type=int, default=1000)
    p.add_argument("--eval-min-length", type=int, default=30)
    p.add_argument("--eval-max-length", type=int, default=510)
    a = p.parse_args()
    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    records = {}
    with open(a.records, newline="") as h:
        for r in csv.DictReader(h, delimiter="\t"):
            r["length"] = int(r["length"])
            r["quality_score"] = float(r["quality_score"])
            records[r["id"]] = r
    uf = UnionFind(records)
    with open(a.clusters) as h:
        for line in h:
            rep, member = line.rstrip("\n").split("\t")[:2]
            if rep in records and member in records:
                uf.union(rep, member)
    by_pdb, by_seq = defaultdict(list), defaultdict(list)
    for pid, row in records.items():
        by_pdb[row["pdb_id"]].append(pid)
        by_seq[row["aa"]].append(pid)
    for groups in (by_pdb.values(), by_seq.values()):
        for members in groups:
            first = members[0]
            for member in members[1:]:
                uf.union(first, member)

    groups = defaultdict(list)
    for pid in records:
        groups[uf.find(pid)].append(pid)
    eligible = []
    for root, members in groups.items():
        choices = [pid for pid in members if a.eval_min_length <= records[pid]["length"] <= a.eval_max_length]
        if choices:
            representative = min(choices, key=lambda pid: (-records[pid]["quality_score"], pid))
            eligible.append((root, representative))
    if len(eligible) < a.valid_groups + a.test_groups:
        raise SystemExit(f"Only {len(eligible)} eligible groups; need {a.valid_groups + a.test_groups}")
    eligible.sort()
    random.Random(a.seed).shuffle(eligible)
    test_groups = {root for root, _ in eligible[:a.test_groups]}
    valid_groups = {root for root, _ in eligible[a.test_groups:a.test_groups+a.valid_groups]}
    # Emit one evaluation representative per held-out group. Keep other members excluded.
    test_ids = [rep for root, rep in eligible[:a.test_groups]]
    valid_ids = [rep for root, rep in eligible[a.test_groups:a.test_groups+a.valid_groups]]
    heldout = test_groups | valid_groups

    # De-duplicate exact sequences in training using the highest-confidence, then ID rule.
    train_by_seq = {}
    for root, members in groups.items():
        if root in heldout:
            continue
        for pid in members:
            seq = records[pid]["aa"]
            old = train_by_seq.get(seq)
            if old is None or (-records[pid]["quality_score"], pid) < (-records[old]["quality_score"], old):
                train_by_seq[seq] = pid
    train_ids = sorted(train_by_seq.values())
    split_sets = {"train": train_ids, "valid": valid_ids, "test": test_ids}
    for name, ids in split_sets.items():
        (out / f"{name}_ids.txt").write_text("".join(f"{pid}\n" for pid in ids))

    with open(out / "cluster_membership.tsv", "w") as h:
        h.write("id\tmerged_group\tsplit\tselected_representative\n")
        rep_to_split = {pid: name for name in ("valid", "test") for pid in split_sets[name]}
        for root, members in sorted(groups.items()):
            split = "heldout_excluded" if root in heldout else "train"
            for pid in sorted(members):
                h.write(f"{pid}\t{root}\t{rep_to_split.get(pid, split)}\t{int(pid in rep_to_split)}\n")
    selected = set(train_ids + valid_ids + test_ids)
    with open(out / "excluded.tsv", "w") as h:
        h.write("id\treason\n")
        for root, members in sorted(groups.items()):
            for pid in sorted(members):
                if pid not in selected:
                    reason = "heldout_group_nonrepresentative" if root in heldout else "exact_aa_duplicate"
                    h.write(f"{pid}\t{reason}\n")
    with open(out / "split_counts.tsv", "w") as h:
        h.write(f"candidate_records\t{len(records)}\nmerged_groups\t{len(groups)}\n")
        h.write(f"train\t{len(train_ids)}\nvalid\t{len(valid_ids)}\ntest\t{len(test_ids)}\n")
        h.write(f"evaluation_groups\t{len(eligible)}\nseed\t{a.seed}\n")


if __name__ == "__main__":
    main()
