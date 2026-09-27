"""Deterministically remove train hits and replace valid/test sequence conflicts."""

from __future__ import annotations

import argparse
import csv
import json
import random
from collections import defaultdict
from pathlib import Path


def read_ids(path: Path) -> list[str]:
    return path.read_text().splitlines()


def read_hits(path: Path) -> list[tuple[str, str]]:
    if not path.exists():
        return []
    hits = []
    with path.open() as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if len(row) >= 8:
                try:
                    identity, qcov, tcov = float(row[2]), float(row[4]), float(row[5])
                except ValueError:
                    continue
                if identity >= 0.30 and qcov >= 0.80 and tcov >= 0.80:
                    hits.append((row[0], row[1]))
    return hits


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--records", required=True)
    p.add_argument("--manifest-dir", required=True)
    p.add_argument("--audit-dir", required=True)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--valid-count", type=int, default=1000)
    p.add_argument("--test-count", type=int, default=1000)
    p.add_argument("--min-length", type=int, default=30)
    p.add_argument("--max-length", type=int, default=510)
    a = p.parse_args()
    manifests, audit = Path(a.manifest_dir), Path(a.audit_dir)
    state_path = manifests / "repair_state.json"
    membership_path = manifests / "cluster_membership.tsv"

    id_group: dict[str, str] = {}
    initial_split: dict[str, str] = {}
    group_members: dict[str, list[str]] = defaultdict(list)
    with membership_path.open() as h:
        for row in csv.DictReader(h, delimiter="\t"):
            pid, group = row["id"], row["merged_group"]
            id_group[pid] = group
            if row["split"] == "train":
                initial_split[group] = "train"
            elif row["selected_representative"] == "1" and row["split"] in {"valid", "test"}:
                initial_split[group] = row["split"]
            group_members[group].append(pid)

    if state_path.exists():
        state = json.loads(state_path.read_text())
    else:
        initial_eval = {g: s for g, s in initial_split.items() if s in {"valid", "test"}}
        initial_train = sorted(g for g, s in initial_split.items() if s == "train")
        state = {
            "seed": a.seed,
            "initial_eval_groups": initial_eval,
            "initial_train_groups_order": initial_train,
            "blocked_groups": [],
            "audit_rounds": [],
        }
        random.Random(a.seed).shuffle(state["initial_train_groups_order"])

    with open(a.records, newline="") as h:
        records = {r["id"]: r for r in csv.DictReader(h, delimiter="\t")}
    train = set(read_ids(manifests / "train_ids.txt"))
    valid = read_ids(manifests / "valid_ids.txt")
    test = read_ids(manifests / "test_ids.txt")
    valid_set, test_set = set(valid), set(test)
    blocked = set(state["blocked_groups"])

    train_hits = read_hits(audit / "valid_train.tsv") + read_hits(audit / "test_train.tsv")
    remove_train_groups = {id_group[target] for _, target in train_hits if target in id_group and target in train}
    for group in remove_train_groups:
        blocked.add(group)
        train.difference_update(group_members[group])

    eval_hits = read_hits(audit / "valid_test.tsv")
    remove_test_groups = {id_group[target] for _, target in eval_hits if target in id_group and target in test_set}
    if not train_hits and not eval_hits:
        state["audit_rounds"].append({"train_hits": 0, "valid_test_hits": 0, "train_records_removed": 0})
        state["blocked_groups"] = sorted(blocked)
        state["status"] = "sequence_audit_passed"
        state_path.write_text(json.dumps(state, indent=2) + "\n")
        print("No cross-split sequence hits at identity >=0.30 and bidirectional coverage >=0.80")
        return

    removed_train_count = len(remove_train_groups)
    if remove_test_groups:
        test = [pid for pid in test if id_group.get(pid) not in remove_test_groups]
        blocked.update(remove_test_groups)

    current_eval_groups = {id_group[x] for x in valid + test if x in id_group}
    initial_eval_groups = set(state["initial_eval_groups"])
    # Replace removed test groups only with unused original training groups. A later
    # complete audit either accepts the candidate or removes it deterministically.
    if len(test) < a.test_count:
        selected = set(valid + test)
        for group in state["initial_train_groups_order"]:
            if group in current_eval_groups or group in blocked:
                continue
            members = group_members[group]
            choices = [pid for pid in members if a.min_length <= int(records[pid]["length"]) <= a.max_length]
            if not choices:
                continue
            rep = min(choices, key=lambda pid: (-float(records[pid]["quality_score"]), pid))
            if rep in selected:
                continue
            test.append(rep)
            selected.add(rep)
            current_eval_groups.add(group)
            train.difference_update(members)
            if len(test) == a.test_count:
                break
    if len(test) != a.test_count:
        raise SystemExit(f"Could replenish only {len(test)} of {a.test_count} test groups")

    valid = sorted(valid)
    test = sorted(test)
    train = sorted(train)
    (manifests / "train_ids.txt").write_text("".join(f"{x}\n" for x in train))
    (manifests / "valid_ids.txt").write_text("".join(f"{x}\n" for x in valid))
    (manifests / "test_ids.txt").write_text("".join(f"{x}\n" for x in test))

    updated_valid_set, updated_test_set = set(valid), set(test)
    eval_by_group = {id_group[x]: ("valid" if x in updated_valid_set else "test") for x in valid + test}
    current_train_groups = {id_group[x] for x in train}
    with membership_path.open("w") as h:
        h.write("id\tmerged_group\tsplit\tselected_representative\n")
        for group in sorted(group_members):
            for pid in sorted(group_members[group]):
                if group in eval_by_group:
                    split = eval_by_group[group] if pid in updated_valid_set | updated_test_set else "heldout_excluded"
                    selected_flag = int(pid in updated_valid_set | updated_test_set)
                elif group in blocked:
                    split, selected_flag = "sequence_audit_removed", 0
                elif group in current_train_groups:
                    split, selected_flag = "train", 0
                elif group in initial_eval_groups:
                    split, selected_flag = "heldout_excluded", 0
                else:
                    split, selected_flag = "train_excluded", 0
                h.write(f"{pid}\t{group}\t{split}\t{selected_flag}\n")

    selected_ids = set(train + valid + test)
    with (manifests / "excluded.tsv").open("w") as h:
        h.write("id\treason\n")
        for pid in sorted(records):
            if pid in selected_ids:
                continue
            group = id_group[pid]
            if group in blocked:
                reason = "sequence_audit_conflict_group_removed"
            elif group in initial_eval_groups or group in eval_by_group:
                reason = "heldout_group_nonrepresentative"
            else:
                reason = "exact_aa_duplicate_or_group_member"
            h.write(f"{pid}\t{reason}\n")

    state["blocked_groups"] = sorted(blocked)
    state["audit_rounds"].append({
        "train_hit_rows": len(train_hits),
        "valid_test_hit_rows": len(eval_hits),
        "train_groups_removed": removed_train_count,
        "train_records_removed": sum(len(group_members[g]) for g in remove_train_groups),
        "test_groups_replaced": len(remove_test_groups),
        "train_records_after": len(train),
        "valid_records_after": len(valid),
        "test_records_after": len(test),
    })
    state["status"] = "repair_pending_reaudit"
    state_path.write_text(json.dumps(state, indent=2) + "\n")
    with (manifests / "split_counts.tsv").open("w") as h:
        h.write(f"candidate_records\t{len(records)}\nmerged_groups\t{len(group_members)}\n")
        h.write(f"train\t{len(train)}\nvalid\t{len(valid)}\ntest\t{len(test)}\n")
        h.write(f"sequence_audit_blocked_groups\t{len(blocked)}\nseed\t{a.seed}\n")
    print(state["audit_rounds"][-1])


if __name__ == "__main__":
    main()
