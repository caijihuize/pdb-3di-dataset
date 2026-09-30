"""Translate stable structure IDs to numeric Foldseek DB keys."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


MODEL_ID = re.compile(r"^(?P<pdb>[^_]{4})_MODEL_(?P<model>\d+)_(?P<chain>.+)$")


def normalize_id(value: str) -> str:
    """Normalize only the case-insensitive PDB entry part of a Foldseek ID."""
    pdb_id, separator, rest = value.partition("_")
    if not separator or len(pdb_id) != 4 or not rest:
        raise ValueError(f"Invalid Foldseek structure ID: {value}")
    return f"{pdb_id.upper()}_{rest}"


def components(value: str) -> tuple[str, str, str]:
    match = MODEL_ID.match(value)
    if match:
        return match.group("pdb").upper(), match.group("model"), match.group("chain")
    pdb_id, _, chain = value.partition("_")
    return pdb_id.upper(), "", chain


def header_keys(header_db: Path, header_index: Path, wanted: set[str]) -> dict[str, str]:
    """Read exact, model-aware IDs from the Foldseek header database."""
    found = {}
    with header_db.open("rb") as headers, header_index.open() as index:
        for line in index:
            key, offset, length = line.split()[:3]
            headers.seek(int(offset))
            raw_id = headers.read(int(length)).split(None, 1)[0].rstrip(b"\0").decode(
                "utf-8", errors="strict"
            )
            stable_id = normalize_id(raw_id)
            if stable_id in wanted:
                if stable_id in found:
                    raise ValueError(f"Duplicate Foldseek header ID: {stable_id}")
                found[stable_id] = key
    return found


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--header-db", type=Path, required=True)
    p.add_argument("--header-index", type=Path, required=True)
    p.add_argument("--ids", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--mapping", type=Path, required=True)
    p.add_argument("--split", required=True)
    a = p.parse_args()
    ordered = [line.strip() for line in a.ids.read_text().splitlines() if line.strip()]
    if len(ordered) != len(set(ordered)):
        raise ValueError(f"Duplicate IDs in {a.ids}")
    wanted = set(ordered)
    found = header_keys(a.header_db, a.header_index, wanted)
    missing = sorted(wanted - set(found))
    if missing:
        raise SystemExit(
            f"{len(missing)} IDs not found in Foldseek header database; examples: {missing[:5]}"
        )
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text("".join(f"{found[x]}\n" for x in ordered))
    a.mapping.parent.mkdir(parents=True, exist_ok=True)
    with a.mapping.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["row_index", "split", "id", "pdb_id", "model_id", "chain_id",
                         "foldseek_db_key"])
        for row_index, stable_id in enumerate(ordered):
            pdb_id, model_id, chain_id = components(stable_id)
            writer.writerow([row_index, a.split, stable_id, pdb_id, model_id, chain_id,
                             found[stable_id]])


if __name__ == "__main__":
    main()
