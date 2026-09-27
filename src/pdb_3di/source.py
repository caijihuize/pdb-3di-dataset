"""Freeze an experimental-protein PDB ID inventory and download its mmCIF files."""

from __future__ import annotations

import argparse
import concurrent.futures
import csv
import gzip
import hashlib
import json
import os
import time
from datetime import datetime, timezone
from pathlib import Path

import requests

SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
FILES = "https://files.rcsb.org/download"


def discover(out: Path, page_size: int) -> None:
    out.mkdir(parents=True, exist_ok=True)
    query = {
        "query": {"type": "terminal", "service": "text", "parameters": {
            "attribute": "rcsb_entry_info.polymer_entity_count_protein",
            "operator": "greater", "value": 0}},
        "return_type": "entry",
        "request_options": {"results_content_type": ["experimental"],
                            "paginate": {"start": 0, "rows": page_size},
                            "sort": [{"sort_by": "rcsb_entry_container_identifiers.entry_id",
                                      "direction": "asc"}]},
    }
    session = requests.Session()
    ids: list[str] = []
    total = None
    for start in range(0, 10**9, page_size):
        query["request_options"]["paginate"]["start"] = start
        for attempt in range(5):
            try:
                response = session.post(SEARCH, json=query, timeout=120)
                response.raise_for_status()
                payload = response.json()
                break
            except (requests.RequestException, ValueError):
                if attempt == 4:
                    raise
                time.sleep(2 ** attempt)
        if total is None:
            total = payload["total_count"]
        elif total != payload["total_count"]:
            raise RuntimeError("PDB holdings changed during pagination; rerun discovery")
        ids.extend(row["identifier"].upper() for row in payload["result_set"])
        print(f"discovered {len(ids)}/{total}", flush=True)
        if len(ids) >= total:
            break
    if len(ids) != total or len(ids) != len(set(ids)):
        raise RuntimeError(f"Incomplete/duplicate ID inventory: {len(ids)} rows, {total} expected")
    ids.sort()
    (out / "source_ids.txt").write_text("".join(f"{x}\n" for x in ids))
    query["request_options"]["paginate"] = {"start": 0, "rows": page_size}
    (out / "source_query.json").write_text(json.dumps({
        "retrieved_utc": datetime.now(timezone.utc).isoformat(), "count": len(ids),
        "query": query, "coordinate_format": "PDBx/mmCIF.gz",
        "download_url_template": FILES + "/{pdb_id}.cif.gz",
    }, indent=2) + "\n")


def location(root: Path, pdb_id: str) -> Path:
    return root / pdb_id[1:3].lower() / f"{pdb_id.lower()}.cif.gz"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fetch_one(pdb_id: str, root: Path) -> tuple[str, str, int, str]:
    target = location(root, pdb_id)
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.is_file():
        try:
            with gzip.open(target, "rb") as handle:
                if not handle.read(64).startswith(b"data_"):
                    raise ValueError("invalid mmCIF header")
                for _ in iter(lambda: handle.read(1024 * 1024), b""):
                    pass
            return pdb_id, "ok", target.stat().st_size, sha256(target)
        except (OSError, ValueError):
            target.unlink()
    part = target.with_suffix(target.suffix + f".{os.getpid()}.part")
    for attempt in range(5):
        try:
            with requests.get(f"{FILES}/{pdb_id}.cif.gz", stream=True, timeout=120) as response:
                response.raise_for_status()
                with part.open("wb") as handle:
                    for chunk in response.iter_content(1024 * 1024):
                        handle.write(chunk)
            with gzip.open(part, "rb") as handle:
                if not handle.read(64).startswith(b"data_"):
                    raise ValueError("invalid mmCIF header")
                for _ in iter(lambda: handle.read(1024 * 1024), b""):
                    pass
            part.replace(target)
            return pdb_id, "ok", target.stat().st_size, sha256(target)
        except (requests.RequestException, OSError, ValueError) as exc:
            part.unlink(missing_ok=True)
            if attempt == 4:
                return pdb_id, f"error:{type(exc).__name__}:{exc}", 0, ""
            time.sleep(2 ** attempt)
    raise AssertionError("unreachable")


def download(ids_file: Path, root: Path, out: Path, workers: int, limit: int,
             shard_index: int, shard_count: int) -> None:
    ids = [x.strip().upper() for x in ids_file.read_text().splitlines() if x.strip()]
    if limit:
        ids = ids[:limit]
    if not 0 <= shard_index < shard_count:
        raise ValueError("shard_index must be in [0, shard_count)")
    ids = [x for i, x in enumerate(ids) if i % shard_count == shard_index]
    out.parent.mkdir(parents=True, exist_ok=True)
    failures = 0
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool, out.open("w") as handle:
        handle.write("pdb_id\tstatus\tbytes\tsha256\n")
        for index, (pdb_id, status, size, checksum) in enumerate(
            pool.map(lambda x: fetch_one(x, root), ids), 1
        ):
            handle.write(f"{pdb_id}\t{status}\t{size}\t{checksum}\n")
            failures += status != "ok"
            if index % 1000 == 0:
                handle.flush()
                print(f"downloaded/checked {index}/{len(ids)}", flush=True)
    print(f"wrote {out}")
    if failures:
        raise SystemExit(f"{failures} downloads failed; rerun the shard to retry")


def verify(ids_file: Path, root: Path, logs: Path, report: Path) -> None:
    ids = set(ids_file.read_text().splitlines())
    checked = {}
    failures = []
    for log in sorted(logs.glob("downloads_*.tsv")):
        with log.open(newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                key = row["pdb_id"]
                if key in checked:
                    failures.append(f"duplicate log ID {key}")
                checked[key] = row
                if row["status"] != "ok" or not location(root, key).is_file():
                    failures.append(f"{key}: {row['status']}")
    missing = sorted(ids - checked.keys())
    extra = sorted(checked.keys() - ids)
    result = {"source_ids": len(ids), "logged_ids": len(checked),
              "missing_count": len(missing), "extra_count": len(extra),
              "failures_count": len(failures), "missing_examples": missing[:20],
              "extra_examples": extra[:20], "failure_examples": failures[:20],
              "status": "passed" if not (missing or extra or failures) else "failed"}
    report.parent.mkdir(parents=True, exist_ok=True)
    report.write_text(json.dumps(result, indent=2) + "\n")
    print(result)
    if result["status"] != "passed":
        raise SystemExit(1)


def main() -> None:
    p = argparse.ArgumentParser()
    sub = p.add_subparsers(dest="command", required=True)
    d = sub.add_parser("discover")
    d.add_argument("--out-dir", type=Path, required=True)
    d.add_argument("--page-size", type=int, default=10000)
    f = sub.add_parser("download")
    f.add_argument("--ids", type=Path, required=True)
    f.add_argument("--raw-dir", type=Path, required=True)
    f.add_argument("--log", type=Path, required=True)
    f.add_argument("--workers", type=int, default=8)
    f.add_argument("--limit", type=int, default=0, help="First N IDs for a pilot")
    f.add_argument("--shard-index", type=int, default=0)
    f.add_argument("--shard-count", type=int, default=1)
    v = sub.add_parser("verify")
    v.add_argument("--ids", type=Path, required=True)
    v.add_argument("--raw-dir", type=Path, required=True)
    v.add_argument("--log-dir", type=Path, required=True)
    v.add_argument("--report", type=Path, required=True)
    a = p.parse_args()
    if a.command == "discover":
        discover(a.out_dir, a.page_size)
    elif a.command == "download":
        download(a.ids, a.raw_dir, a.log, a.workers, a.limit,
                 a.shard_index, a.shard_count)
    else:
        verify(a.ids, a.raw_dir, a.log_dir, a.report)


if __name__ == "__main__":
    main()
