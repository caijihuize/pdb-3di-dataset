"""Split FASTA into bounded target shards for uncapped cross-set searches."""

from __future__ import annotations

import argparse
from pathlib import Path


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--records-per-shard", type=int, default=50000)
    a = p.parse_args()
    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    shard = -1
    count = 0
    handle = None
    try:
        with open(a.input) as source:
            for line in source:
                if line.startswith(">"):
                    if count % a.records_per_shard == 0:
                        if handle:
                            handle.close()
                        shard += 1
                        handle = (out / f"target_{shard:03d}.fasta").open("w")
                    count += 1
                if handle is None:
                    raise ValueError("FASTA must begin with a header")
                handle.write(line)
    finally:
        if handle:
            handle.close()
    print(f"records={count} shards={shard + 1} records_per_shard={a.records_per_shard}")


if __name__ == "__main__":
    main()
