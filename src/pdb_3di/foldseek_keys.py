"""Translate stable AlphaFold IDs to numeric Foldseek DB keys for createsubdb."""

from __future__ import annotations

import argparse
from pathlib import Path


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--lookup", required=True)
    p.add_argument("--ids", required=True)
    p.add_argument("--output", required=True)
    a = p.parse_args()
    wanted = {line.strip().upper() for line in open(a.ids) if line.strip()}
    found = {}
    with open(a.lookup) as h:
        for line in h:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            ident = parts[1].upper()
            if ident in wanted:
                found[ident] = parts[0]
    missing = sorted(wanted - found.keys())
    if missing:
        raise SystemExit(f"{len(missing)} IDs not found in Foldseek lookup; examples: {missing[:5]}")
    Path(a.output).write_text("".join(f"{found[x]}\n" for x in sorted(wanted)))


if __name__ == "__main__":
    main()
