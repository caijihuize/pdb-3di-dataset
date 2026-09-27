"""Stage a verified Hugging Face dataset release without raw PDB archive files."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from pathlib import Path


def load_json(path: Path) -> dict:
    if not path.is_file():
        raise RuntimeError(f"Required report is missing: {path}")
    return json.loads(path.read_text())


def copy(source: Path, target: Path) -> None:
    if not source.is_file():
        raise RuntimeError(f"Required release file is missing: {source}")
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, target)


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def card(counts: dict[str, int], source_count: int) -> str:
    return f"""---
pretty_name: Experimental PDB AA–3Di Dataset
license: cc0-1.0
task_categories:
  - text-generation
tags:
  - biology
  - protein
  - protein-structure
configs:
  - config_name: default
    data_files:
      - split: train
        path: data/train-00000-of-00001.parquet
      - split: validation
        path: data/validation-00000-of-00001.parquet
      - split: test
        path: data/test-00000-of-00001.parquet
---

# Experimental PDB AA–3Di Dataset

Paired amino acid (AA) and Foldseek 3Di sequences derived from experimentally
determined protein structures in the wwPDB PDB Core Archive. Version `pdb_v1`
was built from a frozen RCSB Search API inventory of {source_count:,} released
experimental PDB entries containing protein polymers. This is the *source
inventory*, before method, quality and chain filtering.

| Split | Records | Evaluation mmCIF files |
| --- | ---: | ---: |
| Train | {counts['train']:,} | Not included |
| Validation | {counts['valid']:,} | {counts['valid']:,} |
| Test | {counts['test']:,} | {counts['test']:,} |

Each Parquet row contains `id`, `sequence_3di`, and `sequence_aa`. The record ID
has the form `<PDB_ID>_<Foldseek chain ID>`. ID lists in `ids/` follow Parquet
row order. The full deposited mmCIF entries for validation and test are in
`structures/`; the target chain for a record is identified by its ID and the
metadata in `provenance/records.tsv`. Other chains may be present in the file.

```python
from datasets import load_dataset

ds = load_dataset("caijihuize/pdb-3di-dataset")
row = ds["test"][0]
print(row["id"], row["sequence_3di"], row["sequence_aa"])
```

## Method

The coordinate source is the wwPDB PDB Core Archive, accessed through RCSB.
The primary version uses X-ray entries with reported resolution at most 3.0 Å
and electron microscopy entries with reported resolution at most 4.0 Å.
Only Foldseek-observed protein residues are included in paired sequences;
unresolved residues of a deposited polymer are not added as AA-only positions.
NMR and integrative/hybrid structures are outside the primary version.

Protein chains must have at least 30 observed residues, at most 20% unknown AA,
matching AA/3Di lengths, and valid Foldseek 3Di symbols. MMseqs2 clustering
uses 30% identity and 80% bidirectional coverage. Chains from the same PDB
entry and exact AA duplicates stay in the same split group. Validation and
test each select one representative from 1,000 groups, with evaluation length
30–510 residues. A final all-member sequence search and deterministic repair
found no cross-split alignment meeting the published 30% identity and 80%
bidirectional coverage threshold. Foldseek structural neighbors are reported
in `reports/` and do not control split membership.

## Reproducibility and limitations

The inventory query and date, per-file checksums, method and chain metadata,
group membership, exclusions, tool versions, audit summaries and validation
results are in `provenance/` and `reports/`. Verify the release with
`sha256sum -c checksums/SHA256SUMS` after downloading all files.

PDB entries may represent only a domain or an engineered construct and may
contain missing coordinates. Experimental structures have method-dependent
quality and are not necessarily full-length natural proteins. Sequence
isolation under the reported search settings does not establish fold novelty
or absence of exposure during pretraining. The raw PDB archive and training
structure coordinates are not mirrored here.

PDB archive data files are available under CC0 1.0. Cite the original PDB
entries and their structure publications where possible, as well as the
wwPDB/RCSB archive and Foldseek method. Source and code:
<https://www.wwpdb.org/>, <https://www.rcsb.org/>,
<https://github.com/caijihuize/pdb-3di-dataset>.
"""


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--root", type=Path, required=True)
    p.add_argument("--out-dir", type=Path, required=True)
    a = p.parse_args()
    root = a.root.resolve()
    output = a.out_dir.resolve()
    source = load_json(root / "reports/pdb_v1/source_validation.json")
    validation = load_json(root / "reports/pdb_v1/validation.json")
    summary = load_json(root / "reports/pdb_v1/summary.json")
    if source.get("status") != "passed" or validation.get("status") != "passed":
        raise RuntimeError("Source or release validation has not passed")
    counts = summary["split_counts"]
    if counts != validation["split_counts"] or counts["valid"] != 1000 or counts["test"] != 1000:
        raise RuntimeError("Release split counts are incomplete or inconsistent")
    if any(result.get("threshold_violations") != 0
           for result in summary["sequence_audit"].values()):
        raise RuntimeError("Sequence audit is missing or has threshold violations")
    if any(result.get("status") == "missing"
           for result in summary["structure_audit"].values()):
        raise RuntimeError("Structural audit is missing")
    stage = output.with_name(output.name + ".staging")
    if stage.exists():
        shutil.rmtree(stage)
    stage.mkdir(parents=True)
    processed = root / "data/processed/pdb_v1"
    manifests = root / "manifests/pdb_v1"
    reports = root / "reports/pdb_v1"
    for split in ("train", "validation", "test"):
        copy(processed / f"{split}.parquet",
             stage / "data" / f"{split}-00000-of-00001.parquet")
    for split in ("train", "valid", "test"):
        copy(manifests / f"{split}_ids.txt", stage / "ids" / f"{split}_ids.txt")
        copy(manifests / f"{split}_row_index.tsv",
             stage / "provenance" / f"{split}_row_index.tsv")
    for split in ("validation", "test"):
        ids = (manifests / f"{'valid' if split == 'validation' else split}_ids.txt").read_text().splitlines()
        directory = processed / "structures" / split
        files = sorted(directory.glob("*.cif.gz"))
        if len(files) != len(ids):
            raise RuntimeError(f"{split}: {len(files)} structures for {len(ids)} records")
        for path in files:
            copy(path, stage / "structures" / split / path.name)
    for name in ("source_ids.txt", "source_query.json", "cluster_membership.tsv",
                 "excluded.tsv", "split_counts.tsv"):
        copy(manifests / name, stage / "provenance" / name)
    for path in sorted(manifests.glob("downloads_*.tsv")):
        copy(path, stage / "provenance" / "downloads" / path.name)
    for name in ("records.tsv", "filter_counts.tsv", "structure_metadata.tsv"):
        copy(root / "data/interim/candidates" / name, stage / "provenance" / name)
    for name in ("source_validation.json", "validation.json", "summary.json"):
        copy(reports / name, stage / "reports" / name)
    (stage / "README.md").write_text(card(counts, source["source_ids"]))
    checksums = stage / "checksums/SHA256SUMS"
    checksums.parent.mkdir(parents=True, exist_ok=True)
    with checksums.open("w") as handle:
        for path in sorted(stage.rglob("*")):
            if path.is_file() and path != checksums:
                handle.write(f"{digest(path)}  {path.relative_to(stage)}\n")
    if output.exists():
        shutil.rmtree(output)
    stage.rename(output)
    print(f"Prepared {output}")


if __name__ == "__main__":
    main()
