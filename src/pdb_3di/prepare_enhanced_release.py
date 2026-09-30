"""Stage validated enhanced versions for an additive Hugging Face upload."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from pathlib import Path

from pdb_3di.export_versions import SPLITS, VERSIONS


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    validation = json.loads((args.root / "reports/pdb_v2/validation.json").read_text())
    if validation.get("status") != "passed":
        raise RuntimeError("Enhanced version validation has not passed")
    stage = args.output.with_name(args.output.name + ".staging")
    if stage.exists():
        shutil.rmtree(stage)
    stage.mkdir(parents=True)
    descriptions = {
        "pdb_v1_1": "The pdb_v1 split with chain quality, date, 3Di complexity and SIFTS annotations.",
        "pdb_v2_sequence": "UniProt-aware 30% sequence-isolated split with inverse component-size sampling metadata.",
        "pdb_v2_structure": "Sequence isolation plus Foldseek TM-score and CATH topology connected components.",
        "pdb_v2_time": "Chronological release split with sequence- and structure-novelty labels.",
    }
    index = ["# Enhanced experimental PDB AA–3Di versions", "", "Validated additive versions of the dataset:", ""]
    for version in VERSIONS:
        counts = validation["versions"][version]
        index.append(f"- `{version}`: {descriptions[version]} "
                     f"({counts['train']:,}/{counts['valid']:,}/{counts['test']:,} train/validation/test)")
        target = stage / version
        (target / "data").mkdir(parents=True)
        (target / "ids").mkdir(parents=True)
        for split in SPLITS:
            parquet_name = "validation" if split == "valid" else split
            shutil.copy2(args.root / "data/processed" / version / f"{parquet_name}.parquet",
                         target / "data" / f"{parquet_name}.parquet")
            release_split = "validation" if split == "valid" else split
            shutil.copy2(args.root / "manifests" / version / f"{split}_ids.txt",
                         target / "ids" / f"{release_split}_ids.txt")
            shutil.copy2(args.root / "manifests" / version / f"{split}_structure_ids.tsv",
                         target / "ids" / f"{release_split}_structure_ids.tsv")
        shutil.copy2(args.root / "data/processed" / version / "schema.txt", target / "schema.txt")
        for name in ("cluster_membership.tsv", "split_counts.tsv", "record_annotations.tsv",
                     "novelty.tsv", "VERSION.txt"):
            source = args.root / "manifests" / version / name
            if source.is_file():
                (target / "provenance").mkdir(exist_ok=True)
                shutil.copy2(source, target / "provenance" / name)
        (target / "README.md").write_text(
            f"# {version}\n\n{descriptions[version]}\n\n"
            "Parquet row order matches the files in `ids/`. Each structure mapping records "
            "the PDB entry, model, author chain and exact Foldseek database key.\n")
    (stage / "README.md").write_text("\n".join(index) + "\n")
    (stage / "provenance").mkdir()
    shutil.copy2(args.root / "configs/pdb_v2.yaml", stage / "provenance/pdb_v2.yaml")
    shutil.copy2(args.root / "reports/pdb_v2/validation.json", stage / "provenance/validation.json")
    shutil.copy2(args.root / "data/raw/mappings/SHA256SUMS", stage / "provenance/mapping_SHA256SUMS")
    checksums = stage / "SHA256SUMS"
    with checksums.open("w") as handle:
        for path in sorted(stage.rglob("*")):
            if path.is_file() and path != checksums:
                handle.write(f"{digest(path)}  {path.relative_to(stage)}\n")
    if args.output.exists():
        shutil.rmtree(args.output)
    stage.rename(args.output)


if __name__ == "__main__":
    main()
