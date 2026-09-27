"""Read experimental mmCIF metadata and pair Foldseek AA/3Di chains."""

from __future__ import annotations

import argparse
import concurrent.futures
import csv
import gzip
import re
from collections import Counter
from pathlib import Path

from Bio.PDB.MMCIF2Dict import MMCIF2Dict

AA = set("ACDEFGHIKLMNPQRSTVWY")
DI = set("algvsredtipkfqnymhwc")


def first(values, default=""):
    if values is None:
        return default
    return values[0] if isinstance(values, list) else values


def best_float(values) -> float | None:
    if values is None:
        return None
    if not isinstance(values, list):
        values = [values]
    parsed = []
    for value in values:
        try:
            parsed.append(float(value))
        except (TypeError, ValueError):
            pass
    return min(parsed) if parsed else None


def sequence_length(seq: str) -> int:
    return len(re.sub(r"\s+", "", re.sub(r"\([^)]*\)", "X", seq)))


def metadata(path: Path) -> dict:
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        data = MMCIF2Dict(handle)
    method = ";".join(data.get("_exptl.method", []))
    if "X-RAY DIFFRACTION" in method:
        category = "xray"
        resolution = best_float(data.get("_refine.ls_d_res_high"))
    elif "ELECTRON MICROSCOPY" in method:
        category = "em"
        resolution = best_float(data.get("_em_3d_reconstruction.resolution"))
    else:
        category = "other"
        resolution = None
    entity_length = {}
    for entity, poly_type, seq in zip(data.get("_entity_poly.entity_id", []),
                                      data.get("_entity_poly.type", []),
                                      data.get("_entity_poly.pdbx_seq_one_letter_code_can", [])):
        if poly_type == "polypeptide(L)":
            entity_length[entity] = sequence_length(seq)
    auth_entity = {}
    for auth, entity, atom, group in zip(data.get("_atom_site.auth_asym_id", []),
                                         data.get("_atom_site.label_entity_id", []),
                                         data.get("_atom_site.label_atom_id", []),
                                         data.get("_atom_site.group_PDB", [])):
        if atom == "CA" and group == "ATOM" and entity in entity_length:
            auth_entity[auth.upper()] = entity
    return {"pdb_id": path.name[:4].upper(), "method": method, "category": category,
            "resolution": resolution, "entity_lengths": entity_length,
            "auth_entity": auth_entity}


def safe_metadata(path: Path) -> tuple[dict | None, tuple | None]:
    try:
        return metadata(path), None
    except Exception as exc:
        return None, (path.name, type(exc).__name__, str(exc).replace("\t", " "))


def fasta(path: Path) -> dict[str, str]:
    out = {}
    key, chunks = None, []
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                if key is not None:
                    if key in out:
                        raise ValueError(f"Duplicate Foldseek ID: {key}")
                    out[key] = "".join(chunks)
                key = line[1:].split()[0].upper()
                chunks = []
            elif line:
                chunks.append(line)
    if key is not None:
        if key in out:
            raise ValueError(f"Duplicate Foldseek ID: {key}")
        out[key] = "".join(chunks)
    return out


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--raw-dir", type=Path, required=True)
    p.add_argument("--aa", type=Path, required=True)
    p.add_argument("--di", type=Path, required=True)
    p.add_argument("--out-dir", type=Path, required=True)
    p.add_argument("--min-length", type=int, default=30)
    p.add_argument("--max-unknown-fraction", type=float, default=0.20)
    p.add_argument("--max-xray-resolution", type=float, default=3.0)
    p.add_argument("--max-em-resolution", type=float, default=4.0)
    p.add_argument("--workers", type=int, default=4)
    a = p.parse_args()
    a.out_dir.mkdir(parents=True, exist_ok=True)
    aa, di = fasta(a.aa), fasta(a.di)
    meta = {}
    metadata_errors = []
    paths = sorted(a.raw_dir.rglob("*.cif.gz"))
    with concurrent.futures.ProcessPoolExecutor(max_workers=a.workers) as pool:
        for i, (row, error) in enumerate(pool.map(safe_metadata, paths, chunksize=16), 1):
            if row is not None:
                meta[row["pdb_id"]] = row
            if error is not None:
                metadata_errors.append(error)
            if i % 10000 == 0:
                print(f"parsed metadata {i}/{len(paths)}", flush=True)
    with (a.out_dir / "structure_metadata.tsv").open("w", newline="") as handle:
        w = csv.writer(handle, delimiter="\t", lineterminator="\n")
        w.writerow(["pdb_id", "method", "category", "resolution"])
        for pid, row in sorted(meta.items()):
            w.writerow([pid, row["method"], row["category"], row["resolution"] or ""])
    with (a.out_dir / "metadata_errors.tsv").open("w", newline="") as handle:
        w = csv.writer(handle, delimiter="\t", lineterminator="\n")
        w.writerow(["file", "error_type", "message"])
        w.writerows(metadata_errors)
    excluded = []
    kept = []
    for key in sorted(set(aa) | set(di)):
        pdb_id, _, chain = key.partition("_")
        info = meta.get(pdb_id)
        reason = None
        if key not in aa:
            reason = "missing_aa"
        elif key not in di:
            reason = "missing_3di"
        elif not chain or info is None:
            reason = "missing_structure_metadata"
        else:
            aaseq, diseq = aa[key].upper(), di[key].lower()
            resolution = info["resolution"]
            if len(aaseq) != len(diseq):
                reason = "length_mismatch"
            elif len(aaseq) < a.min_length:
                reason = "too_short"
            elif not aaseq or sum(x not in AA for x in aaseq) / len(aaseq) > a.max_unknown_fraction:
                reason = "too_many_unknown_aa"
            elif any(x not in DI for x in diseq):
                reason = "invalid_3di"
            elif info["category"] == "other":
                reason = "method_excluded"
            elif resolution is None or resolution <= 0:
                reason = "missing_resolution"
            elif info["category"] == "xray" and resolution > a.max_xray_resolution:
                reason = "low_xray_resolution"
            elif info["category"] == "em" and resolution > a.max_em_resolution:
                reason = "low_em_resolution"
            if reason is None:
                entity = info["auth_entity"].get(chain, "")
                full_len = info["entity_lengths"].get(entity, 0)
                coverage = min(1.0, len(aaseq) / full_len) if full_len else None
                kept.append((key, pdb_id, chain, entity, info["method"], info["category"],
                             resolution, full_len, coverage, aaseq, diseq))
        if reason:
            excluded.append((key, reason))
    with (a.out_dir / "records.tsv").open("w", newline="") as handle, \
         (a.out_dir / "aa.filtered.fasta").open("w") as aaout, \
         (a.out_dir / "3di.filtered.fasta").open("w") as diout:
        w = csv.writer(handle, delimiter="\t", lineterminator="\n")
        w.writerow(["id", "pdb_id", "chain_id", "entity_id", "method", "method_category",
                    "resolution", "deposited_length", "coverage", "length", "quality_score", "aa", "3di"])
        for key, pdb_id, chain, entity, method, category, resolution, full_len, coverage, aaseq, diseq in kept:
            w.writerow([key, pdb_id, chain, entity, method, category, f"{resolution:.3f}",
                        full_len or "", f"{coverage:.4f}" if coverage is not None else "",
                        len(aaseq), f"{1 / resolution:.6f}", aaseq, diseq])
            aaout.write(f">{key}\n{aaseq}\n")
            diout.write(f">{key}\n{diseq}\n")
    with (a.out_dir / "excluded.tsv").open("w", newline="") as handle:
        w = csv.writer(handle, delimiter="\t", lineterminator="\n")
        w.writerow(["id", "reason"])
        w.writerows(excluded)
    counts = Counter(reason for _, reason in excluded)
    (a.out_dir / "filter_counts.tsv").write_text(
        f"downloaded_structures\t{len(meta)}\nfoldseek_aa_chains\t{len(aa)}\n"
        f"foldseek_3di_chains\t{len(di)}\nkept\t{len(kept)}\nexcluded\t{len(excluded)}\n" +
        "".join(f"{name}\t{n}\n" for name, n in sorted(counts.items())))
    print(f"structures={len(meta)} kept_chains={len(kept)} excluded_chains={len(excluded)} metadata_errors={len(metadata_errors)}")


if __name__ == "__main__":
    main()
