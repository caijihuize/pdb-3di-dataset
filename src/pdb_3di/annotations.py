"""Create chain-level quality, provenance and external mapping annotations."""

from __future__ import annotations

import argparse
import concurrent.futures
import csv
import gzip
import math
from collections import Counter, defaultdict
from pathlib import Path

from Bio.PDB.MMCIF2Dict import MMCIF2Dict

from pdb_3di.foldseek_keys import components


MISSING = {"", ".", "?"}
FIELDS = [
    "id", "pdb_id", "model_id", "chain_id", "deposition_date", "release_date",
    "r_work", "r_free", "assembly_count", "assembly_descriptions",
    "observed_ca_residues", "complete_backbone_residues", "backbone_complete_fraction",
    "internal_missing_residues", "terminal_missing_residues", "mean_ca_occupancy",
    "min_ca_occupancy", "mean_ca_b_factor", "alternate_location_residues",
    "three_di_entropy_bits", "three_di_max_token_fraction", "three_di_unique_tokens",
    "native_uniprot_ids", "sifts_uniprot_ids", "sifts_pfam_ids", "sifts_cath_ids",
]


def values(data: dict, name: str) -> list[str]:
    value = data.get(name, [])
    return value if isinstance(value, list) else [value]


def scalar(data: dict, name: str) -> str:
    for value in values(data, name):
        if value not in MISSING:
            return value
    return ""


def number(data: dict, name: str) -> str:
    value = scalar(data, name)
    try:
        return f"{float(value):.4f}"
    except ValueError:
        return ""


def sequence_stats(sequence: str) -> tuple[str, str, str]:
    counts = Counter(sequence)
    length = len(sequence)
    entropy = -sum((n / length) * math.log2(n / length) for n in counts.values())
    return f"{entropy:.6f}", f"{max(counts.values()) / length:.6f}", str(len(counts))


def requested_model(model_id: str) -> str:
    return model_id or "1"


def parse_one(task: tuple[Path, list[dict]]) -> tuple[list[dict], tuple | None]:
    path, targets = task
    try:
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
            data = MMCIF2Dict(handle)
        wanted = {(requested_model(row["model_id"]), row["chain_id"]) for row in targets}
        residues = defaultdict(lambda: {"atoms": set(), "alts": set(), "ca_occ": [],
                                        "ca_b": [], "label_seq": None})
        columns = [
            values(data, "_atom_site.pdbx_PDB_model_num"),
            values(data, "_atom_site.auth_asym_id"),
            values(data, "_atom_site.label_seq_id"),
            values(data, "_atom_site.auth_seq_id"),
            values(data, "_atom_site.pdbx_PDB_ins_code"),
            values(data, "_atom_site.label_atom_id"),
            values(data, "_atom_site.label_alt_id"),
            values(data, "_atom_site.occupancy"),
            values(data, "_atom_site.B_iso_or_equiv"),
            values(data, "_atom_site.group_PDB"),
        ]
        for model, chain, label_seq, auth_seq, ins, atom, alt, occ, bfac, group in zip(*columns):
            if group != "ATOM" or (model, chain) not in wanted:
                continue
            key = (model, chain, label_seq if label_seq not in MISSING else f"{auth_seq}:{ins}")
            residue = residues[key]
            residue["atoms"].add(atom)
            if alt not in MISSING:
                residue["alts"].add(alt)
            if label_seq not in MISSING:
                try:
                    residue["label_seq"] = int(label_seq)
                except ValueError:
                    pass
            if atom == "CA":
                try:
                    residue["ca_occ"].append(float(occ))
                except ValueError:
                    pass
                try:
                    residue["ca_b"].append(float(bfac))
                except ValueError:
                    pass

        native = defaultdict(set)
        refs = {}
        for ref_id, db_name, accession in zip(values(data, "_struct_ref.id"),
                                               values(data, "_struct_ref.db_name"),
                                               values(data, "_struct_ref.pdbx_db_accession")):
            if db_name in {"UNP", "UNIPROT"} and accession not in MISSING:
                refs[ref_id] = accession
        for ref_id, accession, strands in zip(values(data, "_struct_ref_seq.ref_id"),
                                               values(data, "_struct_ref_seq.pdbx_db_accession"),
                                               values(data, "_struct_ref_seq.pdbx_strand_id")):
            accession = accession if accession not in MISSING else refs.get(ref_id, "")
            if accession:
                for chain in strands.replace(";", ",").split(","):
                    native[chain.strip()].add(accession)

        revisions = sorted(x for x in values(data, "_pdbx_audit_revision_history.revision_date")
                           if x not in MISSING)
        entry = {
            "deposition_date": scalar(data, "_pdbx_database_status.recvd_initial_deposition_date"),
            "release_date": revisions[0] if revisions else "",
            "r_work": number(data, "_refine.ls_R_factor_R_work"),
            "r_free": number(data, "_refine.ls_R_factor_R_free"),
            "assembly_count": str(len(set(x for x in values(data, "_pdbx_struct_assembly.id")
                                          if x not in MISSING))),
            "assembly_descriptions": ";".join(sorted(set(
                x for x in values(data, "_pdbx_struct_assembly.oligomeric_details") if x not in MISSING))),
        }
        output = []
        for target in targets:
            model = requested_model(target["model_id"])
            chain_rows = [r for (m, c, _), r in residues.items() if m == model and c == target["chain_id"]]
            ca_rows = [r for r in chain_rows if "CA" in r["atoms"]]
            complete = sum({"N", "CA", "C", "O"}.issubset(r["atoms"]) for r in chain_rows)
            positions = sorted({r["label_seq"] for r in ca_rows if r["label_seq"] is not None})
            deposited = target["deposited_length"]
            internal = (positions[-1] - positions[0] + 1 - len(positions)) if positions else ""
            terminal = ((positions[0] - 1) + max(0, deposited - positions[-1])) if positions and deposited else ""
            occupancies = [max(r["ca_occ"]) for r in ca_rows if r["ca_occ"]]
            bfactors = [sum(r["ca_b"]) / len(r["ca_b"]) for r in ca_rows if r["ca_b"]]
            row = dict(target)
            row.update(entry)
            row.update({
                "observed_ca_residues": str(len(ca_rows)),
                "complete_backbone_residues": str(complete),
                "backbone_complete_fraction": f"{complete / len(ca_rows):.6f}" if ca_rows else "",
                "internal_missing_residues": str(internal),
                "terminal_missing_residues": str(terminal),
                "mean_ca_occupancy": f"{sum(occupancies) / len(occupancies):.6f}" if occupancies else "",
                "min_ca_occupancy": f"{min(occupancies):.6f}" if occupancies else "",
                "mean_ca_b_factor": f"{sum(bfactors) / len(bfactors):.6f}" if bfactors else "",
                "alternate_location_residues": str(sum(bool(r["alts"]) for r in ca_rows)),
                "native_uniprot_ids": ";".join(sorted(native[target["chain_id"]])),
            })
            output.append(row)
        return output, None
    except Exception as exc:
        return [], (path.name, type(exc).__name__, str(exc).replace("\t", " ").replace("\n", " "))


def csv_rows(path: Path):
    with path.open(newline="", encoding="utf-8-sig") as handle:
        lines = (line for line in handle if not line.startswith("#"))
        yield from csv.DictReader(lines)


def external_mappings(directory: Path, wanted: set[tuple[str, str]]) -> dict:
    output = defaultdict(lambda: {"sifts_uniprot_ids": set(), "sifts_pfam_ids": set(),
                                  "sifts_cath_ids": set()})
    specs = [
        ("pdb_chain_uniprot.csv", "sifts_uniprot_ids", ("SP_PRIMARY", "UNIPROT")),
        ("pdb_chain_pfam.csv", "sifts_pfam_ids", ("PFAM_ID",)),
        ("pdb_chain_cath_uniprot.csv", "sifts_cath_ids", ("CATH_ID",)),
    ]
    for filename, field, candidates in specs:
        path = directory / filename
        if not path.is_file():
            continue
        for row in csv_rows(path):
            normalized = {k.upper(): v for k, v in row.items() if k is not None}
            key = (normalized.get("PDB", "").upper(), normalized.get("CHAIN", ""))
            if key not in wanted:
                continue
            value = next((normalized.get(name, "") for name in candidates if normalized.get(name)), "")
            if value not in MISSING:
                output[key][field].add(value)
    return output


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--records", type=Path, required=True)
    parser.add_argument("--raw-dir", type=Path, required=True)
    parser.add_argument("--sifts-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--errors", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()
    grouped = defaultdict(list)
    wanted = set()
    with args.records.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            pdb_id, model_id, chain_id = components(row["id"])
            entropy, max_fraction, unique = sequence_stats(row["3di"])
            target = {"id": row["id"], "pdb_id": pdb_id, "model_id": model_id,
                      "chain_id": chain_id, "deposited_length": int(row["deposited_length"] or 0),
                      "three_di_entropy_bits": entropy,
                      "three_di_max_token_fraction": max_fraction,
                      "three_di_unique_tokens": unique}
            grouped[pdb_id].append(target)
            wanted.add((pdb_id, chain_id))
    mappings = external_mappings(args.sifts_dir, wanted)
    tasks = []
    for pdb_id, targets in sorted(grouped.items()):
        path = args.raw_dir / pdb_id[1:3].lower() / f"{pdb_id.lower()}.cif.gz"
        tasks.append((path, targets))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.errors.parent.mkdir(parents=True, exist_ok=True)
    part = args.output.with_suffix(args.output.suffix + ".part")
    errors = []
    count = 0
    with part.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        with concurrent.futures.ProcessPoolExecutor(max_workers=args.workers) as pool:
            for index, (rows, error) in enumerate(pool.map(parse_one, tasks, chunksize=1), 1):
                if error:
                    errors.append(error)
                for row in rows:
                    mapped = mappings[(row["pdb_id"], row["chain_id"])]
                    for field in ("sifts_uniprot_ids", "sifts_pfam_ids", "sifts_cath_ids"):
                        row[field] = ";".join(sorted(mapped[field]))
                    writer.writerow({field: row.get(field, "") for field in FIELDS})
                    count += 1
                if index % 10000 == 0:
                    print(f"annotated structures {index}/{len(tasks)} records={count}", flush=True)
    part.replace(args.output)
    with args.errors.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["file", "error_type", "message"])
        writer.writerows(errors)
    print(f"records={count} errors={len(errors)}")


if __name__ == "__main__":
    main()
