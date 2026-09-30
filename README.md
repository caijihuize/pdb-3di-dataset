# Experimental PDB AA–3Di dataset

This project builds paired amino acid (AA) and Foldseek 3Di sequences from
experimentally determined protein structures in the wwPDB PDB Core Archive.
It follows the split and audit design of the sibling `sp-dataset` project.
The source inventory is obtained from RCSB's Search API and frozen in
`manifests/pdb_v1/source_ids.txt` with the exact query and retrieval date.
Coordinates are downloaded as PDBx/mmCIF.gz.

## Scope of version 1

- Released experimental PDB entries with at least one protein polymer entity.
- X-ray structures with reported resolution at most 3.0 Å and electron
  microscopy structures with reported resolution at most 4.0 Å.
- Protein chains with at least 30 Foldseek-observed residues, at most 20% unknown
  amino acids, equal AA/3Di lengths, and valid 3Di symbols.
- NMR and integrative/hybrid models are not in the primary version 1 dataset.
- One record is one Foldseek chain, identified as `<PDB_ID>_<chain>`. Entries
  with multiple coordinate models use `<PDB_ID>_MODEL_<model>_<chain>`.
  `records.tsv` retains the PDB entry, author chain ID, mapped polymer entity,
  method, resolution, sequence coverage and both sequences.
- The evaluation mmCIF is the full deposited entry. Its target chain is specified
  by the record ID and metadata; multiple proteins may occur in one file.

The AA sequence is the chain sequence extracted by Foldseek from the coordinates,
not the full deposited SEQRES sequence. Missing residues therefore do not create
AA-only positions. The stored coverage is informative; complex entries may have
missing or ambiguous entity mappings and then have empty coverage.

## Build

Create the environment with `mamba env create -f environment.yml`, then run
from this repository root:

```bash
bash scripts/run.sh discover
bash scripts/run.sh download
bash scripts/run.sh features
bash scripts/run.sh candidates
bash scripts/run.sh cluster
bash scripts/run.sh split
bash scripts/finalize.sh
```

Set `PDB3DI_PYTHON`, `FOLDSEEK_BIN`, `MMSEQS_BIN`, `THREADS` and
`DOWNLOAD_WORKERS` when needed. `DOWNLOAD_LIMIT=N` limits the download stage
to the first N IDs for a pilot. Stages are restartable, and the download stage
checks existing gzip files and records SHA-256 checksums in `downloads.tsv`.
The full build requires substantial storage, network transfer and batch compute.

For the current Slurm cluster, `submit_download.sh` downloads the frozen ID
list in 16 array shards, with at most two active shards. `submit_build.sh`
verifies every downloaded entry before feature extraction and clustering;
`submit_finalize.sh` performs audits, export and validation. Submit with
`afterok` dependencies, for example:

```bash
DOWNLOAD_JOB=$(sbatch --parsable scripts/submit_download.sh)
BUILD_JOB=$(sbatch --parsable --dependency="afterok:${DOWNLOAD_JOB}" scripts/submit_build.sh)
sbatch --dependency="afterok:${BUILD_JOB}" scripts/submit_finalize.sh
```

The stage scripts also run interactively for small pilots. Do not run
`finalize.sh` on a partial source inventory with release-sized validation/test
settings.

## Release output

`data/processed/pdb_v1/` contains `train.parquet`, `validation.parquet` and
`test.parquet`, with columns `id`, `sequence_3di`, `sequence_aa`, plus paired
FASTA files and validation/test mmCIF structures. `manifests/pdb_v1/` contains
source IDs, checksums, split IDs, group membership and exclusion reasons.
`reports/pdb_v1/` contains validation and audit summaries.
The published `ids/` directory contains one stable ID per Parquet row plus
`*_structure_ids.tsv` mappings with the PDB entry, model, author chain and
build-specific Foldseek database key. Stable IDs should be used to reference
records; the numeric Foldseek keys are included for exact build provenance.
After all checks pass, `scripts/finalize.sh` also stages a Hugging Face dataset
repository in `release/pdb_v1/`, including a dataset card, three Parquet splits,
evaluation mmCIFs, provenance, reports and SHA-256 checksums. It does not include
the downloaded source archive or intermediate databases. Publish the staged
folder with `hf upload caijihuize/pdb-3di-dataset release/pdb_v1 . --repo-type dataset`
after checking its contents and target repository settings.
For this build, `scripts/submit_publish_hf.sh` can run after the finalization
job succeeds. It verifies every staged file against `checksums/SHA256SUMS`,
then creates the public dataset repository if needed and uploads the release.

MMseqs2 clustering uses 30% identity and 80% bidirectional coverage. PDB-entry
members and exact AA duplicates are kept in the same group. Each evaluation
split selects one chain from each of 1,000 groups. A full cross-split MMseqs2
search repairs detected threshold conflicts before release. Foldseek structural
neighbors are reported and do not determine split membership.

The PDB archive is updated regularly. Reproducibility requires preserving the
source query, retrieval date, entry ID list, downloaded file checksums, tool
versions, parameters and final split IDs. A live archive inventory is not a
permanent historical snapshot. Historical wwPDB snapshots are an alternative
when exact coordinate version pinning is required before download.

## Source and reuse

PDB Core Archive: <https://www.wwpdb.org/>. Search and file APIs:
<https://search.rcsb.org/> and
<https://www.rcsb.org/docs/programmatic-access/file-download-services>.
The source data and original code/documentation have their respective terms;
see `LICENSE` for this project's code license. Core split/audit routines were
adapted from the sibling MIT-licensed `sp-dataset` project.
