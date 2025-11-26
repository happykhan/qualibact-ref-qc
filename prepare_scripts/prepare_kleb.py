#!/usr/bin/env python3

"""
Generate sample sheets for the Klebsiella reference QC pipeline.

The script scans the FASTA directories inside Genomes_20251008 (or a custom
root) and emits CSV files containing the `sample_id` and absolute FASTA path.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path
from typing import Iterable

SPECIES_FASTA_DIRS = {
    "Kquasipneumoniae": "Kquasipneumoniae_fastas",
    "Kvariicola": "Kvariicola_fastas",
}

VALID_SUFFIXES = (
    ".fa",
    ".fa.gz",
    ".fasta",
    ".fasta.gz",
    ".fna",
    ".fna.gz",
)


def infer_sample_id(fasta_path: Path) -> str:
    """Derive the sample identifier from the FASTA filename."""
    name = fasta_path.name
    if name.lower().endswith(".gz"):
        name = name[:-3]

    lower_name = name.lower()
    for ext in (".fasta", ".fna", ".fa"):
        if lower_name.endswith(ext):
            name = name[: -len(ext)]
            break
    return name


def iter_fastas(fasta_dir: Path) -> Iterable[Path]:
    """Yield FASTA files inside fasta_dir, filtering by allowed suffixes."""
    for candidate in sorted(fasta_dir.iterdir()):
        if not candidate.is_file():
            continue
        if not candidate.name.lower().endswith(VALID_SUFFIXES):
            continue
        yield candidate


def write_samplesheet(fasta_dir: Path, output_path: Path) -> int:
    """Write the CSV samplesheet and return the number of rows produced."""
    fastas = list(iter_fastas(fasta_dir))
    if not fastas:
        raise ValueError(f"No FASTA files found in {fasta_dir}")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["sample_id", "fasta"])
        writer.writeheader()
        for fasta in fastas:
            writer.writerow(
                {
                    "sample_id": infer_sample_id(fasta),
                    "fasta": str(fasta.resolve()),
                }
            )
    return len(fastas)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Create Nextflow sample sheets for Klebsiella references."
    )
    parser.add_argument(
        "--genomes-root",
        default="Genomes_20251008",
        help="Root directory holding FASTA folders (default: %(default)s).",
    )
    parser.add_argument(
        "--out-dir",
        default=None,
        help="Optional output directory for the CSV files (defaults to genomes root).",
    )
    args = parser.parse_args()

    genomes_root = Path(args.genomes_root).expanduser().resolve()
    out_dir = Path(args.out_dir).expanduser().resolve() if args.out_dir else genomes_root

    if not genomes_root.exists():
        parser.error(f"Genomes root not found: {genomes_root}")

    written = 0
    for species, subdir in SPECIES_FASTA_DIRS.items():
        fasta_dir = genomes_root / subdir
        if not fasta_dir.exists():
            print(
                f"[WARN] Missing FASTA directory for {species}: {fasta_dir}",
                file=sys.stderr,
            )
            continue

        output_csv = out_dir / f"{species}_samplesheet.csv"
        try:
            count = write_samplesheet(fasta_dir, output_csv)
            print(f"Wrote {count} rows to {output_csv}")
            written += 1
        except ValueError as exc:
            print(f"[WARN] {exc}", file=sys.stderr)

    if written == 0:
        sys.exit("No sample sheets were created. Check the input directories.")


if __name__ == "__main__":
    main()
