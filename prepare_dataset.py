"""
Download NCBI reference genomes for a species and build a samplesheet.

The script wraps the NCBI Datasets CLI (`datasets`) so we do not have to hand-roll
API requests. Provide a species name (taxon) and it will:
1. Download all reference assemblies for that species (optionally filtered by
   assembly level).
2. Copy the FASTA files into a directory named after the species.
3. Emit a samplesheet (`samplesheet.csv`) with `sample_id` and `assembly_path`
   columns that can be used as input for the Nextflow workflow.

Dependencies are provided via Pixi; see the project `pixi.toml`.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import shutil
import subprocess
import tempfile
import zipfile
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence


DATASETS_CMD = ["pixi", "run", "datasets"]


@dataclass
class AssemblyRecord:
    sample_id: str
    source_fasta: Path
    dest_fasta: Path


def slugify(value: str) -> str:
    """Return a filesystem-friendly slug."""
    value = re.sub(r"[^A-Za-z0-9._-]+", "_", value.strip())
    return value.strip("_") or "sample"


def run_command(cmd: Sequence[str]) -> None:
    """Run a command and raise a helpful error on failure."""
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"Command failed ({result.returncode}): {' '.join(cmd)}\n"
            f"stdout:\n{result.stdout}\n\nstderr:\n{result.stderr}"
        )


def check_datasets_cli_available() -> None:
    if shutil.which(DATASETS_CMD) is None:
        raise RuntimeError(
            "NCBI Datasets CLI ('datasets') is not available. "
            "Install it with `pixi install` (see pixi.toml) or ensure it is on your PATH."
        )


def ensure_output_dir(target: Path, force: bool) -> None:
    if target.exists():
        if any(target.iterdir()) and not force:
            raise FileExistsError(
                f"Output directory {target} already exists and is not empty. "
                "Use --force to overwrite."
            )
        if force:
            shutil.rmtree(target)
    target.mkdir(parents=True, exist_ok=True)


def download_package(
    species: str, assembly_levels: Iterable[str], reference_only: bool, tmpdir: Path
) -> Path:
    """
    Use NCBI Datasets CLI to download genomes for a taxon.

    Returns the path to the unzipped package root.
    """
    zip_path = tmpdir / "ncbi_dataset.zip"
    cmd = [
        DATASETS_CMD*,
        "download",
        "genome",
        "taxon",
        species,
        "--include",
        "genome",
        "--filename",
        str(zip_path),
    ]
    if reference_only:
        cmd.append("--reference")
    if assembly_levels:
        cmd += ["--assembly-level", ",".join(assembly_levels)]
    run_command(cmd)

    extract_root = tmpdir / "package"
    with zipfile.ZipFile(zip_path, "r") as zf:
        zf.extractall(extract_root)
    return extract_root


def load_metadata(package_root: Path) -> Dict[str, dict]:
    """
    Parse assembly_data_report.jsonl (if present) to get metadata keyed by accession.
    """
    report = package_root / "ncbi_dataset" / "data" / "assembly_data_report.jsonl"
    metadata: Dict[str, dict] = {}
    if not report.exists():
        return metadata

    with report.open() as handle:
        for line in handle:
            if not line.strip():
                continue
            try:
                entry = json.loads(line)
            except json.JSONDecodeError:
                continue
            accession = entry.get("assembly_accession")
            if accession:
                metadata[accession] = entry
    return metadata


def derive_sample_id(accession: str, meta: Dict[str, dict]) -> str:
    """Use accession plus strain (if present) to build a readable sample_id."""
    entry = meta.get(accession, {})
    strain = (
        entry.get("organism", {})
        .get("infraspecific_names", {})
        .get("strain")
    )
    if strain:
        return slugify(f"{accession}_{strain}")
    return accession


def find_genome_fasta(assembly_dir: Path) -> Optional[Path]:
    """
    Locate the genomic FASTA within an assembly directory.
    """
    candidates = list(
        assembly_dir.glob("*_genomic.fna")
    ) or list(assembly_dir.glob("*_genomic.fasta")) or list(
        assembly_dir.glob("*.fna")
    ) or list(
        assembly_dir.glob("*.fasta")
    )
    if not candidates:
        return None
    return candidates[0]


def collect_assemblies(package_root: Path, dest_dir: Path) -> List[AssemblyRecord]:
    data_dir = package_root / "ncbi_dataset" / "data"
    metadata = load_metadata(package_root)
    records: List[AssemblyRecord] = []

    for assembly_dir in sorted(p for p in data_dir.iterdir() if p.is_dir()):
        accession = assembly_dir.name
        fasta = find_genome_fasta(assembly_dir)
        if not fasta:
            continue
        sample_id = derive_sample_id(accession, metadata)
        dest_fasta = dest_dir / f"{accession}.fna"
        records.append(AssemblyRecord(sample_id, fasta, dest_fasta))
    return records


def copy_assemblies(records: Sequence[AssemblyRecord]) -> None:
    for record in records:
        record.dest_fasta.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(record.source_fasta, record.dest_fasta)


def write_samplesheet(records: Sequence[AssemblyRecord], path: Path) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["sample_id", "assembly_path"])
        for record in records:
            writer.writerow([record.sample_id, record.dest_fasta.resolve()])


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Download NCBI reference genomes for a species and build a samplesheet."
    )
    parser.add_argument(
        "species",
        help="Species/taxon name (e.g. 'Salmonella enterica').",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path.cwd(),
        help="Directory where the species folder and samplesheet will be written.",
    )
    parser.add_argument(
        "--assembly-levels",
        default="complete,chromosome",
        help="Comma-separated assembly levels to request (e.g. complete,chromosome,scaffold).",
    )
    parser.add_argument(
        "--include-non-reference",
        action="store_true",
        help="Include non-reference assemblies (reference assemblies only by default).",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite an existing species directory.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    assembly_levels = [
        level.strip()
        for level in args.assembly_levels.split(",")
        if level.strip()
    ]

    species_slug = slugify(args.species)
    species_dir = args.outdir / species_slug
    assemblies_dir = species_dir / "assemblies"
    samplesheet_path = species_dir / "samplesheet.csv"

    check_datasets_cli_available()
    ensure_output_dir(species_dir, args.force)

    with tempfile.TemporaryDirectory() as tmp:
        tmpdir = Path(tmp)
        package_root = download_package(
            args.species,
            assembly_levels,
            reference_only=not args.include_non_reference,
            tmpdir=tmpdir,
        )
        records = collect_assemblies(package_root, assemblies_dir)
        if not records:
            raise RuntimeError("No assemblies found in the downloaded package.")
        copy_assemblies(records)
        write_samplesheet(records, samplesheet_path)

    print(f"Assemblies written to: {assemblies_dir.resolve()}")
    print(f"Samplesheet written to: {samplesheet_path.resolve()}")


if __name__ == "__main__":
    main()
